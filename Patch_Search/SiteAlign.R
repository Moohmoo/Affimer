#!/usr/local/bin/Rscript
library(bio3d)
library(argparser, quietly=TRUE)
#
#
#
MAXALIGN=200 # safety parameter
#patchsearch_AA=c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","HSD","HSE","HSP")

#
#
#
p=arg_parser("A program to superimpose a patch onto proteins")
p=add_argument(p, "alfile",  help="input align file")
p=add_argument(p, "peptide",  help="peptide chain in query")
p=add_argument(p, "--outprefix", default="target",  help="file prefix")
p=add_argument(p, "--minscore", default="0.",  help="min. score")
p=add_argument(p, "--maxal", default=MAXALIGN,  help="maximum number of alignments")
p=add_argument(p, "--apofile", default="",  help="apo pdb file")
p=add_argument(p, "--column", default="10",help="score column number")

argv = parse_args(p)
alfile=argv$alfile
prefix=argv$outprefix
minscore=as.double(argv$minscore)
maxal=as.integer(argv$maxal)
pepchain=argv$peptide
apofile=argv$apofile
ncolumn=as.integer(argv$column)+1
if (!file.exists(alfile)) {
   message("no align file")
   q()
}


L=readLines(alfile)
Lindex=grep(">",L)
n=length(L)
Lindex=append(Lindex,n+1)
n=length(Lindex)
noal=which(Lindex[2:(n+1)]-Lindex[1:n]==1)
if (length(noal)>0) Lindex=Lindex[-noal]
if (maxal>0) nbal=min(maxal,MAXALIGN)
nbal=min(nbal, length(Lindex))
#nbal=length(Lindex)
Lindex=Lindex[1:(nbal-1)]
cat("nbal=",nbal,"\n")

al_num=1
num_target=0
apo=NULL
if (file.exists(apofile)) apo=suppressWarnings(bio3d::read.pdb(apofile))

for (k in Lindex) {
    num_target=num_target+1
    x=unlist(strsplit(L[k]," "))
    query_file=x[2]
    cat("query_file",query_file, pepchain, "\n")
    #!! si pepchain dans nom de fichier (patch adhoc)
    pchain=pepchain
    if (pepchain==":") {
     pos=gregexpr(":",query_file)[[1]][1]
     pchain=substring(query_file,pos+1,pos+1)
    }
    #!!
    query=suppressWarnings(bio3d::read.pdb(query_file))
    # 
    target_file=x[3]
    cat("target_file=", target_file, "\n")
    cl_len=as.integer(x[7])
    cat(ncolumn, cl_len, x[ncolumn], minscore,"\n")
    if (as.double(x[ncolumn])<minscore) next
    target=suppressWarnings(bio3d::read.pdb(target_file))
    target=trim.pdb(target, string="protein")
    x=unlist(strsplit(L[k+1]," "))
    lq=unlist(strsplit(x,split=":"))
    n=length(lq)
    numq=as.integer(lq[seq(1,n,4)])
    insertq=lq[seq(2,n,4)]
    namq=lq[seq(3,n,4)]
    chq=lq[seq(4,n,4)]
    query_chains=unique(query$atom$chain)
    if (pchain=="?") {
       pchain=setdiff(query_chains, chq)
       cat("pchain=",pchain,"\n")
    }
    query.ind=atom.select(query, resno=numq, insert=insertq, elety=namq, chain=chq)
    x=unlist(strsplit(L[k+2]," "))
    lt=unlist(strsplit(x,split=":"))
    n=length(lt)
    numt=as.integer(lt[seq(1,n,4)])
    insertt=lt[seq(2,n,4)]
    namt=lt[seq(3,n,4)]
    cht=lt[seq(4,n,4)]
    chain=cht[1] ## suppose une seule chaine
    #ind=grep(substring(target_file,1,4),patchlist)
    #true_patch_files=patchlist[ind]
    #ind=grep(paste(chain,"-patch",sep=""), true_patch_files)
    #true_patch_files=true_patch_files[ind]
    
    sel=NULL
    #super_len=cl_len
    if (length(numt)!=length(numq)) {
       message("inequal length in alignements!")
    }
    super_len=length(numt)
    for (i in 1:super_len) {
    	x=atom.select(target, resno=numt[i], insert=insertt[i], elety=namt[i], chain=cht[i])
	sel=c(sel, x$xyz)
    }
    target.ind=sel
    sel=NULL
    for (i in 1:super_len) {
    	x=atom.select(query, resno=numq[i], insert=insertq[i], elety=namq[i], chain=chq[i])
	sel=c(sel, x$xyz)
    }
    query.ind=sel

    ligand=trim.pdb(query, chain=pchain)
    cat(super_len, "target.ind=", length(target.ind), "mobile.inds", length(query.ind), "\n")
    result=fit.xyz(fixed=target$xyz, mobile=query$xyz, fixed.inds=target.ind, mobile.inds=query.ind, verbose=FALSE)
    query.moved=query
    query.moved$xyz=result
    hitrmsd=bio3d::rmsd(query.moved$xyz, b=target$xyz, a.inds=query.ind, b.inds=target.ind)
    cat("-----------> bs rmsd=", hitrmsd,"\n")
    ligand.moved=trim.pdb(query.moved, chain=pchain)
    ind=which(ligand.moved$atom$elesy%in%c("a","b","A","C"))
    ligand.moved$atom$elesy[ind]="C"
    pos=min(gregexpr("\\.",target_file)[[1]])
    movedfile=paste(prefix,"-",num_target,".",substring(target_file,pos+1),sep="")
    unlink(movedfile)
    target=suppressWarnings(bio3d::read.pdb(target_file))
    target=trim.pdb(target, string="protein", chain=cht)
    result=bio3d::binding.site(target, ligand.moved, cutoff=1.4, byres=FALSE)
    nclash=length(result$inds$atom)
    cat("nclash=", nclash, "nligand=", nrow(ligand.moved$atom), "\n")
    #if (abs(nrow(ligand.moved$atom)-nclash)>=10) {
    if (nclash<=3) {
       cat("creating ",movedfile,"\n")
       ligand.moved$atom$alt=NA # patch 
       write.pdb(pdb=ligand.moved, file = movedfile)
       if (!is.null(apo)) write.pdb(pdb=apo, file = movedfile, append=TRUE)
    }
}
#
