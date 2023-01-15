patchsearch_HOME=Sys.getenv("patchsearch_HOME")
source(paste(patchsearch_HOME,"/NS3lib.R",sep=""))
# default values
patchsearch_PRECISION=1.
patchsearch_MAXHITS=1
patchsearch_ALFILE=""
patchsearch_MAXDELTADIST=1.0
patchsearch_DELTADIST=1
patchsearch_NEIDIST=15
patchsearch_SCORE=mapping_dist_sum2
patchsearch_BCMIN=0.0
#
#
#
library(argparser, quietly=TRUE)
# by default, cliques based on Calpha
#types=c("a","b")
#types=c("a","A")
#types="a"
types=c("a","A","N","O")
#
source(paste(patchsearch_HOME,"/NS3lib.R",sep=""))

output_results=function(id1,id2, prec, M, alen, clen, rms, fcov, fdist,fscore, outfile=outfile) {
   if (outfile=="" | !file.exists(outfile)) {
      cat("query target precision query_len match_len  clique_len rmsd coverage meandistorsion score\n",file=outfile,append=TRUE)
   } 
   cat(id1,id2, prec, M, alen, clen, rms, fcov, fdist, fscore, "\n",file=outfile,append=TRUE)
}

#
#
#
p=arg_parser("A program to search a protein surface patch against a proteins surface")
p=add_argument(p, "query", help="pdb file")
p=add_argument(p, "target",  help="pdb file")
p=add_argument(p, "--qchain", default="",  help="query chain")
p=add_argument(p, "--tchain", default="",  help="target chain")
p=add_argument(p, "--verbose",flag=TRUE, help="details about computations")
p=add_argument(p, "--outfile", default="", help="output file")
p=add_argument(p, "--alfile", default=patchsearch_ALFILE, help="alignment output file")
p=add_argument(p, "--maxhits", default=patchsearch_MAXHITS, help="maximum number of hits")
p=add_argument(p, "--prec", default=patchsearch_PRECISION, help="")
p=add_argument(p, "--clprec", default=patchsearch_PRECISION, help="")
p=add_argument(p, "--flex", flag=TRUE, help="always this method")
p=add_argument(p, "--surface", flag=TRUE, help="consider surface atoms only")
p=add_argument(p, "--resi",  default="", help="select target residues")

argv = parse_args(p)
queryfile=argv$query
targetfile=argv$target
qchain=argv$qchain
qchain=unlist(strsplit(qchain,split=","))
if (length(qchain)==0) qchain=""
tchain=argv$tchain
tchain=unlist(strsplit(tchain,split=","))
if (length(tchain)==0) tchain=""
alignfile=argv$alfile
outfile=argv$outfile
verbose=argv$verbose
precision=as.double(argv$prec)
if (precision<=0) precision=patchsearch_NEIDIST
cl_prec=as.double(argv$clprec)
flex=argv$flex
maxhits=as.integer(argv$maxhits)
surface=argv$surface
resi=argv$resi
resi=unlist(strsplit(resi,split=","))
resi=as.integer(resi)

#unlink(alignfile)

score_thresh=0
deep=FALSE
#
if (verbose) cat("begin PS2...\n")
#
# read the query pdb file
#
if (!file.exists(queryfile)) {
       mess=paste(queryfile, "no file")
       message(mess)
       q()
}

cat("surface=",surface,"\n")
query=read_pdb(queryfile, chain=qchain, access=surface)#FALSE=all atoms
#if (!surface) query$atom$b=99
cat("query chains=", unique(query$ch),"\n")
#
# test
#
Y=query$coord
M=query$natoms
# patch_maxsize en test
patch_maxsize=patchsearch_NEIDIST
#
# read the target pdb file
#
print(targetfile)
if (!file.exists(targetfile)) {
       mess=paste(targetfile, "is not available")
       message(mess)
       q()
}

cat("surface=",surface,"\n")
target=read_pdb(targetfile, chain=tchain, access=surface, resi=resi)
#if (!surface) target$atom$b=99

cat("target chains=", unique(target$chains),"\n")
if (1) {
    cat("query composition:")
    print(sort(table(query$type[query$access>0])))
    cat("target composition:")
    print(sort(table(target$type[target$access>0])))
}

if (sum(query$type%in%types)<patchsearch_MINCLIQUEATOMS | sum(target$type%in%types)<patchsearch_MINCLIQUEATOMS) {
   types=c("a","A","N","O")
}
cat("types=", types, "\n")

N=target$natoms
if (N==0) {
       message ("empty target: no atoms")
       q()
}

if (M>N) {
        message ("query larger than target")
        #output_nohit("query larger than target ", queryfile, targetfile)
        #q()
	M=min(M,N)
	N=max(M,N)
}
#
# search the query on all the chains of target
#

   qcliques=patch_search(query, target, deltadist=patchsearch_DELTADIST, maxdist=patchsearch_MAXDELTADIST, neidist=patch_maxsize, precision=precision, thresh=cl_prec, merge=TRUE, ptypes=types, verbose=verbose, flex=flex, maxcliques=maxhits, score_function=patchsearch_SCORE)

if (qcliques$ok==TOO_BIG) { # can not occur now !!
   mess="no result,graph too big or empty"
   message(mess)
   q() #no output
}

if (qcliques$ok<0) {
    mess="no valid clique"
    message(mess)
    imax=1
    output_results(queryfile, targetfile, qcliques$ok, M, qcliques$alen[imax],  qcliques$clen[imax], qcliques$rms[imax], qcliques$cov[imax], qcliques$dist[imax], qcliques$score[imax], outfile=outfile)
    q()
}
#
# output the nb_output best matches
#
    #imax=which.max(qcliques$score)
    #nb_output=1
    nb_output=length(qcliques$target)
    if (maxhits>0) nb_output=min(nb_output, maxhits)
    for (imax in 1:nb_output) {
	I1=qcliques$target[[imax]] # target indices
	I2=qcliques$query[[imax]]  # query indices
  	if (qcliques$bc[imax]>patchsearch_BCMIN & length(I1)>6) {
	  output_results(queryfile, targetfile, precision, M, qcliques$alen[imax],qcliques$clen[imax], qcliques$rms[imax], qcliques$cov[imax], qcliques$dist[imax], qcliques$score[imax],outfile=outfile)
	  if (alignfile!="") {
	     coma = apply(query$coord[I2,], 2, mean)
	     com =  apply(query$coord, 2, mean)
	     dcom = sqrt(sum((coma-com)**2))
	     if (dcom >= 10) next
	     ind=1:length(I1)
	     #ind=which(target$access[I1]>0 & query$access[I2]>0) # select atom on surface
	      I1=I1[ind]
	      I2=I2[ind]
	     query_res=paste(query$resnum[I2],":",query$aname[I2],":",query$chains[I2], sep="")
	     target_res=paste(target$resnum[I1],":",target$aname[I1],":",target$chains[I1], sep="")
	     cat(">>",queryfile, targetfile, precision, qcliques$len[1], qcliques$alen[imax],  qcliques$clen[imax], qcliques$rms[imax], qcliques$cov[imax],  qcliques$dist[imax], qcliques$score[imax], "\n", file=alignfile, append=TRUE)
	     cat(query_res, "\n", file=alignfile, append=TRUE)
	     cat(target_res, "\n", file=alignfile, append=TRUE)
	  }
      }
    }


