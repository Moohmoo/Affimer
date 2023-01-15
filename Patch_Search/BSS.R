# 
# search against binding site bank
#
patchsearch_MAXHITS=1
patchsearch_PRECISION=1.
patchsearch_NEIDIST=15
patchsearch_HOME=Sys.getenv("patchsearch_HOME")
library(argparser, quietly=TRUE)
#
#
#
p=arg_parser("A program to search a protein surface patch against a proteins surface")
p=add_argument(p, "query", help="pdb file")
p=add_argument(p, "bank",  help="pdb directory or a single file")
p=add_argument(p, "--tchain", default="",  help="target chain")
p=add_argument(p, "--bchain", default="",  help="bs chain")
p=add_argument(p, "--verbose",flag=TRUE, help="details about computations")
p=add_argument(p, "--outfile", default="", help="output file")
p=add_argument(p, "--alfile", default="", help="alignment output file")
p=add_argument(p, "--precision", default=patchsearch_PRECISION, help="")
p=add_argument(p, "--maxhits", default=patchsearch_MAXHITS, help="maximum number of hits")
p=add_argument(p, "--prefix", default="",  help="file prefix")
p=add_argument(p, "--minscore", default="0",  help="score threshold")
p=add_argument(p, "--align", flag=TRUE,  help="align or not")
p=add_argument(p, "--resi", default="",  help="select target residue")
p=add_argument(p, "--col", default="10",  help="column to sort")

argv = parse_args(p)
query=argv$query
bank=argv$bank
prefix=argv$prefix
outfile=argv$outfile
if (outfile=="") outfile=paste(prefix,".out",sep="")
alignfile=argv$alfile
align=argv$align
if (alignfile=="" & align) alignfile=paste(prefix,".al",sep="")
scorefile=paste(prefix,".out",sep="")
verbose=argv$verbose
precision=as.double(argv$prec)
minscore=as.double(argv$minscore)
resi=argv$resi
column=as.integer(argv$col)

if (nchar(resi)>0) resi=paste(" --resi ", resi, sep="")

if (precision<=0) precision=patchsearch_NEIDIST
maxhits=as.integer(argv$maxhits)

tchain=""
if (nchar(argv$tchain)>0) {
   tchain=paste(" --tchain ",argv$tchain,sep="")
}

bchain=""
if (nchar(argv$bchain)>0) {
   bchain=paste(" --qchain ",argv$bchain,sep="")
}

apofile=paste(prefix,"_apo.pdb",sep="")

if (!file.exists(apofile)) {
   typedfile=paste(prefix,"_t.pdb",sep="")
   command=paste("Rscript ", patchsearch_HOME,"/SimpleAtomTyping.R ",query, " -o ", typedfile, sep="")
   print(command)
   system(command)
   command=paste("Rscript ", patchsearch_HOME,"/AnnotateSurface.R ",typedfile," -o ",apofile, tchain, " --add calpha ",sep="")
   print(command)
   system(command)
   unlink(typedfile)
}
#
#
if (dir.exists(bank)) { # if bank is directory of bs files
   bslist=dir(bank)
} else if (file.exists(bank)) { # if bank is a bs file
   pos=gregexpr("/",bank)
   last=length(pos[[1]])
   bslist=substring(bank,first=pos[[1]][last]+1)
   bank=substring(bank,first=1,last=pos[[1]][last])
   if (bank=="") bank="./"
}
n=length(bslist)
if (!file.exists(outfile)) cat("bs query precision bslen alen clen rmsd coverage meandist score\n", file=outfile, append=TRUE)
for (bs in bslist) {
   cat("bs=",bs,"\n")
   id=toupper(substring(bs,1,4))
   bsfile=paste(bank,"/",bs,sep="")
   cat("bsfile=",bsfile,"\n")
   if (alignfile=="") {
      command=paste("Rscript ",patchsearch_HOME,"/PS3.R ", bsfile, bchain," ", apofile, tchain, " --surface", " -p ", precision, " --maxhits ",maxhits," ", resi, " --flex  -o ", outfile, sep="")
   } else {
      command=paste("Rscript ",patchsearch_HOME,"/PS3.R ", bsfile, bchain," ", apofile, tchain, " --surface", " -p ", precision, " --maxhits ",maxhits," ", resi, " --flex --alfile ", alignfile," -o ", outfile, sep="")
   }
   print(command)
   system(command)
}

if (align) {
   command=paste("Rscript ", patchsearch_HOME,"/sortResults.R ", outfile, " ", alignfile, " ", column,sep="")
   print(command)
   system(command)
   system(command)

   command=paste("Rscript ", patchsearch_HOME,"/SiteAlign.R ", alignfile," --outprefix ", prefix, sep="")
   print(command)
   system(command)
}

# unlink(apofile)
