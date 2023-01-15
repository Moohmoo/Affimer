library(argparser, quietly=TRUE)
p=arg_parser("A program to proceed results from an alignment file")
p=add_argument(p, "score", help="score file")
p=add_argument(p, "alfile", help="align file")
p=add_argument(p, "column", help="score column number")

argv = parse_args(p)
scorefile=argv$score
alignfile=argv$alfile
nscore=as.integer(argv$column)+1
outfile=paste(basename(alignfile),".index",sep="")
unlink(outfile)
if (!file.exists(alignfile)) {
   message("no align file")
   q()
}
#
# Pvalues estimation
#
D = read.table(scorefile, header=TRUE)
score = D$score
m = mean(score)
s = sd(score)
scale = sqrt(6)/pi*s
loc = m + scale*digamma(1) # digamma(1)=-Euler
sc = (score - loc)/scale
#faux pval = exp(-exp(sc))
pval = 1 - exp(-exp(-sc))
#
#
# sort alfile and add Pvalues 
#
L=readLines(alignfile)
Lindex=grep(">",L)
n=length(L)
Lindex=append(Lindex,n+1)
n=length(Lindex)
noal=which(Lindex[2:(n+1)]-Lindex[1:n]==1)
if (length(noal)>0) Lindex=Lindex[-noal]
#if (maxal>0) nbal=min(maxal,MAXALIGN)
#nbal=min(nbal, length(Lindex))
nbal=length(Lindex)
Lindex=Lindex[-nbal]
score=c()
num_target=0
for (k in Lindex) {
    num_target=num_target+1
    x=unlist(strsplit(L[k]," "))
    if (nscore>0) score[num_target]=as.double(x[nscore])
    y=unlist(strsplit(L[k+1]," "))
    if (nscore==0) score[num_target]=length(grep("CA",y))
    L[k]=paste(L[k]," ", pval[num_target], sep="")
    cat(num_target, basename(x[2]), score[num_target], pval[num_target], "\n",file=outfile, append=TRUE)
}

o=order(score,decreasing=TRUE)
Lindex=Lindex[o]
unlink(alignfile)
for (k in Lindex) {
    cat(L[k],"\n", file=alignfile, append=TRUE) 
    cat(L[k+1],"\n", file=alignfile, append=TRUE)
    cat(L[k+2],"\n", file=alignfile, append=TRUE)
}
