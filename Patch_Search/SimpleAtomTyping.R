#!/usr/local/bin/Rscript
library(bio3d)
library(argparser, quietly=TRUE)

p=arg_parser("Atom typing")
p=add_argument(p, "pdb", help="pdb file")
p=add_argument(p, "--outfile", default="", help="output pdb file")
argv = parse_args(p)

#
# pdb typing
#

patchsearch_type=function(pdbfile) {
    pdb=suppressWarnings(bio3d::read.pdb(pdbfile))
    if (any(is.na(pdb$atom$elesy))) {
       pdb$atom$elesy=substring(pdb$atom$elety,1,1) # si pas de elesy
    }
    pdb$atom$elesy[pdb$calpha]="a"
    pdb$atom$elesy[pdb$atom$elety=="CB"]="b"
    pdb$atom$insert[is.na(pdb$atom$insert)]=""

    resno=unique(pdb$atom[,c("resno", "insert", "chain")])
    for (k in 1:nrow(resno))  {
      resn=resno[k,1]
      ins=resno[k,2]
      ch=resno[k,3]
      resid=pdb$atom$resid[pdb$atom$resno==resn & pdb$atom$insert==ins & pdb$atom$chain==ch][1]# insert!
      ch=pdb$atom$chain[pdb$atom$resno==resn & pdb$atom$insert==ins & pdb$atom$chain==ch][1]
      if (resid=="HIS") {
        inds=bio3d::atom.select(pdb, resno=resn, insert=ins, chain=ch, elety=c("CG","CD2","CE1"))
	if (length(inds$atom)>0) {
	   pdb$atom$elesy[inds$atom]="A"
	}
      }
      if (resid=="PHE") {
        inds=bio3d::atom.select(pdb, resno=resn, insert=ins, chain=ch, elety=c("CG","CD1","CD2","CE1","CE2","CZ"))
	if (length(inds$atom)>0) {
	   pdb$atom$elesy[inds$atom]="A"
	}
      }
      if (resid=="TYR") {
        inds=bio3d::atom.select(pdb, resno=resn, insert=ins, chain=ch, elety=c("CG","CD1","CD2","CE1","CE2","CZ"))
        cat("TYR", ch, resn, ins,"\n")
        cat(inds$atom,"\n")
	if (length(inds$atom)>0) {
	   pdb$atom$elesy[inds$atom]="A"
	}
      }
      if (resid=="TRP") {
       	   # replacement of first cycle
	 inds=bio3d::atom.select(pdb, resno=resn, insert=ins, chain=ch, elety=c("CG","CD1","CD2","CE2"))
	 if (length(inds$atom)>0) {
	    pdb$atom$elesy[inds$atom]="A"
	 }
       	   # replacement of second cycle
       	 inds=bio3d::atom.select(pdb, resno=resn, insert=ins, chain=ch, elety=c("CD2","CE2","CE3","CZ2","CH2","CZ3"))
	 if (length(inds$atom)>0) {
	     pdb$atom$elesy[inds$atom]="A"
	 }
       }
   }
   pdb
}


pdbfile=argv$pdb
outfile=argv$outfile
pdb=patchsearch_type(pdbfile)

if (outfile=="") outfile=pdbfile
print(outfile)
write.pdb(pdb, file=outfile)


