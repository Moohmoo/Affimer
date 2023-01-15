#!/usr/local/bin/Rscript
# Rscript getBindingSite.R ~/3D/PDB_2017/g3/pdb1g3j.ent.gz  1g3jA:B+pep.xyz --qchain A --lchain B --add calpha --cutoff 10 --ligand --surface --output xyz

ACC=10

library(bio3d)
library(argparser)
library(vanddraabe)

patchsearch_HOME=Sys.getenv("patchsearch_HOME")
source(paste(patchsearch_HOME,"/NS3lib.R",sep=""))


p=arg_parser("A basic program to extract a protein surface patch bound to a ligand")
p=add_argument(p, "pdbfile", help="input pdb file")
p=add_argument(p, "outfile", default=NA, help="output pdb file")
p=add_argument(p, "--lchain",  default=NA, help="ligand chain selection")
p=add_argument(p, "--qchain", default=NA,help="query chain selection")
p=add_argument(p, "--cutoff", default="5",help="distance cutoff")
p=add_argument(p, "--add", default="", help="add complementary atoms:\"\",\"calpha\", \"cbeta\",\"cab\",\"backbone\", \"residue\", \"deep\"")
p=add_argument(p, "--ligand", flag=TRUE, help="add ligand to outfile")
p=add_argument(p, "--surface", flag=TRUE, help="keep surface atoms only")
p=add_argument(p, "--output", default="pdb", help="output format: pdb | xyz")
p=add_argument(p, "--radius", default=1.4, help="output pdb surface")

argv = parse_args(p)
pdbfile=argv$pdbfile
outfile=argv$outfile
add.ligand=argv$ligand
surface=argv$surface
output=argv$output
radius=argv$radius
#
# just to see the chains in pdb
#
pdb=patchsearch_type(pdbfile) 
pdb=bio3d::trim.pdb(pdb, string="protein")
chainlist=unique(pdb$atom$chain)

#
# process arguments
#
chain=argv$qchain
if (is.na(chain)) chain=chainlist else {
   chain=unlist(strsplit(chain,split=","))
   chain=intersect(chain, chainlist)
}

ligand=argv$lchain
if (is.na(ligand)) {
    message("no ligand chain")
    q()
}
ligand=unlist(strsplit(ligand,split=","))
ligchain=ligand
# ligand is a residu or a list of residu (not fully tested)
ligres=NULL
if (all(nchar(ligand)>1)) {
   ligchain=substring(ligand,1,1)
   ligres=as.integer(substring(ligand,first=2))
}
if (ligchain%in%chain) chain=setdiff(chain, ligchain)
#
#
#
cutoff=as.double(argv$cutoff)
add=argv$add
add_values=c("", "calpha","cbeta","cab","backbone", "residue","deep")
add=intersect(add, add_values)
if (length(add)==0) {
   message("add what?")
   q()
}

fullpdb=pdb
ligand.pdb=bio3d::trim.pdb(pdb, chain=ligchain, string="protein")
pdb=bio3d::trim.pdb(pdb, chain=chain, string="protein",resno=ligres)
pdb=surf_extract(pdb, chain, add=add, probeRadius=radius, minacc=ACC)
selbs=bio3d::binding.site(pdb, ligand.pdb, cutoff=cutoff, byres=FALSE)

if (length(selbs$inds$atom)==0) {
   message("no binding sites")
   q()
}

keep=selbs$inds$atom
if (surface) {
   message("extract surface")
   keep=intersect(selbs$inds$atom,which(pdb$atom$b>=ACC))
   pdb$atom$b[]=0
   pdb$atom$b[keep]=99
}

pdb=add.atoms(pdb, selatom=pdb$atom$eleno[keep], chain, add) # annote with 99 added atoms
keep=which(pdb$atom$b>=ACC)

# equivalent to options --extract of PS2
#if (surface) {
   eleno=pdb$atom$eleno[keep]
   chain=pdb$atom$chain[keep]
   ind=bio3d::atom.select(pdb,  eleno=eleno,  chain=chain)
   pdb=bio3d::trim.pdb(pdb, inds=ind)
#}



# add the ligand to the patch file for script SiteAlign.R
if (add.ligand) {
   ligand.pdb$atom$b[]=0
   pdb=suppressWarnings(bio3d::cat.pdb(pdb, ligand.pdb, rechain=FALSE, renumber=FALSE))
}
if (output=="pdb")  write.pdb(pdb, file=outfile)
if (output=="xyz") {
   print(pdb)
   D=pdb2xyz (pdb, access=surface)
   write.table(D,quote=FALSE, row=FALSE, file=outfile)
}
