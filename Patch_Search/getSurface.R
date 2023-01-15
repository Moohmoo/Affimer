library(argparser, quietly=TRUE)
library(bio3d)
library(vanddraabe)
ACC=10

patchsearch_HOME=Sys.getenv("patchsearch_HOME")
source(paste(patchsearch_HOME,"/NS3lib.R",sep=""))


p = arg_parser("A program to extract a protein surface")
p = add_argument(p, "pdbfile", help="input pdb file")
p = add_argument(p, "outfile", help="output pdb file")
p = add_argument(p, "--chain", default=NA,help="protein chain")
p = add_argument(p, "--add", default="", help="add complementary atoms:\"\",\"calpha\",  \"cbeta\", \"cab\", \"backbone\", \"residue\", \"deep\"")
p = add_argument(p, "--cutoff", default="5",help="distance cutoff for buried atoms")
p = add_argument(p, "--radius", default=1.4, help="output pdb surface")
p = add_argument(p, "--output", default="pdb", help="output format: pdb | xyz")

#
# process arguments
#
argv = parse_args(p)
pdbfile=argv$pdbfile
outfile=argv$outfile
chain=argv$chain
add.ligand=argv$ligand
radius=argv$radius
output=argv$output

add=argv$add
add_values=c("", "calpha","backbone", "cbeta", "cab", "residue","deep")
add=intersect(add, add_values)
if (length(add)==0) {
   message("add what?")
   q()
}
cutoff=as.double(argv$cutoff)

#
pdb=patchsearch_type(pdbfile) 
pdb=bio3d::trim.pdb(pdb, string="protein")
if (nrow(pdb$atom)==0) {
	message(paste(pdbfile," is not a protein"))
	q()
}

# chain=user argument; chains=surf_extract argument; ch=pdb chains
ch=unique(pdb$atom$chain)
if (is.na(chain)) {
   chains=ch
} else {
   chain=unlist(strsplit(chain,split=","))
   chain=intersect(chain,ch)
   if (length(chain)==0) {
      message("chains not available")
      q()
   }
   chains=chain
}

pdb=bio3d::trim(pdb, chain=chains)
surf.pdb=surf_extract(pdb, chains, add=add, probeRadius=radius,minacc=ACC)
eleno=surf.pdb$atom$eleno[surf.pdb$atom$b>=ACC]
chain=surf.pdb$atom$chain[surf.pdb$atom$b>=ACC]
ind=bio3d::atom.select(surf.pdb,  eleno=eleno,  chain=chain)
surf.pdb=bio3d::trim.pdb(surf.pdb, inds=ind)
if (output=="pdb") {
   bio3d::write.pdb(surf.pdb, file=outfile)
}
if (output=="xyz") {
   D=pdb2xyz(surf.pdb, access=TRUE)
   write.table(D,quote=FALSE, row=FALSE, file=outfile)
}

