library(argparser, quietly=TRUE)
library(bio3d)
library(vanddraabe)
ACC=0

SASA=function(pdb, probeRadius = radius) {
    result=FreeSASA.diff(atoms.oi = pdb$atom, probeRadius = radius)
    acc=result$SASA.prot
    x=strsplit(result$uniq.atom.ids,split="_")
    for (i in 1:length(acc)) {
        resno=as.integer(x[[i]][2])
        chain=x[[i]][3]
        eleno=as.integer(x[[i]][5])
	ind=atom.select(pdb, resno=resno, chain=chain, eleno=eleno)
	#ind=atom.select(pdb,  eleno=eleno)
	pdb$atom[ind$atom,"b"]=acc[i]
    }
    pdb
}

surf_annotate=function(pdb, chains, add="calpha", minacc=ACC, probeRadius = radius) {
    cat("surf_annotate...\n")
    pdb$atom$b[]=0 # otherwise all chain not in chains remains acc.
    pdb=chain_surf_annotate(pdb, chains[1], minacc, add, probeRadius)
    
    if (length(chains)>1) {
	for (i in 2:length(chains)) {
	    pdb=chain_surf_annotate(pdb, chains[i], minacc, add, probeRadius)
	    # surf.pdb=bio3d::cat.pdb(surf.pdb, tmp.pdb, renumber=FALSE, rechain=FALSE)
	    # renumber=FALSE doesn't work when chain break
	 }
    }
    return(pdb)
}

chain_surf_annotate=function(fullpdb, chain, minacc, add, probeRadius) {
    cat("chain_surf_extract=", chain, "\n")
    selchain=bio3d::atom.select(fullpdb, string="protein", chain=chain)
    pdb=bio3d::trim.pdb(fullpdb,string="protein",chain=chain)
    modified.pdb=SASA(pdb, probeRadius)# b-factors replaced by acc
    modified.pdb$atom$b[modified.pdb$atom$b<minacc]=0
    fullpdb$atom[selchain$atom,]=modified.pdb$atom
    selatom=fullpdb$atom$eleno[fullpdb$atom$b>=minacc] # seleted indices of surface atoms (eleno)
    ind=c()# 
    if (add=="deep") {
       surf.pdb=bio3d::trim(fullpdb, inds=ind)
       ind=atoms.add.layer(fullpdb, surf.pdb, chain)
    }
    if (add=="calpha") ind=add.atoms(fullpdb, selatom, chain, string="calpha")
    if (add=="cbeta") ind=add.atoms(fullpdb, selatom, chain, string="", elety="CB")
    if (add=="cab") ind=add.atoms(fullpdb, selatom, chain, string="", elety=c("CA","CB"))
    if (add=="backbone") ind=add.atoms(fullpdb, selatom, chain, string="backbone")
    if (add=="residue") ind=add.atoms(fullpdb, selatom, chain, string="protein")
    fullpdb$atom$b[ind]=99
    cat("len select=", sum(fullpdb$atom$chain==chain),sum(fullpdb$atom$b>=minacc & fullpdb$atom$chain==chain),"\n")
    return(fullpdb)
}


# eleno+chain=selected atoms (surface or patch atoms)
add.atoms=function(pdb, eleno, chain, string="calpha", elety=NULL) {
    sel=bio3d::atom.select(pdb, eleno=eleno, chain=chain)
    res=unique(pdb$atom[sel$atom,c("resno", "insert", "chain")])
    selresn=res[,1]
    selins=res[,2]
    selch=res[,3]
 
    #selres=unique(pdb$atom$resno[sel$atom])
    sel=bio3d::atom.select(pdb, string=string, elety=elety, chain=selch, resno=selresn, insert=selch)
    sel$atom
}

atoms.add.layer=function(pdb, surf.pdb, chain, string="calpha") {
    ind=bio3d::atom.select(pdb, chain=chain)
    result=bio3d::binding.site(pdb, a.inds=ind, b=surf.pdb,  cutoff=cutoff, byres=FALSE)
    #bio3d::trim.pdb(pdb,result$inds, string=string)
    result$inds
}

p = arg_parser("A program to extract a protein surface")
p = add_argument(p, "pdbfile", help="input pdb file")
p = add_argument(p, "--outfile", default="", help="output pdb file")
p = add_argument(p, "--tchain", default=NA,help="protein chain")
p = add_argument(p, "--lchain", default="",help="ligand chain")
p = add_argument(p, "--add", default="", help="add complementary atoms:\"\",\"calpha\", \"backbone\", \"residue\", \"deep\"")
p=add_argument(p, "--cutoff", default="5",help="distance cutoff for buried atoms")
p=add_argument(p, "--surface", flag=TRUE, help="output pdb surface")
p=add_argument(p, "--radius", default=1.4, help="output pdb surface")
argv = parse_args(p)
pdbfile=argv$pdbfile
outfile=argv$outfile
if (outfile=="") outfile=pdbfile
chain=argv$tchain
ligand=argv$lchain
surface=argv$surface
radius=argv$radius

add=argv$add
add_values=c("", "calpha","cbeta", "cab", "backbone", "residue","deep")
add=intersect(add, add_values)
if (length(add)==0) {
   message("add what?")
   q()
}
cutoff=as.double(argv$cutoff)

# main
pdb=suppressWarnings(bio3d::read.pdb(pdbfile))
pdb=bio3d::trim.pdb(pdb, string="protein")
if (nrow(pdb$atom)==0) {
	message(paste(pdbfile," is not a protein"))
	q()
}

pdb$atom$insert[is.na(pdb$atom$insert)]=""

# chain=user argument; chains=surf_extract argument; ch=pdb chains
ch=unique(pdb$atom$chain)
if (is.na(chain)) {
   chains=ch
} else {
   chain=unlist(strsplit(chain,split=","))
   if (!all(chain%in%ch)) {
      message("chain not available")
      q()
   }
   chains=chain
}
if (ligand%in%chains) chains=setdiff(chains,ligand)
print(chains)

surf.pdb=surf_annotate(pdb, chains, add=add, probeRadius=radius, minacc=ACC)
ligsel=bio3d::atom.select(surf.pdb, chain=ligand)
surf.pdb$atom$b[ligsel$atom]=99 # force to keep peptide complete
if (surface) {
   eleno=surf.pdb$atom$eleno[surf.pdb$atom$b>=ACC]
   chain=surf.pdb$atom$chain[surf.pdb$atom$b>=ACC]
   ind=bio3d::atom.select(surf.pdb,  eleno=eleno,  chain=chain)
   surf.pdb=bio3d::trim.pdb(surf.pdb, inds=ind)
}

bio3d::write.pdb(surf.pdb, file=outfile)

