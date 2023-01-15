library(Rcpp)
patchsearch_HOME=Sys.getenv("patchsearch_HOME")
cppsource=paste(patchsearch_HOME, "/NS3lib.cpp",sep="")
sourceCpp(cppsource, cacheDir=paste(patchsearch_HOME,"/Tmp",sep=""))

library(bio3d)
patchsearch_MAXDELTADIST=3
patchsearch_MINDIST=0
patchsearch_MINCLIQUE=6
patchsearch_MINCLIQUEATOMS=20
patchsearch_BCMIN=0
patchsearch_INTERCLIQUE=4
patchsearch_NBNEI=4
patchsearch_MAXEDGES=3000000
EMPTY_GRAPH=-3
NO_CLIQUE=-2
TOO_BIG=-1

#patchsearch_AA=c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","HSD","HSE","HSP")
#patchsearch_codes=0:22
#names(patchsearch_codes)=patchsearch_AA

#' Compute BC score between X and Y
#'
#' @param X matrix Nx3
#' @param Y matrix Nx3
#' @export
BCscore<-function(X,Y) {
    X=scale(X,scale=FALSE)
    Y=scale(Y,scale=FALSE)
    if (det(t(X)%*%X)<1 | det(t(Y)%*%Y)<1) return(0)
    det(t(X)%*%Y)/(sqrt(det(t(X)%*%X))*sqrt(det(t(Y)%*%Y)))
}


#' Compute rmsd between X and Y
#'
#' @param X matrix Nx3
#' @param Y matrix Nx3
#' @export
patchsearch_rmsd=function(X,Y) {
   r=bio3d::rmsd(as.vector(t(X)),as.vector(t(Y)),fit=TRUE)
   r
}

select.chain=function(resno) {
    resno=unique(resno)
    n=length(resno)
    d=diff(resno)
    keep=NULL
    if (d[1]<=1) keep=c(keep,resno[1])
    for (k in 2:(n-1)) {
    	if (d[k-1]<=1 | d[k]<=1)
	   keep=c(keep,resno[k])
    }
    if (d[n-1]<=1) keep=c(keep,resno[n])
    return(keep)
}



#' read a PDB file 
#'
#' @param pdbfile  pdb file
#' @export 
read_pdb <- function(pdbfile, chain="", reduce=FALSE, access=FALSE, resi=integer(0)) {
	pdb=suppressWarnings(bio3d::read.pdb(pdbfile, rm.insert=FALSE))
	if (length(resi)>0) pdb=bio3d::trim.pdb(pdb, resno=resi)
	cat("trim to accessible atoms", access, nrow(pdb$atom),"-->")
	if (access) {
	   pdb=bio3d::trim.pdb(pdb, eleno=pdb$atom$eleno[pdb$atom$b>0])
	}
	cat(nrow(pdb$atom),"\n")
#	pdb=bio3d::trim.pdb(pdb, string="protein", resid=patchsearch_AA)
	pdb=bio3d::trim.pdb(pdb, string="protein")
	if (chain[1]!="") {
#	   pdb=bio3d::trim.pdb(pdb, string="protein", resid=patchsearch_AA, chain=chain)
	   pdb=bio3d::trim.pdb(pdb, string="protein", chain=chain)
	}
	cat(nrow(pdb$atom),"\n")
	if (reduce) {
   	   resno=select.chain(pdb$atom$resno) #??
   	   pdb=bio3d::trim.pdb(pdb, resno=resno)
	}
	cat("trim to accessible atoms",nrow(pdb$atom),"\n")
	x=pdb$atom$x
	y=pdb$atom$y
	z=pdb$atom$z
	access=pdb$atom$b
	type=gsub(" ","",pdb$atom$elesy)
	n=length(type)
	#resname=patchsearch_codes[pdb$atom$resid]
	resname=aa321(pdb$atom$resid)
	ch=pdb$atom$chain
	ch[is.na(ch)]=" "
	pdb$atom$chain=ch
	# modif du 4/7/18 (sinon non bijection atomnum et resnum:chain:atomname)
	pdb$atom$insert[is.na(pdb$atom$insert)]=""
	resnum=paste(pdb$atom$resno,":",pdb$atom$insert,sep="")
	#resnum=pdb$atom$resno
	#
    	list(coord=cbind(x,y,z), type=type, access=access, num=pdb$atom$eleno, resno=pdb$atom$resno, resnum=resnum, insert=pdb$atom$insert, resname=resname, natoms=n, aname=pdb$atom$elety, chains=ch)
}
    
## version du 19/12/18
extract_cliques<-function(gp) {
        #Cl=igraph::largest_cliques(gp)
	size=igraph::clique.number(gp)
	cat("extract_cliques-->", igraph::clique.number(gp),"\n")
        Cl=igraph::maximal.cliques(gp,min=size-1)
	#Cl=cliques(gp, min=6,max=6)
	if (length(Cl)==0) return(list(clique=NULL,size=0))
	##
	for (i in 1:length(Cl)) Cl[[i]]=as.integer(names(Cl[[i]]))
	##
	list(clique=Cl, size=max(mapply(length,Cl)))
}	

union_clique <- function(C1,C2,N) {
  I1=(C1-1)%%N+1
  I2=(C2-1)%%N+1
  J1=(C1-1)%/%N+1
  J2=(C2-1)%/%N+1
  I=which(!I2%in%I1 & !J2%in%J1)
  base::union(C1,C2[I])
}
#
# sans BC
#
merge_cliques<-function(clusters,N,common=patchsearch_INTERCLIQUE) {
   cat("------> merge_cliques", length(clusters),"\n")
   stop=FALSE
   while (!stop) {
         stop=TRUE
   	 n=length(clusters)
	 i=1
  	 while (i<= n-1) {
	   j=i+1
	   while (j<= n) {
	     if (length(base::intersect(clusters[[i]],clusters[[j]]))>=common) {
	     	#clusters[[i]]=union(clusters[[i]],clusters[[j]]) # may create forks
	     	C=union_clique(clusters[[i]],clusters[[j]],N)
 		clusters[[i]]=C
		clusters[[j]]=c()
		n=n-1
		j=j-1
		stop=FALSE
	     }
	     j=j+1
     	    }
	    i=i+1
	  }
   }
   cat("------> merge_cliques", length(clusters),"\n")
   return(clusters)
}
#
# sort clusters
#
sort_clusters=function(clusters,N, X, XRes, Y, YRes, score_function=mapping_dist_sum2, thresh=0) {
    nbclique=length(clusters)
    if (verbose) cat("nb of clique clusters=",nbclique,"\n")
    bc=alen0=double(nbclique)
    for (i in 1:nbclique) {
    	C=clusters[[i]]
   	size=length(C)
   	I=(C-1)%%N+1 # target indices
   	J=(C-1)%/%N+1# query indices
   	s=sign(BCscore(X[I,], Y[J,]))
	bc[i]=s*matching_score(C,X,XRes,Y,YRes, score_function, thresh) 
	alen0[i]=length(C)
    }
    clusters=clusters[bc>=0]
    bc=bc[bc>=0]
    o=order(bc,decreasing=TRUE)
    clusters=clusters[o]
    clusters
}
#
#
#
remove_redundant_clusters = function(clusters, common=patchsearch_INTERCLIQUE) {
   n=length(clusters)
   i=1
   while (i<= n-1) {
   	 j=i+1
	 while (j<= n) {
	       lmin=min(length(clusters[[i]]),length(clusters[[j]]),common)
	       #lmin=common
	       if (length(base::intersect(clusters[[i]],clusters[[j]]))>=lmin) {
		  clusters[[j]]=c()
		  n=n-1
		  j=j-1
	       }
	       j=j+1
     	 }
	 i=i+1
   }
   return(clusters)
}

#
#
#
build_graph=function(types, deltadist, X, XProp, XRes, Xaccess, Y, YProp, YRes, Yaccess, mindist, maxdist, thresh=0.9, verbose=FALSE) {
	N=nrow(X)
	J=which(YProp%in%types)
	I=which(XProp%in%types)
	Jaccess=which(Yaccess>0)
	Iaccess=which(Xaccess>0)
	I=intersect(I, Iaccess)
	J=intersect(J, Jaccess)
	cat ("build_graph::",length(I),length(J),length(intersect(I, Iaccess)),length(intersect(J, Jaccess)),"\n")
	Mmotif=length(J)
	Nmotif=length(I)
	if (Mmotif<patchsearch_MINCLIQUEATOMS | Nmotif<patchsearch_MINCLIQUEATOMS) {
   	   types=c("a","A","C","O","N")
	   J=which(YProp%in%types)
	   I=which(XProp%in%types)
	   I=intersect(I, Iaccess)
	   J=intersect(J, Jaccess)
	   cat ("+build_graph::",length(I),length(J),length(intersect(I, Iaccess)),length(intersect(J, Jaccess)),"\n")
	}
	V=vertexc(XProp[I], YProp[J])
	if (verbose) cat("vertices:",length(I), length(J),length(V)/2,"\n")
	if (all(V==0) | length(V)<=1) return(NULL)
	V=cbind(J[V[,1]],I[V[,2]]) # atom ids 1..M, 1..N
	#nV=nrow(V)
	E=buildGraphc(X, Y, V, thresh, mindist, maxdist)
	while (nrow(E)>patchsearch_MAXEDGES){
	   cat("edges:",nrow(E),"\n")
	   maxdist=maxdist/2.;
	   E=buildGraphc(X, Y, V, thresh, mindist, maxdist)
	}
	cat("edges:",nrow(E),"\n")
	if (nrow(E)==0) {
	   message("no edge...")
	   return(igraph::make_empty_graph())
	}
	P=cbind((E[,3]-1)*N+E[,1],(E[,4]-1)*N+E[,2])
	mode(P)="character" #vertex label and not id to avoid huge graphs with non connected nodes
        gp=igraph::graph.edgelist(P, directed=FALSE)
	if (verbose) cat("graph:vertices:",igraph::vcount(gp), "graph.edges:", igraph::ecount(gp),"\n")
	return(gp)
}

clique_patchsearch=function(g, X, Y, minclique=patchsearch_MINCLIQUE, bcmin=patchsearch_BCMIN) {
	N=nrow(X);M=nrow(Y)
	cl=extract_cliques(g)
	##
	if (cl$size<minclique) {
	   cat("cl$size<minclique\n")
	   return(NULL)
	}
	nbclique=length(cl$clique)
	size=bc=double(nbclique)
	for (i in 1:nbclique) {
	    C=cl$clique[[i]]
	    size[i]=length(C)
	    I=(C-1)%%N+1 # target indices
	    J=(C-1)%/%N+1# query indices
	    bc[i]=BCscore(X[I,], Y[J,])
	}
	cl$clique=cl$clique[bc>=bcmin]
	#
	if (length(cl$clique)==0) return(NULL)
	bc=bc[bc>=bcmin]
	o=order(bc,decreasing=TRUE)
	cl$bc=bc[o]
	cl$clique=cl$clique[o]
	cl$size=size[o]
	return(cl)
}


filter_greedy=function(C, X, XRes, Y, YRes, score_function, thresh) {
    CI=(C-1)%%N+1 # target indices
    CJ=(C-1)%/%N+1# query indices
    W=score_function(X, XRes, CI, CI, Y, YRes, CJ, CJ, thresh)
    #maxs=sum(W)
    ind=which.min(W)
    while (W[ind]<=0) {
    	  CI=CI[-ind];CJ=CJ[-ind];C=C[-ind]
    	  W=score_function(X, XRes, CI, CI, Y, YRes, CJ, CJ,thresh)
   	  ind=which.min(W)
    }
    C
}

matching_score=function(C, X, XRes, Y, YRes, score_function, thresh=0) {
    CI=(C-1)%%N+1 # target indices
    CJ=(C-1)%/%N+1# query indices
    W=score_function(X, XRes, CI, CI, Y, YRes, CJ, CJ,thresh)
    sum(W)
}

max_bipartite_score=function(C, K, X, XRes, Y, YRes, gp, score_function, thresh) {
	C0=C
  	CI=(C-1)%%N+1 # target indices
   	CJ=(C-1)%/%N+1# query indices
       	KI=(K-1)%%N+1
	KJ=(K-1)%/%N+1
       	W=score_function(X, XRes,KI, CI, Y,YRes, KJ, CJ, thresh)
	W[W<=0]=-1
	#
       	igraph::E(gp)$weight=W
	igraph::V(gp)$type=igraph::bipartite_mapping(gp)$type
      	result=igraph::max_bipartite_match(gp)
	x=result$matching[igraph::V(gp)$type]
	x=x[!is.na(x)]
	jnd=as.integer(substring(x,first=2))
	ind=as.integer(substring(names(x),first=2))
	#W[W<=0]=0 # inutile
	list(C=(jnd-1)*N+ind, weight=W, score=result$matching_weight)
}


max_bipartite_flex=function(C, K, X, XRes, Y, YRes, thresh, verbose=FALSE, score_function=mapping_dist_sum2) {
	Nref=min(nrow(X), nrow(Y))
       	KI=(K-1)%%N+1
	KJ=(K-1)%/%N+1
      	P=cbind(paste("q",KJ,sep=""),paste("t",KI,sep=""))
       	gp=igraph::graph.edgelist(P, directed=FALSE)
	stop=FALSE
	niter=0
	Cprev=C
	C0=C
	    stop=FALSE
	    niter=0
	    score_init=matching_score(C,X, XRes, Y, YRes,score_function,thresh)
	    if (verbose) cat("C0=",length(C),score_init, score_init/length(C),"\n")
	    while(!stop & niter<20) {
		result=max_bipartite_score(Cprev, K, X, XRes, Y, YRes, gp, score_function, thresh) #f(C1=result$C,C0=Cprev)
		nC=length(result$C)
		C=filter_greedy(result$C,X,XRes,Y,YRes,score_function, thresh) # C1bar
		if (verbose) cat("C0=",length(Cprev),"C1=",length(result$C),"C1b=",length(C),"C0.C1=",length(intersect(result$C,Cprev)),"C0.C1b=",length(intersect(C,Cprev)),"f(C1,C0)=",result$score,"\n")
		if (all(Cprev%in%C) & all(C%in%Cprev)) {stop=TRUE}
		if (verbose) cat("f(C1)=", matching_score(result$C,X, XRes,Y,YRes,score_function,thresh),"f(C1b)=", matching_score(C,X,XRes,Y,YRes,score_function,thresh),"\n")
		Cprev=C
		niter=niter+1
	    }
	    C=union(C0,C) # hope that clique comes first
	    nC=length(C)
	    score=matching_score(C,X, XRes,Y,YRes, score_function,thresh)
	    dist=thresh*(1-1/nC**2*score)
	    if (stop & verbose) cat ("convergence\n")
	list(C=C, score=score/nC, dist=dist, coverage=nC/Nref)
}



# search for matchings between query pdb and a surface pdb files
patch_search=function(query, target, deltadist, maxdist, neidist, ptypes=patchsearch_types, merge=TRUE, greedy=TRUE, verbose=FALSE, precision=precision, thresh=cl_prec, flex=TRUE, maxcliques=1, score_function) {
    X=target$coord
    Y=query$coord
    XProp=target$type
    YProp=query$type
    cat("XProp=", unique(target$type), "\n")
    cat("YProp=", unique(query$type), "\n")
    
    XRes=target$resno #
    YRes=query$resno # 
    Xaccess=target$access
    Yaccess=query$access
    N=nrow(X)
    M=nrow(Y)
    cat("patch_search::M,N=",M,N,"\n")
    
    #
    # construction du graphe produit
    #
    cliques=NULL
    cat(ptypes, deltadist, nrow(X), nrow(Y), patchsearch_MINDIST, neidist, thresh,"\n")
    ## test!!
    ## clique_types="a"
    ##
    if (M==0 | N==0) {
	 if (verbose) cat("empty pdb...\n")
   	 return(list(target=NULL, query=NULL, len=M, alen=0, rmsd=100, cov=0, score=0, dist=0, ok=EMPTY_GRAPH))
    }
    graph=build_graph(ptypes, deltadist, X, XProp, XRes, Xaccess, Y, YProp, YRes, Yaccess, patchsearch_MINDIST, neidist, thresh, verbose)
    if (length(igraph::E(graph))==0) {# from make_empty_graph()
	 if (verbose) cat("empty graph...\n")
   	 return(list(target=NULL, query=NULL, len=M, alen=0, rmsd=100, cov=0, score=0, dist=0, ok=EMPTY_GRAPH))
    }
     if (verbose) cat("clique search...\n")
#
# recherche des cliques
#
    cliques=clique_patchsearch(graph, X, Y, patchsearch_MINCLIQUE, patchsearch_BCMIN)
    if (verbose) cat("number of cliques=",length(cliques$clique),"\n")
    if (is.null(cliques)) {
   	return(list(target=NULL, query=NULL, len=M, alen=0, rmsd=100, cov=0, score=0, dist=100, ok=NO_CLIQUE))
    }
    cliques=cliques$clique
#
# merging cliques
#
    if (verbose) cat("clique merge...\n")
    clusters=merge_cliques(cliques,N, patchsearch_INTERCLIQUE)
    nbclique=length(clusters)
    clusters=sort_clusters(cliques,N,X,XRes,Y,YRes,score_function) 
    nbclique=length(clusters)
    if (nbclique==0) {
	message("no valid clique after merging")
   	return(list(target=NULL, query=NULL, len=M, alen=0, rmsd=100, cov=0, score=0, dist=100, ok=NO_CLIQUE))
    }
    if (1) {
       cat("nb of selected clique clusters=",nbclique,"\n")
    }
#
# enriching the maxcliques best valid cliques 
#
nbclique=min(maxcliques,nbclique)

#thresh=precision
if (greedy==TRUE) {
##
    V=vertexc(XProp, YProp)
    cov=double(nbclique)
    dist=double(nbclique)
    score=double(nbclique)
    clen=double(nbclique)
    for (ic in 1:nbclique) {
       C=clusters[[ic]]
       nbefore= length(C)
       cat("enlarging clique ", ic, "\n")
    #
    # calcul de K=ensemble de links à tester pour enrichir clique choisie
    # 1. ensemble des voisins de C
    #
       #K=getConnectedNeighborsc(C, V, X, Y, thresh, patchsearch_NBNEI)
       #K=union(K,C)
       K=(V[,1]-1)*N+V[,2]
       K=union(K,C)
       clen[ic]=length(C)
       result=max_bipartite_flex(C, K, X, XRes, Y, YRes, thresh=precision, verbose, score_function)
       Cplus=result$C
       clusters[[ic]]=Cplus
       cov[ic]=result$coverage
       dist[ic]=result$dist
       score[ic]=result$score
    }
}
#
# remove redundant clusters
#
if (0) {
cat("------> remove_redundant_clusters", length(clusters),"\n")
o=order(unlist(lapply(clusters, length)), decreasing=TRUE)
clusters=clusters[o]
clusters=remove_redundant_clusters(clusters, 6)
cat("------> remove_redundant_clusters", length(clusters),"\n")
}
#
# returns matched atoms as indices of query (Iquery) and target (Itarget)
#
   nclique=length(clusters)
   Itarget=list()
   Iquery=list()
   rms=double(nbclique)
   alen=integer(nbclique)
   bc=double(nbclique)
#
#
#
   for (i in 1:nbclique) {
       I=(clusters[[i]]-1)%%N+1 # target indices
       J=(clusters[[i]]-1)%/%N+1# query indices
       alen[i]=length(J)
       rms[i]=patchsearch_rmsd(X[I,],Y[J,])
       bc[i]=BCscore(X[I,],Y[J,])
       #cat("clique ic=",rms[i],bc[i],"\n")
       Itarget[[i]]=I
       Iquery[[i]]=J
   }
   o=order(score, decreasing=TRUE)
   Itarget=Itarget[o]
   Iquery=Iquery[o]
   alen=alen[o]
   clen=clen[o]
   rms=rms[o]
   cov=cov[o]
   score=score[o]
   dist=dist[o]
   return(list(target=Itarget, query=Iquery, len=M, alen=alen, clen=clen, rmsd=rms, cov=cov, score=score, dist=dist, bc=bc, ok=0))
}

