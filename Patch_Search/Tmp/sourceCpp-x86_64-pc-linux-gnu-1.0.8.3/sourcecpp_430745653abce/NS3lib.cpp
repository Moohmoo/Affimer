// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <ctime>
#include <math.h>
using namespace Rcpp;
const double EPS =  1.e-8;
const double MAXITER =  200;


// [[Rcpp::export]]
double distloc(NumericMatrix X, int ind1, int ind2) {
  int i;
  double x,y,sum=0;
  for (i=0;i<3;i++) {
    x=X(ind1-1,i);
    y=X(ind2-1,i);
    sum+=(x-y)*(x-y);
  }
  return sqrt(sum);
}
// [[Rcpp::export]]
NumericVector mapping_dist_sum(NumericMatrix X, IntegerVector XRes, IntegerVector I, IntegerVector CI, NumericMatrix Y, IntegerVector YRes, IntegerVector J, IntegerVector CJ, double thresh) {
  double  score, sc, d1, d2, deltadist;
  int i,j,n;

  n=I.size();
  //printf("mapping_dist_sum %d %d\n", I.size(), CI.size());
  if (thresh<=0.) thresh=15.; //patchsearch_NEIDIST
  NumericVector scores(n);
  for (i=0; i<I.size();i++) {
    score=0.;
    for (j=0; j<CI.size();j++) {
       if ((CJ(j)==J(i) && CI(j)!=I(i)) || (CI(j)==I(i) && CJ(j)!=J(i))) {
        score=0.;
        break;
      } else {
	d1=distloc(Y,CJ(j),J(i));
	d2=distloc(X,CI(j),I(i));
	deltadist=fabs(d1-d2);
	sc=1-deltadist/thresh;
      }
      //if (deltadist>thresh) sc=0.;
      score=score+sc;
    }
    //printf("scores[%d]=%lf\n", i, score);
    scores[i]=score/CI.size();
  }
  return scores;
}


// [[Rcpp::export]]
NumericVector mapping_dist_sum2(NumericMatrix X, IntegerVector XRes, IntegerVector I, IntegerVector CI, NumericMatrix Y, IntegerVector YRes, IntegerVector J, IntegerVector CJ, double thresh) {
  double  score, sc, d1, d2, deltadist;
  int i,j,n;

  n=I.size();
  if (thresh<=0.) thresh=15.; //patchsearch_NEIDIST
  NumericVector scores(n);
  for (i=0; i<I.size();i++) {
    score=0.;
    for (j=0; j<CI.size();j++) {
       if ((CJ(j)==J(i) && CI(j)!=I(i)) || (CI(j)==I(i) && CJ(j)!=J(i))) {
        score=0.;
        break;
      } else {
	d1=distloc(Y,CJ(j),J(i));
	d2=distloc(X,CI(j),I(i));
	deltadist=fabs(d1-d2);
	sc=1-deltadist/thresh;
      }
      //if (deltadist>thresh) sc=0.;
      score=score+sc;
    }
    //printf("scores[%d]=%lf\n", i, score);
    scores[i]=score;
  }
  return scores;
}

// all nodes in V connected to at least nbnei in C
// [[Rcpp::export]]
IntegerVector getConnectedNeighborsc(IntegerVector C, IntegerMatrix V, NumericMatrix X, NumericMatrix Y, double thresh, int nbnei) {
  int i,j,inC,count=0;
  int nC=C.size();
  IntegerVector K;
  IntegerVector I(nC), J(nC);
  int N=X.nrow();

  for (int i=0; i < nC; i++) {
    I(i)=(C(i)-1)%N+1;
    J(i)=(int)(C(i)-1)/N+1;
  }

  double deltadist, d1, d2;
  for (j=0; j<V.nrow(); j++) {
    count=0;
    inC=0;
    for (i=0; i<nC; i++) {
      if (J(i)!=V(j,0) || I(i)!=V(j,1)) {
	d1=distloc(Y,J(i),V(j,0));
	d2=distloc(X,I(i),V(j,1));
	deltadist=fabs(d1-d2);
	//deltadist=fabs(d1-d2)/std::min(d1,d2);
	if (deltadist<=thresh) count++;
      }
      else inC=1; // j is in C
    }
    if (count>=nbnei && !inC) K.push_back((V(j,0)-1)*N+V(j,1));
  }
  return K;
}

// (V(v,0),V(v,1)): link number v (query, target) (!:base 1)
//' Constructs vertex table V
//'
//' @param XProp
//' @param YProp
//' @export
// [[Rcpp::export]]
IntegerMatrix vertexc(CharacterVector XProp,  CharacterVector YProp) {
  int i,ip,v=0;
  int N=XProp.size();
  int M=YProp.size();
  IntegerMatrix V(N*M,2);
  for(i=0; i<M; i++)
    for (ip=0; ip<N; ip++)
      if (strcmp(XProp(ip),YProp(i))==0) { 
	V(v,0)=i+1;
	V(v,1)=ip+1;
	v++;
      }
  if (v<=1) {
    fprintf(stderr,"no product graph...\n");
    return V(Range(0,0),_);
  }
  return V(Range(0,v-1),_);
}

//' Constructs product graph: absolute delta_dist
//'
//' @param X matrix Nx3
//' @param Y matrix Nx3
//' @param V matrix Nvx2
//' @export
// [[Rcpp::export]]
IntegerMatrix buildGraphc(NumericMatrix X, NumericMatrix Y, IntegerMatrix V,  double thresh, double mindist, double maxdist) {
  int i,j, nv, e=0;
  double d,d1,d2;
  std::vector<int> Etmp;

  nv=V.nrow();
  for(i=0; i<nv-1; i++) {
    for (j=i+1; j<nv; j++) {
      //printf("%d %d %d %d\n",V(i,0),V(j,0),V(i,1),V(j,1));
      d1=distloc(Y,V(i,0),V(j,0)); // dist between atoms i and j in Y
      d2=distloc(X,V(i,1),V(j,1)); // dist between atoms i and j in X
      d=fabs(d1-d2);
      if (d<=thresh && d1>=mindist && d2>=mindist && d1<=maxdist && d2<=maxdist) {
	      Etmp.push_back(V(i,1));
	      Etmp.push_back(V(j,1));
	      Etmp.push_back(V(i,0));
	      Etmp.push_back(V(j,0));
	      e+=1;
	}
    }
  }
  int n=Etmp.size();
  IntegerMatrix E(e,4);
  for (i=0,j=0; i<n-3; i+=4,j++) {
    E(j,0)=Etmp[i];
    E(j,1)=Etmp[i+1];
    E(j,2)=Etmp[i+2];
    E(j,3)=Etmp[i+3];
  }
  Etmp.clear();
return E;
}

//' Constructs product graph: absolute delta_dist
//'
//' @param X matrix Nx3
//' @param Y matrix Nx3
//' @param V matrix Nvx2
//' @export
// [[Rcpp::export]]
IntegerMatrix buildGraph2c(NumericMatrix X, IntegerVector XRes, NumericMatrix Y, IntegerVector YRes, IntegerMatrix V,  double thresh, double mindist, double maxdist) {
  int i,j, tres, qres, nv, e=0;
  double d,d1,d2;
  std::vector<int> Etmp;
  
  nv=V.nrow();
  for(i=0; i<nv-1; i++) {
    for (j=i+1; j<nv; j++) {
      //printf("%d %d %d %d\n",V(i,0),V(j,0),V(i,1),V(j,1));
      d1=distloc(Y,V(i,0),V(j,0)); // dist between atoms i and j in Y
      d2=distloc(X,V(i,1),V(j,1)); // dist between atoms i and j in X
      d=fabs(d1-d2);
      qres=(YRes(V(i,0)-1)==YRes(V(j,0)-1));
      tres=(XRes(V(i,1)-1)==XRes(V(j,1)-1));
      if (d<=thresh && d1>=mindist && d2>=mindist && d1<=maxdist && d2<=maxdist && !(qres || tres)) {
	      Etmp.push_back(V(i,1));
	      Etmp.push_back(V(j,1));
	      Etmp.push_back(V(i,0));
	      Etmp.push_back(V(j,0));
	      e+=1;
	}
    }
  }
  int n=Etmp.size();
  IntegerMatrix E(e,4);
  for (i=0,j=0; i<n-3; i+=4,j++) {
    E(j,0)=Etmp[i];
    E(j,1)=Etmp[i+1];
    E(j,2)=Etmp[i+2];
    E(j,3)=Etmp[i+3];
  }
  Etmp.clear();
return E;
}



#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// distloc
double distloc(NumericMatrix X, int ind1, int ind2);
RcppExport SEXP sourceCpp_7_distloc(SEXP XSEXP, SEXP ind1SEXP, SEXP ind2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type ind1(ind1SEXP);
    Rcpp::traits::input_parameter< int >::type ind2(ind2SEXP);
    rcpp_result_gen = Rcpp::wrap(distloc(X, ind1, ind2));
    return rcpp_result_gen;
END_RCPP
}
// mapping_dist_sum
NumericVector mapping_dist_sum(NumericMatrix X, IntegerVector XRes, IntegerVector I, IntegerVector CI, NumericMatrix Y, IntegerVector YRes, IntegerVector J, IntegerVector CJ, double thresh);
RcppExport SEXP sourceCpp_7_mapping_dist_sum(SEXP XSEXP, SEXP XResSEXP, SEXP ISEXP, SEXP CISEXP, SEXP YSEXP, SEXP YResSEXP, SEXP JSEXP, SEXP CJSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XRes(XResSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type I(ISEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type CI(CISEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YRes(YResSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type J(JSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type CJ(CJSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(mapping_dist_sum(X, XRes, I, CI, Y, YRes, J, CJ, thresh));
    return rcpp_result_gen;
END_RCPP
}
// mapping_dist_sum2
NumericVector mapping_dist_sum2(NumericMatrix X, IntegerVector XRes, IntegerVector I, IntegerVector CI, NumericMatrix Y, IntegerVector YRes, IntegerVector J, IntegerVector CJ, double thresh);
RcppExport SEXP sourceCpp_7_mapping_dist_sum2(SEXP XSEXP, SEXP XResSEXP, SEXP ISEXP, SEXP CISEXP, SEXP YSEXP, SEXP YResSEXP, SEXP JSEXP, SEXP CJSEXP, SEXP threshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XRes(XResSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type I(ISEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type CI(CISEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YRes(YResSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type J(JSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type CJ(CJSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    rcpp_result_gen = Rcpp::wrap(mapping_dist_sum2(X, XRes, I, CI, Y, YRes, J, CJ, thresh));
    return rcpp_result_gen;
END_RCPP
}
// getConnectedNeighborsc
IntegerVector getConnectedNeighborsc(IntegerVector C, IntegerMatrix V, NumericMatrix X, NumericMatrix Y, double thresh, int nbnei);
RcppExport SEXP sourceCpp_7_getConnectedNeighborsc(SEXP CSEXP, SEXP VSEXP, SEXP XSEXP, SEXP YSEXP, SEXP threshSEXP, SEXP nbneiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int >::type nbnei(nbneiSEXP);
    rcpp_result_gen = Rcpp::wrap(getConnectedNeighborsc(C, V, X, Y, thresh, nbnei));
    return rcpp_result_gen;
END_RCPP
}
// vertexc
IntegerMatrix vertexc(CharacterVector XProp, CharacterVector YProp);
RcppExport SEXP sourceCpp_7_vertexc(SEXP XPropSEXP, SEXP YPropSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type XProp(XPropSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type YProp(YPropSEXP);
    rcpp_result_gen = Rcpp::wrap(vertexc(XProp, YProp));
    return rcpp_result_gen;
END_RCPP
}
// buildGraphc
IntegerMatrix buildGraphc(NumericMatrix X, NumericMatrix Y, IntegerMatrix V, double thresh, double mindist, double maxdist);
RcppExport SEXP sourceCpp_7_buildGraphc(SEXP XSEXP, SEXP YSEXP, SEXP VSEXP, SEXP threshSEXP, SEXP mindistSEXP, SEXP maxdistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< double >::type mindist(mindistSEXP);
    Rcpp::traits::input_parameter< double >::type maxdist(maxdistSEXP);
    rcpp_result_gen = Rcpp::wrap(buildGraphc(X, Y, V, thresh, mindist, maxdist));
    return rcpp_result_gen;
END_RCPP
}
// buildGraph2c
IntegerMatrix buildGraph2c(NumericMatrix X, IntegerVector XRes, NumericMatrix Y, IntegerVector YRes, IntegerMatrix V, double thresh, double mindist, double maxdist);
RcppExport SEXP sourceCpp_7_buildGraph2c(SEXP XSEXP, SEXP XResSEXP, SEXP YSEXP, SEXP YResSEXP, SEXP VSEXP, SEXP threshSEXP, SEXP mindistSEXP, SEXP maxdistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type XRes(XResSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type YRes(YResSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< double >::type mindist(mindistSEXP);
    Rcpp::traits::input_parameter< double >::type maxdist(maxdistSEXP);
    rcpp_result_gen = Rcpp::wrap(buildGraph2c(X, XRes, Y, YRes, V, thresh, mindist, maxdist));
    return rcpp_result_gen;
END_RCPP
}
