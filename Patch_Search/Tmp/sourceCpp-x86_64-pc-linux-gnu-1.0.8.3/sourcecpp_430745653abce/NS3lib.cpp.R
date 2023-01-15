`.sourceCpp_7_DLLInfo` <- dyn.load('/home/concorde/bochinski/Documents/Recherche_Pattern/Patch_Search/Tmp/sourceCpp-x86_64-pc-linux-gnu-1.0.8.3/sourcecpp_430745653abce/sourceCpp_13.so')

distloc <- Rcpp:::sourceCppFunction(function(X, ind1, ind2) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_distloc')
mapping_dist_sum <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_mapping_dist_sum')
mapping_dist_sum2 <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_mapping_dist_sum2')
getConnectedNeighborsc <- Rcpp:::sourceCppFunction(function(C, V, X, Y, thresh, nbnei) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_getConnectedNeighborsc')
vertexc <- Rcpp:::sourceCppFunction(function(XProp, YProp) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_vertexc')
buildGraphc <- Rcpp:::sourceCppFunction(function(X, Y, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_buildGraphc')
buildGraph2c <- Rcpp:::sourceCppFunction(function(X, XRes, Y, YRes, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_buildGraph2c')

rm(`.sourceCpp_7_DLLInfo`)
