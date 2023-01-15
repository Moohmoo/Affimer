`.sourceCpp_2_DLLInfo` <- dyn.load('/mnt/c/Users/oussa/Downloads/Master/Projet Long/Mohamed/Patch_Search/Tmp/sourceCpp-x86_64-pc-linux-gnu-1.0.9/sourcecpp_2fca2f2be7ff/sourceCpp_3.so')

distloc <- Rcpp:::sourceCppFunction(function(X, ind1, ind2) {}, FALSE, `.sourceCpp_2_DLLInfo`, 'sourceCpp_2_distloc')
mapping_dist_sum <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_2_DLLInfo`, 'sourceCpp_2_mapping_dist_sum')
mapping_dist_sum2 <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_2_DLLInfo`, 'sourceCpp_2_mapping_dist_sum2')
getConnectedNeighborsc <- Rcpp:::sourceCppFunction(function(C, V, X, Y, thresh, nbnei) {}, FALSE, `.sourceCpp_2_DLLInfo`, 'sourceCpp_2_getConnectedNeighborsc')
vertexc <- Rcpp:::sourceCppFunction(function(XProp, YProp) {}, FALSE, `.sourceCpp_2_DLLInfo`, 'sourceCpp_2_vertexc')
buildGraphc <- Rcpp:::sourceCppFunction(function(X, Y, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_2_DLLInfo`, 'sourceCpp_2_buildGraphc')
buildGraph2c <- Rcpp:::sourceCppFunction(function(X, XRes, Y, YRes, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_2_DLLInfo`, 'sourceCpp_2_buildGraph2c')

rm(`.sourceCpp_2_DLLInfo`)
