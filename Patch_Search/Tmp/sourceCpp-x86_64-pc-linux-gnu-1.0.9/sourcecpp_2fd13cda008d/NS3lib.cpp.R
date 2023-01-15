`.sourceCpp_4_DLLInfo` <- dyn.load('/mnt/c/Users/oussa/Downloads/Master/Projet Long/Mohamed/Patch_Search/Tmp/sourceCpp-x86_64-pc-linux-gnu-1.0.9/sourcecpp_2fd13cda008d/sourceCpp_9.so')

distloc <- Rcpp:::sourceCppFunction(function(X, ind1, ind2) {}, FALSE, `.sourceCpp_4_DLLInfo`, 'sourceCpp_4_distloc')
mapping_dist_sum <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_4_DLLInfo`, 'sourceCpp_4_mapping_dist_sum')
mapping_dist_sum2 <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_4_DLLInfo`, 'sourceCpp_4_mapping_dist_sum2')
getConnectedNeighborsc <- Rcpp:::sourceCppFunction(function(C, V, X, Y, thresh, nbnei) {}, FALSE, `.sourceCpp_4_DLLInfo`, 'sourceCpp_4_getConnectedNeighborsc')
vertexc <- Rcpp:::sourceCppFunction(function(XProp, YProp) {}, FALSE, `.sourceCpp_4_DLLInfo`, 'sourceCpp_4_vertexc')
buildGraphc <- Rcpp:::sourceCppFunction(function(X, Y, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_4_DLLInfo`, 'sourceCpp_4_buildGraphc')
buildGraph2c <- Rcpp:::sourceCppFunction(function(X, XRes, Y, YRes, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_4_DLLInfo`, 'sourceCpp_4_buildGraph2c')

rm(`.sourceCpp_4_DLLInfo`)
