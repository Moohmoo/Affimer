`.sourceCpp_5_DLLInfo` <- dyn.load('/home/concorde/bochinski/Documents/Projet/Patch_Search/Tmp/sourceCpp-x86_64-pc-linux-gnu-1.0.8.3/sourcecpp_8a1597b40b994/sourceCpp_6.so')

distloc <- Rcpp:::sourceCppFunction(function(X, ind1, ind2) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_distloc')
mapping_dist_sum <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_mapping_dist_sum')
mapping_dist_sum2 <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_mapping_dist_sum2')
getConnectedNeighborsc <- Rcpp:::sourceCppFunction(function(C, V, X, Y, thresh, nbnei) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_getConnectedNeighborsc')
vertexc <- Rcpp:::sourceCppFunction(function(XProp, YProp) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_vertexc')
buildGraphc <- Rcpp:::sourceCppFunction(function(X, Y, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_buildGraphc')
buildGraph2c <- Rcpp:::sourceCppFunction(function(X, XRes, Y, YRes, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_buildGraph2c')

rm(`.sourceCpp_5_DLLInfo`)
