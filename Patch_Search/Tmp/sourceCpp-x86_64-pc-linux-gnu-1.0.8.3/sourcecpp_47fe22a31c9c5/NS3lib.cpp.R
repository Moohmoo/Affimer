`.sourceCpp_3_DLLInfo` <- dyn.load('/home/concorde/bochinski/Documents/Projet/PS3.1/Tmp/sourceCpp-x86_64-pc-linux-gnu-1.0.8.3/sourcecpp_47fe22a31c9c5/sourceCpp_4.so')

distloc <- Rcpp:::sourceCppFunction(function(X, ind1, ind2) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_distloc')
mapping_dist_sum <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_mapping_dist_sum')
mapping_dist_sum2 <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_mapping_dist_sum2')
getConnectedNeighborsc <- Rcpp:::sourceCppFunction(function(C, V, X, Y, thresh, nbnei) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_getConnectedNeighborsc')
vertexc <- Rcpp:::sourceCppFunction(function(XProp, YProp) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_vertexc')
buildGraphc <- Rcpp:::sourceCppFunction(function(X, Y, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_buildGraphc')
buildGraph2c <- Rcpp:::sourceCppFunction(function(X, XRes, Y, YRes, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_buildGraph2c')

rm(`.sourceCpp_3_DLLInfo`)
