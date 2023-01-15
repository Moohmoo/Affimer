`.sourceCpp_1_DLLInfo` <- dyn.load('/home/concorde/bochinski/Téléchargements/PS3.1/Tmp/sourceCpp-x86_64-pc-linux-gnu-1.0.8.3/sourcecpp_ff1413473576/sourceCpp_2.so')

distloc <- Rcpp:::sourceCppFunction(function(X, ind1, ind2) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_distloc')
mapping_dist_sum <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_mapping_dist_sum')
mapping_dist_sum2 <- Rcpp:::sourceCppFunction(function(X, XRes, I, CI, Y, YRes, J, CJ, thresh) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_mapping_dist_sum2')
getConnectedNeighborsc <- Rcpp:::sourceCppFunction(function(C, V, X, Y, thresh, nbnei) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_getConnectedNeighborsc')
vertexc <- Rcpp:::sourceCppFunction(function(XProp, YProp) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_vertexc')
buildGraphc <- Rcpp:::sourceCppFunction(function(X, Y, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_buildGraphc')
buildGraph2c <- Rcpp:::sourceCppFunction(function(X, XRes, Y, YRes, V, thresh, mindist, maxdist) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_buildGraph2c')

rm(`.sourceCpp_1_DLLInfo`)
