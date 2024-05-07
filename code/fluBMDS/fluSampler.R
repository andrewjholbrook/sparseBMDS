# input arguments
args = commandArgs(trailingOnly=TRUE)
print(args)
virus <- as.character(args[1]) # options = h1, h3, vic, yam
band <- as.integer(args[2])
iteration <- as.integer(args[3])
burnin <- as.double(args[4])
landmark <- (as.character(args[5]) == "TRUE")

# import libraries
library(MassiveMDS)
library(mvtnorm)
library(truncnorm)
library(tidyverse)

# necessary functions
source("flufunc.R")

# read tree file
readbeast <- function(virus, priorRootSampleSize = 0.001) {
  # Read BEAST tree
  tree <- ape::read.nexus(paste0(virus, "_small_sample.trees"))
  N <- length(tree[[1]]$tip.label)
  
  # Get tree VCV conditional on root
  treeLength <- sum(tree[[1]]$edge.length)
  treeVcv <- caper::VCV.array(tree[[1]]) / treeLength
  class(treeVcv) <- "matrix"
  
  # Integrate out fully conjugate root [??]
  treeVcv <- treeVcv + matrix(1 / priorRootSampleSize, ncol = N, nrow = N)
  
  return(list(treeVcv = treeVcv))
}
tree_info <- readbeast(virus = virus)

# read distance matrix with all viral info
fluCombi_Deff <- read.delim("fluCombi_Deff.txt", header = TRUE, row.names = 1)

# select appropriate subset based on virus 
if (virus == "h1"){
  distMat <- fluCombi_Deff[1:1370, 1:1370]
} else if (virus == "h3"){
  distMat <- fluCombi_Deff[1371:2759, 1371:2759]
} else if (virus == "vic"){
  distMat <- fluCombi_Deff[2760:4152, 2760:4152]
} else {
  distMat <- fluCombi_Deff[4153:5392, 4153:5392]
}

colnames(distMat) <- rownames(distMat)

# randomize data 
set.seed(666)
rand <- sample(1:nrow(distMat), nrow(distMat))
tree_info$treeVcvRand <- tree_info$treeVcv[rand, rand]
distMatRand <- distMat[rand, rand]

# run sampler
output <- hmcsampler(n_iter = iteration, burnIn = burnin, data = distMatRand, 
                     virus = virus, beast = tree_info, latentDimension = 6, 
                     StepSize = 0.01, StepSizeSigma = 0.01, 
                     mdsPrecision = 5, sparse_band = band, 
                     targetAccept = 0.65, targetAcceptSigma = 0.44, thin = 100,
                     landmarks = landmark)

# create distance matrix
# size <- nrow(distMat)
# distArray <- array(NA, dim = c(length(output$samples), size, size))
# for (i in 1:dim(distArray)[1]){
#   distArray[i, , ] <- as.matrix(dist(output$samples[[i]]))
# }

# calculate miness 
# miness <- 100000
# for (k in 2:size){
#   alless <- coda::effectiveSize(distArray[, k, 1:(k - 1)])
#   ess <- min(alless)
#   if (ess < miness){
#     miness <- ess
#   }
# }
# print(miness)

flu_comp <- tibble(virus = virus, band = band, landmark = landmark, 
                   time = output$time[3])

write_csv(flu_comp, append = TRUE, file = "fluResults.txt")
