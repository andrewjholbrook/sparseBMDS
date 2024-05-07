library(ape)
library(caper)

### necessary functions
source("flufunc.R")

# read tree file
readbeast <- function(virus, priorRootSampleSize = 0.001) {
  # Read BEAST tree
  tree <- ape::read.nexus(paste0(virus, "_small_sample.trees"))
  N <- length(tree$PAUP_1$tip.label)
  
  # Get tree VCV conditional on root
  treeLength <- sum(tree$PAUP_1$edge.length)
  treeVcv <- caper::VCV.array(tree$PAUP_1) / treeLength
  class(treeVcv) <- "matrix"
  
  # Integrate out fully conjugate root [??]
  treeVcv <- treeVcv + matrix(1 / priorRootSampleSize, ncol = N, nrow = N)
  
  return(list(treeVcv = treeVcv))
}

### application 
# read tree file for h1n1
h1n1_tree_info <- readbeast("h1")

# permute treePrec
set.seed(666)
rand <- sample(1:1370, 1370)
h1n1_tree_info$treeVcvRand <- h1n1_tree_info$treeVcv[rand, rand]

# read distance matrix
#load("~/Downloads/distMat.Rdata")
fluCombi_Deff <- read.delim("fluCombi_Deff.txt", header = TRUE,
                            row.names = 1)
#h1n1_dist <- distMat[1:1370, 1:1370]
h1n1_dist <- fluCombi_Deff[1:1370, 1:1370]
colnames(h1n1_dist) <- rownames(h1n1_dist)
h1n1_distRand <- h1n1_dist[rand, rand]

# run sampler
hmc_sparse <- hmcsampler(n_iter = 10000, burnIn = 100, data = h1n1_distRand, 
                         beast = h1n1_tree_info, latentDimension = 6, 
                         StepSize = 10^-3, StepSizeSigma = 10^-3, 
                         mdsPrecision = 5, sparse_band = 20, virus = "h1.2",
                         targetAccept = 0.65, targetAcceptSigma = 0.44)
# read in files
chain <- read.csv("fluOutput/chain_otherh1B20.csv", 
                  col.names = c("time", "target", "prec"))
cov <- read.csv("fluOutput/latent_covh1B20.csv", header = FALSE)
coord <- read.csv("fluOutput/latent_coordh1B20.csv", header = FALSE)

# plots
plot(chain$target, type = "l"); coda::effectiveSize(chain$target)
plot(chain$prec, type = "l"); coda::effectiveSize(chain$prec)
plot(density(rowSums(cov) / 14.2))

dist_array <- array(0, dim = c(nrow(coord), 100, 100))
for (i in 1:dim(dist_array)[1]){
  coordmat <- matrix(coord[i, ], ncol = 6)
  sub <- sample(1:nrow(coordmat), 100)
  dist_array[i, , ] <- as.matrix(dist(coordmat[sub, ]))
}

miness_temp <- 110000
#medess <- c()
for (k in 2:100){
  alless <- coda::effectiveSize(dist_array[, k, 1:(k - 1)])
  ess <- min(alless)
  print(ess)
  if (ess < miness_temp){
    miness_temp <- ess
  }
  #medess[k - 1] <- median(alless)
}

maxess <- 110000
#medess <- c()
for (k in 2:100){
  alless <- coda::effectiveSize(dist_array[, k, 1:(k - 1)])
  ess <- max(alless)
  print(ess)
  if (ess > maxess){
    maxess <- ess
  }
  #medess[k - 1] <- median(alless)
}

