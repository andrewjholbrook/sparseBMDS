# import libraries
library(MassiveMDS)
library(mvtnorm)
library(tidyverse)

# set seed 
set.seed(12345)

# inputs
items <- 10000
dimension <- 2
ll_time <- c()
grad_time <- c()
tau <- 100

# simulate data 
Z_ex <- rmvnorm(items, mean = rep(0, dimension), sigma = diag(dimension))
D_ex <- as.matrix(dist(Z_ex))

# evaluate computational time
for (i in 1:items){
  print(i)
  # create engine
  engine_eval <- MassiveMDS::createEngine(embeddingDimension = dimension, 
                                          locationCount = items, 
                                          truncation = TRUE, tbb = 0, simd = 0, 
                                          gpu = 0, single = 0, bandwidth = i)
  engine_eval <- MassiveMDS::setPairwiseData(engine_eval, data = D_ex)
  engine_eval <- MassiveMDS::setPrecision(engine_eval, precision = tau)
  engine_eval <- MassiveMDS::updateLocations(engine_eval, locations = Z_ex)
  
  # likelihood
  stime <- proc.time()
  MassiveMDS::getLogLikelihood(engine_eval, landmarks = FALSE)
  endtime <- proc.time() - stime
  ll_time[i] <- endtime[3]
  
  # gradient
  stime <- proc.time()
  MassiveMDS::getGradient(engine_eval, landmarks = FALSE)
  endtime <- proc.time() - stime
  grad_time[i] <- endtime[3]
}

# results
comp_df <- tibble(N = items, band.no = seq(1, items),
                  likelihood = ll_time,
                  gradient = grad_time)

write_csv(comp_df, append = TRUE, file = "~/sBMDS_results/rawllgrad_n10000B.txt")