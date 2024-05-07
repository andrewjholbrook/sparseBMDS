# import libraries
library(MassiveMDS)
library(mvtnorm)
library(tidyverse)

# set seed 
set.seed(12345)

# input arguments
args = commandArgs(trailingOnly=TRUE)
print(args)
item <- as.double(args[1])
band <- as.integer(args[2])
landmark <- (as.character(args[3]) == "TRUE")

# pre-specified arguments
dimension <- 2
tau <- 100 # precision

# evaluate computational time for different item numbers at fixed bands
# simulate data 
Z_ex <- rmvnorm(item, mean = rep(0, dimension), sigma = diag(dimension))
D_ex <- as.matrix(dist(Z_ex))
  
# create engine
engine_eval <- MassiveMDS::createEngine(embeddingDimension = dimension, 
                                        locationCount = item, 
                                        truncation = TRUE, tbb = 0, simd = 0, 
                                        gpu = 0, single = 0, bandwidth = band)
engine_eval <- MassiveMDS::setPairwiseData(engine_eval, data = D_ex)
engine_eval <- MassiveMDS::setPrecision(engine_eval, precision = tau)
engine_eval <- MassiveMDS::updateLocations(engine_eval, locations = Z_ex)
  
# likelihood
stime <- proc.time()
MassiveMDS::getLogLikelihood(engine_eval, landmark)
endtime <- proc.time() - stime
ll_time <- endtime[3]
  
# gradient
stime <- proc.time()
MassiveMDS::getGradient(engine_eval, landmark)
endtime <- proc.time() - stime
grad_time <- endtime[3]

# results
comp_df <- tibble(N = item, band.no = band, 
                  likelihood = ll_time,
                  gradient = grad_time)

write_csv(comp_df, append = TRUE, file = "code/txt_simdata/rawllgrad_PUB.txt")




