# input arguments
args = commandArgs(trailingOnly=TRUE)
print(args)
item <- as.double(args[1])
method <- as.integer(args[2])
band <- as.integer(args[3])

# import libraries 
library(MassiveMDS)
library(mvtnorm)
library(truncnorm)
library(plyr)
library(tidyverse)

# source functions
source("~/sparseBMDS/code/Rscripts/sBMDS_MassiveMDS_functions.R")
#source("~/sBMDS_scripts/sBMDS_MassiveMDS_functions.R")

# simulate data
set.seed(12345)
Z_ex <- rmvnorm(item, mean = rep(0, 2), sigma = diag(2))
D_ex <- as.matrix(dist(Z_ex))

# adding noise to distance matrix
noise <- matrix(rtruncnorm(item^2, mean = 0, sd = .1), item, item) 
noise <- (noise + t(noise)) / 2
diag(noise) <- 0
D_noise <- noise + D_ex

# run simulation & time eval

if (method == "FMH"){
  # Full likelihood Metropolis-Hastings (FMH)
  start_time <- proc.time()
  sim_result <- sbmds_metropolis(maxIts = 100000, dims = 2, data = D_noise, 
                                 bandwidth = item, precision = 100, 
                                 targetAccept = 0.238, stepSize = .1, thin = 100)
  end_time <- proc.time() - start_time
} else if (method == "FHMC"){
  # Full likelihood Hamiltonian Monte Carlo (FHMC)
  start_time <- proc.time()
  sim_result <- sbmds_nonadapt_sigma_hmc(maxIts = 100000, dims = 2, data = D_noise, 
                                         bandwidth = item, precision = 100, 
                                         targetAccept = 0.65, stepSize = .1, thin = 100)
  end_time <- proc.time() - start_time
} else if (method == "SMH"){
  # Sparse likelihood MH (SMH)
  start_time <- proc.time()
  sim_result <- sbmds_metropolis(maxIts = 100000, dims = 2, data = D_noise, 
                                 bandwidth = band, precision = 100, 
                                 targetAccept = 0.238, stepSize = .1, thin = 100)
  end_time <- proc.time() - start_time
} else if (method == "SSHMC") {
  # Sparse likelihood & sparse gradient HMC 
  start_time <- proc.time()
  sim_result <- sbmds_nonadapt_sigma_hmc(maxIts = 100000, dims = 2, data = D_noise, 
                                         bandwidth = band, precision = 100, 
                                         targetAccept = 0.65, stepSize = .1, thin = 100)
  end_time <- proc.time() - start_time
} else {
  # Full likelihood & sparse gradient HMC 
  start_time <- proc.time()
  sim_result <- fsbmds_nonadapt_sigma_hmc(maxIts = 100000, dims = 2, data = D_noise, 
                                         bandwidth = band, precision = 100, 
                                         targetAccept = 0.65, stepSize = .1, thin = 100)
  end_time <- proc.time() - start_time
}
 
# calculate the distance matrix

simdist_mat <- array(0, dim = c(1000, item, item))

for (j in 2:1001){
  simdist_mat[j - 1, , ] <- as.matrix(dist(sim_result[j, , ]))
}

# calculate minimum ess
miness <- 100000
medess <- c()
for (k in 1:item){
  alless <- coda::effectiveSize(simdist_mat[-(1:200), -k, k])
  ess <- min(alless) # minimum ess
  if (ess < miness){
    miness <- ess
  }
  medess[k] <- median(alless) # median ess 
}

# efficiency 

mineff <- miness / end_time[3] # effective sample size per second
medeff <- mean(medess) / end_time[3]

# make a dataframe & write to csv file

df_effcomp <- tibble(item = item, method = method, time = end_time[3], 
                     miness = miness, medess = mean(medess), minratio = mineff,
                     medratio = medeff)

write_csv(df_effcomp, append = TRUE, file = "efficiency_results.txt")
