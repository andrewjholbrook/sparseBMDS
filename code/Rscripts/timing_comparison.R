# input arguments
args = commandArgs(trailingOnly=TRUE)
print(args)
items <- as.double(args[1])
method <- as.character(args[2])
band <- as.integer(args[3])
iteration <- as.integer(args[4])
burnin <- as.double(args[5])
landmark <- (as.character(args[6]) == "TRUE")

# import libraries
library(MassiveMDS)
library(mvtnorm)
library(truncnorm)
library(tidyverse)

# source functions
#source("~/sparseBMDS/code/Rscripts/sBMDS_MassiveMDS_functions.R")
source("~/sBMDS_scripts/sBMDS_MassiveMDS_functions.R")

# simulate data
set.seed(12345)
latent <- rmvnorm(items, mean = rep(0, 2), sigma = diag(2))
D_latent <- as.matrix(dist(latent))

# adding noise to distance matrix
noise <- matrix(0, nrow = items, ncol = items)
noise[upper.tri(noise)] <- rtruncnorm(items * (items - 1) / 2, mean = 0, 
                                      sd = .2, a = 0)
noise <- noise + t(noise)
D_noise <- noise + D_latent

miness <- 0

# run simulation & time eval
#while (miness < 100){
  if (method == "FMH"){
    # Full likelihood Metropolis-Hastings (FMH)
    #start_time <- proc.time()
    sim_result <- bmds_metropolis(maxIts = iteration, dims = 2, data = D_noise, 
                                  bandwidth = items, precision = 25, 
                                  landmarks = landmark, targetAccept = 0.238, 
                                  stepSize = .01, thin = 100, burnin = burnin)
    #end_time <- proc.time() - start_time
  } else if (method == "FHMC"){
    # Full likelihood Hamiltonian Monte Carlo (FHMC)
    #start_time <- proc.time()
    sim_result <- bmds_nonadapt_sigma_hmc(maxIts = iteration, dims = 2, 
                                          data = D_noise, bandwidth = item, 
                                          precision = 25, landmarks = landmark,
                                          targetAccept = 0.65, stepSize = .01, 
                                          thin = 100, burnin = burnin)
    #end_time <- proc.time() - start_time
  } else if (method == "SMH"){
    # Sparse likelihood MH (SMH)
    #start_time <- proc.time()
    sim_result <- bmds_metropolis(maxIts = iteration, dims = 2, data = D_noise, 
                                  bandwidth = band, precision = 25, 
                                  landmarks = landmark, targetAccept = 0.238, 
                                  stepSize = .02, thin = 100, burnin = burnin)
    #end_time <- proc.time() - start_time
  } else {
    # Sparse likelihood & sparse gradient HMC 
    #start_time <- proc.time()
    sim_result <- bmds_nonadapt_sigma_hmc(maxIts = iteration, dims = 2, 
                                          data = D_noise, bandwidth = band, 
                                          precision = 25, landmarks = landmark,
                                          targetAccept = 0.65, stepSize = .02, 
                                          thin = 100, burnin = burnin)
    #end_time <- proc.time() - start_time
  }
  
  # calculate the distance matrix
  
  dist_sim <- array(0, dim = c(dim(sim_result$latent)[1], items, items))
  
  for (j in 1:dim(sim_result$latent)[1]){
    dist_sim[j, , ] <- as.matrix(dist(sim_result$latent[j, , ]))
  }
  
  # calculate minimum ess, only need to look at lower-triangle 
  miness_temp <- 110000
  medess <- c()
  for (k in 2:items){
    alless <- coda::effectiveSize(dist_sim[, k, 1:(k - 1)])
    ess <- min(alless)
    if (ess < miness_temp){
      miness_temp <- ess
    }
    medess[k - 1] <- median(alless)
  }
  miness <- miness_temp
  print(miness)
  # add more iterations
  #iteration <- iteration + 1000 
#}

# make a dataframe & write to csv file
int_min <- tibble(items = items, method = method, landmarks = landmark, 
                  time = sim_result$time[3], #iter = iteration - 1000, 
                  miness = miness, medess = mean(medess))

#write_csv(int_min, append = TRUE, file = "~/sBMDS_results/min_iter.txt")
write_csv(int_min, append = TRUE, file = "code/txt_simdata/time_bands.txt")
