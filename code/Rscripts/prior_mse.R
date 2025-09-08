# define arguments
args = commandArgs(trailingOnly=TRUE)
print(args)

items <- as.integer(args[1])
band <- as.double(args[2])
sd_opt <- as.double(args[3])
sd_prior <- as.double(args[4])
iter <- as.double(args[5])
thinner <- as.double(args[6])
burnin <- as.double(args[7])

# import libraries
library(MassiveMDS)
library(invgamma)
library(mvtnorm)
#library(tidyverse)
#library(dplyr)
library(truncnorm)
library(readr)

# load necessary functions
source("~/sBMDS_scripts/sBMDS_MassiveMDS_functions.R")
#source("code/Rscripts/sBMDS_MassiveMDS_functions.R")

# replace target & grad function to include new prior
s_target <- function(locations, engine, sigmasq, landmarks) {
  # update engine to calculate likelihood at location & sigmasq of input
  engine <- MassiveMDS::setPrecision(engine, precision = 1 / sigmasq)
  engine <- MassiveMDS::updateLocations(engine, locations)
  N <- engine$locationCount
  band.no <- engine$bandwidth
  m_full <- N * (N - 1 ) / 2
  m_sub <- N * band.no - band.no * (band.no + 1) / 2
  output <- MassiveMDS::getLogLikelihood(engine, landmarks) +
    # independent, Gaussian prior for theta centered at 0 & sd = sd_prior
    (m_sub / m_full) * 
    (sum(mvtnorm::dmvnorm(locations, mean = rep(0, engine$embeddingDimension), 
                         sigma = diag(rep(sd_prior^2, engine$embeddingDimension)),
                         log = TRUE)) +
    # inverse-gamma prior for sigma^2
    dinvgamma(sigmasq, shape = 1, rate = 1, log = TRUE))
  return(output)
}

sbmds_grad_sigma <- function(locations, engine, sigmasq, landmarks){
  # update engine to calculate gradient at location of input
  engine <- MassiveMDS::setPrecision(engine, precision = 1 / sigmasq)
  engine <- MassiveMDS::updateLocations(engine, locations)
  N <- engine$locationCount
  band.no <- engine$bandwidth
  m_full <- N * (N - 1 ) / 2
  m_sub <- N * band.no - band.no * (band.no + 1) / 2
  gradient <- MassiveMDS::getGradient(engine, landmarks) - 
    (m_sub / m_full) * (locations / sd_prior^2) # from the MVN prior
  return(gradient)
}

# set seed
set.seed(12345)

# simulate data
Z.sim <- mvtnorm::rmvnorm(items, mean = rep(0, 2), sigma = diag(2)) 
D.sim <- as.matrix(dist(Z.sim)) # n x n matrix of distances
noise <- matrix(0, nrow = items, ncol = items)
noise[upper.tri(noise)] <- rtruncnorm(items * (items - 1) / 2, mean = 0, 
                                      sd = sd_opt, a = 0)
noise <- noise + t(noise)
D.noise <- noise + D.sim

# run simulation

# sim <- sbmds_nonadapt_sigma_hmc(maxIts = iter, dims = 2, data = D.noise, 
#                                 precision = 1/(sd_opt^2), bandwidth = band, 
#                                 targetAccept = 0.65, stepSize = .1, thin = thinner)

sim <- bmds_hmc(maxIts = iter, dims = 2, data = D.noise, bandwidth = band,
                landmarks = TRUE, 
                targetAccept = .65, targetAccept_Sigma = 0.44, 
                stepSize = .01, stepSizeSigma = .01, 
                thin = thinner, burnin = burnin)

# calculate mse and relative error for latent variables and sigma respectively
error <- c()

for (t in 1:(dim(sim$latent)[1])){
  if (items > 1000){
    subset <- sample(1:items, 1000)
    error[t] <- mean((as.matrix(dist(sim$latent[t, subset, ])) - 
                        as.matrix(dist(Z.sim[subset, ])))^2)
  } else {
    error[t] <- mean((as.matrix(dist(sim$latent[t, , ])) - D.sim)^2)
  }
}

#rel_sigma <- c()
#for (t in 1:length(sim$sigma)){
#subset <- sample(1000)
#rel_sigma[t] <- abs(sim$sigma[t] - sd_opt) / (sd_opt)
#}

#mse_band_hmc <-  mean(error) #mean(error[-(1:burnin)])

error_results <- data.frame(sd = sd_prior, locations = items, band.no = band, 
                            method = "LM", mse = mean(error), 
                            rel_sigma = median(sim$sigma),
                            AR = sim$AR, ARS = sim$ARS)

write_csv(error_results, append = TRUE, 
          file = "~/sBMDS_results/mse_prior_weighted.txt")
