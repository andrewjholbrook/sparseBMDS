# define arguments
args = commandArgs(trailingOnly=TRUE)
print(args)

items <- as.integer(args[1])
band <- as.double(args[2])
sd_opt <- as.double(args[3])
iter <- as.double(args[4])
thinner <- as.double(args[5])
burnin <- as.double(args[6])

# oFile <- paste("errresults/", "sd_", sd_opt, "band#_", band, ".txt", sep="")

# import libraries
library(MassiveMDS)
library(invgamma)
library(mvtnorm)
library(tidyverse)
library(truncnorm)

# load necessary functions
#source("~/sBMDS_scripts/sBMDS_MassiveMDS_functions.R")
source("code/Rscripts/sBMDS_MassiveMDS_functions.R")

# set seed
set.seed(12345)

#mse_band_hmc <- c()

# simulate data
Z.sim <- mvtnorm::rmvnorm(items, mean = rep(0, 2), sigma = diag(2)) 
D.sim <- as.matrix(dist(Z.sim)) # n x n matrix of distances
# noise <- matrix(rtruncnorm(items^2, mean = 0, sd = sd_opt, a = 0), items, items) 
# noise <- (noise + t(noise)) / 2
# diag(noise) <- 0
noise <- matrix(0, nrow = items, ncol = items)
noise[upper.tri(noise)] <- rtruncnorm(items * (items - 1) / 2, mean = 0, 
                                      sd = sd_opt, a = 0)
noise <- (noise + t(noise))
D.noise <- noise + D.sim

# run simulation

# sim <- sbmds_nonadapt_sigma_hmc(maxIts = iter, dims = 2, data = D.noise, 
#                                 precision = 1/(sd_opt^2), bandwidth = band, 
#                                 targetAccept = 0.65, stepSize = .1, thin = thinner)

sim <- bmds_hmc(maxIts = iter, dims = 2, data = D.noise, bandwidth = band,
                targetAccept = .65, targetAccept_Sigma = 0.44, 
                stepSize = .01, stepSizeSigma = .001, 
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

rel_sigma <- c()
for (t in 1:length(sim$sigma)){
  #subset <- sample(1000)
  rel_sigma[t] <- abs(sim$sigma[t] - sd_opt) / (sd_opt)
}

#mse_band_hmc <-  mean(error) #mean(error[-(1:burnin)])

error_results <- tibble(sd = sd_opt, locations = items, band.no = band, 
                        mse = mean(error), rel_sigma = median(rel_sigma),
                        AR = sim$AR, ARS = sim$ARS)

#write.csv(df_results, file = oFile, sep = ",", row.names = FALSE)

write_csv(error_results, append = TRUE, file = "~/sBMDS_results/mse_results_wsigma.txt")
