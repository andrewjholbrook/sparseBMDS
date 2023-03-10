# define arguments
args = commandArgs(trailingOnly=TRUE)
print(args)

items <- as.double(args[1])
band <- as.integer(args[2])
sd_opt <- as.double(args[3])
iter <- as.double(args[4])
burnin <- 20000

# oFile <- paste("errresults/", "sd_", sd_opt, "band#_", band, ".txt", sep="")

# import libraries
library(MassiveMDS)
library(invgamma)
library(mvtnorm)
library(truncnorm)

# load necessary functions
#source("~/sBMDS_scripts/sBMDS_MassiveMDS_functions.R")
source("code/Rscripts/sBMDS_MassiveMDS_functions.R")

# set seed
set.seed(12345)

mse_band_hmc <- c()

# simulate data
Z.sim <- mvtnorm::rmvnorm(items, mean = rep(0, 2), sigma = diag(2)) 
D.sim <- as.matrix(dist(Z.sim)) # n x n matrix of distances
noise <- matrix(rtruncnorm(items^2, mean = 0, sd = sd_opt), items, items) 
noise <- (noise + t(noise)) / 2
diag(noise) <- 0
D.noise <- noise + D.sim

# run simulation

for (s in 1:band){
  sim <- sbmds_nonadapt_sigma_hmc(maxIts = iter, dims = 2, data = D.noise, 
                                  precision = 1/(sd_opt^2), bandwidth = s, 
                                  targetAccept = 0.238, stepSize = .1)
  error <- c()
  for (t in 1:iter){
    error[t] <- mean((as.matrix(dist(sim[t, , ])) - D.sim)^2)
  }
  mse_band_hmc[s] <- mean(error[-(1:burnin)])
}

error_results <- tibble(sd = sd_opt, locations = items, band.no = seq(band), 
                        mse = mse_band_hmc)

#write.csv(df_results, file = oFile, sep = ",", row.names = FALSE)

write_csv(error_results, append = TRUE, file = "~/sBMDS_results/mse_results.txt")
