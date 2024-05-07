# define arguments to be changed

args = commandArgs(trailingOnly=TRUE)
print(args)
item <- as.double(args[1])
sigma2 <- as.integer(args[2])
# oFile <- paste("errresults/", "sd_", sd_opt, "band#_", band, ".txt", sep="")

# import libraries
library(tidyverse)
library(invgamma)
library(mvtnorm)
library(truncnorm)
#library(gridExtra)

# load necessary functions
source("sBMDS_functions.R")

# set seed

# make functions for evaluating computational time
set.seed(12345)

comp_time <- function(n, dims, sigmasq){
  ll_time <- c()
  grad_time <- c()
  Z_ex <- rmvnorm(n, mean = rep(0, dims), sigma = diag(dims))
  D_ex <- as.matrix(dist(Z_ex, diag = TRUE, upper = TRUE))
  for (i in 1:n){
    # likelihood
    stime <- proc.time()
    sll_results <- sbmds_ll(latent = Z_ex, D = D_ex, 
                            sigmasq, band.no = i)
    endtime <- proc.time() - stime
    ll_time[i] <- endtime[3]
    
    # gradient
    stime <- proc.time()
    sgrad_results <- sbmds_grad(latent = Z_ex, D = D_ex, 
                                sigmasq = sigma2, band.no = i)
    endtime <- proc.time() - stime
    grad_time[i] <- endtime[3]
  }
  return(list(ll_time, grad_time))
}

# results

comp_results <- comp_time(n = item, dims = 2, sigmasq = sigma2)

comp_df <- tibble(N = item, band.no = seq(1, item), sigmasq = sqrt(sigma2),
                     ll_time = comp_results[[1]],
                     grad_time = comp_results[[2]])

#write.csv(df_results, file = oFile, sep = ",", row.names = FALSE)

write_csv(comp_df, append = TRUE, file = "~/sBMDS_results/raw_like_grad_eval.txt")