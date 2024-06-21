# Sparse BMDS (sBMDS)

The folder `code` contains five different folders. We briefly describe each one and highlight important files:

1. `Rscript`:
`sBMDS_MassiveMDS_functions.R`: functions to compute sBMDS posterior distribution and MCMC sampling using the `MassiveMDS` package. This R package efficiently evaluates the sBMDS log-likelihoods and gradients and integraes directly into BEAST.
`sBMDS_functions.R`: functions to compute sBMDS log-likelihood, gradient, posterior, and MCMC algorithms directly in R
`MassiveMDS_mse.R`: calculate MSE and relative error for latent variables and sigma respectively
`MassiveMDS_rawllgrad.R`: give the computational time to evaluate the raw log-likelihood and gradient in the `MassiveMDS` package. The `MassiveMDS` package 
`R_rawllgrad_eval.R`: same as above but uses the `sBMDS_functions.R` script in R

2. `data_graphs`: simulation output used for the plots in the sBMDS manuscript
   
3. `fluBMDS`: sBMDS application
   
4. `results_graphs`: the plots in the sBMDS manuscript
   
5. `rmd_files`: code used to generate the plots 



