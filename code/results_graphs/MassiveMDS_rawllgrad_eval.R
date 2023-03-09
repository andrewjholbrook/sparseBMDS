# import libraries

library(MassiveMDS)
library(mvtnorm)
library(tidyverse)

# set seed 
set.seed(12345)

# inputs
items <- 10
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
  MassiveMDS::getLogLikelihood(engine_eval)
  endtime <- proc.time() - stime
  ll_time[i] <- endtime[3]
  
  # gradient
  stime <- proc.time()
  MassiveMDS::getGradient(engine_eval)
  endtime <- proc.time() - stime
  grad_time[i] <- endtime[3]
}

# results
comp_df <- tibble(N = items, band.no = seq(1, items), 
                  likelihood = ll_time,
                  gradient = grad_time)

#write.csv(df_results, file = oFile, sep = ",", row.names = FALSE)

write_csv(comp_df, append = TRUE, file = "code/txt_simdata/MassiveMDS_rawllgrad_eval.txt")

# plots
c1 <- comp_df %>% 
  gather(., key = "raw type", value = "elapsed time", c(likelihood, gradient)) %>%
  ggplot(., aes(band.no, `elapsed time`, color = `raw type`)) + 
  geom_point() +
  labs(title = "Raw Likelihood & Gradient Evaluations",
       x = "Number of Bands", y = "Computation Time (sec)",
       color = "") +
  theme_bw() +
  scale_color_manual(values = c("blue", "orange"),
                     labels = c("gradient", "likelihood"))

### heatmap data for likelihood

tau <- 100
items <- 10000
heatmap_ll <- matrix(NA, nrow = items, ncol = items)
heatmap_grad <- matrix(NA, nrow = items, ncol = items)

for (i in 1:items){
  print(i)
  # simulate data
  Z_i <- rmvnorm(i, mean = rep(0, dimension), sigma = diag(dimension))
  D_i <- as.matrix(dist(Z_i))
  for (j in 1:i){
    # create engine
    eng <- MassiveMDS::createEngine(embeddingDimension = dimension, 
                                            locationCount = i, 
                                            truncation = TRUE, tbb = 0, simd = 0, 
                                            gpu = 0, single = 0, bandwidth = j)
    eng <- MassiveMDS::setPairwiseData(eng, data = D_i)
    eng <- MassiveMDS::setPrecision(eng, precision = tau)
    eng <- MassiveMDS::updateLocations(eng, locations = Z_i)
    
    # compute likelihood time
    time <- proc.time()
    MassiveMDS::getLogLikelihood(eng)
    endtime <- proc.time() - time
    heatmap_ll[j, i] <- endtime[3]
    
    # compute gradient time
    time <- proc.time()
    MassiveMDS::getGradient(eng)
    endtime <- proc.time() - time
    heatmap_grad[j, i] <- endtime[3]
  }
}

## heatmap results
heatmap(heatmap_grad[1100:1150, 1100:1150], Rowv = NA, Colv = NA, 
        xlab = "Number of Locations", ylab = "Number of Bands", 
        main = "Gradient")

heatmap(heatmap_ll[1100:1150, 1100:1150], Rowv = NA, Colv = NA, 
        xlab = "Number of Locations", ylab = "Number of Bands",
        main = "Likelihood", na.rm = TRUE)

