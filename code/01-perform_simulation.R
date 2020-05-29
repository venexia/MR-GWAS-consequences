remove(list = ls())
graphics.off()

# Specify paths

source("code/specify_paths.R", echo = TRUE)

# Source functions

source("code/fn-simulation.R", echo = TRUE)

# Run simulation for positive covariate and positive treatment

c_pos_t_pos <- simulation(name = "c_pos_t_pos",
                          seed = 1,
                          N = 500000,
                          snps = 200,
                          alpha1 = 10,
                          alpha2 = 0.3,
                          beta1 = log(2),
                          beta2 = log(1),
                          treat = 2,
                          c_mean = 10,
                          c_sd = 2)

# Run simulation for positive covariate and negative treatment

c_pos_t_neg <- simulation(name = "c_pos_t_neg",
                          seed = 1,
                          N = 500000,
                          snps = 200,
                          alpha1 = 10,
                          alpha2 = 0.3,
                          beta1 = log(2),
                          beta2 = log(1),
                          treat = -2,
                          c_mean = 10,
                          c_sd = 2)

# Run simulation for negativecovariate and positive treatment

c_neg_t_pos <- simulation(name = "c_neg_t_pos",
                          seed = 12345,
                          N = 500000,
                          snps = 200,
                          alpha1 = 10,
                          alpha2 = 0.3,
                          beta1 = log(2),
                          beta2 = log(1),
                          treat = 2,
                          c_mean = -10,
                          c_sd = 2)

# Run simulation for negative covariate and negative treatment

c_neg_t_neg <- simulation(name = "c_neg_t_neg",
                          seed = 12345,
                          N = 500000,
                          snps = 200,
                          alpha1 = 10,
                          alpha2 = 0.3,
                          beta1 = log(2),
                          beta2 = log(1),
                          treat = -2,
                          c_mean = -10,
                          c_sd = 2)