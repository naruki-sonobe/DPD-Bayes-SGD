#######################################
##    Load Necessary Libraries
#######################################
# install.packages("doParallel")
# install.packages("foreach")
# install.packages("ggplot2")
# install.packages("tictoc")
# install.packages("dqrng")

library(doParallel)  # For parallel computing
library(foreach)     # For 'foreach' loops
library(ggplot2)     # For plotting
library(tictoc)      # For timing
library(dqrng)       # For generating random numbers (dqrexp, dqrnorm)

# Clear the workspace
rm(list = ls())

#######################################
##        Helper Functions
#######################################

# Function to generate data with contamination
generate_data <- function(n = 100, contamination_rate = 0.01) {
  n_main <- ceiling(n * (1 - contamination_rate))  # Number of samples from the main distribution
  n_outliers <- n - n_main                         # Number of contaminated (outlier) samples
  
  x_main <- dqrnorm(n = n_main, mean = 0, sd = 1)
  x_outliers <- dqrnorm(n = n_outliers, mean = 10, sd = 0.1)
  
  return(list(
    n_main = n_main,
    n_outliers = n_outliers,
    x_main = x_main,
    x_outliers = x_outliers,
    x_all = c(x_main, x_outliers),
    contamination_rate = contamination_rate
  ))
}

# Calculate the density of N(mu, sigma^2) for a given theta = (mu, sigma)
calc_density <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  dnorm(x = x, mean = mu, sd = sigma)
}

# Calculate the score function
calc_score <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  cbind(
    (x - mu) / (sigma^2),
    ((x - mu)^2) / (sigma^3) - 1 / sigma
  )
}

# Regularization term (used in Exact LLB)
calc_regularizer <- function(theta, beta) {
  sigma <- theta[2]
  c(
    0,
    -(2 * pi)^(-beta / 2) * beta * (1 + beta)^(-3 / 2) * sigma^(-(beta + 1))
  )
}

# Generate samples Y ~ N(mu, sigma^2)
generate_y <- function(theta, m = 10) {
  mu <- theta[1]
  sigma <- theta[2]
  dqrnorm(n = m, mean = mu, sd = sigma)
}

#######################################
##        Main Part
#######################################

# Sample size
sample_size <- 1000
# Contamination rate
contamination_rate <- 0.05

# Generate data
data_list <- generate_data(n = sample_size, contamination_rate = contamination_rate)
observed_data <- data_list$x_all  # Observed data vector

# Iteration parameters
max_iter <- 500

# Learning rate parameters
step_size_init <- 1
step_decay_rate <- 0.7
step_decay_interval <- 25

# Initial parameter values
theta_init <- c(median(observed_data), mad(observed_data))

# Other parameters
beta_param <- 0.5
num_draws_y <- 100    # Number of Y samples per iteration
num_sim <- 10000      # Number of repetitions for parallel simulation

##################################################
##  1. Proposed Method (LLB + Monte Carlo approximation)
##################################################

# Start parallel processing
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

tictoc::tic("Proposed Method & Exact LLB - Total Time")

theta_estimates_proposed <- foreach(s = 1:num_sim, 
                                    .combine = rbind, 
                                    .packages = "dqrng") %dopar% {
                                      # Initialize parameters
                                      theta_current <- theta_init
                                      step_size <- step_size_init
                                      
                                      # Generate weighted exponential samples
                                      u <- dqrexp(length(observed_data), 1)
                                      weights <- length(observed_data) * u / sum(u)
                                      
                                      iter_count <- 0
                                      while (iter_count < max_iter) {
                                        # Generate Y samples
                                        y_sample <- generate_y(theta_current, m = num_draws_y)
                                        
                                        # Calculate density and score for both X and Y
                                        dens_x <- calc_density(theta_current, observed_data)
                                        score_x <- calc_score(theta_current, observed_data)
                                        
                                        dens_y <- calc_density(theta_current, y_sample)
                                        score_y <- calc_score(theta_current, y_sample)
                                        
                                        # Calculate stochastic gradient
                                        stoc_grad <- -colMeans(weights * dens_x^beta_param * score_x) +
                                          colMeans(dens_y^beta_param * score_y)
                                        
                                        # Update step size
                                        if (iter_count %% step_decay_interval == 0) {
                                          step_size <- step_size * step_decay_rate
                                        }
                                        
                                        # Update parameters
                                        theta_current <- theta_current - step_size * stoc_grad
                                        iter_count <- iter_count + 1
                                      }
                                      # Return parameter estimates (combined with rbind)
                                      theta_current
                                    }

stopCluster(cl)

##################################################
##  2. Exact LLB
##################################################

# Start parallel processing
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

theta_estimates_llb <- foreach(s = 1:num_sim, 
                               .combine = rbind, 
                               .packages = "dqrng") %dopar% {
                                 theta_current <- theta_init
                                 step_size <- step_size_init
                                 
                                 # Generate weighted exponential samples
                                 u <- dqrexp(length(observed_data), 1)
                                 weights <- length(observed_data) * u / sum(u)
                                 
                                 iter_count <- 0
                                 while (iter_count < max_iter) {
                                   dens_x <- calc_density(theta_current, observed_data)
                                   score_x <- calc_score(theta_current, observed_data)
                                   
                                   # Calculate gradient with regularization for Exact LLB
                                   grad_exact <- -colMeans(weights * dens_x^beta_param * score_x) + 
                                     calc_regularizer(theta_current, beta_param)
                                   
                                   # Update step size
                                   if (iter_count %% step_decay_interval == 0) {
                                     step_size <- step_size * step_decay_rate
                                   }
                                   
                                   # Update parameters
                                   theta_current <- theta_current - step_size * grad_exact
                                   iter_count <- iter_count + 1
                                 }
                                 theta_current
                               }

stopCluster(cl)

tictoc::toc()  # End timing

##############################################
##  Collect and summarize results
##############################################

# For Proposed method
means_proposed <- colMeans(theta_estimates_proposed)
vars_proposed  <- apply(theta_estimates_proposed, 2, var)

# For Exact LLB
means_llb <- colMeans(theta_estimates_llb)
vars_llb  <- apply(theta_estimates_llb, 2, var)

# Display results
cat("\nPosterior Means:\n")
cat("Proposed - mu:", means_proposed[1], "sigma:", means_proposed[2], "\n")
cat("LLB      - mu:", means_llb[1], "sigma:", means_llb[2], "\n")

cat("\nPosterior Variances:\n")
cat("Proposed - mu:", vars_proposed[1], "sigma:", vars_proposed[2], "\n")
cat("LLB      - mu:", vars_llb[1], "sigma:", vars_llb[2], "\n\n")

##############################
##  Visualize Posterior Distributions
##############################

hist_location <- ggplot() +
  stat_bin(mapping = aes(x = theta_estimates_proposed[,1], fill = "Proposed", alpha = 0.5), binwidth = 0.02) +
  stat_bin(mapping = aes(x = theta_estimates_llb[,1], fill = "LLB", alpha = 0.5), binwidth = 0.02) +
  labs(x = expression(mu), y = "Frequency", fill = "Method") +
  guides(alpha = "none") +
  theme_classic(base_size = 14) +
  theme(
    aspect.ratio = 0.6,
    legend.position = c(0.8, 0.85),
    legend.justification = c(0, 1),
    legend.margin = margin(2, 2, 2, 2),
    legend.key.size = unit(1, "lines"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "lines")
  ) +
  coord_cartesian(
    xlim = c(
      min(c(theta_estimates_proposed[,1], theta_estimates_llb[,1])),
      max(c(theta_estimates_proposed[,1], theta_estimates_llb[,1]))
    )
  )

hist_scale <- ggplot() +
  stat_bin(mapping = aes(x = theta_estimates_proposed[,2], fill = "Proposed", alpha = 0.5), binwidth = 0.02) +
  stat_bin(mapping = aes(x = theta_estimates_llb[,2], fill = "LLB", alpha = 0.5), binwidth = 0.02) +
  labs(x = expression(sigma), y = "Frequency", fill = "Method") +
  guides(alpha = "none") +
  theme_classic(base_size = 14) +
  theme(
    aspect.ratio = 0.6,
    legend.position = c(0.8, 0.85),
    legend.justification = c(0, 1),
    legend.margin = margin(2, 2, 2, 2),
    legend.key.size = unit(1, "lines"),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "lines")
  ) +
  coord_cartesian(
    xlim = c(
      min(c(theta_estimates_proposed[,2], theta_estimates_llb[,2])),
      max(c(theta_estimates_proposed[,2], theta_estimates_llb[,2]))
    )
  )

# Display histograms
print(hist_location)
print(hist_scale)