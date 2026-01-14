#######################################
##         Load Necessary Libraries
#######################################
library(doParallel)
library(tictoc)
library(actuar)
library(ggplot2)
rm(list = ls())

#######################################
##           Helper Functions
#######################################

# Function to generate data (Inverse Gaussian for main data and Normal for outliers)
generate_data <- function(n = 1000, contamination_rate = 0.01) {
  n_main <- ceiling(n * (1 - contamination_rate))
  n_outliers <- n - n_main
  
  x_main <- rinvgauss(n = n_main, mean = 1, shape = 3)
  x_outliers <- rnorm(n = n_outliers, mean = 10, sd = 0.1)
  
  return(list(
    n_main = n_main,
    n_outliers = n_outliers,
    x_main = x_main,
    x_outliers = x_outliers,
    x_all = c(x_main, x_outliers),
    contamination_rate = contamination_rate
  ))
}

# Function to calculate the density of the Inverse Gaussian distribution
calc_density <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  dinvgauss(x = x, mean = mu, shape = sigma)
}

# Function to calculate the score function for the Inverse Gaussian distribution
calc_score <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  cbind(sigma * (x - mu) / (mu^3),
        1 / (2 * sigma) - (x - mu)^2 / (2 * mu^2 * x))
}

# Function to generate Y samples from the Inverse Gaussian distribution
generate_y <- function(theta, m = 10) {
  mu <- theta[1]
  sigma <- theta[2]
  rinvgauss(n = m, mean = mu, shape = sigma)
}

#######################################
##           Main Computation
#######################################

# Set sample size and contamination rate
sample_size <- 1000
contamination_rate <- 0.05

# Generate data
data_list <- generate_data(n = sample_size, contamination_rate = contamination_rate)
observed_data <- data_list$x_all

# Optimization parameters
beta_param <- 0.5
num_draws_y <- 100
max_iter <- 1000
step_size_init <- 1
step_decay_rate <- 0.7
step_decay_interval <- 25

# Initial parameter values (MLE-based initialization)
theta_init <- c(mean(observed_data), 1 / mean(1 / observed_data - 1 / mean(observed_data)))

# Number of simulation iterations
num_sim <- 1000

tic()
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

# LLB (Monte Carlo approximation) parameter estimation
theta_estimates_proposed <- foreach(s = 1:num_sim, .combine = rbind, .packages = "actuar") %dopar% {
  theta_current <- theta_init
  step_size <- step_size_init
  
  # Generate weighted exponential random variables
  u <- rexp(sample_size, 1)
  weights <- sample_size * u / sum(u)
  
  iter_count <- 0
  while (iter_count < max_iter) {
    prev_theta <- theta_current
    Y <- generate_y(theta_current, m = num_draws_y)
    
    density_X <- calc_density(theta_current, observed_data)
    score_X   <- calc_score(theta_current, observed_data)
    
    density_Y <- calc_density(theta_current, Y)
    score_Y   <- calc_score(theta_current, Y)
    
    stoc_grad <- - colMeans(weights * density_X^beta_param * score_X) +
      colMeans(density_Y^beta_param * score_Y)
    
    if (iter_count %% step_decay_interval == 0) {
      step_size <- step_size * step_decay_rate
    }
    
    theta_current <- theta_current - step_size * stoc_grad
    
    iter_count <- iter_count + 1
  }
  theta_current
}

stopCluster(cl)
toc()

# Compute the average estimated parameters
theta_est_proposed <- colMeans(theta_estimates_proposed)
mu_hat    <- theta_est_proposed[1]
sigma_hat <- theta_est_proposed[2]

# Prepare data for visualization
observed_data_df <- data.frame(x = observed_data)
xx <- seq(min(observed_data), max(observed_data), length.out = 200)
density_df <- data.frame(x = xx, 
                         density = dinvgauss(xx, mean = mu_hat, shape = sigma_hat))

# Plot histogram and estimated density
dens_plot <- ggplot(observed_data_df, aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 30, 
                 fill = "lightblue", color = "black", alpha = 0.5) +
  geom_line(data = density_df, aes(x = x, y = density), 
            color = "magenta", linewidth = 1.2) +
  labs(title = element_blank(),
       x = "x", y = "Density") +
  theme_classic(base_size = 14)

print(dens_plot)