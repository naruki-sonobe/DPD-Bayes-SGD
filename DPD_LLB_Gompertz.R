#######################################
##         Load Necessary Libraries
#######################################
library(doParallel)
library(tictoc)
library(VGAM)
library(ggplot2)
rm(list = ls())

#######################################
##          Helper Functions
#######################################

# Function to generate data using the Gompertz distribution for main data
# and Normal distribution for outliers
generate_data <- function(n = 1000, contamination_rate = 0.01) {
  n_main <- ceiling(n * (1 - contamination_rate))
  n_outliers <- n - n_main
  
  x_main <- rgompertz(n = n_main, scale = 1, shape = 0.1)
  x_outliers <- rnorm(n = n_outliers, mean = 10, sd = 0.1)
  
  return(list(
    n_main = n_main,
    n_outliers = n_outliers,
    x_main = x_main,
    x_outliers = x_outliers,
    observed_data = append(x_main, x_outliers),
    contamination_rate = contamination_rate
  ))
}

# Function to calculate the density of the Gompertz distribution
calc_density <- function(theta, x) {
  scale <- theta[1]
  shape <- theta[2]
  dgompertz(x = x, scale = scale, shape = shape)
}

# Function to calculate the score function for the Gompertz distribution
calc_score <- function(theta, x) {
  scale <- theta[1]
  shape <- theta[2]
  coef <- (1 - exp(scale * x)) / (scale^2) + x * exp(scale * x) / scale
  cbind(x - shape * coef, 1 / shape + (1 - exp(scale * x)) / scale)
}

# Function to generate Y samples from the Gompertz distribution
generate_y <- function(theta, m = 10) {
  scale <- theta[1]
  shape <- theta[2]
  rgompertz(n = m, scale = scale, shape = shape)
}

# Function to compute the MLE for the Gompertz distribution parameters
estimate_mle <- function(x) {
  f <- function(w) {
    weights <- 1 / mean(1 - exp(w * x))
    v1 <- (1 - exp(w * x)) / w
    v2 <- x * exp(w * x)
    sum(x) + weights * (sum(v1) + sum(v2))
  }
  scale_est <- uniroot(f, c(0.01, 10))$root
  shape_est <- -scale_est / mean(1 - exp(scale_est * x))
  return(c(scale_est, shape_est))
}

#######################################
##          Main Computation
#######################################

# Set sample size and contamination rate
sample_size <- 1000
contamination_rate <- 0.05

# Generate data
data_list <- generate_data(n = sample_size, contamination_rate = contamination_rate)
observed_data <- data_list$observed_data

# Optimization parameters
beta_param <- 0.5
num_draws_y <- 100
max_iter <- 1000
step_size_init <- 0.5
step_decay_rate <- 0.7
step_decay_interval <- 25

# Initial parameter values (obtained via MLE)
theta_init <- estimate_mle(observed_data)

# Number of simulation iterations
num_sim <- 1000

tictoc::tic()
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

# LLB (Monte Carlo approximation) parameter estimation
# The optimization is performed in the log-parameter space.
THETA <- foreach(s = 1:num_sim, .packages = "VGAM", .combine = "rbind") %dopar% {
  # 'theta' is in the log-scale so that exp(theta) > 0
  theta <- theta_init
  step_size <- step_size_init
  
  # Generate weighted exponential random variables
  u <- rexp(sample_size, 1)
  weights <- sample_size * u / sum(u)
  
  iter_count <- 0
  while (iter_count < max_iter) {
    prev_theta <- theta
    # Evaluate functions at the exponentiated parameter values
    Y <- generate_y(exp(theta), m = num_draws_y)
    
    density_X <- calc_density(exp(theta), observed_data)
    score_X   <- calc_score(exp(theta), observed_data)
    
    density_Y <- calc_density(exp(theta), Y)
    score_Y   <- calc_score(exp(theta), Y)
    
    # Compute stochastic gradient
    stoc_grad <- - colMeans(weights * density_X^beta_param * score_X) +
      colMeans(density_Y^beta_param * score_Y)
    
    if (iter_count %% step_decay_interval == 0) {
      step_size <- step_size * step_decay_rate
    }
    
    # Update in log-scale (hence the multiplication by exp(theta))
    theta <- theta - step_size * stoc_grad * exp(theta)
    
    iter_count <- iter_count + 1
  }
  # Return the estimated parameters in the original scale
  t(exp(theta))
}
stopCluster(cl)
tictoc::toc()

# Compute the average estimated parameters over simulations
theta_est <- colMeans(THETA)
scale_est <- theta_est[1]
shape_est <- theta_est[2]

# Prepare data for visualization
observed_data_df <- data.frame(x = observed_data)
xx <- seq(min(observed_data), max(observed_data), length.out = 200)
density_df <- data.frame(x = xx, density = dgompertz(xx, scale = scale_est, shape = shape_est))

# Plot histogram and estimated density curve
dens_plot <- ggplot(observed_data_df, aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 30,
                 fill = "lightblue", color = "black", alpha = 0.6) +
  geom_line(data = density_df, aes(x = x, y = density),
            color = "magenta", linewidth = 1.5) +
  labs(title = element_blank(),
       x = "x", y = "Density") +
  theme_classic(base_size = 14)

print(dens_plot)