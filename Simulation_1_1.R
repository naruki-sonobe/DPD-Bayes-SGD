#######################################
##   Load Necessary Libraries
#######################################
# install.packages("doParallel")
# install.packages("foreach")
# install.packages("ggplot2")
# install.packages("tictoc")
# install.packages("dqrng")
# install.packages("tidyr")
# install.packages("dplyr") # For data summarization

library(doParallel)  # For parallel computing
library(foreach)     # For 'foreach' loops
library(ggplot2)     # For plotting
library(tictoc)      # For timing
library(dqrng)       # For generating random numbers (dqrexp, dqrnorm)
library(tidyr)       # For data manipulation (pivot_longer)
library(dplyr)       # For data summarization

# Clear the workspace
rm(list = ls())

#######################################
##         Helper Functions
#######################################

# Function to generate data with contamination
generate_data <- function(n = 100, contamination_rate = 0.01) {
  n_main <- ceiling(n * (1 - contamination_rate))
  n_outliers <- n - n_main
  
  x_main <- dqrnorm(n = n_main, mean = 0, sd = 1)
  x_outliers <- dqrnorm(n = n_outliers, mean = 10, sd = 0.1)
  
  return(c(x_main, x_outliers))
}

# Calculate the density of N(mu, sigma^2)
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
##       Simulation Parameters
#######################################

# General parameters
sample_size <- 1000
contamination_rate <- 0.05
max_iter <- 500
beta_param <- 0.5
num_draws_y <- 10
num_sim <- 10000 # Posterior samples per repetition

# Number of repetitions for the entire simulation
N_repetitions <- 20

# Learning rate parameters
step_size_init <- 0.1
step_decay_rate <- 0.7
step_decay_interval <- 25

##################################################
##  Main Simulation Loop
##  (Run the entire experiment N_repetitions times)
##################################################

# Start parallel processing for the outer loop
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

tictoc::tic("Total Simulation Time")

# Use foreach to parallelize the N_repetitions
all_results <- foreach(
  rep = 1:N_repetitions, 
  .combine = rbind, 
  .packages = c("dqrng", "doParallel")
) %dopar% {
  
  # 1. Generate a new dataset for this repetition
  observed_data <- generate_data(n = sample_size, contamination_rate = contamination_rate)
  theta_init <- c(median(observed_data), mad(observed_data))
  
  # --- Proposed Method ---
  theta_estimates_proposed <- foreach(s = 1:num_sim, .combine = rbind) %do% {
    theta_current <- theta_init
    step_size <- step_size_init
    u <- dqrexp(length(observed_data), 1)
    weights <- length(observed_data) * u / sum(u)
    
    for (iter_count in 1:max_iter) {
      y_sample <- generate_y(theta_current, m = num_draws_y)
      dens_x <- calc_density(theta_current, observed_data)
      score_x <- calc_score(theta_current, observed_data)
      dens_y <- calc_density(theta_current, y_sample)
      score_y <- calc_score(theta_current, y_sample)
      stoc_grad <- -colMeans(weights * dens_x^beta_param * score_x) + colMeans(dens_y^beta_param * score_y)
      if (iter_count %% step_decay_interval == 0) step_size <- step_size * step_decay_rate
      theta_current <- theta_current - step_size * stoc_grad
    }
    theta_current
  }
  
  # --- Exact LLB Method ---
  theta_estimates_llb <- foreach(s = 1:num_sim, .combine = rbind) %do% {
    theta_current <- theta_init
    step_size <- step_size_init
    u <- dqrexp(length(observed_data), 1)
    weights <- length(observed_data) * u / sum(u)
    
    for (iter_count in 1:max_iter) {
      dens_x <- calc_density(theta_current, observed_data)
      score_x <- calc_score(theta_current, observed_data)
      grad_exact <- -colMeans(weights * dens_x^beta_param * score_x) + calc_regularizer(theta_current, beta_param)
      if (iter_count %% step_decay_interval == 0) step_size <- step_size * step_decay_rate
      theta_current <- theta_current - step_size * grad_exact
    }
    theta_current
  }
  
  # 2. Combine results into a data frame for this repetition
  data.frame(
    repetition = rep,
    method = rep(c("LLB with SGD", "LLB"), each = num_sim),
    mu = c(theta_estimates_proposed[, 1], theta_estimates_llb[, 1]),
    sigma = c(theta_estimates_proposed[, 2], theta_estimates_llb[, 2])
  )
}

stopCluster(cl)
tictoc::toc()

#################################################################
##  Calculate and Display Averages of Posterior Statistics
#################################################################

all_summary_stats <- all_results %>%
  group_by(repetition, method) %>%
  summarise(
    posterior_mean_mu = mean(mu, na.rm = TRUE),
    posterior_var_mu = var(mu, na.rm = TRUE),
    posterior_mean_sigma = mean(sigma, na.rm = TRUE),
    posterior_var_sigma = var(sigma, na.rm = TRUE),
    .groups = 'drop'
  )

final_summary <- all_summary_stats %>%
  group_by(method) %>%
  summarise(
    avg_posterior_mean_mu = mean(posterior_mean_mu, na.rm = TRUE),
    avg_posterior_var_mu = mean(posterior_var_mu, na.rm = TRUE),
    avg_posterior_mean_sigma = mean(posterior_mean_sigma, na.rm = TRUE),
    avg_posterior_var_sigma = mean(posterior_var_sigma, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n--- Average of Posterior Statistics across", N_repetitions, "Repetitions ---\n")
print(final_summary)


##############################################
##  Prepare Data for Plotting
##############################################

# Convert the results from wide to long format for easier plotting with ggplot2
plot_data <- all_results %>%
  pivot_longer(
    cols = c("mu", "sigma"),
    names_to = "parameter",
    values_to = "value"
  )

# Ensure 'parameter' is a factor to control facet order
plot_data$parameter <- factor(plot_data$parameter, levels = c("mu", "sigma"))

##################################################
##  Visualize Posterior Distributions with KDE
##################################################

kde_plot <- ggplot(plot_data, aes(x = value, color = method)) +
  # Add the geom_density layer. 
  # The 'group' aesthetic is crucial: it tells ggplot to draw a separate line
  # for each combination of method and repetition.
  # 'alpha' is set to a low value to make the lines transparent, showing overlap.
  geom_density(aes(group = interaction(method, repetition)), alpha = 0.2, trim = TRUE) +
  
  # Use facet_wrap to create separate panels for mu and sigma.
  # 'scales = "free"' allows the x and y axes to adapt to the data in each panel.
  # 'labeller' is used to parse Greek letters for facet titles.
  facet_wrap(~ parameter, scales = "free", labeller = label_parsed) +
  
  # Add labels and a title
  labs(
    title = element_blank(),
    x = "Parameter Value",
    y = "Density",
    color = "Method" 
  ) +
  
  # Apply a clean theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "gray90", color = "gray90"),
    strip.text = element_text(face = "bold", size = 12)
  )

# Display the plot
print(kde_plot)