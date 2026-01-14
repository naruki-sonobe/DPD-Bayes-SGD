###################################
##      Load Required Libraries
###################################

library(doParallel) # For parallel computation in 'foreach'
library(dqrng)      # For random number generation
library(ggplot2)    # For plotting
library(tidyr)      # For data manipulation
library(dplyr)      # For data summarization
library(mcmcse)     # For calculating Effective Sample Size (ESS)
library(tictoc)     # For timing the script

# Clear the workspace
rm(list = ls())

###################################
##          Helper Functions
###################################

# --- Functions for Data Generation and DPD LLB ---
generate_data <- function(n = 100, contamination_rate = 0.01) {
  n_main <- ceiling(n * (1 - contamination_rate))
  n_outliers <- n - n_main
  x_main <- dqrnorm(n = n_main, mean = 0, sd = 1)
  x_outliers <- dqrnorm(n = n_outliers, mean = 10, sd = 0.1)
  return(list(
    x_all = c(x_main, x_outliers)
  ))
}
calc_density <- function(theta, x) {
  dnorm(x = x, mean = theta[1], sd = theta[2])
}
calc_score <- function(theta, x) {
  mu <- theta[1]; sigma <- theta[2]
  cbind((x - mu) / (sigma^2), ((x - mu)^2) / (sigma^3) - 1 / sigma)
}
generate_y <- function(theta, m = 10) {
  dqrnorm(n = m, mean = theta[1], sd = theta[2])
}

# --- Functions for Standard Bayes MH Algorithm ---

# Log-prior for Standard Bayes: p(mu, sigma) ∝ 1/sigma
log_prior_standard <- function(theta) {
  sigma <- theta[2]
  if (sigma <= 0) {
    return(-Inf) # Return -Inf if sigma is non-positive
  }
  return(-log(sigma))
}

# Log-likelihood for Standard Bayes: Normal distribution
log_likelihood_standard <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  if (sigma <= 0) {
    return(-Inf)
  }
  sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

# Log-posterior for Standard Bayes
log_posterior_standard <- function(theta, x) {
  log_prior_standard(theta) + log_likelihood_standard(theta, x)
}


###################################
##       Simulation Parameters
###################################
n <- 1000
contamination_rate <- 0.05
alpha <- 0.5
max_iter <- 500
step_size_init <- 1
step_decay_rate <- 0.7
step_decay_interval <- 25
m <- 100
num_sim <- 1000 # Posterior samples

# Number of repetitions for the entire simulation
N_repetitions <- 20 


##################################################
##            Main Simulation Loop
##################################################
tictoc::tic()

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

# The foreach loop now returns a list of lists, so .combine is removed
loop_results <- foreach(
  rep = 1:N_repetitions,
  .packages = c("doParallel", "dqrng", "stats", "mcmcse") # Added mcmcse
) %dopar% {
  
  # --- 1. Generate Data for this repetition ---
  data_list <- generate_data(n = n, contamination_rate = contamination_rate)
  D <- data_list$x_all
  theta_init <- c(median(D), mad(D))
  
  # --- 2-a. DPD LLB Posterior (Stochastic Optimization) ---
  theta_estimates_llb <- foreach(s = 1:num_sim, .combine = rbind) %do% {
    theta_current <- c(theta_init[1], log(theta_init[2]))
    step_size <- step_size_init
    u <- dqrexp(length(D), 1)
    weights <- length(D) * u / sum(u)
    for (iter_count in 1:max_iter) {
      y_sample <- generate_y(c(theta_current[1], exp(theta_current[2])), m = m)
      dens_x <- calc_density(c(theta_current[1], exp(theta_current[2])), D)
      score_x <- calc_score(c(theta_current[1], exp(theta_current[2])), D)
      dens_y <- calc_density(c(theta_current[1], exp(theta_current[2])), y_sample)
      score_y <- calc_score(c(theta_current[1], exp(theta_current[2])), y_sample)
      stoc_grad <- -colMeans(weights * dens_x^alpha * score_x) + colMeans(dens_y^alpha * score_y)
      if (iter_count %% step_decay_interval == 0) step_size <- step_size * step_decay_rate
      theta_current <- theta_current - step_size * stoc_grad * c(1, exp(theta_current[2]))
    }
    c(theta_current[1], exp(theta_current[2]))
  }
  
  # --- 2-b. Metropolis-Hastings (Standard Bayes) ---
  thin <- 50
  num_total_mh <- num_sim * thin
  theta_mh_samples <- matrix(0, nrow = num_sim, ncol = 2)
  colnames(theta_mh_samples) <- c("mu", "sigma") # Name columns for clarity
  theta_current <- theta_init
  
  # Standard deviation for the proposal distribution (may need tuning)
  step_size_mh <- c(1/20, 1/20) 
  
  current_log_post <- log_posterior_standard(theta_current, D)
  
  for (s in 1:num_total_mh) {
    # Propose new parameters
    theta_proposal <- theta_current + rnorm(2, mean = 0, sd = step_size_mh)
    
    # Check if proposal is in the valid parameter space (sigma > 0)
    if (theta_proposal[2] > 0) {
      proposal_log_post <- log_posterior_standard(theta_proposal, D)
      
      # Calculate acceptance ratio
      log_ratio <- proposal_log_post - current_log_post
      
      if (log(dqrunif(1)) < log_ratio) {
        theta_current <- theta_proposal
        current_log_post <- proposal_log_post
      }
    }
    # If theta_proposal[2] <= 0, the proposal is automatically rejected
    
    # Store the sample after thinning
    if (s %% thin == 0) {
      theta_mh_samples[s / thin, ] <- theta_current
    }
  }
  
  # --- 3. Calculate ESS for the MH samples ---
  ess_results <- c(
    mu_ess = mcmcse::ess(theta_mh_samples[, "mu"]),
    sigma_ess = mcmcse::ess(theta_mh_samples[, "sigma"])
  )
  
  # --- 4. Combine and Return Results for this Repetition ---
  # Combine samples into a data.frame
  samples_df <- data.frame(
    repetition = rep,
    method = rep(c("DPD LLB Posterior", "Standard Bayes (MH)"), each = num_sim),
    mu = c(theta_estimates_llb[, 1], theta_mh_samples[, "mu"]),
    sigma = c(theta_estimates_llb[, 2], theta_mh_samples[, "sigma"])
  )
  
  # Return both samples and ESS as a list
  return(list(samples = samples_df, ess = ess_results))
}
stopCluster(cl)

#################################################################
##           Process and Summarize Simulation Results
#################################################################

# Extract samples and ESS results from the loop output
all_results_df <- dplyr::bind_rows(lapply(loop_results, function(x) x$samples))
all_ess_df <- dplyr::bind_rows(lapply(loop_results, function(x) as.data.frame(t(x$ess))))

# --- Display Average ESS ---
avg_ess <- colMeans(all_ess_df, na.rm = TRUE)
cat("\n--- Average Effective Sample Size (ESS) for MH method across", N_repetitions, "Repetitions ---\n")
print(avg_ess)

# --- Calculate and Display Averages of Posterior Statistics ---
all_summary_stats <- all_results_df %>%
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
##      Prepare Data for Plotting
##############################################

plot_data <- all_results_df %>%
  pivot_longer(cols = c("mu", "sigma"), names_to = "parameter", values_to = "value")
plot_data$parameter <- factor(plot_data$parameter, levels = c("mu", "sigma"))

##################################################
##  Visualize Posterior Distributions with KDE
##################################################

kde_plot <- ggplot(plot_data, aes(x = value, color = method)) +
  geom_density(aes(group = interaction(method, repetition)), alpha = 0.2, trim = TRUE) +
  facet_wrap(~ parameter, scales = "free", labeller = label_parsed) +
  labs(
    title = "Comparison of Posterior Distributions",
    x = "Parameter Value",
    y = "Density",
    color = "Method"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "gray90", color = "gray90"),
    strip.text = element_text(face = "bold", size = 12)
  )

print(kde_plot)
tictoc::toc()