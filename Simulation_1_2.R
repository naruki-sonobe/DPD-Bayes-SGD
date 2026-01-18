###################################
##      Load Required Libraries
###################################

library(cmdstanr)   # For CmdStan-based MCMC
library(doParallel) # For parallel computation in 'foreach'
library(dqrng)      # For random number generation
library(ggplot2)    # For plotting
library(mcmcse)     # For ESS calculation
library(posterior)  # For R-hat calculation
library(tidyr)      # For data manipulation
library(dplyr)      # For data summarization

# Clear the workspace
rm(list = ls())

# If cmdstan is not installed, run the following line once
# check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
# install_cmdstan()

###################################
##         Helper Functions
###################################

generate_data <- function(n = 100, contamination_rate = 0.01) {
  n_main <- ceiling(n * (1 - contamination_rate))
  n_outliers <- n - n_main
  x_main <- dqrnorm(n = n_main, mean = 0, sd = 1)
  x_outliers <- dqrnorm(n = n_outliers, mean = 10, sd = 0.1)
  return(list(
    n_main = n_main, n_outliers = n_outliers, x_main = x_main,
    x_outliers = x_outliers, x_all = c(x_main, x_outliers),
    contamination_rate = contamination_rate
  ))
}
calc_density <- function(theta, x) {
  dnorm(x = x, mean = theta[1], sd = theta[2])
}
calc_score <- function(theta, x) {
  mu <- theta[1]; sigma <- theta[2]
  cbind((x - mu) / (sigma^2), ((x - mu)^2) / (sigma^3) - 1 / sigma)
}
calc_regularizer <- function(theta, alpha) {
  sigma <- theta[2]
  c(0, -(2 * pi)^(-alpha / 2) * alpha * (1 + alpha)^(-3 / 2) * sigma^(-(alpha + 1)))
}
generate_y <- function(theta, m = 10) {
  dqrnorm(n = m, mean = theta[1], sd = theta[2])
}
dpd_function <- function(para, x_data, alpha) {
  mu <- para[1]; sigma <- para[2]
  if (sigma <= 0) return(Inf) 
  ww <- dnorm(x_data, mu, sigma)^alpha
  val <- -ww / alpha + (2 * pi * sigma^2)^(-alpha / 2) * (1 + alpha)^(-3 / 2)
  return(val)
}
log_prior <- function(theta) {
  sigma <- theta[2]
  if (sigma <= 0) return(-Inf) 
  return(-log(sigma))
}
log_likelihood <- function(theta, x, alpha) {
  sum(dpd_function(theta, x, alpha))
}
log_posterior <- function(theta, x, w, alpha) {
  if(theta[2] <= 0) return(-Inf)
  log_prior(theta) - w * log_likelihood(theta, x, alpha)
}


###################################
##       Simulation Parameters
###################################
n <- 1000
contamination_rate <- 0.05
alpha <- 0.5
max_iter <- 500
step_size_init <- 0.1
step_decay_rate <- 0.7
step_decay_interval <- 25
m <- 10
num_sim_llb <- 10000 
thin <- 20
num_total_mh <- thin * num_sim_llb

# Number of repetitions for the entire simulation
N_repetitions <- 20

##################################################
##  Compile Stan Model from External File
##################################################
stan_model_file <- "DPD_GB_Normal.stan"
stan_model <- cmdstan_model(stan_model_file)

##################################################
##  Main Simulation Loop
##################################################
cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

all_results <- foreach(
  rep = 1:N_repetitions,
  .combine = rbind,
  .packages = c("doParallel", "dqrng", "cmdstanr", "posterior", "stats") 
) %dopar% {
  
  # --- 1. Generate Data for this repetition ---
  data_list <- generate_data(n = n, contamination_rate = contamination_rate)
  D <- data_list$x_all
  theta_init <- c(median(D), mad(D))
  
  # --- 2-a. LLB (Stochastic Optimization) ---
  theta_estimates_llb <- foreach(s = 1:num_sim_llb, .combine = rbind) %do% {
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
  
  # --- 2-b. Metropolis-Hastings (MH) ---
  theta_dpd <- theta_init
  eta <- step_size_init
  for (iter_count in 1:max_iter) {
    grad <- -colMeans(calc_density(theta_dpd, D)^alpha * calc_score(theta_dpd, D)) + calc_regularizer(theta_dpd, alpha)
    if (iter_count %% step_decay_interval == 0) eta <- eta * step_decay_rate
    theta_dpd <- theta_dpd - eta * grad
  }
  ep <- 0.01; est <- theta_dpd; R <- 2
  dDPD <- sapply(1:R, function(r) { ee <- rep(0, R); ee[r] <- ep; (dpd_function(est + ee, D, alpha) - dpd_function(est - ee, D, alpha)) / (2 * ep) })
  ddDPD <- array(sapply(1:R, function(r1) sapply(1:R, function(r2) {
    ee1 <- ee2 <- rep(0, R); ee1[r1] <- ep; ee2[r2] <- ep
    (dpd_function(est + ee1 + ee2, D, alpha) + dpd_function(est - ee1 - ee2, D, alpha) - dpd_function(est + ee1 - ee2, D, alpha) - dpd_function(est - ee1 + ee2, D, alpha)) / (4 * ep^2)
  })), dim = c(n, R, R))
  I_mat <- t(dDPD) %*% dDPD / n
  J_mat <- apply(ddDPD, c(2, 3), mean)
  opt_w <- sum(diag(J_mat %*% solve(I_mat) %*% J_mat)) / sum(diag(J_mat))
  
  theta_mh_samples <- matrix(0, num_sim_llb, 2)
  theta_current <- theta_init
  step_size_mh <- 1 / 20
  for (s in 1:num_total_mh) {
    theta_proposal <- theta_current + rnorm(2, 0, step_size_mh)
    if (theta_proposal[2] > 0) {
      log_ratio <- log_posterior(theta_proposal, D, opt_w, alpha) - log_posterior(theta_current, D, opt_w, alpha)
      if (log_ratio > log(dqrunif(1))) theta_current <- theta_proposal
    }
    if (s %% thin == 0) theta_mh_samples[s / thin, ] <- theta_current
  }
  
  # --- 2-c. HMC (CmdStan-based MCMC) ---
  D_stan <- list(
    N = length(data_list$x_all),
    x = data_list$x_all,
    theta0 = theta_init,
    beta = alpha,
    w = opt_w
  )
  stan_fit <- stan_model$sample(data = D_stan, iter_sampling = num_sim_llb, chains = 1, refresh = 0, show_messages = FALSE)
  theta_cmdstan <- as_draws_matrix(stan_fit$draws(variables = c("mu", "sigma")))
  
  # --- 3. Combine and Return Results for this Repetition ---
  data.frame(
    repetition = rep,
    method = rep(c("LLB with SGD", "MH", "HMC"), each = num_sim_llb),
    mu = c(theta_estimates_llb[, 1], theta_mh_samples[, 1], theta_cmdstan[, 1]),
    sigma = c(theta_estimates_llb[, 2], theta_mh_samples[, 2], theta_cmdstan[, 2])
  )
}
stopCluster(cl)

#################################################################
##  Calculate and Display Averages of Posterior Statistics
#################################################################

# 1. 
all_summary_stats <- all_results %>%
  group_by(repetition, method) %>%
  summarise(
    posterior_mean_mu = mean(mu, na.rm = TRUE),
    posterior_var_mu = var(mu, na.rm = TRUE),
    posterior_mean_sigma = mean(sigma, na.rm = TRUE),
    posterior_var_sigma = var(sigma, na.rm = TRUE),
    .groups = 'drop'
  )

# 2. 
final_summary <- all_summary_stats %>%
  group_by(method) %>%
  summarise(
    avg_posterior_mean_mu = mean(posterior_mean_mu, na.rm = TRUE),
    avg_posterior_var_mu = mean(posterior_var_mu, na.rm = TRUE),
    avg_posterior_mean_sigma = mean(posterior_mean_sigma, na.rm = TRUE),
    avg_posterior_var_sigma = mean(posterior_var_sigma, na.rm = TRUE),
    .groups = 'drop'
  )

# results
cat("\n--- Average of Posterior Statistics across", N_repetitions, "Repetitions ---\n")
print(final_summary)

##############################################
##  Prepare Data for Plotting
##############################################
plot_data <- all_results %>%
  pivot_longer(cols = c("mu", "sigma"), names_to = "parameter", values_to = "value")
plot_data$parameter <- factor(plot_data$parameter, levels = c("mu", "sigma"))

##################################################
##  Visualize Posterior Distributions with KDE
##################################################
kde_plot <- ggplot(plot_data, aes(x = value, color = method)) +
  geom_density(aes(group = interaction(method, repetition)), alpha = 0.2, trim = TRUE) +
  facet_wrap(~ parameter, scales = "free", labeller = label_parsed) +
  labs(
    title = element_blank(),
    x = "Parameter Value", y = "Density", color = "Method"
  ) +
  scale_color_brewer(palette = "Set1") + # Use a color-blind friendly palette
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "gray90", color = "gray90"),
    strip.text = element_text(face = "bold", size = 12)
  )

print(kde_plot)
