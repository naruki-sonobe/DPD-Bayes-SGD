###################################
##       Load Required Libraries
###################################
# If needed, uncomment and install:
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("doParallel")
# install.packages("dqrng")
# install.packages("ggplot2")
# install.packages("mcmcse")
# install.packages("posterior")

library(cmdstanr)   # For CmdStan-based MCMC
library(doParallel) # For parallel computation in 'foreach'
library(dqrng)      # For random number generation (dqrexp, dqrnorm, dqrunif)
library(ggplot2)    # For plotting histograms
library(mcmcse)     # For ESS calculation (ess())
library(posterior)  # For R-hat calculation (rhat())

# Clear the workspace
rm(list = ls())

###################################
##         Helper Functions
###################################

# 1. Data Generation with Contamination
generate_data <- function(n = 100, contamination_rate = 0.01) {
  n_main <- ceiling(n * (1 - contamination_rate))
  n_outliers <- n - n_main
  
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

# 2. Calculate the density of N(mu, sigma^2) for a given theta = (mu, sigma)
calc_density <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  dnorm(x = x, mean = mu, sd = sigma)
}

# 3. Calculate the score function
calc_score <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  cbind(
    (x - mu) / (sigma^2),
    ((x - mu)^2) / (sigma^3) - 1 / sigma
  )
}

# 4. Calculate the regularization term (used in LLB)
calc_regularizer <- function(theta, alpha) {
  # Only the sigma part is used in this regularizer
  sigma <- theta[2]
  c(
    0,
    -(2 * pi)^(-alpha / 2) * alpha * (1 + alpha)^(-3 / 2) * sigma^(-(alpha + 1))
  )
}

# 5. Generate Y ~ N(mu, sigma^2)
generate_y <- function(theta, m = 10) {
  mu <- theta[1]
  sigma <- theta[2]
  dqrnorm(n = m, mean = mu, sd = sigma)
}

# 6. Define DPD function
dpd_function <- function(para, alpha) {
  mu <- para[1]
  sigma <- para[2]
  ww <- dnorm(D, mu, sigma)^alpha
  val <- -ww / alpha + (2 * pi * sigma^2)^(-alpha / 2) * (1 + alpha)^(-3 / 2)
  return(val)
}

# 7. Define log-prior (improper prior ~ 1/sigma)
log_prior <- function(theta) {
  sigma <- theta[2]
  return(-log(sigma))
}

# 8. Define log-likelihood (based on DPD)
log_likelihood <- function(theta, x) {
  sum(dpd_function(theta, alpha))
}

# 9. Define log-posterior
log_posterior <- function(theta, x, w) {
  # Posterior ~ prior * (likelihood^w)
  log_prior(theta) - w * log_likelihood(theta, x)
}

###################################
##        Main Computations
###################################

## 1. Generate Data
n <- 1000
contamination_rate <- 0.05
data_list <- generate_data(n = n, contamination_rate = contamination_rate)
D <- data_list$x_all  # Observed data

## 2. LLB settings
max_iter <- 500         # Number of iterations
step_size_init <- 1
step_decay_rate <- 0.7
step_decay_interval <- 25

theta_init <- c(median(D), mad(D))  # Initial guess for (mu, sigma)
alpha <- 0.5
m <- 100       # Sample size of Y in each iteration
num_sim_llb <- 10000  # Repetitions for LLB

##--------------------------------------
## 2-a. LLB (Stochastic Optimization)
##--------------------------------------

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

# LLB uses log(sigma) parameterization internally
# so we transform: (mu, log(sigma)) and update in that space
theta_estimates_llb <- foreach(s = 1:num_sim_llb, 
                               .combine = rbind, 
                               .packages = "dqrng") %dopar% {
                                 theta_current <- c(theta_init[1], log(theta_init[2])) 
                                 step_size <- step_size_init
                                 
                                 # Generate weighted exponential samples
                                 u <- dqrexp(length(D), 1)
                                 weights <- length(D) * u / sum(u)
                                 
                                 iter_count <- 0
                                 while (iter_count < max_iter) {
                                   # Generate Y samples with the current (mu, exp(log_sigma))
                                   y_sample <- generate_y(c(theta_current[1], exp(theta_current[2])), m = m)
                                   
                                   # Calculate density and score for X
                                   dens_x <- calc_density(c(theta_current[1], exp(theta_current[2])), D)
                                   score_x <- calc_score(c(theta_current[1], exp(theta_current[2])), D)
                                   
                                   # Calculate density and score for Y
                                   dens_y <- calc_density(c(theta_current[1], exp(theta_current[2])), y_sample)
                                   score_y <- calc_score(c(theta_current[1], exp(theta_current[2])), y_sample)
                                   
                                   # Stochastic gradient
                                   stoc_grad <- -colMeans(weights * dens_x^alpha * score_x) +
                                     colMeans(dens_y^alpha * score_y)
                                   
                                   # Decay the step size
                                   if (iter_count %% step_decay_interval == 0) {
                                     step_size <- step_size * step_decay_rate
                                   }
                                   
                                   # Update (mu, log_sigma)
                                   # The gradient w.r.t log_sigma is multiplied by exp(log_sigma)
                                   theta_current <- theta_current - step_size * stoc_grad * c(1, exp(theta_current[2]))
                                   iter_count <- iter_count + 1
                                 }
                                 
                                 # Convert back from (mu, log_sigma) to (mu, sigma)
                                 c(theta_current[1], exp(theta_current[2]))
                               }

stopCluster(cl)

##--------------------------------------
## 2-b. Metropolis-Hastings (MH)
##--------------------------------------

theta_mh_samples <- matrix(0, num_sim_llb, 2)
theta_current <- theta_init
step_size_mh <- 1 / 20
accept_count <- 0
thin <- 20

# First, obtain DPD estimator (for initialization)
theta_dpd <- theta_init
max_iter_dpd <- 500
eta <- step_size_init

iter_count <- 0
while (iter_count < max_iter_dpd) {
  dens_x <- calc_density(theta_dpd, D)
  score_x <- calc_score(theta_dpd, D)
  
  # Gradient for DPD (includes small regularizer calc_regularizer)
  grad <- -colMeans(dens_x^alpha * score_x) + calc_regularizer(theta_dpd, alpha)
  
  # Decay the step size
  if (iter_count %% step_decay_interval == 0) {
    eta <- eta * step_decay_rate
  }
  
  theta_dpd <- theta_dpd - eta * grad
  iter_count <- iter_count + 1
}

# Estimate w (opt_w) for MH using empirical approximation of the Fisher matrix
ep <- 0.01
est <- theta_dpd
R <- 2

# First derivative approximation
dDPD <- matrix(NA, n, R)
for (r in 1:R) {
  ee <- rep(0, R)
  ee[r] <- ep
  dDPD[, r] <- (dpd_function(est + ee, alpha) - dpd_function(est - ee, alpha)) / (2 * ep)
}

# Second derivative approximation
ddDPD <- array(NA, c(n, R, R))
for (r1 in 1:R) {
  for (r2 in 1:R) {
    ee1 <- ee2 <- rep(0, R)
    ee1[r1] <- ep
    ee2[r2] <- ep
    ddDPD[, r1, r2] <- (
      dpd_function(est + ee1 + ee2, alpha) +
        dpd_function(est - ee1 - ee2, alpha) -
        dpd_function(est + ee1 - ee2, alpha) -
        dpd_function(est - ee1 + ee2, alpha)
    ) / (4 * ep^2)
  }
}

# Calculate I & J matrices
I_mat <- 0
for (j in 1:n) {
  I_mat <- I_mat + dDPD[j, ] %*% t(dDPD[j, ]) / n
}
J_mat <- apply(ddDPD, c(2, 3), mean)

opt_w <- sum(diag(J_mat %*% solve(I_mat) %*% J_mat)) / sum(diag(J_mat))

# Run MH
num_total_mh <- thin * num_sim_llb
for (s in 1:num_total_mh) {
  # Propose new theta
  theta_proposal <- theta_current + rnorm(2, 0, step_size_mh)
  
  # Calculate log acceptance ratio
  log_ratio <- log_posterior(theta_proposal, D, opt_w) - 
    log_posterior(theta_current, D, opt_w)
  
  # Accept/reject
  if (log_ratio > log(dqrunif(1))) {
    theta_current <- theta_proposal
    accept_count <- accept_count + 1
  }
  
  # Thinning
  if (s %% thin == 0) {
    idx <- s / thin
    theta_mh_samples[idx, ] <- theta_current
  }
}

# MH convergence diagnosis
accept_ratio <- accept_count / (thin * num_sim_llb)
ess_mh <- ess(theta_mh_samples)  # from mcmcse
rhat_mh <- apply(theta_mh_samples, 2, rhat) # from posterior

##--------------------------------------
## 2-c. CmdStan-based MCMC
##--------------------------------------

# Compile Stan model (adjust file path as needed)
stan_model_file <- "/absolute/path/to/DPD_GB_Normal.stan"
stan_model <- cmdstan_model(stan_model_file)

# Data for Stan
D_stan <- list(
  N = length(data_list$x_all),
  x = data_list$x_all,
  theta0 = theta_init,
  beta = alpha,
  w = opt_w
)

# Sampling
stan_fit <- fit_model$sample(
  data = D_stan,
  iter_sampling = num_sim_llb,
  thin = 4
)

# Extract samples
theta_cmdstan <- posterior::as_draws_matrix(stan_fit$draws(variables = c("mu", "sigma")))

# Convergence diagnostics (CmdStan)
ess_stan <- ess(theta_cmdstan)  # from mcmcse or posterior
rhat_stan <- apply(theta_cmdstan, 2, rhat)

###################################
##       Summarize Results
###################################

# For Proposed (LLB) - (theta_estimates_llb)
llb_means <- colMeans(theta_estimates_llb)
llb_vars  <- apply(theta_estimates_llb, 2, var)

# For MH - (theta_mh_samples)
mh_means <- colMeans(theta_mh_samples)
mh_vars  <- apply(theta_mh_samples, 2, var)

# For Stan - (theta_cmdstan)
stan_means <- colMeans(theta_cmdstan)
stan_vars  <- apply(theta_cmdstan, 2, var)

cat("\nPosterior Means:\n")
cat("LLB-Proposed - mu:", llb_means[1], " sigma:", llb_means[2], "\n")
cat("MH           - mu:", mh_means[1],  " sigma:", mh_means[2],  "\n")
cat("Stan         - mu:", stan_means[1]," sigma:", stan_means[2],"\n")

cat("\nPosterior Variances:\n")
cat("LLB-Proposed - mu:", llb_vars[1],  " sigma:", llb_vars[2],  "\n")
cat("MH           - mu:", mh_vars[1],   " sigma:", mh_vars[2],   "\n")
cat("Stan         - mu:", stan_vars[1], " sigma:", stan_vars[2], "\n\n")

###################################
##       Plot Histograms
###################################

# Histogram of the data
hist_data <- ggplot() +
  stat_bin(aes(x = D, fill = "Data", alpha = 0.5), bins = 30) +
  guides(fill = "none", alpha = "none") +
  labs(x = "Value", y = "Frequency") +
  theme_classic()

# Location parameter
hist_mu <- ggplot() +
  stat_bin(aes(x = theta_estimates_llb[, 1], fill = "Proposed", alpha = 0.5), binwidth = 0.02) +
  stat_bin(aes(x = theta_mh_samples[, 1], fill = "MH", alpha = 0.5), binwidth = 0.02) +
  stat_bin(aes(x = theta_cmdstan[, 1], fill = "Stan", alpha = 0.5), binwidth = 0.02) +
  labs(x = expression(mu), y = "Frequency", fill = "Method") +
  guides(alpha = "none") +
  coord_cartesian(
    xlim = c(
      min(c(theta_estimates_llb[,1], theta_mh_samples[,1], theta_cmdstan[,1])),
      max(c(theta_estimates_llb[,1], theta_mh_samples[,1], theta_cmdstan[,1]))
    )
  ) +
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
  )

# Scale parameter
hist_sigma <- ggplot() +
  stat_bin(aes(x = theta_estimates_llb[, 2], fill = "Proposed", alpha = 0.5), binwidth = 0.02) +
  stat_bin(aes(x = theta_mh_samples[, 2], fill = "MH", alpha = 0.5), binwidth = 0.02) +
  stat_bin(aes(x = theta_cmdstan[, 2], fill = "Stan", alpha = 0.5), binwidth = 0.02) +
  labs(x = expression(sigma), y = "Frequency", fill = "Method") +
  guides(alpha = "none") +
  coord_cartesian(
    xlim = c(
      min(c(theta_estimates_llb[,2], theta_mh_samples[,2], theta_cmdstan[,2])),
      max(c(theta_estimates_llb[,2], theta_mh_samples[,2], theta_cmdstan[,2]))
    )
  ) +
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
  )

# Print plots
print(hist_data)
print(hist_mu)
print(hist_sigma)