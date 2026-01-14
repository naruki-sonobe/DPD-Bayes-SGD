###################################
##  Libraries
###################################
library(dqrng)      # RNG
library(ggplot2)    # Plot
library(mvtnorm)    # dmvnorm
library(dplyr)      # Summaries
library(tidyr)      # Data reshape
library(gridExtra)  # Arrange plots

# Clear the workspace
rm(list = ls())

###################################
##  DPD for 2D Normal (mean only)
##  - Target: theta = (mu1, mu2)
##  - Known covariance: Sigma = I2
###################################

# Generate 2D data (optionally with 5% outliers)
generate_data_2d <- function(n = 800, contamination_rate = 0.05,
                             mu_true = c(0, 0),
                             mu_out  = c(6, 6)) {
  n_main <- ceiling(n * (1 - contamination_rate))
  n_out  <- n - n_main
  # Inliers ~ N(mu_true, I2)
  xin1 <- dqrnorm(n_main, mean = mu_true[1], sd = 1)
  xin2 <- dqrnorm(n_main, mean = mu_true[2], sd = 1)
  Xin  <- cbind(xin1, xin2)
  # Outliers ~ N(mu_out, 0.1^2 I2)
  if (n_out > 0) {
    xout1 <- dqrnorm(n_out, mean = mu_out[1], sd = 0.1)
    xout2 <- dqrnorm(n_out, mean = mu_out[2], sd = 0.1)
    Xout  <- cbind(xout1, xout2)
    Xall  <- rbind(Xin, Xout)
  } else {
    Xall  <- Xin
  }
  list(X = Xall, n = n, n_main = n_main, n_out = n_out,
       contamination_rate = contamination_rate)
}

# Multivariate normal density with Sigma = I2
dens2 <- function(x, mu) {
  # x: N x 2, mu: length-2
  mvtnorm::dmvnorm(x, mean = mu, sigma = diag(2))
}

# Score with respect to mu when Sigma = I:  s(x; mu) = x - mu
score_mu <- function(x, mu) {
  sweep(x, 2, mu, "-")  # N x 2 matrix
}

# Generate y ~ N(mu, I2), m samples
gen_y <- function(mu, m = 20) {
  cbind(dqrnorm(m, mean = mu[1], sd = 1),
        dqrnorm(m, mean = mu[2], sd = 1))
}

# DPD contribution per observation (vector length N)
# For dimension d=2 and Sigma = I:
#   ∫ f^{1+α} = (2π)^(-α d/2) |I|^(-α/2) (1+α)^(-d/2)
#             = (2π)^(-α) (1+α)^(-1)
dpd_vec <- function(mu, X, alpha) {
  # negative log pseudo-lik piece per obs (up to constant), matching univariate style
  wx  <- dens2(X, mu)^alpha
  cst <- (2 * pi)^(-alpha) * (1 + alpha)^(-2)
  -wx / alpha + cst
}

# Log-prior for mu (weakly informative normal, independent)
log_prior_mu <- function(mu, sd0 = 10) {
  sum(dnorm(mu, mean = 0, sd = sd0, log = TRUE))
}

# Log-likelihood-like sum for DPD (sum over obs of dpd piece)
log_like_dpd <- function(mu, X, alpha) {
  sum(dpd_vec(mu, X, alpha))
}

# DPD-based posterior (power-posterior style)
log_post_mu <- function(mu, X, w, alpha, sd0 = 10) {
  log_prior_mu(mu, sd0 = sd0) - w * log_like_dpd(mu, X, alpha)
}

###################################
##  LLB with SGD (2D mean)
##   - Follows your univariate structure
###################################
llb_sgd_sampler_mu <- function(X, alpha = 0.5,
                               num_samples = 4000,
                               max_iter = 400,
                               step_size_init = 0.15,
                               step_decay_rate = 0.7,
                               step_decay_interval = 25,
                               m = 20) {
  n <- nrow(X)
  # random Dirichlet-like weights via normalized Exp(1)
  llb_draw <- function(mu_init) {
    mu_cur <- mu_init
    step   <- step_size_init
    # one set of Exponential weights per sample (like in your code)
    u  <- dqrng::dqrexp(n, 1)
    wN <- n * u / sum(u)
    for (t in 1:max_iter) {
      Y      <- gen_y(mu_cur, m = m)
      dens_x <- dens2(X,  mu_cur)
      sc_x   <- score_mu(X, mu_cur)             # N x 2
      dens_y <- dens2(Y,  mu_cur)
      sc_y   <- score_mu(Y, mu_cur)             # m x 2
      # stochastic gradient (2-vector)
      # - E_{weighted data}[ f(x)^α s(x) ] + E_y[ f(y)^α s(y) ]
      g1 <- colMeans(wN * (dens_x^alpha) * sc_x)
      g2 <- colMeans((dens_y^alpha) * sc_y)
      stoc_grad <- -g1 + g2
      # decay
      if (t %% step_decay_interval == 0) step <- step * step_decay_rate
      mu_cur <- mu_cur - step * stoc_grad
    }
    mu_cur
  }
  
  # initialize at robust center
  mu0 <- apply(as.matrix(X), 2, median) 
  if (any(is.na(mu0))) mu0 <- colMeans(X)
  
  # draw num_samples iid "LLB estimates"
  out <- t(replicate(num_samples, llb_draw(mu0)))
  colnames(out) <- c("mu1", "mu2")
  out
}

###################################
##  MH sampler (2D mean)
###################################
mh_sampler_mu <- function(X, alpha = 0.5,
                          w = 1.0,
                          num_draws = 4000,
                          thin = 50,
                          step_sd = 0.10,
                          burn = 1000,
                          prior_sd = 10) {
  total_iter <- burn + num_draws * thin
  mu_cur <- colMeans(X)  # start near data center
  keep   <- matrix(NA_real_, nrow = num_draws, ncol = 2)
  k <- 0
  for (s in 1:total_iter) {
    mu_prop <- mu_cur + rnorm(2, 0, step_sd)
    loga <- log_post_mu(mu_prop, X, w, alpha, sd0 = prior_sd) -
      log_post_mu(mu_cur,  X, w, alpha, sd0 = prior_sd)
    if (loga > log(runif(1))) mu_cur <- mu_prop
    if (s > burn && (s - burn) %% thin == 0) {
      k <- k + 1
      keep[k, ] <- mu_cur
    }
  }
  colnames(keep) <- c("mu1", "mu2")
  keep
}

###################################
##  Weight w (same spirit as your code)
##  - Numerical I and J with finite diff on DPD vector
###################################
opt_weight_mu <- function(mu_hat, X, alpha, eps = 1e-2) {
  n  <- nrow(X)
  R  <- 2
  # First derivatives (n x R)
  dDPD <- sapply(1:R, function(r) {
    e <- rep(0, R); e[r] <- eps
    (dpd_vec(mu_hat + e, X, alpha) - dpd_vec(mu_hat - e, X, alpha)) / (2 * eps)
  })
  # Second derivatives per obs (n x R x R)
  ddDPD <- array(NA_real_, dim = c(n, R, R))
  for (r1 in 1:R) for (r2 in 1:R) {
    e1 <- rep(0, R); e2 <- rep(0, R)
    e1[r1] <- eps; e2[r2] <- eps
    dd <- dpd_vec(mu_hat + e1 + e2, X, alpha) +
      dpd_vec(mu_hat - e1 - e2, X, alpha) -
      dpd_vec(mu_hat + e1 - e2, X, alpha) -
      dpd_vec(mu_hat - e1 + e2, X, alpha)
    ddDPD[, r1, r2] <- dd / (4 * eps^2)
  }
  I_mat <- t(dDPD) %*% dDPD / n
  J_mat <- apply(ddDPD, c(2, 3), mean)  # average Hessian
  # same formula as your univariate code
  as.numeric(sum(diag(J_mat %*% solve(I_mat) %*% J_mat)) / sum(diag(J_mat)))
}

###################################
##  Credible ellipse (95%)
###################################
ellipse_points <- function(mu_hat, Sigma_hat, level = 0.95, npt = 200) {
  # returns data.frame with columns x,y
  csq <- qchisq(level, df = 2)
  eig <- eigen(Sigma_hat, symmetric = TRUE)
  A   <- eig$vectors %*% diag(sqrt(eig$values * csq)) %*% t(eig$vectors)
  t   <- seq(0, 2 * pi, length.out = npt)
  circle <- cbind(cos(t), sin(t))
  E <- sweep(circle %*% t(A), 2, mu_hat, "+")
  data.frame(x = E[,1], y = E[,2])
}

###################################
##  One experiment: return samples & ellipse
###################################
run_once <- function(contamination_rate = 0.00,
                     n = 800,
                     alpha = 0.5,
                     num_llb = 4000,
                     num_mh  = 4000,
                     seed = 1234) {
  dqset.seed(seed)
  dat <- generate_data_2d(n = n, contamination_rate = contamination_rate)
  X   <- dat$X
  
  # Robust pilot for w: do a short LLB to get mu_hat
  llb_pilot <- llb_sgd_sampler_mu(X, alpha = alpha,
                                  num_samples = 200,
                                  max_iter = 200,
                                  step_size_init = 0.2,
                                  m = 10)
  mu_hat <- colMeans(llb_pilot)
  
  # Compute w (same spirit as your code)
  wopt <- opt_weight_mu(mu_hat, X, alpha)
  
  # Full LLB samples
  llb_samps <- llb_sgd_sampler_mu(X, alpha = alpha,
                                  num_samples = num_llb,
                                  max_iter = 400,
                                  step_size_init = 0.15,
                                  m = 20)
  
  # MH samples (DPD posterior with weight wopt)
  mh_samps <- mh_sampler_mu(X, alpha = alpha, w = wopt,
                            num_draws = num_mh, thin = 50,
                            step_sd = 0.12, burn = 1500)
  
  # Summaries & 95% credible ellipses
  llb_mean <- colMeans(llb_samps)
  mh_mean  <- colMeans(mh_samps)
  llb_cov  <- stats::cov(llb_samps)
  mh_cov   <- stats::cov(mh_samps)
  
  ell_llb <- ellipse_points(llb_mean, llb_cov, level = 0.95)
  ell_mh  <- ellipse_points(mh_mean,  mh_cov,  level = 0.95)
  
  list(
    X = X, alpha = alpha, wopt = wopt,
    llb = list(samples = llb_samps, mean = llb_mean, cov = llb_cov, ellipse = ell_llb),
    mh  = list(samples = mh_samps,  mean = mh_mean,  cov = mh_cov,  ellipse = ell_mh)
  )
}

###################################
##  Run: clean vs contaminated
###################################
alpha <- 0.5
set.seed(1)

res_clean <- run_once(contamination_rate = 0.00, n = 800, alpha = alpha,
                      num_llb = 3000, num_mh = 3000, seed = 1)
res_cont  <- run_once(contamination_rate = 0.05, n = 800, alpha = alpha,
                      num_llb = 3000, num_mh = 3000, seed = 2)

###################################
##  Credible ellipse (95%)  
###################################
ellipse_points <- function(mu_hat, Sigma_hat, level = 0.95, npt = 200) {
  csq <- qchisq(level, df = 2)
  eig <- eigen(Sigma_hat, symmetric = TRUE)
  A   <- eig$vectors %*% diag(sqrt(eig$values * csq)) %*% t(eig$vectors)
  t   <- seq(0, 2 * pi, length.out = npt)
  circle <- cbind(cos(t), sin(t))
  E <- sweep(circle %*% t(A), 2, mu_hat, "+")
  E <- rbind(E, E[1, , drop = FALSE])  
  data.frame(x = E[,1], y = E[,2])
}

###################################
##  Multi-level ellipses (50% & 95%)
###################################
build_ellipses <- function(res, levels = c(0.50, 0.95)) {
  # LLB
  llb_mean <- res$llb$mean; llb_cov <- res$llb$cov
  ell_llb <- do.call(rbind, lapply(levels, function(lv) {
    e <- ellipse_points(llb_mean, llb_cov, level = lv)
    cbind(e, method = "LLB with SGD", level = paste0(round(lv*100), "%"))
  }))
  # MH
  mh_mean <- res$mh$mean; mh_cov <- res$mh$cov
  ell_mh <- do.call(rbind, lapply(levels, function(lv) {
    e <- ellipse_points(mh_mean, mh_cov, level = lv)
    cbind(e, method = "MH", level = paste0(round(lv*100), "%"))
  }))
  rbind(ell_llb, ell_mh)
}

###################################
##  Plot helper
###################################
plot_scene_overlay <- function(res, title_txt = "Clean", show_points = FALSE,
                               levels = c(0.50, 0.95)) {
  X <- as.data.frame(res$X); colnames(X) <- c("x","y")
  mu_llb <- as.data.frame(t(res$llb$mean)); colnames(mu_llb) <- c("x","y"); mu_llb$method <- "LLB with SGD"
  mu_mh  <- as.data.frame(t(res$mh$mean));  colnames(mu_mh)  <- c("x","y"); mu_mh$method  <- "MH"
  mu_df  <- rbind(mu_llb, mu_mh)
  
  ell_df <- build_ellipses(res, levels)
  ell_df$level <- factor(ell_df$level, levels = paste0(sort(levels)*100, "%"))
  
  p <- ggplot() +
    { if (show_points) geom_point(data = X, aes(x = x, y = y), alpha = 0.12, size = 0.8) } +
    
    ## 95%
    geom_polygon(
      data = subset(ell_df, level == paste0(round(100*max(levels)), "%")),
      aes(x = x, y = y, fill = method, color = method, linetype = level),
      alpha = 0.15, linewidth = 0.7
    ) +
    
    ## 50%
    geom_polygon(
      data = subset(ell_df, level == paste0(round(100*min(levels)), "%")),
      aes(x = x, y = y, fill = method, color = method, linetype = level),
      alpha = 0.30, linewidth = 0.7
    ) +
    
    geom_point(data = mu_df, aes(x = x, y = y, shape = method), size = 3) +
    geom_point(aes(x = 0, y = 0), shape = 4, size = 3, stroke = 1.1) +
    coord_equal() +
    labs(
      title = title_txt,
      x = expression(mu[1]), y = expression(mu[2]),
      fill = "Method", color = "Method", shape = "Method", linetype = "Credible level"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")
  
  # p <- p + scale_fill_manual(values = c("LLB with SGD"="#1f77b4","MH"="#d62728")) +
  #          scale_color_manual(values = c("LLB with SGD"="#1f77b4","MH"="#d62728"))
  
  p
}

p1 <- plot_scene_overlay(res_clean, title_txt = "No outliers (0%)",
                         show_points = FALSE, levels = c(0.50, 0.95))
p2 <- plot_scene_overlay(res_cont,  title_txt = "With outliers (5%)",
                         show_points = FALSE, levels = c(0.50, 0.95))

leg <- cowplot::get_legend(
  p1 + theme(
    legend.position = "top",
    legend.box.spacing = unit(2, "pt"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )
)
p1_nl <- p1 + theme(legend.position = "none",
                    plot.margin = margin(t = 2, r = 6, b = 6, l = 6))
p2_nl <- p2 + theme(legend.position = "none",
                    plot.margin = margin(t = 2, r = 6, b = 6, l = 6))

cowplot::plot_grid(
  leg,
  cowplot::plot_grid(p1_nl, p2_nl, ncol = 2, align = "hv"),
  ncol = 1,
  rel_heights = c(0.08, 1)  
)
