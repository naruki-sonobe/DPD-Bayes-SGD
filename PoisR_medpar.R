# ============ 1. Load Libraries and Data ============
library(COUNT)
library(robustbase)
library(doParallel)
library(foreach)
library(dqrng)
library(mvnfast)
library(tictoc)
library(MCMCpack) # Standard Bayesian estimation

rm(list = ls())

# Load data
data("medpar", package = "COUNT")

# Response and design matrix: los ~ hmo + white + age80 + died + type
y <- medpar$los
X <- model.matrix(~ hmo + white + age80 + died + type, data = medpar)

# ============ 2. Required Functions ============
# Clip eta for numerical stability (avoid exp overflow/underflow)

poissonDensity <- function(theta, y, X) {
  eta <- as.vector(X %*% theta)
  eta <- pmin(pmax(eta, -700), 700)
  lambda <- exp(eta)
  lambda[lambda <= 0] <- 1e-12
  lambda[is.infinite(lambda)] <- .Machine$double.xmax
  dpois(x = y, lambda = lambda)
}

poissonScore <- function(theta, y, X) {
  eta <- as.vector(X %*% theta)
  eta <- pmin(pmax(eta, -700), 700)
  lambda <- exp(eta)
  lambda[lambda <= 0] <- 1e-12
  lambda[is.infinite(lambda)] <- .Machine$double.xmax
  diff <- y - lambda
  sweep(X, 1, diff, FUN = "*")
}

generateAuxiliaryResponses <- function(theta, X, numAuxiliary = 10) {
  eta <- as.vector(X %*% theta)
  eta <- pmin(pmax(eta, -700), 700)
  lambda <- exp(eta)
  lambda[lambda <= 0] <- 1e-12
  lambda[is.infinite(lambda)] <- .Machine$double.xmax
  nObs <- nrow(X)
  matrix(rpois(nObs * numAuxiliary, lambda), nrow = nObs, ncol = numAuxiliary)
}

runParallelOptimization <- function() {
  cores <- getOption("mc.cores", detectCores())
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  results <- tryCatch({
    foreach(replicateIndex = 1:numReplicates,
            .combine = rbind,
            .packages = c("dqrng", "mvnfast"),
            .export = c("n", "p", "X", "y", "initialTheta", "initialLearningRate",
                        "maxIterations", "learningRateDecay", "decayInterval", "alpha",
                        "numAuxiliarySamples", "generateAuxiliaryResponses",
                        "poissonDensity", "poissonScore")) %dopar% {
                          
                          currentTheta <- initialTheta
                          currentLearningRate <- initialLearningRate
                          
                          # Dirichlet(1,...,1)-equivalent weights via exponential draws
                          u <- dqrexp(n, 1)
                          weights <- n * u / sum(u)
                          
                          iter <- 0
                          while (iter < maxIterations) {
                            densityY <- poissonDensity(currentTheta, y, X)  # length n
                            scoreY   <- poissonScore(currentTheta, y, X)    # n x p
                            
                            # Auxiliary samples from current model
                            Z <- generateAuxiliaryResponses(currentTheta, X, numAuxiliarySamples)
                            
                            # Monte Carlo estimate of auxiliary term
                            rTheta <- matrix(0, n, p)
                            for (k in 1:numAuxiliarySamples) {
                              densZk  <- poissonDensity(currentTheta, Z[, k], X)
                              scoreZk <- poissonScore(currentTheta, Z[, k], X)
                              rTheta  <- rTheta + (densZk^alpha) * scoreZk / numAuxiliarySamples
                            }
                            
                            # Stochastic gradient (DPD objective)
                            grad_obs <- colMeans(weights * (densityY^alpha) * scoreY)
                            grad_aux <- colMeans(weights * rTheta)
                            stochasticGradient <- -grad_obs + grad_aux
                            
                            # Step-size decay
                            if (iter > 0 && (iter %% decayInterval == 0)) {
                              currentLearningRate <- currentLearningRate * learningRateDecay
                            }
                            
                            currentTheta <- as.numeric(currentTheta - currentLearningRate * stochasticGradient)
                            iter <- iter + 1
                          }
                          
                          t(currentTheta)
                        }
  }, error = function(e) {
    stopCluster(cl); stop(e)
  })
  
  stopCluster(cl)
  results
}

# ============ 3. Set Parameters ============
n <- nrow(X)
p <- ncol(X)

numReplicates <- 1000  # bootstrap reps
maxIterations <- 500   # SGD iterations

# Initialize with standard Poisson GLM
glm.fit <- glm(los ~ hmo + white + age80 + died + type,
               data = medpar, family = poisson)
initialTheta <- coef(glm.fit)

# SGD/DPD hyperparameters
initialLearningRate <- 0.01
numAuxiliarySamples <- 100
learningRateDecay   <- 0.7
decayInterval       <- 25
alpha               <- 0.5  # DPD robustness level

# ============ 4. Run Estimations (SGD and MCMC) ============
tic("Robust Method with SGD (LLB + SGD)")
thetaEstimates_SGD <- runParallelOptimization()
toc()

tic("Standard Bayesian MCMC (Poisson)")
mcmc.fit <- MCMCpoisson(
  los ~ hmo + white + age80 + died + type,
  data    = medpar,
  mcmc    = 10000,
  thin    = 10,
  burnin  = 10000,
  verbose = 0
)
toc()

# ============ 5. Summaries and Comparison ============

# LLB + SGD summary
results_LLB_SGD <- data.frame(
  Mean     = colMeans(thetaEstimates_SGD),
  Variance = apply(thetaEstimates_SGD, 2, var),
  CI_Lower = apply(thetaEstimates_SGD, 2, quantile, probs = 0.025),
  CI_Upper = apply(thetaEstimates_SGD, 2, quantile, probs = 0.975)
)
results_LLB_SGD$CI_Width <- results_LLB_SGD$CI_Upper - results_LLB_SGD$CI_Lower
rownames(results_LLB_SGD) <- colnames(X)

# MCMC summary
mcmc_summary <- summary(mcmc.fit)
results_MCMC <- data.frame(
  Mean     = mcmc_summary$statistics[, "Mean"],
  Variance = (mcmc_summary$statistics[, "SD"])^2,
  CI_Lower = mcmc_summary$quantiles[, "2.5%"],
  CI_Upper = mcmc_summary$quantiles[, "97.5%"]
)
results_MCMC$CI_Width <- results_MCMC$CI_Upper - results_MCMC$CI_Lower

cat("\n--- LLB with SGD Results (COUNT::medpar, y = los) ---\n")
print(results_LLB_SGD)

cat("\n--- Standard Bayesian MCMC Results (COUNT::medpar, y = los) ---\n")
print(results_MCMC)

# Focus on the bulk: exclude top 5% largest los
idx_mid <- which(y <= quantile(y, 0.95))
Xmid <- X[idx_mid, , drop = FALSE]
ymid <- y[idx_mid]

# Stable log(mean(exp(.))) utility
logmeanexp <- function(a) { m <- max(a); m + log(mean(exp(a - m))) }

# Bulk LPD for MCMC
B_mcmc <- as.matrix(mcmc.fit)
lpd_mcmc_mid <- sapply(seq_len(nrow(B_mcmc)), function(s) {
  mu <- as.vector(exp(Xmid %*% B_mcmc[s, ]))
  sum(dpois(ymid, mu, log = TRUE))
})
LPD_mcmc_mid <- logmeanexp(lpd_mcmc_mid)

# Bulk LPD for LLB
B_llb <- thetaEstimates_SGD
lpd_llb_mid <- sapply(seq_len(nrow(B_llb)), function(s) {
  mu <- as.vector(exp(Xmid %*% B_llb[s, ]))
  sum(dpois(ymid, mu, log = TRUE))
})
LPD_llb_mid <- logmeanexp(lpd_llb_mid)

cat("\nBulk LPD (<=95% quantile of los):\n")
print(c(LPD_llb_mid = LPD_llb_mid,
        LPD_mcmc_mid = LPD_mcmc_mid,
        diff = LPD_llb_mid - LPD_mcmc_mid))
