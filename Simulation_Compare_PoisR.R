# Load necessary libraries
library(doParallel)
library(dqrng)
library(monomvn)   # for mvrnorm
library(mvnfast)
library(tictoc)
rm(list = ls())

# ============ Functions ============

# Function to generate the response variable for Poisson regression.
# It computes lambda = exp(X %*% trueTheta) and then draws y ~ Poisson(lambda).
generateResponse <- function(n, X, trueTheta) {
  lambda <- exp(X %*% trueTheta)
  y <- rpois(n, lambda)
  return(y)
}

# Function to compute the Poisson probability (density) for given theta, response y, and design matrix X.
poissonDensity <- function(theta, y, X) {
  lambda <- exp(X %*% theta)
  return(dpois(x = y, lambda = lambda))
}

# Function to compute the Poisson score (gradient component) for given theta, response y, and design matrix X.
poissonScore <- function(theta, y, X) {
  lambda <- exp(X %*% theta)
  diff <- y - lambda
  return(sweep(X, 1, diff, FUN = "*"))
}

# Function to generate a single pseudo response given theta and design matrix X.
generatePseudoResponse <- function(theta, X) {
  lambda <- exp(X %*% theta)
  return(rpois(n = nrow(X), lambda = lambda))
}

# Function to generate auxiliary responses.
# It returns a matrix of dimension (n x numAuxiliary) where each row has auxiliary draws.
generateAuxiliaryResponses <- function(theta, X, numAuxiliary = 10) {
  lambda <- exp(X %*% theta)
  nObs <- nrow(X)
  Z <- matrix(rpois(nObs * numAuxiliary, lambda), nrow = nObs, ncol = numAuxiliary)
  return(Z)
}

# Function to perform parallel optimization using either Stochastic Gradient Descent (SGD)
# or Full Batch Gradient Descent (FBGD) for the Poisson regression model.
runParallelOptimization <- function(algorithmType) {
  cores <- getOption("mc.cores", detectCores())
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Run replicates in parallel
  results <- tryCatch({
    foreach(replicateIndex = 1:numReplicates,
            .combine = rbind,
            .packages = c("dqrng", "mvnfast"),
            .export = c("n", "p", "X", "y", "initialTheta", "initialLearningRate", "initialStepSize", 
                        "maxIterations", "learningRateDecay", "decayInterval", "alpha", 
                        "numAuxiliarySamples", "gridPoints", "generatePseudoResponse", "generateAuxiliaryResponses",
                        "poissonDensity", "poissonScore")) %dopar% {
                          currentTheta <- initialTheta
                          currentLearningRate <- initialLearningRate
                          stepSize <- initialStepSize
                          
                          # Generate weights using exponential random variables
                          u <- dqrexp(n, 1)
                          weights <- n * u / sum(u)
                          iter <- 0
                          
                          while (iter < maxIterations) {
                            previousTheta <- currentTheta
                            
                            if (algorithmType == "SGD") {
                              # For SGD, generate auxiliary responses using Monte Carlo (m draws)
                              Z <- generateAuxiliaryResponses(currentTheta, X, numAuxiliarySamples)
                            } else {  
                              # For FBGD, create a grid of values (0, 1, ..., gridPoints - 1)
                              gridValues <- 0:(gridPoints - 1)
                              Z <- matrix(rep(gridValues, n), nrow = n, ncol = gridPoints)
                            }
                            
                            # Compute density and score for the observed responses
                            densityY <- poissonDensity(currentTheta, y, X)
                            scoreY <- poissonScore(currentTheta, y, X)
                            
                            # Compute the auxiliary component rTheta
                            rTheta <- matrix(0, n, p + 1)
                            if (algorithmType == "SGD") {
                              for (k in 1:numAuxiliarySamples) {
                                rTheta <- rTheta + (1 / numAuxiliarySamples) * 
                                  (poissonDensity(currentTheta, Z[, k], X)^alpha) * 
                                  poissonScore(currentTheta, Z[, k], X)
                              }
                            } else {  # FBGD
                              for (k in 1:gridPoints) {
                                rTheta <- rTheta + (poissonDensity(currentTheta, Z[, k], X)^(alpha + 1)) * 
                                  poissonScore(currentTheta, Z[, k], X)
                              }
                            }
                            
                            # Update theta using the computed gradient
                            if (algorithmType == "SGD") {
                              stochasticGradient <- - colMeans(weights * (poissonDensity(currentTheta, y, X)^alpha) * 
                                                                 poissonScore(currentTheta, y, X)) +
                                colMeans(weights * rTheta)
                              if (iter %% decayInterval == 0) {
                                currentLearningRate <- currentLearningRate * learningRateDecay
                              }
                              currentTheta <- currentTheta - currentLearningRate * stochasticGradient
                            } else {
                              fullBatchGradient <- - colMeans(weights * (poissonDensity(currentTheta, y, X)^alpha) * 
                                                                poissonScore(currentTheta, y, X)) +
                                colMeans(weights * rTheta)
                              currentTheta <- currentTheta - stepSize * fullBatchGradient
                            }
                            
                            iter <- iter + 1
                          }
                          t(currentTheta)
                        }
  }, error = function(e) {
    stopCluster(cl)
    stop(e)
  })
  
  stopCluster(cl)
  return(results)
}

# ============ Main Simulation ============

# Set the number of observations and the dimensions of the predictor (excluding the intercept)
n <- 300
parameterDimensions <- c(2, 5, 10, 20)
numSimulations <- 100

# Initialize matrices to store results:
# - mseResults: Mean Squared Error of the parameter estimates.
# - avgCIWidth: Average width of the 95% confidence intervals.
# - coverageProb: Coverage probability of the confidence intervals.
mseResults <- matrix(0, nrow = 2, ncol = length(parameterDimensions))
avgCIWidth   <- matrix(0, nrow = 2, ncol = length(parameterDimensions))
coverageProb <- matrix(0, nrow = 2, ncol = length(parameterDimensions))

for (i in seq_along(parameterDimensions)) {
  # Initialize a matrix to count coverage for each coefficient (rows: algorithm, columns: coefficients)
  coverageCount <- matrix(0, nrow = 2, ncol = parameterDimensions[i] + 1)
  
  for (sim in 1:numSimulations) {
    p <- parameterDimensions[i]
    
    # Generate design matrix with intercept
    X <- cbind(1, mvrnorm(n, rep(0, p), diag(p)))
    # Draw true parameter vector from Uniform(0, 1/4) (length: p+1)
    trueTheta <- dqrunif(p + 1, 0, 1/4)
    
    # Generate response variable using the true parameters
    y <- generateResponse(n = n, X = X, trueTheta = trueTheta)
    
    # Set optimization parameters
    numReplicates <- 1000
    maxIterations <- 500
    # Obtain initial parameter estimates using a Poisson GLM (using predictors without the intercept column)
    initialTheta <- as.numeric(glm(y ~ ., data = data.frame(y, X[, -1, drop = FALSE]), 
                                   family = poisson)$coefficients)
    
    # Parameters for SGD
    initialLearningRate <- 1
    numAuxiliarySamples <- n  # number of auxiliary samples for SGD
    learningRateDecay <- 0.7
    decayInterval <- 25
    alpha <- 0.5
    
    # Parameters for FBGD
    initialStepSize <- mean(initialLearningRate * learningRateDecay^(0:((maxIterations / decayInterval) - 1)))
    gridPoints <- n  # number of grid points for FBGD
    
    # Run SGD optimization
    tic("SGD")
    thetaEstimates_SGD <- runParallelOptimization("SGD")
    sgdTime <- toc(quiet = TRUE)$callback_msg  # (Optional) capture timing message
    
    mseResults[1, i] <- mseResults[1, i] + mean((colMeans(thetaEstimates_SGD) - trueTheta)^2)
    
    CI_SGD <- apply(thetaEstimates_SGD, 2, function(x) quantile(x, c(0.025, 0.975)))
    for (j in 1:ncol(CI_SGD)) {
      if (CI_SGD[1, j] <= trueTheta[j] && trueTheta[j] <= CI_SGD[2, j]) {
        coverageCount[1, j] <- coverageCount[1, j] + 1
      }
    }
    avgCIWidth[1, i] <- avgCIWidth[1, i] + mean(CI_SGD[2, ] - CI_SGD[1, ])
    coverageProb[1, i] <- mean(coverageCount[1, ])
    toc(log = FALSE)  # End timing for SGD
    
    # Run FBGD optimization
    tic("FBGD")
    thetaEstimates_FBGD <- runParallelOptimization("FBGD")
    fbgdTime <- toc(quiet = TRUE)$callback_msg  # (Optional) capture timing message
    
    mseResults[2, i] <- mseResults[2, i] + mean((colMeans(thetaEstimates_FBGD) - trueTheta)^2)
    
    CI_FBGD <- apply(thetaEstimates_FBGD, 2, function(x) quantile(x, c(0.025, 0.975)))
    for (j in 1:ncol(CI_FBGD)) {
      if (CI_FBGD[1, j] <= trueTheta[j] && trueTheta[j] <= CI_FBGD[2, j]) {
        coverageCount[2, j] <- coverageCount[2, j] + 1
      }
    }
    avgCIWidth[2, i] <- avgCIWidth[2, i] + mean(CI_FBGD[2, ] - CI_FBGD[1, ])
    coverageProb[2, i] <- mean(coverageCount[2, ])
    toc(log = FALSE)  # End timing for FBGD
  }
  
  # Average the results over the number of simulation replications
  mseResults[, i] <- mseResults[, i] / numSimulations
  avgCIWidth[, i] <- avgCIWidth[, i] / numSimulations
  coverageProb[, i] <- coverageProb[, i] / numSimulations
}

# (Optional) Print final results
print("Mean Squared Errors (MSE):")
print(mseResults)
print("Average Confidence Interval Widths:")
print(avgCIWidth)
print("Coverage Probabilities:")
print(coverageProb)