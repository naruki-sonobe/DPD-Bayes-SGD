# Load necessary libraries
library(doParallel)
library(tictoc)
library(mvnfast)
library(parallel)  # For detectCores()

# Clear the workspace
rm(list = ls())

# ============ Functions ============

# Function to generate synthetic data with contamination
generateData <- function(n = 100, p = 2, contaminationRate = 0.01) {
  nNonContaminated <- ceiling(n * (1 - contaminationRate))
  nContaminated <- n - nNonContaminated
  
  nonContaminatedData <- rmvn(n = nNonContaminated, mu = rep(0, p), sigma = diag(p))
  contaminatedData    <- rmvn(n = nContaminated,    mu = rep(10, p), sigma = 0.01 * diag(p))
  
  dataMatrix <- rbind(nonContaminatedData, contaminatedData)
  
  return(list(nNonContaminated = nNonContaminated,
              nContaminated = nContaminated,
              nonContaminatedData = nonContaminatedData,
              contaminatedData = contaminatedData,
              dataMatrix = dataMatrix,
              contaminationRate = contaminationRate))
}

# Function to compute the multivariate normal density at points in X
mvnormDensity <- function(theta, X) {
  p <- ncol(X)
  return(dmvn(X = X, mu = theta, sigma = diag(p)))
}

# Function to center the rows of matrix X by subtracting theta
centerData <- function(theta, X) {
  return(sweep(X, 2, theta, FUN = "-"))
}

# Function to generate auxiliary samples from a multivariate normal distribution
generateY <- function(theta, m = 10) {
  p <- length(theta)
  return(rmvn(n = m, mu = theta, sigma = diag(p)))
}

# Function to run parallel optimization (either SGD or Full Batch Gradient Descent)
runParallelOptimization <- function(algorithmType) {
  cores <- getOption("mc.cores", detectCores())
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  results <- tryCatch({
    foreach(rep = 1:numReplicates,
            .combine = rbind,
            .packages = "mvnfast",
            .export = c("n", "p", "initialTheta", "initialLearningRate",
                        "initialStepSize", "dataMatrix", "maxIterations",
                        "learningRateDecay", "decayInterval", "beta", "batchSize",
                        "gridBound", "gridPoints", "generateY", "mvnormDensity", "centerData")
    ) %dopar% {
      currentTheta <- initialTheta
      currentLearningRate <- initialLearningRate
      stepSize <- initialStepSize
      
      # Generate weights using exponential random variables
      u <- rexp(n, 1)
      weights <- n * u / sum(u)
      iter <- 0
      
      # Optimization loop
      while (iter < maxIterations) {
        if (algorithmType == "SGD") {
          # Generate a batch of auxiliary samples for SGD
          Y <- generateY(currentTheta, m = batchSize)
        } else {
          # Generate grid points for Full Batch Gradient Descent
          gridSeq <- seq(-gridBound, gridBound, length.out = gridPoints)
          Y <- as.matrix(expand.grid(rep(list(gridSeq), p)))
        }
        
        densityX <- mvnormDensity(currentTheta, dataMatrix)
        centeredX <- centerData(currentTheta, dataMatrix)
        
        densityY <- mvnormDensity(currentTheta, Y)
        centeredY <- centerData(currentTheta, Y)
        
        if (algorithmType == "SGD") {
          # Compute stochastic gradient
          stochasticGradient <- -colMeans(weights * densityX^beta * centeredX) +
            colMeans(densityY^beta * centeredY)
          # Decay the learning rate at specified intervals
          if (iter %% decayInterval == 0) {
            currentLearningRate <- currentLearningRate * learningRateDecay
          }
          currentTheta <- currentTheta - currentLearningRate * stochasticGradient
        } else {
          # Compute full batch gradient
          gradient <- -colMeans(weights * centeredX * densityX^beta) +
            colMeans((2 * gridBound)^p * centeredY * densityY^(beta + 1))
          currentTheta <- currentTheta - stepSize * gradient
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

# ============ Main Computation ============

# Set sample size and contamination rate
n <- 100
dimensions <- c(2, 3)
contaminationRate <- 0.05

# Matrices to store execution times and mean squared errors
executionTimes <- matrix(0, nrow = length(dimensions), ncol = 2)
meanSquaredErrors <- matrix(0, nrow = length(dimensions), ncol = 2)

# Loop over different sample dimensions
for(i in seq_along(dimensions)) {
  p <- dimensions[i]
  
  # Generate data for the current dimension
  dataList <- generateData(n = n, p = p, contaminationRate = contaminationRate)
  dataMatrix <- dataList$dataMatrix
  
  # Set simulation parameters
  numReplicates <- 1000
  maxIterations <- 500
  initialTheta <- apply(dataMatrix, 2, median)
  
  # Parameters for SGD
  initialLearningRate <- 1
  learningRateDecay <- 0.7
  decayInterval <- 25
  beta <- 0.5
  batchSize <- 10
  
  # Parameters for Full Batch Gradient Descent (FBGD)
  initialStepSize <- mean(initialLearningRate * learningRateDecay^(0:((maxIterations / decayInterval) - 1)))
  gridPoints <- 10
  gridBound <- 2
  
  # Run SGD optimization in parallel and time it
  tic("SGD")
  thetaResults_SGD <- runParallelOptimization("SGD")
  sgdTime <- toc(quiet = TRUE)$callback_msg
  meanSquaredErrors[i, 1] <- mean(colMeans(thetaResults_SGD)^2)
  
  # Run FBGD optimization in parallel and time it
  tic("FBGD")
  thetaResults_FBGD <- runParallelOptimization("FBGD")
  fbgdTime <- toc(quiet = TRUE)$callback_msg
  meanSquaredErrors[i, 2] <- mean(colMeans(thetaResults_FBGD)^2)
  
  # Record execution time messages
  executionTimes[i, ] <- c(sgdTime, fbgdTime)
}

# (Optional) Print the results
print("Mean Squared Errors (by column):")
print(meanSquaredErrors)
print("Execution Times:")
print(executionTimes)