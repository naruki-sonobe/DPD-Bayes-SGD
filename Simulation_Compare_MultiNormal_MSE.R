# Load necessary libraries
library(doParallel)
library(tictoc)
library(mvnfast)
library(parallel)  # for detectCores()

# Clear workspace
rm(list = ls())

# ============ Functions ============

# Function to generate synthetic data with contamination
generateData <- function(n = 100, p = 2, contaminationRate = 0.01) {
  nNonContaminated <- ceiling(n * (1 - contaminationRate))
  nContaminated <- n - nNonContaminated
  
  nonContaminatedData <- rmvn(n = nNonContaminated, mu = rep(0, p), sigma = diag(p))
  contaminatedData    <- rmvn(n = nContaminated,    mu = rep(10, p), sigma = 0.01 * diag(p))
  
  combinedData <- rbind(nonContaminatedData, contaminatedData)
  
  return(list(nNonContaminated = nNonContaminated,
              nContaminated = nContaminated,
              nonContaminatedData = nonContaminatedData,
              contaminatedData = contaminatedData,
              dataMatrix = combinedData,
              contaminationRate = contaminationRate))
}

# Function to compute the multivariate normal density at points in X
mvnDensity <- function(theta, X) {
  p <- ncol(X)
  return(dmvn(X = X, mu = theta, sigma = diag(p)))
}

# Function to center the rows of matrix X by subtracting theta from each row
centerData <- function(theta, X) {
  return(sweep(X, 2, theta, FUN = "-"))
}

# Function to generate auxiliary samples from a multivariate normal distribution
generateAuxiliarySamples <- function(theta, m = 10) {
  p <- length(theta)
  return(rmvn(n = m, mu = theta, sigma = diag(p)))
}

# Parallel optimization function that runs either SGD or Full Batch Gradient Descent (FBGD)
runParallelOptimization <- function(algorithmType) {
  cores <- getOption("mc.cores", detectCores())
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # The foreach loop will run numReplicates independent replicates in parallel.
  results <- tryCatch({
    foreach(rep = 1:numReplicates, 
            .combine = rbind, 
            .packages = "mvnfast",
            .export = c("n", "p", "initialTheta", "initialLearningRate", "initialStepSize",
                        "dataMatrix", "maxIterations", "learningRateDecay", "decayInterval",
                        "beta", "batchSize", "gridBound", "gridPoints",
                        "generateAuxiliarySamples", "mvnDensity", "centerData")) %dopar% {
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
                              auxiliarySamples <- generateAuxiliarySamples(currentTheta, m = batchSize)
                            } else {  # For FBGD
                              gridSequence <- seq(-gridBound, gridBound, length.out = gridPoints)
                              # Create a grid of points for all dimensions
                              auxiliarySamples <- as.matrix(expand.grid(rep(list(gridSequence), p)))
                            }
                            
                            densityX <- mvnDensity(currentTheta, dataMatrix)
                            centeredX <- centerData(currentTheta, dataMatrix)
                            
                            densityY <- mvnDensity(currentTheta, auxiliarySamples)
                            centeredY <- centerData(currentTheta, auxiliarySamples)
                            
                            if (algorithmType == "SGD") {
                              # Compute stochastic gradient
                              stochasticGradient <- -colMeans(weights * (densityX^beta) * centeredX) +
                                colMeans((densityY^beta) * centeredY)
                              # Decay the learning rate at specified intervals
                              if (iter %% decayInterval == 0) {
                                currentLearningRate <- currentLearningRate * learningRateDecay
                              }
                              currentTheta <- currentTheta - currentLearningRate * stochasticGradient
                            } else {  # FBGD
                              fullBatchGradient <- -colMeans(weights * centeredX * (densityX^beta)) +
                                colMeans(((2 * gridBound)^p) * centeredY * (densityY^(beta + 1)))
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

# ============ Main Computation ============

# Set basic parameters
n <- 100
dimensions <- c(2, 3)
contaminationRate <- 0.05
numSimulations <- 1000  # Number of simulation repetitions

# Initialize matrices to store results
mseResults <- matrix(0, nrow = length(dimensions), ncol = 2)
varianceResults <- matrix(0, nrow = length(dimensions), ncol = 2) #<-- ADDED: Matrix for variance

# Loop over each sample dimension
for (i in seq_along(dimensions)) {
  for (sim in 1:numSimulations) {
    p <- dimensions[i]
    
    # Generate data for the current simulation and dimension
    dataSet <- generateData(n = n, p = p, contaminationRate = contaminationRate)
    dataMatrix <- dataSet$dataMatrix
    
    # Set parameters for parallel replicates and optimization iterations
    numReplicates <- 1000     # Number of replicates for the parallel optimization
    maxIterations <- 500      # Maximum number of iterations per replicate
    initialTheta <- apply(dataMatrix, 2, median)
    
    # Parameters for SGD
    initialLearningRate <- 1
    learningRateDecay <- 0.7
    decayInterval <- 25
    beta <- 0.5
    batchSize <- 10
    
    # Parameters for FBGD
    initialStepSize <- mean(initialLearningRate * learningRateDecay^(0:((maxIterations / decayInterval) - 1)))
    gridPoints <- 10
    gridBound <- 2
    
    # Run SGD optimization and time it
    tic("SGD")
    thetaResult_SGD <- runParallelOptimization("SGD")
    sgdTime <- toc(quiet = TRUE)$callback_msg  # (Optional) time message
    mseResults[i, 1] <- mseResults[i, 1] + mean(colMeans(thetaResult_SGD)^2)
    # ADDED: Calculate and store the average variance for SGD parameters
    varianceResults[i, 1] <- varianceResults[i, 1] + mean(apply(thetaResult_SGD, 2, var))
    
    # Run FBGD optimization and time it
    tic("FBGD")
    thetaResult_FBGD <- runParallelOptimization("FBGD")
    fbgdTime <- toc(quiet = TRUE)$callback_msg  # (Optional) time message
    mseResults[i, 2] <- mseResults[i, 2] + mean(colMeans(thetaResult_FBGD)^2)
    # ADDED: Calculate and store the average variance for FBGD parameters
    varianceResults[i, 2] <- varianceResults[i, 2] + mean(apply(thetaResult_FBGD, 2, var))
  }
  # Average the results over the number of simulations
  mseResults[i, ] <- mseResults[i, ] / numSimulations
  varianceResults[i, ] <- varianceResults[i, ] / numSimulations 
}

# Set column names for clarity in the output
colnames(mseResults) <- c("SGD", "FBGD")
colnames(varianceResults) <- c("SGD", "FBGD")
rownames(mseResults) <- paste("Dim=", dimensions)
rownames(varianceResults) <- paste("Dim=", dimensions)

# Print final results
print("Mean Squared Errors (MSE):")
print(mseResults)

print("Average Posterior Variances:")
print(varianceResults)