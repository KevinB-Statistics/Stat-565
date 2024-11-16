############################################################################
# Load required library
library(MASS)
library(mclust)
library(knitr)
library(kableExtra)
############################################################################
set.seed(123)
############################################################################
# SIMULATING DATA FOR EM
############################################################################
# sim.GMMdata <- function(mu, sigma, pi) {
# FUNCTION DETAILS.
# Input:
#   List of 4.
# *Element 1 of the list.
#   n: scalar, sample size.
# *Element 2 of the list.
#   List of 3.
#     mu1: 1 x 2 vector, mean of first component of GMM.
#     mu2: 1 x 2 vector, mean of second component of GMM.
#     mu3: 1 x 2 vector, mean of third component of GMM.
# * Element 3 of the list.
#   List of 3.
#     Sigma1: 1 x 3 vector, covariance matrix of the first component of GMM.
#       Elements ordered as: var(x1), var(x2), cov(x1,x2).
#       *Note that this is not a correlation matrix.
#     Sigma2: 1 x 3 vector, covariance matrix of the second component of GMM.
#       Elements ordered as: same as for Sigma1
#     Sigma3: 1 x 3 vector, covariance matrix of the third component of GMM.
#       Elements ordered as: same as for Sigma1
# * Element 4 of the list.
#   Vector of 1 x 3.
#     pi = (pi1,pi2,pi3)
#     pi1: scalar, mixing probability for the first component of GMM.
#     pi2: scalar, mixing probability for the second component of GMM.
#     pi3: scalar, mixing probability for the third component of GMM.
# Returns:
#   Vector of 1 x n simulated data.
# # }
############################################################################
# Found on stackexchange to visualize the data
density.hist <- function(df, x=NULL, y=NULL) {
  
  require(ggplot2)
  require(gridExtra)
  require(devtools)
  
  htop <- ggplot(data = df, aes_string(x = x)) + 
    geom_histogram(aes(y=..density..), fill = "white", color = "black", bins=100) + 
    stat_density(colour = "blue", geom="line", size = 1, position="identity", show.legend=FALSE) +
    theme_bw() + theme(axis.title.x = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  blank <- ggplot() + geom_point(aes(1,1), colour="white") +
    theme(axis.ticks=element_blank(), panel.background=element_blank(), panel.grid=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), 
          axis.title.y=element_blank())
  
  scatter <- ggplot(data=df, aes_string(x=x, y=y)) + 
    geom_point(size = 0.6) + stat_ellipse(type = "norm", linetype = 2, color="green",size=1) +
    stat_ellipse(type = "t",color="green",size=1) +
    theme_bw() + labs(x=x, y=y)
  
  hright <- ggplot(data=df, aes_string(x=x)) + 
    geom_histogram(aes(y=..density..), fill = "white", color = "black", bins=100) + 
    stat_density(colour = "red", geom="line", size = 1, position="identity", show.legend=FALSE) +
    coord_flip() + theme_bw() + 
    theme(axis.title.y = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  grid.arrange(htop, blank, scatter, hright, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  
}
# simulating data using mvrnorm
generate_gmm_data <- function(input_list) {
  # Initializing variabile
  n <- input_list[[1]]
  mu1 <- input_list[[2]]$mu1
  mu2 <- input_list[[2]]$mu2
  mu3 <- input_list[[2]]$mu3
  
  sigma1_elements <- input_list[[3]]$Sigma1
  sigma2_elements <- input_list[[3]]$Sigma2
  sigma3_elements <- input_list[[3]]$Sigma3
  
  pi <- input_list[[4]]
  
  # Construct covariance matrices from the input vectors
  sigma1 <- matrix(c(sigma1_elements[1], sigma1_elements[3], 
                     sigma1_elements[3], sigma1_elements[2]), nrow = 2)
  
  sigma2 <- matrix(c(sigma2_elements[1], sigma2_elements[3], 
                     sigma2_elements[3], sigma2_elements[2]), nrow = 2)
  
  sigma3 <- matrix(c(sigma3_elements[1], sigma3_elements[3], 
                     sigma3_elements[3], sigma3_elements[2]), nrow = 2)
  
  # Determine number of samples for each component based on mixing probabilities
  n1 <- round(n * pi[1])
  n2 <- round(n * pi[2])
  n3 <- n - (n1 + n2)  # Remaining for the third component
  
  # Generate data for each component
  data1 <- mvrnorm(n = n1, mu = mu1, Sigma = sigma1)
  data2 <- mvrnorm(n = n2, mu = mu2, Sigma = sigma2)
  data3 <- mvrnorm(n = n3, mu = mu3, Sigma = sigma3)
  
  # Combine the data into one dataset and label
  dataset <- rbind(data1, data2, data3)
  labels <- factor(c(rep(1, n1), rep(2, n2), rep(3, n3)))
  
  # Create a data frame with labeled data
  df <- data.frame(X1 = dataset[, 1], X2 = dataset[, 2], Component = labels)
  shuffle_index <- sample(1:n, n)
  df <- df[shuffle_index, ]
  
  return(df)
}

# Example input list with random values
input_list <- list(
  300,  # Sample size
  list(mu1 = c(2, 3), mu2 = c(-2, -1), mu3 = c(4, -3)),  # Means
  list(
    Sigma1 = c(1, 1, 0.8),  # Covariance matrix elements for component 1
    Sigma2 = c(1.5, 1.5, -0.6),  # Covariance matrix elements for component 2
    Sigma3 = c(2, 2, 0.4)  # Covariance matrix elements for component 3
  ),
  c(0.3, 0.4, 0.3)  # Mixing probabilities
)

# Generate the data
simulated_data <- generate_gmm_data(input_list)
View(simulated_data)

# Visualize using the density.hist function
density.hist(simulated_data, x = "X1", y = "X2")
###############################################################################
# OUR SECTION FOR ASSIGNMENT 6
############################################################################
# FUNCTION DETAILS.
# Input:
#   List of two.
#   * Element 1 of the list.
#       epsilon: scalar, convergence tolerance.
#   * Element 2 of the list.
#       x: Vector of 1 x n simulated data fro sim.GMMdata.R
# 
# Returns:
#   List of 3.
#   Maximum likelihood estimates.
#   
#   * Element 2 of the list.
#       List of 3.
#         mu1: 1 x 2 vector, mle of the mean of first component of GMM.
#         mu2: 1 x 2 vector, mle of the mean of second component of GMM.
#         mu3: 1 x 2 vector, mle of the mean of third component of GMM.
#         
#   * Element 3 of the list.
#       List of 3.
#         Sigma1: 1 x 3 vector, mle of the covariance matrix of the first component of GMM.
#         Elements ordered as: var(x1), var(x2), cov(x1,x2).
#         
#         Sigma2: 1 x 3 vector, mle of the covariance matrix of the second component of GMM.
#         Elements ordered as: same as for Sigma1
#         
#         Sigma3: 1 x 3 vector, mle of the covariance matrix of the third component of GMM.
#         Elements ordered as: same as for Sigma1
#     
#     
#     * Element 4 of the list.
#         Vector of 1 x 3.
#         pi = (pi1,pi2,pi3)
#           pi1: scalar, mle of the mixing probability for the first component of GMM.
#           pi2: scalar, mle of the mixing probability for the second component of GMM.
#           pi3: scalar, mle of the mixing probability for the third component of GMM.
#   
############################################################################
# TEST FUNCTION TO DOUBLE CHECK RESULTS
############################################################################
# Fit a GMM using the mclust package
em_model_test <- Mclust(simulated_data[, c("X1", "X2")])

# Display the summary of the fitted model
summary(em_model_test)

# Plot the fitted GMM components
plot(em_model_test, what = "classification")

extract_mles <- function(model) {
  # Extract the MLEs from the model
  means_mle <- model$parameters$mean
  covariances_mle <- model$parameters$variance$sigma
  mixing_proportions_mle <- model$parameters$pro
  
  # Construct the output list with vectors
  result_list <- list(
    means = list(
      mu1 = if (!is.null(means_mle)) as.vector(means_mle[, 1]) else NA,
      mu2 = if (length(means_mle) > 1) as.vector(means_mle[, 2]) else NA,
      mu3 = if (length(means_mle) > 2) as.vector(means_mle[, 3]) else NA
    ),
    covariances = list(
      Sigma1 = if (dim(covariances_mle)[3] >= 1) as.vector(c(covariances_mle[1, 1, 1], covariances_mle[2, 2, 1], covariances_mle[1, 2, 1])) else NA,
      Sigma2 = if (dim(covariances_mle)[3] >= 2) as.vector(c(covariances_mle[1, 1, 2], covariances_mle[2, 2, 2], covariances_mle[1, 2, 2])) else NA,
      Sigma3 = if (dim(covariances_mle)[3] >= 3) as.vector(c(covariances_mle[1, 1, 3], covariances_mle[2, 2, 3], covariances_mle[1, 2, 3])) else NA
    ),
    mixing_proportions = if (!is.null(mixing_proportions_mle)) as.vector(mixing_proportions_mle) else NA
  )
  print(result_list)
  
  return(result_list)
}
# Using function
mle_result <- extract_mles(em_model_test)
print(mle_result)
View(mle_result)

# Function to print table
print_mle_table <- function(mle_result) {
  # Extract data for kable formatting
  means_table <- do.call(rbind, lapply(mle_result$means, function(x) if (is.vector(x)) x else NA))
  rownames(means_table) <- c("mu1", "mu2", "mu3")
  colnames(means_table) <- c("X1", "X2")
  
  covariances_table <- do.call(rbind, lapply(mle_result$covariances, function(x) if (is.vector(x)) x else NA))
  rownames(covariances_table) <- c("Sigma1", "Sigma2", "Sigma3")
  colnames(covariances_table) <- c("var(X1)", "var(X2)", "cov(X1, X2)")
  
  mixing_proportions_table <- data.frame(
    Component = c("pi1", "pi2", "pi3"),
    Mixing_Proportion = mle_result$mixing_proportions
  )
  
  # Print the tables using kable
  cat("MLE for each component's parameter using EM:")
  print(kable(means_table, format = "markdown", caption = "Means (MLE) for Each Component"))
  print(kable(covariances_table, format = "markdown", caption = "Covariance Matrices (MLE) for Each Component"))
  print(kable(mixing_proportions_table, format = "markdown", caption = "Mixing Proportions (MLE)"))
}
# Example of using the function
print_mle_table(mle_result)
############################################################################
# OUR ACTUAL FUNCTION
################################################################################
simGMMData <- function(mu, sigma, pi, n){  
  
  #--------------------------------------------------
  # @description
  # Simulates data from a Gaussian Mixture Model
  #--------------------------------------------------
  # @args in 
  # mu: List of lists containing means for x and y for 3 normal distributions
  # sigma List of lists containing 2 variances then covariance for 3 normal distributions
  # variance must be greater than 0, covariance^2 must be less than sigma1*sigma2
  # pi: selection probabilities for distributions.  Must sum to 1
  # n: the desired length of simulated data
  #--------------------------------------------------  
  # @args out 
  # lists of simulated data
  # outputx: x values simulated from the GMM model
  # outputy: y values simulated from the GMM model
  #--------------------------------------------------  
  # @calls
  # None 
  #-------------------------------------------------
  # @lastChange 2024-11-15
  #--------------------------------------------------
  # @author Brett Flerchinger
  #---------------------------------------------------
  # begin  
  #-------------------------------------------------  
  
  #Error Checking-----------------------------------------
  
  
  if( any(sapply(sigma, is.na)) | length(sigma) != 9){
    stop("Some sigma has an NA or sigma is incorrect length")
  }
  
  if(sigma[1] <= 0){
    stop("Error: bad value for sigma_1 in first distribution")
  }
  
  if(sigma[2] <= 0){
    stop("Error: bad value for sigma_2 in first distribution")
  }
  
  if(sigma[4] <= 0){
    stop("Error: bad value for sigma_1 in second distribution")
  }
  
  if(sigma[5] <= 0){
    stop("Error: bad value for sigma_2 in second distribution")
  }
  
  if(sigma[7] <= 0){
    stop("Error: bad value for sigma_1 in third distribution")
  }
  
  if(sigma[8] <= 0){
    stop("Error: bad value for sigma_2 in third distribution")
  }
  
  if(sigma[2]^2 > (sigma[1]*sigma[3]) ){
    stop("Error: covariance is too high for distribution 1 for given variance values")
  }
  
  if(sigma[6]^2 > (sigma[4]*sigma[5]) ){
    stop("Error: covariance is too high for distribution 1 for given variance values")
  }
  
  if(sigma[9]^2 > (sigma[7]*sigma[8]) ){
    stop("Error: covariance is too high for distribution 1 for given variance values")
  }
  
  if( any(sapply(pi, is.na)) | length(pi) != 3){
    stop("pi has an NA or is incorrect length")
  }
  
  if(abs(1-pi[1]-pi[2]-pi[3]) > .00000001){
    pi[3] = 1-pi[1]-pi[2]
    print("pi values do not sum to 1.  pi_3 set to:")
    print(pi[3])
  }
  
  if(is.na(mu[1]) | is.na(mu[2])){
    stop("mu1 has an NA")
  }
  
  if(is.na(mu[3]) | is.na(mu[4])){
    stop("mu2 has an NA")
  }
  
  if(is.na(mu[5]) | is.na(mu[6])){
    stop("mu3 has an NA")
  }
  
  
  
  #Input Converstion and intialization --------------------------------------------
  mu1 = c(mu[1],mu[2]) #break down mu input
  mu2 = c(mu[3],mu[4])
  mu3 = c(mu[5],mu[6])
  
  sigMatrix1 = matrix(c(sigma[1], sigma[3],sigma[3],sigma[2]), nrow = 2) #convert sigma to matrix
  sigMatrix2 = matrix(c(sigma[4], sigma[6],sigma[6],sigma[5]), nrow = 2)
  sigMatrix3 = matrix(c(sigma[7], sigma[9],sigma[9],sigma[8]), nrow = 2)
  
  randy <- runif(n) #randomly generate
  outputx <- rep(NA, n) #will store x values  
  outputy <- rep(NA, n) #will store y values  
  
  
  sorted <- ifelse(randy < pi[1], 1,
                   ifelse(randy <= (pi[2]+pi[1]), 2, 3)) # assign bivariate distributions
  
  # Main code functionality -----------------------------------------------------
  for(i in 1:n){
    
    if(sorted[i] == 1){
      Z <- matrix(rnorm(2), nrow = 2, ncol = 1) #create Z values
      L <- chol(sigMatrix1) #cholesky to apply variance structure
      X <- mu1 + L %*% Z #convert Z to desired bivariate normal
      outputx[i] = X[1]
      outputy[i] = X[2]
    }#end if
    
    if(sorted[i] == 2){
      Z <- matrix(rnorm(2), nrow = 2, ncol = 1)
      L <- chol(sigMatrix2)
      outputx[i] = X[1]
      outputy[i] = X[2]
    }#end if
    
    if(sorted[i] == 3){
      Z <- matrix(rnorm(2), nrow = 2, ncol = 1)
      L <- chol(sigMatrix3)
      X <- mu3 + L %*% Z
      outputx[i] = X[1]
      outputy[i] = X[2]
    }#end if
    
  }#end for loop
  
  return(list(outputx = outputx, outputy = outputy))
  
}#end function
###################################################################

mu1 = c(0,5) # set mu values
mu2 = c(.2, -7)
mu3 = c(.1,-.01)
mu = c(mu1, mu2, mu3)

pi = c(.7,.2,.1) #set selection probabilities 

sigma1 = c(3,2,1.5) #set sigma matrix values
sigma2 = c(.5,.3,-.35)
sigma3 = c(.2,5,0)

n = 1000

sigma = c(sigma1,sigma2,sigma3)

output <- simGMMData(mu, sigma, pi, n)
output <- as.data.frame(output)
View(output)
#Initialization
set.seed(123)  # For reproducibility
n <- nrow(output)  # Number of data points
k <- 3  # Number of components

# Initial mixing coefficients
pi <- rep(1/k, k)

# Initial mean vectors (randomly selected from data)
mu <- output[sample(1:n, k), ]

# Initial covariance matrices (identity matrices)
Sigma <- lapply(1:k, function(i) diag(ncol(output)))

# Convergence criteria
tolerance <- 1e-6
max_iter <- 100
log_likelihood <- numeric(max_iter)

#E-step

for(i in 1:3)  # This line initializes the for loop to begin running
               # through the each of the three parameters
  
  for(j in 1:n)
    
    pdf_gmm <- (1/sqrt(((2*pi)^2)*det(sigma[i])))*exp^(-0.5*t(x[j]-mu[i])%*%((1/det(sigma[i]))*sigma[i])%*%(x[j]-mu[i]))  # This line defines the PDF of the
                                                                                                                          # bivariate normal distribution.

    delta <- pdf_gmm*pi_gmm[i]/sum(pdf_gmm*pi_gmm[i])  # This line defines the probability that the data point
                                               # x_j belongs to the ith component.
    
    pi[i] <- mean(delta)  # This line updates the mixture probabilities.
    
    mu[i] <- sum(delta*x[j])/sum(delta)  # This line updates the means.
    
    sigma[i] <- sum(delta*(x[j]-mu[i])%*%t(x[j]-mu[i]))/sum(delta)  # This line updates the covariance matrices.

#M-step


#Convergence Check

