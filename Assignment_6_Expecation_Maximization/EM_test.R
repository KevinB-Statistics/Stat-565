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

#HI