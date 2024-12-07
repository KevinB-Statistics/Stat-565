#' EM Algorithm for Gaussian mixture Model (GMM)
#'
#' @description
#' Implements the Expectation-Maximization (EM) algorithm to estimate parameters of a Gaussian Mixture Model (GMM) with three bivariate components.
#'
#' @details The primary input, `x` is a list containing 2 elements:
#'
#' - Element 1: `epsilon` as a scalar value representing convergence tolerance
#'
#' - Element 2: `x` (numeric matrix) of 1 x n simulated data from sim.GMMdata.R
#'
#' @param x A list containing 2 elements as described in the details section
#' @param tolerance Numeric value specifying the convergence tolerance for the EM (epsilon, default = 1e-6)
#' @param max_iter an integer specifying maximum number of iterations (default = 1e4)
#'
#' @returns
#' A list containing 3 elements containing Maximum likelihood estimates:
#'
#' - Element 1 of the list.
#' \itemize{
#'         \item mu1: 1 x 2 vector, mle of the mean of first component of GMM.
#'         \item mu2: 1 x 2 vector, mle of the mean of second component of GMM.
#'         \item mu3: 1 x 2 vector, mle of the mean of third component of GMM.
#'         }
#'
#' - Element 2 of the list.
#' \itemize{
#'         \item Sigma1: 1 x 3 vector, mle of the covariance matrix of the first component of GMM. Elements ordered as: var(x1), var(x2), cov(x1,x2).
#'         \item Sigma2: 1 x 3 vector, mle of the covariance matrix of the second component of GMM. Elements ordered as: same as for Sigma1
#'         \item Sigma3: 1 x 3 vector, mle of the covariance matrix of the third component of GMM.  Elements ordered as: same as for Sigma1
#'         }
#'
#' - Element 3 of the list.
#' \itemize{
#'         \item Vector of 1 x 3. pi = (pi1,pi2,pi3)
#'         \item pi1: scalar, mle of the mixing probability for the first component of GMM.
#'         \item pi2: scalar, mle of the mixing probability for the second component of GMM.
#'         \item pi3: scalar, mle of the mixing probability for the third component of GMM.
#'         }
#'
#' @examples
#' # Create input data for function
#' # 'x' is a list with:
#' #  1. A tolerance value (epsilon)
#' #  2. Simulated data matrix
#'
#' # Simulate data
#' dmatrix <- matrix(rnorm(200, 0, 1), nrow = 100, ncol = 2)
#' data <- list(epsilon = 1e-6, dmatrix)
#'
#' # Run the EM algorithm with default parameters
#' result_default <- EM(x = data)
#' print(result_default)
#'
#' # Run the EM algorithm with custom parameters
#' result_custom <- EM(x = data, tolerance = 1e-4, max_iter = 5000)
#' print(result_custom)
#'
#'
#' ##########################################################################
#' @export
EM <- function(x, tolerance = 1e-6, max_iter = 1e4){
  
  if(typeof(x) != "list"){
    stop("Data must be a list")
  }
  
  if(length(x) != 2){
    stop("Length of data must be 2 for bivariate components")
  }
  
  x <- as.matrix(as.data.frame(x))
  
  n <- nrow(x)  # Number of data points
  k <- 3  # Number of components
  
  # Initial mixing coefficients
  pi_gmm <- rep(1/k, k)
  
  # Initial mean vectors (randomly selected from data)
  mu <- x[sample(1:n, k), ]
  mu <- as.matrix(mu)
  
  # Initial covariance matrices (identity matrices)
  sigma <- lapply(1:k, function(i) diag(ncol(x)))
  
  log_likelihood <- numeric(max_iter)
  
  # Function for getting density of GMM of bivariate 
  dMVN <- function(x, mu, sigma) {
    if (ncol(x) != 2){
      stop("Not bivariate data")
    }
    det_sigma <- det(sigma)
    if(det_sigma <= 0) {
      det_sigma <- 1e-6  # Avoids issues
    }
    inv_sigma <- solve(sigma)
    const <- (2 * pi)^(-1) * det_sigma^(-0.5)
    x_minus_mu <- sweep(x, 2, mu, "-")  # x - mu
    dist_sq <- rowSums((x_minus_mu %*% inv_sigma) * x_minus_mu)
    exp_term <- exp(-0.5 * dist_sq)
    return(const * exp_term)
  }
  
  for (iter in 1:max_iter) {
    # E-step: Calculate mixture probabilities
    delta <- matrix(0, nrow = n, ncol = k)
    for (j in 1:3) {
      delta[, j] <- pi_gmm[j] * dMVN(x, as.matrix(mu[j, ]), sigma[[j]])
    }
    log_likelihood[iter] <- sum(log(rowSums(delta)))
    
    delta <- delta / rowSums(delta)  # This line defines the probability that the data point
                                     # x_j belongs to the ith component. Normalized so each row sums to 1 
  
    # M-step: Update parameters
    N_k <- colSums(delta)
    pi_gmm <- N_k / n
    
    for (j in 1:k) {
      
      mu[j, ] <- colSums(delta[, j] * x) / N_k[j]  # This line updates the means.
      
      x_centered <- sweep(x, 2, mu[j, ], "-")  # This line centers the data.
      
      sigma[[j]] <- (t(x_centered * delta[, j]) %*% x_centered) / N_k[j]  # This line updates the covariance matrices.
      
    }
    
    # Check for convergence to update mu, sigma, and pi using log-likelihood:
    if (log_likelihood[2]!=0 && abs(log_likelihood[iter] - log_likelihood[iter-1]) <= tolerance){
      log_likelihood <- log_likelihood[1:iter]
      break
    }

  }
  
  mu_final <- list(mu_1 = mu[1,], mu_2 = mu[2,], mu_3 = mu[3,])
    
  sigma_final <- list(c(var_x1 = sigma[[1]][1,1], var_x2 = sigma[[1]][2,2], cov_x1x2 = sigma[[1]][1,2]),
                      c(var_x1 = sigma[[2]][1,1], var_x2 = sigma[[2]][2,2], cov_x1x2 = sigma[[2]][1,2]),
                      c(var_x1 = sigma[[3]][1,1], var_x2 = sigma[[3]][2,2], cov_x1x2 = sigma[[3]][1,2]))
    
  pi_gmm_final <- list(pi_1 = pi_gmm[1], pi_2 = pi_gmm[2], pi_3 = pi_gmm[3])
  
  return(list(mu_final, sigma_final, pi_gmm_final))  # This line returns the outputs as a list.

}
