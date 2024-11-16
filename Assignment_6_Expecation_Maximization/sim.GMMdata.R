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