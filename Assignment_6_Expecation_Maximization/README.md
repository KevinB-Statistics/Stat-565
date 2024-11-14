# Stat-565
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
