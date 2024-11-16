library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

source("sim.GMMdata.R")

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

list1 <- output[[1]]
list2 <- output[[2]]


cat("Head of x:\n")
print(head(list1))

cat("Head of y:\n")
print(head(list2))