################################
# Simulating Multi-Normal Data #
################################

# mu is the mean vector; 
# Sigma is the Var-Cov matrix. 


simulate.MN=function(mu,Sigma,seed=NULL) {
  
  n=length(mu)
  
  eigen.Sigma=eigen(Sigma)
  P=eigen.Sigma$vectors%*%diag(sqrt(eigen.Sigma$values))
  
  set.seed(seed)
  e=rnorm(n)
  
  x=P%*%e+mu
  return(x)
}

