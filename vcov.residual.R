#########################################################
# Estimation of Variance-Covariance Matrix of Residuals #
#########################################################

# If bandwidth=h, for any i, Omiga_hat[i,j]=0 if j<i-h or j>i+h. 
# If the error process is weakly stationary, for any i,j,s, cov(e[i],e[i+s])=cov(e[j],e[j+s]). 


vcov.residual=function(e,bandwidth=NULL,stationary=TRUE) {
  
  library(Matrix)
  
  n=length(e)
  
  if(stationary==TRUE) {
    
    if(is.null(bandwidth)) {
      bandwidth=n-1
      Sigma=Diagonal(1,e[1]*e[n])
    } else {
      bandwidth=min(bandwidth,n-1)
      Sigma=Diagonal(n-bandwidth,mean(e[1:(n-bandwidth)]*e[(bandwidth+1):n]))
    }
    
    for(h in (bandwidth-1):1) {
      Sigma=rbind(cbind(0,Sigma),0)+Diagonal(n-h,mean(e[1:(n-h)]*e[(h+1):n]))
    }
    Sigma=rbind(cbind(0,Sigma),0)
    Sigma=Sigma+t(Sigma)+Diagonal(n,mean(e^2))
    
  } else {
    
    Sigma=e%*%t(e)
    
  }
  
  return(Sigma)
}


