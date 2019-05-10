#########################################################
# Estimation of Variance-Covariance Matrix of Residuals #
#########################################################

# If the error process is weakly stationary, for any i,j,s, cov(e[i],e[i+s])=cov(e[j],e[j+s]). 
# If bandwidth=h, for any i, Omiga_hat[i,j]=0 if j<i-h or j>i+h. 
# bandwidth="MA" gives automatic order of a stationary process. 


vcov.residual=function(e,stationary=TRUE,bandwidth=NULL,alpha=0.05) {
  
  library(Matrix)
  
  n=length(e)
  
  if(stationary==TRUE) {
    
    ACF=acf(e,lag.max=n-1,type="covariance",plot=FALSE,demean=FALSE)
    ACF=data.frame(h=ACF$lag,gamma=ACF$acf)
    
    if(is.null(bandwidth)) {
      
      Sigma=toeplitz(ACF$gamma)
      
    } else {
      
      if(bandwidth=="MA") {
        
        ACF$rho=ACF$gamma/ACF$gamma[1]
        ACF$cutoff=abs(sqrt(n)*ACF$rho)<=qnorm(1-alpha/2)
        q=min(ACF$h[ACF$cutoff])-1
        Sigma=toeplitz(sparseVector(ACF$gamma[1:(1+q)],i=1:(1+q),length=n))
        
      } else {
        
        q=floor(max(min(bandwidth,n),0))
        Sigma=toeplitz(sparseVector(ACF$gamma[1:(1+q)],i=1:(1+q),length=n))
        
      }
    }
    
  } else {
    
    if(is.null(bandwidth)) {
      
      Sigma=e%*%t(e)
      
    } else {
      
      q=floor(max(min(bandwidth,n),0))
      
      if(q==0) {
        
        Sigma=Diagonal(n,e^2)
        
      } else {
        
        Sigma=Matrix(0,n,n,sparse=TRUE)
        for(i in 1:(n-1)) {
          for(j in (i+1):min(i+q,n)) {
            Sigma[i,j]=e[i]*e[j]
          }
        }
        Sigma=Diagonal(n,e^2)+Sigma+t(Sigma)
        
      }
    }
  }
  
  return(Sigma)
}


