#########################################################
# Estimation of Variance-Covariance Matrix of Residuals #
#########################################################

# If bandwidth=h, for any i, Omiga_hat[i,j]=0 if j<i-h or j>i+h. 
# If the error process is weakly stationary, for any i,j,s, cov(e[i],e[i+s])=cov(e[j],e[j+s]). 


library(Matrix)

vcov_residual=function(e,bandwidth=NULL,stationary=TRUE) {
  
  T_n=length(e)
  if(is.null(bandwidth)) {bandwidth=T_n-1}
  
  Omiga=Matrix(0,nrow=T_n,ncol=T_n,sparse=TRUE)
  if(stationary==TRUE) {
    
    for(h in 0:bandwidth) {
      gamma_h=sum(e[1:(T_n-h)]*e[(h+1):T_n])/(T_n-h)
      Omiga=Omiga+
        rbind(cbind(matrix(0,nrow=T_n-h,ncol=h),
                    diag(gamma_h,T_n-h)),
              matrix(0,nrow=h,ncol=T_n))
    }
    
  } else {
    
    for(i in 1:T_n) {
      for(j in i:min(i+bandwidth,T_n)) {Omiga[i,j]=e[i]*e[j]}
    }
    
  }
  Omiga=Omiga+t(Omiga-Diagonal(T_n,diag(Omiga)))
  
  return(Omiga)
}


