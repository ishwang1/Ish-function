###########################
# DGP of an ARIMA process #
###########################

# u_t~ARIMA(p,d,q): D^d(u_t)=L^(1:p)(u_t)%*%phi+e_t+L^(1:q)(e_t)%*%theta, 
# e_t~WN(0,sigma^2). 

# By default, if phi and theta are NULL and d=0, returning white noise process; 
# if length(phi)=p, theta is null and d=0, returning AR(p) process;
# if length(theta)=q, phi is null and d=0, returning MA(q) process; 
# if length(phi)=p, length(theta)=q and d=0, returning ARMA(p,q) process;
# if length(phi)=p, length(theta)=q and d>0, returning ARIMA(p,d,q) process. 


DGP.ARIMA=function(n,phi=NULL,theta=NULL,d=0,sigma=1,drop=0,seed=NULL) {
  
  T_n=n+drop
  p=length(phi);q=length(theta)
  rev_phi=rev(phi);rev_theta=rev(theta)
  
  set.seed(seed=seed)
  e=rnorm(T_n+q,sd=sigma)
  
  u=NULL
  if(p==0) {
    if(q==0) {  # WN
      u=e
    } else {  # MA
      for(t in 1:T_n) {
        u=c(u,
            e[q+t]+rev_theta%*%e[1:q+t-1])
      }
    }
  } else {
    u_0=rep(0,p)
    if(q==0) {  # AR
      for(t in 1:T_n) {
        u=c(u,
            rev_phi%*%c(u_0,u)[length(u)+1:p]+e[q+t])
      }
    } else {  # ARMA
      for(t in 1:T_n) {
        u=c(u,
            rev_phi%*%c(u_0,u)[length(u)+1:p]+e[q+t]+rev_theta%*%e[1:q+t-1])
      }
    }
  }
  
  if(d!=0) {  # ARIMA
    for(i in 1:d) {
      u=cumsum(u)
    }
  }
  
  y=u[(drop+1):length(u)]
  return(y)
}


