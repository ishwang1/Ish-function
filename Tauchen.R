Tauchen=function(epsilon_cond,n=9,m=3,rho=0,alpha=0,sigma=1) {
  
  mu=alpha/(1-rho)
  mean_cond=alpha+rho*epsilon_cond
  
  k1=mu-m*sigma
  kn=mu+m*sigma
  
  epsilon=seq(from=k1,to=kn,length.out=n)
  w=(kn-k1)/(n-1)/2
  
  p1=pnorm(k1+w,mean=mean_cond,sd=sigma)
  p=pnorm(epsilon[2:(n-1)]+w,mean=mean_cond,sd=sigma)-pnorm(epsilon[2:(n-1)]-w,mean=mean_cond,sd=sigma)
  pn=1-pnorm(kn-w,mean=mean_cond,sd=sigma)
  
  return(data.frame(epsilon=epsilon,CondProb=c(p1,p,pn)))
}


