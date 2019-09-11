arreg=function(y,order.max=10,HQ.c=2) {
  
  library(foreach)
  
  n=length(y)
  LY=NULL
  for (p in 1:order.max) {
    LY=cbind(LY,c(rep(NA,p),y[-((n-p+1):n)]))
  }
  
  
  
  arreg.order=function(p) {
    if(p>0) {
      y_tilde=y[(p+1):n]
      X=LY[(p+1):n,1:p,drop=FALSE]
      XX=t(X)%*%X
      if(min(eigen(XX)$value)>=.Machine$double.eps) {
        alpha=as.vector(solve(XX)%*%t(X)%*%y_tilde)
        sigma2=mean((y_tilde-as.vector(X%*%alpha))^2)
        AIC=log(sigma2)+2*p/(n-p)
        BIC=log(sigma2)+p*log(n-p)/(n-p)
        HQIC=log(sigma2)+HQ.c*p*log(log(n-p))/(n-p)
        return(c(alpha,rep(NA,order.max-p),sigma2,AIC,BIC,HQIC))
      } else {
        return(rep(NA,order.max+3))
      }
    } else {
      sigma2=mean(y^2)
      return(c(rep(NA,order.max),sigma2,rep(log(sigma2),3)))
    }
  }
  
  result=foreach(p=0:order.max,.combine=rbind) %do% {arreg.order(p)}
  rownames(result)=paste("AR(",0:order.max,")",sep="")
  colnames(result)=c(paste("a",1:order.max,sep=""),"sigma2","AIC","BIC","HQIC")
  
  return(result)
}

