LS.group=function(y,X,group) {
  
  library(Matrix)
  
  K=ncol(X);n=length(y)
  
  # OLS
  
  XX=solve(t(X)%*%X)
  beta_OLS=as.vector(XX%*%t(X)%*%y)
  names(beta_OLS)=colnames(X)
  
  e=y-as.vector(X%*%beta_OLS)
  sigma2_hat=aggregate(e^2,by=list(group),FUN=mean)
  sigma2_hat=plyr::join(data.frame(Group.1=group),sigma2_hat)
  Sigma_hat=Diagonal(x=sigma2_hat$x)
  
  se.OLS=sqrt(diag(XX%*%t(X)%*%Sigma_hat%*%X%*%XX))
  t.OLS=beta_OLS/se.OLS
  p.OLS=pnorm(-abs(t.OLS))*2
  
  # FGLS
  
  Omega=Sigma_hat/mean(e^2)
  XSX=solve(t(X)%*%solve(Omega)%*%X)
  beta_GLS=as.vector(XSX%*%t(X)%*%solve(Omega)%*%y)
  names(beta_GLS)=colnames(X)
  
  u=y-as.vector(X%*%beta_GLS)
  se.GLS=sqrt(sum(u^2)/(n-K)*diag(XSX))
  t.GLS=beta_GLS/se.GLS
  p.GLS=pnorm(-abs(t.GLS))*2
  
  return(list(OLS=cbind(coef=beta_OLS,se=se.OLS,t=t.OLS,`p-value`=p.OLS),
              FGLS=cbind(coef=beta_GLS,se=se.GLS,t=t.GLS,`p-value`=p.GLS)))
}


