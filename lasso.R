###################################################################
# LASSO Estimation by Using Conic Quadratic Optimization in MOSEK #
###################################################################

library(Rmosek)
library(reshape2)
library(ggplot2)


# LASSO: min 1/n*(||y-X%*%b1-Z%*%b0||_2)^2+lambda*sum(tau*|b1|)

# X is the matrix for the regressors to be penalized; 
# Z is the matrix for the regressors to be always kept; by default, the intercept term. 

# method: "P" for Plasso, "S" for Slasso, and "A" for Alasso. 
# Plasso: tau=1; 
# Slasso: tau=colsd(X); 
# Alasso: tau=|beta1_ols|^(-1). 

# If verbose=0, MOSEK will run quietly. 

lasso=function(y,X,Z=1,lambda,method="P",rtol=1e-6,verbose=0) {
  
  n=length(y)
  k1=ncol(X)
  
  if(is.null(Z)) {
    k0=0
  } else {
    if(length(Z)==1) {
      if(Z==1) {Z=matrix(1,n,1)}
    }
    k0=ncol(Z)
  }
  
  # parameter: 
  # b1+,b1-,b0,v,t,(t-1)/2,(t+1)/2
  # v=y-X%*%b1-Z%*%b0
  # b1=(b1+)-(b1-)
  
  # objective function: 
  # min lambda*sum(tau*|b1|)+1/n*t
  P=list(sense="min")
  if(method=="P") {
    P$c=c(rep(lambda,2*k1),rep(0,k0+n),1/n,0,0)
  } else if(method=="S") {
    P$c=c(rep(lambda*apply(X,2,sd),2),rep(0,k0+n),1/n,0,0)
  } else if(method=="A") {
    XZ=cbind(X,Z)
    P$c=c(rep(lambda*abs(as.vector(solve(t(XZ)%*%XZ)%*%t(XZ)%*%y)[1:k1])^(-1),2),
          rep(0,k0+n),1/n,0,0)
  }
  
  # domain:
  # (b1+),(b1-)>=0
  # t>=0, i.e. (t-1)/2>=-0.5, (t+1)/2>=0.5
  P$bx=rbind(blx=c(rep(0,2*k1),rep(-Inf,k0+n),0,-0.5,0.5),
             bux=rep(Inf,2*k1+k0+n+3))
  
  # linear constraint: 
  # X%*%(b1+)-X%*%(b1-)+Z%*%b0+v=y
  # -0.5*t+(t-1)/2=-0.5
  # -0.5*t+(t+1)/2=0.5
  P$A=rbind(cbind(X,-X,Z,Diagonal(n,1),Matrix(0,n,3)),
            cbind(Matrix(0,2,2*k1+k0+n),rbind(c(-0.5,1,0),
                                              c(-0.5,0,1))))
  P$bc=rbind(blc=c(y,-0.5,0.5),
             buc=c(y,-0.5,0.5))
  
  # quadratic cone constraint:
  # (||v||_2)^2<=t, i.e. (||v||_2)^2+[(t-1)/2]^2<=[(t+1)/2]^2
  # therefore, (t+1)/2>=sqrt{(||v||_2)^2+[(t-1)/2]^2}
  P$cones=matrix(list("QUAD",2*k1+k0+c(n+3,1:n,n+2)),
                 nrow=2,ncol=1)
  rownames(P$cones)=c("type","sub")
  
  # optimization: 
  P$dparam$intpnt_nl_tol_rel_gap=rtol
  result=mosek(P,opts=list(verbose=verbose))
  
  status=result$sol$itr$solsta
  estimate=result$sol$itr$xx
  beta1_hat=estimate[1:k1]-estimate[k1+1:k1]
  residual=estimate[2*k1+k0+1:n]
  
  if(k0==0) {
    return(list(status=status,
                coef_X=beta1_hat,
                residual=residual))
  } else {
    beta0_hat=estimate[2*k1+1:k0]
    return(list(status=status,
                coef_X=beta1_hat,
                coef_Z=beta0_hat,
                residual=residual))
  }
}


# lambda will be set in seq(from=0,to=lambda.to,by=lambda.by). 
# label is the name vector containing coefficients of X. 
# legend.position: "none", "left", "right", "bottom", "top". 

lasso.plot=function(y,X,Z=1,lambda.to,lambda.by,method="P",rtol=1e-6,
                    label=NULL,legend.position="right") {
  
  k1=ncol(X)
  if(is.null(label)) {
    label=paste("beta",1:k1,sep="")
  } else if(length(label)==1) {
    label=paste(label,1:k1,sep="")
  }
  
  beta_hat=NULL
  for(lambda in seq(0,lambda.to,lambda.by)) {
    beta_hat=rbind(beta_hat,
                   t(c(lambda,lasso(y,X,Z,lambda,method,rtol)$coef_X)))
  }
  beta_hat=as.data.frame(beta_hat);colnames(beta_hat)=c("lambda",label)
  
  p=ggplot(melt(beta_hat,id.vars="lambda",variable.name="coef",value.name="est"),
           aes(x=lambda,y=est,color=coef))+
    geom_line()+theme(legend.position=legend.position)
  
  return(list(estimate=beta_hat,plot=p))
}


