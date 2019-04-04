empirical_test=function(t,x,alpha=0.05,type="two") {
  
  library(ggplot2)
  g=ggplot(data.frame(x),aes(x))+stat_density(fill="lightyellow",alpha=0.5)+
    labs(x=NULL,y=NULL)+theme(legend.position="bottom")
  
  if(type=="two") {
    CI.lower=quantile(x,alpha/2)
    CI.upper=quantile(x,1-alpha/2)
    p=2*min(c(mean(x>=t),mean(x<=t)))
    vline=data.frame(Value=c("t","CI.lower","CI.upper"),x=c(t,CI.lower,CI.upper))
    g=g+geom_vline(data=vline,aes(xintercept=x,color=Value),size=1)
  }
  if(type=="right") {
    CI.lower=min(x)
    CI.upper=quantile(x,1-alpha)
    p=mean(x>=t)
    vline=data.frame(Value=c("t","CI.upper"),x=c(t,CI.upper))
    g=g+geom_vline(data=vline,aes(xintercept=x,color=Value),size=1)
  }
  if(type=="left") {
    CI.lower=quantile(x,alpha)
    CI.upper=max(x)
    p=mean(x<=t)
    vline=data.frame(Value=c("t","CI.lower"),x=c(t,CI.lower))
    g=g+geom_vline(data=vline,aes(xintercept=x,color=Value),size=1)
  }
  
  result=data.frame(t,p.value=p,median=median(x),CI.lower,CI.upper)
  rownames(result)=NULL
  
  return(list(test=result,plot=g))
}

