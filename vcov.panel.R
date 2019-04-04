vcov.panel=function(X) {
  
  m=nrow(X);n=ncol(X)
  mu=mean(X)
  
  S=matrix(NA,nrow=m,ncol=n)
  for (j in 0:(n-1)) {
    for (i in 0:(m-1)) {
      A=X[1:(m-i),1:(n-j)];B=X[(1+i):m,(1+j):n]
      S[i+1,j+1]=sum((A-mu)*(B-mu))/((m-i)*(n-j)-1)
    }
  }
  S[m,n]=NA
  
  return(S)
}

