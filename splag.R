#######################################
# Generating Weight Matrix from Edges #
#######################################

# edge should be a n*2 matrix, i.e. each row is an edge denoted by two note id codes; 
# strength is the weight vector of each edge, equal by default; 
# n is the number of nodes, i.e. nodes should have id codes from 1 to n; 
# If edges do not have directions and you only have half of edges, mirror should be TRUE; 
# keep rowNorm be TRUE for row normality. 

# If degree is TRUE, the result will contain the degree vector. 
# By default, if sparse is TRUE, the generated weight matrix will be a sparse matrix. 


library(Matrix)

WeightMatrix=function(edge,strength=1,n=NULL,mirror=FALSE,rowNorm=TRUE,
                      degree=FALSE,sparse=TRUE) {
  
  if(is.null(n)) {
    n=max(edge[1],edge[2])
  } else if(n<max(edge[1],edge[2])) {
    stop("The dimension of weight matrix should be no less than the number of nodes. ")
  }
  
  W=Matrix(0,nrow=n,ncol=n,sparse=sparse)
  if(length(strength)==1) {
    for(i in 1:nrow(edge)) {
      W[edge[i,1],edge[i,2]]=1
    }
  } else if(length(strength)==nrow(edge)) {
    for(i in 1:nrow(edge)) {
      W[edge[i,1],edge[i,2]]=strength[i]
    }
  } else {
    stop("The dimension of 'strength' vector should either be 1 or equal to the number of edges. ")
  }
  
  if(mirror==TRUE) {W=W+t(W)}
  
  if(rowNorm==TRUE) {
    D=rowSums(W)
    norm_factor=D
    norm_factor[norm_factor!=0]=1/norm_factor[norm_factor!=0]
    W=Diagonal(n,norm_factor)%*%W
    if(degree==TRUE) {
      return(list(W=W,degree=D))
    } else {
      return(W)
    }
  } else {
    if(degree==TRUE) {
      return(list(W=W,degree=rowSums(W)))
    } else {
      return(W)
    }
  }
  
}




##############################################
# Generating Spatial Lag when Outcome has NA #
##############################################


splag=function(W,x) {
  
  na_index=which(is.na(x))
  W1=W[,-na_index]
  
  norm_factor=rowSums(W1)
  norm_factor[norm_factor!=0]=1/norm_factor[norm_factor!=0]
  W1=Diagonal(length(norm_factor),norm_factor)%*%W1
  
  Lx=as.vector(W1%*%x[-na_index])
  return(Lx)
}




