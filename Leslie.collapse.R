#Leslie.collapse
#AUTHOR: Richard A. Hinrichsen
#CONTACT: rich@hinrichsenenvironmental.com
#This code is based on the following publication
#Hinrichsen, R. A. (2023). Aggregation of Leslie matrix models with 
#application to ten diverse animal species. Population Ecology, 1–21. 
#https://doi.org/10.1002/1438-390X.12149

library("expm")

is.naturalnumber=function(x, tol = .Machine$double.eps^0.5){
  test=(x > tol) & (abs(x - round(x)) < tol)
  return(test)
}
 
#Aggregate a Leslie matrix 
#inputs
#A is a Leslie matrix
#m is the dimensionality of the desired aggregated matrix
#weights are the least squares weights
#and are equal to the stable age distribution
#outputs
#B is the aggregated matrix
#W is the weight matrix
#Ak is the original (or expanded) Leslie matrix raised to the k power
#S is the partitioning matrix
#n is the size of the original (or expanded) Leslie matrix
#m is the size of the aggregated matrix
#k is n/m
#EFF is the effectiveness of aggregation
Leslie.collapse = function(A,m){
 #Check whether the matrix is a Leslie matrix
 if(!is.leslie(A))stop("A is not a Leslie matrix")
 if(!is.naturalnumber(m))stop("m must be a natural number")
 n=dim(A)[1]
 if(n%%m!=0){
   A=expand.leslie(A,m)
 } #if
 n=dim(A)[1]
 k=n/m
 eigA=eigen(A)
 iii=(Re(eigA$values)==max(Re(eigA$values)))
 v=abs(Re(eigA$vectors[,iii]))
 w=v #weights are equal to the stable age distribution
 Ak=A%^%k
 e = rep(1,k)
 Id = diag(1,m)
 S = kronecker(Id,t(e))
 W = diag(w)
 T = W%*%t(S)%*%solve(S%*%W%*%t(S))
 B = S%*%Ak%*%T
 EFF=norm(B%*%S%*%sqrt(W),type="F")^2/norm(S%*%Ak%*%sqrt(W),type="F")^2
 return(list(B=B,W=W,Ak=Ak,S=S,n=n,m=m,k=k,EFF=EFF))
}


#disaggregation step
expand.leslie=function(A,m){
if(!is.leslie(A))return("A must be a Leslie matrix")
if(!is.naturalnumber(m))stop("m must be a natural number")
n=nrow(A)
A.expanded=matrix(0,nrow=n*m,ncol=n*m)
#fill in first row of A.expanded
for(ii in 1:n){
  A.expanded[1,ii*m]=A[1,ii]
}
#next fill in the subdiagonal
for(ii in 2:(n*m)){
 A.expanded[ii,ii-1]=1
}
for(ii in 1:(n-1)){
 A.expanded[m*ii+1,m*ii]=A[ii+1,ii]}
 return(A.expanded)
}

#check whether a matrix is a Leslie matrix
#survival probabilities must be positive
#at least one fertility rate must be positive
#otherwise stable age distribution is not positive
is.leslie=function(A){
 EPS = .Machine$double.eps^0.5
 if(is.vector(A)){
  if(length(A)>1){return(FALSE)}
  if(length(A)==1){if(A<=0)return(FALSE);A=as.matrix(A)};
 }
 dims=dim(A)
 if(nrow(A)!=ncol(A))return(FALSE)
 n=ncol(A)
 if((n==1)&(A[1,1]>0))return(TRUE)
 if(is.na(sum(A)))return(FALSE)
 fertility.rates=A[1,]
 if(n>2){survival.probs=diag(A[2:n,1:(n-1)])}
 if(n==2){survival.probs=A[2,1]}

 Anew=matrix(0,nrow=n,ncol=n)
 if(n>2)Anew[2:n,1:(n-1)]=diag(survival.probs)
 if(n==2){Anew[2,1]=survival.probs}
 Anew[1,]=fertility.rates

 error=c(A)-c(Anew)
 error.squared=sum(error*error)
 test1=(sqrt(error.squared)<EPS)
 iii=(survival.probs>0)&(survival.probs<=1)
 test2=(sum(iii)==(n-1))
 iii=(fertility.rates>=0)
 test3=(sum(iii)==n)
 test4=sum(fertility.rates)>0
 return(test1&test2&test3&test4)
}

