#Leslie.collapse
#AUTHOR: Richard A. Hinrichsen
#CONTACT: rich@hinrichsenenvironmental.com
library("expm")


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
 #first check whether the matrix is a Leslie matrix
 if(!is.leslie(A))stop("A is not a Leslie matrix")
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
n=dim(A)[1]
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
is.leslie=function(A){
EPS=.00001
if(is.na(sum(A)))return(FALSE)
m=dim(A)[1]
if(m==1)return(TRUE)
birth.rates=A[1,]
if(m>2){survival.rates=diag(A[2:m,1:(m-1)])}
if(m==2){survival.rates=A[2,1]}

Anew=matrix(0,nrow=m,ncol=m)
if(m>2)Anew[2:m,1:(m-1)]=diag(survival.rates)
if(m==2){Anew[2,1]=survival.rates}
Anew[1,]=birth.rates

error=c(A)-c(Anew)
error.squared=t(error)%*%error
test1=(sqrt(error.squared)<EPS)
iii=(survival.rates>=0)&(survival.rates<=1)
test2=(sum(iii)==(m-1))
iii=(birth.rates>=0)
test3=(sum(iii)==m)
#return(list(test1=test1,test2=test2,test3=test3))
return(is.logical(test1&test2&test3))
}

