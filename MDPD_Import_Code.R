#First install these packages: MASS, mvtnorm, tclust, robustbase, matlib
#If you are using linux environment like me you have to use:
# "sudo apt install libcurl4-openssl-dev" in the terminal before installing matlib
library(MASS)
library(mvtnorm)
library(tclust)
library(robustbase)
library(matlib)

#First we will define 3 functions
#Function-1: Returns euclidean distance
distance<-function(x,y){
  z = x-y
  d = sqrt(sum(z^2))
  d
}

#Function-2: Another distance- l_infinity
Distance<-function(x,y){
  z = x-y
  D = max(abs(z))
  D
}

#Density Power Divergence Estimate
DPD_Est<-function(X,n,d,alpha,epsilon){
  Mustart = matrix(0,d,1)
  for (i in 1:d) {
    Mustart[i] = median(X[,i])
  }
  Sigmastart = matrix(0,d,d)
  for(i in 1:d){
    for(j in 1:d){
      if(i==j){
        Sigmastart[i,j] = (1.4826*median(abs(X[,i]-median(X[,i]))))^2
      }
      else{
        Sigmastart[i,j] = (1.4826^2)*median((X[,i]-median(X[,i]))*(X[,j]-median(X[,j])))
      }
    }
  }
  Muhat=Mustart
  Sigmahat=Sigmastart
  Sigmahatnew=matrix(0,d,d)
  iteration=0
  repeat
  {
    W=c()
    for(i in 1:n)
    {
      W[i]=exp(-(alpha/2)*t(X[i,]-Muhat)%*%ginv(Sigmahat)%*%(X[i,]-Muhat))
    }
    w=sum(W)
    sum=rep(0,times=d)
    for(i in 1:n)
    {
      sum=sum+(W[i]*X[i,])
    }
    Muhatnew=sum/w
    
    
    if(det(Sigmahat)<10^(-d))
    {
      Sigmahatnew=Sigmahat
    }
    else
    {
      c=alpha/((alpha+1)^((d+2)/2))
      sum=matrix(0,d,d)
      for(i in 1:n)
      {
        S=(X[i,]-Muhat)%*%t(X[i,]-Muhat)
        sum=sum+(S*W[i])
      }
      Sigmahatnew=sum/(w-c)
    }
    error1=Distance(Muhat,Muhatnew)
    error2=Distance(Sigmahat,Sigmahatnew)
    Muhat=Muhatnew
    Sigmahat=Sigmahatnew
    iteration=iteration+1
    if(error1<epsilon && error2<epsilon)
    {
      break
    }
  }
  DPD_Est=cbind(Muhat,Sigmahat)
  DPD_Est
}

#Clustering Code:
n = 1000
d = 4
thresold = 10^(-5)
k = 3
alpha = 0.3
tr = 0.15
epsilon = 0.01
a = 5
mu=matrix(0,k,d)
mu[1,]=rep(0,times=d)
mu[2,]=rep(a,times=d)
mu[3,]=rep(-a,times=d)
Sigma=3*diag(d)
U=runif(n,0,1)
Actual=c()
p=c(0,0.33,0.66,1)
X=matrix(0,n,d)
for(i in 1:n)
{
  for(j in 1:k)
  {
    if(p[j]<U[i] && U[i]<p[j+1])
    {
      X[i,]=mvrnorm(1,mu[j,],Sigma)
      Actual[i]=j
    }
  }
}
#Initialization: One issue is instead of picking a point from the dataset authors have generated
# a point from the dataset as starting point, Additionally not sure where they added the random
#noise
mu0=matrix(0,k,d)
mu0[1,]=mvrnorm(1,mu[1,],diag(d))
mu0[2,]=mvrnorm(1,mu[2,],diag(d))
mu0[3,]=mvrnorm(1,mu[3,],diag(d))

proportion0=c(1/3,1/3,1/3)
cluster0=c() 
#Construction of Clusters Using Maximum Likelihood Principle
for(i in 1:n)
{
  vector=c()
  for(v in 1:k)
  {
    vector[v]=log(proportion0[v])+log(dmvnorm(X[i,],mu0[v,],diag(d)))
  }
  for(v in 1:k)
  {
    if(max(vector)==vector[v])
    {
      cluster0[i]=v
    }
  }
}

clusterinitial=cluster0

cluster=c()
cluster=clusterinitial
iteration=10
I=1
#Update Step
while(I<iteration+1)
{
  count=c()
  sum=rep(0,times=k)
  for(i in 1:n)
  {
    for(v in 1:k)
    {
      if(cluster[i]==v)
      {
        sum[v]=sum[v]+1
      }
    }
  }
  count=sum
  proportion=count/n #Updated the proportion
  Estimate=matrix(0,k*d,d+1) #Because we are placing the mu and sigma side by side hence d+1
  for(v in 1:k)
  {
    Z=matrix(0,count[v],d)
    l=1
    for(i in 1:n)
    {
      if(cluster[i]==v)
      {
        Z[l,]=X[i,]
        l=l+1
      }
    }
    clusternew=c()
    Estimate[(((v-1)*d)+1):(v*d),]=DPD_Est(Z,count[v],d,alpha,epsilon)
  }
  SS=Estimate[(((v-1)*d)+1):(v*d),2:(d+1)] #Why we are not checking for other v's
  SS1=eigen(SS)$vectors
  SS2=diag(eigen(SS)$values)
  SS3=diag(0,d,d)
  for(f in 1:d)
  {
    if(SS2[f,f]>0.1)
    {
      SS3[f,f]=SS2[f,f]
    }
    else
    {
      SS3[f,f]=0.1
    }
  }
  SSnew=SS1%*%SS3%*%t(SS1)
  Estimate[(((v-1)*d)+1):(v*d),2:(d+1)]=SSnew
  for(i in 1:n)
  {
    vector=c()
    for(v in 1:k)
    {
      vector[v]=proportion[v]*dmvnorm(X[i,],Estimate[(((v-1)*d)+1):(v*d),1],Estimate[(((v-1)*d)+1):(v*d),2:(d+1)])
    }
    for(v in 1:k)
    {
      if(max(vector)==vector[v])
      {
        clusternew[i]=v
      }
    }
  }
  cluster=clusternew
  print(I)  
  I=I+1
}
Clusternew=c()
for(i in 1:n)
{
  v=cluster[i]
  likelihood=proportion[v]*dmvnorm(X[i,],Estimate[(((v-1)*d)+1):(v*d),1],Estimate[(((v-1)*d)+1):(v*d),2:(d+1)])
  if(likelihood<thresold)
  {
    Clusternew[i]=0
  }
  else
  {
    Clusternew[i]=cluster[i]
  }
}