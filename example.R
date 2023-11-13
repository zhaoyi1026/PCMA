##################################
# Mediation analysis
# Multiple exposure and multiple mediator

# Example code
##################################

#################################
library("mvtnorm")

rm(list=ls())

source("functions.R")
#################################

#################################
# parameter setting

p<-5
q<-10

#----------------------------
# X parameters

# projection matrix Phi
set.seed(50)
Phi.mat0<-matrix(rnorm(p*p,mean=0,sd=1),nrow=p,ncol=p)
Phi.mat<-qr.Q(qr(Phi.mat0))
Phi<-matrix(NA,nrow(Phi.mat),ncol(Phi.mat))
for(j in 1:ncol(Phi.mat))
{
  if(Phi.mat[which.max(abs(Phi.mat[,j])),j]>0)
  {
    Phi[,j]<-Phi.mat[,j]
  }else
  {
    Phi[,j]<--Phi.mat[,j]
  }
}

if(p>20)
{
  sigma2.X<-sort(exp(c(seq(3,1,length.out=10),seq(-3,-4,length.out=p-10))),decreasing=TRUE) 
}else
{
  sigma2.X<-sort(exp(seq(3,-4,length.out=p)),decreasing=TRUE)
}

# covariance matrix of X
Sigma.X<-Phi%*%diag(sigma2.X)%*%t(Phi)
#----------------------------

#----------------------------
# M parameters

# projection matrix Psi
set.seed(100)
Psi.mat0<-matrix(rnorm(q*q,mean=0,sd=1),nrow=q,ncol=q)
Psi.mat<-qr.Q(qr(Psi.mat0))
Psi<-matrix(NA,nrow(Psi.mat),ncol(Psi.mat))
for(j in 1:ncol(Psi.mat))
{
  if(Psi.mat[which.max(abs(Psi.mat[,j])),j]>0)
  {
    Psi[,j]<-Psi.mat[,j]
  }else
  {
    Psi[,j]<--Psi.mat[,j]
  }
}

# t.X->t.M parameter
alpha0<-c(2,2)
alpha<-matrix(0,nrow=p,ncol=q)
rownames(alpha)<-paste0("CX",1:p)
colnames(alpha)<-paste0("CM",1:q)
if(length(alpha0)==1)
{
  alpha[1,1]<-alpha0
}else
{
  alpha[1:length(alpha0),1:length(alpha0)]<-diag(alpha0) 
}

# M model error parameter
sigma.M<-rep(1,q)
#----------------------------

#----------------------------
# Y parameters

beta0<-c(2,1)
beta<-c(beta0,rep(0,q-length(beta0)))

gamma0<-c(1,-1)
gamma<-c(gamma0,rep(0,p-length(gamma0)))

sigma.Y<-1
#----------------------------

#----------------------------
IE<-alpha%*%diag(beta)
rownames(IE)<-paste0("CX",1:p)
colnames(IE)<-paste0("CM",1:q)
#----------------------------
#################################

#################################
# Generate data
n<-100

set.seed(100)
X<-rmvnorm(n,mean=rep(0,p),sigma=Sigma.X)
X.t<-X%*%Phi

e.M<-rmvnorm(n,mean=rep(0,q),sigma=diag(sigma.M^2))
M.t<-matrix(NA,nrow=n,ncol=q)
if(p==q)
{
  for(k in 1:p)
  {
    M.t[,k]<-X.t[,k]*alpha[k,k]+e.M[,k]
  }
}else
  if(p<q)
  {
    for(k in 1:p)
    {
      M.t[,k]<-X.t[,k]*alpha[k,k]+e.M[,k]
    }
    M.t[,(p+1):q]<-e.M[,(p+1):q]
  }else
  {
    for(k in 1:q)
    {
      M.t[,k]<-X.t[,k]*alpha[k,k]+e.M[,k]
    }
  }
M<-M.t%*%t(Psi)

e.Y<-rnorm(n,mean=0,sd=sigma.Y)
Y<-X.t%*%gamma+M.t%*%beta+e.Y
#################################

#################################
# run proposed method

# method parameters
boot<-TRUE
sims<-1000
boot.ci.type<-c("perc")
conf.level<-0.95
verbose<-TRUE

stop.crt<-c("IE")      # choose the # of components using IE significance
nD<-NULL

re<-NULL
try(re<-HDEM.loglike.opt(X,M,Y,stop.crt=stop.crt,conf.level=conf.level,boot=boot,sims=sims,boot.ci.type=boot.ci.type[1],trace=FALSE,verbose=TRUE))

# correlation of estimated projections and true projections
t(re$Phi)%*%Phi
t(re$Psi)%*%Psi

# inference of the effects
re$coef.inference
#################################
save.image("example.RData")





