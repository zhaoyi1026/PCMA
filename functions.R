###################################
# Mediation analysis with HD exposures and HD mediators
###################################

###################################
# install.packages("quadprog")
# install.packages("ROI")

# library("quadprog")
# library("ROI")

library("MASS")       # general inverse of a matrix
# library("rootSolve")
###################################

###################################
# Given X, M, Y and projections and estimate alpha, beta, and gamma
HDEM.coef<-function(X,M,Y,phi,psi)
{
  n<-length(Y)
  
  x.t<-X%*%phi
  m.t<-M%*%psi
  
  dtmp<-data.frame(X=x.t,M=m.t,Y=Y)
  
  fit.m<-lm(M~X,data=dtmp)
  fit.y<-lm(Y~X+M,data=dtmp)
  
  x.t<-scale(X%*%phi,center=TRUE,scale=FALSE)
  m.t<-scale(M%*%psi,center=TRUE,scale=FALSE)
  y<-scale(Y,center=TRUE,scale=FALSE)
  
  Hx<-diag(rep(1,n))-x.t%*%solve(t(x.t)%*%x.t)%*%t(x.t)
  Hm<-diag(rep(1,n))-m.t%*%solve(t(m.t)%*%m.t)%*%t(m.t)
  
  alpha.est<-c(solve(t(x.t)%*%x.t)%*%(t(x.t)%*%m.t))
  beta.est<-c(solve(t(m.t)%*%Hx%*%m.t)%*%(t(m.t)%*%Hx%*%y))
  gamma.est<-c(solve(t(x.t)%*%Hm%*%x.t)%*%(t(x.t)%*%Hm%*%y))
  
  coef.est<-list(alpha=alpha.est,beta=beta.est,gamma=gamma.est,IE=alpha.est*beta.est)
  
  return(coef.est)
}
###################################

###################################
# Maximize the joint likelihood function

# Joint (negative) likelihood function
HDEM.loglike<-function(X,M,Y,phi,psi,alpha,beta,gamma,sigma2,tau2)
{
  n<-length(Y)
  p<-ncol(X)
  q<-ncol(M)
  
  x.t<-X%*%phi
  m.t<-M%*%psi
  
  s1<-sum((c(m.t-x.t*alpha))^2)/sigma2
  s2<-sum((c(Y-x.t*gamma-m.t*beta))^2)/tau2
  s3<-n*log(sigma2)+n*log(tau2)
  
  return(s1+s2+s3)
}

# Estimate by minimizing the negative joint loglikelihood function
HDEM.loglike.D1.base<-function(X,M,Y,max.itr=1000,tol=1e-4,trace=FALSE,phi0=NULL,psi0=NULL)
{
  # X: n by p exposure matrix
  # M: n by q mediator matrix
  # Y: n by 1 outcome vector
  
  n<-length(Y)
  p<-ncol(X)
  q<-ncol(M)
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",1:p)
  }
  if(is.null(colnames(M))==TRUE)
  {
    colnames(M)<-paste0("M",1:q)
  }
  
  #-----------------------
  # initial values
  # set initial values
  if(is.null(phi0))
  {
    phi0<-rep(1/sqrt(p),p)
  }
  if(is.null(psi0))
  {
    psi0<-rep(1/sqrt(q),q)
  }
  
  alpha0<-0
  beta0<-0
  gamma0<-0
  sigma20<-1
  tau20<-1
  
  if(trace)
  {
    coef.trace<-NULL
    phi.trace<-NULL
    psi.trace<-NULL
    var.trace<-NULL
    loglike<-NULL
  }
  #-----------------------
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update (alpha,beta,gamma)
    coef.tmp<-HDEM.coef(X,M,Y,phi0,psi0)
    alpha.new<-coef.tmp$alpha
    beta.new<-coef.tmp$beta
    gamma.new<-coef.tmp$gamma
    
    # update sigma2 and tau2
    sigma2.new<-sum((M%*%psi0-(X%*%phi0)*alpha.new)^2)/n
    tau2.new<-sum((Y-(X%*%phi0)*gamma.new-(M%*%psi0)*beta.new)^2)/n
    
    # update phi
    phi.lambda.func<-function(lambda)
    {
      U<-t(X)%*%((alpha.new/sigma2.new)*M%*%psi0+(gamma.new/tau2.new)*Y-(beta.new*gamma.new/tau2.new)*M%*%psi0)
      
      otmp<-rep(NA,length(lambda))
      for(i in 1:length(lambda))
      {
        inv.tmp<-ginv((alpha.new^2/sigma2.new+gamma.new^2/tau2.new)*t(X)%*%X+lambda[i]*diag(rep(1,p)))
        otmp[i]<-c(t(U)%*%inv.tmp%*%inv.tmp%*%U)-1
      }
      
      return(otmp)
    }
    re.tmp<-NULL
    try(re.tmp<-uniroot(phi.lambda.func,lower=0,upper=1e6),silent=TRUE)
    if(is.null(re.tmp)==FALSE)
    {
      phi.lambda.val<-re.tmp$root
    }else
    {
      phi.lambda.val<-0
    }
    phi.new<-c(ginv((alpha.new^2/sigma2.new+gamma.new^2/tau2.new)*t(X)%*%X+phi.lambda.val*diag(rep(1,p)))%*%(t(X)%*%((alpha.new/sigma2.new)*M%*%psi0+(gamma.new/tau2.new)*Y-(beta.new*gamma.new/tau2.new)*M%*%psi0)))
    if(sqrt(sum(phi.new^2))!=1)
    {
      phi.new<-phi.new/sqrt(sum(phi.new^2))
    }
    
    # update psi
    psi.lambda.func<-function(lambda)
    {
      V<-t(M)%*%((alpha.new/sigma2.new-beta.new*gamma.new/tau2.new)*X%*%phi.new+beta.new/tau2.new*Y)
      
      otmp<-rep(NA,length(lambda))
      for(i in 1:length(lambda))
      {
        inv.tmp<-solve((1/sigma2.new+beta.new^2/tau2.new)*t(M)%*%M+lambda[i]*diag(rep(1,q))) 
        otmp[i]<-c(t(V)%*%inv.tmp%*%inv.tmp%*%V)-1
      }
      
      return(otmp)
    }
    re.tmp<-NULL
    try(re.tmp<-psi.lambda.val<-uniroot(psi.lambda.func,lower=0,upper=1e6),silent=TRUE)
    if(is.null(re.tmp)==FALSE)
    {
      psi.lambda.val<-re.tmp$root
    }else
    {
      psi.lambda.val<-0
    }
    psi.new<-c(ginv((1/sigma2.new+beta.new^2/tau2.new)*t(M)%*%M+psi.lambda.val*diag(rep(1,q)))%*%(t(M)%*%((alpha.new/sigma2.new-beta.new*gamma.new/tau2.new)*X%*%phi.new+beta.new/tau2.new*Y)))
    if(sqrt(sum(psi.new^2))!=1)
    {
      psi.new<-psi.new/sqrt(sum(psi.new^2))
    }
    
    if(trace)
    {
      coef.trace<-cbind(coef.trace,c(alpha.new,beta.new,gamma.new))
      phi.trace<-cbind(phi.trace,phi.new)
      psi.trace<-cbind(psi.trace,psi.new)
      var.trace<-cbind(var.trace,c(sigma2.new,tau2.new))
      loglike<-c(loglike,HDEM.loglike(X,M,Y,phi.new,psi.new,alpha.new,beta.new,gamma.new,sigma2.new,tau2.new))
    }
    
    diff<-max(abs(alpha0-alpha.new),abs(beta0-beta.new),abs(gamma0-gamma.new))
    # print(diff)
    
    alpha0<-alpha.new
    beta0<-beta.new
    gamma0<-gamma.new
    sigma20<-sigma2.new
    tau20<-tau2.new
    phi0<-phi.new
    psi0<-psi.new
  }
  
  if(phi.new[which.max(abs(phi.new))]<0)
  {
    phi.new<--phi.new
  }
  if(psi.new[which.max(abs(psi.new))]<0)
  {
    psi.new<--psi.new
  }
  
  if(trace)
  {
    rownames(coef.trace)<-c("alpha","beta","gamma")
    rownames(var.trace)<-c("sigma2","tau2")
    
    re<-list(alpha=alpha.new,beta=beta.new,gamma=gamma.new,phi=phi.new,psi=psi.new,sigma2=sigma2.new,tau2=tau2.new,
             logLik=HDEM.loglike(X,M,Y,phi.new,psi.new,alpha.new,beta.new,gamma.new,sigma2.new,tau2.new),
             convergence=(s<max.itr),nitr=s,
             coef.trace=coef.trace,phi.trace=phi.trace,psi.trace=psi.trace,var.trace=var.trace,logLik.trace=loglike)
  }else
  {
    re<-list(alpha=alpha.new,beta=beta.new,gamma=gamma.new,phi=phi.new,psi=psi.new,sigma2=sigma2.new,tau2=tau2.new,
             logLik=HDEM.loglike(X,M,Y,phi.new,psi.new,alpha.new,beta.new,gamma.new,sigma2.new,tau2.new),
             convergence=(s<max.itr),nitr=s)
  }
  
  return(re)
}
# try several initial value and optimize over the objective function
HDEM.loglike.D1<-function(X,M,Y,max.itr=1000,tol=1e-4,trace=FALSE,phi0.mat=NULL,psi0.mat=NULL,ninitial=NULL,seed=100)
{
  # X: n by p exposure matrix
  # M: n by q mediator matrix
  # Y: n by 1 outcome vector
  
  n<-length(Y)
  p<-ncol(X)
  q<-ncol(M)
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",1:p)
  }
  if(is.null(colnames(M))==TRUE)
  {
    colnames(M)<-paste0("M",1:q)
  }
  
  #-----------------------
  # initial values
  # set initial values
  if(is.null(phi0.mat))
  {
    set.seed(seed)
    phi.tmp<-matrix(rnorm((max(p,q)+1+5)*p,mean=0,sd=1),nrow=p)
    phi0.mat<-apply(phi.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
    #--------------------------------
  }
  if(is.null(psi0.mat))
  {
    set.seed(seed)
    psi.tmp<-matrix(rnorm((max(p,q)+1+5)*q,mean=0,sd=1),nrow=q)
    psi0.mat<-apply(psi.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
    #--------------------------------
  }
  if(is.null(ninitial))
  {
    ninitial<-min(c(ncol(phi0.mat),ncol(psi0.mat),10))
  }else
  {
    if(ninitial>min(ncol(phi0.mat),ncol(psi0.mat)))
    {
      ninitial<-min(ncol(phi0.mat),ncol(psi0.mat))
    }
  }
  set.seed(seed)
  phi0.mat<-matrix(phi0.mat[,sort(sample(1:ncol(phi0.mat),ninitial,replace=FALSE))],ncol=ninitial)
  psi0.mat<-matrix(psi0.mat[,sort(sample(1:ncol(psi0.mat),ninitial,replace=FALSE))],ncol=ninitial)
  
  re.tmp<-vector("list",length=ninitial)
  obj<-rep(NA,ninitial)
  for(kk in 1:ninitial)
  {
    try(re.tmp[[kk]]<-HDEM.loglike.D1.base(X,M,Y,max.itr=max.itr,tol=tol,trace=trace,phi0=phi0.mat[,kk],psi0=psi0.mat[,kk]))
    
    if(is.null(re.tmp[[kk]])==FALSE)
    {
      obj[kk]<-re.tmp[[kk]]$logLik
    }
  }
  
  opt.idx<-which.min(obj)[1]
  re<-re.tmp[[opt.idx]]
  
  return(re)
}

# Second and higher direction
HDEM.loglike.Dk<-function(X,M,Y,Phi0=NULL,Psi0=NULL,beta0=NULL,gamma0=NULL,max.itr=1000,tol=1e-4,trace=FALSE,phi0.mat=NULL,psi0.mat=NULL,ninitial=NULL,seed=100)
{
  if(is.null(Phi0)|is.null(Psi0))
  {
    return(HDEM.loglike.D1(X,M,Y,max.itr=max.itr,tol=tol,trace=trace,phi0.mat=phi0.mat,psi0.mat=psi0.mat,ninitial=ninitial,seed=seed))
  }else
  {
    # X: n by p exposure matrix
    # M: n by q mediator matrix
    # Y: n by 1 outcome vector
    
    n<-length(Y)
    p<-ncol(X)
    q<-ncol(M)
    
    if(is.null(colnames(X))==TRUE)
    {
      colnames(X)<-paste0("X",1:p)
    }
    if(is.null(colnames(M))==TRUE)
    {
      colnames(M)<-paste0("M",1:q)
    }
    
    p0<-ncol(Phi0)
    q0<-ncol(Psi0)
    # remove identified components of X
    Xtmp<-X-X%*%(Phi0%*%t(Phi0))
    # remove identified components of M
    Mtmp<-M-M%*%(Psi0%*%t(Psi0))
    # remove identified effect from X and M in Y
    Ytmp<-Y-(X%*%Phi0)%*%gamma0-(M%*%Psi0)%*%beta0
    
    # re<-HDEM.loglike.D1(Xtmp,Mtmp,Y,max.itr=max.itr,tol=tol,trace=trace,phi0.mat=phi0.mat,psi0.mat=psi0.mat,ninitial=ninitial,seed=seed)
    re<-HDEM.loglike.D1(Xtmp,Mtmp,Ytmp,max.itr=max.itr,tol=tol,trace=trace,phi0.mat=phi0.mat,psi0.mat=psi0.mat,ninitial=ninitial,seed=seed)
    
    re$orthogonality<-list(phi=c(c(t(re$phi)%*%Phi0)),psi=c(t(re$psi)%*%Psi0))
    
    return(re)
  }
}
###################################

###################################
# Inference based on bootstrap
HDEM.inf<-function(X,M,Y,phi,psi,boot=TRUE,sims=1000,boot.ci.type=c("perc","SE"),conf.level=0.95,verbose=TRUE)
{
  # X: n by p exposure matrix
  # M: n by q mediator matrix
  # Y: n by 1 outcome vector
  
  n<-length(Y)
  p<-ncol(X)
  q<-ncol(M)
  
  if(boot)
  {
    coef.boot<-matrix(NA,sims,5)
    colnames(coef.boot)<-c("alpha","beta","gamma","IE","DE")
    
    for(b in 1:sims)
    {
      set.seed(100+b)
      idx.tmp<-sample(1:n,n,replace=TRUE)
      
      Xtmp<-X[idx.tmp,]
      Mtmp<-M[idx.tmp,]
      Ytmp<-matrix(Y[idx.tmp,],ncol=1)
      
      otmp<-HDEM.coef(Xtmp,Mtmp,Ytmp,phi,psi)
      coef.boot[b,]<-c(otmp$alpha,otmp$beta,otmp$gamma,otmp$IE,otmp$gamma)
      
      if(verbose)
      {
        print(paste0("Bootstrap sample ",b))
      }
    }
    
    coef.out<-matrix(NA,5,6)
    colnames(coef.out)<-c("Estimate","SE","statistic","pvalue","LB","UB")
    rownames(coef.out)<-colnames(coef.boot)
    coef.out[,1]<-apply(coef.boot,2,mean,na.rm=TRUE)
    coef.out[,2]<-apply(coef.boot,2,sd,na.rm=TRUE)
    coef.out[,3]<-coef.out[,1]/coef.out[,2]
    coef.out[,4]<-(1-pnorm(abs(coef.out[,3])))*2
    if(boot.ci.type[1]=="perc")
    {
      coef.out[,c(5,6)]<-t(apply(coef.boot,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE))
    }
    if(boot.ci.type[1]=="SE")
    {
      coef.out[,5]<-coef.out[,1]-qnorm(1-(1-conf.level)/2)*coef.out[,2]
      coef.out[,6]<-coef.out[,1]+qnorm(1-(1-conf.level)/2)*coef.out[,2]
    }
    
    re<-list(Inference=coef.out,boot.output=coef.boot)
    return(re)
  }
}

# Inference based on asymptotics
HDEM.inf.asmp<-function(X,M,Y,phi,psi,conf.level=0.95)
{
  # X: n by p exposure matrix
  # M: n by q mediator matrix
  # Y: n by 1 outcome vector
  
  n<-length(Y)
  p<-ncol(X)
  q<-ncol(M)
  
  re.coef<-HDEM.coef(X,M,Y,phi,psi)
  
  # variance
  sigma2<-mean((c(M%*%psi)-c(X%*%phi)*re.coef$alpha)^2)
  tau2<-mean((c(Y)-c(X%*%phi)*re.coef$gamma-c(M%*%psi)*re.coef$beta)^2)
  
  # Fisher information matrix of alpha, beta, gamma
  fs.info<-matrix(0,3,3)
  colnames(fs.info)=rownames(fs.info)<-c("alpha","beta","gamma")
  fs.info[1,1]<-c(t(X%*%phi)%*%(X%*%phi))/sigma2
  fs.info[2,2]<-c(t(M%*%psi)%*%(M%*%psi))/tau2
  fs.info[3,3]<-c(t(X%*%phi)%*%(X%*%phi))/tau2
  fs.info[2,3]<-c(t(M%*%psi)%*%(X%*%phi))/tau2
  fs.info[3,2]<-c(t(X%*%phi)%*%(M%*%psi))/tau2
  coef.cov<-ginv(fs.info)
  colnames(coef.cov)=rownames(coef.cov)<-c("alpha","beta","gamma")
  
  coef.out<-matrix(NA,5,6)
  colnames(coef.out)<-c("Estimate","SE","statistic","pvalue","LB","UB")
  rownames(coef.out)<-c("alpha","beta","gamma","IE","DE")
  coef.out[,1]<-c(re.coef$alpha,re.coef$beta,re.coef$gamma,re.coef$IE,re.coef$gamma)
  coef.out["alpha",2]<-sqrt(coef.cov[1,1])
  coef.out["beta",2]<-sqrt(coef.cov[2,2])
  coef.out["gamma",2]<-sqrt(coef.cov[3,3])
  coef.out["IE",2]<-sqrt((coef.out["alpha",2]*coef.out["beta",1])^2+(coef.out["alpha",1]*coef.out["beta",2])^2)
  coef.out["DE",2]<-sqrt(coef.cov[3,3])
  coef.out[,3]<-coef.out[,1]/coef.out[,2]
  coef.out[,4]<-(1-pnorm(abs(coef.out[,3])))*2
  coef.out[,5]<-coef.out[,1]-qnorm(1-(1-conf.level)/2)*coef.out[,2]
  coef.out[,6]<-coef.out[,1]+qnorm(1-(1-conf.level)/2)*coef.out[,2]
  
  re<-list(Inference=coef.out,coef.cov=coef.cov)
  return(re)
}
###################################

###################################
# Maximize the joint likelihood function

# Decide the number of projections based on significance of IE or prespecify
HDEM.loglike.opt<-function(X,M,Y,stop.crt=c("nD","IE"),nD=NULL,conf.level=0.95,boot=TRUE,sims=1000,boot.ci.type=c("perc","SE"),
                           max.itr=1000,tol=1e-4,trace=FALSE,phi0.mat=NULL,psi0.mat=NULL,ninitial=NULL,seed=100,verbose=TRUE)
{
  # X: n by p exposure matrix
  # M: n by q mediator matrix
  # Y: n by 1 outcome vector
  
  n<-length(Y)
  p<-ncol(X)
  q<-ncol(M)
  
  if(stop.crt[1]=="nD"&is.null(nD))
  {
    stop.crt<-"IE"
  }
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",1:p)
  }
  if(is.null(colnames(M))==TRUE)
  {
    colnames(M)<-paste0("M",1:q)
  }
  
  #--------------------------------------------
  # First direction
  tm1<-system.time(re1<-HDEM.loglike.D1(X,M,Y,max.itr=max.itr,tol=tol,trace=trace,phi0.mat=phi0.mat,psi0.mat=psi0.mat,ninitial=ninitial,seed=seed))
  
  # bootstrap inference
  re1.inf<-HDEM.inf(X,M,Y,phi=re1$phi,psi=re1$psi,boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,verbose=FALSE)
  
  Phi.est<-matrix(re1$phi,ncol=1)
  Psi.est<-matrix(re1$psi,ncol=1)
  beta.est<-c(re1$beta)
  gamma.est<-c(re1$gamma)
  
  cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
  rownames(cp.time)<-c("user","system","elapsed")
  
  if(verbose)
  {
    print(paste0("Component ",ncol(Phi.est)))
  }
  #--------------------------------------------
  
  if(stop.crt[1]=="nD")
  {
    re.inference<-vector("list",nD)
    names(re.inference)<-paste0("D",1:nD)
    re.inference[[1]]<-re1.inf$Inference
    
    if(nD>1)
    {
      for(j in 2:nD)
      {
        re.tmp<-NULL
        
        try(tm.tmp<-system.time(re.tmp<-HDEM.loglike.Dk(X,M,Y,Phi0=Phi.est,Psi0=Psi.est,beta0=beta.est,gamma0=gamma.est,max.itr=max.itr,tol=tol,
                                                        trace=trace,phi0.mat=phi0.mat,psi0.mat=psi0.mat,ninitial=ninitial,seed=seed)))
        
        if(is.null(re.tmp)==FALSE)
        {
          Xtmp<-X-X%*%(Phi.est%*%t(Phi.est))
          Mtmp<-M-M%*%(Psi.est%*%t(Psi.est))
          Ytmp<-Y-(X%*%Phi.est)%*%gamma.est-(M%*%Psi.est)%*%beta.est
          re.tmp.inf<-HDEM.inf(Xtmp,Mtmp,Ytmp,phi=re.tmp$phi,psi=re.tmp$psi,boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,verbose=FALSE)
          re.inference[[j]]<-re.tmp.inf$Inference
          
          Phi.est<-cbind(Phi.est,re.tmp$phi)
          Psi.est<-cbind(Psi.est,re.tmp$psi)
          beta.est<-c(beta.est,re.tmp$beta)
          gamma.est<-c(gamma.est,re.tmp$gamma)
          
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
          if(verbose)
          {
            print(paste0("Component ",ncol(Phi.est)))
          }
        }else
        {
          break
        }
      }
    }
    
    colnames(Phi.est)=colnames(Psi.est)<-paste0("D",1:ncol(Phi.est))
    rownames(Phi.est)<-paste0("V",1:p)
    rownames(Psi.est)<-paste0("V",1:q)
    
    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)<-c(paste0("D",1:ncol(Phi.est)),"Total")
    
    if(ncol(Phi.est)>1)
    {
      Phi.orthog<-t(Phi.est)%*%Phi.est
    }else
    {
      Phi.orthog<-1
    }
    if(ncol(Psi.est)>1)
    {
      Psi.orthog<-t(Psi.est)%*%Psi.est
    }else
    {
      Psi.orthog<-1
    }
  }
  if(stop.crt[1]=="IE")
  {
    nD<-1
    
    re.inference<-vector("list",nD)
    names(re.inference)<-paste0("D",1:nD)
    re.inference[[1]]<-re1.inf$Inference
    
    IE.sig.tmp<-as.numeric(re1.inf$Inference["IE","LB"]*re1.inf$Inference["IE","UB"]>0)
    while(IE.sig.tmp==1)
    {
      re.tmp<-NULL
      try(tm.tmp<-system.time(re.tmp<-HDEM.loglike.Dk(X,M,Y,Phi0=Phi.est,Psi0=Psi.est,beta0=beta.est,gamma0=gamma.est,max.itr=max.itr,tol=tol,
                                                      trace=trace,phi0.mat=phi0.mat,psi0.mat=psi0.mat,ninitial=ninitial,seed=seed)))
      
      if(is.null(re.tmp)==FALSE)
      {
        nD<-nD+1
        
        Xtmp<-X-X%*%(Phi.est%*%t(Phi.est))
        Mtmp<-M-M%*%(Psi.est%*%t(Psi.est))
        Ytmp<-Y-(X%*%Phi.est)%*%gamma.est-(M%*%Psi.est)%*%beta.est
        re.tmp.inf<-HDEM.inf(Xtmp,Mtmp,Ytmp,phi=re.tmp$phi,psi=re.tmp$psi,boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,verbose=FALSE)
        
        IE.sig.tmp<-as.numeric(re.tmp.inf$Inference["IE","LB"]*re.tmp.inf$Inference["IE","UB"]>0)
        
        if(IE.sig.tmp==1)
        {
          Phi.est<-cbind(Phi.est,re.tmp$phi)
          Psi.est<-cbind(Psi.est,re.tmp$psi)
          beta.est<-c(beta.est,re.tmp$beta)
          gamma.est<-c(gamma.est,re.tmp$gamma)
          
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
          
          re.inference[[nD]]<-re.tmp.inf$Inference
          
          if(verbose)
          {
            print(paste0("Component ",ncol(Phi.est)))
          }
        }
      }else
      {
        break
      }
    }
    
    names(re.inference)<-paste0("D",1:ncol(Phi.est))
    
    colnames(Phi.est)=colnames(Psi.est)<-paste0("D",1:ncol(Phi.est))
    rownames(Phi.est)<-paste0("V",1:p)
    rownames(Psi.est)<-paste0("V",1:q)
    
    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)<-c(paste0("D",1:ncol(Phi.est)),"Total")
    
    if(ncol(Phi.est)>1)
    {
      Phi.orthog<-t(Phi.est)%*%Phi.est
    }else
    {
      Phi.orthog<-1
    }
    if(ncol(Psi.est)>1)
    {
      Psi.orthog<-t(Psi.est)%*%Psi.est
    }else
    {
      Psi.orthog<-1
    }
  }
  
  re<-list(coef.inference=re.inference,Phi=Phi.est,Psi=Psi.est,Phi.orthogonality=Phi.orthog,Psi.orthogonality=Psi.orthog,time=cp.time)
  
  return(re)
}
###################################

















###################################
# Estimate phi and psi by maximizing |IE|
# Not finished
# HDEM.IEmax<-function(X,M,Y)
# {
#   n<-length(Y)
#   p<-ncol(X)
#   q<-ncol(M)
#   
#   obj.tmp<-function(theta)
#   {
#     phi<-theta[1:p]
#     psi<-theta[(p+1):(p+q)]
#     
#     otmp<-HDEM.coef(X,M,Y,phi,psi)
#     
#     return(abs(otmp$IE))
#   }
#   const.tmp<-function(theta)
#   {
#     phi<-theta[1:p]
#     psi<-theta[(p+1):(p+q)]
#     
#     phi.norm<-sqrt(sum(phi^2))
#     psi.norm<-sqrt(sum(psi^2))
#     
#     return((phi.norm-1)^2+(psi.norm-1)^2)
#   }
#   const.phi<-function(theta)
#   {
#     phi<-theta[1:p]
#     
#     phi.norm2<-sum(phi^2)
#     
#     return(phi.norm2)
#   }
#   const.psi<-function(theta)
#   {
#     psi<-theta[(p+1):(p+q)]
#     
#     psi.norm2<-sum(psi^2)
#     
#     return(psi.norm2)
#   }
#   
#   # prob<-OP(objective=F_objective(obj.tmp,n=p+q),
#   #          constraints=F_constraint(list(const.phi,const.psi),dir=c("==","=="),rhs=c(1,1)),
#   #          bounds=V_bound(ld=-1,ud=1,nobj=p+q),
#   #          maximum=TRUE)
#   # 
#   # otmp<-ROI_solve(prob,solver="nlminb",start=rep(0,p+q))
#   
#   
#   
#   obj.phi<-function(phi)
#   {
#     otmp<-HDEM.coef(X,M,Y,phi,psi)
#     
#     return(abs(otmp$IE))
#   }
#   # prob.phi<-OP(objective=F_objective(obj.phi,n=p),
#   #              constraints=F_constraint(const.phi,dir="==",rhs=1),
#   #              bounds=V_bound(ld=-1,ud=1,nobj=p),
#   #              maximum=TRUE)
# }
###################################








