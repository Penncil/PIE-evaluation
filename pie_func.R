library(numDeriv)
library(MASS)

################### FUNCTIONS ######################
expit <- function(t) exp(t)/(1+exp(t))
logit <- function(prob) log(prob/(1-prob)) 

## observed negative log likelihood
## t = (alpha1,alpha0,beta0,beta1)
nloglik=function(t)
{
  tem=t[3]+t[4]*x
  p=(1-t[2])+(t[1]+t[2]-1)*exp(tem)/(1+exp(tem))
  loglikelihood=s*log(p)+(1-s)*log(1-p)
  return(-sum(loglikelihood))
}

## prior = uniform distribution.
dprior.u=function(x,a,b){
  as.numeric(x>=a & x<=b)/(b-a)
}

NaiveVsPIE.uni.norm <- function(s, x, adjust=NULL){
  n <- length(s)
  ## naive estimator
  ## observed negative log likelihood, ignoring errors
  ## t = (beta0,beta1)
  nloglik.naive = function(t)
  {
    tem=t[1]+t[2]*x
    p=exp(tem)/(1+exp(tem))
    loglikelihood=s*log(p)+(1-s)*log(1-p)
    return(-sum(loglikelihood))
  }
  est.naive = nlminb(c(0,0),nloglik.naive,lower=c(-5,-5),upper=c(5,5))$par
  cov.est.naive = ginv(hessian(nloglik.naive, est.naive))
  
  ## PIE estimator of beta
  norm.int.lik.u=function(t)
  {
    ngrid=10
    a1.grid = seq(low1,up1,length = ngrid)
    a2.grid = seq(low2,up2,length = ngrid)
    
    prod.lik=array(dim=c(ngrid,ngrid))
    for(i in 1:ngrid){
      for(j in 1:ngrid){
        tem=t[1]+t[2]*x
        p=(1-a2.grid[j])+(a1.grid[i]+a2.grid[j]-1)*exp(tem)/(1+exp(tem))
        #### NUMERICAL ISSUE, need to time a number (e.g., 1.4) to avoid -inf after multiplication
        #### As the no. of sample size increases and the no. of prevalence decreases, 
        #### the size of adjustment needed is larger.
        #### Practically, may need to try a few choices.
        likelihood=(adjust*p)^s*(adjust-adjust*p)^(1-s) 
        prod.lik[i,j] = prod(likelihood)*dprior.u(a1.grid[i],low1,up1)*dprior.u(a2.grid[j],low2,up2)
      }
    }
    return(-log(sum(prod.lik)))
  }
  
  ## prior distributions for sensitivity and specificity both uniform
  t.s <- proc.time()
  est.int.u= nlminb(c(0,0),norm.int.lik.u,lower=c(-5,-5),upper=c(5,5))$par
  t.run <- proc.time()-t.s
  
  cov.est.int.u = ginv(hessian(norm.int.lik.u, est.int.u))
  
  nloglik.gs <- function(t)
  {
    tem=t[1]+t[2]*x
    p=(1-alpha0)+(alpha0+alpha1-1)*exp(tem)/(1+exp(tem))
    loglikelihood=s*log(p)+(1-s)*log(1-p)
    return(-sum(loglikelihood))
  }
  est.gs = nlminb(c(0,0),nloglik.gs,lower=c(-5,-5),upper=c(5,5))$par
  cov.est.gs = ginv(hessian(nloglik.gs, est.gs))
  
  list(out = c(est.int.u,cov.est.int.u[1,1], cov.est.int.u[2,2], cov.est.int.u[1,2],
               est.naive, cov.est.naive[1,1], cov.est.naive[2,2], cov.est.naive[1,2], 
               est.gs,cov.est.gs[1,1], cov.est.gs[2,2], cov.est.gs[1,2]),
       run.time = t.run[3])
}

## prior=logit normal distribution.

dprior=function(x,mu,sigma){
  y = (x-0.5)*2
  dlogitnorm(y, logit(mu),sigma)
}

PIE.logitnorm <- function(s, x, adjust=NULL){
  n <- length(s)
  ## t = (beta0,beta1)

  ## PIE estimator of beta
  norm.int.lik=function(t)
  {
    ngrid=100
    a1.grid = seq(0.501,0.999,length = ngrid)
    a2.grid = seq(0.9,0.999,length = 20)
    
    prod.lik=array(dim=c(ngrid,20))
    for(i in 1:ngrid){
      for(j in 1:20){
        tem=t[1]+t[2]*x
        p=(1-a2.grid[j])+(a1.grid[i]+a2.grid[j]-1)*exp(tem)/(1+exp(tem))
        #### NUMERICAL ISSUE, need to time a number (e.g., 1.4) to avoid -inf after multiplication
        #### As the no. of sample size increases and the no. of prevalence decreases, 
        #### the size of adjustment needed is larger.
        #### Practically, may need to try a few choices.
        likelihood=(adjust*p)^s*(adjust-adjust*p)^(1-s) 
        prod.lik[i,j] = prod(likelihood)*dprior(a1.grid[i],mu1,var1)*dprior(a2.grid[j],mu2,var2)
      }
    }
    return(-log(sum(prod.lik)))
  }
  
  est.int= nlminb(c(0,0),norm.int.lik,lower=c(-5,-5),upper=c(5,5))$par
  #cov.est.int = ginv(hessian(norm.int.lik, est.int))
  
  list(out = c(est.int)) #,cov.est.int[1,1], cov.est.int[2,2], cov.est.int[1,2]))
}

# prior=beta distribution
dprior.beta=function(x, a,b){
  dbeta(x,a,b)
}

PIE.beta <- function(s, x, adjust=NULL){
  n <- length(s)
  ## t = (beta0,beta1)
  ## PIE estimator of beta
  int.lik.beta=function(t)
  {
    ngrid=100
    a1.grid = seq(0.2,0.999,length = ngrid)
    a2.grid = seq(0.95,0.999,length = 10)
    
    prod.lik=array(dim=c(ngrid,10))
    for(i in 1:ngrid){
      for(j in 1:10){
        tem=t[1]+t[2]*x
        p=(1-a2.grid[j])+(a1.grid[i]+a2.grid[j]-1)*exp(tem)/(1+exp(tem))
        #### NUMERICAL ISSUE, need to time a number (e.g., 1.4) to avoid -inf after multiplication
        #### As the no. of sample size increases and the no. of prevalence decreases, 
        #### the size of adjustment needed is larger.
        #### Practically, may need to try a few choices.
        likelihood=(adjust*p)^s*(adjust-adjust*p)^(1-s) 
        prod.lik[i,j] = prod(likelihood)*dprior.beta(a1.grid[i],a,b)*dprior.beta(a2.grid[j],a2,b2)
      }
    }
    return(-log(sum(prod.lik)))
  }
  ## prior distributions for sensitivity and specificity both uniform
  est.int= nlminb(c(0,0),int.lik.beta,lower=c(-5,-5),upper=c(5,5))$par
  list(out = c(est.int))
}