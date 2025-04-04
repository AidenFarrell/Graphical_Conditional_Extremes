# CDF, quantile function, Jacobian for transformation to delta-Laplace margins.

#------------------------------------------------------------------------------------------------------------------
# pdlaplace: univariate distribution function


pdlaplace<-function(z,mu,sigma,delta,lower.T=T)
{
  
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  nz<-length(z)
  result<-numeric(nz)
  if(length(delta)==1){delta<-rep(delta,nz)}
  if(lower.T==T){
    result[z<mu]<- 0.5*pgamma((((mu-z)/sigma)[z<mu])^delta[z<mu],shape=1/delta[z<mu],scale=1,lower.tail = F)
  result[z>=mu]<- 0.5+0.5*pgamma((((z-mu)/sigma)[z>=mu])^delta[z>=mu],shape=1/delta[z>=mu],scale=1)
  }else{
    
    result[z<mu]<- 1-0.5*pgamma((((mu-z)/sigma)[z<mu])^delta[z<mu],shape=1/delta[z<mu],scale=1,lower.tail = F)
    result[z>=mu]<- 0.5*pgamma((((z-mu)/sigma)[z>=mu])^delta[z>=mu],shape=1/delta[z>=mu],scale=1,lower.tail=F)
  }
 
  return(result)
}

#------------------------------------------------------------------------------------------------------------------
# qdlaplace: univariate quantile function

qdlaplace<-function(q,mu,sigma,delta)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  nq<-length(q)
  result<-numeric(nq)
  if(length(mu)==1){mu<-rep(mu,nq)}
  if(length(sigma)==1){sigma<-rep(sigma,nq)}
  if(length(delta)==1){delta<-rep(delta,nq)}
  result[q<0.5]<- mu[q<0.5] - sigma[q<0.5]*qgamma(1-2*q[q<0.5],shape=1/delta[q<0.5],scale=1)^(1/delta[q<0.5])
  result[q>=0.5]<- mu[q>=0.5] + sigma[q>=0.5]*qgamma(2*q[q>=0.5]-1,shape=1/delta[q>=0.5],scale=1)^(1/delta[q>=0.5])
  return(result)
}


#------------------------------------------------------------------------------------------------------------------
# rdlaplace: univariate random number generation

rdlaplace<-function(n,mu,sigma,delta)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  return(qdlaplace(runif(n),mu=mu,delta=delta,sigma=sigma))
}


#------------------------------------------------------------------------------------------------------------------
# ddlaplace: univariate density function

ddlaplace<-function(z,mu,sigma,delta,log=FALSE)
{
  k=sqrt(gamma(1/delta)/gamma(3/delta))
  sigma=k*sigma
  ld<- -abs((z-mu)/(sigma))^delta +log(delta)-log(2*sigma)-lgamma(1/delta)
  if(log==TRUE){return(ld)}
  else{return(exp(ld))}
}

#------------------------------------------------------------------------------------------------------------------
# dmvdlaplace: multivariate density function with Gaussian copula


# with faster MVN computation: requires matrix Sigma AND its cholesky factorization, chol(Sigma), which results in faster
# computation when used in apply()

dmvdlaplace<-function(z,mu,sigmad,SigmaChol,Sigma,delta,log=FALSE)
{
  z_list <- lapply(apply(z, 2, list), function(x){x[[1]]})
  un <- pmap(.l = list(z = z_list, mu = mu, sigma = sigmad, delta = delta), .f = pdlaplace, lower.T = F)
  zn <- pmap(.l = list(p = un, mean = mu, sd = sqrt(diag(Sigma))), .f = qnorm, lower.tail = F)
  # zn<-qnorm(pdlaplace(z,mu=mu,sigma=sigmad,delta=delta,lower.T=F),mean=mu,sd=sqrt(diag(Sigma)),lower.tail =F)
  ld1<-sum(dmvn(do.call(cbind, zn),mu=mu,Sigma = Sigma,log=TRUE))
  ld2<-sum(do.call(c, pmap(.l = list(z = z_list, mu = mu, sigma = sigmad, delta = delta), .f = ddlaplace, log = TRUE)))
  ld3<-sum(do.call(c, pmap(.l = list(x = zn, mean = mu, sd = sqrt(diag(Sigma))), .f = dnorm, log = TRUE)))
  if(log==TRUE){return(ld1+ld2-ld3)}
  else{return(exp(ld1+ld2-ld3))}
}

#------------------------------------------------------------------------------------------------------------------
# rmvdlaplace: multivariate random number generation with Gaussian copula

rmvdlaplace<-function(n,dim,mu,sigmad,Sigma,delta)
{
  if(dim(Sigma)[1]!=dim){stop("dim and Sigma do not match")}
  NU<-t(apply(rmvn(n,mu=rep(0,length=dim),Sigma=Sigma),1,pnorm,sd=sqrt(diag(Sigma))))
  return(t(apply(NU,1,qdlaplace,mu=mu,sigma=sigmad,delta=delta)))
}

#------------------------------------------------------------------------------------------------------------------
## Function to evaluate the log-likelihood of the delta laplace distribution for fitting
llh_ddlaplce <- function(par, z, negative = FALSE){
  ## Extract parameters
  mu <- par[1]
  sigma <- par[2]
  delta <- par[3]
  
  if(sigma <= 0 | delta <= 0){
    return(-10^10*((-1)^negative))
  }
  
  llh <- sum(ddlaplace(z, mu, sigma, delta, log = TRUE))
  
  if(is.infinite(llh)){
    return(-10^10*((-1)^negative))
  }
  
  return(llh*((-1)^negative))
}


## Fit the log-likelihood for the model using optim
fit_ddlaplce <- function(par, data){
  ## checks
  if(!is.numeric(data)){stop("data must be a vector")}
  if(length(par) != 3){stop("invalid number of parameters")}
  
  ## Fit the model
  fit <- try(optim(par = par, fn = llh_ddlaplce, z = data, negative = TRUE,
                   control = list(maxit = 1e+9), method = "Nelder-Mead"),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from fit_agg")
    out <- list()
    out$par <- matrix(NA, nrow = 3, ncol = 1)
    rownames(out$par) <- c("loc", "scale", "shape")
    out$llh <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in fit_agg")
    out <- list()
    out$par <- matrix(NA, nrow = 3, ncol = 1)
    rownames(out$par) <- c("loc", "scale", "shape")
    out$llh <- NA
    out$convergence <- NA
  }
  else{
    ## Extract the output
    out <- list()
    out$par <- as.matrix(fit$par)
    rownames(out$par) <- c("loc", "scale", "shape")
    out$llh <- -fit$value
    out$convergence <- fit$convergence 
  }
  return(out)
}
