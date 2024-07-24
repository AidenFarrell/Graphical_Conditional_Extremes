################################################################################
## Reading in required packages
required_pckgs <- c("Matrix", "sparseMVN")
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################

#Functions for the AGGD
#Parameterisation as given by https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8959693

## Univariate functions

## Probability density function
dagg <- function(x, loc, scale_1, scale_2, shape, log = FALSE){
  
  ## checks
  if(!is.numeric(x)){stop("x must be a vector")}
  if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("scale and shape must be a positive real number")
  }
  
  ## set up the parameters
  n <- length(q)
  z <- numeric(n)
  
  ## Get the density on the log scale 
  l_const <- log(shape) - log(scale_1 + scale_2) - lgamma(1/shape)
  z[x < loc] <- l_const - ((-x[x < loc] + loc)/(scale_1))^(shape)
  z[x >= loc] <- l_const - ((x[x >= loc] - loc)/(scale_2))^(shape)
  
  ## exponentiate the result if required
  if(!log){
    z <- exp(z)
  }
  return(z)
}

## Cumulative distribution function
pagg <- function(q, loc, scale_1, scale_2, shape, log = FALSE){
  
  ## checks
  if(!is.numeric(q)){stop("q must be a vector")}
  if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("scale and shape must be a positive real number")
  }
  
  ## set up the parameters
  n <- length(q)
  z <- numeric(n)
  
  ## get the result
  z[q < loc] <- (scale_1/(scale_1 + scale_2))*(1 - pgamma(q = ((loc - q[q < loc])/scale_1)^(shape), shape = 1/shape, rate = 1))
  z[q >= loc] <- (scale_1 + scale_2*pgamma(q = ((q[q >= loc] - loc)/scale_2)^(shape), shape = 1/shape, rate = 1))/(scale_1 + scale_2)  
  
  ## log the result if needed
  if(log){
    z <- log(z)
  }
  
  ## return the output
  return(z)
}

## Qunatile function
qagg <- function(p, loc, scale_1, scale_2, shape){
  if(!is.numeric(p)){stop("p must be a vector")}
  if(any(p < 0) | any(p > 1)){stop("p must be in the region [0,1]")}
  if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("scale and shape must be a positive real number")
  }
  
  n <- length(p)
  z <- numeric(n)
  
  z[p < (scale_1)/(scale_1 + scale_2)] <- loc - scale_1*qgamma(1 - p[p < (scale_1)/(scale_1 + scale_2)]*(scale_1 + scale_2)/scale_1, shape = 1/shape, rate = 1)^(1/shape)
  z[p >= (scale_1)/(scale_1 + scale_2)] <- loc + scale_2*qgamma((p[p >= (scale_1)/(scale_1 + scale_2)]*(scale_1 + scale_2) - scale_1)/scale_2, shape = 1/shape, rate = 1)^(1/shape) 
  
  return(z)
}

## Random generation
ragg <- function(n, loc, scale_1, scale_2, shape){
  if(length(n) > 1 | n <= 0 | n%%1 != 0){
    stop("n must be a single positive integer")
  }
  if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("scale must be a positive real number")
  }
  
  p <- runif(n)
  z <- qagg(p, loc, scale_1, scale_2, shape)
  return(z)
}

## log-likelihood for the model
llh_agg <- function(x, par, negative = FALSE){
  
  ## extract parameter estimates
  loc <- par[1]
  scale_1 <- par[2]
  scale_2 <- par[3]
  shape <- par[4]
  
  ## check parameters
  if(scale_1 <= 0 | scale_2 <= 0 | shape <= 0){
    return((-10^10)*(-1)^negative)
  }
  else{
    ## calculate log-likelihood
    z <- sum(dagg(x, loc, scale_1, scale_2, shape, log = TRUE))
    return(z*(-1)^negative)
  }
}

## Fit the log-likelihood for the model using optim
fit_agg <- function(par, data){
  ## checks
  if(!is.vector(data)){stop("data must be a vector")}
  if(length(par) != 4){stop("invalid number of parameters")}
  
  ## Fit the model
  fit <- try(optim(par = par, fn = llh_agg, x = data, negative = TRUE,
               control = list(maxit = 1e+9), method = "Nelder-Mead"),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from fit_agg")
    out <- list()
    out$par <- matrix(NA, nrow = 4, ncol = 1)
    rownames(out$par) <- c("loc", "scale_1", "scale_2", "shape")
    out$llh <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in fit_agg")
    out <- list()
    out$par <- matrix(NA, nrow = 4, ncol = 1)
    rownames(out$par) <- c("loc", "scale_1", "scale_2", "shape")
    out$llh <- NA
    out$convergence <- NA
  }
  else{
    ## Extract the output
    out <- list()
    out$par <- as.matrix(fit$par)
    rownames(out$par) <- c("loc", "scale_1", "scale_2", "shape")
    out$llh <- -fit$value
    out$convergence <- fit$convergence 
  }
  return(out)
}

################################################################################
## Functions for the multivariate AGG using a Gaussian Copula

## Probability density function
dmvagg <- function(data, loc, scale_1, scale_2, shape, Gamma, log = FALSE){
  
  ## check data is in correct format
  if(!is.matrix(data)){stop("data must be a matrix")}
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  ## checks on the inputs
  if(!is.numeric(loc) | !is.numeric(scale_1) | !is.numeric(scale_2) | !is.numeric(shape)){
    stop("parameters must be vectors")
  }
  else if(length(loc) != d | length(scale_1) != d | length(scale_2) != d | length(shape) != d){
    stop("parameters need to have the same dimension of data")
  }
  else if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("Invalid parameters provided")
  }
  if(!missing(Gamma)){
    if(!(class(Gamma) != "dgCMatrix" | class(Gamma) != "dsCMatrix")){
      stop("Gamma must be sparseMatrix")
    }
  }
  
  l_f_z <- sapply(1:d, function(i){dagg(x = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i], log = TRUE)})
  
  Q_F_z <- sapply(1:d, function(i){qnorm(pagg(q = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i]))})
  if(missing(Gamma)){
    warning("No precision matrix has been provided. The empirical precision matrix will be used")
    Gamma <- as(solve(cor(Q_F_z)), "sparseMatrix") 
  }
  l_mvnorm <- dmvn.sparse(x = Q_F_z, mu = rep(0, d), CH = Cholesky(Gamma), log = TRUE)
  l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
  
  z <- l_mvnorm + apply((l_f_z - l_dnorm), 1, sum)
  if(!log){
    z <- exp(z)
  }
  return(z)
}

## Random generation
rmvagg <- function(n, loc, scale_1, scale_2, shape, Gamma){
  #Checks on parameters
  if(!(class(Gamma) != "dgCMatrix" | class(Gamma) != "dsCMatrix")){
    stop("Gamma must be sparseMatrix")
  }
  else if(length(n) > 1 | n <= 0 | n%%1 != 0){
    stop("n must be a single positive integer")
  }
  else if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("Invalid parameters provided")
  }
  else{
    d <- dim(Gamma)[1]
    if(length(loc) == 1){loc <- rep(loc, d)}
    if(length(scale_1) == 1){scale_1 <- rep(scale_1, d)}
    if(length(scale_2) == 1){scale_2 <- rep(scale_2, d)}
    if(length(shape) == 1){shape <- rep(shape, d)}
    
    p <- apply(rmvn.sparse(n = n, mu = rep(0, d), CH = Cholesky(Gamma)), 1, pnorm)
    z <- sapply(1:d, function(i){qagg(p = p[i,], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i])})
    return(z)
  }
}