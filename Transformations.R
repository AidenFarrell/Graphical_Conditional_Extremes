################################################################################
## Reading in required packages
required_pckgs <- c("LaplacesDemon")
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in required scripts
source("threshold_selection_paper/thresh_qq_metric.R")

################################################################################
## Coles and Tawn (1991) semi-parametric approach for transformation to uniform margins
semiparametric_unif_transform <- function(data, par){
  if(!is.numeric(data)){
    stop("data must be a vector")
  }
  if(!is.numeric(par) | length(par) != 3){
    stop("par must be a vector of length three")
  }
  ## Obtain infromation from the data
  n <- length(data)
  
  ## Extract parameters
  u <- par[1]
  scale <- par[2]
  shape <- par[3]
  
  ## Index of excesses
  index <- which(data > u)
  
  ## Below the threshold we use the empirical CDF
  z <- rank(data)/(n + 1)
  
  ## Above the threhsold we assume a GPD tail
  pu <- length(index)/n
  if(abs(shape) < 1e-6){
    z[index] <- 1 - pu*exp(-exp(-(data[index] - u)/scale))
  }else{
    z[index] <- 1 - pu*pmax(0, 1 + shape*(data[index] - u)/scale)^(-1/shape)
  }
  return(z)
}

## Coles and Tawn (1991) semi-parametric approach for transformation from uniform
## margins to original data margins
semiparametric_gpd_transform <- function(p, par, distFun){
  if(!is.numeric(p) | any(p <= 0) | any(p >= 1)){
    stop("p must be a vector of probabilties in the region(0,1)")
  }
  if(!any(class(distFun) == "ecdf")){
    stop("distFun must have class ecdf")
  }
  
  # extract MLEs
  u <- par[1]
  scale <- par[2]
  shape <- par[3]
  qu <- par[4]
  
  ## Index of excesses
  index <- which(p > qu)
  
  ## Below the threshold
  z <- quantile(distFun, p)
  ## Above the threshold
  if(abs(shape) <= 1e-6){
    z[index] <- u + scale*log((1 - qu)/(1 - p[index]))
  }
  else{
    z[index] <- u + (scale/shape)*(((1 - qu)/(1 - p[index]))^(shape) - 1)
  }
  return(z)
}


## Transform onto standard Laplace margins
## Uses the semi-parametric method of Coles & Tawn (1991)
## Threshold is selected using Murphy et al. (2024)

X_to_Laplace <- function(x, q = seq(0.55, 0.99, by = 0.01), k = 200, m = 500){
  ## Checks on the inputs
  if(!is.numeric(x)){
    stop("x must be a vector")
  }
  if(!is.numeric(q)){
    stop("q must be a vector of quantiles to be tested")
  }
  if(min(q) <= 0.5 | max(q) > 0.999){
    stop("q must be between 0.5 and 0.999 exclusive")
  }
  if(k <= 0 | k %% 1 != 0){
    stop("Number of bootstrapped samples must be a positive integer")
  }
  if(m <= 0 | m %% 1 != 0){
    stop("Number of equally spaced probabilities must be a positive integer")
  }
  
  ## Perform the threshold selection
  u_poss <- quantile(x, q)
  u_out <- thresh_qq_metric(data = x, thresh = u_poss, k = k, m = m)
  
  ## Get the output from the threshold selection
  u_star <- u_out$thresh
  qu_star <- q[which(u_poss == u_star)]
  scale_star <- u_out$par[1]
  shape_star <- u_out$par[2]
  
  ## Transform the data onto Laplace margins
  y <- plaplace(semiparametric_unif_transform(data = x, u = u_star, par = par_star))
  
  ## Output
  out <- list(data = list(X = x, Y = y),
              par = list(u = u_star,
                         qu = qu_star,
                         scale = scale_star,
                         shape = shape_star))
  class(out) <- "Laplace_Transform"
  return(out)
}