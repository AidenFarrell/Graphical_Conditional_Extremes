#read in required file
required_pckgs <- c("parallel", "purrr")
lapply(required_pckgs, require, character.only = TRUE)
source("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/Threshold_Selection/_gpd.R")

################################################################################
#empirical cdf function
EmpCDF <- function(x, denom = c("n", "n+1")){
  if(!is.null(dim(x))){
    stop("x must be a one dimensional vector")
  }
  denom <- match.arg(denom)
  #length of original data
  n <- length(x)
  if(denom == "n+1"){
    n <- n+1
  }
  
  #Calculate probabilities
  counts <- as.vector(table(x))
  probs <- cumsum(counts)/n
  
  #output
  return(cbind(sort(unique(x)), probs))
  
}

loglike_gpd <- function(data, u, par = c(10, 0.1), negative = FALSE){
  
  #make sure all data is above the threshold
  data <- data[data > u]
  #set parameters
  sig_u <- par[1]
  xi  <- par[2]
  
  #check all parameter values
  if(any(sig_u <= 0)){
    return(-10^10*(-1)^negative)
  }
  # check all xi below UEP (if it exists, adjusting for rounding)
  else if(xi < 0 & max(data) >= u - sig_u/xi){
    return(-10^10*(-1)^negative)
  }
  else{
    if(abs(xi) > 1e-06 & any(1 + xi*(data - u)/sig_u <= 0)) {
      return(-10^10*(-1)^negative)
    }
    else if(abs(xi) > 1e-06 & any(1 + xi*(data - u - 1)/sig_u <= 0)) {
      return(-10^10*(-1)^negative)
    }
    else{
      z <- sum(dgpd(x = data, shape = xi, scale = sig_u, mu = u, log = TRUE))
    }
  }
  return(z*(-1)^negative)
}

mle_gpd <- function(data, u, par = c(10, 0.1), llh_val = TRUE, hessian = FALSE){
  #set parameters
  if (length(par) != 2) {
    stop("invalid number of paramters provided")
  }
  sig_init <- par[1]
  xi_init <- par[2]
  
  #make sure no data is below the threhsold
  data <- data[data > u]
  
  #Check parameter values
  stopifnot(all(sig_init > 0))
  
  if((xi_init < 0) & any((u - sig_init/xi_init) <= max(data))){
    stop("Invalid starting point. Data above upper end point of distribution.")
  }
  
  #Calculate the negative log-likelihood
  optout <- optim(par = as.vector(par), fn = loglike_gpd,
                  data = data, u = u, negative = TRUE, 
                  method = "Nelder-Mead", hessian = hessian)
  
  if(!llh_val & !hessian){
    out <- optout$par
  }
  else{
    out <- list(mle = optout$par)
  }
  if(llh_val){
    out$nllh <- optout$value
  }
  if(hessian){
    out$hessian <- optout$hessian
  }
  return(out)
}

dq1_metric <- function(x, y){
  mean(abs(x - y))
}

#filter data
filter_data <- function(data, thresholds){
  data %>%
    mutate(v = thresholds) %>% #mutates the data with the thresholds
    filter(y > v) #filters such that the rounded data is greater than the threshold
}

################################################################################
#Try and parallise the above code
boot_mles_and_dist <- function(data, u, init_mles, start_par, p, B){
  #create a list of non-parametrically bootstrapped samples
  boot_samples <- replicate(n = B, 
                            expr = sample(x = data, size = length(data), replace = TRUE),
                            simplify = FALSE)
  
  #initialise starting parameters for bootstrapping
  if(init_mles[2] <= 0){
    par_boot <- rbind(sapply(X = boot_samples, FUN = mean), rep(start_par[2], times = B))
    #convert above to a list 
    par_boot <- lapply(X = seq_len(ncol(par_boot)), FUN = function(j){par_boot[,j]})
  }
  else{
    par_boot <- rep(list(start_par), B)
  }
  
  #get the mles of the bootstrapped samples
  boot_mles <- pmap(.f = mle_gpd,
                    .l = list(data = boot_samples,
                              par = par_boot),
                    u = u,
                    llh_val = FALSE)
  
  #extract mles as vectors
  sig_boot_mles <- sapply(X = boot_mles, FUN = function(vec){vec[1]})
  xi_boot_mles <- sapply(X = boot_mles, FUN = function(vec){vec[2]})
  
  #calculate the theoretical quantiles using the bootstrapped MLES
  tq <- pmap(.f = qgpd,
             .l = list(scale = sig_boot_mles, shape = xi_boot_mles),
             p = p,
             mu = u)
  
  #calculate the empirical quantiles of the bootstrapped samples using the asymmetrically spaced probabilities
  #Need to calculate the percentiles that we are going to test at
  eq <- pmap(.f = quantile,
             .l = list(x = boot_samples),
             probs = p)
  
  #obtain the distance metrics
  dq1 <- sapply(X = pmap(.l = list(x = eq, y = tq), .f = dq1_metric), FUN = mean, na.rm = TRUE)
  return(dq1)
}

thrs_distance_parallel <- function(data, u, par = c(10, 0.1), B = 500, m = 501, l = 5){
  #checks on the inputs
  if(length(dim(data)) != 0){
    stop("Data must be a vector")
  }
  if(length(dim(u)) != 0){
    stop("u to be tested needs to be a vector needs to be a vector")
  }
  if(length(par) != 2){
    stop("Incorrect numebr of parameters to be estimated provided")
  }
  if(B <= 0 | B%%1 != 0){
    stop("Number of bootstrapped samples must be a positive integer")
  }
  if(m <= 0 | m%%1 != 0){
    stop("Number of equally spaced probabilities must be a positive integer")
  }
  
  #ensure l number of unique values in the bulk and tail
  nu_au <- sapply(X = u, FUN = function(x, u){length(unique(x[x > u]))}, x = data)
  if(all(nu_au < l)){
    stop("Not enough unique exceedances in the tail")
  }
  if(any(nu_au < l)){
    index <- which(nu_au < l)
    u <- u[-index]
  }
  
  #fit the models to get the point estimates
  data_au <- lapply(X = u, FUN = function(x, y){sort(y[y > x])}, y = data)
  init_mles <- mcmapply(FUN = mle_gpd,
                        data = data_au,
                        u = u,
                        MoreArgs = list(par = as.vector(par),
                                        llh_val = FALSE),
                        mc.cores = detectCores(),
                        SIMPLIFY = FALSE)
  
  #setting up the equally spaced probabilities
  p <- ppoints(m)
  
  #calculate the distance metrics
  dists <- mcmapply(FUN = boot_mles_and_dist,
                    u = u,
                    init_mles = init_mles,
                    data = data_au,
                    MoreArgs = list(p = p,
                                    B = B,
                                    start_par = par),
                    mc.cores = detectCores(),
                    SIMPLIFY = FALSE)
  dists <- sapply(X = dists, FUN = mean, na.rm = TRUE)
  
  #obtain the choice
  ind <- which.min(dists)
  choice <- u[ind]
  xi <- init_mles[[ind]][2]
  sig <- init_mles[[ind]][1]
  out <- list(u = choice,
              par = c(sig, xi),
              all_dists = data.frame(u = u, dq1 = dists))
  return(out)
}