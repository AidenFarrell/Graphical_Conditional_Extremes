#Functions for the AGGD
#Parameterisation as given by https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8959693
daggd <- function(x, loc, scale_1, scale_2, shape, log = FALSE){
  
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

paggd <- function(q, loc, scale_1, scale_2, shape, log = FALSE){
  
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

qaggd <- function(p, loc, scale_1, scale_2, shape){
  if(!is.null(dim(p))){
    stop("q must be a vector")
  }
  if(any(p < 0) | any(p > 1)){
    stop("p must be in the region [0,1]")
  }
  if(any(scale_1 <= 0) | any(scale_2 <= 0)){
    stop("scale must be a positive real number")
  }
  
  n <- length(p)
  z <- numeric(n)
  
  z[p < (scale_1)/(scale_1 + scale_2)] <- loc - scale_1*qgamma(1 - p[p < (scale_1)/(scale_1 + scale_2)]*(scale_1 + scale_2)/scale_1, shape = 1/shape, rate = 1)^(1/shape)
  z[p >= (scale_1)/(scale_1 + scale_2)] <- loc + scale_2*qgamma((p[p >= (scale_1)/(scale_1 + scale_2)]*(scale_1 + scale_2) - scale_1)/scale_2, shape = 1/shape, rate = 1)^(1/shape) 
  
  return(z)
}

raggd <- function(n, loc, scale_1, scale_2, shape){
  if(length(n) > 1 | n%%1 != 0){
    stop("n must be a single positive integer")
  }
  if(any(scale_1 <= 0) | any(scale_2 <= 0)){
    stop("scale must be a positive real number")
  }
  
  p <- runif(n)
  z <- qaggd(p, loc, scale_1, scale_2, shape)
  return(z)
}

llh_aggd <- function(x, par, negative = FALSE){
  
  loc <- par[1]
  scale_1 <- par[2]
  scale_2 <- par[3]
  shape <- par[4]
  
  if(scale_1 <= 0 | scale_2 <= 0 | shape <= 0){
    return((-10^10)*(-1)^negative)
  }
  else{
    z <- sum(daggd(x, loc, scale_1, scale_2, shape, log = TRUE))
    return(z*(-1)^negative)
  }
}

fit_aggd <- function(par, data){
  fit <- optim(par = par, fn = llh_aggd, x = data, negative = TRUE,
               control = list(maxit = 1e+9), method = "Nelder-Mead")
  
  out <- list()
  out$par$main <- as.matrix(fit$par)
  rownames(out$par$main) <- c("loc", "scale_1", "scale_2", "shape")
  
  out$llh <- -fit$value
  out$convergence <- fit$convergence
  
  return(out)
}

mean_aggd <- function(loc, scale_1, scale_2, shape){
  if(any(c(scale_1, scale_2, shape) <= 0)){
    stop("Invalid parameters provided")
  }
  z <- loc - (scale_1 - scale_2)*exp(lgamma(2/shape) - lgamma(1/shape))
  return(z)
}

sd_aggd <- function(loc, scale_1, scale_2, shape){
  if(any(c(scale_1, scale_2, shape) <= 0)){
    stop("Invalid parameters provided")
  }
  z <- (scale_1^3 + scale_2^3)/(scale_1 + scale_2)*exp(lgamma(3/shape) - lgamma(1/shape)) - 
    ((scale_1 - scale_2)*exp(lgamma(2/shape) - lgamma(1/shape)))^2
  return(sqrt(z))
}

## Functions for the multivariate AGG using a Gaussian Copula
dmvagg <- function(data, loc, scale_1, scale_2, shape, log = FALSE){
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  if(length(loc) != d | length(scale_1) != d | length(scale_2) != d | length(shape) != d){
    stop("parameters need to have the same dimension of data")
  }
  else if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("Invalid parameters provided")
  }
  else{
    l_f_z <- sapply(1:d, function(i){daggd(x = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i], log = TRUE)})
    
    Q_F_z <- sapply(1:d, function(i){qnorm(paggd(q = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i]))})
    l_mvnorm <- mvtnorm::dmvnorm(x = Q_F_z, mean = rep(0, d), sigma = cor(Q_F_z), log = TRUE)
    l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
    
    z <- l_mvnorm + apply((l_f_z - l_dnorm), 1, sum)
    if(!log){
      z <- exp(z)
    }
    return(z)
  }
}

rmvagg <- function(n, dim, loc, scale_1, scale_2, shape, Sigma){
  #Checks on parameters
  if(dim(Sigma)[1] != dim){
    stop("dim and Sigma do not match")
  }
  else if(dim(Sigma)[1] != dim | dim(Sigma)[2] != dim){
    stop("Sigma is not a dxd matrix")
  }
  else if(any(diag(Sigma) < 0) | min(eigen(Sigma)$values) <= 0 | !isTRUE(all.equal(Sigma, t(Sigma)))){
    stop("Sigma is not a valid covaraince matrix")
  }
  else if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("Invalid parameters provided")
  }
  else{
    if(length(loc) == 1){loc <- rep(loc, dim)}
    if(length(scale_1) == 1){scale_1 <- rep(scale_1, dim)}
    if(length(scale_2) == 1){scale_2 <- rep(scale_2, dim)}
    if(length(shape) == 1){shape <- rep(shape, dim)}
    
    p <- t(apply(mvtnorm::rmvnorm(n = n, mean = rep(0, dim), sigma = Sigma), 1, pnorm))
    z <- sapply(1:dim, function(i){qaggd(p = p[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i])})
    return(z)
  }
}

llh_mvagg <- function(par, data, negative = FALSE){
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  #extract the parameters
  loc <- par[1:d]
  scale_1 <- par[(d + 1):(2*d)]
  scale_2 <- par[(2*d + 1):(3*d)]
  shape <- par[(3*d + 1):(4*d)]
  
  #check the parameter values
  if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    return((-10^10)*(-1)^negative)
  }
  else{
    llh <- tryCatch(sum(dmvagg(data = data, loc = loc, scale_1 = scale_1, scale_2 = scale_2, shape = shape, log = TRUE)),
                    error = function(e){-10^10})
    return(llh*(-1)^negative) 
  }
}

fit_mvagg <- function(par, data){
  fit <- optim(par = par, fn = llh_mvagg, data = data,
               negative = TRUE, hessian = FALSE,
               method = "BFGS", control = list(maxit = 1e+6))
  
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  #extract the MLEs
  loc <- fit$par[1:d]
  scale_1 <- fit$par[(d + 1):(2*d)]
  scale_2 <- fit$par[(2*d + 1):(3*d)]
  shape <- fit$par[(3*d + 1):(4*d)]
  
  #Get the correlation matrix from the model
  Q_F_Z <- sapply(1:d, function(i){
    qnorm(paggd(q = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i]))
  })
  Sigma <- cor(Q_F_Z)
  
  #Make the output
  out <- list()
  out$par$main <- rbind(loc, scale_1, scale_2, shape)
  out$par$Sigma <- Sigma
  rownames(out$par$main) <- c("loc", "scale_1", "scale_2", "shape")
  
  out$llh <- -fit$value
  out$convergence <- fit$convergence
  return(out)
}

## above but assume everything is independent
dmvagg_indep <- function(data, loc, scale_1, scale_2, shape, log = FALSE){
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  if(length(loc) != d | length(scale_1) != d | length(scale_2) != d | length(shape) != d){
    stop("parameters need to have the same dimension of data")
  }
  else if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    stop("Invalid parameters provided")
  }
  else{
    l_f_z <- sapply(1:d, function(i){daggd(x = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i], log = TRUE)})
    
    Q_F_z <- sapply(1:d, function(i){qnorm(paggd(q = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i]))})
    l_mvnorm <- mvtnorm::dmvnorm(x = Q_F_z, mean = rep(0, d), sigma = diag(1, d), log = TRUE)
    l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
    
    z <- l_mvnorm + apply((l_f_z - l_dnorm), 1, sum)
    if(!log){
      z <- exp(z)
    }
    return(z)
  }
}

llh_mvagg_indep <- function(par, data, negative = FALSE){
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  #extract the parameters
  loc <- par[1:d]
  scale_1 <- par[(d + 1):(2*d)]
  scale_2 <- par[(2*d + 1):(3*d)]
  shape <- par[(3*d + 1):(4*d)]
  
  #check the parameter values
  if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    return((-10^10)*(-1)^negative)
  }
  else{
    llh <- tryCatch(sum(dmvagg_indep(data = data, loc = loc, scale_1 = scale_1, scale_2 = scale_2, shape = shape, log = TRUE)),
                    error = function(e){-10^10})
    return(llh*(-1)^negative) 
  }
}

fit_mvagg_indep <- function(par, data){
  fit <- optim(par = par, fn = llh_mvagg_indep, data = data,
               negative = TRUE, hessian = FALSE,
               method = "BFGS", control = list(maxit = 1e+6))
  
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  #extract the MLEs
  loc <- fit$par[1:d]
  scale_1 <- fit$par[(d + 1):(2*d)]
  scale_2 <- fit$par[(2*d + 1):(3*d)]
  shape <- fit$par[(3*d + 1):(4*d)]
  
  Sigma <- diag(1, d)
  
  #Make the output
  out <- list()
  out$par$main <- rbind(loc, scale_1, scale_2, shape)
  out$par$Sigma <- Sigma
  rownames(out$par$main) <- c("loc", "scale_1", "scale_2", "shape")
  
  out$llh <- -fit$value
  out$convergence <- fit$convergence
  return(out)
}

################################################################################
llh_mvagg_graph <- function(par, data, cliques, seps, negative = FALSE){
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  #extract the parameters
  loc <- par[1:d]
  scale_1 <- par[(d + 1):(2*d)]
  scale_2 <- par[(2*d + 1):(3*d)]
  shape <- par[(3*d + 1):(4*d)]
  
  #check the parameter values
  if(any(scale_1 <= 0) | any(scale_2 <= 0) | any(shape <= 0)){
    return((-10^10)*(-1)^negative)
  }
  else{
    ## Transform data onto joint standard Gaussian margins
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i]))}))
    
    ## correlation matrix to be used
    Sigma <- as.matrix(cor(Q_F_z))
    
    l_f_z <- sapply(1:d, function(i){daggd(x = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i], log = TRUE)})
    l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
    
    l_mvnorm_cliques <- sapply(cliques, function(i){sum(mvtnorm::dmvnorm(x = as.matrix(Q_F_z[,i]), mean = rep(0, length(i)), sigma = Sigma[i,i], log = TRUE))})
    l_mvnorm_seps  <- sapply(seps, function(i){sum(mvtnorm::dmvnorm(x = as.matrix(Q_F_z[,i]), mean = rep(0, length(i)), sigma = as.matrix(Sigma[i,i]), log = TRUE))})
    
    llh_cliques <- sum(l_mvnorm_cliques) + sum(sapply(cliques, function(i){sum((l_f_z - l_dnorm)[,i])}))
    llh_seps <- sum(l_mvnorm_seps) + sum(sapply(seps, function(i){sum((l_f_z - l_dnorm)[,i])}))
    llh <- llh_cliques - llh_seps
    return(llh*(-1)^negative) 
  }
}

fit_mvagg_graph <- function(par, data, cliques, seps){
  
  ## if there are no separators then the cliques should not be connected
  ## and so they should be treated as independent
  
  if(is_empty(seps)){
    #extract info from data
    n <- dim(data)[1]
    d <- dim(data)[2]
    
    #extract the MLEs
    loc <- par[1:d]
    scale_1 <- par[(d + 1):(2*d)]
    scale_2 <- par[(2*d + 1):(3*d)]
    shape <- par[(3*d + 1):(4*d)]
    
    ## fit the model
    n_cliques <- length(cliques)
    fit <- vector("list", length = n_cliques)
    for(i in 1:n_cliques){
      index <- cliques[[i]]
      fit[[i]] <- fit_mvagg(data = data[,index],
                            par = c(loc[index], scale_1[index], scale_2[index], shape[index]))
    }
    
    loc <- do.call(c, lapply(fit, function(x){x$par$main[1,]}))
    scale_1 <- do.call(c, lapply(fit, function(x){x$par$main[2,]}))
    scale_2 <- do.call(c, lapply(fit, function(x){x$par$main[3,]}))
    shape <- do.call(c, lapply(fit, function(x){x$par$main[4,]}))
    
    Sigma <- matrix(data = 0, nrow = d, ncol = d)
    for(i in 1:n_cliques){
      Sig_elements <- permutations(n = length(cliques[[i]]), r = 2, v = cliques[[i]], repeats.allowed = TRUE)
      Sigma[Sig_elements] <- fit[[i]]$par$Sigma
    }
    
    #Make the output
    out <- list()
    out$par$main <- rbind(loc, scale_1, scale_2, shape)
    out$par$Sigma <- Sigma
    rownames(out$par$main) <- c("loc", "scale_1", "scale_2", "shape")
    
    out$llh <- sapply(fit, function(x){x$llh})
    out$convergence <- sapply(fit, function(x){x$convergence})
    return(out)
  }
  else{
    fit <- optim(par = par, fn = llh_mvagg_graph, data = data,
                 cliques = cliques, seps = seps,
                 negative = TRUE, hessian = FALSE,
                 method = "BFGS", control = list(maxit = 1e+9))
    
    #extract info from data
    n <- dim(data)[1]
    d <- dim(data)[2]
    
    #extract the MLEs
    loc <- fit$par[1:d]
    scale_1 <- fit$par[(d + 1):(2*d)]
    scale_2 <- fit$par[(2*d + 1):(3*d)]
    shape <- fit$par[(3*d + 1):(4*d)]
    
    #Get the correlation matrix from the model
    Q_F_Z <- sapply(1:d, function(i){
      qnorm(paggd(q = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i]))
    })
    Sigma <- as.matrix(cor(Q_F_Z))
    
    #Make the output
    out <- list()
    out$par$main <- rbind(loc, scale_1, scale_2, shape)
    out$par$Sigma <- Sigma
    rownames(out$par$main) <- c("loc", "scale_1", "scale_2", "shape")
    
    out$llh <- -fit$value
    out$convergence <- fit$convergence
    return(out) 
  }
}

################################################################################
## Update the univaraite AGG functions to be in line with this paper
## https://pdf.sciencedirectassets.com/271506/1-s2.0-S0957417422X00134/1-s2.0-S0957417422008429/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEMj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJGMEQCIAM%2FHP5DWlj3aAZkL4PVCbEVGSlVnjr9n6Jw4P%2FRxaeOAiB1SB2Uy4teIqj1dryxC7TeHFYZdzIkBDFLQkNQFTFVTyqyBQgREAUaDDA1OTAwMzU0Njg2NSIMa5CeYyb5EGXdMbsoKo8Fiovwg3We%2F2T5XkjaShuhNtXI7PU1N6n9Vrmy1VZ8%2BBoZyhgm%2BSRxG3runqmJHUd%2FS37xh5%2FVp38RRe4niPOt3cCRfNBuwIaYkj6xUbU4AIuiiV1socKs7sILfiSoPXCqtvUQ%2B7iQssvISNj6VyHGsJrjZGOaoss1qVxHfMhg6qPoxTj7aRCHvocgOsQb2h9yvMPso1KWOLPFHnINEjaHL9dKsL3OHdR5kqUjQxIwoMqLUP9oPoa8MX5Gx7LSVCguGWdiFRI%2BmPNsBEQbKIiJSnbeetzkgqtJLj6oE9dAzdCtNC%2FXaCCNKlT%2BCE1TI8KvvzCSgRhOKcfw3Q%2FMgU922sgYVGvwJhfklSVAWFmdL7z%2BGC1FqOmKrigyD5pAbtLY%2BShPYOp%2B1l3UMxPdRbyqKHYHbHU43G5j%2BcTgloUd9rtnzqCLgnvwTtP9sn%2FfmFf0%2FCw048BqV4hdEnA9%2F1SjVdn%2BFnE8BrjcITOR6rSjTqz728%2FsYKn5IsDKRaHBGveukAyOmTvd88AXr6zUI%2Fj9bkSYvF2QKEymZZJlusDN2JX9MhB%2BEwboKOPfkX0tSeFjw0sNSnc9%2BU%2FCUvHY47Es%2FFdY88fST9RSQ9hCARuVdvTIWuPdkfr2QYf%2FYLq4D1dDtKEcaZan0Pzea98b%2B4ru7sUxFLdM2bDhysiXb9uUhi52YFxh5xsExe2Auc76dPhcwfV%2Fk5tUGsPVzJ%2BaOaE%2BzUbPuSOi%2F9%2Bgs%2BeJk5qC9eHStGG7D81555jDHY3ysZDk1V4Zb8NmdhRMQfygj7nHmy%2Fayg6IIwaS4%2FwZeqmD%2FG2sRXUDKduGlvYK%2Bb9oEMvJnLFRbURCTAaYHoUiFs%2FzimzHkruxZ%2F0sKmHuzjiPKjDUyYixBjqyASdnKDwA7MinCqvnO3G6x3Xp17bWwFR3%2BjRQScY8gGERyKAXAyS%2FWeqEhPB1HqwwlKChHYfx01aDVas2gyBd7UeD5g0n%2BAYdUeaW%2BTGRttEQA2r1NIGb80F%2FtePSDYU7pobtXhSWqXkoMFRDbg237xCe0Rz%2BuLIW8U%2BFsYXigI9mInmHUdJJoDvvvaYibsEDrsd92Rb9RNEvdYnwDdC0zyvLsdWiHq3%2F0NfzbHx6w5w%2BJAo%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240419T092127Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYYFR7YFGS%2F20240419%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=90074fc66fcfd0b8e19e81fc23c7e6cd20dbf3a1140a5123c641c04ae2131f48&hash=238f6152e3ac556c7144dfa8c2c5faa0b2eb03c750525ffa79a540a75ad1d516&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0957417422008429&tid=spdf-31763a1c-2f71-497d-8f09-905c89f6d9e1&sid=cfc82d59929d494cb25b52a50c3067fb6766gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0104585657030b5d5556&rr=876bd3120c1e0747&cc=gb

# daggd <- function(x, loc, scale_1, scale_2, shape, log = FALSE){
#   if(!is.null(dim(x))){
#     stop("x must be a vector")
#   }
#   if(any(scale_1 <= 0) | any(scale_2 <= 0)){
#     stop("scale must be a positive real number")
#   }
#   
#   n <- length(x)
#   z <- rep(NA, n)
#   index_low <- which(x < loc)
#   
#   l_A_shape <- (shape/2)*(lgamma(3/shape) - lgamma(1/shape))
#   l_const <- log(shape) + l_A_shape/shape - log(scale_1 + scale_2) - lgamma(1/shape)
#   
#   index_low <- which(x < loc)
#   z[index_low] <- l_const - exp(l_A_shape)*(((-x[index_low] + loc)/(scale_1))^(shape))
#   z[-index_low] <- l_const - exp(l_A_shape)*(((x[-index_low] - loc)/(scale_2))^(shape))
#   if(!log){
#     z <- exp(z)
#   }
#   return(z)
# }
# 
# paggd <- function(q, loc, scale_1, scale_2, shape, log = FALSE){
#   if(!is.null(dim(q))){
#     stop("q must be a vector")
#   }
#   if(any(scale_1 <= 0) | any(scale_2 <= 0)){
#     stop("scale must be a positive real number")
#   }
#   
#   l_A_shape <- (shape/2)*(lgamma(3/shape) - lgamma(1/shape))
#   n <- length(q)
#   if(n == 1){
#     if(q < loc){
#       z <- (scale_1/(scale_1 + scale_2))*(1 - pgamma(q = exp(l_A_shape)*(((loc - q)/scale_1)^(shape)), shape = 1/shape, rate = 1))
#     }
#     else{
#       z <- (scale_1 + scale_2*pgamma(q = exp(l_A_shape)*(((q - loc)/scale_2)^(shape)), shape = 1/shape, rate = 1))/(scale_1 + scale_2) 
#     }
#   }
#   else{
#     z <- rep(NA, n)
#     index_low <- which(q < loc)
#     
#     z[index_low] <- (scale_1/(scale_1 + scale_2))*(1 - pgamma(q = exp(l_A_shape)*(((loc - q[index_low])/scale_1)^(shape)), shape = 1/shape, rate = 1))
#     z[-index_low] <- (scale_1 + scale_2*pgamma(q = exp(l_A_shape)*(((q[-index_low] - loc)/scale_2)^(shape)), shape = 1/shape, rate = 1))/(scale_1 + scale_2) 
#   }
#   
#   if(log){
#     z <- log(z)
#   }
#   return(z)
# }
# 
# qaggd <- function(p, loc, scale_1, scale_2, shape){
#   if(!is.null(dim(p))){
#     stop("q must be a vector")
#   }
#   if(any(p < 0) | any(p > 1)){
#     stop("p must be in the region [0,1]")
#   }
#   if(any(scale_1 <= 0) | any(scale_2 <= 0)){
#     stop("scale must be a positive real number")
#   }
#   
#   l_A_shape <- (shape/2)*(lgamma(3/shape) - lgamma(1/shape))
#   
#   n <- length(p)
#   z <- rep(NA, n)
#   index_low <- which(p < (scale_1)/(scale_1 + scale_2))
#   
#   if(is_empty(index_low)){
#     z <-  loc + scale_2*qgamma((p*(scale_1 + scale_2) - scale_1)/scale_2, shape = 1/shape, rate = 1)^(1/shape)
#   }
#   else{
#     z[index_low] <- loc - scale_1*(qgamma(1 - p[index_low]*(scale_1 + scale_2)/scale_1, shape = 1/shape, rate = 1)/exp(l_A_shape))^(1/shape)
#     z[-index_low] <- loc + scale_2*(qgamma((p[-index_low]*(scale_1 + scale_2) - scale_1)/scale_2, shape = 1/shape, rate = 1)/exp(l_A_shape))^(1/shape) 
#   }
#   
#   return(z)
# }
