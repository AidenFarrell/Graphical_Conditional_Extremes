################################################################################
#functions re cliques from Engleke Github

#' Get Cliques and Separators of a graph
#'
#' Finds all cliques, separators, and (recursively) separators of separators
#' in a graph.
#'
#' @param graph An \[`igraph::graph`\] object
#' @return A list of vertex sets that represent the cliques and (recursive)
#' separators of `graph`, ordered such that separators come before cliques they separate.
#'
#' @keywords internal
get_cliques_and_separators <- function(graph, sortIntoLayers = FALSE, includeSingletons = FALSE){
  # start with maximal cliques of graph
  newCliques <- lapply(igraph::max_cliques(graph), as.numeric)
  allCliques <- newCliques
  
  while(length(newCliques) > 0){
    # compute all separators
    separators <- get_separators(allCliques, includeSingletons)
    # check which are actually new cliques
    newCliques <- setdiff(separators, allCliques)
    # add to list of cliques
    allCliques <- c(newCliques, allCliques)
  }
  if(!sortIntoLayers){
    return(allCliques)
  }
  layers <- sort_cliques_and_separators(allCliques)
  return(layers)
}

#' Order Cliques
#'
#' Orders the cliques in a connected decomposable graph so that they fulfill the running intersection property.
#' @keywords internal
order_cliques <- function(cliques) {
  n <- length(cliques)
  ret <- list()
  for (i in 1:n) {
    foundNextClique <- FALSE
    for (j in seq_along(cliques)) {
      candidate <- cliques[[j]]
      rest <- Reduce(union, cliques[-j], c())
      separator <- intersect(candidate, rest)
      sepInClique <- sapply(cliques[-j], function(C) all(separator %in% C))
      if (length(sepInClique) == 0 || any(sepInClique)) {
        # add clique to return list
        # ret[[i]] <- candidate
        ret <- c(list(candidate), ret)
        # remove clique from input list
        cliques[j] <- NULL
        foundNextClique <- TRUE
        break
      }
    }
    if (!foundNextClique) {
      stop("Graph not decomposable or not connected!")
    }
  }
  return(ret)
}

# Get all separators (non-recursive) between a set of cliques
get_separators <- function(cliques, includeSingletons = FALSE){
  separators <- list()
  for(i in seq_along(cliques)){
    for(j in seq_len(i-1)){
      sep <- sort(intersect(cliques[[i]], cliques[[j]]))
      if(length(sep)>1 || (includeSingletons && length(sep) == 1)){
        separators <- c(separators, list(sep))
      }
    }
  }
  separators <- unique(separators)
  return(separators)
}

################################################################################
#Coles and Tawn (1991) semi-parametric approach for transformation to uniform margins
semiparametric_unif_transform <- function(data, u, scale, shape){
  #get length of data
  n <- length(data)
  
  #get index
  index <- which(data > u)
  
  #set up empty vector for output
  z <- rep(NA, n)
  
  #get below the threshold
  z1 <- rank(data)/(n + 1)
  z[-index] <- z1[-index]
  
  #get above the threhsold
  pu <- length(index)/n
  if(abs(shape) < 1e-6){
    z2 <- 1 - pu*exp(-exp(-(data - u)/scale))
  }else{
    z2 <- 1 - pu*pmax(0, 1 + shape*(data - u)/scale)^(-1/shape)
  }
  z[index] <- z2[index]
  
  return(z)
}

semiparametric_gpd_transform <- function(p, u_unif, gpd_par, distFun){
  
  # extract MLEs
  u <- gpd_par[1]
  scale <- gpd_par[2]
  shape <- gpd_par[3]
  
  #index
  index <- which(p <= u_unif)
  
  #get the output
  if(is_empty(index)){
    if(abs(shape) <= 1e-6){
      z <- u + scale*log((1 - u_unif)/(1 - p))
    }
    else{
      z <- u + (scale/shape)*(((1 - u_unif)/(1 - p))^(shape) - 1)
    }
  }
  else{
    #set up empty vector for the output
    n <- length(p)
    z <- rep(NA, n)
    
    #below the threshold
    z[index] <- quantile(distFun, p[index])
    
    #above the threshold
    if(abs(shape) <= 1e-6){
      z[-index] <- u_orig + scale*log((1 - u_unif)/(1 - p[-index]))
    }
    else{
      z[-index] <- u_orig + (scale/shape)*(((1 - u_unif)/(1 - p[-index]))^(shape) - 1)
    }
  }
  return(z)
}

################################################################################
## Take a some data X and transform onto Laplace margins
## Uses the semi-parametric method of Coles & Tawn (1991) to do this
## Threshold is selected using Murphy et al. (2024)

X_to_Laplace <- function(x, q = seq(0.55, 0.99, by = 0.01), k = 200, m = 500){
  ## Checks on the inputs
  if(!is.numeric(q)){stop("q must be a vector of probabilities")}
  if(min(q) <= 0.5 | max(q) > 0.999){stop("q must be in the interval (0.5, 0.999)")}
  
  ## Perform the threshold selection
  u_poss <- quantile(x, q)
  u_out <- thresh_qq_metric(data = x, thresh = u_poss, k = k, m = m)
  
  ## Get the output from the threshold selection
  u_star <- u_out$thresh
  qu_star <- q[which(u_poss == u_star)]
  scale_star <- u_out$par[1]
  shape_star <- u_out$par[2]
  
  ## Transform the data onto Laplace margins
  y <- qlaplace(semiparametric_unif_transform(data = x, u = u_star, scale = scale_star, shape = shape_star))
  
  ## Output
  out <- list(X = x,
              Y = y,
              u = u_star,
              qu = qu_star,
              scale = scale_star,
              scale = shape_star)
  return(out)
}

################################################################################
## Functions to fit the conditional extremes model using a graphical structure
Cond_extremes_graph <- function(data, cond, graph = NA, 
                                constrain = TRUE, q = c(0,1), v = 10, aLow = -1, 
                                maxit = 1e+6, start = c(0.1, 0.1), nOptim = 1){
  
  theCall <- match.call()
  # if (!inherits(x, "migpd")) 
  #   stop("you need to use an object created by migpd")
  # margins <- list(casefold(margins), p2q = switch(casefold(margins), 
  #                                                 gumbel = function(p) -log(-log(p)), laplace = function(p) ifelse(p < 
  #                                                                                                                    0.5, log(2 * p), -log(2 * (1 - p)))), q2p = switch(casefold(margins), 
  #                                                                                                                                                                       gumbel = function(q) exp(-exp(-q)), laplace = function(q) ifelse(q < 
  #                                                                                                                                                                                                                                          0, exp(q)/2, 1 - 0.5 * exp(-q))))
  # x <- mexTransform(x, margins = margins, method = marTransform, 
  #                   r = referenceMargin)
  # x$referenceMargin <- referenceMargin
  # if (margins[[1]] == "gumbel" & constrain) {
  #   warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
  #   constrain <- FALSE
  # }
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(is.character(cond)){ 
    cond <- match(cond, dimnames(x$transformed)[[2]])
  }
  # if(missing(dqu)){
  #   message("Assuming same quantile for dependence thesholding as was used \n     to fit corresponding marginal model...\n")
  #   dqu <- x$mqu[which]
  # }
  # dth <- quantile(x$transformed[, which], dqu)
  dependent <- (1:(dim(data)[[2]]))[-cond]
  # if(length(dqu) < length(dependent)){
  #   dqu <- rep(dqu, length = length(dependent))
  #   aLow <- if else(margins[[1]] == "gumbel", 10^(-10), -1 + 10^(-10))
  # }
  if(missing(start)){
    start <- matrix(rep(c(0.1, 0.1), length(dependent)), nrow = length(dependent))
  }
  if(length(start) != 2 * length(dependent)){
    stop("start should be of type 'mex' or be a vector of length 2, or be a matrix with 2 rows and ncol equal to the number of dependence models to be estimated")
  }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    comps_final <- sapply(1:(d-1), list)
  }
  else{
    #Remove the cond node from the graph
    adj_matrix <- as.matrix(as_adj(graph))
    adj_matrix <- adj_matrix[-cond, -cond]
    graph_z <- graph_from_adjacency_matrix(adj_matrix)
    V(graph_z)$id <- 1:(d-1)
    
    #get the components of the new graph
    comps_z <- components(graph_z)
    #separate the nodes into their components
    comps_memb <- unique(comps_z$membership)
    comps_final <- lapply(comps_memb, function(x){V(graph_z)$id[which(comps_z$membership == x)]})
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- data[,cond]
  ydep <- data[,-cond]
  
  #fit the profile log-likelihood
  res <- lapply(comps_final, function(x){
    qfun(yex = yex, ydep = as.matrix(ydep[,x]), comps = x, 
         constrain = constrain, aLow = aLow, q = q, v = v,
         maxit = maxit, start = start[x,], nOptim = nOptim)})
  out <- list()
  out$par$a <- matrix(sapply(res, function(x){x$par$a}), nrow = 1)
  out$par$b <- matrix(sapply(res, function(x){x$par$b}), nrow = 1)
  out$par$mu <- matrix(sapply(res, function(x){x$par$mu}), nrow = 1)
  out$par$Sigma <- lapply(res, function(x){matrix(x$par$Sigma, ncol(x$par$Sigma))})
  
  colnames(out$par$a) <- colnames(out$par$b) <- colnames(out$par$mu) <-
    sapply(dependent, function(x){paste0("Column", x)})
  
  out$comps <- lapply(comps_final, function(x){dependent[x]})
  
  out$loglike <- sapply(res, function(x){-x$value})
  out$convergence <- sapply(res, function(x){-x$convergence})
  
  out$Z <- do.call(cbind, lapply(res, function(x){x$Z}))
  
  return(out)
  #   #get alpha and beta parameters and the residuals
  #   pl_par <- do.call(cbind, lapply(pl_fit, function(x){t(matrix(x$par, ncol = 2))}))
  #   a <- lapply(comps_final, function(x){pl_par[1,x]})
  #   b <- lapply(comps_final, function(x){pl_par[2,x]})
  #   alpha_yi <- lapply(a, function(x){sapply(x, function(z){yex*z})})
  #   beta_yi <- lapply(b, function(x){sapply(x, function(z){yex^z})})
  #   Z <- pmap(.l = list(a = alpha_yi, b = beta_yi,
  #                       y = lapply(comps_final, function(x){ydep[,x]})),   
  #             .f = function(a, b, y){(y - a)/b})
  #   
  #   #get mu and sigma parameters
  #   mu <- lapply(Z, function(x){apply(x, 2, mean)})
  #   sigma <- lapply(Z, function(x){cov(x)})
  #   sigma_hat <- diag(0, length(do.call(c, comps_final)))
  #   sigma_locations <- lapply(comps_final, function(x){permutations(n = length(x), r = 2, v = x, repeats.allowed = TRUE)})
  #   for(i in 1:length(comps_final)){
  #     sigma_hat[sigma_locations[[i]]] <- sigma[[i]]
  #   }
  #   
  #   #organise the output
  #   
  #   #alpha, beta and mu parameters
  #   pl_par <- as.data.frame(rbind(pl_par, do.call(c, mu)))
  #   rownames(pl_par) <- c("a", "b", "m")
  #   colnames(pl_par) <- sapply((1:d)[-which], function(x){paste0("Column", x)})
  #   pl_par_by_comp <- lapply(comps_final, function(x){pl_par[,x]})
  #   
  #   #Sigma parameter
  #   if(is_igraph(graph)){
  #     #Label using node id rather than the new id above
  #     V(graph_z)$id <- (1:d)[-which]
  #     #get the components of this new graph
  #     comps_z <- components(graph_z)
  #     #separate the members
  #     comps_memb <- unique(comps_z$membership)
  #     if(length(comps_memb) != 1){
  #       comps_final <- lapply(comps_memb, function(x){V(graph_z)$id[which(comps_z$membership == x)]})
  #     }
  #   }
  #   
  #   #Get the residuals together
  #   Z <- data.frame(do.call(cbind, Z))
  #   colnames(Z) <- sapply((1:d)[-which], function(x){paste0("Column", x)})
  #   
  #   #Get the output
  #   out <- list(par = list(main_by_comp = pl_par_by_comp,
  #                          Sigma_by_comp = lapply(sigma, function(x){matrix(x, ncol(x))}),
  #                          Sigma_full = sigma_hat),
  #               convergence = max(sapply(pl_fit, function(x){x$convergence})),
  #               components = comps_final,
  #               llh = sapply(pl_fit, function(x){-x$value}),
  #               Z = Z)
  #   
  #   return(out)
  # }
}

################################################################################
## Above but with the profile log-likelihood using cliques and separators
Cond_extremes_graph <- function(data, cond = 1, graph = NA, 
                                constrain = TRUE, q = c(0,1), v = 10, aLow = -1, 
                                maxit = 1e+6, start = c(0.1, 0.1), nOptim = 1){
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  dependent <- (1:d)[-cond]
  
  ## check the conditioning random variable is valid
  if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(cond%%1 != 0 | cond <= 0 | cond > d){
    stop("cond must be a single positie integer")
  }
  
  ## get the starting parameters
  if(missing(start)){
    start <- c(0.1, 0.1)
  }
  else if(!is.numeric(start)){
    stop("start must be a vector")
  }
  else if(length(start) > 2){
    stop("start must be a vector")
  }
  else if(any(abs(start)) > 1){
    stop("Initial starting values are outside the parameter sapce")
  }
  if(length(start) == 2){
    start <- matrix(rep(start, d-1), ncol = 2)
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- data[,cond]
  ydep <- data[,-cond]
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## Fit the model
    res <- lapply(1:(d-1), function(i){
      qfun_indep(yex = yex, ydep = as.matrix(ydep[,i]), 
                 constrain = constrain, aLow = aLow, q = q, v = v,
                 maxit = maxit, start = start[i,], nOptim = nOptim)})
    
    #get the output
    out <- list()
    out$par$a <- matrix(do.call(c, lapply(res, function(x){x$par$a})), nrow = 1)
    out$par$b <- matrix(do.call(c, lapply(res, function(x){x$par$b})), nrow = 1)
    out$par$mu <- matrix(do.call(c, lapply(res, function(x){x$par$mu})), nrow = 1)
    
    Sigma_out <- matrix(0, ncol = d-1, nrow = d-1)
    for(i in 1:length(cliques)){
      Sigma_elements <- permutations(n = length(cliques[[i]]), r = 2, v = cliques[[i]], repeats.allowed = TRUE)
      Sigma_out[Sigma_elements] <-  res[[i]]$par$Sigma
    }
    out$par$Sigma <- Sigma_out
    
    colnames(out$par$a) <- colnames(out$par$b) <- colnames(out$par$mu) <-
      sapply(dependent, function(x){paste0("Column", x)})
    
    out$loglike <- sapply(res, function(x){-x$value})
    out$convergence <- sapply(res, function(x){x$convergence})
    out$Z <- do.call(cbind, lapply(res, function(x){x$Z}))
  }
  else{
    #Graph is provided so we need to figure out which edges are present
    graph_cond <- delete.vertices(graph, cond)
    graphs_ind <- lapply(1:d, delete.vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    #get the cliques and the separators in the graph
    cliques <- order_cliques(get_cliques_and_separators(graph_cond))
    #now get the separators but we only include singletons a vertex only appears once in a clique of length two
    seps <- get_separators(cliques, includeSingletons = inc_sings)
    #check there is no cross over. If there is remove from the clique set
    index <- which(cliques %in% intersect(cliques, seps))
    if(length(index) != 0){
      cliques <- cliques[-index]
    }
    
    #sort the cliques and separators into ascending order
    cliques <- order_cliques(cliques)
    
    cliques <- lapply(cliques, function(x){sort(x)})
    seps <- lapply(seps, function(x){sort(x)})
  }
  
  #fit the profile log-likelihood
  if(is_empty(seps)){
    res <- lapply(cliques, function(x){
      qfun(yex = yex, ydep = as.matrix(ydep[,x]), 
           constrain = constrain, aLow = aLow, q = q, v = v,
           maxit = maxit, start = start[x,], nOptim = nOptim)}) 
    
    return(out)
  }
  else{
    Sigma_zero <- matrix(0, d-1, d-1)
    Sig_non_zero <- upper.tri(Sigma_zero, diag = TRUE)
    # Sig_non_zero <- as_edgelist(graph_cond)
    # Sig_non_zero <- rbind(Sig_non_zero, 
    #                       matrix(rep(1:(d-1), 2), ncol = 2))
    
    res <- qfun_new(yex = yex, ydep = as.matrix(ydep), comps = 1:(d-1), 
                    constrain = constrain, aLow = aLow, q = q, v = v,
                    cliques = cliques, seps = seps, Sig_non_zero = Sig_non_zero,
                    maxit = maxit, start = start, nOptim = nOptim) 
    
    #get the output
    out <- list()
    out$par$a <- matrix(res$par$a, nrow = 1)
    out$par$b <- matrix(res$par$b, nrow = 1)
    out$par$mu <- matrix(res$par$mu, nrow = 1)
    out$par$Sigma <- matrix(res$par$Sigma, nrow = d - length(cond), ncol = d-length(cond))
    
    colnames(out$par$a) <- colnames(out$par$b) <- colnames(out$par$mu) <-
      sapply(dependent, function(x){paste0("Column", x)})
    
    out$loglike <- -res$value
    out$convergence <- res$convergence
    out$Z <- res$Z
    
    return(out)
  }
}

qfun_indep <- function(yex, ydep, constrain, q, v, aLow, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, constrain, q, v, aLow){
    #get the starting parameters
    d <- ncol(ydep)
    a <- param[1:d]
    b <- param[(d + 1):length(param)]
    
    #get the value of the profile log_likelihood for fixed a and b
    res <- ProfileLogLik_a_b_indep(yex, ydep, a, b, constrain, q, v, aLow)
    res$profLik
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   yex = yex, ydep = ydep,
                   constrain = constrain, aLow = aLow, q = q, v = v, method = "BFGS"),
             silent = TRUE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       yex = yex, ydep = ydep,
                       constrain = constrain, aLow = aLow, v = v, q = q, method = "BFGS"),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d <- ncol(ydep)
    alpha_hat <- fit$par[1:d]
    beta_hat <- fit$par[(d + 1):length(fit$par)]
    
    #obtain the residuals
    a_yi_hat <- sapply(alpha_hat, function(x){yex*x})
    b_yi_hat <- sapply(beta_hat, function(x){yex^x})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    out <- list()
    out$par <- list(a = alpha_hat, b = beta_hat, mu = apply(Z, 2, mean), Sigma = cov(Z))
    out$value <- fit$value
    out$convergence <- fit$convergence
    out$Z <- Z
  }
  return(out)
}

#function for the profile log-likelihood
ProfileLogLik_a_b_indep <- function(yex, ydep, a, b, constrain = TRUE, q, v = 10, aLow = -1){
  
  #get the residuals
  a_yi <- sapply(a, function(x){yex*x})
  b_yi <- sapply(b, function(x){yex^x})
  Z <- (ydep - a_yi)/b_yi
  
  #get the mean and covariance matrix of the residuals
  mu <- apply(Z, 2, mean)
  sigma <- cov(Z)
  
  #get the negative log-likelihood for the given parameters
  nllh <- negloglike_indep(yex, ydep, a_yi, b_yi, a, b, mu, sigma, constrain, q, v, aLow)
  res <- list(profLik = nllh, mu = mu, sigma = sigma)
  return(res)
}

negloglike_indep <- function(yex, ydep, a_yi, b_yi, a, b, m, sig, constrain, q, v, aLow){
  
  #conditions on the parameters
  Delta <- 10^40
  delta <- 10^(-10)
  
  if(any(is.na(m)) | any(is.nan(m))){
    res <- Delta
  }
  else if(any(a < aLow[1]) | any(a > 1 - delta) | any(b > 1 - delta)){
    res <- Delta
  }
  else if(any(is.nan(sig)) | any(is.infinite(sig))){
    res <- Delta
  }
  else if(is.complex(sig) | any(diag(sig) <= 0) | any(eigen(sig)$values <= 0) | !isTRUE(all.equal(sig, t(sig)))){
    res <- Delta
  }
  else{
    if(length(a) == 1){
      #get the mean and sd matrix of the 1-dimensional Gaussian
      mu <- lapply(apply(a_yi + b_yi*m, 1, list), function(y){y[[1]]})
      sigma <- apply(b_yi, 1, function(x){x*sig*x}, simplify = FALSE)
      ydep_list <- sapply(ydep, list)
    }
    else{
      #get the mean and covariance matrix of the d-dimensional Gaussian
      mu <- lapply(apply(a_yi + t(apply(b_yi, 1, function(x){x*m})), 1, list), function(y){y[[1]]})
      sigma <- apply(b_yi, 1, function(x){diag(x)%*%sig%*%diag(x)}, simplify = FALSE) 
      ydep_list <- lapply(apply(ydep, 1, list), function(y){y[[1]]})
    }
    
    #Calculate the log-likelihood
    res <- pmap(.l = list(x = ydep_list, mean = mu, sigma = sigma),
                .f = dmvnorm, log = TRUE) 
    res <- -Reduce("+", res)
    
    #check the value is valid
    if(is.infinite(res)){
      if(res < 0) {
        res <- -Delta
      }
      else{ 
        res <- Delta
      }
      warning("Infinite value of Q in mexDependence")
    }
    #Checks if the constraints are satisfied to reduce the parameter space
    else if(constrain){ 
      zpos <- apply(ydep - yex, 2, quantile, probs = q, simplify = FALSE)
      z <- apply((ydep - a_yi)/b_yi, 2, quantile, probs = q, simplify = FALSE)
      zneg <- apply(ydep + yex, 2, quantile, probs = q, simplify = FALSE)
      
      constraints_sat <- sapply(pmap(.l = list(a = a, b = b, z1 = z, zpos1 = zpos, zneg1 = zneg),
                                     .f = function(a, b, z1, zpos1, zneg1){
                                       pmap(.l = list(a = a, b = b, z = z1, zpos = zpos1, zneg = zneg1), 
                                            .f = texmex:::ConstraintsAreSatisfied, v)}), function(x){all(x == TRUE)})
      if(!all(constraints_sat)){
        res <- Delta
      }
    }
  }
  res
}

#Might be able to get rid of the comps argument
qfun_new <- function(yex, ydep, comps, constrain, q, v, aLow, cliques, seps, Sig_non_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, comps, constrain, q, v, aLow, cliques, seps, Sig_non_zero){
    #get the starting parameters
    n_comps <- length(comps)
    a <- param[1:n_comps]
    b <- param[(n_comps + 1):length(param)]
    
    #get the value of the profile log_likelihood for fixed a and b
    res <- ProfileLogLik_a_b_new(yex, ydep, a, b, constrain, q, v, aLow, cliques, seps, Sig_non_zero)
    res$profLik
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   yex = yex, ydep = ydep, comps = comps, cliques = cliques, seps = seps, Sig_non_zero = Sig_non_zero,
                   constrain = constrain, aLow = aLow, q = q, v = v, method = "BFGS"),
             silent = TRUE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       yex = yex, ydep = ydep, comps = comps, cliques = cliques, seps = seps,
                       constrain = constrain, aLow = aLow, v = v, q = q, method = "BFGS"),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    alpha_hat <- fit$par[1:length(comps)]
    beta_hat <- fit$par[-(1:length(alpha_hat))][1:length(comps)]
    
    #obtain the residuals
    a_yi_hat <- sapply(alpha_hat, function(x){yex*x})
    b_yi_hat <- sapply(beta_hat, function(x){yex^x})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    sigma_Z <- cov(Z)
    sigma <- matrix(0, nrow = ncol(Z), ncol = ncol(Z))
    sigma[Sig_non_zero] <- sigma_Z[Sig_non_zero]
    sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
    
    out <- list()
    out$par <- list(a = alpha_hat, b = beta_hat, mu = apply(Z, 2, mean), Sigma = sigma)
    out$value <- fit$value
    out$convergence <- fit$convergence
    out$Z <- Z
  }
  return(out)
}

#function for the profile log-likelihood
ProfileLogLik_a_b_new <- function(yex, ydep, a, b, constrain = TRUE, q, v = 10, aLow = -1, cliques, seps, Sig_non_zero){
  
  #get the residuals
  a_yi <- sapply(a, function(x){yex*x})
  b_yi <- sapply(b, function(x){yex^x})
  Z <- (ydep - a_yi)/b_yi
  
  #get the mean and covariance matrix of the residuals
  mu <- apply(Z, 2, mean)
  sigma_Z <- cov(Z)
  sigma <- matrix(0, nrow = ncol(Z), ncol = ncol(Z))
  sigma[Sig_non_zero] <- sigma_Z[Sig_non_zero]
  sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
  
  #get the negative log-likelihood for the given parameters
  nllh <- negloglike_new(yex, ydep, a_yi, b_yi, a, b, mu, sigma, constrain, q, v, aLow, cliques, seps)
  res <- list(profLik = nllh, mu = mu, sigma = sigma)
  return(res)
}

negloglike_new <- function(yex, ydep, a_yi, b_yi, a, b, m, sig, constrain, q, v, aLow, cliques, seps){
  
  #conditions on the parameters
  Delta <- 10^40
  delta <- 10^(-10)
  
  if(any(is.na(m)) | any(is.nan(m))){
    res <- Delta
  }
  else if(any(a < aLow[1]) | any(a > 1 - delta) | any(b > 1 - delta)){
    res <- Delta
  }
  else if(any(is.nan(sig)) | any(is.infinite(sig))){
    res <- Delta
  }
  else if(is.complex(sig) | any(diag(sig) <= 0) | any(eigen(sig)$values <= 0) | !isTRUE(all.equal(sig, t(sig)))){
    res <- Delta
  }
  else{
    #get the mean and covaraince matrix of the d-dimensional Guassian
    mu <- lapply(apply(a_yi + t(apply(b_yi, 1, function(x){x*m})), 1, list), function(y){y[[1]]})
    sigma <- apply(b_yi, 1, function(x){diag(x)%*%sig%*%diag(x)}, simplify = FALSE) 
    ydep_list <- lapply(apply(ydep, 1, list), function(y){y[[1]]})
    
    #Separate out into the cliques and the separator(s)
    mu_cliques <- lapply(cliques, function(x){lapply(mu, function(y){y[x]})})
    sigma_cliques <- lapply(cliques, function(x){lapply(sigma, function(y){as.matrix(y[x,x])})})
    ydep_cliques <- lapply(cliques, function(x){lapply(ydep_list, function(y){y[x]})})
    
    mu_seps <- lapply(seps, function(x){lapply(mu, function(y){y[x]})})
    sigma_seps <- lapply(seps, function(x){lapply(sigma, function(y){as.matrix(y[x,x])})})
    ydep_seps <- lapply(seps, function(x){lapply(ydep_list, function(y){y[x]})})
    
    #Calculate the log-likelihood for the cliques and separator(s)
    res_cliques <- lapply(1:length(cliques), function(x){
      pmap(.l = list(x = ydep_cliques[[x]], mean = mu_cliques[[x]], sigma = sigma_cliques[[x]]), 
           .f = dmvnorm, log = TRUE)})
    res_cliques <- sum(sapply(1:length(cliques), function(x){Reduce("+", res_cliques[[x]])}))
    
    res_seps <- lapply(1:length(seps), function(x){
      pmap(.l = list(x = ydep_seps[[x]], mean = mu_seps[[x]], sigma = sigma_seps[[x]]), 
           .f = dmvnorm, log = TRUE)})
    res_seps <- sum(sapply(1:length(seps), function(x){Reduce("+", res_seps[[x]])}))
    
    #Calculate the total negative log-likelihood 
    res <- -(res_cliques - res_seps)
    
    #check the value is valid
    if(is.infinite(res)){
      if(res < 0) {
        res <- -Delta
      }
      else{ 
        res <- Delta
      }
      warning("Infinite value of Q in mexDependence")
    }
    #Checks if the constraints are satisfied to reduce the parameter space
    else if(constrain){ 
      zpos <- apply(ydep - yex, 2, quantile, probs = q, simplify = FALSE)
      z <- apply((ydep - a_yi)/b_yi, 2, quantile, probs = q, simplify = FALSE)
      zneg <- apply(ydep + yex, 2, quantile, probs = q, simplify = FALSE)
      
      constraints_sat <- do.call(rbind, pmap(.l = list(a = a, b = b, z = z, zpos1 = zpos, zneg1 = zneg),
                                             .f = function(a, b, z1, zpos1, zneg1){sapply(pmap(.l = list(a = a, b = b, z = z1, zpos = zpos1, zneg = zneg1), 
                                                                                               .f = texmex:::ConstraintsAreSatisfied, v), function(x){x})}))
      if(!all(constraints_sat)){
        res <- Delta
      }
    }
  }
  res
}

qfun_HT <- function(yex, ydep, constrain, q, v, aLow, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, constrain, q, v, aLow){
    ## Extract the starting parameters
    a <- param[1]
    b <- param[2]
    
    ## Get the residuals
    a_yi <- yex*a
    b_yi <- yex^b
    Z <- (ydep - a_yi)/b_yi
    
    ## Obtain the mean and covariance matrix of the residuals
    mu <- mean(Z)
    sigma <- as.numeric(var(Z))
    
    ## Conditions on the parameters
    Delta <- 10^40
    delta <- 10^(-10)
    if(any(a < aLow[1]) | any(a > 1 - delta) | any(b > 1 - delta)){
      res <- Delta
    }
    else{
      mu <- a_yi + mu*b_yi
      sigma <- sigma*(b_yi^2)
      res <- -sum(dnorm(ydep, mean = mu, sd = sqrt(sigma), log = TRUE))
    }
    
    if(is.infinite(res)){
      if(res < 0) {
        res <- -Delta
      }
      else{ 
        res <- Delta
      }
      warning("Infinite value of Q in mexDependence")
    }
    
    #Checks if the constraints are satisfied to reduce the parameter space
    else if(constrain){ 
      zpos <- quantile(ydep - yex, probs = q)
      z <- quantile(Z, probs = q)
      zneg <- quantile(ydep + yex, probs = q)
      
      constraints_sat <- texmex:::ConstraintsAreSatisfied(a, b, z, zpos, zneg, v)
      
      # constraints_sat <- sapply(pmap(.l = list(a = a, b = b, z1 = z, zpos1 = zpos, zneg1 = zneg),
      #                                .f = function(a, b, z1, zpos1, zneg1){
      #                                  pmap(.l = list(a = a, b = b, z = z1, zpos = zpos1, zneg = zneg1), 
      #                                       .f = texmex:::ConstraintsAreSatisfied, v)}), function(x){all(x == TRUE)})
      if(!all(constraints_sat)){
        res <- Delta
      }
    }
    res
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   yex = yex, ydep = ydep,
                   constrain = constrain, aLow = aLow, q = q, v = v, method = "BFGS"),
             silent = TRUE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       yex = yex, ydep = ydep,
                       constrain = constrain, aLow = aLow, v = v, q = q, method = "BFGS"),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    ## Extract MLEs of alpha and beta
    a_hat <- fit$par[1]
    b_hat <- fit$par[2]
    
    ## Obtain the residuals
    Z <- (ydep - yex*a_hat)/(yex^b_hat)
    
    ## Organise the output
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = mean(Z), Sigma = var(Z))
    out$value <- fit$value
    out$convergence <- fit$convergence
    out$Z <- Z
  }
  return(out)
}

################################################################################
## Conditional extremes model assuming the residuals follow a (multivariate) asymmetric generalised Gaussian distribution
Cond_extremes_MVAGGD <- function(data, cond, graph = NA, start,
                                 constrain = TRUE, q = c(0,1), v = 10, aLow = -1,
                                 maxit = 1e+6, nOptim = 1){
  
  theCall <- match.call()
  # if (!inherits(x, "migpd")) 
  #   stop("you need to use an object created by migpd")
  # margins <- list(casefold(margins), p2q = switch(casefold(margins), 
  #                                                 gumbel = function(p) -log(-log(p)), laplace = function(p) ifelse(p < 
  #                                                                                                                    0.5, log(2 * p), -log(2 * (1 - p)))), q2p = switch(casefold(margins), 
  #                                                                                                                                                                       gumbel = function(q) exp(-exp(-q)), laplace = function(q) ifelse(q < 
  #                                                                                                                                                                                                                                          0, exp(q)/2, 1 - 0.5 * exp(-q))))
  # x <- mexTransform(x, margins = margins, method = marTransform, 
  #                   r = referenceMargin)
  # x$referenceMargin <- referenceMargin
  # if (margins[[1]] == "gumbel" & constrain) {
  #   warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
  #   constrain <- FALSE
  # }
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(is.character(cond)){ 
    cond <- match(cond, dimnames(x$transformed)[[2]])
  }
  # if(missing(dqu)){
  #   message("Assuming same quantile for dependence thesholding as was used \n     to fit corresponding marginal model...\n")
  #   dqu <- x$mqu[cond]
  # }
  # dth <- quantile(x$transformed[, cond], dqu)
  dependent <- (1:(dim(data)[[2]]))[-cond]
  # if(length(dqu) < length(dependent)){
  #   dqu <- rep(dqu, length = length(dependent))
  #   aLow <- if else(margins[[1]] == "gumbel", 10^(-10), -1 + 10^(-10))
  # }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graph_cond <- delete_vertices(graph, cond)
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    #get the cliques and the separators in the graph
    cliques <- order_cliques(get_cliques_and_separators(graph_cond))
    #now get the separators but we only include singletons a vertex only appears once in a clique of length two
    seps <- get_separators(cliques, includeSingletons = inc_sings)
    #check there is no cross over. If there is remove from the clique set
    index <- which(cliques %in% intersect(cliques, seps))
    if(length(index) != 0){
      cliques <- cliques[-index]
    }
    
    #sort the cliques and separators into ascending order
    cliques <- order_cliques(cliques)
    
    cliques <- lapply(cliques, function(x){sort(x)})
    seps <- lapply(seps, function(x){sort(x)})
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  n_comps <- length(cliques)
  if(is_empty(seps)){
    
    #fit the model
    res <- lapply(1:n_comps, function(i){qfun_MVAGGD(yex = yex, ydep = as.matrix(ydep[,cliques[[i]]]), comps = cliques[[i]], 
                                                     constrain = constrain, aLow = aLow, q = q, v = v,
                                                     maxit = maxit, start = start[cliques[[i]],], nOptim = nOptim)
    })
    
    ## Get the output from the model
    out <- list()
    
    par_out <- rbind(do.call(c, lapply(res, function(x){x$par$a})), 
                     do.call(c, lapply(res, function(x){x$par$b})), 
                     do.call(c, lapply(res, function(x){x$par$mu})),
                     do.call(c, lapply(res, function(x){x$par$sigma_1})),
                     do.call(c, lapply(res, function(x){x$par$sigma_2})),
                     do.call(c, lapply(res, function(x){x$par$shape})))
    par_out <- as.data.frame(par_out)
    rownames(par_out) <- c("alpha", "beta", "mu", "sigma_1", "sigma_2", "shape")
    colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
    
    Sigma_out <- diag(1, d-1)
    for(i in 1:n_comps){
      if(length(cliques[[i]]) > 1){
        Sig_entries <- permutations(n = length(cliques[[i]]), r = 2, v = min(cliques[[i]]):max(cliques[[i]]), repeats.allowed = TRUE)
        Sigma_out[Sig_entries] <- res[[i]]$par$Sigma
      }
    }
    Sigma_out[lower.tri(Sigma_out)] <- t(Sigma_out)[lower.tri(Sigma_out)]
    out$par <- list(main = par_out, Sigma = Sigma_out)
    
    if(n_comps == 1){
      out_llh <- data.frame(llh = -res[[1]]$value)
      colnames(out_llh) <- c("Total")
    }
    else{
      out_llh <- c(sapply(res, function(x){-x$value}), sum(sapply(res, function(x){-x$value})))
      out_llh <- t(as.data.frame(out_llh))
      colnames(out_llh) <- c(sapply(1:n_comps, function(i){paste0("Comp", i)}), "Total") 
    }
    rownames(out_llh) <- c("llh")
    out$loglike <- out_llh
    
    out$convergence <- sapply(res, function(x){x$convergence})
    
    out_Z <- do.call(cbind, lapply(res, function(x){x$Z}))
    out_Z <- as.data.frame(out_Z)
    colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
    out$Z <- out_Z
  }
  else{
    res <- qfun_MVAGGD_graph(yex = yex, ydep = as.matrix(ydep), 
                             comps = 1:(d-1), cliques = cliques, seps = seps,
                             constrain = constrain, aLow = aLow, q = q, v = v,
                             maxit = maxit, start = start, nOptim = nOptim)
    
    ## Get the output from the model
    out <- list()
    
    par_out <- rbind(res$par$a, res$par$b, res$par$mu, res$par$sigma_1, res$par$sigma_2, res$par$shape)
    par_out <- as.data.frame(par_out)
    rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
    colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
    
    Sigma_out <- res$par$Sigma
    
    out$par <- list(main = par_out, Sigma = Sigma_out)
    
    out_llh <- data.frame(llh = -res$value)
    colnames(out_llh) <- c("Total")
    rownames(out_llh) <- c("llh")
    out$loglike <- out_llh
    
    out$convergence <- res$convergence
    
    out_Z <- res$Z
    out_Z <- as.data.frame(out_Z)
    colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
    out$Z <- out_Z
  }
  return(out)
}

qfun_MVAGGD <- function(yex, ydep, comps, constrain, q, v, aLow, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, comps, constrain, q, v, aLow, negative = FALSE){
    #get the starting parameters
    d_comps <- length(comps)
    n_data <- length(yex)
    
    a <- param[1:d_comps]
    b <- param[(d_comps + 1):(2*d_comps)]
    mu <- param[(2*d_comps + 1):(3*d_comps)]
    sigma_1 <- param[(3*d_comps + 1):(4*d_comps)]
    sigma_2 <- param[(4*d_comps + 1):(5*d_comps)]
    shape <- param[(5*d_comps + 1):(6*d_comps)]
    
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else if(any(abs(a) > 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      res <- sum(dmvagg(data = z, loc = mu, scale_1 = sigma_1, scale_2 = sigma_2,
                        shape = shape, log = TRUE)) +
        sum(log(1/b_yi))
      
      #check the value is valid
      if(is.infinite(res)){
        if(res < 0) {
          res <- return((-10^10)*(-1)^negative)
        }
        else{ 
          res <- return((-10^10)*(-1)^negative)
        }
        warning("Infinite value of Q in mexDependence")
      }
      #Checks if the constraints are satisfied to reduce the parameter space
      else if(constrain){ 
        zpos <- apply(apply(ydep, 2, function(y){y - yex}), 2, quantile, probs = q, simplify = FALSE)
        z <- apply((ydep - a_yi)/b_yi, 2, quantile, probs = q, simplify = FALSE)
        zneg <- apply(apply(ydep, 2, function(y){y + yex}), 2, quantile, probs = q, simplify = FALSE)
        
        constraints_sat <- sapply(pmap(.l = list(a = a, b = b, z1 = z, zpos1 = zpos, zneg1 = zneg),
                                       .f = function(a, b, z1, zpos1, zneg1){
                                         pmap(.l = list(a = a, b = b, z = z1, zpos = zpos1, zneg = zneg1), 
                                              .f = texmex:::ConstraintsAreSatisfied, v)}), function(x){all(x == TRUE)})
        if(!all(constraints_sat)){
          res <- return((-10^10)*(-1)^negative)
        }
      }
    }
    
    return(res*(-1)^negative)
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   yex = yex, ydep = ydep, comps = comps,
                   constrain = constrain, aLow = aLow, q = q, v = v, 
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       yex = yex, ydep = ydep, comps = comps,
                       constrain = constrain, aLow = aLow, v = v, q = q,
                       negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta and AGG
    d_comps <- length(comps)
    a_hat <- fit$par[1:d_comps]
    b_hat <- fit$par[(d_comps + 1):(2*d_comps)]
    mu_hat <- fit$par[(2*d_comps + 1):(3*d_comps)]
    sigma_1_hat <- fit$par[(3*d_comps + 1):(4*d_comps)]
    sigma_2_hat <- fit$par[(4*d_comps + 1):(5*d_comps)]
    shape_hat <- fit$par[(5*d_comps + 1):(6*d_comps)]
    
    #obtain the residuals
    a_yi_hat <- sapply(a_hat, function(a){yex*a})
    b_yi_hat <- sapply(b_hat, function(b){yex^b})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## get the correlation matrix
    Q_F_Z <- sapply(1:d_comps, function(i){
      qnorm(paggd(q = Z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))
    })
    Sigma_hat <- cor(Q_F_Z)
    
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Sigma = Sigma_hat)
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

qfun_MVAGGD_graph <- function(yex, ydep, comps, cliques, seps, constrain, q, v, aLow, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, comps, cliques, seps, constrain, q, v, aLow, negative = FALSE){
    #get the starting parameters
    d_comps <- length(comps)
    n_data <- length(yex)
    
    a <- param[1:d_comps]
    b <- param[(d_comps + 1):(2*d_comps)]
    mu <- param[(2*d_comps + 1):(3*d_comps)]
    sigma_1 <- param[(3*d_comps + 1):(4*d_comps)]
    sigma_2 <- param[(4*d_comps + 1):(5*d_comps)]
    shape <- param[(5*d_comps + 1):(6*d_comps)]
    
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else if(any(abs(a) > 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      res <- llh_dmvagg_graph(data = z, loc = mu, scale_1 = sigma_1, scale_2 = sigma_2,
                              shape = shape, cliques = cliques, seps = seps) +
        sum(log(1/b_yi))
      
      #check the value is valid
      if(is.infinite(res)){
        if(res < 0) {
          res <- return((-10^10)*(-1)^negative)
        }
        else{ 
          res <- return((-10^10)*(-1)^negative)
        }
        warning("Infinite value of Q in mexDependence")
      }
      #Checks if the constraints are satisfied to reduce the parameter space
      else if(constrain){ 
        zpos <- apply(apply(ydep, 2, function(y){y - yex}), 2, quantile, probs = q, simplify = FALSE)
        z <- apply((ydep - a_yi)/b_yi, 2, quantile, probs = q, simplify = FALSE)
        zneg <- apply(apply(ydep, 2, function(y){y + yex}), 2, quantile, probs = q, simplify = FALSE)
        
        constraints_sat <- sapply(pmap(.l = list(a = a, b = b, z1 = z, zpos1 = zpos, zneg1 = zneg),
                                       .f = function(a, b, z1, zpos1, zneg1){
                                         pmap(.l = list(a = a, b = b, z = z1, zpos = zpos1, zneg = zneg1), 
                                              .f = texmex:::ConstraintsAreSatisfied, v)}), function(x){all(x == TRUE)})
        if(!all(constraints_sat)){
          res <- return((-10^10)*(-1)^negative)
        }
      }
    }
    return(res*(-1)^negative)
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   yex = yex, ydep = ydep, comps = comps,
                   constrain = constrain, aLow = aLow, q = q, v = v,
                   cliques = cliques, seps = seps,
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       yex = yex, ydep = ydep, comps = comps,
                       constrain = constrain, aLow = aLow, v = v, q = q,
                       negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d_comps <- length(comps)
    a_hat <- fit$par[1:d_comps]
    b_hat <- fit$par[(d_comps + 1):(2*d_comps)]
    mu_hat <- fit$par[(2*d_comps + 1):(3*d_comps)]
    sigma_1_hat <- fit$par[(3*d_comps + 1):(4*d_comps)]
    sigma_2_hat <- fit$par[(4*d_comps + 1):(5*d_comps)]
    shape_hat <- fit$par[(5*d_comps + 1):(6*d_comps)]
    
    #obtain the residuals
    a_yi_hat <- sapply(a_hat, function(a){yex*a})
    b_yi_hat <- sapply(b_hat, function(b){yex^b})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## get the correlation matrix
    Q_F_Z <- sapply(1:d_comps, function(i){
      qnorm(paggd(q = Z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))
    })
    Sigma_hat <- cor(Q_F_Z)
    
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Sigma = Sigma_hat)
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

llh_dmvagg_graph <- function(data, loc, scale_1, scale_2, shape, cliques, seps, negative = FALSE){
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  ## Transform data onto joint standard Gaussian margins
  Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i]))}))
  
  ## correlation matrix to be used
  Sigma <- as.matrix(cor(Q_F_z))
  
  l_f_z <- sapply(1:d, function(i){daggd(x = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i], log = TRUE)})
  l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
  
  l_mvnorm_cliques <- sapply(cliques, function(i){sum(mvtnorm::dmvnorm(x = as.matrix(Q_F_z[,i]), mean = rep(0, length(i)), sigma = Sigma[i,i], log = TRUE))})
  l_mvnorm_seps  <- sapply(seps, function(i){sum(mvtnorm::dmvnorm(x = as.matrix(Q_F_z[,i]), mean = rep(0, length(i)), sigma = as.matrix(Sigma[i,i]), log = TRUE))})
  
  llh_cliques <- sum(l_mvnorm_cliques + sapply(cliques, function(i){sum((l_f_z - l_dnorm)[,i])}))
  llh_seps <- sum(l_mvnorm_seps + sapply(seps, function(i){sum((l_f_z - l_dnorm)[,i])}))
  llh <- llh_cliques - llh_seps
  return(llh*(-1)^negative)
}

################################################################################
## Fit the Conditional extremes model with graphical structure and MVAGG residuals
llh_dmvagg_graph_Gamma <- function(data, loc, scale_1, scale_2, shape, Gamma_zero, cliques, seps, negative = FALSE){
  #extract info from data
  n <- dim(data)[1]
  d <- dim(data)[2]
  
  ## Transform data onto joint standard Gaussian margins
  Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i]))}))
  Sigma_start <- cor(Q_F_z)
  
  if(any(is.infinite(Q_F_z)) | any(is.na(Q_F_z)) | any(is.nan(Q_F_z))){
    return((-10^10)*(-1)^negative)
  }
  else if(any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
    return((-10^10)*(-1)^negative)
  }
  else{
    if(is_empty(Gamma_zero)){
      Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, nobs = n, penalize.diagonal = FALSE, thr = 1e-9))
    }
    else{
      Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, nobs = n, zero = as.matrix(Gamma_zero), penalize.diagonal = FALSE, thr = 1e-9))
    }
    Sigma <- Lasso_est$w
    
    ## Likelihood calculation
    l_f_z <- sapply(1:d, function(i){daggd(x = data[,i], loc = loc[i], scale_1 = scale_1[i], scale_2 = scale_2[i], shape = shape[i], log = TRUE)})
    l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
    
    l_mvnorm_cliques <- sapply(cliques, function(i){sum(mvtnorm::dmvnorm(x = as.matrix(Q_F_z[,i]), mean = rep(0, length(i)), sigma = as.matrix(Sigma[i,i]), log = TRUE))})
    llh <- sum(l_mvnorm_cliques + sapply(cliques, function(i){sum((l_f_z - l_dnorm)[,i])}))
    
    if(!is_empty(seps)){
      l_mvnorm_seps  <- sapply(seps, function(i){sum(mvtnorm::dmvnorm(x = as.matrix(Q_F_z[,i]), mean = rep(0, length(i)), sigma = as.matrix(Sigma[i,i]), log = TRUE))})
      
      llh_seps <- sum(l_mvnorm_seps + sapply(seps, function(i){sum((l_f_z - l_dnorm)[,i])}))
      llh <- llh - llh_seps 
    }
    return(llh*(-1)^negative)
  }
}

## Conditional extremes model assuming the residuals follow a (multivariate) asymmetric generalised Gaussian distribution
## Gamma is estimated 
Cond_extremes_MVAGGD_Gamma <- function(data, cond, graph = NA, start,
                                       maxit = 1e+6, nOptim = 1){
  
  theCall <- match.call()
  # if (!inherits(x, "migpd")) 
  #   stop("you need to use an object created by migpd")
  # margins <- list(casefold(margins), p2q = switch(casefold(margins), 
  #                                                 gumbel = function(p) -log(-log(p)), laplace = function(p) ifelse(p < 
  #                                                                                                                    0.5, log(2 * p), -log(2 * (1 - p)))), q2p = switch(casefold(margins), 
  #                                                                                                                                                                       gumbel = function(q) exp(-exp(-q)), laplace = function(q) ifelse(q < 
  #                                                                                                                                                                                                                                          0, exp(q)/2, 1 - 0.5 * exp(-q))))
  # x <- mexTransform(x, margins = margins, method = marTransform, 
  #                   r = referenceMargin)
  # x$referenceMargin <- referenceMargin
  # if (margins[[1]] == "gumbel" & constrain) {
  #   warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
  #   constrain <- FALSE
  # }
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(is.character(cond)){ 
    cond <- match(cond, dimnames(x$transformed)[[2]])
  }
  # if(missing(dqu)){
  #   message("Assuming same quantile for dependence thesholding as was used \n     to fit corresponding marginal model...\n")
  #   dqu <- x$mqu[cond]
  # }
  # dth <- quantile(x$transformed[, cond], dqu)
  dependent <- (1:(dim(data)[[2]]))[-cond]
  # if(length(dqu) < length(dependent)){
  #   dqu <- rep(dqu, length = length(dependent))
  #   aLow <- if else(margins[[1]] == "gumbel", 10^(-10), -1 + 10^(-10))
  # }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    non_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graph_cond <- delete.vertices(graph, cond)
    graphs_ind <- lapply(1:d, delete.vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    #get the cliques and the separators in the graph
    cliques <- order_cliques(get_cliques_and_separators(graph_cond))
    #now get the separators but we only include singletons a vertex only appears once in a clique of length two
    seps <- get_separators(cliques, includeSingletons = inc_sings)
    #check there is no cross over. If there is remove from the clique set
    index <- which(cliques %in% intersect(cliques, seps))
    if(length(index) != 0){
      cliques <- cliques[-index]
    }
    
    #sort the cliques and separators into ascending order
    cliques <- order_cliques(cliques)
    
    cliques <- lapply(cliques, function(x){sort(x)})
    seps <- lapply(seps, function(x){sort(x)})
    
    ## Get the non_edges in the graph
    all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
    edges_in_graph <- data.frame(as_edgelist(graph_cond))
    all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
    non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  ## fit the model
  res <- qfun_MVAGGD_graph_Gamma(yex = yex, ydep = as.matrix(ydep), 
                                 comps = 1:(d-1), cliques = cliques, seps = seps, Gamma_zero = non_edges,
                                 maxit = maxit, start = c(start), nOptim = nOptim)
  
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(res$par$a, res$par$b, res$par$mu, res$par$sigma_1, res$par$sigma_2, res$par$shape)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out_llh <- data.frame(llh = -res$value)
  colnames(out_llh) <- c("Total")
  rownames(out_llh) <- c("llh")
  out$loglike <- out_llh
  
  out$convergence <- res$convergence
  
  out_Z <- res$Z
  out_Z <- as.data.frame(out_Z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
}

qfun_MVAGGD_graph_Gamma <- function(yex, ydep, comps, cliques, seps, Gamma_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, comps, cliques, seps, Gamma_zero, negative = FALSE){
    #get the starting parameters
    d_comps <- length(comps)
    n_data <- length(yex)
    a <- param[1:d_comps]
    b <- param[(d_comps + 1):(2*d_comps)]
    mu <- param[(2*d_comps + 1):(3*d_comps)]
    sigma_1 <- param[(3*d_comps + 1):(4*d_comps)]
    sigma_2 <- param[(4*d_comps + 1):(5*d_comps)]
    shape <- param[(5*d_comps + 1):(6*d_comps)]
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0) | any(abs(a) > 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      res <- llh_dmvagg_graph_Gamma(data = z, loc = mu, scale_1 = sigma_1, scale_2 = sigma_2,
                                    shape = shape, Gamma_zero = Gamma_zero, cliques = cliques, seps = seps) -
        sum(log(b_yi))
      
      #check the value is valid
      if(is.infinite(res)){
        if(res < 0) {
          res <- return((-10^10)*(-1)^negative)
        }
        else{ 
          res <- return((-10^10)*(-1)^negative)
        }
        warning("Infinite value of Q in mexDependence")
      }
    }
    return(res*(-1)^negative)
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   yex = yex, ydep = ydep, comps = comps,
                   cliques = cliques, seps = seps, Gamma_zero = Gamma_zero,
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = length(comps)),
                    b = rep(NA, times = length(comps)),
                    mu = rep(NA, times = length(comps)),
                    Sigma = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       yex = yex, ydep = ydep, comps = comps, Gamma_non_zero = Gamma_non_zero,
                       negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = length(comps)),
                        b = rep(NA, times = length(comps)),
                        mu = rep(NA, times = length(comps)),
                        Sigma = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d_comps <- length(comps)
    a_hat <- fit$par[1:d_comps]
    b_hat <- fit$par[(d_comps + 1):(2*d_comps)]
    mu_hat <- fit$par[(2*d_comps + 1):(3*d_comps)]
    sigma_1_hat <- fit$par[(3*d_comps + 1):(4*d_comps)]
    sigma_2_hat <- fit$par[(4*d_comps + 1):(5*d_comps)]
    shape_hat <- fit$par[(5*d_comps + 1):(6*d_comps)]
    
    #obtain the residuals
    a_yi_hat <- sapply(a_hat, function(a){yex*a})
    b_yi_hat <- sapply(b_hat, function(b){yex^b})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d_comps, function(i){qnorm(paggd(q = Z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    if(is_empty(Gamma_zero)){
      Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, thr = 1e-9))
    }
    else{
      Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, zero = as.matrix(Gamma_zero), penalize.diagonal = FALSE, thr = 1e-9))
    }
    Gamma_hat <- Lasso_est$wi
    
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

################################################################################
## As above but now we use the two-step method
## Take the fitted residuals as an input and then fit the MVAGG to the data with some graphical structure
## If IID or Saturated models are used use fit_agg or fit_mvagg

Cond_extremes_MVAGGD_Gamma_TS <- function(data, cond, graph = NA, start,
                                          constrain = TRUE, q = c(0,1), v = 10, aLow = -1,
                                          maxit = 1e+6, nOptim = 1){
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  dependent <- (1:(dim(data)[[2]]))[-cond]
  
  ## First step where we fit the original conditional extremes model
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  res_OG <- lapply(1:(d-1), function(i){
    qfun(yex = yex, ydep = as.matrix(ydep[,i]), 
         constrain = constrain, aLow = aLow, q = q, v = v,
         maxit = maxit, start = start[i,1:2], nOptim = nOptim)})
  
  ## Get the necessary output
  z <- sapply(res_OG, function(x){x$Z})
  a_hat <- sapply(res_OG, function(x){x$par$a})
  b_hat <- sapply(res_OG, function(x){x$par$b})
  
  ## Now fit the graphical model to Z|i
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    non_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graph_cond <- delete_vertices(graph, cond)
    graphs_ind <- lapply(1:d, delete.vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    #get the cliques and the separators in the graph
    cliques <- order_cliques(get_cliques_and_separators(graph_cond))
    #now get the separators but we only include singletons a vertex only appears once in a clique of length two
    seps <- get_separators(cliques, includeSingletons = inc_sings)
    #check there is no cross over. If there is remove from the clique set
    index <- which(cliques %in% intersect(cliques, seps))
    if(length(index) != 0){
      cliques <- cliques[-index]
    }
    
    #sort the cliques and separators into ascending order
    cliques <- order_cliques(cliques)
    
    cliques <- lapply(cliques, function(x){sort(x)})
    seps <- lapply(seps, function(x){sort(x)})
    
    ## Get the non_edges in the graph
    all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
    edges_in_graph <- data.frame(as_edgelist(graph_cond))
    all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
    non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
  }
  
  ## fit the model
  res <- qfun_MVAGGD_graph_Gamma_TS(z = z, comps = 1:(d-1), 
                                    cliques = cliques, seps = seps, Gamma_zero = non_edges,
                                    maxit = maxit, start = c(start[,-c(1:2)]), nOptim = nOptim)
  
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(a_hat, b_hat, res$par$mu, res$par$sigma_1, res$par$sigma_2, res$par$shape)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out_llh <- data.frame(llh = -res$value)
  colnames(out_llh) <- c("Total")
  rownames(out_llh) <- c("llh")
  out$loglike <- out_llh
  
  out$convergence <- res$convergence
  
  out_Z <- as.data.frame(z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
}

qfun_MVAGGD_graph_Gamma_TS <- function(z, comps, cliques, seps, Gamma_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, comps, cliques, seps, Gamma_zero, negative = FALSE){
    #get the starting parameters
    d_comps <- length(comps)
    mu <- param[1:d_comps]
    sigma_1 <- param[(d_comps + 1):(2*d_comps)]
    sigma_2 <- param[(2*d_comps + 1):(3*d_comps)]
    shape <- param[(3*d_comps + 1):(4*d_comps)]
    
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      res <- llh_dmvagg_graph_Gamma(data = z, loc = mu, scale_1 = sigma_1, scale_2 = sigma_2,
                                    shape = shape, Gamma_zero = Gamma_zero, cliques = cliques, seps = seps)
      
      #check the value is valid
      if(is.infinite(res)){
        if(res < 0) {
          res <- return((-10^10)*(-1)^negative)
        }
        else{ 
          res <- return((-10^10)*(-1)^negative)
        }
        warning("Infinite value of Q in mexDependence")
      }
    }
    return(res*(-1)^negative)
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, comps = comps,
                   cliques = cliques, seps = seps, Gamma_zero = Gamma_zero,
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = length(comps)),
                    sigma_1 = rep(NA, times = length(comps)),
                    sigma_2 = rep(NA, times = length(comps)),
                    shape = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = length(comps)),
                    sigma_1 = rep(NA, times = length(comps)),
                    sigma_2 = rep(NA, times = length(comps)),
                    shape = diag(NA, length(comps)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, comps = comps,
                       cliques = cliques, seps = seps, Gamma_zero = Gamma_zero,
                       negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = length(comps)),
                        sigma_1 = rep(NA, times = length(comps)),
                        sigma_2 = rep(NA, times = length(comps)),
                        shape = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = length(comps)),
                        sigma_1 = rep(NA, times = length(comps)),
                        sigma_2 = rep(NA, times = length(comps)),
                        shape = diag(NA, length(comps)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d_comps <- length(comps)
    mu_hat <- fit$par[1:d_comps]
    sigma_1_hat <- fit$par[(d_comps + 1):(2*d_comps)]
    sigma_2_hat <- fit$par[(2*d_comps + 1):(3*d_comps)]
    shape_hat <- fit$par[(3*d_comps + 1):(4*d_comps)]
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d_comps, function(i){qnorm(paggd(q = z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    if(is_empty(Gamma_zero)){
      Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, thr = 1e-9))
    }
    else{
      Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, zero = as.matrix(Gamma_zero), penalize.diagonal = FALSE, thr = 1e-9))
    }
    Gamma_hat <- Lasso_est$wi
    
    out <- list()
    out$par <- list(mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

################################################################################
## Fit the two-step model
## Take the fitted residuals as an input and then fit the MVAGG to the data with some graphical structure
## If IID or Saturated models are used use fit_agg or fit_mvagg

## under the graphical model
fit_MVAGG_graph_two_step <- function(Z, cond, graph, start = NULL){
  
  #Get dimension of the data
  dim_data <- dim(Z)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2] + 1
    n <- dim_data[1]
  }
  dependent <- (1:d)[-cond]
  
  #check the graphical structure
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graph_cond <- delete_vertices(graph, cond)
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    #get the cliques and the separators in the graph
    cliques <- order_cliques(get_cliques_and_separators(graph_cond))
    #now get the separators but we only include singletons a vertex only appears once in a clique of length two
    seps <- get_separators(cliques, includeSingletons = inc_sings)
    #check there is no cross over. If there is remove from the clique set
    index <- which(cliques %in% intersect(cliques, seps))
    if(length(index) != 0){
      cliques <- cliques[-index]
    }
    
    #sort the cliques and separators into ascending order
    cliques <- order_cliques(cliques)
    
    cliques <- lapply(cliques, function(x){sort(x)})
    seps <- lapply(seps, function(x){sort(x)})
  }
  
  ## get the starting parameters
  if(is.null(start)){
    start <- cbind(apply(Z, 2, quantile, probs = 0.5),
                   apply(Z, 2, var),
                   apply(Z, 2, var),
                   rep(1.5, d-1))
  }
  
  if(is_empty(seps)){
    fit_res <- vector("list", length(cliques))
    for(i in 1:length(cliques)){
      if(length(cliques[[i]]) == 1){
        Sigma_start <- 1
      }
      else{
        Sigma_start <- cor(Z[,cliques[[i]]])
        Sigma_start <- Sigma_start[upper.tri(Sigma_start)] 
      }
      fit_res[[i]] <- fit_mvagg(data = as.matrix(Z[,cliques[[i]]]),
                                par = c(c(start[cliques[[i]],]), Sigma_start))
      # }
    }
    
    #get the output
    out <- list()
    
    if(length(fit_res) == 1){
      out$par$main <- fit_res[[1]]$par$main
    }
    else{
      out$par$main <- do.call(cbind, lapply(fit_res, function(x){x$par$main}))
    }
    rownames(out$par$main) <- c("loc", "scale_1", "scale_2", "shape")
    colnames(out$par$main) <- sapply(dependent, function(x){paste("Column", x)})
    
    out$loglike <- sapply(fit_res, function(x){x$llh})
    out$convergence <- sapply(fit_res, function(x){x$convergence})
    
    out$par$Sigma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:length(cliques)){
      Sig_entries <- permutations(n = length(cliques[[i]]), r = 2, v = cliques[[i]], repeats.allowed = TRUE)
      out$par$Sigma[Sig_entries] <- fit_res[[i]]$par$Sigma
    }
    return(out)
  }
  else{
    ## Get the non-zero elements of Sigma
    Sigma_Z <- cor(Z)
    Sig_non_zero <- upper.tri(Sigma_Z)
    ## Get the starting parameters
    Sigma_start <- Sigma_Z[Sig_non_zero]
    
    ## Fit the model
    fit_res <- fit_mvagg_graph(data = as.matrix(Z), 
                               par = c(c(start), Sigma_start),
                               cliques = cliques, 
                               seps = seps,
                               Sig_non_zero = Sig_non_zero)
    
    ## Get the output
    out <- list()
    
    out$par$main <- fit_res$par$main
    rownames(out$par$main) <- c("loc", "scale_1", "scale_2", "shape")
    colnames(out$par$main) <- sapply(dependent, function(x){paste("Column", x)})
    
    out$par$Sigma <- fit_res$par$Sigma
    
    out$loglike <- fit_res$llh
    out$convergence <- fit_res$convergence
    
    return(out)
  }
}

################################################################################
Cond_extremes_MVAGGD_Gamma_new <- function(data, cond, graph = NA, start,
                                           maxit = 1e+6, nOptim = 1){
  
  theCall <- match.call()
  # if (!inherits(x, "migpd")) 
  #   stop("you need to use an object created by migpd")
  # margins <- list(casefold(margins), p2q = switch(casefold(margins), 
  #                                                 gumbel = function(p) -log(-log(p)), laplace = function(p) ifelse(p < 
  #                                                                                                                    0.5, log(2 * p), -log(2 * (1 - p)))), q2p = switch(casefold(margins), 
  #                                                                                                                                                                       gumbel = function(q) exp(-exp(-q)), laplace = function(q) ifelse(q < 
  #                                                                                                                                                                                                                                          0, exp(q)/2, 1 - 0.5 * exp(-q))))
  # x <- mexTransform(x, margins = margins, method = marTransform, 
  #                   r = referenceMargin)
  # x$referenceMargin <- referenceMargin
  # if (margins[[1]] == "gumbel" & constrain) {
  #   warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
  #   constrain <- FALSE
  # }
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(is.character(cond)){ 
    cond <- match(cond, dimnames(x$transformed)[[2]])
  }
  # if(missing(dqu)){
  #   message("Assuming same quantile for dependence thesholding as was used \n     to fit corresponding marginal model...\n")
  #   dqu <- x$mqu[cond]
  # }
  # dth <- quantile(x$transformed[, cond], dqu)
  dependent <- (1:(dim(data)[[2]]))[-cond]
  # if(length(dqu) < length(dependent)){
  #   dqu <- rep(dqu, length = length(dependent))
  #   aLow <- if else(margins[[1]] == "gumbel", 10^(-10), -1 + 10^(-10))
  # }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    non_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    
    ## Now figure out the cliques and separators
    graph_cond <- delete_vertices(graph, cond)
    comps <- components(graph_cond)
    
    if(comps$no != 1){
      ## Get the separate components
      v_comps <- groups(comps)
      g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
      
      ## Get the non_edges in the graph(s)
      non_edges <- vector("list", comps$no)
      for(i in 1:comps$no){
        if(length(V(g_comps[[i]])) <= 2){
          non_edges[[i]] <- matrix(NA, nrow = 0, ncol = 2)
        }
        else{
          all_edges <- data.frame(combinations(n = length(V(g_comps[[i]])), r = 2, v = V(g_comps[[i]])))
          edges_in_comp <- data.frame(as_edgelist(g_comps[[i]]))
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_comp)
          non_edges[[i]] <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        }
      }
      
      ## Get the cliques and separators
      for(i in 1:comps$no){
        V(g_comps[[i]])$name <- v_comps[[i]]
      }
      cliques <- lapply(g_comps, function(x){max_cliques(x)})
      seps <- lapply(cliques, function(x){get_separators(x, includeSingletons = inc_sings)})
      
    }
    else{
      #get the cliques and the separators in the graph
      cliques <- max_cliques(graph_cond)
      #now get the separators but we only include singletons a vertex only appears once in a clique of length two
      seps <- get_separators(cliques, includeSingletons = inc_sings)
      
      ## Get the non_edges in the graph
      all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
      edges_in_graph <- data.frame(as_edgelist(graph_cond))
      all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
      non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
    }
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  ## fit the models and get the output
  if(is_empty(seps) & length(cliques) == d-1){
    ## Make this function and do them independently of one another
    res_1 <- lapply(1:(d-1), function(i){
      qfun_MVAGGD_Gamma_indep(yex = yex, ydep = as.matrix(ydep[,i]),
                              maxit = maxit, start = c(start[i,]), nOptim = nOptim)})
    
    res <- list()
    res$par$a <- do.call(c, lapply(res_1, function(x){x$par$a}))
    res$par$b <- do.call(c, lapply(res_1, function(x){x$par$b}))
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- diag(1, d-1)
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
    res$Z <- do.call(cbind, lapply(res_1, function(x){x$Z}))
  }
  else if(is_empty(seps) & length(cliques) == 1){
    res <- qfun_MVAGGD_Gamma_full(yex = yex, ydep = as.matrix(ydep), 
                                  maxit = maxit, start = c(start), nOptim = nOptim)
  }
  else if(!is_connected(graph_cond)){
    res_1 <- vector("list", comps$no)
    for(i in 1:comps$no){
      if(is_empty(seps[[i]])){
        res_1[[i]] <- qfun_MVAGGD_Gamma_full(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]),
                                             maxit = maxit, start = c(start[v_comps[[i]],]), nOptim = nOptim)
      }
      else{
        res_1[[i]] <- qfun_MVAGGD_Gamma_graph(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]), Gamma_zero = non_edges[[i]],
                                              maxit = maxit, start = c(start[v_comps[[i]],]), nOptim = nOptim)
      }
    }
    
    res <- list()
    res$par$a <- do.call(c, lapply(res_1, function(x){x$par$a}))
    res$par$b <- do.call(c, lapply(res_1, function(x){x$par$b}))
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:comps$no){
      if(length(v_comps[[i]]) == 1){
        clique_edges <- matrix(rep(v_comps[[i]], 2), ncol = 2)
      }
      else{
        clique_edges <- permutations(n = length(v_comps[[i]]), r = 2, v = v_comps[[i]], repeats.allowed = TRUE)
      }
      res$par$Gamma[clique_edges] <- res_1[[i]]$par$Gamma
    }
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
    res$Z <- do.call(cbind, lapply(res_1, function(x){x$Z}))
  }
  else{
    res <- qfun_MVAGGD_Gamma_graph(yex = yex, ydep = as.matrix(ydep), Gamma_zero = non_edges,
                                   maxit = maxit, start = c(start), nOptim = nOptim)
  }
  
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(res$par$a, res$par$b, res$par$mu, res$par$sigma_1, res$par$sigma_2, res$par$shape)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out$loglike <- -sum(res$value)
  
  out$convergence <- max(res$convergence)
  
  out_Z <- res$Z
  out_Z <- as.data.frame(out_Z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
}

Cond_extremes_MVAGGD_Gamma_decomp <- function(data, cond, graph = NA, start,
                                              maxit = 1e+6, nOptim = 1){
  
  theCall <- match.call()
  # if (!inherits(x, "migpd")) 
  #   stop("you need to use an object created by migpd")
  # margins <- list(casefold(margins), p2q = switch(casefold(margins), 
  #                                                 gumbel = function(p) -log(-log(p)), laplace = function(p) ifelse(p < 
  #                                                                                                                    0.5, log(2 * p), -log(2 * (1 - p)))), q2p = switch(casefold(margins), 
  #                                                                                                                                                                       gumbel = function(q) exp(-exp(-q)), laplace = function(q) ifelse(q < 
  #                                                                                                                                                                                                                                          0, exp(q)/2, 1 - 0.5 * exp(-q))))
  # x <- mexTransform(x, margins = margins, method = marTransform, 
  #                   r = referenceMargin)
  # x$referenceMargin <- referenceMargin
  # if (margins[[1]] == "gumbel" & constrain) {
  #   warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
  #   constrain <- FALSE
  # }
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(is.character(cond)){ 
    cond <- match(cond, dimnames(x$transformed)[[2]])
  }
  # if(missing(dqu)){
  #   message("Assuming same quantile for dependence thesholding as was used \n     to fit corresponding marginal model...\n")
  #   dqu <- x$mqu[cond]
  # }
  # dth <- quantile(x$transformed[, cond], dqu)
  dependent <- (1:(dim(data)[[2]]))[-cond]
  # if(length(dqu) < length(dependent)){
  #   dqu <- rep(dqu, length = length(dependent))
  #   aLow <- if else(margins[[1]] == "gumbel", 10^(-10), -1 + 10^(-10))
  # }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    non_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    
    ## Now figure out the cliques and separators
    graph_cond <- delete_vertices(graph, cond)
    comps <- components(graph_cond)
    
    if(comps$no != 1){
      ## Get the separate components
      v_comps <- groups(comps)
      g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
      
      ## Get the non_edges in the graph(s)
      non_edges <- vector("list", comps$no)
      for(i in 1:comps$no){
        if(length(V(g_comps[[i]])) <= 2){
          non_edges[[i]] <- matrix(NA, nrow = 0, ncol = 2)
        }
        else{
          all_edges <- data.frame(combinations(n = length(V(g_comps[[i]])), r = 2, v = V(g_comps[[i]])))
          edges_in_comp <- data.frame(as_edgelist(g_comps[[i]]))
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_comp)
          non_edges[[i]] <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        }
      }
      
      ## Get the cliques and separators
      for(i in 1:comps$no){
        V(g_comps[[i]])$name <- v_comps[[i]]
      }
      cliques <- lapply(g_comps, function(x){max_cliques(x)})
      seps <- lapply(cliques, function(x){get_separators(x, includeSingletons = inc_sings)})
      
    }
    else{
      #get the cliques and the separators in the graph
      cliques <- max_cliques(graph_cond)
      #now get the separators but we only include singletons a vertex only appears once in a clique of length two
      seps <- get_separators(cliques, includeSingletons = inc_sings)
      
      ## Get the non_edges in the graph
      all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
      edges_in_graph <- data.frame(as_edgelist(graph_cond))
      all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
      non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
    }
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  ## fit the models and get the output
  if(is_empty(seps) & length(cliques) == d-1){
    ## Make this function and do them independently of one another
    res_1 <- lapply(1:(d-1), function(i){
      qfun_MVAGGD_Gamma_indep(yex = yex, ydep = as.matrix(ydep[,i]),
                              maxit = maxit, start = c(start[i,]), nOptim = nOptim)})
    
    res <- list()
    res$par$a <- do.call(c, lapply(res_1, function(x){x$par$a}))
    res$par$b <- do.call(c, lapply(res_1, function(x){x$par$b}))
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- diag(1, d-1)
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
    res$Z <- do.call(cbind, lapply(res_1, function(x){x$Z}))
  }
  else if(is_empty(seps) & length(cliques) == 1){
    res <- qfun_MVAGGD_Gamma_full(yex = yex, ydep = as.matrix(ydep), 
                                  maxit = maxit, start = c(start), nOptim = nOptim)
  }
  else if(!is_connected(graph_cond)){
    res_1 <- vector("list", comps$no)
    for(i in 1:comps$no){
      if(is_empty(seps[[i]])){
        res_1[[i]] <- qfun_MVAGGD_Gamma_full(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]),
                                             maxit = maxit, start = c(start[v_comps[[i]],]), nOptim = nOptim)
      }
      else{
        res_1[[i]] <- qfun_MVAGGD_Gamma_graph_decomp(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]),
                                                     cliques = cliques[[i]], seps = seps[[i]], Gamma_zero = non_edges[[i]],
                                                     maxit = maxit, start = c(start[v_comps[[i]],]), nOptim = nOptim)
      }
    }
    
    res <- list()
    res$par$a <- do.call(c, lapply(res_1, function(x){x$par$a}))
    res$par$b <- do.call(c, lapply(res_1, function(x){x$par$b}))
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:comps$no){
      if(length(v_comps[[i]]) == 1){
        clique_edges <- matrix(rep(v_comps[[i]], 2), ncol = 2)
      }
      else{
        clique_edges <- permutations(n = length(v_comps[[i]]), r = 2, v = v_comps[[i]], repeats.allowed = TRUE)
      }
      res$par$Gamma[clique_edges] <- res_1[[i]]$par$Gamma
    }
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
    res$Z <- do.call(cbind, lapply(res_1, function(x){x$Z}))
  }
  else{
    res <- qfun_MVAGGD_Gamma_graph_decomp(yex = yex, ydep = as.matrix(ydep),
                                          cliques = cliques, seps = seps, Gamma_zero = non_edges,
                                          maxit = maxit, start = c(start), nOptim = nOptim)
  }
  
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(res$par$a, res$par$b, res$par$mu, res$par$sigma_1, res$par$sigma_2, res$par$shape)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out$loglike <- -sum(res$value)
  
  out$convergence <- max(res$convergence)
  
  out_Z <- res$Z
  out_Z <- as.data.frame(out_Z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
}

qfun_MVAGGD_Gamma_indep <- function(yex, ydep, maxit, start, nOptim){
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, negative = FALSE){
    #get the starting parameters
    d <- ncol(ydep)
    n <- length(yex)
    a <- param[1:d]
    b <- param[(d + 1):(2*d)]
    mu <- param[(2*d + 1):(3*d)]
    sigma_1 <- param[(3*d + 1):(4*d)]
    sigma_2 <- param[(4*d + 1):(5*d)]
    shape <- param[(5*d + 1):(6*d)]
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0) | any(abs(a) > 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## get the Z data
      b_yi <- yex^b
      z <- (ydep - yex*a)/b_yi
      res <- sum(daggd(x = z, loc = mu, scale_1 = sigma_1, scale_2 = sigma_2, shape = shape, log = TRUE)) - sum(log(b_yi))
      
      #check the value is valid
      if(is.infinite(res)){
        res <- return((-10^10)*(-1)^negative)
        warning("Infinite value of Q in mexDependence")
      }
      return(res*(-1)^negative)
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                   control = list(maxit = maxit),
                   negative = TRUE, method = "Nelder-Mead", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
  }
  else if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    a_hat <- fit$par[1]
    b_hat <- fit$par[2]
    
    #obtain the residuals
    Z <- (ydep - yex*a_hat)/(yex^b_hat)
    
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = fit$par[3], sigma_1 = fit$par[4], sigma_2 = fit$par[5], shape = fit$par[6])
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    out <- list()
  }
  return(out)
}

qfun_MVAGGD_Gamma_graph <- function(yex, ydep, Gamma_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, Gamma_zero, negative = FALSE){
    #get the starting parameters
    d <- ncol(ydep)
    n <- length(yex)
    a <- param[1:d]
    b <- param[(d + 1):(2*d)]
    mu <- param[(2*d + 1):(3*d)]
    sigma_1 <- param[(3*d + 1):(4*d)]
    sigma_2 <- param[(4*d + 1):(5*d)]
    shape <- param[(5*d + 1):(6*d)]
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0) | any(abs(a) > 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## get the Z data
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      ## get this on standard Gaussian margins
      Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i]))}))
      Sigma_start <- cor(Q_F_z)
      if(any(is.infinite(Q_F_z)) | any(is.na(Q_F_z)) | any(is.nan(Q_F_z)) |
         any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
        return((-10^10)*(-1)^negative)
      }
      else if(any(diag(Sigma_start) <= 0) | any(eigen(Sigma_start)$values <= 0)){
        return((-10^10)*(-1)^negative)
      }
      else{
        ## Calculate the log likelihood function
        l_f_z <- sapply(1:d, function(i){daggd(x = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i], log = TRUE)})
        l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
        
        Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = FALSE, zero = Gamma_zero, thr = 1e-9))
        chol_Gamma <- chol(Lasso_est$wi)
        y <- chol_Gamma%*%t(Q_F_z)
        l_mvnorm <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
        
        res <- l_mvnorm + sum(l_f_z) - sum(l_dnorm) - sum(log(b_yi))
        
        #check the value is valid
        if(is.infinite(res)){
          res <- return((-10^10)*(-1)^negative)
          warning("Infinite value of Q in mexDependence")
        }
        return(res*(-1)^negative)
      }
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                   control = list(maxit = maxit), Gamma_zero = Gamma_zero,
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
  }
  else if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d <- ncol(ydep)
    n <- length(yex)
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[(d + 1):(2*d)]
    mu_hat <- fit$par[(2*d + 1):(3*d)]
    sigma_1_hat <- fit$par[(3*d + 1):(4*d)]
    sigma_2_hat <- fit$par[(4*d + 1):(5*d)]
    shape_hat <- fit$par[(5*d + 1):(6*d)]
    
    #obtain the residuals
    a_yi_hat <- sapply(a_hat, function(a){yex*a})
    b_yi_hat <- sapply(b_hat, function(b){yex^b})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = Z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    
    Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, zero = Gamma_zero, thr = 1e-9))
    Gamma_hat <- Lasso_est$wi
    
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    out <- list()
  }
  return(out)
}

qfun_MVAGGD_Gamma_graph_decomp <- function(yex, ydep, cliques, seps, Gamma_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, cliques, seps, Gamma_zero, negative = FALSE){
    #get the starting parameters
    d <- ncol(ydep)
    n <- length(yex)
    a <- param[1:d]
    b <- param[(d + 1):(2*d)]
    mu <- param[(2*d + 1):(3*d)]
    sigma_1 <- param[(3*d + 1):(4*d)]
    sigma_2 <- param[(4*d + 1):(5*d)]
    shape <- param[(5*d + 1):(6*d)]
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0) | any(abs(a) > 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## get the Z data
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      ## get this on standard Gaussian margins
      Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i]))}))
      Sigma_start <- cor(Q_F_z)
      if(any(is.infinite(Q_F_z)) | any(is.na(Q_F_z)) | any(is.nan(Q_F_z)) |
         any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
        return((-10^10)*(-1)^negative)
      }
      else if(any(diag(Sigma_start) <= 0) | any(eigen(Sigma_start)$values <= 0)){
        return((-10^10)*(-1)^negative)
      }
      else{
        l_f_z <- sapply(1:d, function(i){daggd(x = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i], log = TRUE)})
        l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
        l_diff <- l_f_z - l_dnorm
        
        Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, zero = Gamma_zero, thr = 1e-9))
        chol_Gamma <- chol(Lasso_est$wi)
        chol_Gamma_cliques <- lapply(cliques, function(i){as.matrix(chol_Gamma[i,i])})
        chol_Gamma_seps <- lapply(seps, function(i){as.matrix(chol_Gamma[i,i])})
        
        n_cliques <- length(cliques)
        n_seps <- length(seps)
        y_cliques <- lapply(1:n_cliques, function(i){chol_Gamma_cliques[[i]]%*%t(Q_F_z[,cliques[[i]]])})
        y_seps <- lapply(1:n_seps, function(i){chol_Gamma_seps[[i]]%*%t(Q_F_z[,seps[[i]]])})
        
        l_mvnorm_cliques <- sapply(1:n_cliques, function(i){
          -n*length(cliques[[i]])*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma_cliques[[i]]))) - sum(y_cliques[[i]]^2)/2
        })
        
        l_mvnorm_seps <- sapply(1:n_seps, function(i){
          -n*length(seps[[i]])*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma_seps[[i]]))) - sum(y_seps[[i]]^2)/2
        })
        
        res <- sum(l_mvnorm_cliques) - sum(l_mvnorm_seps) + sum(l_diff) - sum(log(b_yi))
        
        #check the value is valid
        if(is.infinite(res)){
          res <- return((-10^10)*(-1)^negative)
          warning("Infinite value of Q in mexDependence")
        }
        return(res*(-1)^negative)
      }
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                   control = list(maxit = maxit),
                   cliques = cliques, seps = seps, Gamma_zero = Gamma_zero,
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
  }
  else if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d <- ncol(ydep)
    n <- length(yex)
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[(d + 1):(2*d)]
    mu_hat <- fit$par[(2*d + 1):(3*d)]
    sigma_1_hat <- fit$par[(3*d + 1):(4*d)]
    sigma_2_hat <- fit$par[(4*d + 1):(5*d)]
    shape_hat <- fit$par[(5*d + 1):(6*d)]
    
    #obtain the residuals
    a_yi_hat <- sapply(a_hat, function(a){yex*a})
    b_yi_hat <- sapply(b_hat, function(b){yex^b})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = Z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    
    Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, zero = Gamma_zero, thr = 1e-9))
    Gamma_hat <- Lasso_est$wi
    
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    out <- list()
  }
  return(out)
}

qfun_MVAGGD_Gamma_full <- function(yex, ydep, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, negative = FALSE){
    #get the starting parameters
    d <- ncol(ydep)
    n <- length(yex)
    a <- param[1:d]
    b <- param[(d + 1):(2*d)]
    mu <- param[(2*d + 1):(3*d)]
    sigma_1 <- param[(3*d + 1):(4*d)]
    sigma_2 <- param[(4*d + 1):(5*d)]
    shape <- param[(5*d + 1):(6*d)]
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0) | any(abs(a) > 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## get the Z data
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      ## get this on standard Gaussian margins
      Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i]))}))
      Sigma_start <- cor(Q_F_z)
      if(any(is.infinite(Q_F_z)) | any(is.na(Q_F_z)) | any(is.nan(Q_F_z)) |
         any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
        return((-10^10)*(-1)^negative)
      }
      else if(any(diag(Sigma_start) <= 0) | any(eigen(Sigma_start)$values <= 0)){
        return((-10^10)*(-1)^negative)
      }
      else{
        l_f_z <- sapply(1:d, function(i){daggd(x = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i], log = TRUE)})
        l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
        
        Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = FALSE, thr = 1e-9))
        chol_Gamma <- chol(Lasso_est$wi)
        y <- chol_Gamma%*%t(Q_F_z)
        l_mvnorm <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
        
        res <- l_mvnorm + sum(l_f_z) - sum(l_dnorm) - sum(log(b_yi))
        
        #check the value is valid
        if(is.infinite(res)){
          res <- return((-10^10)*(-1)^negative)
          warning("Infinite value of Q in mexDependence")
        }
        return(res*(-1)^negative)
      }
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                   control = list(maxit = maxit),
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
  }
  else if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d <- ncol(ydep)
    n <- length(yex)
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[(d + 1):(2*d)]
    mu_hat <- fit$par[(2*d + 1):(3*d)]
    sigma_1_hat <- fit$par[(3*d + 1):(4*d)]
    sigma_2_hat <- fit$par[(4*d + 1):(5*d)]
    shape_hat <- fit$par[(5*d + 1):(6*d)]
    
    #obtain the residuals
    a_yi_hat <- sapply(a_hat, function(a){yex*a})
    b_yi_hat <- sapply(b_hat, function(b){yex^b})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = Z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    
    Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, thr = 1e-9))
    Gamma_hat <- Lasso_est$wi
    
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    out <- list()
  }
  return(out)
}

Cond_extremes_MVAGGD_Gamma_TS_new <- function(data, cond, graph = NA, start,
                                              constrain = TRUE, q = c(0,1), v = 10, aLow = -1,
                                              maxit = 1e+6, nOptim = 1){
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  dependent <- (1:(dim(data)[[2]]))[-cond]
  
  ## First step where we fit the original conditional extremes model
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  res_OG <- lapply(1:(d-1), function(i){
    qfun_HT(yex = yex, ydep = as.matrix(ydep[,i]),
            constrain = constrain, aLow = aLow, q = q, v = v,
            maxit = maxit, start = start[i,1:2], nOptim = nOptim)})
  
  ## Get the necessary output
  z <- sapply(res_OG, function(x){x$Z})
  a_hat <- sapply(res_OG, function(x){x$par$a})
  b_hat <- sapply(res_OG, function(x){x$par$b})
  
  ## Now fit the graphical model to Z|i
  ## determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    non_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    
    ## Now figure out the cliques and separators
    graph_cond <- delete_vertices(graph, cond)
    comps <- components(graph_cond)
    
    if(comps$no != 1){
      ## Get the separate components
      v_comps <- groups(comps)
      g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
      
      ## Get the non_edges in the graph(s)
      non_edges <- vector("list", comps$no)
      for(i in 1:comps$no){
        if(length(V(g_comps[[i]])) <= 2){
          non_edges[[i]] <- matrix(NA, nrow = 0, ncol = 2)
        }
        else{
          all_edges <- data.frame(combinations(n = length(V(g_comps[[i]])), r = 2, v = V(g_comps[[i]])))
          edges_in_comp <- data.frame(as_edgelist(g_comps[[i]]))
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_comp)
          non_edges[[i]] <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        }
      }
      
      ## Get the cliques and separators
      for(i in 1:comps$no){
        V(g_comps[[i]])$name <- v_comps[[i]]
      }
      cliques <- lapply(g_comps, function(x){max_cliques(x)})
      seps <- lapply(cliques, function(x){get_separators(x, includeSingletons = inc_sings)})
      
    }
    else{
      #get the cliques and the separators in the graph
      cliques <- max_cliques(graph_cond)
      #now get the separators but we only include singletons a vertex only appears once in a clique of length two
      seps <- get_separators(cliques, includeSingletons = inc_sings)
      
      ## Get the non_edges in the graph
      all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
      edges_in_graph <- data.frame(as_edgelist(graph_cond))
      all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
      non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
    }
  }
  
  n_cliques <- length(cliques)
  ## fit the models and get the output
  if(is_empty(seps) & n_cliques == d-1){
    ## Make this function and do them independently of one another
    res_1 <- lapply(1:(d-1), function(i){
      qfun_MVAGGD_Gamma_indep_TS(z = as.matrix(z[,i]), maxit = maxit, start = c(start[i,-c(1:2)]), nOptim = nOptim)})
    
    res <- list()
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- diag(1, d-1)
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
  }
  else if(is_empty(seps) & n_cliques == 1){
    res <- qfun_MVAGGD_Gamma_full_TS(z = z, maxit = maxit, start = c(start[,-c(1:2)]), nOptim = nOptim)
  }
  else if(!is_connected(graph_cond)){
    res_1 <- vector("list", comps$no)
    for(i in 1:comps$no){
      if(is_empty(seps[[i]])){
        res_1[[i]] <- qfun_MVAGGD_Gamma_full_TS(z = as.matrix(z[,v_comps[[i]]]), 
                                                maxit = maxit, start = c(start[v_comps[[i]],-c(1:2)]), nOptim = nOptim)
      }
      else{
        res_1[[i]] <- qfun_MVAGGD_Gamma_graph_TS(z = z[,v_comps[[i]]], Gamma_zero = non_edges[[i]],
                                                 maxit = maxit, start = c(start[v_comps[[i]],-c(1:2)]), nOptim = nOptim)
      }
    }
    
    res <- list()
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:comps$no){
      if(length(v_comps[[i]]) == 1){
        clique_edges <- matrix(rep(v_comps[[i]], 2), ncol = 2)
      }
      else{
        clique_edges <- permutations(n = length(v_comps[[i]]), r = 2, v = v_comps[[i]], repeats.allowed = TRUE)
      }
      res$par$Gamma[clique_edges] <- res_1[[i]]$par$Gamma
    }
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
  }
  else{
    res <- qfun_MVAGGD_Gamma_graph_TS(z = z, Gamma_zero = non_edges,
                                      maxit = maxit, start = c(start[,-c(1:2)]), nOptim = nOptim)
  }
  
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(a_hat, b_hat, res$par$mu, res$par$sigma_1, res$par$sigma_2, res$par$shape)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out$loglike <- -sum(res$value)
  
  out$convergence <- max(res$convergence)
  
  out_Z <- as.data.frame(z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
} 

Cond_extremes_MVAGGD_Gamma_TS_decomp <- function(data, cond, graph = NA, start,
                                                 constrain = TRUE, q = c(0,1), v = 10, aLow = -1,
                                                 maxit = 1e+6, nOptim = 1){
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  dependent <- (1:(dim(data)[[2]]))[-cond]
  
  ## First step where we fit the original conditional extremes model
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  res_OG <- lapply(1:(d-1), function(i){
    qfun_HT(yex = yex, ydep = as.matrix(ydep[,i]),
            constrain = constrain, aLow = aLow, q = q, v = v,
            maxit = maxit, start = start[i,1:2], nOptim = nOptim)})
  
  ## Get the necessary output
  z <- sapply(res_OG, function(x){x$Z})
  a_hat <- sapply(res_OG, function(x){x$par$a})
  b_hat <- sapply(res_OG, function(x){x$par$b})
  
  ## Now fit the graphical model to Z|i
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    non_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    
    ## Now figure out the cliques and separators
    graph_cond <- delete_vertices(graph, cond)
    comps <- components(graph_cond)
    
    if(comps$no != 1){
      ## Get the separate components
      v_comps <- groups(comps)
      g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
      
      ## Get the non_edges in the graph(s)
      non_edges <- vector("list", comps$no)
      for(i in 1:comps$no){
        if(length(V(g_comps[[i]])) <= 2){
          non_edges[[i]] <- matrix(NA, nrow = 0, ncol = 2)
        }
        else{
          all_edges <- data.frame(combinations(n = length(V(g_comps[[i]])), r = 2, v = V(g_comps[[i]])))
          edges_in_comp <- data.frame(as_edgelist(g_comps[[i]]))
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_comp)
          non_edges[[i]] <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        }
      }
      
      ## Get the cliques and separators
      for(i in 1:comps$no){
        V(g_comps[[i]])$name <- v_comps[[i]]
      }
      cliques <- lapply(g_comps, function(x){max_cliques(x)})
      seps <- lapply(cliques, function(x){get_separators(x, includeSingletons = inc_sings)})
      
    }
    else{
      #get the cliques and the separators in the graph
      cliques <- max_cliques(graph_cond)
      #now get the separators but we only include singletons a vertex only appears once in a clique of length two
      seps <- get_separators(cliques, includeSingletons = inc_sings)
      
      ## Get the non_edges in the graph
      all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
      edges_in_graph <- data.frame(as_edgelist(graph_cond))
      all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
      non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
    }
  }
  
  n_cliques <- length(cliques)
  ## fit the models and get the output
  if(is_empty(seps) & n_cliques == d-1){
    ## Make this function and do them independently of one another
    res_1 <- lapply(1:(d-1), function(i){
      qfun_MVAGGD_Gamma_indep_TS(z = as.matrix(z[,i]), maxit = maxit, start = c(start[i,-c(1:2)]), nOptim = nOptim)})
    
    res <- list()
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- diag(1, d-1)
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
  }
  else if(is_empty(seps) & n_cliques == 1){
    res <- qfun_MVAGGD_Gamma_full_TS(z = z, maxit = maxit, start = c(start[,-c(1:2)]), nOptim = nOptim)
  }
  else if(!is_connected(graph_cond)){
    res_1 <- vector("list", comps$no)
    for(i in 1:comps$no){
      if(is_empty(seps[[i]])){
        res_1[[i]] <- qfun_MVAGGD_Gamma_full_TS(z = as.matrix(z[,v_comps[[i]]]), 
                                                maxit = maxit, start = c(start[v_comps[[i]],-c(1:2)]), nOptim = nOptim)
      }
      else{
        res_1[[i]] <- qfun_MVAGGD_Gamma_graph_TS_decomp(z = z[,v_comps[[i]]], Gamma_zero = non_edges[[i]],
                                                        cliques = cliques[[i]], seps = seps[[i]],
                                                        maxit = maxit, start = c(start[v_comps[[i]],-c(1:2)]), nOptim = nOptim)
      }
    }
    
    res <- list()
    res$par$a <- do.call(c, lapply(res_1, function(x){x$par$a}))
    res$par$b <- do.call(c, lapply(res_1, function(x){x$par$b}))
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:comps$no){
      if(length(v_comps[[i]]) == 1){
        clique_edges <- matrix(rep(v_comps[[i]], 2), ncol = 2)
      }
      else{
        clique_edges <- permutations(n = length(v_comps[[i]]), r = 2, v = v_comps[[i]], repeats.allowed = TRUE)
      }
      res$par$Gamma[clique_edges] <- res_1[[i]]$par$Gamma
    }
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
  }
  else{
    res <- qfun_MVAGGD_Gamma_graph_TS_decomp(z = z, Gamma_zero = non_edges, cliques = cliques, seps = seps,
                                             maxit = maxit, start = c(start[,-c(1:2)]), nOptim = nOptim)
  }
  
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(a_hat, b_hat, res$par$mu, res$par$sigma_1, res$par$sigma_2, res$par$shape)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out$loglike <- -sum(res$value)
  
  out$convergence <- max(res$convergence)
  
  out_Z <- as.data.frame(z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
} 

qfun_MVAGGD_Gamma_indep_TS <- function(z, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, negative = FALSE){
    #get the starting parameters
    sigma_1 <- param[2]
    sigma_2 <- param[3]
    shape <- param[4]
    
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      res <- sum(daggd(x = z, loc = param[1], scale_1 = sigma_1, scale_2 = sigma_2, shape = shape, log = TRUE))
      
      #check the value is valid
      if(is.infinite(res)){
        res <- return((-10^10)*(-1)^negative)
        warning("Infinite value of Q in mexDependence")
      }
      return(res*(-1)^negative)
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, negative = TRUE, method = "Nelder-Mead", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, negative = TRUE, method = "Nelder-Mead", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of AGG
    out <- list()
    out$par <- list(mu = fit$par[1], sigma_1 = fit$par[2], sigma_2 = fit$par[3], shape = fit$par[4])
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

qfun_MVAGGD_Gamma_graph_TS <- function(z, Gamma_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, Gamma_zero, negative = FALSE){
    #get the starting parameters
    d <- ncol(z)
    n <- nrow(z)
    mu <- param[1:d]
    sigma_1 <- param[(d + 1):(2*d)]
    sigma_2 <- param[(2*d + 1):(3*d)]
    shape <- param[(3*d + 1):(4*d)]
    
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i]))}))
      Sigma_start <- cor(Q_F_z)
      if(any(is.infinite(Q_F_z)) | any(is.na(Q_F_z)) | any(is.nan(Q_F_z)) |
         any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
        return((-10^10)*(-1)^negative)
      }
      else if(any(diag(Sigma_start) <= 0) | any(eigen(Sigma_start)$values <= 0)){
        return((-10^10)*(-1)^negative)
      }
      else{
        l_f_z <- sapply(1:d, function(i){daggd(x = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i], log = TRUE)})
        l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
        
        Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = FALSE, thr = 1e-9, zero = Gamma_zero))
        chol_Gamma <- chol(Lasso_est$wi)
        y <- chol_Gamma%*%t(Q_F_z)
        l_mvnorm <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
        
        res <- l_mvnorm + sum(l_f_z) - sum(l_dnorm)
        
        #check the value is valid
        #check the value is valid
        if(is.infinite(res)){
          res <- return((-10^10)*(-1)^negative)
          warning("Infinite value of Q in mexDependence")
        }
        return(res*(-1)^negative)
      }
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, Gamma_zero = Gamma_zero, negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, Gamma_zero = Gamma_zero, negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    sigma_1_hat <- fit$par[(d + 1):(2*d)]
    sigma_2_hat <- fit$par[(2*d + 1):(3*d)]
    shape_hat <- fit$par[(3*d + 1):(4*d)]
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, thr = 1e-9, zero = Gamma_zero))
    Gamma_hat <- Lasso_est$wi
    
    out <- list()
    out$par <- list(mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

qfun_MVAGGD_Gamma_graph_TS_decomp <- function(z, Gamma_zero, cliques, seps, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, Gamma_zero, cliques, seps, negative = FALSE){
    #get the starting parameters
    d <- ncol(z)
    n <- nrow(z)
    mu <- param[1:d]
    sigma_1 <- param[(d + 1):(2*d)]
    sigma_2 <- param[(2*d + 1):(3*d)]
    shape <- param[(3*d + 1):(4*d)]
    
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i]))}))
      Sigma_start <- cor(Q_F_z)
      if(any(is.infinite(Q_F_z)) | any(is.na(Q_F_z)) | any(is.nan(Q_F_z)) |
         any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
        return((-10^10)*(-1)^negative)
      }
      else if(any(diag(Sigma_start) <= 0) | any(eigen(Sigma_start)$values <= 0)){
        return((-10^10)*(-1)^negative)
      }
      else{
        l_f_z <- sapply(1:d, function(i){daggd(x = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i], log = TRUE)})
        l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
        l_diff <- l_f_z - l_dnorm
        
        Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, zero = Gamma_zero, thr = 1e-9))
        chol_Gamma <- chol(Lasso_est$wi)
        chol_Gamma_cliques <- lapply(cliques, function(i){as.matrix(chol_Gamma[i,i])})
        chol_Gamma_seps <- lapply(seps, function(i){as.matrix(chol_Gamma[i,i])})
        
        n_cliques <- length(cliques)
        n_seps <- length(seps)
        y_cliques <- lapply(1:n_cliques, function(i){chol_Gamma_cliques[[i]]%*%t(Q_F_z[,cliques[[i]]])})
        y_seps <- lapply(1:n_seps, function(i){chol_Gamma_seps[[i]]%*%t(Q_F_z[,seps[[i]]])})
        
        l_mvnorm_cliques <- sapply(1:n_cliques, function(i){
          -n*length(cliques[[i]])*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma_cliques[[i]]))) - sum(y_cliques[[i]]^2)/2
        })
        
        l_mvnorm_seps <- sapply(1:n_seps, function(i){
          -n*length(seps[[i]])*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma_seps[[i]]))) - sum(y_seps[[i]]^2)/2
        })
        
        res <- sum(l_mvnorm_cliques) - sum(l_mvnorm_seps) + sum(l_diff)
        
        #check the value is valid
        if(is.infinite(res)){
          res <- return((-10^10)*(-1)^negative)
          warning("Infinite value of Q in mexDependence")
        }
        return(res*(-1)^negative)
      }
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, Gamma_zero = Gamma_zero, cliques = cliques, seps = seps,
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, Gamma_zero = Gamma_zero, cliques = cliques, seps = seps,
                       negative = TRUE, method = "Nelder-Mead", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    sigma_1_hat <- fit$par[(d + 1):(2*d)]
    sigma_2_hat <- fit$par[(2*d + 1):(3*d)]
    shape_hat <- fit$par[(3*d + 1):(4*d)]
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, thr = 1e-9, zero = Gamma_zero))
    Gamma_hat <- Lasso_est$wi
    
    out <- list()
    out$par <- list(mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

qfun_MVAGGD_Gamma_full_TS <- function(z, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, negative = FALSE){
    #get the starting parameters
    d <- ncol(z)
    n <- nrow(z)
    mu <- param[1:d]
    sigma_1 <- param[(d + 1):(2*d)]
    sigma_2 <- param[(2*d + 1):(3*d)]
    shape <- param[(3*d + 1):(4*d)]
    
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i]))}))
      Sigma_start <- cor(Q_F_z)
      if(any(is.infinite(Q_F_z)) | any(is.na(Q_F_z)) | any(is.nan(Q_F_z)) |
         any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
        return((-10^10)*(-1)^negative)
      }
      else if(any(diag(Sigma_start) <= 0) | any(eigen(Sigma_start)$values <= 0)){
        return((-10^10)*(-1)^negative)
      }
      else{
        l_f_z <- sapply(1:d, function(i){daggd(x = z[,i], loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i], log = TRUE)})
        l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
        
        Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = FALSE, thr = 1e-9))
        chol_Gamma <- chol(Lasso_est$wi)
        y <- chol_Gamma%*%t(Q_F_z)
        l_mvnorm <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
        
        res <- l_mvnorm + sum(l_f_z) - sum(l_dnorm)
        
        #check the value is valid
        if(is.infinite(res)){
          res <- return((-10^10)*(-1)^negative)
          warning("Infinite value of Q in mexDependence")
        }
        return(res*(-1)^negative)
      }
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, negative = TRUE, method = "Nelder-Mead", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    sigma_1_hat <- fit$par[(d + 1):(2*d)]
    sigma_2_hat <- fit$par[(2*d + 1):(3*d)]
    shape_hat <- fit$par[(3*d + 1):(4*d)]
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(paggd(q = z[,i], loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, thr = 1e-9))
    Gamma_hat <- Lasso_est$wi
    
    out <- list()
    out$par <- list(mu = mu_hat, sigma_1 = sigma_1_hat, sigma_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

################################################################################
## Improved functions to fit the one-step method when we assume that the residuals
## are multivariate Gaussian

Cond_extremes_graph_new <- function(data, cond, graph = NA, 
                                    constrain = TRUE, q = c(0,1), v = 10, aLow = -1, 
                                    maxit = 1e+6, start = c(0.1, 0.1), nOptim = 1){
  
  theCall <- match.call()
  # if (!inherits(x, "migpd")) 
  #   stop("you need to use an object created by migpd")
  # margins <- list(casefold(margins), p2q = switch(casefold(margins), 
  #                                                 gumbel = function(p) -log(-log(p)), laplace = function(p) ifelse(p < 
  #                                                                                                                    0.5, log(2 * p), -log(2 * (1 - p)))), q2p = switch(casefold(margins), 
  #                                                                                                                                                                       gumbel = function(q) exp(-exp(-q)), laplace = function(q) ifelse(q < 
  #                                                                                                                                                                                                                                          0, exp(q)/2, 1 - 0.5 * exp(-q))))
  # x <- mexTransform(x, margins = margins, method = marTransform, 
  #                   r = referenceMargin)
  # x$referenceMargin <- referenceMargin
  # if (margins[[1]] == "gumbel" & constrain) {
  #   warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
  #   constrain <- FALSE
  # }
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(is.character(cond)){ 
    cond <- match(cond, dimnames(x$transformed)[[2]])
  }
  # if(missing(dqu)){
  #   message("Assuming same quantile for dependence thesholding as was used \n     to fit corresponding marginal model...\n")
  #   dqu <- x$mqu[cond]
  # }
  # dth <- quantile(x$transformed[, cond], dqu)
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  dependent <- (1:(dim(data)[[2]]))[-cond]
  # if(length(dqu) < length(dependent)){
  #   dqu <- rep(dqu, length = length(dependent))
  #   aLow <- if else(margins[[1]] == "gumbel", 10^(-10), -1 + 10^(-10))
  # }
  if(missing(start)){
    start <- c(0.1, 0.1)
  }
  else if(inherits(start, "mex")){
    start <- start$dependence$coefficients[1:2,]
  }
  if(length(start) == 2){
    start <- matrix(rep(start, d-1), ncol = 2, byrow = TRUE)
  }
  
  ## determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    non_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    
    ## Now figure out the cliques and separators
    graph_cond <- delete_vertices(graph, cond)
    comps <- components(graph_cond)
    
    if(comps$no != 1){
      ## Get the separate components
      v_comps <- groups(comps)
      g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
      
      ## Get the non_edges in the graph(s)
      non_edges <- vector("list", comps$no)
      for(i in 1:comps$no){
        if(length(V(g_comps[[i]])) <= 2){
          non_edges[[i]] <- matrix(NA, nrow = 0, ncol = 2)
        }
        else{
          all_edges <- data.frame(combinations(n = length(V(g_comps[[i]])), r = 2, v = V(g_comps[[i]])))
          edges_in_comp <- data.frame(as_edgelist(g_comps[[i]]))
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_comp)
          non_edges[[i]] <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        }
      }
      
      ## Get the cliques and separators
      for(i in 1:comps$no){
        V(g_comps[[i]])$name <- v_comps[[i]]
      }
      cliques <- lapply(g_comps, function(x){max_cliques(x)})
      seps <- lapply(cliques, function(x){get_separators(x, includeSingletons = inc_sings)})
      
    }
    else{
      #get the cliques and the separators in the graph
      cliques <- max_cliques(graph_cond)
      #now get the separators but we only include singletons a vertex only appears once in a clique of length two
      seps <- get_separators(cliques, includeSingletons = inc_sings)
      
      ## Get the non_edges in the graph
      all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
      edges_in_graph <- data.frame(as_edgelist(graph_cond))
      all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
      non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
    } 
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- data[,cond]
  ydep <- data[,-cond]
  
  #fit the log-likelihood
  n_cliques <- length(cliques)
  res <- list()
  if(is_empty(seps) & n_cliques == d-1){
    res_1 <- lapply(cliques, function(i){
      qfun_HT(yex = yex, ydep = as.matrix(ydep[,i]), 
              constrain = constrain, q = q, v = v, aLow = aLow,
              maxit = maxit, start = start[i,], nOptim = nOptim)}) 
    
    #get the output
    res$par$a <- sapply(res_1, function(x){x$par$a})
    res$par$b <- sapply(res_1, function(x){x$par$b})
    res$par$mu <- sapply(res_1, function(x){x$par$mu})
    res$par$Gamma <- diag(1/sapply(res_1, function(x){x$par$Sigma}))
    
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
    res$Z <- sapply(res_1, function(x){x$Z})
  }
  else if(is_empty(seps) & length(cliques) == 1){
    res <- qfun_MVN_Gamma_full(yex = yex, ydep = as.matrix(ydep),
                               maxit = maxit, start = c(start), nOptim = nOptim)
  }
  else if(!is_connected(graph_cond)){
    res_1 <- vector("list", comps$no)
    for(i in 1:comps$no){
      if(is_empty(seps[[i]])){
        res_1[[i]] <- qfun_MVN_Gamma_full(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]),
                                          maxit = maxit, 
                                          start = c(start[v_comps[[i]],]), 
                                          nOptim = nOptim)
      }
      else{
        res_1[[i]] <- qfun_MVN_Gamma_graph(yex = yex, ydep = ydep[,v_comps[[i]]], 
                                           Gamma_zero = non_edges[[i]],
                                           maxit = maxit, 
                                           start = c(start[v_comps[[i]],]), 
                                           nOptim = nOptim)
      }
    }
    
    res$par$a <- do.call(c, lapply(res_1, function(x){x$par$a}))
    res$par$b <- do.call(c, lapply(res_1, function(x){x$par$b}))
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$Gamma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:comps$no){
      if(length(v_comps[[i]]) == 1){
        clique_edges <- matrix(rep(v_comps[[i]], 2), ncol = 2)
      }
      else{
        clique_edges <- permutations(n = length(v_comps[[i]]), r = 2, v = v_comps[[i]], repeats.allowed = TRUE)
      }
      res$par$Gamma[clique_edges] <- res_1[[i]]$par$Gamma
    }
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
    res$Z <- do.call(cbind, lapply(res_1, function(x){x$Z}))
  }
  else{
    res <- qfun_MVN_Gamma_graph(yex = yex, ydep = ydep, Gamma_zero = non_edges,
                                maxit = maxit, start = c(start), nOptim = nOptim) 
  }
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(res$par$a, res$par$b, res$par$mu)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "mu")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out$loglike <- -sum(res$value)
  
  out$convergence <- max(res$convergence)
  
  out_Z <- as.data.frame(res$Z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
}

qfun_MVN_Gamma_graph <- function(yex, ydep, Gamma_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, Gamma_zero, negative = FALSE){
    ## Obtain parameters
    d <- ncol(ydep)
    n <- length(yex)
    a <- param[1:d]
    b <- param[(d + 1):(2*d)]
    
    #conditions on the parameters
    if(any(abs(a) >= 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## Obtain the residuals
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      Sigma_start <- cov(z)
      if(any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
        return((-10^10)*(-1)^negative)
      }
      else if(any(diag(Sigma_start) <= 0) | any(eigen(Sigma_start)$values <= 0)){
        return((-10^10)*(-1)^negative)
      }
      else{
        Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = TRUE, thr = 1e-9, zero = Gamma_zero))
        chol_Gamma <- chol(Lasso_est$wi)
        mu <- matrix(rep(apply(z, 2, mean), n), ncol = d, byrow = TRUE)
        y <- chol_Gamma%*%t(z - mu)
        res <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2 - sum(log(b_yi))
        
        if(is.infinite(res)){
          warning("Infinite value of Q in mexDependence")
          return((-10^10)*(-1)^negative)
        }
      }
    }
    return(res*(-1)^negative)
  }
  
  ## Fit the model
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep, Gamma_zero = Gamma_zero,
                   negative = TRUE,  control = list(maxit = maxit), method = "BFGS"),
             silent = TRUE)
  
  d <- ncol(ydep)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = d),
                    b = rep(NA, times = d),
                    mu = rep(NA, times = d),
                    Sigma = diag(NA, d))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(fit$convergence != 0){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = d),
                    b = rep(NA, times = d),
                    mu = rep(NA, times = d),
                    Sigma = diag(NA, d))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep, Gamma_zero = Gamma_zero,
                       negative = TRUE,  control = list(maxit = maxit), method = "BFGS"),
                 silent = TRUE)
      
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = d),
                        b = rep(NA, times = d),
                        mu = rep(NA, times = d),
                        Sigma = diag(NA, d))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
      }
      else if(fit$convergence != 0){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = d),
                        b = rep(NA, times = d),
                        mu = rep(NA, times = d),
                        Sigma = diag(NA, d))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
      }
    }
  }
  
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[(d + 1):(2*d)]
    
    #obtain the residuals
    a_yi_hat <- sapply(a_hat, function(x){yex*x})
    b_yi_hat <- sapply(b_hat, function(x){yex^x})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    mu_hat <- apply(Z, 2, mean)
    Lasso_est <- suppressWarnings(glasso(s = cov(Z), rho = 0, penalize.diagonal = TRUE, thr = 1e-9, zero = Gamma_zero))
    
    ## Organise the output
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = mu_hat, Gamma = Lasso_est$wi)
    out$value <- fit$value
    out$convergence <- fit$convergence
    out$Z <- Z
  }
  return(out)
}

qfun_MVN_Gamma_full <- function(yex, ydep, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, negative = FALSE){
    ## Obtain parameters
    d <- ncol(ydep)
    n <- length(yex)
    a <- param[1:d]
    b <- param[(d + 1):(2*d)]
    
    #conditions on the parameters
    if(any(abs(a) >= 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## Obtain the residuals
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      ## Now maximise the Z
      mu_z <- matrix(rep(apply(z, 2, mean), n), ncol = d, byrow = TRUE)
      chol_Gamma <- chol(solve(cov(z)))
      
      y <- chol_Gamma%*%t(z - mu_z)
      llh_z <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
      
      res <- (llh_z - sum(log(b_yi)))
      if(is.infinite(res)){
        warning("Infinite value of Q in mexDependence")
        return((-10^10)*(-1)^negative)
      }
    }
    return(res*(-1)^negative)
  }
  
  ## Fit the model
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                   negative = TRUE,  control = list(maxit = maxit), method = "BFGS"),
             silent = TRUE)
  
  d <- ncol(ydep)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = d),
                    b = rep(NA, times = d),
                    mu = rep(NA, times = d),
                    Sigma = diag(NA, d))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(fit$convergence != 0){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, times = d),
                    b = rep(NA, times = d),
                    mu = rep(NA, times = d),
                    Sigma = diag(NA, d))
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                       negative = TRUE,  control = list(maxit = maxit), method = "BFGS"),
                 silent = TRUE)
      
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = d),
                        b = rep(NA, times = d),
                        mu = rep(NA, times = d),
                        Sigma = diag(NA, d))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
      }
      else if(fit$convergence != 0){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = rep(NA, times = d),
                        b = rep(NA, times = d),
                        mu = rep(NA, times = d),
                        Sigma = diag(NA, d))
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
      }
    }
  }
  
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[-(1:d)]
    
    #obtain the residuals
    a_yi_hat <- sapply(a_hat, function(x){yex*x})
    b_yi_hat <- sapply(b_hat, function(x){yex^x})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## Organise the output
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = apply(Z, 2, mean), Gamma = solve(cov(Z)))
    out$value <- fit$value
    out$convergence <- fit$convergence
    out$Z <- Z
  }
  return(out)
}

################################################################################
## Two-step MVN model
## Fits the original H+T model
## Then fits a MVN to the fitted residuals
## Purely for comparative purposes and there is no expectation this model will perform well
## as the marginal distributions of the fitted residuals are AGG

Cond_extremes_MVN_TS <- function(data, cond, graph = NA, start = c(0.1, 0.1),
                                 constrain = TRUE, q = c(0,1), v = 10, aLow = -1,
                                 maxit = 1e+6, nOptim = 1){
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  if(missing(start)){
    start <- c(0.1, 0.1)
  }
  else if(inherits(start, "mex")){
    start <- start$dependence$coefficients[1:2,]
  }
  if(length(start) == 2){
    start <- matrix(rep(start, d-1), ncol = 2, byrow = TRUE)
  }
  
  dependent <- (1:(dim(data)[[2]]))[-cond]
  
  ## First step where we fit the original conditional extremes model
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  res_OG <- lapply(1:(d-1), function(i){
    qfun_HT(yex = yex, ydep = as.matrix(ydep[,i]),
            constrain = constrain, aLow = aLow, q = q, v = v,
            maxit = maxit, start = start[i,], nOptim = nOptim)})
  
  ## Get the necessary output
  z <- sapply(res_OG, function(x){x$Z})
  a_hat <- sapply(res_OG, function(x){x$par$a})
  b_hat <- sapply(res_OG, function(x){x$par$b})
  
  ## Now fit the graphical model to Z|i
  ## determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    edges <- matrix(rep(1:(d-1), 2), ncol = 2)
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    
    ## Now figure out the cliques and separators
    graph_cond <- delete_vertices(graph, cond)
    comps <- components(graph_cond)
    
    if(comps$no != 1){
      ## Get the separate components
      v_comps <- groups(comps)
      g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
      
      ## Get the non_edges in the graph(s)
      edges <- vector("list", comps$no)
      for(i in 1:comps$no){
        if(length(V(g_comps[[i]])) <= 1){
          edges[[i]] <- matrix(rep(V(g_comps[[i]]), 2), ncol = 2)
        }
        else{
          all_edges <- data.frame(combinations(n = length(V(g_comps[[i]])), r = 2, v = V(g_comps[[i]]), repeats.allowed = TRUE))
          edges_in_comp <- data.frame(rbind(as_edgelist(g_comps[[i]]), matrix(rep(V(g_comps[[i]]), 2), ncol = 2)))
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_comp)
          edges[[i]] <- as.matrix(all_edges[which(all_edges$exists == TRUE), 1:2])
        }
      }
      
      ## Get the cliques and separators
      for(i in 1:comps$no){
        V(g_comps[[i]])$name <- v_comps[[i]]
      }
      cliques <- lapply(g_comps, function(x){max_cliques(x)})
      seps <- lapply(cliques, function(x){get_separators(x, includeSingletons = inc_sings)})
      
    }
    else{
      #get the cliques and the separators in the graph
      cliques <- max_cliques(graph_cond)
      #now get the separators but we only include singletons a vertex only appears once in a clique of length two
      seps <- get_separators(cliques, includeSingletons = inc_sings)
      
      ## Get the non_edges in the graph
      all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1), repeats.allowed = TRUE))
      edges_in_graph <- data.frame(rbind(as_edgelist(graph_cond), matrix(rep(1:(d-1), 2), ncol = 2)))
      all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
      edges <- as.matrix(all_edges[which(all_edges$exists == TRUE), 1:2])
    } 
  }
  
  ## fit the models and get the output
  mu_start <- apply(z, 2, mean)
  Gamma_start <- solve(cov(z))
  Gamma_chol <- chol(Gamma_start)
  
  res <- list()
  n_cliques <- length(cliques)
  if(is_empty(seps) & n_cliques == d-1){
    ## Make this function and do them independently of one another
    res$par$mu <- mu_start
    res$par$Gamma <- diag(1/apply(z, 2, var), d-1)
    res$value <- -sum(dmvnorm(x = z, mean = mu_start, sigma = solve(res$par$Gamma), log = TRUE))
    res$convergence <- 0
  }
  else if(is_empty(seps) & n_cliques == 1){
    res$par$mu <- mu_start
    res$par$Gamma <- Gamma_start
    res$value <- -sum(dmvnorm(x = z, mean = mu_start, sigma = solve(res$par$Gamma), log = TRUE))
    res$convergence <- 0
  }
  else if(!is_connected(graph_cond)){
    res_1 <- vector("list", comps$no)
    for(i in 1:comps$no){
      if(is_empty(seps[[i]])){
        res_1[[i]] <- qfun_MVN_full_TS(z = as.matrix(z[,v_comps[[i]]]),
                                       Gamma_non_zero = edges[[i]],
                                       maxit = maxit, 
                                       start = c(mu_start[v_comps[[i]]], Gamma_chol[edges[[i]]]), 
                                       nOptim = nOptim)
      }
      else{
        res_1[[i]] <- qfun_MVN_graph_TS(z = z[,v_comps[[i]]], 
                                        Gamma_non_zero = edges[[i]],
                                        maxit = maxit, 
                                        start = c(mu_start[v_comps[[i]]], Gamma_chol[edges[[i]]]), 
                                        nOptim = nOptim)
      }
    }
    
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$Gamma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:comps$no){
      if(length(v_comps[[i]]) == 1){
        clique_edges <- matrix(rep(v_comps[[i]], 2), ncol = 2)
      }
      else{
        clique_edges <- permutations(n = length(v_comps[[i]]), r = 2, v = v_comps[[i]], repeats.allowed = TRUE)
      }
      res$par$Gamma[clique_edges] <- res_1[[i]]$par$Gamma
    }
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
  }
  else{
    res <- qfun_MVN_graph_TS(z = z, Gamma_non_zero = edges,
                             maxit = maxit, start = c(mu_start, Gamma_chol[edges]), 
                             nOptim = nOptim)
  }
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(a_hat, b_hat, res$par$mu)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out$loglike <- -sum(res$value)
  
  out$convergence <- max(res$convergence)
  
  out_Z <- as.data.frame(z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
} 

qfun_MVN_graph_TS <- function(z, Gamma_non_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, Gamma_non_zero, negative = FALSE){
    #get the starting parameters
    d <- ncol(z)
    n <- nrow(z)
    mu <- param[1:d]
    Gamma_chol <- matrix(0, nrow = d, ncol = d)
    Gamma_chol[Gamma_non_zero] <- param[-c(1:d)]
    
    if(any(diag(Gamma_chol) <= 0) | any(eigen(Gamma_chol)$values <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      mu_mat <- matrix(rep(mu, n), nrow = n, ncol = d, byrow = TRUE)
      y <- Gamma_chol%*%t(z - mu_mat)
      res <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(Gamma_chol))) - sum(y^2)/2
      
      #check the value is valid
      if(is.infinite(res)){
        if(res < 0) {
          res <- return((-10^10)*(-1)^negative)
        }
        else{ 
          res <- return((-10^10)*(-1)^negative)
        }
        warning("Infinite value of Q in mexDependence")
      }
      return(res*(-1)^negative)
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, Gamma_non_zero = Gamma_non_zero, negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, Gamma_non_zero = Gamma_non_zero, negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    Gamma_chol <- matrix(0, nrow = d, ncol = d)
    Gamma_chol[Gamma_non_zero] <- fit$par[-c(1:d)]
    Gamma_hat <- t(Gamma_chol)%*% Gamma_chol
    
    out <- list()
    out$par <- list(mu = mu_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

qfun_MVN_full_TS <- function(z, Gamma_non_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, Gamma_non_zero, negative = FALSE){
    #get the starting parameters
    d <- ncol(z)
    n <- nrow(z)
    mu <- param[1:d]
    Gamma_chol <- matrix(0, nrow = d, ncol = d)
    Gamma_chol[Gamma_non_zero] <- param[-c(1:d)]
    
    if(any(diag(Gamma_chol) <= 0) | any(eigen(Gamma_chol)$values <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      mu_mat <- matrix(rep(mu, n), nrow = n, ncol = d, byrow = TRUE)
      y <- Gamma_chol%*%t(z - mu_mat)
      res <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(Gamma_chol))) - sum(y^2)/2
      
      #check the value is valid
      if(is.infinite(res)){
        if(res < 0) {
          res <- return((-10^10)*(-1)^negative)
        }
        else{ 
          res <- return((-10^10)*(-1)^negative)
        }
        warning("Infinite value of Q in mexDependence")
      }
      return(res*(-1)^negative)
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, Gamma_non_zero = Gamma_non_zero, 
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, Gamma_non_zero = Gamma_non_zero,
                       negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    Gamma_chol <- matrix(0, nrow = d, ncol = d)
    Gamma_chol[Gamma_non_zero] <- fit$par[-c(1:d)]
    Gamma_hat <- t(Gamma_chol)%*% Gamma_chol
    
    out <- list()
    out$par <- list(mu = mu_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

Cond_extremes_MVN_TS_Z <- function(z, a, b, cond, graph = NA, maxit = 1e+6, nOptim = 1){
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  
  #get information from the data
  dim_data <- dim(z)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 1 columns")
  }
  else{
    d <- dim_data[2] + 1
    n <- dim_data[1]
  }
  
  dependent <- (1:d)[-cond]
  
  ## Now fit the graphical model to Z|i
  ## determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    edges <- matrix(rep(1:(d-1), 2), ncol = 2)
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    
    ## Now figure out the cliques and separators
    graph_cond <- delete_vertices(graph, cond)
    comps <- components(graph_cond)
    
    if(comps$no != 1){
      ## Get the separate components
      v_comps <- groups(comps)
      g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
      
      ## Get the non_edges in the graph(s)
      edges <- vector("list", comps$no)
      for(i in 1:comps$no){
        if(length(V(g_comps[[i]])) <= 1){
          edges[[i]] <- matrix(rep(V(g_comps[[i]]), 2), ncol = 2)
        }
        else{
          all_edges <- data.frame(combinations(n = length(V(g_comps[[i]])), r = 2, v = V(g_comps[[i]]), repeats.allowed = TRUE))
          edges_in_comp <- data.frame(rbind(as_edgelist(g_comps[[i]]), matrix(rep(V(g_comps[[i]]), 2), ncol = 2)))
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_comp)
          edges[[i]] <- as.matrix(all_edges[which(all_edges$exists == TRUE), 1:2])
        }
      }
      
      ## Get the cliques and separators
      for(i in 1:comps$no){
        V(g_comps[[i]])$name <- v_comps[[i]]
      }
      cliques <- lapply(g_comps, function(x){max_cliques(x)})
      seps <- lapply(cliques, function(x){get_separators(x, includeSingletons = inc_sings)})
      
    }
    else{
      #get the cliques and the separators in the graph
      cliques <- max_cliques(graph_cond)
      #now get the separators but we only include singletons a vertex only appears once in a clique of length two
      seps <- get_separators(cliques, includeSingletons = inc_sings)
      
      ## Get the non_edges in the graph
      all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1), repeats.allowed = TRUE))
      edges_in_graph <- data.frame(rbind(as_edgelist(graph_cond), matrix(rep(1:(d-1), 2), ncol = 2)))
      all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
      edges <- as.matrix(all_edges[which(all_edges$exists == TRUE), 1:2])
    } 
  }
  
  ## fit the models and get the output
  mu_start <- apply(z, 2, mean)
  Gamma_start <- solve(cov(z))
  Gamma_chol <- chol(Gamma_start)
  
  res <- list()
  n_cliques <- length(cliques)
  if(is_empty(seps) & n_cliques == d-1){
    ## Make this function and do them independently of one another
    res$par$mu <- mu_start
    res$par$Gamma <- diag(1/apply(z, 2, var), d-1)
    res$value <- -sum(dmvnorm(x = z, mean = mu_start, sigma = solve(res$par$Gamma), log = TRUE))
    res$convergence <- 0
  }
  else if(is_empty(seps) & n_cliques == 1){
    res$par$mu <- mu_start
    res$par$Gamma <- Gamma_start
    res$value <- -sum(dmvnorm(x = z, mean = mu_start, sigma = solve(res$par$Gamma), log = TRUE))
    res$convergence <- 0
  }
  else if(!is_connected(graph_cond)){
    res_1 <- vector("list", comps$no)
    for(i in 1:comps$no){
      if(is_empty(seps[[i]])){
        res_1[[i]] <- qfun_MVN_full_TS(z = as.matrix(z[,v_comps[[i]]]),
                                       Gamma_non_zero = edges[[i]],
                                       maxit = maxit, 
                                       start = c(mu_start[v_comps[[i]]], Gamma_chol[edges[[i]]]), 
                                       nOptim = nOptim)
      }
      else{
        res_1[[i]] <- qfun_MVN_graph_TS(z = z[,v_comps[[i]]], 
                                        Gamma_non_zero = edges[[i]],
                                        maxit = maxit, 
                                        start = c(mu_start[v_comps[[i]]], Gamma_chol[edges[[i]]]), 
                                        nOptim = nOptim)
      }
    }
    
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$Gamma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:comps$no){
      if(length(v_comps[[i]]) == 1){
        clique_edges <- matrix(rep(v_comps[[i]], 2), ncol = 2)
      }
      else{
        clique_edges <- permutations(n = length(v_comps[[i]]), r = 2, v = v_comps[[i]], repeats.allowed = TRUE)
      }
      res$par$Gamma[clique_edges] <- res_1[[i]]$par$Gamma
    }
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
  }
  else{
    res <- qfun_MVN_graph_TS(z = z, Gamma_non_zero = edges,
                             maxit = maxit, start = c(mu_start, Gamma_chol[edges]), 
                             nOptim = nOptim)
  }
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(a, b, res$par$mu)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out$loglike <- -sum(res$value)
  
  out$convergence <- max(res$convergence)
  
  out_Z <- as.data.frame(z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
} 

qfun_MVN_graph_TS <- function(z, Gamma_non_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, Gamma_non_zero, negative = FALSE){
    #get the starting parameters
    d <- ncol(z)
    n <- nrow(z)
    mu <- param[1:d]
    Gamma_chol <- matrix(0, nrow = d, ncol = d)
    Gamma_chol[Gamma_non_zero] <- param[-c(1:d)]
    
    if(any(diag(Gamma_chol) <= 0) | any(eigen(Gamma_chol)$values <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      mu_mat <- matrix(rep(mu, n), nrow = n, ncol = d, byrow = TRUE)
      y <- Gamma_chol%*%t(z - mu_mat)
      res <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(Gamma_chol))) - sum(y^2)/2
      
      #check the value is valid
      if(is.infinite(res)){
        if(res < 0) {
          res <- return((-10^10)*(-1)^negative)
        }
        else{ 
          res <- return((-10^10)*(-1)^negative)
        }
        warning("Infinite value of Q in mexDependence")
      }
      return(res*(-1)^negative)
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, Gamma_non_zero = Gamma_non_zero, negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, Gamma_non_zero = Gamma_non_zero, negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    Gamma_chol <- matrix(0, nrow = d, ncol = d)
    Gamma_chol[Gamma_non_zero] <- fit$par[-c(1:d)]
    Gamma_hat <- t(Gamma_chol)%*% Gamma_chol
    
    out <- list()
    out$par <- list(mu = mu_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

qfun_MVN_full_TS <- function(z, Gamma_non_zero, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, z, Gamma_non_zero, negative = FALSE){
    #get the starting parameters
    d <- ncol(z)
    n <- nrow(z)
    mu <- param[1:d]
    Gamma_chol <- matrix(0, nrow = d, ncol = d)
    Gamma_chol[Gamma_non_zero] <- param[-c(1:d)]
    
    if(any(diag(Gamma_chol) <= 0) | any(eigen(Gamma_chol)$values <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      mu_mat <- matrix(rep(mu, n), nrow = n, ncol = d, byrow = TRUE)
      y <- Gamma_chol%*%t(z - mu_mat)
      res <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(Gamma_chol))) - sum(y^2)/2
      
      #check the value is valid
      if(is.infinite(res)){
        if(res < 0) {
          res <- return((-10^10)*(-1)^negative)
        }
        else{ 
          res <- return((-10^10)*(-1)^negative)
        }
        warning("Infinite value of Q in mexDependence")
      }
      return(res*(-1)^negative)
    }
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, Gamma_non_zero = Gamma_non_zero, 
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(mu = rep(NA, times = ncol(z)),
                    sigma_1 = rep(NA, times = ncol(z)),
                    sigma_2 = rep(NA, times = ncol(z)),
                    shape = diag(NA, ncol(z)))
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       z = z, Gamma_non_zero = Gamma_non_zero,
                       negative = TRUE, method = "BFGS", hessian = FALSE),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(mu = rep(NA, times = ncol(z)),
                        sigma_1 = rep(NA, times = ncol(z)),
                        sigma_2 = rep(NA, times = ncol(z)),
                        shape = diag(NA, ncol(z)))
        out$value <- NA
        out$convergence <- NA
        break
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    Gamma_chol <- matrix(0, nrow = d, ncol = d)
    Gamma_chol[Gamma_non_zero] <- fit$par[-c(1:d)]
    Gamma_hat <- t(Gamma_chol)%*% Gamma_chol
    
    out <- list()
    out$par <- list(mu = mu_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  return(out)
}

Cond_extremes_MVAGGD_Gamma_TS_new_Z <- function(z, a, b, cond, graph = NA, start,
                                                maxit = 1e+6, nOptim = 1){
  
  #need to update this section of code
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  
  #get information from the data
  dim_data <- dim(z)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 1 columns")
  }
  else{
    d <- dim_data[2] + 1
    n <- dim_data[1]
  }
  
  dependent <- (1:d)[-cond]
  
  ### Fit the graphical model to Z|i
  ## determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## get the cliques separators and non-edges in the graph
    cliques <- sapply(1:(d - length(cond)), list)
    seps <- list()
    non_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    
    #determine if removing a single nodes results in more than one component
    #If it does it is a separator and we need to include single separators
    graphs_ind <- lapply(1:d, delete_vertices, graph = graph)
    comps_ind <- lapply(graphs_ind, components)
    comps_ind_unique <- sapply(comps_ind, function(x){length(unique(x$membership))})
    if(any(comps_ind_unique != 1)){
      inc_sings <- TRUE
    }
    else{
      inc_sings <- FALSE
    }
    
    ## Now figure out the cliques and separators
    graph_cond <- delete_vertices(graph, cond)
    comps <- components(graph_cond)
    
    if(comps$no != 1){
      ## Get the separate components
      v_comps <- groups(comps)
      g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
      
      ## Get the non_edges in the graph(s)
      non_edges <- vector("list", comps$no)
      for(i in 1:comps$no){
        if(length(V(g_comps[[i]])) <= 2){
          non_edges[[i]] <- matrix(NA, nrow = 0, ncol = 2)
        }
        else{
          all_edges <- data.frame(combinations(n = length(V(g_comps[[i]])), r = 2, v = V(g_comps[[i]])))
          edges_in_comp <- data.frame(as_edgelist(g_comps[[i]]))
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_comp)
          non_edges[[i]] <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        }
      }
      
      ## Get the cliques and separators
      for(i in 1:comps$no){
        V(g_comps[[i]])$name <- v_comps[[i]]
      }
      cliques <- lapply(g_comps, function(x){max_cliques(x)})
      seps <- lapply(cliques, function(x){get_separators(x, includeSingletons = inc_sings)})
      
    }
    else{
      #get the cliques and the separators in the graph
      cliques <- max_cliques(graph_cond)
      #now get the separators but we only include singletons a vertex only appears once in a clique of length two
      seps <- get_separators(cliques, includeSingletons = inc_sings)
      
      ## Get the non_edges in the graph
      all_edges <- data.frame(combinations(n = d - 1, r = 2, v = 1:(d-1)))
      edges_in_graph <- data.frame(as_edgelist(graph_cond))
      all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, edges_in_graph)
      non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
    }
  }
  
  n_cliques <- length(cliques)
  ## fit the models and get the output
  if(is_empty(seps) & n_cliques == d-1){
    ## Make this function and do them independently of one another
    res_1 <- lapply(1:(d-1), function(i){
      qfun_MVAGGD_Gamma_indep_TS(z = as.matrix(z[,i]), maxit = maxit, start = c(start[i,]), nOptim = nOptim)})
    
    res <- list()
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- diag(1, d-1)
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
  }
  else if(is_empty(seps) & n_cliques == 1){
    res <- qfun_MVAGGD_Gamma_full_TS(z = z, maxit = maxit, start = c(start), nOptim = nOptim)
  }
  else if(!is_connected(graph_cond)){
    res_1 <- vector("list", comps$no)
    for(i in 1:comps$no){
      if(is_empty(seps[[i]])){
        res_1[[i]] <- qfun_MVAGGD_Gamma_full_TS(z = as.matrix(z[,v_comps[[i]]]), 
                                                maxit = maxit, start = c(start[v_comps[[i]],]), nOptim = nOptim)
      }
      else{
        res_1[[i]] <- qfun_MVAGGD_Gamma_graph_TS(z = z[,v_comps[[i]]], Gamma_zero = non_edges[[i]],
                                                 maxit = maxit, start = c(start[v_comps[[i]],]), nOptim = nOptim)
      }
    }
    
    res <- list()
    res$par$mu <- do.call(c, lapply(res_1, function(x){x$par$mu}))
    res$par$sigma_1 <- do.call(c, lapply(res_1, function(x){x$par$sigma_1}))
    res$par$sigma_2 <- do.call(c, lapply(res_1, function(x){x$par$sigma_2}))
    res$par$shape <- do.call(c, lapply(res_1, function(x){x$par$shape}))
    res$par$Gamma <- matrix(0, nrow = d-1, ncol = d-1)
    for(i in 1:comps$no){
      if(length(v_comps[[i]]) == 1){
        clique_edges <- matrix(rep(v_comps[[i]], 2), ncol = 2)
      }
      else{
        clique_edges <- permutations(n = length(v_comps[[i]]), r = 2, v = v_comps[[i]], repeats.allowed = TRUE)
      }
      res$par$Gamma[clique_edges] <- res_1[[i]]$par$Gamma
    }
    res$value <- sapply(res_1, function(x){x$value})
    res$convergence <- sapply(res_1, function(x){x$convergence})
  }
  else{
    res <- qfun_MVAGGD_Gamma_graph_TS(z = z, Gamma_zero = non_edges,
                                      maxit = maxit, start = c(start), nOptim = nOptim)
  }
  
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(a, b, res$par$mu, res$par$sigma_1, res$par$sigma_2, res$par$shape)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  Gamma_out <- res$par$Gamma
  
  out$par <- list(main = par_out, Gamma = Gamma_out)
  
  out$loglike <- -sum(res$value)
  
  out$convergence <- max(res$convergence)
  
  out_Z <- as.data.frame(z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
} 


################################################################################

## Extrapolations of X to calculate probabilities such as P[X_{j} > u_{j} | X_{i} > u_{i}]

## Get P[X_{j} > u_{j} | X_{i} > u_{i}] for a data set
p_data_surv_multi <- function(x, cond, u, uncon){
  #get dimension of data
  d <- dim(x)[2]
  
  #get number of times Xi > ui
  nu_x <- length(which(x[,cond] > u[cond]))
  
  #get the probabilities 
  index_uni <- which(sapply(uncon, function(x){length(x) == 1}))
  if(!is_empty(index_uni)){
    combos_uni <- uncon[index_uni]
    cond_in_uni <- sapply(combos_uni, function(x){cond%in%x})
    if(!any(cond_in_uni)){
      index_cond_not_in_uni <- which(cond_in_uni == FALSE)
      for(i in 1:length(index_cond_not_in_uni)){
        combos_uni[[index_cond_not_in_uni[i]]] <- c(combos_uni[[index_cond_not_in_uni[i]]], cond)
      }
    }
    p_uni <- sapply(combos_uni, function(i){length(which(apply(x[,i], 1, function(y){all(y > u[i])})))/nu_x})  
    if(length(index_uni) == length(uncon)){
      return(p_uni)
    }
    else{
      #get the multivariate probabilities
      combos_multi <- uncon[-index_uni]
      cond_in_multi <- sapply(combos_multi, function(x){cond%in%x})
      if(!any(cond_in_multi)){
        index_cond_not_in_multi <- which(cond_in_multi == FALSE)
        for(i in 1:length(index_cond_not_in_multi)){
          combos_multi[[index_cond_not_in_multi[i]]] <- c(combos_multi[[index_cond_not_in_multi[i]]], cond)
        }
      }
      p_multi <- sapply(combos_multi, function(i){
        length(which(apply(x[,i], 1, function(y){
          all(y > u[i])})))/nu_x})
      return(c(p_uni, p_multi))
    } 
  }
  else{
    combos_multi <- uncon
    cond_in_multi <- sapply(combos_multi, function(x){cond%in%x})
    if(!any(cond_in_multi)){
      index_cond_not_in_multi <- which(cond_in_multi == FALSE)
      for(i in 1:length(index_cond_not_in_multi)){
        combos_multi[[index_cond_not_in_multi[i]]] <- c(combos_multi[[index_cond_not_in_multi[i]]], cond)
      }
    }
    p_multi <- sapply(combos_multi, function(i){
      length(which(apply(x[,i], 1, function(y){
        all(y > u[i])})))/nu_x})
    return(p_multi)
  }
}

## Get P[X_{j} > u_{j} | X_{i} > u_{i}] for a data set
p_data_surv_multi_new <- function(x, cond, u, uncon){
  #get dimension of data
  d <- dim(x)[2]
  
  #get number of times Xi > ui
  nu_x <- length(which(x[,cond] > u[cond]))
  
  ## Now get the conditional probabilities
  p <- sapply(uncon, function(i){
    length(which(apply(t(x[,c(cond, i)]) > u[c(cond, i)], 2, all) == TRUE))/nu_x})
  return(p)
}

## Get P[X_{j} < u_{j} | X_{i} > u_{i}] for a data set
p_data_cond_cdf <- function(x, cond, u, uncon){
  #get dimension of data
  d <- dim(x)[2]
  n <- dim(x)[1]
  
  #get number of times Xi > ui
  index_large_cond <- which(x[,cond] > u[cond])
  n_ux <- length(index_large_cond)
  
  ## Now get the conditional probabilities
  p <- sapply(uncon, function(j){
    length(which(apply(t(x[index_large_cond,j]) < u[j], 2, all)))/n_ux})
  return(p)
}

## Get P[X_{j} > u_{j} | X_{i} > u_{i}] if X is multivariate Gaussian
p_mvnorm <- function(u, mu, Sigma, cond, uncond){
  
  ## get info from the data
  d <- length(mu)
  
  #Calculate the prob of the conditioning variable being large
  p_cond <- unname(pmvnorm(upper = u[cond], mean = mu[cond], sigma = Sigma[cond,cond])[1])
  
  ## Calculate the prob of the unconditioned variables being large
  d_uncond <- length(c(cond, uncond))
  uncond_final <- lapply(1:d_uncond, function(r){combinations(n = d_uncond, r = r, v = c(cond, uncond))})
  n_uncond_final <- sapply(uncond_final, nrow)
  
  p_uncond <- sapply(1:d_uncond, function(i){sapply(1:n_uncond_final[[i]], function(j){
    unname(pmvnorm(upper = u[uncond_final[[i]][j,]], mean = mu[uncond_final[[i]][j,]], sigma = Sigma[uncond_final[[i]][j,], uncond_final[[i]][j,]])[1])*((-1)^(i))})})
  
  ## Calcualte the overall probability
  p <- (1 + sum(do.call(c, p_uncond)))/(1 - p_cond)
  return(p)
}

## Function to calculate joint exceedance probabilities
# Prob_calc <- function(data, cond, u_p, B, uncon){
#   ## Number of times conditioning variable exceeding the threshold
#   nu_cond <- sapply(1:B, function(i){
#     length(which(data[[i]][,cond] > u_p[cond]))})
#   
#   index_uni <- which(sapply(uncon, function(x){length(x) == 1}))
#   if(!is_empty(index_uni)){
#     
#     ## get the univariate combinations and add in the conditioning random variable if needed
#     combos_uni <- uncon[index_uni]
#     cond_in_uni <- sapply(combos_uni, function(x){cond%in%x})
#     if(!any(cond_in_uni)){
#       index_cond_not_in_uni <- which(cond_in_uni == FALSE)
#       for(i in 1:length(index_cond_not_in_uni)){
#         combos_uni[[index_cond_not_in_uni[i]]] <- sort(c(combos_uni[[index_cond_not_in_uni[i]]], cond))
#       }
#     }
#     
#     ## calculate conditional probabilities
#     p_uni <- lapply(combos_uni, function(i){sapply(1:B, function(j){
#       length(which(apply(data[[j]][,i], 1, function(y){all(y > u_p[i])})))/nu_cond[j]})})  
#     
#     ## calculate multivariate probs if required
#     if(length(index_uni) == length(uncon)){
#       p_all <- p_uni
#     }
#     else{
#       #get the multivariate probabilities
#       combos_multi <- uncon[-index_uni]
#       cond_in_multi <- sapply(combos_multi, function(x){cond%in%x})
#       if(!any(cond_in_multi)){
#         index_cond_not_in_multi <- which(cond_in_multi == FALSE)
#         for(i in 1:length(index_cond_not_in_multi)){
#           combos_multi[[index_cond_not_in_multi[i]]] <- sort(c(combos_multi[[index_cond_not_in_multi[i]]], cond))
#         }
#       }
#       p_multi <- lapply(combos_multi, function(i){
#         sapply(1:B, function(j){
#           length(which(apply(data[[j]][,i], 1, function(x){
#             all(x > u_p[i])})))/nu_cond[j]})})
#       p_all <- c(p_uni, p_multi)
#     }
#   }
#   else{
#     combos_multi <- uncon
#     cond_in_multi <- sapply(combos_multi, function(x){cond%in%x})
#     if(!any(cond_in_multi)){
#       index_cond_not_in_multi <- which(cond_in_multi == FALSE)
#       for(i in 1:length(index_cond_not_in_multi)){
#         combos_multi[[index_cond_not_in_multi[i]]] <- c(combos_multi[[index_cond_not_in_multi[i]]], cond)
#       }
#     }
#     p_all <- lapply(combos_multi, function(i){
#       sapply(1:B, function(j){
#         length(which(apply(X[[j]][,i], 1, function(x){
#           all(x > u_p[i])})))/nu_cond[j]})})
#   }
#   if(length(u_p[[1]]) > 1){
#     p_all <- lapply(p_all, function(x){apply(x, 1, mean)})
#   }
#   else{
#     p_all <- lapply(p_all, mean)
#   }
#   return(p_all)
# }

Prob_calc <- function(data, cond, u_p, B, uncon){
  ## Number of times conditioning variable exceeding the threshold
  nu_cond <- sapply(1:B, function(i){
    length(which(data[[i]][,cond] > u_p[cond]))})
  
  ## probability for each sample
  p <- sapply(1:B, function(i){sapply(uncon, function(j){
    length(which(apply(t(data[[i]][,c(cond, j)]) > u_p[c(cond, j)], 2, all) == TRUE))/nu_cond[[i]]})})
  
  p <- apply(p, 1, mean)
  return(p)
}

## Get P[X_{j} > u_{j} | X_{i} > u_{i}] by sampling from Z|i 
Prob_Joint_Large_Sample_Z <- function(data, cond = 1, n_excesses, Z, a, b, scale_gpd, shape_gpd, u_orig, qu_orig, u_unif, u_p, B = 1, uncon){
  
  ## Get dimension
  d <- dim(data)[2]
  
  ## make data a list for later
  data_list <- lapply(apply(data, 2, list), function(x){x[[1]]})
  
  ## Step 1 - Simulate Yi| Yi > ti(u_{X_{i}}) from a Laplace distribution
  ui <- replicate(n = B, expr = runif(n = n_excesses, min = u_unif), simplify = FALSE)
  Yi <- lapply(ui, qlaplace)
  Xi <- lapply(ui, function(u){qgpd(p = 1 - (1 - u)/(1 - qu_orig[cond]), 
                                    mu = u_orig[cond], 
                                    scale = scale_gpd[cond], 
                                    shape = shape_gpd[cond])})
  
  #Step 2 - Simulate Z|i from G|i
  samps <- replicate(n = B, expr = sample(1:nrow(Z), n_excesses, replace = TRUE), simplify = FALSE)
  Z_star_sample <- lapply(samps, function(i){Z[i,]})
  
  #Step 2 - Obtain Y_{-i} (Y not including i)
  a_yi_hat_sim <- lapply(Yi, function(x){sapply(a, function(a){x*a})})
  b_yi_hat_sim <- lapply(Yi, function(x){sapply(b, function(b){x^b})})
  Y_not_i <- pmap(.l = list(a = a_yi_hat_sim, b = b_yi_hat_sim, z = Z_star_sample),
                  .f = function(a, b, z){a + b*z})
  
  ## Step 3 - transform Y onto X
  U_not_i <- lapply(Y_not_i, function(y){apply(y, 2, plaplace, simplify = FALSE)})
  
  X_not_i <- lapply(U_not_i, function(u){do.call(cbind, pmap(.l = list(x = u,
                                                                       data = data_list[-cond],
                                                                       qu = qu_orig[-cond],
                                                                       th = u_orig[-cond],
                                                                       sigma = scale_gpd[-cond],
                                                                       xi = shape_gpd[-cond]),
                                                             .f = function(x, data, qu, th, sigma, xi){
                                                               texmex:::revTransform(x = x, data = data, qu = qu, th = th,
                                                                                     sigma = sigma, xi = xi, method = "mixture")}))})
  
  ## combine to get a matrix Y and X
  if(cond == 1){
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x1, x2)})
  }
  else if(cond == d){
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x2, x1)})
  }
  else{
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x2[,(1:(cond - 1))], x1, x2[,(cond:(d - 1))])})
  }
  
  #Step 4 - Calculate probabilities
  p_all <- Prob_calc(data = X, cond = cond, u_p = u_p, B = B, uncon = uncon)
  return(p_all)
}

## Get P[X_{j} > u_{j} | X_{i} > u_{i}] by simulating from Z|i - MVN
Prob_Joint_Large_Simulate_MVN <- function(data, cond = 1, n_excesses, a, b, mu, Sigma, scale_gpd, shape_gpd, u_orig, qu_orig, u_unif, u_p, B = 1, uncon){
  
  ## Get dimension
  d <- dim(data)[2]
  
  ## make data a list for later
  data_list <- lapply(apply(data, 2, list), function(x){x[[1]]})
  
  ## Step 1 - Simulate Yi| Yi > ti(u_{X_{i}}) from a Laplace distribution
  ui <- replicate(n = B, expr = runif(n = n_excesses, min = u_unif), simplify = FALSE)
  Yi <- lapply(ui, qlaplace)
  Xi <- lapply(ui, function(u){qgpd(p = 1 - (1 - u)/(1 - qu_orig[cond]), 
                                    mu = u_orig[cond], 
                                    scale = scale_gpd[cond], 
                                    shape = shape_gpd[cond])})
  
  #Step 2 - Simulate Z|i from G|i
  Z_star_simulate <- replicate(n = B, expr = rmvnorm(n = n_excesses, mean = mu, sigma = Sigma), simplify = FALSE)
  
  #Step 2 - Obtain Y_{-i} (Y not including i)
  a_yi_hat_sim <- lapply(Yi, function(x){sapply(a, function(a){x*a})})
  b_yi_hat_sim <- lapply(Yi, function(x){sapply(b, function(b){x^b})})
  Y_not_i <- pmap(.l = list(a = a_yi_hat_sim, b = b_yi_hat_sim, z = Z_star_simulate),
                  .f = function(a, b, z){a + b*z})
  
  ## Step 3 - transform Y onto X
  U_not_i <- lapply(Y_not_i, function(y){apply(y, 2, plaplace, simplify = FALSE)})
  
  X_not_i <- lapply(U_not_i, function(u){do.call(cbind, pmap(.l = list(x = u,
                                                                       data = data_list[-cond],
                                                                       qu = qu_orig[-cond],
                                                                       th = u_orig[-cond],
                                                                       sigma = scale_gpd[-cond],
                                                                       xi = shape_gpd[-cond]),
                                                             .f = function(x, data, qu, th, sigma, xi){
                                                               texmex:::revTransform(x = x, data = data, qu = qu, th = th,
                                                                                     sigma = sigma, xi = xi, method = "mixture")}))})
  
  ## combine to get a matrix Y and X
  if(cond == 1){
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x1, x2)})
  }
  else if(cond == d){
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x2, x1)})
  }
  else{
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x2[,(1:(cond - 1))], x1, x2[,(cond:(d - 1))])})
  }
  
  #Step 4 - Calculate probabilities
  p_all <- Prob_calc(data = X, cond = cond, u_p = u_p, B = B, uncon = uncon)
  return(p_all)
}

## Get P[X_{j} > u_{j} | X_{i} > u_{i}] by simulating from Z|i - MVAGG
Prob_Joint_Large_Simulate_MVAGGD <- function(data, cond = 1, n_excesses, a, b, loc, scale_1, scale_2, shape, Sigma, scale_gpd, shape_gpd, u_orig, qu_orig, u_unif, u_p, B = 1, uncon){
  ## Get dimension
  d <- dim(data)[2]
  ## make data a list for later
  data_list <- lapply(apply(data, 2, list), function(x){x[[1]]})
  ## Step 1 - Simulate Yi| Yi > ti(u_{X_{i}}) from a Laplace distribution
  ui <- replicate(n = B, expr = runif(n = n_excesses, min = u_unif), simplify = FALSE)
  Yi <- lapply(ui, qlaplace)
  Xi <- lapply(ui, function(u){qgpd(p = 1 - (1 - u)/(1 - qu_orig[cond]), 
                                    mu = u_orig[cond],
                                    scale = scale_gpd[cond],
                                    shape = shape_gpd[cond])})
  #Step 2 - Simulate Z|i from G|i
  Z_star_simulate <- replicate(n = B, 
                               expr = rmvagg(loc = loc, scale_1 = scale_1, scale_2 = scale_2, shape = shape, Sigma = Sigma, n = n_excesses, dim = d - 1), 
                               simplify = FALSE)
  
  #Step 2 - Obtain Y_{-i} (Y not including i)
  a_yi_hat_sim <- lapply(Yi, function(x){sapply(a, function(a){x*a})})
  b_yi_hat_sim <- lapply(Yi, function(x){sapply(b, function(b){x^b})})
  Y_not_i <- pmap(.l = list(a = a_yi_hat_sim, b = b_yi_hat_sim, z = Z_star_simulate),
                  .f = function(a, b, z){a + b*z})
  ## Step 3 - transform Y onto X
  U_not_i <- lapply(Y_not_i, function(y){apply(y, 2, plaplace, simplify = FALSE)})
  X_not_i <- lapply(U_not_i, function(u){do.call(cbind, pmap(.l = list(x = u,
                                                                       data = data_list[-cond],
                                                                       qu = qu_orig[-cond],
                                                                       th = u_orig[-cond],
                                                                       sigma = scale_gpd[-cond],
                                                                       xi = shape_gpd[-cond]),
                                                             .f = function(x, data, qu, th, sigma, xi){
                                                               texmex:::revTransform(x = x, data = data, qu = qu, th = th,
                                                                                     sigma = sigma, xi = xi, method = "mixture")}))})
  ## combine to get a matrix Y and X
  if(cond == 1){
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x1, x2)})
  }
  else if(cond == d){
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x2, x1)})
  }
  else{
    X <- pmap(.l = list(x1 = Xi, x2 = X_not_i), .f = function(x1, x2){cbind(x2[,(1:(cond - 1))], x1, x2[,(cond:(d - 1))])})
  }
  ## Step 4 - Calculate probabilities
  p_all <- Prob_calc(data = X, cond = cond, u_p = u_p, B = B, uncon = uncon)
  return(p_all)
}

################################################################################
## Conditional extremes model but we use a three-step approach

Cond_extremes_MVAGGD_Three_step <- function(data, cond = 1, graph = NA,
                                            constrain = TRUE, q = c(0,1), v = 10, aLow = -1, 
                                            maxit = 1e+6, nOptim = 1,
                                            start_HT = c(0.1, 0.1), start_AGG = c(0, 1, 2, 1.5)){
  
  ## make sure conditioning random variable is fine
  if(missing(cond)){
    message("Missing 'cond'. Conditioning on", dimnames(x$transformed)[[2]][1], "\n")
    cond <- 1
  }
  else if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  
  #get information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  dependent <- (1:(dim(data)[[2]]))[-cond]
  
  ## Step 1
  ## Fit the original conditional multivaraite extremes model to the data to
  ## obtain the fitted residuals and dependence parameters
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  fit_HT <- lapply(1:(d-1), function(i){
    qfun_HT(yex = yex, ydep = as.matrix(ydep[,i]),
            constrain = constrain, aLow = aLow, q = q, v = v,
            maxit = maxit, start = start_HT, nOptim = nOptim)})
  
  ## Get the necessary output
  z <- sapply(fit_HT, function(x){x$Z})
  a_hat <- sapply(fit_HT, function(x){x$par$a})
  b_hat <- sapply(fit_HT, function(x){x$par$b})
  
  ## Step 2
  ## Fit an asymmetric generalised Gaussian distribution to the fitted residuals
  ## Do this marginally
  z_mean <- apply(z, 2, mean)
  fit_AGG <- lapply(1:(d-1), function(i){fit_aggd(par = c(z_mean[i], 1.5, 1.5, 1.5), data = z[,i])})
  
  ## Get the necessary output
  loc_hat <- sapply(fit_AGG, function(x){x$par$main[1]})
  scale_1_hat <- sapply(fit_AGG, function(x){x$par$main[2]})
  scale_2_hat <- sapply(fit_AGG, function(x){x$par$main[3]})
  shape_hat <- sapply(fit_AGG, function(x){x$par$main[4]})
  
  ## Step 3
  ## Determine the covariance matrix from the data
  
  ## First transform the data onto standard Gaussian margins
  Z_Gaussian <- sapply(1:(d-1), function(i){
    qnorm(paggd(q = z[,i],
                loc = loc_hat[i],
                scale_1 = scale_1_hat[i],
                scale_2 = scale_2_hat[i],
                shape = shape_hat[i]))})
  
  ## Determine the non_edges in the graph
  if(is_igraph(graph)){
    all_edges <- data.frame(combinations(n = d-1, r = 2, v = 1:(d-1)))
    graph_edges <- data.frame(as_edgelist(delete_vertices(graph, cond)))
    all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, graph_edges)
    non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2]) 
  }
  else{
    non_edges <- combinations(n = d-1, r = 2, v = 1:(d-1))
  }
  
  ## Now fit the graphical lasso
  if(is_empty(non_edges)){
    fit_glasso <- glasso(s = cor(Z_Gaussian), rho = 0, nobs = n, 
                         thr = 1e-8, maxit = 1e+6, penalize.diagonal = FALSE)
  }
  else{
    fit_glasso <- glasso(s = cor(Z_Gaussian), rho = 0, nobs = n, 
                         zero = non_edges, thr = 1e-8, maxit = 1e+6, penalize.diagonal = FALSE) 
  }
  
  ## Get the output from the model
  out <- list()
  
  par_out <- rbind(a_hat, b_hat, loc_hat, scale_1_hat, scale_2_hat, shape_hat)
  par_out <- as.data.frame(par_out)
  rownames(par_out) <- c("alpha", "beta", "loc", "scale_1", "scale_2", "shape")
  colnames(par_out) <- sapply(dependent, function(i){paste0("Column", i)})
  
  out$par <- list(main = par_out, Sigma = fit_glasso$w, Gamma = fit_glasso$wi)
  
  out_Z <- as.data.frame(z)
  colnames(out_Z) <- sapply(dependent, function(i){paste0("Column", i)})
  out$Z <- out_Z
  return(out)
}

################################################################################
## Simulation of the entire spatial process
## See Richards et al. 2022 Algorithm 1
Sim_Surface_HT <- function(data, data_Laplace, n_sim, q, a, b, Z,
                           qu_gpd, u_gpd, scale_gpd, shape_gpd){
  if(!is.matrix(data)){
    stop("Data must be a matrix is rows corresponding to replicates of the process
         and columns corresponding to the dimension of the process.")
  }
  ## Get information about the data set
  n <- nrow(data)
  d <- ncol(data)
  ## Get the conditioning sites
  cond_sites <- sample(x = 1:d, size = n_sim, replace = TRUE)
  ## Convert the threshold onto Laplace margins and simulate from the conditioning site
  v <- qlaplace(q)
  large_Laplace <- v +  rexp(n = n_sim, rate = 1)
  ## Simulate residuals from the fitted model
  Z_sample <- t(sapply(cond_sites, function(i){
    unlist(unname(Z[[i]][sample(1:nrow(Z[[i]]), 1),]))}))
  ## Convert the simulated data onto Laplace margins
  a_cond_sites <- t(sapply(cond_sites, function(i){a[[i]]}))
  b_cond_sites <- t(sapply(cond_sites, function(i){b[[i]]}))
  Laplace_Samples <- matrix(NA, nrow = n_sim, ncol = d-1)
  for(i in 1:(d-1)){
    Laplace_Samples[,i] <- a_cond_sites[,i]*large_Laplace + (large_Laplace^(b_cond_sites[,i]))*Z_sample[,i]
  }
  ## Now combine Y{-i}|Y_{i} > u_{Y_{i}} and Y{-i}|Y_{i} > u_{Y_{i}}
  Laplace_Samples_Full <- matrix(NA, nrow = n_sim, ncol = d)
  for(i in 1:n_sim){
    r <- cond_sites[i]
    Laplace_Samples_Full[i,r] <- large_Laplace[i]
    if(r == 1){
      Laplace_Samples_Full[i,-1] <- Laplace_Samples[i,]
    }
    else if(r == d){
      Laplace_Samples_Full[i,-d] <- Laplace_Samples[i,]
    }
    else{
      Laplace_Samples_Full[i,(1:(r-1))] <- Laplace_Samples[i,(1:(r-1))]
      Laplace_Samples_Full[i,((r+1):d)] <- Laplace_Samples[i,(r:(d-1))]
    }
  }
  ## Now use importance sampling to resample the simulated data
  ## Idea is to up-weight those close to the boundary and down-weight those in the
  ## centre of the spatial process
  weights <- d/apply(Laplace_Samples_Full, 1, function(x){length(which(x > v))})
  Index <- sample(x = 1:n_sim, size = n, prob = weights/sum(weights))
  ## Now we convert the data onto the original scale
  Final_Laplace_Samples <- Laplace_Samples_Full[Index,]
  ## convert onto original scale
  Final_Uniform_Samples <- apply(Final_Laplace_Samples, 2, plaplace, simplify = FALSE)
  data_list <- lapply(apply(data, 2, list), function(x){x[[1]]})
  Final_Original_Samples <- mcmapply(FUN = texmex:::revTransform,
                                     x = Final_Uniform_Samples,
                                     data = data_list,
                                     qu = qu_gpd,
                                     th = u_gpd,
                                     sigma = scale_gpd,
                                     xi = shape_gpd,
                                     MoreArgs = list(method = "mixture"),
                                     SIMPLIFY = TRUE,
                                     mc.cores = detectCores() - 1)
  ## Above simulated X(s) | max(X(s)) > v
  ## We want to get the unconditioned process
  ## First get the data points were we have no extremes
  Index_No_Extremes <- which(apply(data_Laplace, 1, max) < v)
  Data_Body <- data[Index_No_Extremes,]
  Data_Body_Laplace <- data_Laplace[Index_No_Extremes,]
  ## then get the probability we will draw from the joint tail
  Index_At_Least_One_Extreme <- which(apply(data_Laplace, 1, max) > v)
  p_tail <- length(Index_At_Least_One_Extreme)/(n + 1)
  p_accecpt <- runif(n = n)
  ## Get the index of body/tail
  Index_Body <- sample(x = 1:length(Index_No_Extremes), size = length(which(p_accecpt >= p_tail)), replace = TRUE)
  Index_Tail <- which(p_accecpt < p_tail)
  ## Get the final data sets
  Data_Final_Laplace_Margins <- rbind(Data_Body_Laplace[Index_Body,],
                                      Final_Laplace_Samples[Index_Tail,])
  Data_Final_Original_Margins <- rbind(Data_Body[Index_Body,],
                                       Final_Original_Samples[Index_Tail,])
  return(list(Data_Margins = Data_Final_Original_Margins,
              Laplace_Margins = Data_Final_Laplace_Margins))
}

## As above but we simulate from the fitted model for the residuals rather than
## sampling from the fitted residuals
Sim_Surface_MVAGG <- function(data, data_Laplace, n_sim, q, a, b, loc, scale_1, scale_2, shape, Sigma,
                              qu_gpd, u_gpd, scale_gpd, shape_gpd){
  if(!is.matrix(data)){
    stop("Data must be a matrix is rows corresponding to replicates of the process
         and columns corresponding to the dimension of the process.")
  }
  ## Get information about the data set
  n <- nrow(data)
  d <- ncol(data)
  ## Get the conditioning sites
  cond_sites <- sample(x = 1:d, size = n_sim, replace = TRUE)
  ## Convert the threshold onto Laplace margins and simulate from the conditioning site
  v <- qlaplace(q)
  large_Laplace <- v +  rexp(n = n_sim, rate = 1)
  ## Simulate residuals from the fitted model
  Z <- t(mcmapply(FUN = rmvagg,
                  loc = lapply(cond_sites, function(i){loc[[i]]}),
                  scale_1 = lapply(cond_sites, function(i){scale_1[[i]]}),
                  scale_2 = lapply(cond_sites, function(i){scale_2[[i]]}),
                  shape = lapply(cond_sites, function(i){shape[[i]]}),
                  Sigma = lapply(cond_sites, function(i){Sigma[[i]]}),
                  MoreArgs = list(n = 1, dim = d-1),
                  SIMPLIFY = TRUE,
                  mc.cores = detectCores() - 1))
  ## Convert the simulated data onto Laplace margins
  a_cond_sites <- t(sapply(cond_sites, function(i){a[[i]]}))
  b_cond_sites <- t(sapply(cond_sites, function(i){b[[i]]}))
  Laplace_Samples <- matrix(NA, nrow = n_sim, ncol = d-1)
  for(i in 1:(d-1)){
    Laplace_Samples[,i] <- a_cond_sites[,i]*large_Laplace + (large_Laplace^(b_cond_sites[,i]))*Z[,i]
  }
  ## Now combine Y{-i}|Y_{i} > u_{Y_{i}} and Y{-i}|Y_{i} > u_{Y_{i}}
  Laplace_Samples_Full <- matrix(NA, nrow = n_sim, ncol = d)
  for(i in 1:n_sim){
    r <- cond_sites[i]
    Laplace_Samples_Full[i,r] <- large_Laplace[i]
    if(r == 1){
      Laplace_Samples_Full[i,-1] <- Laplace_Samples[i,]
    }
    else if(r == d){
      Laplace_Samples_Full[i,-d] <- Laplace_Samples[i,]
    }
    else{
      Laplace_Samples_Full[i,(1:(r-1))] <- Laplace_Samples[i,(1:(r-1))]
      Laplace_Samples_Full[i,((r+1):d)] <- Laplace_Samples[i,(r:(d-1))]
    }
  }
  ## Now use importance sampling to resample the simulated data
  ## Idea is to up-weight those close to the boundary and down-weight those in the
  ## centre of the spatial process
  weights <- d/apply(Laplace_Samples_Full, 1, function(x){length(which(x > v))})
  Index <- sample(x = 1:n_sim, size = n, prob = weights/sum(weights))
  ## Now we convert the data onto the original scale
  Final_Laplace_Samples <- Laplace_Samples_Full[Index,]
  ## convert onto original scale
  Final_Uniform_Samples <- apply(Final_Laplace_Samples, 2, plaplace, simplify = FALSE)
  data_list <- lapply(apply(data, 2, list), function(x){x[[1]]})
  Final_Original_Samples <- mcmapply(FUN = texmex:::revTransform,
                                     x = Final_Uniform_Samples,
                                     data = data_list,
                                     qu = qu_gpd,
                                     th = u_gpd,
                                     sigma = scale_gpd,
                                     xi = shape_gpd,
                                     MoreArgs = list(method = "mixture"),
                                     SIMPLIFY = TRUE,
                                     mc.cores = detectCores() - 1)
  ## Above simulated X(s) | max(X(s)) > v
  ## We want to get the unconditioned process
  ## First get the data points were we have no extremes
  Index_No_Extremes <- which(apply(data_Laplace, 1, max) < v)
  Data_Body <- data[Index_No_Extremes,]
  Data_Body_Laplace <- data_Laplace[Index_No_Extremes,]
  ## then get the probability we will draw from the joint tail
  Index_At_Least_One_Extreme <- which(apply(data_Laplace, 1, max) > v)
  p_tail <- length(Index_At_Least_One_Extreme)/(n + 1)
  p_accecpt <- runif(n = n)
  ## Get the index of body/tail
  Index_Body <- sample(x = 1:length(Index_No_Extremes), size = length(which(p_accecpt >= p_tail)), replace = TRUE)
  Index_Tail <- which(p_accecpt < p_tail)
  ## Get the final data sets
  Data_Final_Laplace_Margins <- rbind(Data_Body_Laplace[Index_Body,],
                                      Final_Laplace_Samples[Index_Tail,])
  Data_Final_Original_Margins <- rbind(Data_Body[Index_Body,],
                                       Final_Original_Samples[Index_Tail,])
  return(list(Data_Margins = Data_Final_Original_Margins,
              Laplace_Margins = Data_Final_Laplace_Margins))
}
