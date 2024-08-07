## Functions to fit the CMEVM when we assume that the residuals are multivariate Gaussian
## with some graphical structure
Cond_Extremes_MVN <- function(data, cond, graph = NA, 
                              constrain = TRUE, q = c(0,1), v = 20, aLow = -1, 
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
  
  ## check the conditioning random variable is valid and get dependent random variable
  if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(cond%%1 != 0 | cond <= 0 | cond > d){
    stop("cond must be a single positie integer")
  }
  dependent <- (1:d)[-cond]
  
  ## get the starting parameters
  if(missing(start)){
    start <- c(0.1, 0.1)
  }
  else if(!is.numeric(start)){
    stop("start must be a vector")
  }
  else if(length(start) != 2){
    stop("start must be a vector of length 2")
  }
  else if(any(abs(start) > 1)){
    stop("Initial starting values are outside the parameter sapce")
  }
  if(length(start) == 2){
    start <- matrix(rep(start, d-1), ncol = 2)
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- data[,cond]
  ydep <- data[,-cond]
  
  ## Determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## Fit the independence model
    res <- lapply(1:(d-1), function(i){
      qfun_MVN_indep(yex = yex, ydep = as.matrix(ydep[,i]), 
                     constrain = constrain, aLow = aLow, q = q, v = v,
                     maxit = maxit, start = start[i,], nOptim = nOptim)})
    
    ## Extract the output
    out <- list()
    
    out$par$main <- matrix(data = c(do.call(c, lapply(res, function(x){x$par$a})),
                                    do.call(c, lapply(res, function(x){x$par$b})),
                                    do.call(c, lapply(res, function(x){x$par$mu}))),
                           nrow = 3, ncol = d-1, byrow = TRUE)
    if(any(is.na(out$par$main))){
      out$par$main <- matrix(NA, nrow = 3, ncol = d-1)
    }
    rownames(out$par$main) <- c("a", "b", "mu")
    colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
    
    if(any(is.na(out$par$main))){
      out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = rep(NA, d-1)) 
    }
    else{
      out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = sapply(res, function(x){x$par$Gamma}))
    }
    
    out$loglike <- sum(sapply(res, function(x){-x$value}))
    out$convergence <- max(sapply(res, function(x){x$convergence}))
    out$Z <- do.call(cbind, lapply(res, function(x){x$Z}))
    if(any(is.na(out$Z))){
      out$Z <- matrix(NA, nrow = nrow(out$Z), ncol = ncol(out$Z)) 
    }
    colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
  }
  else{
    ## Graph is provided so we need to determine its structure
    graph_cond <- delete_vertices(graph, cond)
    
    ## Check if the graph is full, if so fit then model instantly
    all_edges <- as.data.frame(combinations(n = d-1, r = 2, v = 1:(d-1)))
    g_edges <- as.data.frame(as_edgelist(graph_cond))
    n_edges_full <- nrow(all_edges)
    n_edges_graph_cond <- nrow(g_edges)
    
    if(n_edges_full == n_edges_graph_cond){
      ## Fit the saturated model
      res <- qfun_MVN_full(yex = yex, ydep = as.matrix(ydep),
                           maxit = maxit, start = c(start), nOptim = nOptim)
      
      ## Extract the output
      out <- list()
      
      out$par$main <- matrix(data = c(res$par$a, res$par$b, res$par$mu),
                             nrow = 3, ncol = d-1, byrow = TRUE)
      rownames(out$par$main) <- c("a", "b", "mu")
      colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
      
      out$par$Gamma <- as(res$par$Gamma, "sparseMatrix")
      
      out$loglike <- -res$value
      out$convergence <- res$convergence
      out$Z <- res$Z
      colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
    }
    else{
      if(is_connected(graph_cond)){
        ## If the graph is connected use the graphical model
        ## determine the edges in the conditional graph
        all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, g_edges)
        non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        
        ## Fit the graphical model
        res <- qfun_MVN_Gamma_graph(yex = yex, ydep = ydep, Gamma_zero = non_edges,
                                    maxit = maxit, start = c(start), nOptim = nOptim)
        
        ## Extract the output
        out <- list()
        
        out$par$main <- matrix(data = c(res$par$a, res$par$b, res$par$mu),
                               nrow = 3, ncol = d-1, byrow = TRUE)
        rownames(out$par$main) <- c("a", "b", "mu")
        colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
        
        out$par$Gamma <- as(res$par$Gamma, "sparseMatrix")
        
        out$loglike <- -res$value
        out$convergence <- res$convergence
        out$Z <- res$Z
        colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
      }
      else{
        ## The graph is disconnected so we can treat each components as
        ## exactly independent
        comps <- components(graph_cond)
        n_comps <- comps$no
        v_comps <- groups(comps)
        g_comps <- lapply(v_comps, function(i){subgraph(graph_cond, i)})
        res_1 <- list()
        for(i in 1:n_comps){
          ## determine if the component is a full graph or not
          n_vertices <- length(V(g_comps[[i]]))
          if(n_vertices == 1){
            res_1[[i]] <- qfun_MVN_indep(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]), 
                                         constrain = TRUE, aLow = aLow, q = q, v = v,
                                         maxit = maxit, start = c(start[v_comps[[i]],]), nOptim = nOptim)
          }
          else{
            all_edges_comp <- as.data.frame(combinations(n = n_vertices, r = 2, v = 1:n_vertices))
            g_edges_comp <- as.data.frame(as_edgelist(g_comps[[i]]))
            ## fit saturated model to this component
            if(nrow(all_edges_comp) == nrow(g_edges_comp)){
              res_1[[i]] <- qfun_MVN_full(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]),
                                          maxit = maxit, 
                                          start = c(start[v_comps[[i]],]), 
                                          nOptim = nOptim)
            }
            else{
              ## fit graphical model to this component
              all_edges_comp$exists <- do.call(paste0, all_edges_comp) %in% do.call(paste0, g_edges_comp)
              non_edges_comp <- as.matrix(all_edges_comp[which(all_edges_comp$exists == FALSE), 1:2])
              res_1[[i]] <- qfun_MVN_graph(yex = yex, ydep = ydep[,v_comps[[i]]], 
                                           Gamma_zero = non_edges_comp,
                                           maxit = maxit, 
                                           start = c(start[v_comps[[i]],]), 
                                           nOptim = nOptim)
            } 
          }
        }
        
        ## Extract the output
        out <- list()
        
        out$par$main <- matrix(data = c(do.call(c, lapply(res_1, function(x){unname(x$par$a)})), 
                                        do.call(c, lapply(res_1, function(x){unname(x$par$b)})),
                                        do.call(c, lapply(res_1, function(x){unname(x$par$mu)}))),
                               nrow = 3, ncol = d-1, byrow = TRUE)
        if(any(is.na(out$par$main))){
          out$par$main <- matrix(NA, nrow = 3, ncol = d-1)
        }
        rownames(out$par$main) <- c("a", "b", "mu")
        colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
        
        Gamma_intrim <- matrix(0, nrow = d-1, ncol = d-1)
        Gamma_entries <- lapply(v_comps, function(x){permutations(n = length(x), r = 2, v = x, repeats.allowed = TRUE)})
        for(i in 1:n_comps){
          Gamma_intrim[Gamma_entries[[i]]] <- res_1[[i]]$par$Gamma
        }
        if(any(is.na(Gamma_intrim))){
          Gamma_intrim <- matrix(NA, nrow = d-1, ncol = d-1)
        }
        out$par$Gamma <- as(Gamma_intrim, "sparseMatrix")
        
        out$loglike <- sum(sapply(res_1, function(x){-x$value}))
        out$convergence <- max(sapply(res_1, function(x){x$convergence}))
        out$Z <- do.call(cbind, lapply(1:n_comps, function(i){res_1[[i]]$Z}))
        if(any(is.na(out$Z))){
          out$Z <- matrix(NA, nrow = nrow(out$Z), ncol = ncol(out$Z)) 
        }
        colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
      }
    }
  }
  class(out) <- "Cond_Extremes_MVN"
  return(out)
}

qfun_MVN_indep <- function(yex, ydep, constrain, q, v, aLow, maxit, start, nOptim){
  
  ## Function for the optim to optimise over
  Qpos <- function(param, yex, ydep, constrain, q, v, aLow, negative = FALSE){
    
    ## Extract starting parameters
    a <- param[1]
    b <- param[2]
    if(a < aLow | abs(a) > 1 | b > 1){
      return((-10^10)*(-1)^negative)
    }
    
    ## Obtain the residuals
    a_yi <- a*yex
    b_yi <- yex^b
    Z <- (ydep - a_yi)/b_yi
    if(any(is.infinite(Z)) | any(is.na(Z)) | any(is.nan(Z))){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## Obtain the mean and covariance matrix of the residuals
      mu_Z <- mean(Z)
      sigma_Z <- as.numeric(var(Z))
      
      ## Fit the profile log-likelihood
      mu <- a_yi + b_yi*mu_Z
      sigma <- sigma_Z*(b_yi^2)
      res <- sum(dnorm(ydep, mean = mu, sd = sqrt(sigma), log = TRUE))
      
      ## Check the value is valid
      if(is.infinite(res)){
        return((-10^10)*(-1)^negative)
        warning("Infinite value of Q in mexDependence")
      }
      ## Checks if the constraints are satisfied to reduce the parameter space
      else if(constrain){ 
        zpos <- quantile(ydep - yex, probs = q)
        z <- quantile(Z, probs = q)
        zneg <- quantile(ydep + yex, probs = q)
        constraints_sat <- texmex:::ConstraintsAreSatisfied(a = a, b = b, z = z, zpos = zpos, zneg = zneg, v = v)
        if(!all(constraints_sat)){
          return((-10^10)*(-1)^negative)
        }
      } 
    }
    ## Obtain the value of the profile log_likelihood for fixed a and b
    return(res*(-1)^negative)
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   yex = yex, ydep = ydep, negative = TRUE,
                   constrain = constrain, aLow = aLow, q = q, v = v, method = "BFGS"),
             silent = TRUE)
  n <- length(yex)
  
  if(nOptim > 1){
    if(inherits(fit, "try-error")){
      next()
    }
    else if(fit$convergence != 0 | fit$value == 1e+10){
      next()
    }
    else{
      for(i in 2:nOptim){
        fit <- try(optim(par = par_start <- fit$par, fn = Qpos, control = list(maxit = maxit),
                         yex = yex, ydep = ydep, negative = TRUE,
                         constrain = constrain, aLow = aLow, v = v, q = q, method = "BFGS"),
                   silent = TRUE)
        if(inherits(fit, "try-error")){
          break()
        }
        else if(fit$convergence != 0 | fit$value == 1e+10){
          break()
        }
      } 
    }
  }
  
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_Extremes_MVN")
    out <- list()
    out$par <- list(a = NA, b = NA, mu = NA, Sigma = NA)
    out$value <- NA
    out$convergence <- NA
    out$Z <- matrix(NA, nrow = n, ncol = 1)
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in Cond_Extremes_MVN")
    out <- list()
    out$par <- list(a = NA, b = NA, mu = NA, Sigma = NA)
    out$value <- NA
    out$convergence <- NA
    out$Z <- matrix(NA, nrow = n, ncol = 1)
  }
  else if(!is.na(fit$par[1])){
    ## Extract MLEs of alpha and beta
    alpha_hat <- fit$par[1]
    beta_hat <- fit$par[2]
    
    ## Obtain the residuals
    a_yi_hat <- yex*alpha_hat
    b_yi_hat <- yex^beta_hat
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## Organise the output
    out <- list()
    out$par <- list(a = alpha_hat, b = beta_hat, mu = mean(Z), Gamma = 1/as.numeric(var(Z)))
    out$value <- fit$value
    out$convergence <- fit$convergence
    out$Z <- Z
  }
  else{
    warning("Unknown error in Cond_Extremes_MVN")
    out <- list()
    out$par <- list(a = NA, b = NA, mu = NA, Sigma = NA)
    out$value <- NA
    out$convergence <- NA
    out$Z <- matrix(NA, nrow = n, ncol = 1)
  }
  class(out) <- "Cond_Extremes_MVN"
  return(out)
}

qfun_MVN_graph <- function(yex, ydep, Gamma_zero, maxit, start, nOptim){
  
  ## Function for the optim to optimise over
  Qpos <- function(param, yex, ydep, Gamma_zero, negative = FALSE){
    ## Extract the the starting parameters
    d <- ncol(ydep)
    n <- length(yex)
    a <- param[1:d]
    b <- param[(d + 1):(2*d)]
    if(any(abs(a) >= 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## Obtain the residuals
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      ## Use a graphical lasso to fit the graphical structure to the residuals
      Sigma_start <- cov(z)
      if(any(is.infinite(Sigma_start)) | any(is.na(Sigma_start)) | any(is.nan(Sigma_start))){
        return((-10^10)*(-1)^negative)
      }
      else if(any(diag(Sigma_start) <= 0) | any(eigen(Sigma_start)$values <= 0)){
        return((-10^10)*(-1)^negative)
      }
      else{
        ## Calculate the log likelihood function
        Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = TRUE, thr = 1e-9, zero = Gamma_zero))
        chol_Gamma <- chol(Lasso_est$wi)
        mu <- matrix(rep(apply(z, 2, mean), n), ncol = d, byrow = TRUE)
        y <- chol_Gamma%*%t(z - mu)
        res <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2 - sum(log(b_yi))
        
        ## Check the value is valid
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
  
  ## Extract the output
  d <- ncol(ydep)
  n <- nrow(ydep)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Error in optim call from Cond_Extremes_MVN")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep, Gamma_zero = Gamma_zero,
                       negative = TRUE,  control = list(maxit = maxit), method = "BFGS"),
                 silent = TRUE)
      
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_Extremes_MVN")
        out <- list()
        out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
        out$Z <- matrix(NA, nrow = n, ncol = d)
        out$value <- NA
        out$convergence <- NA
      }
      else if(fit$convergence != 0 | fit$value == 1e+10){
        warning("Error in optim call from Cond_Extremes_MVN")
        out <- list()
        out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
        out$Z <- matrix(NA, nrow = n, ncol = d)
        out$value <- NA
        out$convergence <- NA
      }
    }
  }
  else if(!is.na(fit$par[1])){
    ## Extract MLEs
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[(d + 1):(2*d)]
    
    ## Obtain the residuals
    a_yi_hat <- sapply(a_hat, function(x){yex*x})
    b_yi_hat <- sapply(b_hat, function(x){yex^x})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## Extract MLE of Gamma
    Lasso_est <- suppressWarnings(glasso(s = cov(Z), rho = 0, penalize.diagonal = TRUE, thr = 1e-9, zero = Gamma_zero))
    
    ## Organise the output
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, mu = apply(Z, 2, mean), Gamma = Lasso_est$wi)
    out$value <- fit$value
    out$convergence <- fit$convergence
    out$Z <- Z
  }
  else{
    warning("Unknown error in Cond_Extremes_MVN")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  return(out)
}

qfun_MVN_full <- function(yex, ydep, maxit, start, nOptim){
  
  ## Function for the optim to optimise over
  Qpos <- function(param, yex, ydep, negative = FALSE){
    ## Extract the the starting parameters
    d <- ncol(ydep)
    n <- length(yex)
    a <- param[1:d]
    b <- param[(d + 1):(2*d)]
    if(any(abs(a) >= 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## Obtain the residuals
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      ## Now maximise assuming the residuals have a saturated Gaussian model
      mu_z <- matrix(rep(apply(z, 2, mean), n), ncol = d, byrow = TRUE)
      chol_Gamma <- chol(solve(cov(z)))
      y <- chol_Gamma%*%t(z - mu_z)
      llh_z <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
      res <- (llh_z - sum(log(b_yi)))
      
      ## Check the value is valid
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
  
  ## Extract the output
  d <- ncol(ydep)
  n <- length(yex)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_Extremes_MVN")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Error in optim call from Cond_Extremes_MVN")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                       negative = TRUE,  control = list(maxit = maxit), method = "BFGS"),
                 silent = TRUE)
      
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_Extremes_MVN")
        out <- list()
        out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
        out$Z <- matrix(NA, nrow = n, ncol = d)
        out$value <- NA
        out$convergence <- NA
      }
      else if(fit$convergence != 0 | fit$value == 1e+10){
        warning("Error in optim call from Cond_Extremes_MVN")
        out <- list()
        out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
        out$Z <- matrix(NA, nrow = n, ncol = d)
        out$value <- NA
        out$convergence <- NA
      }
    }
  }
  else if(!is.na(fit$par[1])){
    ## Extract MLEs of alpha and beta
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[-(1:d)]
    
    ## Obtain the residuals
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
  else{
    warning("Unknown error in Cond_Extremes_MVN")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), mu = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  return(out)
}
