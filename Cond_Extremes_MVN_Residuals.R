## Functions to fit the CMEVM when we assume that the residuals are multivariate Gaussian
## with some graphical structure
Cond_extremes_graph_new <- function(data, cond, graph = NA, 
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
  
  ## determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## Fit the model
    res <- lapply(1:(d-1), function(i){
      qfun_MVN_indep(yex = yex, ydep = as.matrix(ydep[,i]), 
                     constrain = constrain, aLow = aLow, q = q, v = v,
                     maxit = maxit, start = start[i,], nOptim = nOptim)})
    
    #get the output
    out <- list()
    
    out$par$main <- matrix(data = c(do.call(c, lapply(res, function(x){x$par$a})),
                                    do.call(c, lapply(res, function(x){x$par$b})),
                                    do.call(c, lapply(res, function(x){x$par$mu}))),
                           nrow = 3, ncol = d-1, byrow = TRUE)
    rownames(out$par$main) <- c("a", "b", "mu")
    colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
    
    out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = sapply(res, function(x){x$par$Gamma}))
    
    out$loglike <- sum(sapply(res, function(x){-x$value}))
    out$convergence <- max(sapply(res, function(x){x$convergence}))
    out$Z <- do.call(cbind, lapply(res, function(x){x$Z}))
  }
  else{
    #Graph is provided so we need to figure out the separators and the cliques
    graph_cond <- delete_vertices(graph, cond)
    
    ## check if the graph is full, if so fit the model instantly
    all_edges <- as.data.frame(combinations(n = d-1, r = 2, v = 1:(d-1)))
    g_edges <- as.data.frame(as_edgelist(graph_cond))
    n_edges_full <- nrow(all_edges)
    n_edges_graph_cond <- nrow(g_edges)
    
    if(n_edges_full == n_edges_graph_cond){
      ## Fit the model
      res <- qfun_MVN_full(yex = yex, ydep = as.matrix(ydep),
                           maxit = maxit, start = c(start), nOptim = nOptim)
      
      #get the output
      out <- list()
      
      out$par$main <- matrix(data = c(res$par$a, res$par$b, res$par$mu),
                             nrow = 3, ncol = d-1, byrow = TRUE)
      rownames(out$par$main) <- c("a", "b", "mu")
      colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
      
      lower_tri_elements <- lower.tri(res$par$Gamma, diag = TRUE)
      non_zero_elements <- which(res$par$Gamma*lower_tri_elements != 0, arr.ind = TRUE)
      out$par$Gamma <- sparseMatrix(i = non_zero_elements[,1], j = non_zero_elements[,2], x = c(res$par$Gamma[non_zero_elements]), symmetric = TRUE)
      
      out$loglike <- -res$value
      out$convergence <- res$convergence
      out$Z <- res$Z
    }
    else{
      ## if the graph is connected use the graphical model
      if(is_connected(graph_cond)){
        
        ## determine the edges in the model
        all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, g_edges)
        non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        
        ## Fit the model
        res <- qfun_MVN_Gamma_graph(yex = yex, ydep = ydep, Gamma_zero = non_edges,
                                    maxit = maxit, start = c(start), nOptim = nOptim)
        
        #get the output
        out <- list()
        
        out$par$main <- matrix(data = c(res$par$a, res$par$b, res$par$mu),
                               nrow = 3, ncol = d-1, byrow = TRUE)
        rownames(out$par$main) <- c("a", "b", "mu")
        colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
        
        lower_tri_elements <- lower.tri(res$par$Gamma, diag = TRUE)
        non_zero_elements <- which(res$par$Gamma*lower_tri_elements != 0, arr.ind = TRUE)
        out$par$Gamma <- sparseMatrix(i = non_zero_elements[,1], j = non_zero_elements[,2], x = c(res$par$Gamma[non_zero_elements]), symmetric = TRUE)
        
        out$loglike <- -res$value
        out$convergence <- res$convergence
        out$Z <- res$Z
      }
      else{
        ## the graph is disconnected so we can treat the components as exactly independent
        ## extract the components
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
        
        #get the output
        out <- list()
        
        out$par$main <- matrix(data = c(do.call(c, lapply(res_1, function(x){unname(x$par$a)})), 
                                        do.call(c, lapply(res_1, function(x){unname(x$par$b)})),
                                        do.call(c, lapply(res_1, function(x){unname(x$par$mu)}))),
                               nrow = 3, ncol = d-1, byrow = TRUE)
        rownames(out$par$main) <- c("a", "b", "mu")
        colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
        
        Gamma_intrim <- matrix(0, nrow = d-1, ncol = d-1)
        Gamma_entries <- lapply(v_comps, function(x){permutations(n = length(x), r = 2, v = x, repeats.allowed = TRUE)})
        for(i in 1:n_comps){
          Gamma_intrim[Gamma_entries[[i]]] <- res_1[[i]]$par$Gamma
        }
        lower_tri_elements <- lower.tri(Gamma_intrim, diag = TRUE)
        non_zero_elements <- which(Gamma_intrim*lower_tri_elements != 0, arr.ind = TRUE)
        out$par$Gamma <- sparseMatrix(i = non_zero_elements[,1], j = non_zero_elements[,2], x = c(Gamma_intrim[non_zero_elements]), symmetric = TRUE)
        
        out$loglike <- sum(sapply(res_1, function(x){-x$value}))
        out$convergence <- max(sapply(res_1, function(x){x$convergence}))
        out$Z <- do.call(cbind, lapply(1:n_comps, function(i){res_1[[i]]$Z}))
      }
    }
  }
  class(out) <- "Cond_extremes_MVN"
  return(out)
}

qfun_MVN_indep <- function(yex, ydep, constrain, q, v, aLow, maxit, start, nOptim){
  
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, constrain, q, v, aLow, negative = FALSE){
    
    #get the starting parameters
    a <- param[1]
    b <- param[2]
    
    ## check alpha and beta parameters
    if(a < aLow | abs(a) > 1 | b > 1){
      return((-10^10)*(-1)^negative)
    }
    
    #get the residuals
    a_yi <- a*yex
    b_yi <- yex^b
    Z <- (ydep - a_yi)/b_yi
    
    #get the mean and covariance matrix of the residuals
    mu_Z <- mean(Z)
    sigma_Z <- as.numeric(var(Z))
    
    ## Fit the profile log-likelihood
    mu <- a_yi + b_yi*mu_Z
    sigma <- sigma_Z*(b_yi^2)
    res <- sum(dnorm(ydep, mean = mu, sd = sqrt(sigma), log = TRUE))
    
    #check the value is valid
    if(is.infinite(res)){
      return((-10^10)*(-1)^negative)
      warning("Infinite value of Q in mexDependence")
    }
    #Checks if the constraints are satisfied to reduce the parameter space
    else if(constrain){ 
      zpos <- quantile(ydep - yex, probs = q)
      z <- quantile(Z, probs = q)
      zneg <- quantile(ydep + yex, probs = q)
      constraints_sat <- texmex:::ConstraintsAreSatisfied(a = a, b = b, z = z, zpos = zpos, zneg = zneg, v = v)
      if(!all(constraints_sat)){
        return((-10^10)*(-1)^negative)
      }
    }
    
    #get the value of the profile log_likelihood for fixed a and b
    return(res*(-1)^negative)
  }
  
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   yex = yex, ydep = ydep, negative = TRUE,
                   constrain = constrain, aLow = aLow, q = q, v = v, method = "BFGS"),
             silent = TRUE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_graph")
    out <- list()
    out$par <- list(a = NA, b = NA, mu = NA, Sigma = NA)
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
    return(out)
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_graph")
    out <- list()
    out$par <- list(a = NA, b = NA, mu = NA, Sigma = NA)
    out$value <- NA
    out$convergence <- NA
    out$Z <- NA
    return(out)
  }
  else if(nOptim > 1){
    for(i in 2:nOptim){
      par_start <- fit$par
      par_start[which(par_start < 0)] <- 0.1
      #par_start[which(par_start < 0)] <- mean(par_start[which(par_start > 0)])
      fit <- try(optim(par = par_start, fn = Qpos, control = list(maxit = maxit),
                       yex = yex, ydep = ydep, negative = TRUE,
                       constrain = constrain, aLow = aLow, v = v, q = q, method = "BFGS"),
                 silent = TRUE)
      if(inherits(fit, "try-error")){
        warning("Error in optim call from Cond_extremes_graph")
        out <- list()
        out$par <- list(a = NA, b = NA, mu = NA, Sigma = NA)
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
        return(out)
      }
      else if(fit$convergence != 0){
        warning("Non-convergence in Cond_extremes_graph")
        out <- list()
        out$par <- list(a = NA, b = NA, mu = NA, Sigma = NA)
        out$value <- NA
        out$convergence <- NA
        out$Z <- NA
        return(out)
      }
    }
  }
  #Get the output
  if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    alpha_hat <- fit$par[1]
    beta_hat <- fit$par[2]
    
    #obtain the residuals
    a_yi_hat <- yex*alpha_hat
    b_yi_hat <- yex^beta_hat
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    out <- list()
    out$par <- list(a = alpha_hat, b = beta_hat, mu = mean(Z), Gamma = 1/as.numeric(var(Z)))
    out$value <- fit$value
    out$convergence <- fit$convergence
    out$Z <- Z
  }
  return(out)
}

qfun_MVN_graph <- function(yex, ydep, Gamma_zero, maxit, start, nOptim){
  
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
    return(out)
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
    return(out)
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
        return(out)
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
        return(out)
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

qfun_MVN_full <- function(yex, ydep, maxit, start, nOptim){
  
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
    return(out)
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
    return(out)
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
        return(out)
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
        return(out)
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
