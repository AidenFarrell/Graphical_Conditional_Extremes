################################################################################
## Reading in required scripts
source("~/PhD/Project_2/Graphical_Conditional_Extremes/MVAGG_Functions.R")

################################################################################

## Conditional extremes model under the one-step approach assuming the residuals
## follow a multivariate asymmetric generalised Gaussian distribution
Cond_extremes_MVAGG <- function(data, cond, graph = NA, start,
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
  
  ## check the conditioning random variable is valid and get dependent random variables
  if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(cond%%1 != 0 | cond <= 0 | cond > d){
    stop("cond must be a single positie integer")
  }
  dependent <- (1:d)[-cond]
  
  ## get the starting parameters
  if(missing(start)){
    start <- c(0.1, 0.1, 0, 1.5, 2, 1.5)
  }
  else if(!is.numeric(start)){
    stop("start must be a vector")
  }
  else if(length(start) != 6*(d-1)){
    stop("start must be a vector of length 6(d-1)")
  }
  else if(any(abs(start[,1:2]) > 1)){
    stop("Initial starting values are outside the parameter sapce")
  }
  if(length(start) == 6){
    start <- matrix(rep(start, d-1), ncol = 6, byrow = TRUE)
  }
  
  #separate data into the conditioning and unconditioned random variables
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  ## Determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## Fit the independence model
    res <- lapply(1:(d-1), function(i){
      qfun_MVAGG_indep(yex = yex, ydep = as.matrix(ydep[,i]),
                        maxit = maxit, start = start[i,])})
    
    ## Extract the output
    out <- list()
    
    out$par$main <- matrix(data = c(do.call(c, lapply(res, function(x){x$par$a})),
                                    do.call(c, lapply(res, function(x){x$par$b})),
                                    do.call(c, lapply(res, function(x){x$par$loc})),
                                    do.call(c, lapply(res, function(x){x$par$scale_1})),
                                    do.call(c, lapply(res, function(x){x$par$scale_2})),
                                    do.call(c, lapply(res, function(x){x$par$shape}))),
                           nrow = 6, ncol = d-1, byrow = TRUE)
    if(any(is.na(out$par$main))){
      out$par$main <- matrix(NA, nrow = 6, ncol = d-1)
    }
    rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
    colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
    
    if(any(is.na(out$par$main))){
      out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = rep(NA, d-1)) 
    }
    else{
      out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = rep(1, d-1))
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
    
    ## check if the graph is full, if so fit then model instantly
    all_edges <- as.data.frame(combinations(n = d-1, r = 2, v = 1:(d-1)))
    g_edges <- as.data.frame(as_edgelist(graph_cond))
    n_edges_full <- nrow(all_edges)
    n_edges_graph_cond <- nrow(g_edges)
    
    if(n_edges_full == n_edges_graph_cond){
      ## Fit the saturated model
      res <- qfun_MVAGG_full(yex = yex, ydep = ydep, maxit = maxit, start = c(start))
      
      ## Extract the output
      out <- list()
      
      out$par$main <- matrix(data = c(res$par$a, res$par$b, res$par$loc,
                                      res$par$scale_1, res$par$scale_2, res$par$shape),
                             nrow = 6, ncol = d-1, byrow = TRUE)
      rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
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
        res <- qfun_MVAGG_graph(yex = yex, ydep = ydep, Gamma_zero = non_edges,
                                maxit = maxit, start = c(start))
        
        ## Extract the output
        out <- list()
        
        out$par$main <- matrix(data = c(res$par$a, res$par$b, res$par$loc,
                                        res$par$scale_1, res$par$scale_2, res$par$shape),
                               nrow = 6, ncol = d-1, byrow = TRUE)
        rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
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
            ## Fit independence model to this component
            res_1[[i]] <- qfun_MVAGG_indep(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]), 
                                           maxit = maxit, start = c(start[v_comps[[i]],]))
          }
          else{
            all_edges_comp <- as.data.frame(combinations(n = n_vertices, r = 2, v = 1:n_vertices))
            g_edges_comp <- as.data.frame(as_edgelist(g_comps[[i]]))
            ## Fit saturated model to this component
            if(nrow(all_edges_comp) == nrow(g_edges_comp)){
              res_1[[i]] <- qfun_MVAGG_full(yex = yex, ydep = as.matrix(ydep[,v_comps[[i]]]),
                                            maxit = maxit, start = c(start[v_comps[[i]],]))
            }
            else{
              ## Fit graphical model to this component
              all_edges_comp$exists <- do.call(paste0, all_edges_comp) %in% do.call(paste0, g_edges_comp)
              non_edges_comp <- as.matrix(all_edges_comp[which(all_edges_comp$exists == FALSE), 1:2])
              res_1[[i]] <- qfun_MVAGG_graph(yex = yex, ydep = ydep[,v_comps[[i]]], 
                                           Gamma_zero = non_edges_comp,
                                           maxit = maxit, 
                                           start = c(start[v_comps[[i]],]))
            } 
          }
        }
        
        ## Extract the output
        out <- list()
        
        out$par$main <- matrix(data = c(do.call(c, lapply(res_1, function(x){unname(x$par$a)})), 
                                        do.call(c, lapply(res_1, function(x){unname(x$par$b)})),
                                        do.call(c, lapply(res_1, function(x){unname(x$par$loc)})),
                                        do.call(c, lapply(res_1, function(x){unname(x$par$scale_1)})),
                                        do.call(c, lapply(res_1, function(x){unname(x$par$scale_2)})),
                                        do.call(c, lapply(res_1, function(x){unname(x$par$shape)}))),
                               nrow = 6, ncol = d-1, byrow = TRUE)
        if(any(is.na(out$par$main))){
          out$par$main <- matrix(NA, nrow = 6, ncol = d-1)
        }
        rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
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
  class(out) <- "Cond_extremes_MVAGG"
  return(out)
}

qfun_MVAGG_indep <- function(yex, ydep, maxit, start, nOptim){
  ## Function for the optim to optimise over
  Qpos <- function(param, yex, ydep, negative = FALSE){
    ## Extract starting parameters
    a <- param[1]
    b <- param[2]
    mu <- param[3]
    sigma_1 <- param[4]
    sigma_2 <- param[5]
    shape <- param[6]
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0) | any(abs(a) > 1) | any(b >= 1)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## Obtain the residuals
      b_yi <- yex^b
      z <- c((ydep - yex*a)/b_yi)
      res <- sum(dagg(x = z, loc = mu, scale_1 = sigma_1, scale_2 = sigma_2, shape = shape, log = TRUE)) - sum(log(b_yi))
      
      ## Check the value is valid
      if(is.infinite(res)){
        res <- return((-10^10)*(-1)^negative)
        warning("Infinite value of Q in mexDependence")
      }
      return(res*(-1)^negative)
    }
  }
  
  ## Fit the model
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                   control = list(maxit = maxit),
                   negative = TRUE, method = "Nelder-Mead", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = NA, b = NA, loc = NA, scale_1 = NA, scale_2 = NA, shape = NA)
    out$Z <- NA 
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = NA, b = NA, loc = NA, scale_1 = NA, scale_2 = NA, shape = NA)
    out$Z <- NA 
    out$value <- NA
    out$convergence <- NA
  }
  else if(!is.na(fit$par[1])){
    ## Extract MLEs of alpha and beta
    a_hat <- fit$par[1]
    b_hat <- fit$par[2]
    
    ## Obtain the residuals
    Z <- (ydep - yex*a_hat)/(yex^b_hat)
    
    ## Organise the output
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, loc = fit$par[3], scale_1 = fit$par[4], scale_2 = fit$par[5], shape = fit$par[6])
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    warning("Unknown error in Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = NA, b = NA, loc = NA, scale_1 = NA, scale_2 = NA, shape = NA)
    out$Z <- NA 
    out$value <- NA
    out$convergence <- NA
  }
  return(out)
}

qfun_MVAGG_graph <- function(yex, ydep, Gamma_zero, maxit, start, nOptim){
  
  ## Function for the optim to optimise over
  Qpos <- function(param, yex, ydep, Gamma_zero, negative = FALSE){
    ## Extract the the starting parameters
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
      ## Obtain the residuals
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      ## Transform them onto standard Gaussian margins
      Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(pagg(q = c(z[,i]), loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i]))}))
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
        l_f_z <- sapply(1:d, function(i){dagg(x = c(z[,i]), loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i], log = TRUE)})
        l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
        
        Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = FALSE, zero = Gamma_zero, thr = 1e-9))
        chol_Gamma <- chol(Lasso_est$wi)
        y <- chol_Gamma%*%t(Q_F_z)
        l_mvnorm <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
        
        res <- l_mvnorm + sum(l_f_z) - sum(l_dnorm) - sum(log(b_yi))
        
        ## Check the value is valid
        if(is.infinite(res)){
          res <- return((-10^10)*(-1)^negative)
          warning("Infinite value of Q in mexDependence")
        }
        return(res*(-1)^negative)
      }
    }
  }
  
  ## Fit the model
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                   control = list(maxit = maxit), Gamma_zero = Gamma_zero,
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  
  ## Extract the output
  d <- ncol(ydep)
  n <- nrow(ydep)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, nrow = d, ncol = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else if(!is.na(fit$par[1])){
    ## Extract MLEs from the model
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[(d + 1):(2*d)]
    mu_hat <- fit$par[(2*d + 1):(3*d)]
    sigma_1_hat <- fit$par[(3*d + 1):(4*d)]
    sigma_2_hat <- fit$par[(4*d + 1):(5*d)]
    shape_hat <- fit$par[(5*d + 1):(6*d)]
    
    ## Obtain the residuals
    a_yi_hat <- sapply(a_hat, function(a){yex*a})
    b_yi_hat <- sapply(b_hat, function(b){yex^b})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## Extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(pagg(q = c(Z[,i]), loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, zero = Gamma_zero, thr = 1e-9))
    Gamma_hat <- Lasso_est$wi
    
    ## Organise the output
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, loc = mu_hat, scale_1 = sigma_1_hat, scale_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    warning("Unknown error in Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, nrow = d, ncol = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  return(out)
}

qfun_MVAGG_full <- function(yex, ydep, maxit, start, nOptim){
  
  ## Function for the optim to optimise over
  Qpos <- function(param, yex, ydep, negative = FALSE){
    ## Extract the the starting parameters
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
      ## Obtain the residuals
      a_yi <- sapply(a, function(a){a*yex})
      b_yi <- sapply(b, function(b){yex^b})
      z <- (ydep - a_yi)/b_yi
      
      ## Transform them onto standard Gaussian margins
      Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(pagg(q = c(z[,i]), loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i]))}))
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
        l_f_z <- sapply(1:d, function(i){dagg(x = c(z[,i]), loc = mu[i], scale_1 = sigma_1[i], scale_2 = sigma_2[i], shape = shape[i], log = TRUE)})
        l_dnorm <- sapply(1:d, function(i){dnorm(Q_F_z[,i], log = TRUE)})
        
        chol_Gamma <- chol(solve(Sigma_start))
        y <- chol_Gamma%*%t(Q_F_z)
        l_mvnorm <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
        
        res <- l_mvnorm + sum(l_f_z) - sum(l_dnorm) - sum(log(b_yi))
        
        ## Check the value is valid
        if(is.infinite(res)){
          res <- return((-10^10)*(-1)^negative)
          warning("Infinite value of Q in mexDependence")
        }
        return(res*(-1)^negative)
      }
    }
  }
  
  ## Fit the model
  fit <- try(optim(par = start, fn = Qpos, yex = yex, ydep = ydep,
                   control = list(maxit = maxit),
                   negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  
  ## Extract the output
  d <- ncol(ydep)
  n <- length(yex)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, nrow = d, ncol = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, nrow = d, ncol = d))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else if(!is.na(fit$par[1])){
    ## Extract MLEs from the model
    a_hat <- fit$par[1:d]
    b_hat <- fit$par[(d + 1):(2*d)]
    mu_hat <- fit$par[(2*d + 1):(3*d)]
    sigma_1_hat <- fit$par[(3*d + 1):(4*d)]
    sigma_2_hat <- fit$par[(4*d + 1):(5*d)]
    shape_hat <- fit$par[(5*d + 1):(6*d)]
    
    ## Obtain the residuals
    a_yi_hat <- sapply(a_hat, function(a){yex*a})
    b_yi_hat <- sapply(b_hat, function(b){yex^b})
    Z <- (ydep - a_yi_hat)/b_yi_hat
    
    ## Extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(pagg(q = c(Z[,i]), loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    Gamma_hat <- solve(cor(Q_F_z))
    
    ## Organise the output
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, loc = mu_hat, scale_1 = sigma_1_hat, scale_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$Z <- Z
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    warning("Unknown error in Cond_extremes_MVAGG")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$Z <- matrix(NA, ncol = d, nrow = d)
    out$value <- NA
    out$convergence <- NA
  }
  return(out)
}
