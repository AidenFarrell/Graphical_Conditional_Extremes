################################################################################
## Reading in required scripts
source("~/PhD/Project_2/Graphical_Conditional_Extremes/MVAGG_Functions.R")

################################################################################

## Conditional extremes model under the one-step approach assuming the residuals
## follow a mutltivariate asymmetric generalised Gaussian distribution
Cond_extremes_MVAGGD <- function(data, cond, graph = NA, start,
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
  
  #determine the components in the model
  if(!is_igraph(graph)){
    warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
    
    ## Fit the model
    res <- lapply(1:(d-1), function(i){
      qfun_MVAGGD_indep(yex = yex, ydep = as.matrix(ydep[,i]),
                        maxit = maxit, start = start[i,])})
    
    ## 
    out <- list()
    
    out$par$main <- matrix(data = c(do.call(c, lapply(res, function(x){x$par$a})),
                                    do.call(c, lapply(res, function(x){x$par$b})),
                                    do.call(c, lapply(res, function(x){x$par$loc})),
                                    do.call(c, lapply(res, function(x){x$par$scale_1})),
                                    do.call(c, lapply(res, function(x){x$par$scale_2})),
                                    do.call(c, lapply(res, function(x){x$par$shape}))),
                           nrow = 6, ncol = d-1, byrow = TRUE)
    rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
    colnames(out$par$main) <- sapply(dependent, function(x){paste0("Column", x)})
    
    out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = rep(1, d-1))
    
    out$loglike <- sum(sapply(res, function(x){-x$value}))
    out$convergence <- max(sapply(res, function(x){x$convergence}))
    out$Z <- do.call(cbind, lapply(res, function(x){x$Z}))
    colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
    
    class(out) <- "Cond_extremes_MVAGG"
    return(out)
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

qfun_MVAGGD_indep <- function(yex, ydep, maxit, start, nOptim){
  #function for the optim to optimise over
  Qpos <- function(param, yex, ydep, negative = FALSE){
    #get the starting parameters
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
      ## get the Z data
      b_yi <- yex^b
      z <- c((ydep - yex*a)/b_yi)
      res <- sum(dagg(x = z, loc = mu, scale_1 = sigma_1, scale_2 = sigma_2, shape = shape, log = TRUE)) - sum(log(b_yi))
      
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
    warning("Error in optim call from Cond_extremes_MVAGGD")
    out <- list()
    out$par <- list(a = NA, b = NA, loc = NA, scale_1 = NA, scale_2 = NA, shape = NA)
    out$Z <- NA 
    out$value <- NA
    out$convergence <- NA
    return(out)
  }
  else if(fit$convergence != 0){
    warning("Non-convergence in Cond_extremes_MVAGGD")
    out <- list()
    out$par <- list(a = NA, b = NA, loc = NA, scale_1 = NA, scale_2 = NA, shape = NA)
    out$Z <- NA 
    out$value <- NA
    out$convergence <- NA
    return(out)
  }
  else if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    a_hat <- fit$par[1]
    b_hat <- fit$par[2]
    
    #obtain the residuals
    Z <- (ydep - yex*a_hat)/(yex^b_hat)
    
    out <- list()
    out$par <- list(a = a_hat, b = b_hat, loc = fit$par[3], scale_1 = fit$par[4], scale_2 = fit$par[5], shape = fit$par[6])
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
