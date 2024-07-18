################################################################################
## Reading in required scripts
source("Cond_Extremes_MVN_Residuals.R")
source("MVAGG_Functions.R")

################################################################################
Cond_Extremes_MVAGG_Two_Step <- function(data, cond, graph = NA, 
                                         start_HT = c(0.1, 0.1), start_AGG = c(0, 1, 2, 1.5),
                                         constrain = TRUE, q = c(0,1), v = 10, aLow = -1,
                                         maxit = 1e+6, nOptim = 1){
  
  ## Obtain information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  ## Check the conditioning random variable is valid and get dependent random variables
  if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(cond%%1 != 0 | cond <= 0 | cond > d){
    stop("cond must be a single positie integer")
  }
  dependent <- (1:d)[-cond]
  
  ## Obtain the starting parameters
  if(missing(start_HT)){
    start_HT <- matrix(0.1, nrow = d-1, ncol = 2)
  }
  if(missing(start_AGG)){
    start_AGG <- matrix(rep(c(0, 1, 2, 1.5), d-1), nrow = d-1, ncol = 4, byrow = TRUE)
  }
  else if(!is.numeric(start_HT) | !is.numeric(start_AGG)){
    stop("start_HT and start_AGG must be vectors")
  }
  if(length(start_HT) == 2){
    start_HT <- matrix(rep(start_HT, d-1), ncol = 2, byrow = TRUE)
  }
  if(length(start_AGG) == 4){
    start_AGG <- matrix(rep(start_AGG, d-1), ncol = 4, byrow = TRUE)
  }
  if(length(start_HT) != 2*(d-1) | length(start_AGG) != 4*(d-1)){
    stop("start_HT and start_AGG must be vectors of length 2 or 2(d-1) and 4 or 4(d-1), respectively")
  }
  else if(any(abs(start_HT) > 1)){
    stop("Initial starting values are outside the parameter space for Heffernan and Tawn model")
  }
  else if(any(start_AGG[,-1] <= 0)){
    stop("Initial starting values are outside the parameter space for MVAGG")
  }
  
  ## Separate data into the conditioning and unconditioned random variables
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  ## Step 1
  ## Fit the original conditional multivaraite extremes model to the data to
  ## obtain the fitted residuals and dependence parameters
  res_HT <- lapply(1:(d-1), function(i){
    qfun_MVN_indep(yex = yex, ydep = as.matrix(ydep[,i]),
                   constrain = constrain, aLow = aLow, q = q, v = v,
                   maxit = maxit, start = start_HT[i,], nOptim = nOptim)})
  
  ## Get the necessary output
  z <- as.matrix(sapply(res_HT, function(x){x$Z}))
  a_hat <- sapply(res_HT, function(x){x$par$a})
  b_hat <- sapply(res_HT, function(x){x$par$b})
  
  ## Update the starting parameter for the location parameter as this can be difficult
  ## to initiate
  start_AGG[,1] <- apply(z, 2, mean)
  
  ## If Heffernan and Tawn model has not fit, break the function now
  if(any(is.na(a_hat))){
    warning("Error in optim call from Cond_extremes_MVN")
    out <- list()
    out$par <- list(a = rep(NA, d), b = rep(NA, d), loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = as(matrix(NA, ncol = d, nrow = d), "sparseMatrix"))
    out$Z <- matrix(NA, nrow = n, ncol = d)
    out$value <- NA
    out$convergence <- NA
  }
  else{
    ## Step 2
    ## Model the fitted residuals using a MVAGG distribution
    
    ## Determine the components in the model
    if(!is_igraph(graph)){
      warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID Guassian random variables.")
      
      ## Fit the independence model
      res <- lapply(1:(d-1), function(i){
        qfun_MVAGG_Two_Step_Indep(z = as.matrix(z[,i]), maxit = maxit, start = c(start_AGG[i,]))})
      
      ## Extract the output
      out <- list()
      
      out$par$main <- matrix(data = c(a_hat, b_hat,
                                      do.call(c, lapply(res, function(x){x$par$loc})),
                                      do.call(c, lapply(res, function(x){x$par$scale_1})),
                                      do.call(c, lapply(res, function(x){x$par$scale_2})),
                                      do.call(c, lapply(res, function(x){x$par$shape}))),
                             nrow = 6, ncol = d-1, byrow = TRUE)
      out$Z <- z
      if(any(is.na(out$par$main))){
        out$par$main <- matrix(NA, nrow = 6, ncol = d-1)
        out$Z <- matrix(NA, nrow = nrow(out$Z), ncol = ncol(out$Z)) 
      }
      rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
      colnames(out$par$main) <- colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
      
      if(any(is.na(out$par$main))){
        out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = rep(NA, d-1)) 
      }
      else{
        out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = rep(1, d-1))
      }
      
      out$loglike <- sum(sapply(res, function(x){-x$value}))
      out$convergence <- max(sapply(res, function(x){x$convergence}))
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
        res <- qfun_MVAGG_Two_Step_Full(z = z, start = c(start_AGG), maxit = maxit)
        
        ## Extract the output
        out <- list()
        
        out$par$main <- matrix(data = c(a_hat, b_hat, res$par$loc,
                                        res$par$scale_1, res$par$scale_2, res$par$shape),
                               nrow = 6, ncol = d-1, byrow = TRUE)
        out$Z <- z
        if(any(is.na(out$par$main))){
          out$par$main <- matrix(NA, nrow = 6, ncol = d-1)
          out$Z <- matrix(NA, nrow = nrow(out$Z), ncol = ncol(out$Z)) 
        }
        rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
        colnames(out$par$main) <- colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
        
        out$par$Gamma <- as(res$par$Gamma, "sparseMatrix")
        
        out$loglike <- -res$value
        out$convergence <- res$convergence
      }
      else{
        if(is_connected(graph_cond)){
          ## If the graph is connected use the graphical model
          ## determine the edges in the conditional graph
          all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, g_edges)
          non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
          
          ## Fit the graphical model
          res <- qfun_MVAGG_Two_Step_Graph(z = z, Gamma_zero = non_edges,
                                           maxit = maxit, start = c(start_AGG))
          
          ## Extract the output
          out <- list()
          
          out$par$main <- matrix(data = c(a_hat, b_hat, res$par$loc,
                                          res$par$scale_1, res$par$scale_2, res$par$shape),
                                 nrow = 6, ncol = d-1, byrow = TRUE)
          out$Z <- z
          if(any(is.na(out$par$main))){
            out$par$main <- matrix(NA, nrow = 6, ncol = d-1)
            out$Z <- matrix(NA, nrow = nrow(out$Z), ncol = ncol(out$Z)) 
          }
          rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
          colnames(out$par$main) <- colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
          
          out$par$Gamma <- as(res$par$Gamma, "sparseMatrix")
          
          out$loglike <- -res$value
          out$convergence <- res$convergence
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
              res_1[[i]] <- qfun_MVAGG_Two_Step_Indep(z = as.matrix(z[,v_comps[[i]]]), 
                                                      start = c(start_AGG[v_comps[[i]],]),
                                                      maxit = maxit)
            }
            else{
              all_edges_comp <- as.data.frame(combinations(n = n_vertices, r = 2, v = 1:n_vertices))
              g_edges_comp <- as.data.frame(as_edgelist(g_comps[[i]]))
              ## Fit saturated model to this component
              if(nrow(all_edges_comp) == nrow(g_edges_comp)){
                res_1[[i]] <- qfun_MVAGG_Two_Step_Full(z = as.matrix(z[,v_comps[[i]]]),
                                                       start = c(start_AGG[v_comps[[i]],]),
                                                       maxit = maxit)
              }
              else{
                ## Fit graphical model to this component
                all_edges_comp$exists <- do.call(paste0, all_edges_comp) %in% do.call(paste0, g_edges_comp)
                non_edges_comp <- as.matrix(all_edges_comp[which(all_edges_comp$exists == FALSE), 1:2])
                res_1[[i]] <- qfun_MVAGG_Two_Step_Graph(z = as.matrix(z[,v_comps[[i]]]), 
                                                        Gamma_zero = non_edges_comp,
                                                        start = c(start_AGG[v_comps[[i]],]),
                                                        maxit = maxit)
              } 
            }
          }
          
          ## Extract the output
          out <- list()
          
          out$par$main <- matrix(data = c(a_hat, b_hat,
                                          do.call(c, lapply(res_1, function(x){unname(x$par$loc)})),
                                          do.call(c, lapply(res_1, function(x){unname(x$par$scale_1)})),
                                          do.call(c, lapply(res_1, function(x){unname(x$par$scale_2)})),
                                          do.call(c, lapply(res_1, function(x){unname(x$par$shape)}))),
                                 nrow = 6, ncol = d-1, byrow = TRUE)
          out$Z <- z
          if(any(is.na(out$par$main))){
            out$par$main <- matrix(NA, nrow = 6, ncol = d-1)
            out$Z <- matrix(NA, nrow = nrow(out$Z), ncol = ncol(out$Z)) 
          }
          rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
          colnames(out$par$main) <- colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
          
          Gamma_intrim <- matrix(0, nrow = d-1, ncol = d-1)
          Gamma_entries <- lapply(v_comps, function(x){permutations(n = length(x), r = 2, v = x, repeats.allowed = TRUE)})
          for(i in 1:n_comps){
            if(nrow(Gamma_entries[[i]]) == 1){
              Gamma_intrim[Gamma_entries[[i]]] <- 1
            }
            else{
              Gamma_intrim[Gamma_entries[[i]]] <- res_1[[i]]$par$Gamma 
            }
          }
          if(any(is.na(Gamma_intrim))){
            Gamma_intrim <- matrix(NA, nrow = d-1, ncol = d-1)
          }
          out$par$Gamma <- as(Gamma_intrim, "sparseMatrix")
          
          out$loglike <- sum(sapply(res_1, function(x){-x$value}))
          out$convergence <- max(sapply(res_1, function(x){x$convergence}))
        }
      }
    }
  }
  class(out) <- "Cond_Extremes_MVAGG"
  return(out)
}

qfun_MVAGG_Two_Step_Indep <- function(z, maxit, start){
  
  ## Function for the optim to optimise over
  Qpos <- function(param, z, negative = FALSE){
    ## Extract starting parameters
    sigma_1 <- param[2]
    sigma_2 <- param[3]
    shape <- param[4]
    if(any(sigma_1 <= 0) | any(sigma_2 <= 0) | any(shape <= 0)){
      return((-10^10)*(-1)^negative)
    }
    else{
      ## Fit the model
      res <- sum(dagg(x = c(z), loc = param[1], scale_1 = sigma_1, scale_2 = sigma_2, shape = shape, log = TRUE))
      
      ## Check the value is valid
      if(is.infinite(res)){
        res <- return((-10^10)*(-1)^negative)
        warning("Infinite value of Q in mexDependence")
      }
      return(res*(-1)^negative)
    }
  }
  
  ## Fit the model
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, negative = TRUE, method = "Nelder-Mead", hessian = FALSE),
             silent = FALSE)
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = NA, scale_1 = NA, scale_2 = NA, shape = NA)
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = NA, scale_1 = NA, scale_2 = NA, shape = NA)
    out$value <- NA
    out$convergence <- NA
  }
  else if(!is.na(fit$par[1])){
    ## Extract the MLEs 
    out <- list()
    out$par <- list(loc = fit$par[1], scale_1 = fit$par[2], scale_2 = fit$par[3], shape = fit$par[4])
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    warning("Unknown error in Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = NA, scale_1 = NA, scale_2 = NA, shape = NA)
    out$value <- NA
    out$convergence <- NA
  }
  return(out)
}

qfun_MVAGG_Two_Step_Graph <- function(z, Gamma_zero, maxit, start){
  
  ## Function for the optim to optimise over
  Qpos <- function(param, z, Gamma_zero, negative = FALSE){
    ## Extract the the starting parameters
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
        
        Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = FALSE, thr = 1e-9, zero = Gamma_zero))
        chol_Gamma <- chol(Lasso_est$wi)
        y <- chol_Gamma%*%t(Q_F_z)
        l_mvnorm <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
        
        res <- l_mvnorm + sum(l_f_z) - sum(l_dnorm)
        
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
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, Gamma_zero = Gamma_zero, negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  
  ## Extract the output
  d <- ncol(z)
  
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$value <- NA
    out$convergence <- NA
  }
  else if(!is.na(fit$par[1])){
    ## Extract MLEs from the model
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    sigma_1_hat <- fit$par[(d + 1):(2*d)]
    sigma_2_hat <- fit$par[(2*d + 1):(3*d)]
    shape_hat <- fit$par[(3*d + 1):(4*d)]
    
    ## Extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(pagg(q = c(z[,i]), loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, thr = 1e-9, zero = Gamma_zero))
    Gamma_hat <- Lasso_est$wi
    
    ## Organise the output
    out <- list()
    out$par <- list(loc = mu_hat, scale_1 = sigma_1_hat, scale_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    warning("Unknwon error in Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$value <- NA
    out$convergence <- NA
  }
  return(out)
}

qfun_MVAGG_Two_Step_Full <- function(z, maxit, start){
  
  ## Function for the optim to optimise over
  Qpos <- function(param, z, negative = FALSE){
    ## Extract the the starting parameters
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
        
        # Lasso_est <- suppressWarnings(glasso(s = Sigma_start, rho = 0, penalize.diagonal = FALSE, thr = 1e-9))
        # chol_Gamma <- chol(Lasso_est$wi)
        chol_Gamma <- chol(solve(Sigma_start))
        y <- chol_Gamma%*%t(Q_F_z)
        l_mvnorm <- -n*d*log(2*pi)/2 - n*sum(log(1/diag(chol_Gamma))) - sum(y^2)/2
        
        res <- l_mvnorm + sum(l_f_z) - sum(l_dnorm)
        
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
  fit <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit),
                   z = z, negative = TRUE, method = "BFGS", hessian = FALSE),
             silent = FALSE)
  
  ## Extract the output
  d <- ncol(z)
  
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$value <- NA
    out$convergence <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence in Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$value <- NA
    out$convergence <- NA
  }
  else if(!is.na(fit$par[1])){
    #Extract MLEs of alpha and beta
    d <- ncol(z)
    mu_hat <- fit$par[1:d]
    sigma_1_hat <- fit$par[(d + 1):(2*d)]
    sigma_2_hat <- fit$par[(2*d + 1):(3*d)]
    shape_hat <- fit$par[(3*d + 1):(4*d)]
    
    ## extract MLE of Gamma
    Q_F_z <- as.matrix(sapply(1:d, function(i){qnorm(pagg(q = c(z[,i]), loc = mu_hat[i], scale_1 = sigma_1_hat[i], scale_2 = sigma_2_hat[i], shape = shape_hat[i]))}))
    # Lasso_est <- suppressWarnings(glasso(s = cor(Q_F_z), rho = 0, penalize.diagonal = FALSE, thr = 1e-9))
    # Gamma_hat <- Lasso_est$wi
    Gamma_hat <- solve(cor(Q_F_z))
    
    out <- list()
    out$par <- list(loc = mu_hat, scale_1 = sigma_1_hat, scale_2 = sigma_2_hat, shape = shape_hat, Gamma = Gamma_hat)
    out$value <- fit$value
    out$convergence <- fit$convergence
  }
  else{
    warning("Unknown error in Cond_Extremes_MVAGG")
    out <- list()
    out$par <- list(loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = matrix(NA, ncol = d, nrow = d))
    out$value <- NA
    out$convergence <- NA
  }
  return(out)
}