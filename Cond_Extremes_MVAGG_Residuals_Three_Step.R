################################################################################
## Reading in required scripts
source("MVAGG_Functions.R")
source("Cond_Extremes_MVN_Residuals.R")

################################################################################

Cond_Extremes_MVAGG_Three_Step <- function(data, cond = 1, graph = NA,
                                            constrain = TRUE, q = c(0,1), v = 10, aLow = -1, 
                                            maxit = 1e+6, nOptim = 1,
                                            start_HT, start_AGG){
  
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
  
  ## List for output to be saved
  out <- list()
  
  ## Step 1
  ## Fit the original conditional multivaraite extreme value model to the data to
  ## obtain the fitted residuals and dependence parameters
  yex <- as.matrix(data[,cond])
  ydep <- as.matrix(data[,-cond])
  
  res_HT <- lapply(1:(d-1), function(i){
    qfun_MVN_indep(yex = yex, ydep = as.matrix(ydep[,i]),
                   constrain = constrain, aLow = aLow, q = q, v = v,
                   maxit = maxit, start = start_HT[i,], nOptim = nOptim)})
  
  ## Get the necessary output
  z <- as.matrix(sapply(res_HT, function(x){x$Z}))
  a_hat <- sapply(res_HT, function(x){x$par$a})
  b_hat <- sapply(res_HT, function(x){x$par$b})
  
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
    ## Marginally fit an asymmetric generalised Gaussian distribution to the fitted residuals
    
    ## Update the starting parameter for the location parameter as this can be difficult
    ## to initiate
    start_AGG[,1] <- apply(z, 2, mean)
    res_AGG <- lapply(1:(d-1), function(i){fit_agg(par = start_AGG[i,], data = z[,i])})
    
    ## Get the necessary output
    loc_hat <- sapply(res_AGG, function(x){x$par[1]})
    scale_1_hat <- sapply(res_AGG, function(x){x$par[2]})
    scale_2_hat <- sapply(res_AGG, function(x){x$par[3]})
    shape_hat <- sapply(res_AGG, function(x){x$par[4]})
    
    ## If AGG has not fit to the marginal distributions, break this now
    if(any(is.na(loc_hat))){
      warning("Error in optim call from AGG")
      out <- list()
      out$par <- list(a = rep(NA, d), b = rep(NA, d), loc = rep(NA, d), scale_1 = rep(NA, d), scale_2 = rep(NA, d), shape = rep(NA, d), Gamma = as(matrix(NA, ncol = d, nrow = d), "sparseMatrix"))
      out$Z <- matrix(NA, nrow = n, ncol = d)
      out$value <- NA
      out$convergence <- NA
    }
    else{
      ## Step 3
      ## Determine the covariance matrix from the data
      
      ## First transform the data onto standard Gaussian margins
      Z_Gaussian <- sapply(1:(d-1), function(i){
        qnorm(pagg(q = z[,i],
                   loc = loc_hat[i],
                   scale_1 = scale_1_hat[i],
                   scale_2 = scale_2_hat[i],
                   shape = shape_hat[i]))})
      
      if(!is_igraph(graph)){
        warning("\nNo graphical structure has been provided.\n \nWe assume the residuals are IID AGG random variables.")
        out$par$Gamma <- sparseMatrix(i = 1:(d-1), j = 1:(d-1), x = rep(1, d-1))
      }
      else{
        ## Graph is provided so we need to determine its structure
        graph_cond <- delete_vertices(graph, cond)
        
        ## check if the graph is full to determine the model
        all_edges <- as.data.frame(combinations(n = d-1, r = 2, v = 1:(d-1)))
        g_edges <- as.data.frame(as_edgelist(graph_cond))
        all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, g_edges)
        non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
        
        if(is_empty(non_edges)){
          out$par$Gamma <- as(solve(cor(Z_Gaussian)), "sparseMatrix")
        }
        else{
          ## Use a graphical model
          ## Determine the missing edges in the graph
          fit_glasso <- suppressWarnings(glasso(s = cor(Z_Gaussian), rho = 0, nobs = n,
                                                zero = non_edges, thr = 1e-8, maxit = 1e+6, penalize.diagonal = FALSE)) 
          out$par$Gamma <- as(fit_glasso$wi, "sparseMatrix")
        }
      }
      
      ## Extract the output
      out$par$main <- as.matrix(rbind(a_hat, b_hat, loc_hat, scale_1_hat, scale_2_hat, shape_hat))
      out$Z <- z
      if(any(is.na(out$par$main))){
        out$par$main <- matrix(NA, nrow = 6, ncol = d-1)
        out$Z <- matrix(NA, nrow = nrow(out$Z), ncol = ncol(out$Z)) 
      }
      rownames(out$par$main) <- c("a", "b", "loc", "scale_1", "scale_2", "shape")
      colnames(out$par$main) <- colnames(out$Z) <- sapply(dependent, function(x){paste0("Column", x)})
      
      out$convergence <- max(sapply(res_HT, function(x){x$convergence}),
                             sapply(res_AGG, function(x){x$convergence}))  
      class(out) <- "Cond_Extremes_MVAGG"
      return(out)
    }
  }
}