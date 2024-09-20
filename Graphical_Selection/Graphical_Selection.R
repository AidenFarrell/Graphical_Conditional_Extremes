################################################################################
## Adding required packages
required_pckgs <- c("glasso", "gtools", "igraph", "ppcor", "rlang")
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in required scripts
source("Miscellaneous_Functions/MVAGG_Functions.R")
source("Model_Fitting/Cond_Extremes_MVN_Residuals.R")

################################################################################

## Function to infer the graphical structure of the residuals
Infer_Adj_Matrix_glasso_residuals <- function(data, cond = 1,
                                              constrain = TRUE, q = c(0,1), v = 20, aLow = -1, 
                                              maxit = 1e+6, nOptim = 1,
                                              start_HT, start_AGG, rho = 0.1){
  
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
    start_AGG <- matrix(rep(c(0, 1.5, 2, 1.5), d-1), nrow = d-1, ncol = 4, byrow = TRUE)
  }
  if(!is.numeric(start_HT) | !is.numeric(start_AGG)){
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
  if(any(abs(start_HT) > 1)){
    stop("Initial starting values are outside the parameter space for Heffernan and Tawn model")
  }
  if(any(start_AGG[,-1] <= 0)){
    stop("Initial starting values are outside the parameter space for MVAGG")
  }
  
  ## check rho
  if(length(rho) > 1 | rho < 0){
    stop("rho must be a single positive real number")
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
      
      ## Step 4 
      ## Use glasso to infer the graph
      cor_Z <- cor(Z_Gaussian)
      
      glasso_fit <- glasso(s = cor_Z, rho = rho, nobs = n, thr = 1e-8, maxit = 1e+6, penalize.diagonal = FALSE)
      
      ## Get the adjacency matrix
      adj_matrix <- matrix(0, nrow = d-1, ncol = d-1)
      adj_matrix[which(glasso_fit$wi != 0)] <- 1
      diag(adj_matrix) <- 0
      return(adj_matrix)
    }
  }
}

## Function to infer the graphical structure from the data using a graphical lasso
Infer_Adj_Matrix_glasso_data <- function(data, u, cond = 1, rho = 0.1){
  
  ## Obtain information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  ## Check conditions on the threshold
  if(!is.numeric(u) | length(u) > 1){
    stop("u must be a single real number")
  }
  
  ## Check the conditioning random variable is valid and get dependent random variables
  if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(cond%%1 != 0 | cond <= 0 | cond > d){
    stop("cond must be a single positie integer")
  }
  
  ## check rho
  if(length(rho) > 1 | !is.numeric(rho) | rho < 0){
    stop("rho must be a single positive real number")
  }
  
  ## Use glasso to infer the graph
  data_large <- data[which(data[,cond] > u),]
  
  glasso_fit <- glasso(s = cor(data_large[-cond,-cond]), rho = rho,
                       nobs = n, thr = 1e-8, maxit = 1e+6,
                       penalize.diagonal = FALSE)
  
  ## Get the adjacency matrix
  adj_matrix <- matrix(0, nrow = d-1, ncol = d-1)
  adj_matrix[which(glasso_fit$wi != 0, arr.ind = TRUE)] <- 1
  diag(adj_matrix) <- 0
  return(adj_matrix)
}

## Function to infer the graphical structure from the data using a graphical lasso
Infer_Adj_Matrix_ppcor_data <- function(data, u, cond = 1, alpha = 0.1){
  
  ## Obtain information from the data
  dim_data <- dim(data)
  if(is.null(dim_data)){
    stop("Data must be a matrix with at least d = 2 columns")
  }
  else{
    d <- dim_data[2]
    n <- dim_data[1]
  }
  
  ## Check conditions on the threshold
  if(!is.numeric(u) | length(u) > 1){
    stop("u must be a single real number")
  }
  
  ## Check the conditioning random variable is valid and get dependent random variables
  if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(cond%%1 != 0 | cond <= 0 | cond > d){
    stop("cond must be a single positie integer")
  }
  
  ## check alpha
  if(length(alpha) > 1 | !is.numeric(alpha) | alpha <= 0 | alpha >= 1){
    stop("alpha must be a single real number in the region (0,1)")
  }
  
  ## Use partial correlaiton test to infer the graph
  data_large <- data[which(data[,cond] > u),]
  
  inferred_Edges_i <- which(pcor(data_large[-cond,-cond])$p.value < alpha/2, 
                            arr.ind = TRUE)
  
  ## Get the adjacency matrix
  adj_matrix <- matrix(0, nrow = d-1, ncol = d-1)
  adj_matrix[inferred_Edges_i] <- 1
  diag(adj_matrix) <- 0
  return(adj_matrix)
}

################################################################################

## Combine sub-matrices, obtained by conditioning on each site in turn, into a single
## weighted matrix
Combine_Sub_Adjacency_Matrices <- function(Sub_Matrices){
  if(!is.list(Sub_Matrices)){
    stop("Sub_Matrices must be a list of adjacency matrices")
  }
  if(!all(sapply(Sub_Matrices, is.square.matrix))){
    stop("One or more of the matrices provided is not square")
  }
  
  d <- length(Sub_Matrices)
  if(!all(sapply(Sub_Matrices, dim) == d-1)){
    stop("Matrices provided have different dimensions")
  }
  
  Inferred_Adjacency_Matrix <- matrix(0, nrow = d, ncol = d)
  for(i in 1:(d)){
    Inferred_Adjacency_Matrix[-i,-i] <- 
      Inferred_Adjacency_Matrix[-i,-i] + Sub_Matrices[[i]]
  }
  return(Inferred_Adjacency_Matrix)
}

## Plot the inferred graph when we combine the submatrices together
Plot_Graph_From_Adjacency_Matrix <- function(Adj_Matrix, n_sim = 1, layout, true_graph = NA){
  if(!is.square.matrix(Adj_Matrix)){
    stop("Adj_Matrix must be a square matrix")
  }
  if(!is.symmetric.matrix(Adj_Matrix)){
    stop("Adj_Matrix must be a symmetric matrix")
  }
  if(n_sim <= 0 | n_sim%%1 != 0){
    stop("Number of samples must be a postive integer")
  }
  upper_tri_elements <- upper.tri(Adj_Matrix)
  edges <- which(Adj_Matrix*upper_tri_elements > 0, arr.ind = TRUE)
  edges <- edges[order(edges[,1]),]
  g_out <- graph(t(edges), directed = FALSE)
  
  d <- ncol(Adj_Matrix)
  V(g_out)$label <- NA
  V(g_out)$size <- 15
  E(g_out)$weight <- (Adj_Matrix[edges])/n_sim
  
  if(is_missing(layout)){
    layout <- layout_with_fr(graph)
  }
  if(is_igraph(true_graph)){
    edges_df <- as.data.frame(edges)
    ## Get the true edges from the true graph and reorder them too
    true_edges <- which(as.matrix(as_adjacency_matrix(true_graph))*upper_tri_elements > 0, arr.ind = TRUE)
    true_edges <- true_edges[order(true_edges[,1]),]
    true_edges <- as.data.frame(true_edges)
    
    ## see which edges are true edges and colour these in black, otherwise colour the edges in grey
    edges_df$exists <- do.call(paste0, edges_df) %in% do.call(paste0, true_edges)
    edges_df$colour = "grey"
    edges_df$colour[which(edges_df$exists == TRUE)] <- "black"
    E(g_out)$color <- edges_df$colour
    
    ## plot the graph
    plot(g_out, layout = layout, edge.width = 5*E(g_out)$weight)
  }
  else{
    plot(g_out, layout = layout, edge.width = 5*E(g_out)$weight, edge.color = "grey") 
  }
}
