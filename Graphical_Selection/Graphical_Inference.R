################################################################################
## Adding required packages
required_pckgs <- c("igraph")
t(t(sapply(required_pckgs, require, character.only = TRUE)))
################################################################################
## Set working directory
setwd("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/Graphical_Conditional_Extremes/")

## Reading in required functions
source("Cond_Extremes_MVAGG_Residuals.R")

################################################################################

## Function to infer the graphical structure of the residuals
Infer_Adj_Matrix <- function(data, cond = 1, rho = 0.1, par_start_HT, par_start_AGG){
  
  ## Get information about the data
  dim_data <- dim(data)
  n <- dim_data[1]
  d <- dim_data[2]
  
  ## Step 1
  ## Fit the original Heffernan and Tawn model
  if(is_missing(par_start_HT)){
    par_start_HT <- c(0.1, 0.1)
  }
  else if(length(par_start_HT) != 2){
    stop("par_start_HT must be a vector of length 2")
  }
  else if(any(abs(par_start_HT) > 1)){
    stop("Invalid HT starting parameters")
  }
  fit_HT <- Cond_extremes_graph_new(data = data, cond = cond, graph = NA,
                                    start = par_start_HT, v = ceiling(max(data)) + 1)
  
  ## Step 2
  ## Marginally fit an asymmetric generalised Gaussian distribution to the fitted
  ## residuals
  if(is_missing(par_start_AGG)){
    par_start_AGG <- c(0, 1, 2, 1.5)
  }
  else if(length(par_start_AGG) != 4){
    stop("par_start_HT must be a vector of length 2")
  }
  else if(any(par_start_AGG[-1] <= 0)){
    stop("Invalid AGG starting parameters")
  }
  fits_resid <- apply(fit_HT$Z, 2, fit_aggd, par = par_start_AGG)
  
  ## Step 3
  ## Transform the fitted residuals onto standard Gaussian margins
  Z_Gaussian <- matrix(data = NA, nrow = n, ncol = d-1)
  for(i in 1:(d-1)){
    Z_Gaussian[,i] <- qnorm(paggd(q = fit_HT$Z[,i],
                                  loc = fits_resid[[i]]$par$main[1],
                                  scale_1 = fits_resid[[i]]$par$main[2],
                                  scale_2 = fits_resid[[i]]$par$main[3],
                                  shape = fits_resid[[i]]$par$main[4]))
  }
  
  ## Step 4
  ## Fit the graphical lasso
  fit_glasso <- glasso(cor(Z_Gaussian), rho = rho, nobs = n,
                       thr = 1e-8, penalize.diagonal = FALSE)
  
  ## Step 5
  ## Obtain the inferred sub-graph
  Inferred_Gamma <- fit_glasso$wi
  upper_tri_elements <- upper.tri(Inferred_Gamma)
  Inferred_Edges <- which((abs(Inferred_Gamma) > 0)*upper_tri_elements > 0, arr.ind = TRUE)
  Inferred_Adjacency_Matrix <- matrix(0, nrow = d-1, ncol = d-1)
  Inferred_Adjacency_Matrix[Inferred_Edges] <- 
    Inferred_Adjacency_Matrix[cbind(Inferred_Edges[,2], Inferred_Edges[,1])] <- 1
  return(Inferred_Adjacency_Matrix)
}

## Function to combine the adjacency matrices of the sub-graphs into the adjacency
## matrix for the full graph
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
  for(i in 1:d){
    Inferred_Adjacency_Matrix[-i,-i] <- 
      Inferred_Adjacency_Matrix[-i,-i] + Sub_Matrices[[i]]
  }
  return(Inferred_Adjacency_Matrix)
}

## Obtain the graph from the adacjency matrix above
Graph_From_Adjacency_Matrix <- function(Adj_Matrix){
  if(!is.square.matrix(Adj_Matrix)){
    stop("Adj_Matrix must be a square matrix")
  }
  if(!is.symmetric.matrix(Adj_Matrix)){
    stop("Adj_Matrix must be a symmetric matrix")
  }
  upper_tri_elements <- upper.tri(Adj_Matrix)
  edges <- t(which(Adj_Matrix*upper_tri_elements > 0, arr.ind = TRUE))
  g_out <- graph(edges, directed = FALSE)
  return(g_out)
}

## Plot the inferred sub-graph when we condition on each site
Plot_Graph_From_Adjacency_Matrix_i <- function(Adj_Matrix, cond = 1, n_sim = 1, layout){
  if(!is.square.matrix(Adj_Matrix)){
    stop("Adj_Matrix must be a square matrix")
  }
  if(!is.symmetric.matrix(Adj_Matrix)){
    stop("Adj_Matrix must be a symmetric matrix")
  }
  if(cond <= 0 | cond%%1 != 0){
    stop("cond must be a postive integer")
  }
  if(n_sim <= 0 | n_sim%%1 != 0){
    stop("Number of samples must be a postive integer")
  }
  upper_tri_elements <- upper.tri(Adj_Matrix)
  edges <- which(Adj_Matrix*upper_tri_elements > 0, arr.ind = TRUE)
  g_out <- graph(t(edges), directed = FALSE)
  
  d <- ncol(Adj_Matrix) + 1
  V(g_out)$label <- NA
  V(g_out)$size <- 10
  E(g_out)$weight <- (Adj_Matrix[edges])/n_sim
  
  if(is_missing(layout)){
    layout <- layout_with_fr(graph)
  }
  plot(g_out, layout = layout[-i,], edge.width = 5*E(g_out)$weight, edge.color = "grey50")
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
    plot(g_out, layout = layout, edge.width = 1*E(g_out)$weight)
  }
  else{
    plot(g_out, layout = layout, edge.width = 1*E(g_out)$weight, edge.color = "grey") 
  }
}
