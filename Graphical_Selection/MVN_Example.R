################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "igraph", "graphicalExtremes", "parallel", "purrr")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in R scripts
source("Miscellaneous_Functions/Transformations.R")

## Reading in required functions
source("Graphical_Selection/Graphical_Selection.R")

###############################################################################
## Functions for the Engleke and Hitz model selection method

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

#function to permute the edges when order is not a factor
edges_permute <- function(x){
  set <- c()
  n <- length(x)
  if(n == 0L | n == 1){
    return(rep(NA, 2))
  }
  else{
    for(i in 1:(n - 1)){
      set <- rbind(set, cbind(rep(x[i], rep = n - i), x[-c(1:i)]))
    }
    return(set) 
  }
}

#functions to create a graph with maximum of 3 cliques
Graph_max_cliques_3 <- function(g_start, data, N){
  n_edges_miss <- N*(N-1)/2 - length(E(g_start))
  while(n_edges_miss > 0){
    #get single edges that can be added to the graph
    # singles <- which(as_adjacency_matrix(g_start) == 0, arr.ind = TRUE)
    # singles <- singles[-which(singles[,1] == singles[,2]),]
    #get the cliques
    cliques <- get_cliques_and_separators(g_start)
    #get all the cliques of size two
    cliques_2 <- t(sapply(cliques, function(x){if(length(x) == 2){x}else{rep(NA,2)}}))
    #get the edges from the cliques of size 2
    edges <- lapply(X = 1:N, function(x,y){y[c(which(y[,1] == x), which(y[,2] == x)),]}, y = cliques_2)
    #get the edges in the graph that would make a clique of max size three
    poss_edges <- pmap(.l = list(x = edges, y = 1:N), .f = function(x, y){x[which(x != y, arr.ind = TRUE)]})
    
    #get the final lot of edges that can be added to the graph
    final_edges <- lapply(X = poss_edges, FUN = edges_permute)
    final_edges <- do.call(rbind, final_edges)
    #get rid of na values
    final_edges <- final_edges[-unique(which(is.na(final_edges), arr.ind = TRUE)[,1]),]
    #add in the singles to the possible edges to be added
    # final_edges <- rbind(final_edges, singles)
    # #remove duplicates
    # final_edges <- uniquecombs(final_edges)
    #make into a list
    if(is.null(dim(final_edges))){
      final_edges <- list(c(final_edges))
    }
    else{
      final_edges <- lapply(apply(final_edges, 1, list), function(vec){vec[[1]]}) 
    }
    if(is_empty(final_edges)){
      n_edges_miss <- 0
      Gamma_star <- fits[[which.min(aic_alt)]]
    }
    else{
      #get the set of possible graphs
      poss_graphs <- lapply(X = final_edges, function(x){add_edges(g_start, x, mode = "undirected")})
      #fit the model
      fits <- lapply(X = poss_graphs, FUN = fmpareto_graph_HR, data = data, cens = TRUE)
      #obtain AIC of possible graphs
      log_likes_alt <- pmap(.l = list(graph = poss_graphs, Gamma = fits), .f = loglik_HR, data = data, cens = TRUE)
      
      #Fit the null model
      null_model <- fmpareto_graph_HR(data = data, graph = g_start, cens = TRUE)
      
      #compare teh AIC of null and alternate
      aic_null <- loglik_HR(data = data, graph = g_start, Gamma = null_model, cens = TRUE)
      aic_null <- aic_null[2]
      aic_alt <- sapply(log_likes_alt, function(x){x[2]})
      aic_alt_min <- aic_alt[which.min(aic_alt)]
      
      #set up loop for next iteration
      if(aic_null > aic_alt_min){
        n_edges_miss <- n_edges_miss - 1
        g_start <- poss_graphs[[which.min(aic_alt)]]
      }
      else{
        n_edges_miss <- 0
        aic_alt_min <- aic_null
        Gamma_star <- null_model
      }
    }
  }
  #get the output
  out <- list(g_star = g_start,
              Gamma_star = Gamma_star,
              AIC = aic_alt_min)
  return(out)
}

pairs_edges <- function(mat){
  n <- nrow(mat)
  combos <- combinations(n = n, r = 2, v = 1:n, repeats.allowed = FALSE)
  combos <- lapply(apply(combos, 1, list), function(vec){vec[[1]]})
  out <- lapply(combos, function(x){c(t(mat[x,]))})
  return(out)
}

#does the same as above but not in an iterative sense
Graph_max_cliques_3 <- function(g_start, data, N){
  #m is the number of cliques we want to add
  
  #get single edges that can be added to the graph
  # singles <- which(as_adjacency_matrix(g_start) == 0, arr.ind = TRUE)
  # singles <- singles[-which(singles[,1] == singles[,2]),]
  #get the cliques
  cliques <- get_cliques_and_separators(g_start)
  cliques <- order_cliques(cliques)
  #get all the cliques of size two
  cliques_2 <- t(sapply(cliques, function(x){if(length(x) == 2){x}else{rep(NA,2)}}))
  #get the edges from the cliques of size 2
  edges <- lapply(X = 1:N, function(x,y){y[c(which(y[,1] == x), which(y[,2] == x)),]}, y = cliques_2)
  #get the edges in the graph that would make a clique of max size three
  poss_edges <- pmap(.l = list(x = edges, y = 1:N), .f = function(x, y){x[which(x != y, arr.ind = TRUE)]})
  
  #get the final lot of edges that can be added to the graph
  final_edges <- lapply(X = poss_edges, FUN = edges_permute)
  final_edges <- do.call(rbind, final_edges)
  
  if(!all(is.na(final_edges))){
    #get rid of na values
    final_edges <- final_edges[-unique(which(is.na(final_edges), arr.ind = TRUE)[,1]),]
    #add in the singles to the possible edges to be added
    # final_edges <- rbind(final_edges, singles)
    # #remove duplicates
    # final_edges <- uniquecombs(final_edges)
    #make into a list
    if(is.null(dim(final_edges))){
      final_edges <- list(c(final_edges))
    }
    else{
      final_edges <- pairs_edges(final_edges)
    }
    final_edges_mat <- lapply(final_edges, function(x){matrix(x, ncol = 2, byrow = TRUE)})
    
    # #check if we have any repeats
    # reps <- diag(0, length(final_edges_mat))
    # for(i in 1:length(final_edges_mat)){
    #   edges_i <- final_edges_mat[[i]][order(final_edges_mat[[i]][,1], decreasing = FALSE),]
    #   for(j in 1:length(final_edges_mat)){
    #     if(i >= j){
    #       next()
    #     }
    #     else{
    #       edges_j <- final_edges_mat[[j]][order(final_edges_mat[[j]][,1], decreasing = FALSE),]
    #       if(all(edges_i == edges_j)){
    #         reps[i,j] <- reps[j,i] <- 1
    #       }
    #     }
    #   }
    # }
    
    #get the set of possible graphs
    poss_graphs <- lapply(X = final_edges, function(x){add_edges(g_start, x, mode = "undirected")})
    
    #remove any grpahs that have created more than a three clique
    cliques_poss_graphs <- sapply(poss_graphs, clique_num)
    if(any(cliques_poss_graphs > 3)){
      index <- which(cliques_poss_graphs > 3)
      poss_graphs <- sapply((1:length(poss_graphs))[-index], function(x){poss_graphs[[x]]})
      #this does not work - fix this!
    }
    
    #fit the model
    fits <- lapply(X = poss_graphs, FUN = fmpareto_graph_HR, data = data, method = "vario")
    #obtain AIC of possible graphs
    log_likes_alt <- pmap(.l = list(graph = poss_graphs, Gamma = fits), .f = loglik_HR, data = data)
    
    #compare teh AIC of null and alternate
    aic_alt <- sapply(log_likes_alt, function(x){x[2]})
  }
  
  #get the output
  out <- list(g_star = poss_graphs[[which.min(aic_alt)]],
              Gamma_star = fits[[which.min(aic_alt)]],
              AIC = aic_alt[which.min(aic_alt)])
  return(out)
}

################################################################################

## Create data for the simulation study

## True graph
nv <- d <- 16
connections <- diag(0, nrow = nv)
connections[2:6,1] <- connections[1,2:6] <- 1
connections[c(3,7,8),2] <- connections[2,c(3,7,8)] <- 1
connections[c(2,9,10),3] <- connections[3,c(2,9,10)] <- 1
connections[c(11,12),4] <- connections[4,c(11,12)] <- 1
connections[c(13,14),5] <- connections[5,c(13,14)] <- 1
connections[c(15,16),6] <- connections[6,c(15,16)] <- 1
connections[8,7] <- connections[7,8] <- 1
connections[15,16] <- connections[16,15] <- 1

edge_list <- which(connections == 1, arr.ind = TRUE)
edge_list <- edge_list[edge_list[,2] <= edge_list[,1],]
tree_sim <- graph_from_data_frame(edge_list, directed = FALSE)
layout <- layout_as_tree(tree_sim)

V(tree_sim)$label <- NA
V(tree_sim)$size <- 15

################################################################################
## DO NOT RUN

## number of simulations and data points
n_sim <- 100
n_data <- 1000

## Create the data
seed <- -1034475253
set.seed(seed)

## Precision matrix
simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(tree_sim)), v_sign = -1,
                           v_within = c(0.3, 0.7))
Gamma_true <- simul$omega
Sigma_true <- solve(Gamma_true)
rho_true <- cov2cor(Sigma_true)

## Mean parameter
mu_true <- runif(d, -5, 5)

## Simulate the data
X <- replicate(n = n_sim,
               expr = rmvn.sparse(n = n_data, mu = mu_true, CH = Cholesky(as(solve(rho_true), "sparseMatrix"))),
               simplify = FALSE)

## transform the data onto Laplace margins
X_list_by_data <- lapply(X, function(x){lapply(apply(x, 2, list), function(y){y[[1]]})})
X_list_by_var <- lapply(1:d, function(i){lapply(1:n_sim, function(j){X_list_by_data[[j]][[i]]})})
q_poss <- seq(0.55, 0.999, length.out = 100)
X_to_Y <- lapply(X_list_by_var, function(x){
  mcmapply(x = x,
           FUN = X_to_Laplace,
           MoreArgs = list(q = q_poss),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

## Get the output
u_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$par$u)})})
qu_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$par$qu)})})
scale_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$par$scale)})})
shape_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$par$shape)})})
Y <- lapply(1:n_sim, function(i){sapply(1:d, function(j){X_to_Y[[j]][[i]]$data$Y})})

print("Threshold Selection Complete")

################################################################################
## Save the output
# out <- list(transforms = X_to_Y,
#             par_true = list(n_sim = n_sim, n_data = n_data, graph = tree_sim, Gamma = Gamma_true))
# saveRDS(out, "Data/MVN_D16.RData")

################################################################################
## Read in the data
# out <- readRDS(file = "Data/MVN_D16.RData")
# 
# ## Extract the output
# tree_sim <- out$par_true$graph
# d <- length(V(tree_sim))
# 
# Gamma_true <- out$par_true$Gamma
# 
# n_sim <- out$par_true$n_sim
# n_data <- out$par_true$n_data
# 
# X <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$data$X)})})
# Y <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$data$Y)})})
# 
# u_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$u)})})
# qu_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$qu)})})
# scale_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$scale)})})
# shape_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$shape)})})

################################################################################
## tuning parameters
dqu <- 0.9
rho <- seq(0.2, 0.3, by = 0.01)

## Now we want to subset the data so that each component is large in turn
Y_u <- qlaplace(dqu)
Y_Yi_large <- rep(list(vector("list", d)), n_sim)
for(i in 1:n_sim){
  for(j in 1:d){
    Y_Yi_large[[i]][[j]] <- Y[[i]][which(Y[[i]][,j] > Y_u),] 
  }
}

## Fit the three-step model and infer the graph
n_rho <- length(rho)
Inferred_Subgraphs <- lapply(1:n_rho, function(j){
  lapply(1:n_sim, function(i){
    mcmapply(FUN = Infer_Adj_Matrix_glasso_residuals, 
             data = Y_Yi_large[[i]],
             cond = 1:d,
             MoreArgs = list(rho = rho[j],
                             v = ceiling(max(sapply(Y, max))) + 1),
             SIMPLIFY = FALSE)})})

## Check all the models have converged
## This may require some updates
index <- lapply(Inferred_Subgraphs, function(x){lapply(x, function(y){which(sapply(y, is.matrix) == FALSE)})})
count <- 0
while(any(sapply(index, function(x){!sapply(x, is_empty)}))){
  for(i in 1:n_rho){
    for(j in 1:n_sim){
      if(is_empty(index[[i]][[j]])){
        next()
      }
      else{
        ind <- index[[i]][[j]]
        for(k in 1:length(ind)){
          start_AGG <- c(0, runif(3, 0.5, 3))
          Inferred_Subgraphs[[i]][[j]][[ind[k]]] <- Infer_Adj_Matrix_glasso_residuals(data = Y_Yi_large[[j]][[ind[k]]],
                                                                                      cond = ind[k],
                                                                                      rho = rho[i],
                                                                                      start_AGG = start_AGG,
                                                                                      v = ceiling(max(sapply(Y, max))) + 1)
        }
      }
    }
  } 
  count <- count + 1
  if(count >= 100){
    stop("100 starting values attempted")
  }
  index <- lapply(Inferred_Subgraphs, function(x){lapply(x, function(y){which(sapply(y, is.matrix) == FALSE)})})
  if(!any(sapply(index, function(x){!sapply(x, is_empty)}))){
    print("Model Fitting Complete")
  }
}

################################################################################

## Get the weighted combined graphs for each conditioning variable
Inferred_Adj_Matrices_Weighted <- rep(list(vector("list", n_sim)), n_rho)
for(i in 1:n_rho){
  for(j in 1:n_sim){
    Inferred_Adj_Matrices_Weighted[[i]][[j]] <- Combine_Sub_Adjacency_Matrices(Inferred_Subgraphs[[i]][[j]])
  }
}

## Use majority rule to determine the graphs for each replicate
Inferred_Edges <- lapply(Inferred_Adj_Matrices_Weighted, function(x){lapply(x, function(y){which(y/(d-2) > 0.5, arr.ind = TRUE)})})

## combine the above to a single weighted matrix over the 100 replicates
Inferred_Adj_Matrices <- rep(list(rep(list(matrix(0, nrow = d, ncol = d)), n_sim)), n_rho)
for(i in 1:n_rho){
  for(j in 1:n_sim){
    Inferred_Adj_Matrices[[i]][[j]][Inferred_Edges[[i]][[j]]] <- 1
  } 
}
Inferred_Adj_Matrix_Weighted <- lapply(Inferred_Adj_Matrices, function(x){Reduce("+", x)})

## Inferred graph over the 100 replicates per rho
Inferred_Edges_per_rho <- lapply(Inferred_Adj_Matrix_Weighted, function(x){which(x/n_sim > 0.5, arr.ind = TRUE)})
Inferred_Adj_Matrix_per_rho <- rep(list(matrix(0, ncol = d, nrow = d)), n_rho)
for(i in 1:n_rho){
  Inferred_Adj_Matrix_per_rho[[i]][Inferred_Edges_per_rho[[i]]] <- 1 
}

## Combine the output
Inferred_Edges_Final <- which(Reduce("+", Inferred_Adj_Matrix_per_rho)/n_rho > 0.5, arr.ind = TRUE)
Inferred_Adj_Matrix_Final <- matrix(0, ncol = d, nrow = d)
Inferred_Adj_Matrix_Final[Inferred_Edges_Final] <- 1
Graph_Final <- graph_from_adjacency_matrix(Inferred_Adj_Matrix_Final, mode = "undirected")

print("SCMEVM Graphical Selection Complete")

################################################################################
## Engelke and Hitz model selection

## Transform onto MVP scale
## Threshold will be the largest amongst all GPD fits above

X_Pareto <- lapply(X, function(x){
  data2mpareto(data = x, p = 0.9)
})

## Get and fit the minimum spanning tree
emst_fits <- lapply(X_Pareto, function(x){
  emst(x, cens = TRUE)})
emst_fits_graphs <- lapply(emst_fits, function(x){x$graph})

# Graph_max_cliques_3(g_start = emst_fits_graphs[[1]], data = X_Pareto[[1]], N = d)

## Add edges into the model such that they maximally create a three clique
out <- mcmapply(FUN = Graph_max_cliques_3,
                g_start = emst_fits_graphs, 
                data = X_Pareto, 
                MoreArgs = list(N = d),
                SIMPLIFY = FALSE,
                mc.cores = detectCores() - 2)

print("Engelke and Hitz Graphical Selection Complete")

## Determine the number of edges in each graph
print(table(sapply(out, function(x){length(E(x$g_star))})))

## Obtain the graph
out_subgraphs_EHM <- lapply(out, function(x){igraph::subgraph_from_edges(x$g_star, E(x$g_star), delete.vertices = FALSE)})
adj_EHM <- lapply(out_subgraphs_EHM, as_adjacency_matrix)
adj_total_EHM <- Reduce("+", adj_EHM)
out_graph_EHM <- graph_from_adjacency_matrix(adj_total_EHM, mode = "undirected", weighted = TRUE)

## Add weights to the edges
weights <- c()
for(i in 1:nv){
  row_i <- adj_total_EHM[i,-c(1:i)]
  row_i <- row_i[row_i > 0]
  weights <- c(weights, row_i/n_sim)
}
E(out_graph_EHM)$weight <- weights

## Make higher weighted lines darker and thicker
E(out_graph_EHM)$color <- sapply(round(E(out_graph_EHM)$weight*100), function(x){paste0("grey", abs(100 - x))})

## Format the vertices
V(out_graph_EHM)$label <- NA
V(out_graph_EHM)$size <- 15

################################################################################

## Compare the SCMEVM and EHM methods

## Obtain the correct plot for the SCMEVM
upper_tri_elements <- upper.tri(Inferred_Adj_Matrix_Weighted[[1]])
edges <- which(Inferred_Adj_Matrix_Weighted[[1]]*upper_tri_elements > 0, arr.ind = TRUE)
edges <- edges[order(edges[,1]),]
out_graph_SCMEVM <- make_graph(t(edges), directed = FALSE)

V(out_graph_SCMEVM)$label <- NA
V(out_graph_SCMEVM)$size <- 15
E(out_graph_SCMEVM)$weight <- (Inferred_Adj_Matrix_Weighted[[1]][edges])/n_sim

## Make higher weighted lines darker and thicker
E(out_graph_SCMEVM)$color <- sapply(round(E(out_graph_SCMEVM)$weight*100), function(x){paste0("grey", abs(100 - x))})

## plot the graph
pdf(file = "Images/Simulation_Study/Graphical_Selection/MVN/Inferred_Graphs_Weighted.pdf", width = 10, height = 10)
par(mfrow = c(1, 3), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot.igraph(tree_sim, 
            layout = layout, 
            edge.width = 5,
            edge.color = "black")
plot.igraph(out_graph_SCMEVM, 
            layout = layout, 
            edge.width = 5*E(out_graph_SCMEVM)$weight)
plot.igraph(out_graph_EHM, 
            layout = layout, 
            edge.width = 5*E(out_graph_EHM)$weight)
dev.off()

################################################################################
## Compare the pruned graphs
Edges_Pruned_EHM <- which(adj_total_EHM/n_sim > 0.5, arr.ind = TRUE)
Adj_Matrix_Final_EHM <- matrix(0, ncol = d, nrow = d)
Adj_Matrix_Final_EHM[Edges_Pruned_EHM] <- 1
Graph_Final_EHM <- graph_from_adjacency_matrix(Adj_Matrix_Final_EHM, mode = "undirected")

## Plot the output
V(Graph_Final)$label <- V(Graph_Final_EHM)$label <- NA
V(Graph_Final)$size <- V(Graph_Final_EHM)$size <- 15

pdf(file = "Images/Simulation_Study/Graphical_Selection/MVN/Inferred_Graphs_Final.pdf", width = 10, height = 10)
par(mfrow = c(1, 3), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)

plot.igraph(tree_sim, 
            layout = layout, 
            edge.width = 5, 
            edge.color = "black")

plot.igraph(Graph_Final, 
            layout = layout, 
            edge.width = 5, 
            edge.color = "black")

plot.igraph(Graph_Final_EHM, 
            layout = layout, 
            edge.width = 5, 
            edge.color = "black")

dev.off()
