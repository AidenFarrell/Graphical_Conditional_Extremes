################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("graphicalExtremes", "parallel")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in R scripts
source("Miscellaneous_Functions/Transformations.R")

## Reading in required functions
source("Graphical_Selection/Graphical_Selection.R")

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

par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Simulation_Study/Graphical_Selection/True_Graph.pdf", width = 10, height = 10)
plot.igraph(tree_sim, layout = layout, edge.width = 5, edge.color = "black")
dev.off()

################################################################################
## DO NOT RUN

## number of simulations and data points
n_sim <- 100
n_data <- 100

## Create the data
seed <- -1034475253
set.seed(seed)
Gamma_true <- complete_Gamma(runif(n = nrow(edge_list), 0.5, 1), graph = tree_sim)

## Simulate the data
X <- replicate(n = n_sim,
               expr = rmpareto(n = n_data, model = "HR", par = Gamma_true),
               simplify = FALSE)

## Transform the data onto standard Laplace margins
X_list_by_data <- lapply(X, function(x){
  lapply(apply(x, 2, list), function(y){y[[1]]})})
X_list_by_var <- lapply(1:d, function(i){lapply(1:n_sim, function(j){X_list_by_data[[j]][[i]]})})

X_to_Y <- lapply(1:d, function(i){
  mcmapply(FUN = X_to_Laplace,
           x = X_list_by_var[[i]],
           MoreArgs = list(q = seq(0.55, 0.99, by = 0.01)),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

## Get the output
Y <- lapply(1:n_sim, function(i){sapply(1:d, function(j){
  unname(X_to_Y[[j]][[i]]$data$Y)})})

################################################################################
## Read in the data
out <- readRDS(file = "Data/MVP_D16.RData")

## Extract the output
tree_sim <- out$par_true$graph
d <- length(V(tree_sim))

Gamma_true <- out$par_true$Gamma

n_sim <- out$par_true$n_sim
n_data <- out$par_true$n_data

X <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$data$X)})})
Y <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$data$Y)})})

u_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$u)})})
qu_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$qu)})})
scale_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$scale)})})
shape_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$shape)})})
################################################################################

## tuning parameters
dqu <- 0.7
rho <- seq(from = 0.6, to = 0.7, by = 0.01)

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
for(i in 1:n_rho){
  for(j in 1:n_sim){
    if(is_empty(index[[i]][[j]])){
      next()
    }
    else{
      ind <- index[[i]][[j]]
      for(k in 1:length(ind)){
        Inferred_Subgraphs[[i]][[j]][[ind[k]]] <- Infer_Adj_Matrix_glasso_residuals(data = Y_Yi_large[[j]][[ind[k]]],
                                                                                    cond = ind[k],
                                                                                    rho = rho[i],
                                                                                    start_AGG = c(0, 2, 2, 1.5),
                                                                                    v = ceiling(max(sapply(Y, max))) + 1)
      }
    }
  }
  print(paste0(i, " of ", n_rho, " values of rho tested"))
}

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

## Plot the weighted graphs
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Simulation_Study/Graphical_Selection/Inferred_Graph_Rho_060.pdf", width = 10, height = 10)
Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = Inferred_Adj_Matrix_Weighted[[1]],
                                 n_sim = n_sim,
                                 layout = layout,
                                 true_graph = tree_sim)
dev.off()

## plot the remaining graphs
sapply(Inferred_Adj_Matrix_Weighted, function(x){
  Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = x,
                                   n_sim = n_sim,
                                   layout = layout,
                                   true_graph = tree_sim)
})

## Inferred graph over the 100 replicates per rho
Inferred_Edges_per_rho <- lapply(Inferred_Adj_Matrix_Weighted, function(x){which(x/n_sim > 0.5, arr.ind = TRUE)})
Inferred_Adj_Matrix_per_rho <- rep(list(matrix(0, ncol = d, nrow = d)), n_rho)
for(i in 1:n_rho){
  Inferred_Adj_Matrix_per_rho[[i]][Inferred_Edges_per_rho[[i]]] <- 1 
}

Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = Reduce("+", Inferred_Adj_Matrix_per_rho),
                                 n_sim = n_rho,
                                 layout = layout,
                                 true_graph = tree_sim)

Inferred_Edges_Final <- which(Reduce("+", Inferred_Adj_Matrix_per_rho)/n_rho > 0.5, arr.ind = TRUE)
Inferred_Adj_Matrix_Final <- matrix(0, ncol = d, nrow = d)
Inferred_Adj_Matrix_Final[Inferred_Edges_Final] <- 1
Graph_Final <- graph_from_adjacency_matrix(Inferred_Adj_Matrix_Final, mode = "undirected")

V(Graph_Final)$label <- NA
V(Graph_Final)$size <- 15

par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot.igraph(Graph_Final, layout = layout, edge.width = 5, edge.color = "black")
## We can see we get the correct graph out by using majority rule over the various values of rho

################################################################################
## Infer the graph using graphical lasso on the data

## Transform the data onto standard Gaussian margins
X_Gaussian <- lapply(1:n_sim, function(i){
  sapply(1:d, function(j){
  qnorm(semiparametric_unif_transform(X[[i]][,j], par = c(u_final[[i]][j], scale_final[[i]][j], shape_final[[i]][j])))})})

## tuning parameters
dqu <- 0.7
rho <- seq(from = 0.75, to = 0.85, by = 0.01)

## Now we want to subset the data so that each component is large in turn
X_u <- qnorm(dqu)

## Infer the graph from the graphical lasso on the data
n_rho <- length(rho)
Inferred_Subgraphs <- lapply(1:n_rho, function(j){
  lapply(1:n_sim, function(i){
    mcmapply(FUN = Infer_Adj_Matrix_glasso_data, 
             cond = 1:d,
             MoreArgs = list(data = X_Gaussian[[i]],
                             rho = rho[j],
                             u = X_u),
             SIMPLIFY = FALSE)})})

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

## Plot the weighted graphs
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Simulation_Study/Graphical_Selection/Inferred_Graph_Glasso_Data_Rho_075.pdf", width = 10, height = 10)
Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = Inferred_Adj_Matrix_Weighted[[1]],
                                 n_sim = n_sim,
                                 layout = layout,
                                 true_graph = tree_sim)
dev.off()

## plot the remaining graphs
sapply(Inferred_Adj_Matrix_Weighted, function(x){
  Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = x,
                                   n_sim = n_sim,
                                   layout = layout,
                                   true_graph = tree_sim)
})

## Inferred graph over the 100 replicates per rho
Inferred_Edges_per_rho <- lapply(Inferred_Adj_Matrix_Weighted, function(x){which(x/n_sim > 0.5, arr.ind = TRUE)})
Inferred_Adj_Matrix_per_rho <- rep(list(matrix(0, ncol = d, nrow = d)), n_rho)
for(i in 1:n_rho){
  Inferred_Adj_Matrix_per_rho[[i]][Inferred_Edges_per_rho[[i]]] <- 1 
}

Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = Reduce("+", Inferred_Adj_Matrix_per_rho),
                                 n_sim = n_rho,
                                 layout = layout,
                                 true_graph = tree_sim)

Inferred_Edges_Final <- which(Reduce("+", Inferred_Adj_Matrix_per_rho)/n_rho > 0.5, arr.ind = TRUE)
Inferred_Adj_Matrix_Final <- matrix(0, ncol = d, nrow = d)
Inferred_Adj_Matrix_Final[Inferred_Edges_Final] <- 1
Graph_Final <- graph_from_adjacency_matrix(Inferred_Adj_Matrix_Final, mode = "undirected")

V(Graph_Final)$label <- NA
V(Graph_Final)$size <- 15

par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot.igraph(Graph_Final, layout = layout, edge.width = 5, edge.color = "black")
## We can see we get the correct graph out by using majority rule over the various values of rho

################################################################################
## tuning parameters
dqu <- 0.7
alpha <- 0.1

## Now we want to subset the data so that each component is large in turn
X_u <- qnorm(dqu)

## Infer the graph from the graphical lasso on the data
n_alpha <- length(alpha)
Inferred_Subgraphs <- lapply(1:n_alpha, function(j){
  lapply(1:n_sim, function(i){
    mcmapply(FUN = Infer_Adj_Matrix_ppcor_data, 
             cond = 1:d,
             MoreArgs = list(data = X_Gaussian[[i]],
                             alpha = alpha[j],
                             u = X_u),
             SIMPLIFY = FALSE)})})

## Get the weighted combined graphs for each conditioning variable
Inferred_Adj_Matrices_Weighted <- rep(list(vector("list", n_sim)), n_alpha)
for(i in 1:n_alpha){
  for(j in 1:n_sim){
    Inferred_Adj_Matrices_Weighted[[i]][[j]] <- Combine_Sub_Adjacency_Matrices(Inferred_Subgraphs[[i]][[j]])
  }
}

## Use majority rule to determine the graphs for each replicate
Inferred_Edges <- lapply(Inferred_Adj_Matrices_Weighted, function(x){lapply(x, function(y){which(y/(d-2) > 0.5, arr.ind = TRUE)})})

## combine the above to a single weighted matrix over the 100 replicates
Inferred_Adj_Matrices <- rep(list(rep(list(matrix(0, nrow = d, ncol = d)), n_sim)), n_alpha)
for(i in 1:n_alpha){
  for(j in 1:n_sim){
    Inferred_Adj_Matrices[[i]][[j]][Inferred_Edges[[i]][[j]]] <- 1
  } 
}
Inferred_Adj_Matrix_Weighted <- lapply(Inferred_Adj_Matrices, function(x){Reduce("+", x)})

## Plot the weighted graphs
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Simulation_Study/Graphical_Selection/Inferred_Graph_ppcor_Data_Alpha_01.pdf", width = 10, height = 10)
Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = Inferred_Adj_Matrix_Weighted[[1]],
                                 n_sim = n_sim,
                                 layout = layout,
                                 true_graph = tree_sim)
dev.off()

## plot the remaining graphs
sapply(Inferred_Adj_Matrix_Weighted, function(x){
  Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = x,
                                   n_sim = n_sim,
                                   layout = layout,
                                   true_graph = tree_sim)
})

## Inferred graph over the 100 replicates per rho
Inferred_Edges_per_alpha <- lapply(Inferred_Adj_Matrix_Weighted, function(x){which(x/n_sim > 0.5, arr.ind = TRUE)})
Inferred_Adj_Matrix_per_alpha <- rep(list(matrix(0, ncol = d, nrow = d)), n_alpha)
for(i in 1:n_alpha){
  Inferred_Adj_Matrix_per_alpha[[i]][Inferred_Edges_per_alpha[[i]]] <- 1 
}

Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = Inferred_Adj_Matrix_per_alpha[[1]],
                                 n_sim = n_alpha,
                                 layout = layout,
                                 true_graph = tree_sim)
## Does not infer the graph well
################################################################################
## Save the output
# out <- list(transforms = X_to_Y,
#             par_true = list(n_sim = n_sim, n_data = n_data, graph = tree_sim, Gamma = Gamma_true))
# saveRDS(out, "Data/MVP_D16.RData")
################################################################################
