################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "glasso", "graphicalExtremes", "gridExtra", "gtools", "igraph", "jmuOutlier", "LaplacesDemon", "matrixcalc", "mev", "moments", "rgl", "rlang", "parallel", "pracma", "purrr", "texmex")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Set working directory
setwd("/Users/aidenfarrell/Documents/GitHub/Graphical_Conditional_Extremes/")

## Reading in required functions
source("Cond_Extremes_MVAGG_Residuals.R")
source("General_Functions.R")
source("MVAGG_Functions.R")
source("/Graphical_Selection/Graphical_Inference.R")

## read in threshold selection functions
source("threshold_selection_paper/helper_functions.R")
source("threshold_selection_paper/thresh_qq_metric.R")

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

# pdf(file = "/home/farrel11/PhD/Project_2/Code_for_Paper/Images/Graphical_Selection/True_Graph.pdf",
#     width = 10, height = 10)
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot.igraph(tree_sim, layout = layout, edge.width = 5, edge.color = "black")
# dev.off()

################################################################################
## DO NOT RUN

## Create the data
Gamma_true <- complete_Gamma(runif(n = nrow(edge_list), 0.5, 1), graph = tree_sim)

## number of simulations and data points
n_sim <- 5
n_data <- 100

## tuning parameters
dqu <- 0.7
rho <- 0.65

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
u_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){
  unname(X_to_Y[[j]][[i]]$u)})})
qu_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){
  unname(X_to_Y[[j]][[i]]$qu)})})
scale_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){
  unname(X_to_Y[[j]][[i]]$scale)})})
shape_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){
  unname(X_to_Y[[j]][[i]]$shape)})})
Y <- lapply(1:n_sim, function(i){sapply(1:d, function(j){
  unname(X_to_Y[[j]][[i]]$Y)})})

################################################################################
## Read in the data

out <- readRDS(file = "Data/Graphical_Selection/MVP_D16.Rdata")

## Extract the output
X <- out$data$X
Y <- out$data$Y

n_sim <- out$par_true$n_sim
n_data <- out$par_true$n_data
dqu <- out$par_true$dqu
rho <- out$par_true$rho

edge_list <- out$par_true$edge_list
tree_sim <- out$par_true$graph

Gamma_true <- out$par_true$Gamma

u_final <- out$GPD_par$u
qu_final <- out$GPD_par$qu
scale_final <- out$GPD_par$scale
shape_final <- out$GPD_par$shape
################################################################################

## Now we want to subset the data so that each component is large in turn
Y_u <- qlaplace(dqu)

Y_Yi_large <- rep(list(vector("list", d)), n_sim)
for(i in 1:n_sim){
  for(j in 1:d){
    Y_Yi_large[[i]][[j]] <- Y[[i]][which(Y[[i]][,j] > Y_u),] 
  }
}

## Fit the three-step model and infer the graph
Inferred_Subgraphs <- lapply(1:n_sim, function(i){
  mcmapply(FUN = Infer_Adj_Matrix, 
           data = Y_Yi_large[[i]],
           cond = 1:d,
           MoreArgs = list(rho = rho),
           SIMPLIFY = FALSE)})

## Get the weighted combined graphs for each conditioning variable
Inferred_Adj_Matrices_Weighted <- lapply(Inferred_Subgraphs, Combine_Sub_Adjacency_Matrices)

## Use majority rule to determine the graphs for each replicate
Inferred_Edges <- lapply(Inferred_Adj_Matrices_Weighted, function(x){which(x/(d-2) > 0.5, arr.ind = TRUE)})

## combine the above to a single weighted matrix over the 100 replicates
Inferred_Adj_Matrices <- rep(list(matrix(0, nrow = d, ncol = d)), n_sim)
for(i in 1:n_sim){
  Inferred_Adj_Matrices[[i]][Inferred_Edges[[i]]] <- 1
}
Inferred_Adj_Matrix_Weighted <- Reduce("+", Inferred_Adj_Matrices)

## Plot the weighted graph
# pdf(file = "/home/farrel11/PhD/Project_2/Code_for_Paper/Images/Graphical_Selection/Inferred_Weighted_Graph.pdf",
#     width = 10, height = 10)
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
Plot_Graph_From_Adjacency_Matrix(Adj_Matrix = Inferred_Adj_Matrix_Weighted,
                                 n_sim = n_sim,
                                 layout = layout,
                                 true_graph = tree_sim)
# dev.off()

## Inferred graph over the 100 replicates
Inferred_Edges_Final <- which(Inferred_Adj_Matrix_Weighted/n_sim > 0.5, arr.ind = TRUE)
Inferred_Adj_Matrix_Final <- matrix(0, ncol = d, nrow = d)
Inferred_Adj_Matrix_Final[Inferred_Edges_Final] <- 1
Graph_Final <- graph_from_adjacency_matrix(Inferred_Adj_Matrix_Final, mode = "undirected")

V(Graph_Final)$label <- NA
V(Graph_Final)$size <- 15

par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot.igraph(Graph_Final, layout = layout, edge.width = 5, edge.color = "black")
