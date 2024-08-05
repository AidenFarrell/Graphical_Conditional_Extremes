################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("graphicalExtremes", "parallel")
# install.packages(required_pckgs, dependencies = TRUE)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in functions to transform the data
source("Miscellaneous_Functions/Transformations.R")

## Reading in functions to infer the graph using a glasso on the residuals of the model
source("Graphical_Selection/Graphical_Selection_Glasso_Residuals.R")

## Reading in functions required for model fitting
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")

## Read in functions for prediction
source("Prediction/Conditonal_Probability_Calculations.R")
source("Prediction/Sim_Surfaces.R")

################################################################################
## Plot the data
danube_graph <- graph(t(danube$flow_edges), directed = FALSE)
layout <- cbind(danube$info$PlotCoordX, danube$info$PlotCoordY)

pdf(file = "Images/Danube/Danube_River.pdf", width = 10, height = 10)
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot(danube_graph, layout = layout)
dev.off()

d <- length(V(danube_graph))
flow_connections <- matrix(0, nrow = d, ncol = d)
flow_connections[1,-1] <- 1
flow_connections[2,-c(1:2,13,28:31)] <- 1
flow_connections[3,-c(1:3,13,14:19,28:31)] <- 1
flow_connections[4,-c(1:4,13,14:19,28:31)] <- 1
flow_connections[5,c(6:12,20:22)] <- 1
flow_connections[6,c(7:12,20:22)] <- 1
flow_connections[7,c(8:12,20:22)] <- 1
flow_connections[8,c(9,10,11,12)] <- 1
flow_connections[9,c(10,11,12)] <- 1
flow_connections[10,c(11,12)] <- 1
flow_connections[11,12] <- 1
flow_connections[13,28:31] <- 1
flow_connections[14,15:19] <- 1
flow_connections[15,16:19] <- 1
flow_connections[16,17:19] <- 1
flow_connections[17,18:19] <- 1
flow_connections[18,19] <- 1
flow_connections[20,21:22] <- 1
flow_connections[21,22] <- 1
flow_connections[23,24] <- 1
flow_connections[25,26:27] <- 1
flow_connections[26,27] <- 1
flow_connections[28,29] <- 1
flow_connections[30,c(28,29,31)] <- 1
flow_connections[31,28:29] <- 1

lower_tri_flow <- lower.tri(flow_connections)
flow_connections[lower_tri_flow] <- t(flow_connections)[lower_tri_flow]
is.symmetric.matrix(flow_connections)

################################################################################
## Transform the data onto standard Laplace margins
d <- ncol(danube$data_clustered)
n_data <- nrow(danube$data_clustered)
n_sim <- 1

danube_data_list <- lapply(apply(danube$data_clustered, 2, list), function(x){x[[1]]})
X_to_Y <- mcmapply(FUN = X_to_Laplace,
                   x = danube_data_list,
                   MoreArgs = list(q = seq(0.55, 0.925, by = 0.01)),
                   SIMPLIFY = FALSE,
                   mc.cores = detectCores() - 1)

## Get the output
u_final <- unname(sapply(X_to_Y, function(x){unname(x$par$u)}))
qu_final <- unname(sapply(X_to_Y, function(x){unname(x$par$qu)}))
scale_final <- unname(sapply(X_to_Y, function(x){unname(x$par$scale)}))
shape_final <- unname(sapply(X_to_Y, function(x){unname(x$par$shape)}))
Y <- sapply(X_to_Y, function(x){x$data$Y})

################################################################################
## Infer the graph from the data

## Use partial correlation test

## Need to transform to standard Gaussian margins first
danube_norm <- sapply(1:d, function(i){
  qnorm(semiparametric_unif_transform(danube$data_clustered[,i], par = c(u_final[i], scale_final[i], shape_final[i])))})

## Condition on each variable being large in turn
dqu <- 0.8
u <- qnorm(dqu)
danube_norm_large <- vector("list", d)
for(k in 1:d){
  danube_norm_large[[k]] <- danube_norm[which(danube_norm[,k] > u),]
}

## use partial correlation to infer the graphical structure
alpha <- 0.1
Inferred_Edges_i <- lapply(1:d, function(j){
  which(ppcor::pcor(danube_norm_large[[j]][-j,-j])$p.value < alpha/2, arr.ind = TRUE)})
Adj_Matrices_i <- rep(list(matrix(0, d-1, d-1)), d)
for(i in 1:d){
  Adj_Matrices_i[[i]][Inferred_Edges_i[[i]]] <- 1
  diag(Adj_Matrices_i[[i]]) <- 0
}
Adj_Matrix_Total <- Combine_Sub_Adjacency_Matrices(Adj_Matrices_i)
Inferred_Edges_Final <- which(Adj_Matrix_Total > (d-1)/2, arr.ind = TRUE)
Adj_Matrix_Final <- matrix(0, d, d)
Adj_Matrix_Final[Inferred_Edges_Final] <- 1

Graph_Final <- graph_from_adjacency_matrix(Adj_Matrix_Final, mode = "undirected")
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Danube/Inferred_Graph_ppcor.pdf", width = 10, height = 10)
plot(Graph_Final, layout = layout)
dev.off()

################################################################################
## Infer the graph from the data

## Use glasso on the data
rho <- seq(0.55, 0.85, by = 0.01)
n_rho <- length(rho)
Inferred_Adj_Matrices <- rep(list(rep(list(matrix(0, d-1, d-1)), d)), n_rho)
for(i in 1:n_rho){
  for(j in 1:d){
    glasso_fit <- glasso(s = cor(danube_norm_large[[1]][-1,-1]), rho = rho[i],
                         penalize.diagonal = FALSE, thr = 1e-8, maxit = 1e+6)
    Inferred_edges <- which(abs(glasso_fit$wi) > 0, arr.ind = TRUE)
    Inferred_Adj_Matrices[[i]][[j]][Inferred_edges] <- 1
    diag(Inferred_Adj_Matrices[[i]][[j]]) <- 0
  }
}

## Get the inferred weighted adjacency matrices per rho
Inferred_Adj_Matrices_per_rho <- lapply(Inferred_Adj_Matrices, Combine_Sub_Adjacency_Matrices)

## Prune the edges for each value of rho
Inferred_Edges_per_rho <- lapply(Inferred_Adj_Matrices_per_rho, function(x){which(x > (d-2)/2, arr.ind = TRUE)})

## Get the output as an adjacnecy matrix
Inferred_Sub_Matrices_per_rho <- rep(list(matrix(0, d, d)), n_rho)
for(i in 1:n_rho){
  Inferred_Sub_Matrices_per_rho[[i]][Inferred_Edges_per_rho[[i]]] <- 1
}

## From the pruned graphs per rho, prune again over the various rho
Inferred_Adj_Matrix_Weighted_Final <- Reduce("+", Inferred_Sub_Matrices_per_rho)
Inferred_Edges_Final <- which(Inferred_Adj_Matrix_Weighted_Final > n_rho/2, arr.ind = TRUE)
Inferred_Adj_Matrix_Final <- matrix(0, d, d)
Inferred_Adj_Matrix_Final[Inferred_Edges_Final] <- 1

## Final graph
Graph_Final <- graph_from_adjacency_matrix(Inferred_Adj_Matrix_Final, mode = "undirected")
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Danube/Inferred_Graph_Glasso_Data.pdf", width = 10, height = 10)
plot(Graph_Final, layout = layout)
dev.off()

################################################################################

## Now we want to subset the data so that each component is large in turn
dqu <- 0.8
Y_u <- qlaplace(dqu)

Y_Yi_large <- rep(list(list()), d)
for(i in 1:d){
  Y_Yi_large[[i]] <- Y[which(Y[,i] > Y_u),]
}

## Infer the graphical structure from the residuals of the CMEVM 
## over various values of rho
rho <- seq(0.55, 0.85, by = 0.01)
n_rho <- length(rho)
Inferred_Adj_Matrices <- lapply(rho, function(x){
  mcmapply(FUN = Infer_Adj_Matrix,
           data = Y_Yi_large,
           cond = 1:d,
           MoreArgs = list(rho = x,
                           v = ceiling(max(Y)) + 1),
           mc.cores = detectCores() - 1,
           SIMPLIFY = FALSE)
})

## Get the inferred weighted adjacency matrices per rho
Inferred_Adj_Matrices_per_rho <- lapply(Inferred_Adj_Matrices, Combine_Sub_Adjacency_Matrices)

## Prune the edges for each value of rho
Inferred_Edges_per_rho <- lapply(Inferred_Adj_Matrices_per_rho, function(x){which(x > (d-2)/2, arr.ind = TRUE)})

## Get the output as an adjacnecy matrix
Inferred_Sub_Matrices_per_rho <- rep(list(matrix(0, d, d)), n_rho)
for(i in 1:n_rho){
  Inferred_Sub_Matrices_per_rho[[i]][Inferred_Edges_per_rho[[i]]] <- 1
}

## From the pruned graphs per rho, prune again over the various rho
Inferred_Adj_Matrix_Weighted_Final <- Reduce("+", Inferred_Sub_Matrices_per_rho)
Inferred_Edges_Final <- which(Inferred_Adj_Matrix_Weighted_Final > n_rho/2, arr.ind = TRUE)
Inferred_Adj_Matrix_Final <- matrix(0, d, d)
Inferred_Adj_Matrix_Final[Inferred_Edges_Final] <- 1

## Final graph
Graph_Final <- graph_from_adjacency_matrix(Inferred_Adj_Matrix_Final, mode = "undirected")
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Danube/Inferred_Graph_Glasso_Residuals.pdf", width = 10, height = 10)
plot(Graph_Final, layout = layout)
dev.off()

################################################################################
## Comapre models
fit_graph <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
                      data = Y_Yi_large,
                      cond = 1:d,
                      MoreArgs = list(graph = Graph_Final,
                                      v = ceiling(max(Y)) + 1),
                      SIMPLIFY = FALSE,
                      mc.cores = detectCores() - 1)

fit_full <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
                     data = Y_Yi_large,
                     cond = 1:d,
                     MoreArgs = list(graph = make_full_graph(n = d),
                                     v = ceiling(max(Y)) + 1),
                     SIMPLIFY = FALSE,
                     mc.cores = detectCores() - 1)

sim_graph <- replicate(n = 200,
                       expr = Sim_Surface_MVAGG(n_sim = 10*n_data,
                                                q = 0.8,
                                                transforms = X_to_Y,
                                                CMEVM_fits = fit_graph),
                       simplify = FALSE)

sim_full <- replicate(n = 200,
                      expr =  Sim_Surface_MVAGG(n_sim = 10*n_data,
                                                q = 0.8,
                                                transforms = X_to_Y,
                                                CMEVM_fits = fit_full),
                      simplify = FALSE)

u <- c(0.8, 0.85, 0.9)
eta_data <- array(NA, dim = c(d, d, length(u)))
for(i in 1:d){
  for(j in 1:d){
    tail_data <- mev::taildep(data = cbind(danube$data_clustered[,i], danube$data_clustered[,j]), 
                              u = u,
                              plot = FALSE)
    eta_data[i,j,] <- tail_data$eta[,1]
  }
}

eta_graph <- eta_full <- array(NA, dim = c(d, d, length(u), 200))
for(i in 1:d){
  for(j in 1:d){
    for(k in 1:200){
      tail_graph <- mev::taildep(data = cbind(sim_graph[[k]]$Data_Margins[,i], sim_graph[[k]]$Data_Margins[,j]), 
                                 u = u,
                                 plot = FALSE)
      eta_graph[i,j,,k] <- tail_graph$eta[,1]
      
      tail_full <- mev::taildep(data = cbind(sim_full[[k]]$Data_Margins[,i], sim_full[[k]]$Data_Margins[,j]), 
                                u = u,
                                plot = FALSE)
      eta_full[i,j,,k] <- tail_full$eta[,1] 
    }
  }
}

eta_graph_med <- eta_full_med <- array(NA, dim = c(d, d, length(u)))
for(i in 1:d){
  for(j in 1:d){
    eta_graph_med[i,j,] <- apply(eta_graph[i,j,,], 1, quantile, 0.5)
    eta_full_med[i,j,] <- apply(eta_full[i,j,,], 1, quantile, 0.5)
  }
}

eta_graph_se <- eta_full_se <- array(NA, dim = c(d, d, length(u)))
for(i in 1:d){
  for(j in 1:d){
    eta_graph_se[i,j,] <- apply(eta_graph[i,j,,], 1, sd)/sqrt(200)
    eta_full_se[i,j,] <- apply(eta_full[i,j,,], 1, sd)/sqrt(200)
  }
}

## Sites that are connected (lowest site first)
flow_connected <- which(flow_connections == 1, arr.ind = TRUE)
flow_connected <- flow_connected[which(flow_connected[,1] < flow_connected[,2]),]
flow_connected <- flow_connected[order(flow_connected[,1]),]
flow_connected <- as.data.frame(flow_connected)

pairs <- as.data.frame(combinations(n = d, r = 2, v = 1:d))
pairs$exists <- do.call(paste0, pairs) %in% do.call(paste0, flow_connected)

for(i in 1:length(u)){
  eta_comp <- cbind(pairs, 
                    eta_data[,,i][lower_tri_flow], 
                    eta_graph_med[,,i][lower_tri_flow],
                    eta_full_med[,,i][lower_tri_flow])
  colnames(eta_comp) <- c("Station_1", "Station_2", "Flow_Connected", "Eta_Data", "Eta_Graph", "Eta_Full")
  
  par(mfrow = c(1, 2))
  plot(x = eta_comp$Eta_Data[eta_comp$Flow_Connected == TRUE],
       y = eta_comp$Eta_Graph[eta_comp$Flow_Connected == TRUE],
       col = "blue", pch = 4,
       xlim = c(1/2, 1), ylim = c(1/2, 1),
       xlab = "Empirical", ylab = "Theoretical", main = substitute(eta(u) ~ "- Graph", list(u = u[i])))
  points(x = eta_comp$Eta_Data[eta_comp$Flow_Connected == FALSE],
         y = eta_comp$Eta_Graph[eta_comp$Flow_Connected == FALSE],
         col = "black", pch = 19)
  abline(a = 0, b = 1, col = 2, lty = 2, lwd = 2)
  
  plot(x = eta_comp$Eta_Data[eta_comp$Flow_Connected == TRUE],
       y = eta_comp$Eta_Full[eta_comp$Flow_Connected == TRUE],
       col = "blue", pch = 4,
       xlim = c(1/2, 1), ylim = c(1/2, 1),
       xlab = "Empirical", ylab = "Theoretical", main = substitute(eta(u) ~ "- Full", list(u = u[i])))
points(x = eta_comp$Eta_Data[eta_comp$Flow_Connected == FALSE],
       y = eta_comp$Eta_Full[eta_comp$Flow_Connected == FALSE],
       col = "black", pch = 19)
abline(a = 0, b = 1, col = 2, lty = 2, lwd = 2) 
}

## Look at where the bias is
bias_TS_Graph_boot <- eta_graph_med[,,1] - eta_data[,,1]

x_site <- 1:d
y_site <- 1:d
sites <- expand.grid(x_site, y_site)
sites <- sites[order(sites[,1]),]
sites <- cbind(sites, do.call(paste0, sites[,1:2]) %in% do.call(paste0, flow_connected))
bias_TS_Graph_boot_plot_data <- as.data.frame(cbind(sites, c(bias_TS_Graph_boot)))
colnames(bias_TS_Graph_boot_plot_data) <- c("Conditioning_Site", "Dependent_Site", "Connected", "Value")
ggplot(data = bias_TS_Graph_boot_plot_data, aes(x = Conditioning_Site, y = Dependent_Site, fill = Value, col = Connected)) + 
  scale_x_continuous(breaks=seq(1, d, 1)) + 
  scale_y_continuous(breaks=seq(1, d, 1)) +
  geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
## main bias is in 23 - 27
## Secondary is in 11-12, 20 -22, 16 - 19 and 28 - 30 (basically the mountain regions)