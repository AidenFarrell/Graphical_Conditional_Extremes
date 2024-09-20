################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("ggplot2", "graphicalExtremes", "mev", "parallel")
# install.packages(required_pckgs, dependencies = TRUE)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in functions to transform the data
source("Miscellaneous_Functions/Transformations.R")

## Reading in functions to infer the graph using a glasso on the residuals of the model
source("Graphical_Selection/Graphical_Selection.R")

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

## Get the edges in the graph
danube_edges <- which(as_adjacency_matrix(danube_graph) > 0, arr.ind = TRUE)
danube_edges <- danube_edges[which(danube_edges[,1] < danube_edges[,2]),]
danube_edges <- danube_edges[order(danube_edges[,1]),]
danube_edges <- as.data.frame(danube_edges)

## Get the flow connections for later
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

## Threshold for the tests
dqu <- 0.8
u <- qnorm(dqu)

## use partial correlation to infer the graphical structure
alpha <- 0.1
Adj_Matrices_i <- lapply(1:d, function(j){
  Infer_Adj_Matrix_ppcor_data(data = danube_norm, u = u, cond = j, alpha = alpha)
})
Adj_Matrix_Total <- Combine_Sub_Adjacency_Matrices(Adj_Matrices_i)
Inferred_Edges_Final <- which(Adj_Matrix_Total > (d-1)/2, arr.ind = TRUE)
Adj_Matrix_Final <- matrix(0, d, d)
Adj_Matrix_Final[Inferred_Edges_Final] <- 1

Graph_Final <- graph_from_adjacency_matrix(Adj_Matrix_Final, mode = "undirected")
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Danube/Inferred_Graph/Inferred_Graph_ppcor.pdf", width = 10, height = 10)
plot(Graph_Final, layout = layout)
dev.off()

## Check the percentage of correctly identified edges
Inferred_Edges_Final <- Inferred_Edges_Final[which(Inferred_Edges_Final[,1] < Inferred_Edges_Final[,2]),]
Inferred_Edges_Final <- Inferred_Edges_Final[order(Inferred_Edges_Final[,1]),]
Inferred_Edges_Final <- as.data.frame(Inferred_Edges_Final)
Inferred_Edges_Final$Exists <- do.call(paste0, Inferred_Edges_Final) %in% do.call(paste0, danube_edges)
sum(Inferred_Edges_Final$Exists)/nrow(danube_edges) 

################################################################################
## Infer the graph using a glasso on the data

## Infer the graphical structure over various values of rho
# rho <- seq(0.7, 0.9, by = 0.01)
rho <- 0.9
n_rho <- length(rho)
Inferred_Adj_Matrices <- rep(list(rep(list(matrix(0, d-1, d-1)), d)), n_rho)
for(i in 1:n_rho){
  for(j in 1:d){
    Inferred_Adj_Matrices[[i]][[j]] <-
      Infer_Adj_Matrix_glasso_data(data = danube_norm, u = u, cond = j, rho = rho[i]) 
  }
}

## Get the inferred weighted adjacency matrices per rho
Inferred_Adj_Matrices_per_rho <- lapply(Inferred_Adj_Matrices, Combine_Sub_Adjacency_Matrices)

## Prune the edges for each value of rho
Inferred_Edges_per_rho <- lapply(Inferred_Adj_Matrices_per_rho, function(x){which(x > (d-2)/2, arr.ind = TRUE)})

## Get the output as an adjacency matrix
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
pdf(file = "Images/Danube/Inferred_Graph/Inferred_Graph_Glasso_Data.pdf", width = 10, height = 10)
plot(Graph_Final, layout = layout)
dev.off()

## Check the percentage of correctly identified edges
Inferred_Edges_Final <- Inferred_Edges_Final[which(Inferred_Edges_Final[,1] < Inferred_Edges_Final[,2]),]
Inferred_Edges_Final <- Inferred_Edges_Final[order(Inferred_Edges_Final[,1]),]
Inferred_Edges_Final <- as.data.frame(Inferred_Edges_Final)
Inferred_Edges_Final$Exists <- do.call(paste0, Inferred_Edges_Final) %in% do.call(paste0, danube_edges)
sum(Inferred_Edges_Final$Exists)/nrow(danube_edges) 

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
# rho <- seq(0.7, 0.9, by = 0.01)
rho <- 0.9
n_rho <- length(rho)
Inferred_Adj_Matrices <- lapply(rho, function(x){
  mcmapply(FUN = Infer_Adj_Matrix_glasso_residuals,
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
pdf(file = "Images/Danube/Inferred_Graph/Inferred_Graph_Glasso_Residuals.pdf", width = 10, height = 10)
plot(Graph_Final, layout = layout)
dev.off()

## Check the percentage of correctly identified edges
Inferred_Edges_Final <- Inferred_Edges_Final[which(Inferred_Edges_Final[,1] < Inferred_Edges_Final[,2]),]
Inferred_Edges_Final <- Inferred_Edges_Final[order(Inferred_Edges_Final[,1]),]
Inferred_Edges_Final <- as.data.frame(Inferred_Edges_Final)
Inferred_Edges_Final$Exists <- do.call(paste0, Inferred_Edges_Final) %in% do.call(paste0, danube_edges)
sum(Inferred_Edges_Final$Exists)/nrow(danube_edges) 

################################################################################
## Compare models
fit_graph <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
                      data = Y_Yi_large,
                      cond = 1:d,
                      MoreArgs = list(graph = Graph_Final,
                                      v = ceiling(max(Y)) + 1),
                      SIMPLIFY = FALSE,
                      mc.cores = detectCores() - 1)

sim_graph <- replicate(n = 200,
                       expr = Sim_Surface_MVAGG(n_sim = 20*n_data,
                                                q = 0.8,
                                                transforms = X_to_Y,
                                                CMEVM_fits = fit_graph),
                       simplify = FALSE)

u <- c(0.8, 0.85, 0.9)
n_u <- length(u)
eta_data <- array(NA, dim = c(d, d, length(u)))
for(i in 1:d){
  for(j in 1:d){
    tail_data <- mev::taildep(data = danube$data_clustered[,c(i,j)], u = u, plot = FALSE)
    eta_data[i,j,] <- tail_data$eta[,1]
  }
}

eta_graph <- array(NA, dim = c(d, d, length(u), 200))
for(i in 1:d){
  for(j in 1:d){
    for(k in 1:200){
      tail_graph <- mev::taildep(data = sim_graph[[k]]$Data_Margins[,c(i,j)], u = u, plot = FALSE)
      eta_graph[i,j,,k] <- tail_graph$eta[,1]
    }
  }
}

eta_graph_med <- array(NA, dim = c(d, d, length(u)))
for(i in 1:d){
  for(j in 1:d){
    eta_graph_med[i,j,] <- apply(eta_graph[i,j,,], 1, quantile, 0.5)
  }
}

eta_graph_se <- array(NA, dim = c(d, d, length(u)))
for(i in 1:d){
  for(j in 1:d){
    eta_graph_se[i,j,] <- apply(eta_graph[i,j,,], 1, sd)
  }
}

## plot the output
methods <- c("Three-step - Graphical")
n_methods <- length(methods)

elements <- lapply(apply(combinations(d, r = 2, v = 1:d), 1, list), function(x){x[[1]]})
n_elements <- length(elements)

## eta plot
eta_out_comp <- data.frame(Site_1 = rep(do.call(rbind, elements)[,1], n_methods*n_u),
                           Site_2 = rep(do.call(rbind, elements)[,2], n_methods*n_u))
eta_out_comp$Connected <- do.call(paste0, eta_out_comp[,1:2]) %in% do.call(paste0, flow_connected)
eta_out_comp$u <- rep(rep(u, each = n_elements), n_methods)
eta_out_comp$Method <- rep(methods, each = n_elements*n_u)
eta_out_comp$x <- rep(do.call(c, lapply(1:n_u, function(i){eta_data[,,i][lower_tri_flow]})), n_methods)
eta_out_comp$y <- c(do.call(c, lapply(1:n_u, function(i){eta_graph_med[,,i][lower_tri_flow]})))
eta_out_comp$se <- c(do.call(c, lapply(1:n_u, function(i){eta_graph_se[,,i][lower_tri_flow]})))

eta_out_comp$u <- factor(eta_out_comp$u, levels = u)
eta_out_comp$Method <- factor(eta_out_comp$Method, levels = methods)

label_y <- function(labels) {
  sapply(labels, function(label) {
    substitute(eta(label), list(label = label))
  })
}

ggplot(data = eta_out_comp) + geom_point(aes(x = x, y = y, col = se, shape = Connected, alpha = Connected)) +
  lims(x = c(0.6, 1), y = c(0.6, 1)) +
  labs(x = "Empirical", y = "Theoretical", shape = "Flow-connected", color = "Standard error") +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed", linewidth = 0.5) +
  scale_shape_manual(values = c(1, 2)) +
  scale_colour_gradient(low = "blue", high = "red", breaks = c(0.02, 0.03, 0.04)) +
  scale_alpha_manual(values = c(0.75, 0.4)) +  # Set transparency: 0.3 for FALSE, 1 for TRUE
  theme(legend.position = "top",
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  guides(col = guide_colorbar(title.position = "left", title.vjust = 0.9,
                              barwidth = 20, barheight = 1),  
         shape = guide_legend(title.position = "left"),
         alpha = "none") +
  facet_grid(cols = vars(u),
             labeller = labeller(u = as_labeller(label_y, default = label_parsed)))

## Look at where the bias is
eta_CMEVM_Graph_bias <- eta_graph_med[,,1] - eta_data[,,1]

x_site <- 1:d
y_site <- 1:d
sites <- expand.grid(x_site, y_site)
sites <- sites[order(sites[,1]),]
diag(flow_connections) <- 1
sites <- cbind(sites, as.logical(c(flow_connections)))

eta_CMEVM_Graph_bias_plot_data <- as.data.frame(cbind(sites, c(eta_CMEVM_Graph_bias)))
colnames(eta_CMEVM_Graph_bias_plot_data) <- c("Conditioning_Site", "Dependent_Site", "Connected", "Value")
eta_CMEVM_Graph_bias_plot_data$Conditioning_Site <- factor(eta_CMEVM_Graph_bias_plot_data$Conditioning_Site, levels = 1:d)
eta_CMEVM_Graph_bias_plot_data$Dependent_Site <- factor(eta_CMEVM_Graph_bias_plot_data$Dependent_Site, levels = 1:d)
eta_CMEVM_Graph_bias_plot_data$Value[eta_CMEVM_Graph_bias_plot_data$Conditioning_Site == eta_CMEVM_Graph_bias_plot_data$Dependent_Site] <- 0

ggplot(data = eta_CMEVM_Graph_bias_plot_data, aes(x = Conditioning_Site, y = Dependent_Site, fill = Value, col = Connected)) + 
  scale_x_discrete(breaks = seq(1, d, 1)) + 
  scale_y_discrete(breaks = seq(1, d, 1)) +
  labs(x = "Condtioning Station", y = "Dependent Station", fill = "Bias", col = "Flow-connected") +
  geom_tile(linewidth = 2) + 
  theme(panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "top",
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 0.9,
                               barwidth = 20, barheight = 1)) +  # Adjust the fill legend
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "blue"), breaks = c("FALSE", "TRUE")) +
  scale_fill_gradient2(low = "red", mid = "white", high = "gold", midpoint = 0)