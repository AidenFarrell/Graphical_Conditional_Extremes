################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("ggplot2", "graphicalExtremes", "igraph", "mev", "parallel")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in general functions
source("Miscellaneous_Functions/General_Functions.R")

## Reading in functions to transform the data
source("Miscellaneous_Functions/Transformations.R")

## Reading in functions required for model fitting
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")

## Read in functions for prediction
source("Prediction/Sim_Surfaces.R")

################################################################################
## Plot the upper Danube River basin
danube_graph <- graph(t(danube$flow_edges), directed = FALSE)
layout <- cbind(danube$info$PlotCoordX, danube$info$PlotCoordY)
V(danube_graph)$size <- 15
V(danube_graph)$label.cex <- 2
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot(danube_graph, layout = layout, edge.width = 10)

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

## Euclidean distances between the sites
danube_distances <- as.matrix(dist(x = cbind(danube$info$Long, danube$info$Lat),
                                   diag = TRUE, upper = TRUE)*60*1.852)

## Sites that are connected (lowest site first)
flow_connected <- which(flow_connections == 1, arr.ind = TRUE)
flow_connected <- flow_connected[which(flow_connected[,1] < flow_connected[,2]),]
flow_connected <- flow_connected[order(flow_connected[,1]),]
flow_connected <- as.data.frame(flow_connected)

################################################################################
## Transform data onto standard Laplace margins using the Coles and Tawn Method
d <- ncol(danube$data_clustered)
n_data <- nrow(danube$data_clustered)

danube_list <- lapply(apply(danube$data_clustered, 2, list), function(x){x[[1]]})
X_to_Y <- mcmapply(FUN = X_to_Laplace,
                   x = danube_list,
                   MoreArgs = list(q = seq(0.55, 0.925, by = 0.01)),
                   SIMPLIFY = FALSE,
                   mc.cores = detectCores() - 1)

## Get the output
u_final <- unname(sapply(X_to_Y, function(x){unname(x$par$u)}))
qu_final <- unname(sapply(X_to_Y, function(x){unname(x$par$qu)}))
scale_final <- unname(sapply(X_to_Y, function(x){unname(x$par$scale)}))
shape_final <- unname(sapply(X_to_Y, function(x){unname(x$par$shape)}))
Y <- unname(sapply(X_to_Y, function(x){unname(x$data$Y)}))

## Now we want to subset the data so that each component is large in turn
dqu <- ceiling(max(qu_final)/0.01)*0.01
dqu <- 0.8

Y_u <- apply(Y, 2, quantile, dqu)
Y_Yi_large <- rep(list(list()), d)
for(i in 1:d){
  Y_Yi_large[[i]] <- Y[which(Y[,i] > Y_u[i]),]
}

################################################################################
## Check the GPD fits to ensure they sufficiently capture the tail of the data
ci <- c(0.025, 0.975)
for(i in 1:d){
  n_excesses <- length(which(danube$data_clustered[,i] > u_final[i]))
  p <- (1:n_excesses)/(n_excesses + 1)
  
  excesses_station <- danube$data_clustered[which(danube$data_clustered[,i] > u_final[i]),i]
  q_empirical <- quantile(excesses_station, p)
  q_theoretical <- qgpd(p, shape_final[i], scale = scale_final[i], mu = u_final[i])
  
  ## bootstrapped empirical CI
  n_boot <- 250
  boot_data <- replicate(n = n_boot,
                         expr = danube$data_clustered[sample(x = 1:n_data, size = n_data, replace = TRUE),],
                         simplify = FALSE)
  boot_excesses_station <- lapply(boot_data, function(x){x[which(x[,i] > u_final[i]),i]})
  boot_emp_ci <- unname(apply(unname(sapply(boot_excesses_station, quantile, probs = p)), 1, quantile, probs = ci))
  
  par(mfrow = c(1, 1))
  plot(1, type = "n",
       xlab = "Empirical", ylab = "Theoretical", main = paste0("Station ", i),
       xlim = c(min(excesses_station), max(boot_emp_ci, q_theoretical) + 0.1),
       ylim = c(min(excesses_station), max(boot_emp_ci, q_theoretical) + 0.1))
  polygon(x = c(q_empirical, rev(q_empirical)), y = c(boot_emp_ci[1,], rev(boot_emp_ci[2,])), col = "grey80", border = NA)
  points(x = q_empirical, y = q_theoretical, pch = 19)
  abline(a = 0, b = 1, col = 2, lty = 2, lwd = 2) 
}
## a couple of stations do not estimate the most extreme point well but the model
## performs sufficiently well

################################################################################
## Fit the Engelke + Hitz model to the data
## Use the variogram method for speed purposes
fit_EH_average <- fmpareto_graph_HR(data = danube$data_clustered,
                                    graph = danube_graph,
                                    p = dqu,
                                    method = "vario",
                                    handleCliques = "average")

## Fit the original Heffernan and Tawn model to the data
v <- ceiling(max(Y)) + 1
fit_HT <- mcmapply(FUN = Cond_Extremes_MVN ,
                   data = Y_Yi_large, 
                   cond = 1:d,
                   MoreArgs = list(graph = NA,
                                   v = v),
                   SIMPLIFY = FALSE,
                   mc.cores = detectCores() - 1)

## Fit the conditional multivariate extremes model where we assume the residuals
## have a multivariate asymmetric generalised Gaussian distribution
## Use the three-step method for inference
fit_MVAGG_Three_Step_Graph <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
                                       data = Y_Yi_large, 
                                       cond = 1:d,
                                       MoreArgs = list(graph = danube_graph,
                                                       maxit = 1e+9,
                                                       v = v),
                                       SIMPLIFY = FALSE,
                                       mc.cores = detectCores() - 1)

fit_MVAGG_Three_Step_Full <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
                                      data = Y_Yi_large, 
                                      cond = 1:d,
                                      MoreArgs = list(graph = make_full_graph(n = d),
                                                      maxit = 1e+9,
                                                      v = v),
                                      SIMPLIFY = FALSE,
                                      mc.cores = detectCores() - 1)

################################################################################
## Get a heat map of alpha and beta parameters to assess if the Heffernan and Tawn is consistent
## Conditioning site is the column
a_hat_HT <- sapply(1:d, function(i){Add_NA_Vector(unname(unlist(fit_HT[[i]]$par$main[1,])), i)})
b_hat_HT <- sapply(1:d, function(i){Add_NA_Vector(unname(unlist(fit_HT[[i]]$par$main[2,])), i)})

x_site <- 1:d
y_site <- 1:d
sites <- expand.grid(x_site, y_site)
sites <- sites[order(sites[,1]),]

a_hat_HT_plot <- as.data.frame(cbind(sites, c(a_hat_HT)))
colnames(a_hat_HT_plot) <- c("Conditioning_Site", "Dependent_Site", "Value")
ggplot(data = a_hat_HT_plot, aes(x = Conditioning_Site, y = Dependent_Site, fill = Value)) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high="red", midpoint = 0)

b_hat_HT_plot <- as.data.frame(cbind(sites, c(b_hat_HT)))
colnames(b_hat_HT_plot) <- c("Conditioning_Site", "Dependent_Site", "Value")
ggplot(data = b_hat_HT_plot, aes(x = Conditioning_Site, y = Dependent_Site, fill = Value)) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high="red", midpoint = 0)
## Heffernan and Tawn is not consistent but it's not that bad

################################################################################
## Obtain prediction on the original scale

## Need to simulate surfaces from the models first
n_sim_surfaces <- 200

## Engelke and Hitz model
X_EH_Pareto <- replicate(n = n_sim_surfaces,
                         expr = rmpareto(n = n_data, model = "HR", par = fit_EH_average),
                         simplify = FALSE)
X_EH_Uniform <- lapply(X_EH_Pareto, function(x){apply(x, 2, rank)/(1 + n_data)})
X_EH_Original_Margins <- lapply(1:n_sim_surfaces, function(i){sapply(1:d, function(j){
  texmex:::revTransform(x = X_EH_Uniform[[i]][,j],
                        data = danube$data_clustered[,j],
                        qu = qu_final[j],
                        th = u_final[j],
                        sigma = scale_final[j],
                        xi = shape_final[j],
                        method = "mixture")})})

## Heffernan and Tawn model
X_HT <- replicate(n = n_sim_surfaces,
                  expr = Sim_Surface_HT(n_sim = 10*n_data, q = dqu,
                                        transforms = X_to_Y,
                                        CMEVM_fits = fit_HT),
                  simplify = FALSE)

## Three-step graphical model
X_Three_Step_Graph <- replicate(n = n_sim_surfaces,
                                expr = Sim_Surface_MVAGG(n_sim = 10*n_data, q = dqu,
                                                         transforms = X_to_Y,
                                                         CMEVM_fits = fit_MVAGG_Three_Step_Graph),
                                simplify = FALSE)

## Three_step full model
X_Three_Step_Full <- replicate(n = n_sim_surfaces,
                               expr = Sim_Surface_MVAGG(n_sim = 10*n_data, q = dqu,
                                                        transforms = X_to_Y,
                                                        CMEVM_fits = fit_MVAGG_Three_Step_Full),
                               simplify = FALSE)

## Now look at the eta plots
u <- c(0.8, 0.85, 0.9)
n_u <- length(u)
eta_data <- array(NA, dim = c(d, d, n_u))
for(i in 1:d){
  for(j in 1:d){
    if(j < i){
      next()
    }
    else{
      tail_data <- taildep(data = cbind(danube$data_clustered[,i], danube$data_clustered[,j]),
                           u = u, plot = FALSE)
      eta_data[i,j,] <- eta_data[j,i,] <- tail_data$eta[,1]
    }
  }
}

eta_EH <- eta_HT <- eta_CMEVM_Graph <- eta_CMEVM_Full <- array(NA, dim = c(d, d, n_u, n_sim_surfaces))
eta_EH_med <- eta_HT_med <- eta_CMEVM_Graph_med <- eta_CMEVM_Full_med <- array(NA, dim = c(d, d, n_u))
for(k in 1:n_sim_surfaces){
  for(i in 1:d){
    for(j in 1:d){
      if(j < i){
        next()
      }
      else{
        tail_EH <- taildep(data = cbind(X_EH_Original_Margins[[k]][,i], X_EH_Original_Margins[[k]][,j]),
                           u = u, plot = FALSE)
        eta_EH[i,j,,k] <- eta_EH[j,i,,k] <- tail_EH$eta[,1]
        
        tail_HT <- taildep(data = cbind(X_HT[[k]]$Data_Margins[,i], X_HT[[k]]$Data_Margins[,j]),
                           u = u, plot = FALSE)
        eta_HT[i,j,,k] <- eta_HT[j,i,,k] <- tail_HT$eta[,1]
        
        tail_CMEVM_Graph <- taildep(data = cbind(X_Three_Step_Graph[[k]]$Data_Margins[,i], 
                                                 X_Three_Step_Graph[[k]]$Data_Margins[,j]),
                                    u = u, plot = FALSE)
        eta_CMEVM_Graph[i,j,,k] <- eta_CMEVM_Graph[j,i,,k] <- tail_CMEVM_Graph$eta[,1]
        
        tail_CMEVM_Full <- taildep(data = cbind(X_Three_Step_Full[[k]]$Data_Margins[,i], 
                                                X_Three_Step_Full[[k]]$Data_Margins[,j]),
                                   u = u, plot = FALSE)
        eta_CMEVM_Full[i,j,,k] <- eta_CMEVM_Full[j,i,,k] <- tail_CMEVM_Full$eta[,1]
      }
    }
  } 
}

eta_EH_med <- eta_HT_med <- eta_CMEVM_Graph_med <- eta_CMEVM_Full_med <- array(NA, dim = c(d, d, n_u))
for(i in 1:d){
  for(j in 1:d){
    eta_EH_med[i,j,] <- apply(eta_EH[i,j,,], 1, quantile, 0.5)
    eta_HT_med[i,j,] <- apply(eta_HT[i,j,,], 1, quantile, 0.5)
    eta_CMEVM_Graph_med[i,j,] <- apply(eta_CMEVM_Graph[i,j,,], 1, quantile, 0.5)
    eta_CMEVM_Full_med[i,j,] <- apply(eta_CMEVM_Full[i,j,,], 1, quantile, 0.5)
  }
}

## plot the output
# methods <- c("Engelke & Hitz", "Heffernan and Tawn", "Three-step - Graphical", "Three-step - Saturated")
methods <- c("Engelke & Hitz", "Three-step - Graphical", "Three-step - Saturated")
elements <- lapply(apply(combinations(d, r = 2, v = 1:d), 1, list), function(x){x[[1]]})

n_methods <- length(methods)
n_elements <- length(elements)

eta_out_comp <- data.frame(Site_1 = rep(do.call(rbind, elements)[,1], n_methods*n_u),
                           Site_2 = rep(do.call(rbind, elements)[,2], n_methods*n_u))
eta_out_comp$Connected <- do.call(paste0, eta_out_comp[,1:2]) %in% do.call(paste0, flow_connected)
eta_out_comp$u <- rep(rep(u, each = n_elements), n_methods)
eta_out_comp$Method <- rep(methods, each = n_elements*n_u)
eta_out_comp$x <- rep(c(sapply(1:n_u, function(i){eta_data[,,i][lower_tri_flow]})), n_methods)
eta_out_comp$y <- c(c(sapply(1:n_u, function(i){eta_EH_med[,,i][lower_tri_flow]})),
                    # c(sapply(1:n_u, function(i){eta_HT_med[,,i][lower_tri_flow]})),
                    c(sapply(1:n_u, function(i){eta_CMEVM_Graph_med[,,i][lower_tri_flow]})),
                    c(sapply(1:n_u, function(i){eta_CMEVM_Full_med[,,i][lower_tri_flow]})))

eta_out_comp$u <- factor(eta_out_comp$u, levels = u)
eta_out_comp$Method <- factor(eta_out_comp$Method, levels = methods)

label_y <- function(labels) {
  sapply(labels, function(label) {
    substitute(eta(label), list(label = label))
  })
}

pdf("Images/Danube/Inference_On_Data/Eta_Comparison.pdf", height = 15, width = 15)
ggplot(data = eta_out_comp) + geom_point(aes(x = x, y = y, col = Connected, shape = Connected)) +
  lims(x = c(1/2, 1), y = c(1/2, 1)) +
  labs(x = "Empirical", y = "Theoretical", shape = "Flow-connected", col = "Flow-connected") +
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed", linewidth = 0.5) +
  scale_shape_manual(values = c(19, 4)) +
  scale_color_manual(values = c(alpha(colour = "black", alpha = 0.3), "blue")) +
  theme(legend.position = "top") +
  facet_grid(cols = vars(Method), rows = vars(u),
             labeller = labeller(u = as_labeller(label_y, default = label_parsed)))
dev.off()

## Plot where are the under estimates in the graphical CMEVM coming from
bias_CMEVM_Graph <- eta_CMEVM_Graph_med[,,1] - eta_data[,,1]
x_site <- 1:d
y_site <- 1:d
sites <- expand.grid(x_site, y_site)
sites <- sites[order(sites[,1]),]
sites <- cbind(sites, do.call(paste0, sites[,1:2]) %in% do.call(paste0, flow_connected))
bias_CMEVM_Graph_Plot_Data <- as.data.frame(cbind(sites, c(bias_CMEVM_Graph)))

colnames(bias_CMEVM_Graph_Plot_Data) <- c("Conditioning_Site", "Dependent_Site", "Connected", "Value")
bias_CMEVM_Graph_Plot_Data$Conditioning_Site <- factor(bias_CMEVM_Graph_Plot_Data$Conditioning_Site, levels = 1:d)
bias_CMEVM_Graph_Plot_Data$Dependent_Site <- factor(bias_CMEVM_Graph_Plot_Data$Dependent_Site, levels = 1:d)
bias_CMEVM_Graph_Plot_Data$Value[bias_CMEVM_Graph_Plot_Data$Conditioning_Site == bias_CMEVM_Graph_Plot_Data$Dependent_Site] <- 0

pdf("Images/Danube/Inference_On_Data/Bias_Graphical_CMEVM.pdf", height = 15, width = 15)
ggplot(data = bias_CMEVM_Graph_Plot_Data, aes(x = Conditioning_Site, y = Dependent_Site, fill = Value, col = Connected)) + 
  scale_x_discrete(breaks = seq(1, d, 1)) + 
  scale_y_discrete(breaks = seq(1, d, 1)) +
  labs(x = "Condtioning Station", y = "Dependent Station", fill = "Bias", col = "Flow-connected") +
  geom_tile(size = 2) + 
  theme(
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title.align = 0.5,  # Center the legend title
    legend.text.align = 0.5   # Center the legend text
  ) +
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 0.8)) +  # Adjust the fill legend
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "blue"), breaks = c("FALSE", "TRUE")) +
  scale_fill_gradient2(low = "red", mid = "white", high = "gold", midpoint = 0)
dev.off()