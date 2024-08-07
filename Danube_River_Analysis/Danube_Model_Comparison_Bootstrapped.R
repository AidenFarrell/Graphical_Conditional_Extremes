################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("evd", "ggplot2", "ggpubr", "graphicalExtremes", "igraph", "mev")
required_pckgs <- c("graphicalExtremes", "igraph", "mev", "parallel")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
chibar_emp <- function(X, Y, u){
  if(!is.numeric(X)){
    stop("X must be a vector")
  }
  if(!is.numeric(Y)){
    stop("Y must be a vector")
  }
  n_x <- length(X)
  n_y <- length(Y)
  if(n_x != n_y){
    stop("X and Y are not the same length")
  }
  
  if(!is.numeric(u)){
    stop("u must be a vector of probabilities in interval [0,1] inclusive")
  }
  if(min(u) < 0 | max(u) > 1){
    stop("u must be a vector of probabilities in interval [0,1] inclusive")
  }
  
  U <- rank(X)/(n_x + 1)
  V <- rank(Y)/(n_y + 1)
  
  p_u <- sapply(u, function(v){length(which(U > v))})/n_x
  p_uv <- sapply(u, function(v){length(which(U > v & V > v))})/n_x
  chibar <- 2*log(p_u)/log(p_uv) - 1
  return(chibar)
}

################################################################################
## set seed for script
seed <- 361323376

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
## DO NOT RUN

## Bootstrap the data
n_boot <- 200
n_data <- nrow(danube$data_clustered)
set.seed(seed)
danube_boot <- replicate(n = n_boot, 
                         expr = danube$data_clustered[sample(x = 1:n_data, n_data, replace = TRUE),],
                         simplify = FALSE)

## Transform data onto standard Laplace margins using the Coles and Tawn Method
danube_boot_list <- lapply(danube_boot, function(x){
  lapply(apply(x, 2, list), function(x){x[[1]]})})

X_to_Y <- lapply(1:n_boot, function(i){
  mcmapply(FUN = X_to_Laplace,
           x = danube_boot_list[[i]],
           MoreArgs = list(q = seq(0.55, 0.925, by = 0.01)),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

saveRDS(X_to_Y, file = "Data/Danube_Bootstrapped.RData")
################################################################################
## Read in data
X_to_Y <- readRDS("Data/Danube_Bootstrapped.RData")

n_boot <- 200
n_data <- nrow(danube$data_clustered)

## Get the output
X <- lapply(X_to_Y, function(x){sapply(x, function(y){y$data$X})})
Y <- lapply(X_to_Y, function(x){sapply(x, function(y){y$data$Y})})

u_final <- lapply(X_to_Y, function(x){unname(sapply(x, function(y){unname(y$par$u)}))})
qu_final <- lapply(X_to_Y, function(x){do.call(c, unname(lapply(x, function(y){unname(y$par$qu)})))})
scale_final <- lapply(X_to_Y, function(x){unname(sapply(x, function(y){unname(y$par$scale)}))})
shape_final <- lapply(X_to_Y, function(x){unname(sapply(x, function(y){unname(y$par$shape)}))})

## Now we want to subset the data so that each component is large in turn
dqu <- ceiling(max(sapply(qu_final, max))/0.01)*0.01
dqu <- 0.8

Y_u <- qlaplace(dqu)
Y_Yi_large <- rep(list(vector("list", d)), n_boot)
for(i in 1:n_boot){
  for(j in 1:d){
    Y_Yi_large[[i]][[j]] <- Y[[i]][which(Y[[i]][,j] > Y_u),]
  } 
}

################################################################################
## Fit the Engelke + Hitz model to the data
## Use the variogram method for speed purposes

fit_EH_average <- mcmapply(FUN = fmpareto_graph_HR,
                           data = X,
                           MoreArgs = list(graph = danube_graph,
                                           dqu = dqu,
                                           method = "vario",
                                           handleCliques = "average"),
                           SIMPLIFY = FALSE,
                           mc.cores = detectCores() - 1)

## Fit the original Heffernan and Tawn model to the data
v <- ceiling(max(sapply(Y, max))) + 1
fit_HT <- lapply(1:n_boot, function(i){
  mcmapply(FUN = Cond_Extremes_MVN,
           data = Y_Yi_large[[i]], 
           cond = 1:d,
           MoreArgs = list(graph = NA,
                           v = v),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)
})


## Fit the conditional multivariate extremes model where we assume the residuals
## have a multivariate asymmetric generalised Gaussian distribution
## Use the three-step method for inference
fit_MVAGG_Three_Step_Graph <- lapply(1:n_boot, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
           data = Y_Yi_large[[i]], 
           cond = 1:d,
           MoreArgs = list(graph = danube_graph,
                           maxit = 1e+9,
                           v = v),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)
})

fit_MVAGG_Three_Step_Full <- lapply(1:n_boot, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
           data = Y_Yi_large[[i]], 
           cond = 1:d,
           MoreArgs = list(graph = make_full_graph(n = d),
                           maxit = 1e+9,
                           v = v),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)
})

################################################################################
## Simulate surfaces from the model

## Engelke and Hitz model
X_EH_Pareto <- lapply(fit_EH_average, function(x){
  rmpareto(n = n_data, model = "HR", par = x)})
X_EH_Uniform <- lapply(X_EH_Pareto, function(x){apply(x, 2, rank)/(1 + n_data)})
X_EH_Original_Margins <- lapply(1:n_boot, function(i){sapply(1:d, function(j){
  texmex:::revTransform(x = X_EH_Uniform[[i]][,j],
                        data = X[[i]][,j],
                        qu = qu_final[[i]][j],
                        th = u_final[[i]][j],
                        sigma = scale_final[[i]][j],
                        xi = shape_final[[i]][j],
                        method = "mixture")})})

## Heffernan and Tawn model
X_HT <- mcmapply(Sim_Surface_HT,
                 transforms = X_to_Y,
                 CMEVM_fits = fit_HT,
                 MoreArgs = list(n_sim = 20*n_data, q = dqu),
                 SIMPLIFY = FALSE,
                 mc.cores = detectCores() - 1)

## Three-step graphical model
X_Three_Step_Graph <- mcmapply(Sim_Surface_MVAGG,
                               transforms = X_to_Y,
                               CMEVM_fits = fit_MVAGG_Three_Step_Graph,
                               MoreArgs = list(n_sim = 20*n_data, q = dqu),
                               SIMPLIFY = FALSE,
                               mc.cores = detectCores() - 1)

## Three_step full model
X_Three_Step_Full <- mcmapply(Sim_Surface_MVAGG,
                              transforms = X_to_Y,
                              CMEVM_fits = fit_MVAGG_Three_Step_Full,
                              MoreArgs = list(n_sim = 20*n_data, q = dqu),
                              SIMPLIFY = FALSE,
                              mc.cores = detectCores() - 1)

################################################################################
## Obtain summary statistics from the data
## Look at eta, chi and chibar

u <- c(0.8, 0.85, 0.9)
n_u <- length(u)

eta_boot <- eta_EH <- eta_HT <- eta_CMEVM_Graph <- eta_CMEVM_Full <- array(NA, dim = c(d, d, n_u, n_boot))
eta_boot_med <- eta_EH_med <- eta_HT_med <- eta_CMEVM_Graph_med <- eta_CMEVM_Full_med <- array(NA, dim = c(d, d, n_u))
eta_boot_se <- eta_EH_se <- eta_HT_se <- eta_CMEVM_Graph_se <- eta_CMEVM_Full_se <- array(NA, dim = c(d, d, n_u))

chi_boot <- chi_EH <- chi_HT <- chi_CMEVM_Graph <- chi_CMEVM_Full <- array(NA, dim = c(d, d, n_u, n_boot))
chi_boot_med <- chi_EH_med <- chi_HT_med <- chi_CMEVM_Graph_med <- chi_CMEVM_Full_med <- array(NA, dim = c(d, d, n_u))
chi_boot_se <- chi_EH_se <- chi_HT_se <- chi_CMEVM_Graph_se <- chi_CMEVM_Full_se <- array(NA, dim = c(d, d, n_u))

chibar_boot <- chibar_EH <- chibar_HT <- chibar_CMEVM_Graph <- chibar_CMEVM_Full <- array(NA, dim = c(d, d, n_u, n_boot))
chibar_boot_med <- chibar_EH_med <- chibar_HT_med <- chibar_CMEVM_Graph_med <- chibar_CMEVM_Full_med <- array(NA, dim = c(d, d, n_u))
chibar_boot_se <- chibar_EH_se <- chibar_HT_se <- chibar_CMEVM_Graph_se <- chibar_CMEVM_Full_se <- array(NA, dim = c(d, d, n_u))
for(i in 1:d){
  for(j in 1:d){
    if(j < i){
      next()
    }
    else{
      for(k in 1:n_boot){
        ## eta and chi output
        tail_data <- taildep(data = cbind(X[[k]][,i], X[[k]][,j]), 
                             u = u, plot = FALSE)
        eta_boot[i,j,,k] <- eta_boot[j,i,,k] <- tail_EH$eta[,1]
        chi_boot[i,j,,k] <- chi_boot[j,i,,k] <- tail_EH$chi[,1]
        
        tail_EH <- taildep(data = cbind(X_EH_Original_Margins[[k]][,i], X_EH_Original_Margins[[k]][,j]),
                           u = u, plot = FALSE)
        eta_EH[i,j,,k] <- eta_EH[j,i,,k] <- tail_EH$eta[,1]
        chi_EH[i,j,,k] <- chi_EH[j,i,,k] <- tail_EH$chi[,1]
        
        tail_HT <- taildep(data = cbind(X_HT[[k]]$Data_Margins[,i], X_HT[[k]]$Data_Margins[,j]),
                           u = u, plot = FALSE)
        eta_HT[i,j,,k] <- eta_HT[j,i,,k] <- tail_HT$eta[,1]
        chi_HT[i,j,,k] <- chi_HT[j,i,,k] <- tail_HT$chi[,1]
        
        tail_CMEVM_Graph <- taildep(data = cbind(X_Three_Step_Graph[[k]]$Data_Margins[,i], 
                                                 X_Three_Step_Graph[[k]]$Data_Margins[,j]),
                                    u = u, plot = FALSE)
        eta_CMEVM_Graph[i,j,,k] <- eta_CMEVM_Graph[j,i,,k] <- tail_CMEVM_Graph$eta[,1]
        chi_CMEVM_Graph[i,j,,k] <- chi_CMEVM_Graph[j,i,,k] <- tail_CMEVM_Graph$chi[,1]
        
        tail_CMEVM_Full <- taildep(data = cbind(X_Three_Step_Full[[k]]$Data_Margins[,i], 
                                                X_Three_Step_Full[[k]]$Data_Margins[,j]),
                                   u = u, plot = FALSE)
        eta_CMEVM_Full[i,j,,k] <- eta_CMEVM_Full[j,i,,k] <- tail_CMEVM_Full$eta[,1]
        chi_CMEVM_Full[i,j,,k] <- chi_CMEVM_Full[j,i,,k] <- tail_CMEVM_Full$chi[,1]
        
        ## chibar output
        chibar_boot[i,j,,k] <- chibar_boot[j,i,,k] <- chibar_emp(X = X[[k]][,i], Y = X[[k]][,i], u = u)
        chibar_EH[i,j,,k] <- chibar_EH[j,i,,k] <- chibar_emp(X = X_EH_Original_Margins[[k]][,i], 
                                                             Y = X_EH_Original_Margins[[k]][,j], 
                                                             u = u)
        chibar_HT[i,j,,k] <- chibar_HT[j,i,,k] <- chibar_emp(X = X_HT[[k]]$Data_Margins[,i], 
                                                             Y = X_HT[[k]]$Data_Margins[,j], 
                                                             u = u)
        chibar_CMEVM_Graph[i,j,,k] <- chibar_CMEVM_Graph[j,i,,k] <- chibar_emp(X = X_Three_Step_Full[[k]]$Data_Margins[,i], 
                                                                               Y = X_Three_Step_Full[[k]]$Data_Margins[,j], 
                                                                               u = u)
        chibar_CMEVM_Full[i,j,,k] <- chibar_CMEVM_Full[j,i,,k] <- chibar_emp(X = X_Three_Step_Full[[k]]$Data_Margins[,i], 
                                                                             Y = X_Three_Step_Full[[k]]$Data_Margins[,j], 
                                                                             u = u)
        
      }
      
      ## Eta
      
      ## Get the median estimates over the 200 samples
      eta_boot_med[i,j,] <- eta_boot_med[j,i,] <- apply(eta_boot[i,j,,], 1, quantile, 0.5)
      eta_EH_med[i,j,] <- eta_EH_med[j,i,] <- apply(eta_EH[i,j,,], 1, quantile, 0.5)
      eta_HT_med[i,j,] <- eta_HT_med[j,i,] <- apply(eta_HT[i,j,,], 1, quantile, 0.5)
      eta_CMEVM_Graph_med[i,j,] <- eta_CMEVM_Graph_med[j,i,] <- apply(eta_CMEVM_Graph[i,j,,], 1, quantile, 0.5)
      eta_CMEVM_Full_med[i,j,] <- eta_CMEVM_Full_med[j,i,] <- apply(eta_CMEVM_Full[i,j,,], 1, quantile, 0.5)
      
      ## Get the standard error over the samples
      eta_boot_se[i,j,] <- eta_boot_se[j,i,] <- apply(eta_boot[i,j,,], 1, sd)/sqrt(n_boot)
      eta_EH_se[i,j,] <- eta_EH_se[j,i,] <- apply(eta_EH[i,j,,], 1, sd)/sqrt(n_boot)
      eta_HT_se[i,j,] <- eta_HT_se[j,i,] <- apply(eta_HT[i,j,,], 1, sd)/sqrt(n_boot)
      eta_CMEVM_Graph_se[i,j,] <- eta_CMEVM_Graph_se[j,i,] <- apply(eta_CMEVM_Graph[i,j,,], 1, sd)/sqrt(n_boot)
      eta_CMEVM_Full_se[i,j,] <- eta_CMEVM_Full_se[j,i,] <- apply(eta_CMEVM_Full[i,j,,], 1, sd)/sqrt(n_boot)
      
      ## Chi
      
      ## Get the median estimates over the 200 samples
      chi_boot_med[i,j,] <- chi_boot_med[j,i,] <- apply(chi_boot[i,j,,], 1, quantile, 0.5)
      chi_EH_med[i,j,] <- chi_EH_med[j,i,] <- apply(chi_EH[i,j,,], 1, quantile, 0.5)
      chi_HT_med[i,j,] <- chi_HT_med[j,i,] <- apply(chi_HT[i,j,,], 1, quantile, 0.5)
      chi_CMEVM_Graph_med[i,j,] <- chi_CMEVM_Graph_med[j,i,] <- apply(chi_CMEVM_Graph[i,j,,], 1, quantile, 0.5)
      chi_CMEVM_Full_med[i,j,] <- chi_CMEVM_Full_med[j,i,] <- apply(chi_CMEVM_Full[i,j,,], 1, quantile, 0.5)
      
      ## Get the standard error over the samples
      chi_boot_se[i,j,] <- chi_boot_se[j,i,] <- apply(chi_boot[i,j,,], 1, sd)/sqrt(n_boot)
      chi_EH_se[i,j,] <- chi_EH_se[j,i,] <- apply(chi_EH[i,j,,], 1, sd)/sqrt(n_boot)
      chi_HT_se[i,j,] <- chi_HT_se[j,i,] <- apply(chi_HT[i,j,,], 1, sd)/sqrt(n_boot)
      chi_CMEVM_Graph_se[i,j,] <- chi_CMEVM_Graph_se[j,i,] <- apply(chi_CMEVM_Graph[i,j,,], 1, sd)/sqrt(n_boot)
      chi_CMEVM_Full_se[i,j,] <- chi_CMEVM_Full_se[j,i,] <- apply(chi_CMEVM_Full[i,j,,], 1, sd)/sqrt(n_boot)
      
      ## Chibar
      
      ## Get the median estimates over the 200 samples
      chibar_boot_med[i,j,] <- chibar_boot_med[j,i,] <- apply(chibar_boot[i,j,,], 1, quantile, 0.5)
      chibar_EH_med[i,j,] <- chibar_EH_med[j,i,] <- apply(chibar_EH[i,j,,], 1, quantile, 0.5)
      chibar_HT_med[i,j,] <- chibar_HT_med[j,i,] <- apply(chibar_HT[i,j,,], 1, quantile, 0.5)
      chibar_CMEVM_Graph_med[i,j,] <- chibar_CMEVM_Graph_med[j,i,] <- apply(chibar_CMEVM_Graph[i,j,,], 1, quantile, 0.5)
      chibar_CMEVM_Full_med[i,j,] <- chibar_CMEVM_Full_med[j,i,] <- apply(chibar_CMEVM_Full[i,j,,], 1, quantile, 0.5)
      
      ## Get the standard error over the samples
      chibar_boot_se[i,j,] <- chibar_boot_se[j,i,] <- apply(chibar_boot[i,j,,], 1, sd)/sqrt(n_boot)
      chibar_EH_se[i,j,] <- chibar_EH_se[j,i,] <- apply(chibar_EH[i,j,,], 1, sd)/sqrt(n_boot)
      chibar_HT_se[i,j,] <- chibar_HT_se[j,i,] <- apply(chibar_HT[i,j,,], 1, sd)/sqrt(n_boot)
      chibar_CMEVM_Graph_se[i,j,] <- chibar_CMEVM_Graph_se[j,i,] <- apply(chibar_CMEVM_Graph[i,j,,], 1, sd)/sqrt(n_boot)
      chibar_CMEVM_Full_se[i,j,] <- chibar_CMEVM_Full_se[j,i,] <- apply(chibar_CMEVM_Full[i,j,,], 1, sd)/sqrt(n_boot)
    }
  } 
}
