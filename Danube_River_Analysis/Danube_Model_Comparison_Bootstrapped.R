################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("evd", "ggplot2", "ggpubr", "graphicalExtremes", "igraph", "mev")
required_pckgs <- c("graphicalExtremes", "igraph", "mev", "parallel")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
chibar_emp <- function(data, u){
  ## Obtain information about the data
  if(!is.matrix(data)){
    stop("data must be a matrix with n rows and 2 columns")
  }
  n_data <- nrow(data)
  d <- ncol(data)
  if(d != 2){
    stop("data must be a matrix with n rows and 2 columns")
  }
  
  if(!is.numeric(u)){
    stop("u must be a vector of probabilities in interval [0,1] inclusive")
  }
  if(min(u) < 0 | max(u) > 1){
    stop("u must be a vector of probabilities in interval [0,1] inclusive")
  }
  
  data_uniform <- apply(data, 2, rank)/(n_data + 1)
  mat_u <- lapply(u, function(x){matrix(x, nrow = n_data, ncol = 2)})
  p_uv <- sapply(mat_u, function(x){length(which(apply(data_uniform > x, 1, all)))/(n_data + 1)})
  chibar <- 2*log(1 - u)/log(p_uv) - 1
  return(chibar)
}

chi_multi <- function(data, u){
  ## Obtain information about the data
  if(!is.matrix(data)){
    stop("data must be a matrix with n rows and 2 columns")
  }
  n_data <- nrow(data)
  d <- ncol(data)
  
  if(!is.numeric(u)){
    stop("u must be a vector of probabilities in interval [0,1] inclusive")
  }
  if(min(u) < 0 | max(u) > 1){
    stop("u must be a vector of probabilities in interval [0,1] inclusive")
  }
  
  data_uniform <- apply(data, 2, rank)/(n_data + 1)
  mat_u <- lapply(u, function(x){matrix(x, nrow = n_data, ncol = d)})
  p_surv <- sapply(mat_u, function(x){length(which(apply(data_uniform > x, 1, all)))/(n_data + 1)})
  z <- p_surv/(1 - u)
  return(z)
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
flow_connections[28,29:31] <- 1
flow_connections[29,30:31] <- 1
flow_connections[30,31] <- 1

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

# ## Bootstrap the data
# n_boot <- 200
# n_data <- nrow(danube$data_clustered)
# set.seed(seed)
# danube_boot <- replicate(n = n_boot, 
#                          expr = danube$data_clustered[sample(x = 1:n_data, n_data, replace = TRUE),],
#                          simplify = FALSE)
# 
# ## Transform data onto standard Laplace margins using the Coles and Tawn Method
# danube_boot_list <- lapply(danube_boot, function(x){
#   lapply(apply(x, 2, list), function(x){x[[1]]})})
# 
# X_to_Y <- lapply(1:n_boot, function(i){
#   mcmapply(FUN = X_to_Laplace,
#            x = danube_boot_list[[i]],
#            MoreArgs = list(q = seq(0.55, 0.925, by = 0.01)),
#            SIMPLIFY = FALSE,
#            mc.cores = detectCores() - 1)})
# 
# saveRDS(X_to_Y, file = "Data/Danube_Bootstrapped.RData")
################################################################################
# ## Read in data
X_to_Y <- readRDS("Data/Danube_Bootstrapped.RData")

## Get the output
X <- lapply(X_to_Y, function(x){sapply(x, function(y){y$data$X})})
Y <- lapply(X_to_Y, function(x){sapply(x, function(y){y$data$Y})})

n_boot <- length(X)
n_data <- nrow(X[[1]])

################################################################################

## Now we want to subset the data so that each component is large in turn
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
                           MoreArgs = list(graph = danube_graph, p = dqu,
                                           method = "vario", handleCliques = "average"),
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
## The AGG parametrs may need alternative starting values
## Below will fix this

Index_Three_Step_Graph <- lapply(fit_MVAGG_Three_Step_Graph, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_Three_Step_Graph, length) > 0)){
  for(i in 1:n_boot){
    if(is_empty(Index_Three_Step_Graph[[i]])){
      next()
    }
    else{
      start_par_HT <- c(runif(1, min = 0.1, max = 0.5), runif(1, min = 0.1, max = 0.3))
      start_par_AGG <- c(runif(1, -2, 2), runif(1, 0.5, 5), runif(1, 0.5, 5), runif(1, 0.5, 3))
      
      ind <- Index_Three_Step_Graph[[i]]
      for(j in 1:length(ind)){
        fit_MVAGG_Three_Step_Graph[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large[[i]][[ind[j]]],
                                             cond = ind[j],
                                             graph = danube_graph,
                                             v = v,
                                             start_HT = start_par_HT,
                                             start_AGG = start_par_AGG,
                                             maxit = 1e+9),
              silent = TRUE)
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 100){
    stop("100 starting values attempted")
  }
  
  Index_Three_Step_Graph <- lapply(fit_MVAGG_Three_Step_Graph, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  if(all(sapply(Index_Three_Step_Graph, length) == 0)){
    print("Model Fitting Complete")
  }
}

Index_Three_Step_Full <- lapply(fit_MVAGG_Three_Step_Full, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_Three_Step_Full, length) > 0)){
  for(i in 1:n_boot){
    if(is_empty(Index_Three_Step_Full[[i]])){
      next()
    }
    else{
      start_par_HT <- c(runif(1, min = 0.1, max = 0.5), runif(1, min = 0.1, max = 0.3))
      start_par_AGG <- c(runif(1, -2, 2), runif(1, 0.5, 5), runif(1, 0.5, 5), runif(1, 0.5, 3))
      
      ind <- Index_Three_Step_Full[[i]]
      for(j in 1:length(ind)){
        fit_MVAGG_Three_Step_Full[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large[[i]][[ind[j]]],
                                             cond = ind[j],
                                             graph = make_full_graph(n = d),
                                             v = v,
                                             start_HT = start_par_HT,
                                             start_AGG = start_par_AGG,
                                             maxit = 1e+9),
              silent = TRUE)
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 100){
    stop("100 starting values attempted")
  }
  
  Index_Three_Step_Full <- lapply(fit_MVAGG_Three_Step_Full, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  if(all(sapply(Index_Three_Step_Full, length) == 0)){
    print("Model Fitting Complete")
  }
}

################################################################################
## Simulate surfaces from the model

## Engelke and Hitz model
X_EH_Pareto <- lapply(fit_EH_average, function(x){
  rmpareto(n = n_data, model = "HR", par = x)})

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
        tail_data <- taildep(data = X[[k]][,c(i,j)], u = u, plot = FALSE)
        eta_boot[i,j,,k] <- eta_boot[j,i,,k] <- tail_data$eta[,1]
        chi_boot[i,j,,k] <- chi_boot[j,i,,k] <- tail_data$chi[,1]
        
        tail_EH <- taildep(data = X_EH_Pareto[[k]][,c(i,j)], u = u, plot = FALSE)
        eta_EH[i,j,,k] <- eta_EH[j,i,,k] <- tail_EH$eta[,1]
        chi_EH[i,j,,k] <- chi_EH[j,i,,k] <- tail_EH$chi[,1]
        
        tail_HT <- taildep(data = X_HT[[k]]$Data_Margins[,c(i,j)], u = u, plot = FALSE)
        eta_HT[i,j,,k] <- eta_HT[j,i,,k] <- tail_HT$eta[,1]
        chi_HT[i,j,,k] <- chi_HT[j,i,,k] <- tail_HT$chi[,1]

        tail_CMEVM_Graph <- taildep(data = X_Three_Step_Graph[[k]]$Data_Margins[,c(i,j)], u = u, plot = FALSE)
        eta_CMEVM_Graph[i,j,,k] <- eta_CMEVM_Graph[j,i,,k] <- tail_CMEVM_Graph$eta[,1]
        chi_CMEVM_Graph[i,j,,k] <- chi_CMEVM_Graph[j,i,,k] <- tail_CMEVM_Graph$chi[,1]

        tail_CMEVM_Full <- taildep(data = X_Three_Step_Full[[k]]$Data_Margins[,c(i,j)], u = u, plot = FALSE)
        eta_CMEVM_Full[i,j,,k] <- eta_CMEVM_Full[j,i,,k] <- tail_CMEVM_Full$eta[,1]
        chi_CMEVM_Full[i,j,,k] <- chi_CMEVM_Full[j,i,,k] <- tail_CMEVM_Full$chi[,1]

        ## chibar output
        chibar_boot[i,j,,k] <- chibar_boot[j,i,,k] <- chibar_emp(data = X[[k]][,c(i,j)], u = u)
        chibar_EH[i,j,,k] <- chibar_EH[j,i,,k] <- chibar_emp(data = X_EH_Pareto[[k]][,c(i,j)], u = u)
        chibar_HT[i,j,,k] <- chibar_HT[j,i,,k] <- chibar_emp(data = X_HT[[k]]$Data_Margins[,c(i,j)], u = u)
        chibar_CMEVM_Graph[i,j,,k] <- chibar_CMEVM_Graph[j,i,,k] <- chibar_emp(data = X_Three_Step_Graph[[k]]$Data_Margins[,c(i,j)], u = u)
        chibar_CMEVM_Full[i,j,,k] <- chibar_CMEVM_Full[j,i,,k] <- chibar_emp(data = X_Three_Step_Full[[k]]$Data_Margins[,c(i,j)], u = u)
      }
      
      ## Eta
      
      ## Get the median estimates over the 200 samples
      eta_boot_med[i,j,] <- eta_boot_med[j,i,] <- apply(eta_boot[i,j,,], 1, quantile, 0.5)
      eta_EH_med[i,j,] <- eta_EH_med[j,i,] <- apply(eta_EH[i,j,,], 1, quantile, 0.5)
      eta_HT_med[i,j,] <- eta_HT_med[j,i,] <- apply(eta_HT[i,j,,], 1, quantile, 0.5)
      eta_CMEVM_Graph_med[i,j,] <- eta_CMEVM_Graph_med[j,i,] <- apply(eta_CMEVM_Graph[i,j,,], 1, quantile, 0.5)
      eta_CMEVM_Full_med[i,j,] <- eta_CMEVM_Full_med[j,i,] <- apply(eta_CMEVM_Full[i,j,,], 1, quantile, 0.5)

      # ## Get the standard error over the samples
      eta_boot_se[i,j,] <- eta_boot_se[j,i,] <- apply(eta_boot[i,j,,], 1, sd)
      eta_EH_se[i,j,] <- eta_EH_se[j,i,] <- apply(eta_EH[i,j,,], 1, sd)
      eta_HT_se[i,j,] <- eta_HT_se[j,i,] <- apply(eta_HT[i,j,,], 1, sd)
      eta_CMEVM_Graph_se[i,j,] <- eta_CMEVM_Graph_se[j,i,] <- apply(eta_CMEVM_Graph[i,j,,], 1, sd)
      eta_CMEVM_Full_se[i,j,] <- eta_CMEVM_Full_se[j,i,] <- apply(eta_CMEVM_Full[i,j,,], 1, sd)

      ## Chi

      ## Get the median estimates over the 200 samples
      chi_boot_med[i,j,] <- chi_boot_med[j,i,] <- apply(chi_boot[i,j,,], 1, quantile, 0.5)
      chi_EH_med[i,j,] <- chi_EH_med[j,i,] <- apply(chi_EH[i,j,,], 1, quantile, 0.5)
      chi_HT_med[i,j,] <- chi_HT_med[j,i,] <- apply(chi_HT[i,j,,], 1, quantile, 0.5)
      chi_CMEVM_Graph_med[i,j,] <- chi_CMEVM_Graph_med[j,i,] <- apply(chi_CMEVM_Graph[i,j,,], 1, quantile, 0.5)
      chi_CMEVM_Full_med[i,j,] <- chi_CMEVM_Full_med[j,i,] <- apply(chi_CMEVM_Full[i,j,,], 1, quantile, 0.5)

      ## Get the standard error over the samples
      chi_boot_se[i,j,] <- chi_boot_se[j,i,] <- apply(chi_boot[i,j,,], 1, sd)
      chi_EH_se[i,j,] <- chi_EH_se[j,i,] <- apply(chi_EH[i,j,,], 1, sd)
      chi_HT_se[i,j,] <- chi_HT_se[j,i,] <- apply(chi_HT[i,j,,], 1, sd)
      chi_CMEVM_Graph_se[i,j,] <- chi_CMEVM_Graph_se[j,i,] <- apply(chi_CMEVM_Graph[i,j,,], 1, sd)
      chi_CMEVM_Full_se[i,j,] <- chi_CMEVM_Full_se[j,i,] <- apply(chi_CMEVM_Full[i,j,,], 1, sd)

      ## Chibar

      ## Get the median estimates over the 200 samples
      chibar_boot_med[i,j,] <- chibar_boot_med[j,i,] <- apply(chibar_boot[i,j,,], 1, quantile, 0.5)
      chibar_EH_med[i,j,] <- chibar_EH_med[j,i,] <- apply(chibar_EH[i,j,,], 1, quantile, 0.5)
      chibar_HT_med[i,j,] <- chibar_HT_med[j,i,] <- apply(chibar_HT[i,j,,], 1, quantile, 0.5)
      chibar_CMEVM_Graph_med[i,j,] <- chibar_CMEVM_Graph_med[j,i,] <- apply(chibar_CMEVM_Graph[i,j,,], 1, quantile, 0.5)
      chibar_CMEVM_Full_med[i,j,] <- chibar_CMEVM_Full_med[j,i,] <- apply(chibar_CMEVM_Full[i,j,,], 1, quantile, 0.5)

      ## Get the standard error over the samples
      chibar_boot_se[i,j,] <- chibar_boot_se[j,i,] <- apply(chibar_boot[i,j,,], 1, sd)
      chibar_EH_se[i,j,] <- chibar_EH_se[j,i,] <- apply(chibar_EH[i,j,,], 1, sd)
      chibar_HT_se[i,j,] <- chibar_HT_se[j,i,] <- apply(chibar_HT[i,j,,], 1, sd)
      chibar_CMEVM_Graph_se[i,j,] <- chibar_CMEVM_Graph_se[j,i,] <- apply(chibar_CMEVM_Graph[i,j,,], 1, sd)
      chibar_CMEVM_Full_se[i,j,] <- chibar_CMEVM_Full_se[j,i,] <- apply(chibar_CMEVM_Full[i,j,,], 1, sd)
    }
  } 
}

## plot the output
methods <- c("Engelke & Hitz", "Three-step - Graphical", "Three-step - Saturated")
n_methods <- length(methods)

elements <- lapply(apply(combinations(d, r = 2, v = 1:d), 1, list), function(x){x[[1]]})
n_elements <- length(elements)

## eta plot
eta_out_comp <- data.frame(Site_1 = rep(do.call(rbind, elements)[,1], n_methods*n_u),
                           Site_2 = rep(do.call(rbind, elements)[,2], n_methods*n_u))
eta_out_comp$Connected <- do.call(paste0, eta_out_comp[,1:2]) %in% do.call(paste0, flow_connected)
eta_out_comp$u <- rep(rep(u, each = n_elements), n_methods)
eta_out_comp$Method <- rep(methods, each = n_elements*n_u)
eta_out_comp$x <- rep(do.call(c, lapply(1:n_u, function(i){eta_boot_med[,,i][lower_tri_flow]})), n_methods)
eta_out_comp$y <- c(do.call(c, lapply(1:n_u, function(i){eta_EH_med[,,i][lower_tri_flow]})),
                    do.call(c, lapply(1:n_u, function(i){eta_CMEVM_Graph_med[,,i][lower_tri_flow]})),
                    do.call(c, lapply(1:n_u, function(i){eta_CMEVM_Full_med[,,i][lower_tri_flow]})))

eta_out_comp$se <- c(do.call(c, lapply(1:n_u, function(i){eta_EH_se[,,i][lower_tri_flow]})),
                     do.call(c, lapply(1:n_u, function(i){eta_CMEVM_Graph_se[,,i][lower_tri_flow]})),
                     do.call(c, lapply(1:n_u, function(i){eta_CMEVM_Full_se[,,i][lower_tri_flow]})))

eta_out_comp$u <- factor(eta_out_comp$u, levels = u)
eta_out_comp$Method <- factor(eta_out_comp$Method, levels = methods)

label_y <- function(labels) {
  sapply(labels, function(label) {
    substitute(eta(label), list(label = label))
  })
}

pdf("Images/Danube/Bootstrapped_Ouput/ETA_Comp.pdf", height = 15, width = 15)
ggplot(data = eta_out_comp) + geom_point(aes(x = x, y = y, col = se, shape = Connected)) +
  lims(x = c(1/2, 1), y = c(1/2, 1)) +
  labs(x = "Empirical", y = "Theoretical", shape = "Flow-connected", color = "Standard error") +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed", linewidth = 0.5) +
  scale_shape_manual(values = c(16, 17)) +
  scale_colour_gradient(low = "blue", high = "red",
                        breaks = c(0, 0.02, 0.04, 0.06, 0.08),
                        guide = guide_colorbar(barwidth = 10, barheight = 0.5)) +
  theme(legend.position = "top") +
  guides(col = guide_colorbar(title.position = "left", title.vjust = 0.75)) +  # Adjust the fill legend
  facet_grid(cols = vars(Method), rows = vars(u),
             labeller = labeller(u = as_labeller(label_y, default = label_parsed)))
dev.off()

## chi plot
chi_out_comp <- data.frame(Site_1 = rep(do.call(rbind, elements)[,1], n_methods*n_u),
                           Site_2 = rep(do.call(rbind, elements)[,2], n_methods*n_u))
chi_out_comp$Connected <- do.call(paste0, chi_out_comp[,1:2]) %in% do.call(paste0, flow_connected)
chi_out_comp$u <- rep(rep(u, each = n_elements), n_methods)
chi_out_comp$Method <- rep(methods, each = n_elements*n_u)
chi_out_comp$x <- rep(do.call(c, lapply(1:n_u, function(i){chi_boot_med[,,i][lower_tri_flow]})), n_methods)
chi_out_comp$y <- c(do.call(c, lapply(1:n_u, function(i){chi_EH_med[,,i][lower_tri_flow]})),
                    do.call(c, lapply(1:n_u, function(i){chi_CMEVM_Graph_med[,,i][lower_tri_flow]})),
                    do.call(c, lapply(1:n_u, function(i){chi_CMEVM_Full_med[,,i][lower_tri_flow]})))

chi_out_comp$se <- c(do.call(c, lapply(1:n_u, function(i){chi_EH_se[,,i][lower_tri_flow]})),
                     do.call(c, lapply(1:n_u, function(i){chi_CMEVM_Graph_se[,,i][lower_tri_flow]})),
                     do.call(c, lapply(1:n_u, function(i){chi_CMEVM_Full_se[,,i][lower_tri_flow]})))

chi_out_comp$u <- factor(chi_out_comp$u, levels = u)
chi_out_comp$Method <- factor(chi_out_comp$Method, levels = methods)

label_y <- function(labels) {
  sapply(labels, function(label) {
    substitute(chi(label), list(label = label))
  })
}

pdf("Images/Danube/Bootstrapped_Ouput/Chi_Comp.pdf", height = 15, width = 15)
ggplot(data = chi_out_comp) + geom_point(aes(x = x, y = y, col = se, shape = Connected)) +
  lims(x = c(0, 1), y = c(0, 1)) +
  labs(x = "Empirical", y = "Theoretical", shape = "Flow-connected", color = "Standard error") +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed", linewidth = 0.5) +
  scale_shape_manual(values = c(16, 17)) +
  scale_colour_gradient(low = "blue", high = "red") +
  theme(legend.position = "top") +
  guides(col = guide_colorbar(title.position = "left", title.vjust = 0.75)) +  # Adjust the fill legend
  facet_grid(cols = vars(Method), rows = vars(u),
             labeller = labeller(u = as_labeller(label_y, default = label_parsed)))
dev.off()

## chibar
chibar_out_comp <- data.frame(Site_1 = rep(do.call(rbind, elements)[,1], n_methods*n_u),
                           Site_2 = rep(do.call(rbind, elements)[,2], n_methods*n_u))
chibar_out_comp$Connected <- do.call(paste0, chibar_out_comp[,1:2]) %in% do.call(paste0, flow_connected)
chibar_out_comp$u <- rep(rep(u, each = n_elements), n_methods)
chibar_out_comp$Method <- rep(methods, each = n_elements*n_u)
chibar_out_comp$x <- rep(do.call(c, lapply(1:n_u, function(i){chibar_boot_med[,,i][lower_tri_flow]})), n_methods)
chibar_out_comp$y <- c(do.call(c, lapply(1:n_u, function(i){chibar_EH_med[,,i][lower_tri_flow]})),
                    do.call(c, lapply(1:n_u, function(i){chibar_CMEVM_Graph_med[,,i][lower_tri_flow]})),
                    do.call(c, lapply(1:n_u, function(i){chibar_CMEVM_Full_med[,,i][lower_tri_flow]})))

chibar_out_comp$se <- c(do.call(c, lapply(1:n_u, function(i){chibar_EH_se[,,i][lower_tri_flow]})),
                     do.call(c, lapply(1:n_u, function(i){chibar_CMEVM_Graph_se[,,i][lower_tri_flow]})),
                     do.call(c, lapply(1:n_u, function(i){chibar_CMEVM_Full_se[,,i][lower_tri_flow]})))

chibar_out_comp$u <- factor(chibar_out_comp$u, levels = u)
chibar_out_comp$Method <- factor(chibar_out_comp$Method, levels = methods)

label_y <- function(labels) {
  sapply(labels, function(label) {
    substitute(bar(chi)(label), list(label = label))
  })
}

pdf("Images/Danube/Bootstrapped_Ouput/Chibar_Comp.pdf", height = 15, width = 15)
ggplot(data = chibar_out_comp) + geom_point(aes(x = x, y = y, col = se, shape = Connected)) +
  lims(x = c(0, 1), y = c(0, 1)) +
  labs(x = "Empirical", y = "Theoretical", color = "Standard error", shape = "Flow-connected") +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed", linewidth = 0.5) +
  scale_shape_manual(values = c(16, 17)) +
  scale_colour_gradient(low = "blue", high = "red",
                        breaks = c(0, 0.04, 0.08, 0.12),
                        guide = guide_colorbar(barwidth = 10, barheight = 0.5)) +
  theme(legend.position = "top") +
  guides(col = guide_colorbar(title.position = "left", title.vjust = 0.75)) +  # Adjust the fill legend
  facet_grid(cols = vars(Method), rows = vars(u),
             labeller = labeller(u = as_labeller(label_y, default = label_parsed)))
dev.off()

## Determine where the bias in the graphical model is
eta_CMEVM_Graph_bias <- eta_CMEVM_Graph_med[,,1] - eta_boot_med[,,1]

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

pdf("Images/Danube/Bootstrapped_Ouput/Bias_ETA_80_CMEVM_Graph.pdf", height = 15, width = 15)
ggplot(data = eta_CMEVM_Graph_bias_plot_data, aes(x = Conditioning_Site, y = Dependent_Site, fill = Value, col = Connected)) + 
  scale_x_discrete(breaks = seq(1, d, 1)) + 
  scale_y_discrete(breaks = seq(1, d, 1)) +
  labs(x = "Condtioning Station", y = "Dependent Station", fill = "Bias", col = "Flow-connected") +
  geom_tile(linewidth = 2) + 
  theme(panel.background = element_blank(), axis.ticks = element_blank(), legend.position = "top") +
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 0.75)) +  # Adjust the fill legend
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "blue"), breaks = c("FALSE", "TRUE")) +
  scale_fill_gradient2(low = "red", mid = "white", high = "gold", midpoint = 0)
dev.off()

################################################################################

## Look at multi chi plots
u <- c(0.8, 0.85, 0.9)
n_u <- length(u)
n_combos <- 500

## Get the trivaraite combinations
trivaraite_combos <- combinations(n = d, r = 3, v = 1:d)
r <- sort(sample(x = 1:nrow(trivaraite_combos), size = n_combos, replace = FALSE))
tri_combos <- trivaraite_combos[r,]

## Calculate the value
chi_multi_boot <- chi_multi_EH <- chi_multi_HT <- chi_multi_CMEVM_Graph <- chi_multi_CMEVM_Full <- array(NA, dim = c(n_combos, n_u, n_boot))
chi_multi_boot_med <- chi_multi_EH_med <- chi_multi_HT_med <- chi_multi_CMEVM_Graph_med <- chi_multi_CMEVM_Full_med <- matrix(NA, nrow = n_combos, ncol = n_u)
chi_multi_boot_se <- chi_multi_EH_se <- chi_multi_HT_se <- chi_multi_CMEVM_Graph_se <- chi_multi_CMEVM_Full_se <- matrix(NA, nrow = n_combos, ncol = n_u)

for(i in 1:n_combos){
  index <- tri_combos[i,]
  for(k in 1:n_boot){
    ## chi multi output
    chi_multi_boot[i,,k] <- chi_multi(data = X[[k]][,index], u = u)
    chi_multi_EH[i,,k] <- chi_multi(data = X_EH_Pareto[[k]][,index], u = u)
    chi_multi_HT[i,,k] <- chi_multi(data = X_HT[[k]]$Data_Margins[,index], u = u)
    chi_multi_CMEVM_Graph[i,,k] <- chi_multi(data = X_Three_Step_Graph[[k]]$Data_Margins[,index], u = u)
    chi_multi_CMEVM_Full[i,,k] <- chi_multi(data = X_Three_Step_Full[[k]]$Data_Margins[,index], u = u)
  }
  
  ## Get the median estimates over the 200 samples
  chi_multi_boot_med[i,] <- apply(chi_multi_boot[i,,], 1, quantile, 0.5)
  chi_multi_EH_med[i,] <- apply(chi_multi_EH[i,,], 1, quantile, 0.5)
  chi_multi_HT_med[i,] <- apply(chi_multi_HT[i,,], 1, quantile, 0.5)
  chi_multi_CMEVM_Graph_med[i,] <- apply(chi_multi_CMEVM_Graph[i,,], 1, quantile, 0.5)
  chi_multi_CMEVM_Full_med[i,] <- apply(chi_multi_CMEVM_Full[i,,], 1, quantile, 0.5)
  
  ## Get the standard error over the samples
  chi_multi_boot_se[i,] <- apply(chi_multi_boot[i,,], 1, sd)
  chi_multi_EH_se[i,] <- apply(chi_multi_EH[i,,], 1, sd)
  chi_multi_HT_se[i,] <- apply(chi_multi_HT[i,,], 1, sd)
  chi_multi_CMEVM_Graph_se[i,] <- apply(chi_multi_CMEVM_Graph[i,,], 1, sd)
  chi_multi_CMEVM_Full_se[i,] <- apply(chi_multi_CMEVM_Full[i,,], 1, sd)
}


## plot the output
methods <- c("Engelke & Hitz", "Three-step - Graphical", "Three-step - Saturated")
n_methods <- length(methods)

## eta plot
chi_multi_comp <- data.frame(x = rep(c(chi_multi_boot_med), n_methods))
chi_multi_comp$y <- c(c(chi_multi_EH_med), c(chi_multi_CMEVM_Graph_med), c(chi_multi_CMEVM_Full_med))
chi_multi_comp$se <- c(c(chi_multi_EH_se), c(chi_multi_CMEVM_Graph_se), c(chi_multi_CMEVM_Full_se))
chi_multi_comp$u <- rep(rep(u, each = n_combos), n_methods)
chi_multi_comp$Method <- rep(methods, each = n_combos*n_u)

chi_multi_comp$u <- factor(chi_multi_comp$u, levels = u)
chi_multi_comp$Method <- factor(chi_multi_comp$Method, levels = methods)

label_y <- function(labels) {
  sapply(labels, function(label) {
    substitute(chi[A](label), list(label = label))
  })
}

pdf("Images/Danube/Bootstrapped_Ouput/Chi_Multi_Comp.pdf", height = 15, width = 15)
ggplot(data = chi_multi_comp) + geom_point(aes(x = x, y = y, col = se)) +
  lims(x = c(0, 1), y = c(0, 1)) +
  labs(x = "Empirical", y = "Theoretical", color = "Standard error") +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed", linewidth = 0.5) +
  scale_shape_manual(values = c(16, 17)) +
  scale_colour_gradient(low = "blue", high = "red",
                        breaks = c(0.04, 0.06, 0.08, 0.1)) +
  theme(legend.position = "top") +
  guides(col = guide_colorbar(title.position = "left",
                              title.vjust = 1,
                              barwidth = 10,
                              barheight = 0.5,
                              direction = "horizontal"
  )) +
  facet_grid(cols = vars(Method), rows = vars(u),
             labeller = labeller(u = as_labeller(label_y, default = label_parsed)))
dev.off()




