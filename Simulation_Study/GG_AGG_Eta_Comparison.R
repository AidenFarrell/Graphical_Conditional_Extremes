################################################################################
## Set the seed for the script
seed <- 429934197
set.seed(seed)

################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "ggplot2", "ggpubr", "gtools", "graphicalExtremes", "igraph", "kableExtra", "mev", "parallel", "rlang")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Set working directory

## Reading in general functions
source("Miscellaneous_Functions/MVAGG_Functions.R")
source("Miscellaneous_Functions/General_Functions.R")
source("Miscellaneous_Functions/Transformations.R")
source("Miscellaneous_Functions/Plotting_Functions.R")

## read in threshold selection functions
source("threshold_selection_paper/helper_functions.R")
source("threshold_selection_paper/thresh_qq_metric.R")

## Reading in functions required for model fitting
source("Model_Fitting/Cond_Extremes_MVN_Residuals.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_One_Step.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Two_Step.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")
source("Model_Fitting/Cond_Extremes_MVGG_Residuals_Three_Step.R")

## Read in functions for prediction
source("Prediction/Conditonal_Probability_Calculations.R")
source("Prediction/Sim_Surfaces.R")

################################################################################
## DO NOT RUN

## True graph
d <- 20
prop_edges <- 0.2
all_edges <- combinations(n = d, r = 2, v = 1:d)
g_true <- make_empty_graph()
while(length(V(g_true)) != d && is_connected(g_true) == FALSE){
  edges <- t(all_edges[sample(x = 1:nrow(all_edges), size = nrow(all_edges)*prop_edges, replace = FALSE),])
  g_true <- make_graph(edges = edges, directed = FALSE) 
}
g_true <- make_graph(edges = edges, directed = FALSE)
plot.igraph(g_true)

## Generate a precision matrix with this structure and get Sigma
simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(g_true)), v_sign = -1,
                           v_within = c(0.3, 0.7))
Gamma_true <- simul$omega
Sigma_true <- solve(Gamma_true)
rho_true <- cov2cor(Sigma_true)

## Mean parameter
mu_true <- runif(d, -5, 5)

## number of simulations and data points
n_sim <- 200
n_data <- 5000

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

################################################################################

## Now we want to subset the data so that each component is large in turn
dqu <- 0.9
Y_u <- qlaplace(dqu)

## level i corresponds to data set ith component in the data set given that i is large
## level j corresponds to the jth data set
Y_Yi_large <- rep(list(list()), d)
for(i in 1:d){
  for(j in 1:n_sim){
    Y_Yi_large[[i]][[j]] <- Y[[j]][which(Y[[j]][,i] > Y_u),]
  }
}

################################################################################
## Fit the Three-Step Saturated models using the MVAGG and MVGG for the residual distribution
v <- ceiling(max(sapply(Y, max))) + 1

## MVAGG model
fits_MVAGG <- lapply(1:d, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
           data = Y_Yi_large[[i]],
           MoreArgs = list(graph = make_full_graph(n = d),
                           cond = i,
                           v = v),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

## Catch and issues with model fitting
Index_MVAGG <- lapply(fits_MVAGG, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_MVAGG, length) > 0)){
  for(i in 1:n_boot){
    if(is_empty(Index_MVAGG[[i]])){
      next()
    }
    else{
      start_par_HT <- c(runif(1, min = 0.1, max = 0.5), runif(1, min = 0.1, max = 0.3))
      start_par_AGG <- c(runif(1, -2, 2), runif(1, 0.5, 5), runif(1, 0.5, 5), runif(1, 0.5, 3))
      
      ind <- Index_MVAGG[[i]]
      for(j in 1:length(ind)){
        fits_MVAGG[[i]][[ind[j]]] <-
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
  if(count >= 100){
    stop("100 starting values attempted")
  }
  
  Index_MVAGG <- lapply(fits_MVAGG, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  # if(all(sapply(Index_MVAGG, length) == 0)){
  #   print("Model Fitting Complete")
  # }
}

## MVGG model
fits_MVGG <- lapply(1:d, function(i){
  mcmapply(FUN = Cond_Extremes_MVGG_Three_Step,
           data = Y_Yi_large[[i]],
           MoreArgs = list(graph = make_full_graph(n = d),
                           cond = i,
                           v = v),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

## Catch any issues with model fitting
Index_MVGG <- lapply(fits_MVGG, function(x){which(sapply(x, is.character))})
count <- 0
while(any(sapply(Index_MVGG, length) > 0)){
  for(i in 1:n_boot){
    if(is_empty(Index_MVGG[[i]])){
      next()
    }
    else{
      start_par_HT <- c(runif(1, min = 0.1, max = 0.5), runif(1, min = 0.1, max = 0.3))
      start_par_GG <- c(runif(1, -2, 2), runif(1, 0.5, 5), runif(1, 0.5, 3))
      
      ind <- Index_MVGG[[i]]
      for(j in 1:length(ind)){
        fits_MVGG[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVGG_Three_Step(data = Y_Yi_large[[i]][[ind[j]]],
                                            cond = ind[j],
                                            graph = make_full_graph(n = d),
                                            v = v,
                                            start_HT = start_par_HT,
                                            start_GG = start_par_GG,
                                            maxit = 1e+9),
              silent = TRUE)
      }
    }
  }
  count <- count + 1
  if(count >= 100){
    stop("100 starting values attempted")
  }
  
  Index_MVGG <- lapply(fits_MVGG, function(x){which(sapply(x, is.character))})
  # if(all(sapply(Index_MVGG, length) == 0)){
  #   print("Model Fitting Complete")
  # }
}

Index_MVGG <- lapply(fits_MVGG, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_MVGG, length) > 0)){
  for(i in 1:n_boot){
    if(is_empty(Index_MVGG[[i]])){
      next()
    }
    else{
      start_par_HT <- c(runif(1, min = 0.1, max = 0.5), runif(1, min = 0.1, max = 0.3))
      start_par_GG <- c(runif(1, -2, 2), runif(1, 0.5, 5), runif(1, 0.5, 3))
      
      ind <- Index_MVGG[[i]]
      for(j in 1:length(ind)){
        fits_MVGG[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVGG_Three_Step(data = Y_Yi_large[[i]][[ind[j]]],
                                            cond = ind[j],
                                            graph = make_full_graph(n = d),
                                            v = v,
                                            start_HT = start_par_HT,
                                            start_GG = start_par_GG,
                                            maxit = 1e+9),
              silent = TRUE)
      }
    }
  }
  count <- count + 1
  if(count >= 100){
    stop("100 starting values attempted")
  }
  
  Index_MVGG <- lapply(fits_MVGG, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  # if(all(sapply(Index_MVGG, length) == 0)){
  #   print("Model Fitting Complete")
  # }
}

print("Model Fitting Complete")

################################################################################
## Simulate surfaces from the model

## MVAGG
X_MVAGG <- mcmapply(FUN = Sim_Surface_MVAGG,
                    transforms = lapply(1:n_sim, function(i){lapply(1:d, function(j){X_to_Y[[j]][[i]]})}),
                    CMEVM_fits = lapply(1:n_sim, function(i){lapply(1:d, function(j){fits_MVAGG[[j]][[i]]})}),
                    MoreArgs = list(n_sim = 5*n_data, q = dqu),
                    SIMPLIFY = FALSE,
                    mc.cores = detectCores() - 2)

print("MVAGG Prediction Complete")

## MVGG
X_MVGG <- mcmapply(FUN = Sim_Surface_MVGG,
                   transforms = lapply(1:n_sim, function(i){lapply(1:d, function(j){X_to_Y[[j]][[i]]})}),
                   CMEVM_fits = lapply(1:n_sim, function(i){lapply(1:d, function(j){fits_MVGG[[j]][[i]]})}),
                   MoreArgs = list(n_sim = 5*n_data, q = dqu),
                   SIMPLIFY = FALSE,
                   mc.cores = detectCores() - 2)

print("MVGG Prediction Complete")

################################################################################
## Obtain summary statistics from the data
## Look at eta

u <- c(0.95, 0.99)
n_u <- length(u)

grid <- combinations(d, r = 2, v = 1:d)
grid_list <- lapply(apply(grid, 1, list), function(x){x[[1]]})
n_grid <- length(grid_list)

eta_boot <- eta_MVAGG_boot <- eta_MVGG_boot <- array(NA, dim = c(n_grid, n_u, n_sim))
eta_boot_med <- eta_MVAGG_med <- eta_MVGG_med <- array(NA, dim = c(n_grid, n_u))
eta_boot_se <- eta_MVAGG_se <- eta_MVGG_se <- array(NA, dim = c(n_grid, n_u))

eta_data <- lapply(grid_list, function(i){
  mcmapply(FUN = taildep,
           data = lapply(X, function(x){x[,i]}),
           MoreArgs = list(u = u, 
                           plot = FALSE),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

eta_MVAGG <- lapply(grid_list, function(i){
  mcmapply(FUN = taildep,
           data = lapply(X_MVAGG, function(x){x$Data_Margins[,i]}),
           MoreArgs = list(u = u, 
                           plot = FALSE),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

eta_MVGG <- lapply(grid_list, function(i){
  mcmapply(FUN = taildep,
           data = lapply(X_MVGG, function(x){x$Data_Margins[,i]}),
           MoreArgs = list(u = u, 
                           plot = FALSE),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

## Extract the eta values
for(k in 1:n_grid){
  eta_boot[k,,] <- sapply(eta_data[[k]], function(x){x$eta[,1]})
  eta_MVAGG_boot[k,,] <- sapply(eta_MVAGG[[k]], function(x){x$eta[,1]})
  eta_MVGG_boot[k,,] <- sapply(eta_MVGG[[k]], function(x){x$eta[,1]})
}

## Obtain the mean and sd over the samples
for(i in 1:n_grid){
  for(j in 1:n_u){
    
    ## Data
    eta_boot_med[i,j] <- mean(eta_boot[i,j,])
    eta_boot_se[i,j] <- sd(eta_boot[i,j,])
    
    ## MVAGG
    eta_MVAGG_med[i,j] <- mean(eta_MVAGG_boot[i,j,])
    eta_MVAGG_se[i,j] <- sd(eta_MVAGG_boot[i,j,])
    
    ## MVGG
    eta_MVGG_med[i,j] <- mean(eta_MVGG_boot[i,j,])
    eta_MVGG_se[i,j] <- sd(eta_MVGG_boot[i,j,])
  }
}

## plot the output

## Mehtods
methods <- c("MVAGG", "MVGG")
n_methods <- length(methods)

## Data Frame with the output
eta_out_comp <- data.frame(Site_1 = rep(grid[,1], n_methods*n_u),
                           Site_2 = rep(grid[,2], n_methods*n_u))
eta_out_comp$u <- rep(rep(u, each = n_grid), n_methods)
eta_out_comp$Method <- rep(methods, each = n_grid*n_u)
eta_out_comp$x <- rep(c(eta_boot_med), n_methods)
eta_out_comp$y <- c(c(eta_MVAGG_med), c(eta_MVGG_med))
eta_out_comp$se <- c(c(eta_MVAGG_se), c(eta_MVGG_se))

eta_out_comp$u <- factor(eta_out_comp$u, levels = u)
eta_out_comp$Method <- factor(eta_out_comp$Method, levels = methods)

## Labels
label_y <- function(labels) {
  sapply(labels, function(label) {
    substitute(eta(label), list(label = label))
  })
}

## Plotting paramters
lim_min <- floor(min(eta_out_comp$x, eta_out_comp$y)*10)/10
lim_max <- ceiling(max(eta_out_comp$x, eta_out_comp$y)*10)/10
col_breaks <- seq(from = 0, plyr::round_any(max(eta_out_comp$se), 0.05), by = 0.05)

## Make the plot

pdf("Images/Simulation_Study/GG_AGG_ETA_Comp_MVN.pdf", height = 15, width = 15)
ggplot(data = eta_out_comp) + geom_point(aes(x = x, y = y, col = se)) +
  lims(x = c(lim_min, lim_max), y = c(lim_min, lim_max)) +
  labs(x = "Empirical", y = "Theoretical", color = "Standard error") +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed", linewidth = 0.5) +
  scale_colour_gradient(low = "blue", high = "red",
                        breaks = col_breaks) +
  theme(legend.position = "top",
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  guides(col = guide_colorbar(title.position = "left", title.vjust = 0.9,
                              barwidth = 15),  
         shape = guide_legend(title.position = "left"),
         alpha = "none") +
  facet_grid(cols = vars(Method), rows = vars(u),
             labeller = labeller(u = as_labeller(label_y, default = label_parsed)))
dev.off()

print("Eta Plotting Complete")

################################################################################
