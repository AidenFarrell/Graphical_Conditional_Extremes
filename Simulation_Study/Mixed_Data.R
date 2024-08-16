################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "ggplot2", "graphicalExtremes", "igraph", "parallel")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Algorithm to generate data with mixed dependence structures
mixed_data_generation <- function(n_data, Gamma_MVP, Sigma_MVN){
  
  ## Generate the inner star graph using a MVP distribution
  data_Pareto <- rmpareto(n = n_data, model = "HR", par = Gamma_MVP)
  
  ## Transform the data onto standard Gaussian margins
  data_Gaussian <- qnorm(apply(data_Pareto, 2, rank)/(n_data + 1))
  
  ## Get conditional mean and correlation matrix given X3
  mean_cond <- sapply(Sigma_MVN[1,-1], function(x){x*data_Gaussian[,3]})
  Sigma_cond <- Cond_Sigma(Sigma_MVN, 1)
  rho_cond <- cov2cor(Sigma_cond)
  
  ## Obtain remaining data conditional on X3
  data_MVN <- matrix(NA, nrow = n_data, ncol = 2)
  for(i in 1:n_data){ 
    data_MVN[i,] <- mvtnorm::rmvnorm(n = 1, mean = mean_cond[i,], sigma = rho_cond)
  }
  
  ## return all data
  data_all <- cbind(data_Gaussian, data_MVN)
  return(data_all)
}

## plotting functions for later
boxplot_MLEs <- function(data, methods, y_lab){
  
  ## Extract some information from the data
  d_data <- length(data)
  n_methods <- length(methods)
  p <- ncol(data[[1]])/n_methods
  n <- nrow(data[[1]])
  
  ## get some plotting parameters
  y_min <- floor(min(sapply(data, min, na.rm = TRUE))/0.1)*0.1
  y_max <- ceiling(max(sapply(data, max, na.rm = TRUE))/0.1)*0.1
  
  ## main data to plot
  plot_data <- data.frame(y = do.call(c, data))
  plot_data$Method = rep(rep(rep(methods, each = n), p), d_data)
  plot_data$Method <- factor(plot_data$Method, levels = methods)
  plot_data$Conditioning_Varaible <- rep(1:d_data, each = p*n_sim*n_methods)
  plot_data$Conditioning_Varaible <- factor(plot_data$Conditioning_Varaible, levels = 1:d_data)
  plot_data$Dependent_Variable <- rep(rep(1:p, each = n_sim*n_methods), d_data)
  plot_data$Dependent_Variable <- factor(plot_data$Dependent_Variable, levels = 1:p)
  
  ## Define custom labels for facets
  facet_labels <- setNames(paste0("i = ", 1:d_data), 1:d_data)
  
  ## plot the data
  plot_out <- ggplot(data = plot_data, aes(x = Dependent_Variable, y = y, fill = Method)) + 
    geom_boxplot() +
    theme(legend.position = "top") +
    labs(x = "Depednent Varaible (j)", y = y_lab) +
    facet_grid(cols = vars(Conditioning_Varaible), labeller = labeller(Conditioning_Varaible = facet_labels))
  print(plot_out)
}


boxplot_MLEs_Cov_Mat <- function(data, methods, y_lab){
  
  ## Obtain some information from the data
  if(!is.list(data)){
    stop("data must be a list of parameter estimates to plot")
  }
  d_data <- length(data)
  n_methods = length(methods)
  p <- ncol(data[[1]][[1]])
  n <- nrow(data[[1]][[1]])
  
  ## get some plotting parameters
  x_labels <- apply(combinations(n = d_data, r = 2, v = 1:d_data, repeats.allowed = TRUE), 1, paste0, collapse = "")
  
  ## Construct the plotting data
  plot_data <- data.frame(y = do.call(c, lapply(1:d_data, function(i){
    do.call(c, lapply(1:p, function(j){
      do.call(c, lapply(1:n_methods, function(k){
        data[[i]][[k]][,j]}))}))})))
  
  plot_data$Method = rep(rep(rep(methods, each = n), p), d_data)
  plot_data$Method <- factor(plot_data$Method, levels = methods)
  plot_data$Conditioning_Varaible <- rep(1:d_data, each = p*n*n_methods)
  plot_data$Conditioning_Varaible <- factor(plot_data$Conditioning_Varaible, levels = 1:d_data)
  plot_data$Pair <- rep(rep(x_labels, each = n*n_methods), d_data)
  plot_data$Pair <- factor(plot_data$Pair, levels = x_labels)
  
  ## Define custom labels for facets
  facet_labels <- setNames(paste0("i = ", 1:d_data), 1:d_data)
  
  plot_out <- ggplot(data = plot_data, aes(x = Pair, y = y, col = Method)) + 
    geom_boxplot(outlier.shape = NA) +
    theme(legend.position = "top") +
    labs(x = "Pair", y = y_lab) +
    facet_grid(rows = vars(Conditioning_Varaible), labeller = labeller(Conditioning_Varaible = facet_labels),
               scales = "free_y") +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed", linewidth = 0.25)
  
  return(plot_out)
}

################################################################################
## Set working directory

## Reading in general functions
source("Miscellaneous_Functions/MVAGG_Functions.R")
source("Miscellaneous_Functions/General_Functions.R")
source("Miscellaneous_Functions/Transformations.R")

## Reading in functions required for model fitting
source("Model_Fitting/Cond_Extremes_MVN_Residuals.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_One_Step.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Two_Step.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")

## Read in functions for prediction
source("Prediction/Conditonal_Probability_Calculations.R")
source("Prediction/Sim_Surfaces.R")

################################################################################
## DO NOT RUN

## Set the seed
seed <- 968562880
set.seed(seed)

## AD clique graph
graph_1 <- make_full_graph(n = 3)
plot(graph_1)

## Full graph
edges_full <- rbind(as_edgelist(graph_1), c(3,4), c(3,5), c(4,5))
graph_full <- graph_from_edgelist(edges_full, directed = FALSE)      
plot(graph_full)

## Simulate the inner graph data
d <- length(V(graph_full))
n_data <- 5000

## Get the precision matrix for the inner graph
p_MVP <- length(E(graph_1))
Gamma_MVP <- complete_Gamma(runif(p_MVP, min = 0.5, max = 1), graph = graph_1)

## Get the correlations for the outer edges
simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(graph_1)), v_sign = -1,
                           v_within = c(0.999, 0.99999))
Gamma_MVN <- simul$omega
Sigma_MVN <- solve(Gamma_MVN)
rho_MVN <- cov2cor(Sigma_MVN)

## Generate the data
n_sim = 200
X <- replicate(n = n_sim,
               expr = mixed_data_generation(n_data = n_data,
                                            Gamma_MVP = Gamma_MVP,
                                            Sigma_MVN = Sigma_MVN),
               simplify = FALSE)

## Transform the data onto standard Laplace margins
X_list_by_data <- lapply(X, function(x){lapply(apply(x, 2, list), function(y){y[[1]]})})
X_list_by_var <- lapply(1:d, function(i){lapply(1:n_sim, function(j){X_list_by_data[[j]][[i]]})})
X_to_Y <- lapply(X_list_by_var, function(x){
  mcmapply(x = x,
           FUN = X_to_Laplace,
           MoreArgs = list(q = seq(0.55, 0.99, by = 0.01)),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

## Get the output from the GPD fits
u_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$par$u)})})
qu_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$par$qu)})})
scale_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$par$scale)})})
shape_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$par$shape)})})
Y <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(X_to_Y[[j]][[i]]$data$Y)})})

################################################################################
## Read in data
out <- readRDS("Data/Mixed_Data_D5.RData")

## Obtain the true parameters
graph_full <- out$par_true$graph
d <- length(V(graph_full))
n_sim <- out$par_true$n_sim
n_data <- out$par_true$n_data

Gamma_MVP <- out$par_true$Gamma_MVP
Gamma_MVN <- out$par_true$Gamma_MVN

## Transforms
X_to_Y <- out$transforms

## Obtain the data
X <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$data$X)})})
Y <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$data$Y)})})

## Get the output from the GPD fits
u_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$u)})})
qu_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$qu)})})
scale_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$scale)})})
shape_final <- lapply(1:n_sim, function(i){sapply(1:d, function(j){unname(out$transforms[[j]][[i]]$par$shape)})})

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
## Fit the various models

## Fit the Engelke and Hitz model to the data for comparative purposes
fit_EH <- lapply(X, function(x){
  fmpareto_graph_HR(data = x, graph = graph_full, p = dqu,
                    method = "vario", handleCliques = "average")})

## Heffernan and Tawn model
v <- ceiling(max(sapply(Y, max))) + 1
fit_HT <- lapply(1:d, function(i){
  mcmapply(FUN = Cond_Extremes_MVN,
           data = Y_Yi_large[[i]],
           MoreArgs = list(graph = NA,
                           cond = i,
                           v = v,
                           maxit = 1e+9),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

Z_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){y$Z})})
a_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){unname(y$par$main[1,])})})
b_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){unname(y$par$main[2,])})})

## One-step model fits
## One-step model is sensitive to the starting parameter of the the location parameter in the AGG
## This is heavily correlated with the alpha dependence parameter
## Give an "informed" start for the location parameter in the AGG
loc_start_One_Step <- lapply(1:d, function(j){lapply(1:n_sim, function(k){
  apply((Y_Yi_large[[j]][[k]][,-j] - matrix(0.1*Y_Yi_large[[j]][[k]][,j], nrow = nrow(Y_Yi_large[[j]][[k]]), ncol = d-1))/
          matrix(Y_Yi_large[[j]][[k]][,j]^0.1, nrow = nrow(Y_Yi_large[[j]][[k]]), ncol = d-1), 2, mean)
})})
start_par_One_Step <- lapply(1:d, function(j){lapply(1:n_sim, function(k){
  cbind(rep(0.1, d-1), rep(0.1, d-1), loc_start_One_Step[[j]][[k]],
        rep(1.5, d-1), rep(2, d-1), rep(1.5, d-1))})})

## Graphical
fit_One_Step_Graph <- lapply(1:d, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG,
           data = Y_Yi_large[[i]],
           start = start_par_One_Step[[i]],
           MoreArgs = list(graph = graph_full,
                           cond = i,
                           maxit = 1e+9,
                           nOptim = 5),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

## Graphical
fit_Two_Step_Graph <- lapply(1:d, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG_Two_Step,
           data = Y_Yi_large[[i]],
           MoreArgs = list(graph = graph_full,
                           cond = i,
                           v = v,
                           maxit = 1e+9,
                           nOptim = 2),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

## Three-step model fits
## Independence
fit_Three_Step_Indep <- lapply(1:d, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
           data = Y_Yi_large[[i]],
           MoreArgs = list(graph = NA,
                           cond = i,
                           v = v,
                           maxit = 1e+9),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

## Graphical
fit_Three_Step_Graph <- lapply(1:d, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
           data = Y_Yi_large[[i]],
           MoreArgs = list(graph = graph_full,
                           cond = i,
                           v = v,
                           maxit = 1e+9),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

## Saturated
fit_Three_Step_Full <- lapply(1:d, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
           data = Y_Yi_large[[i]],
           MoreArgs = list(graph = make_full_graph(n = d),
                           cond = i,
                           v = v,
                           maxit = 1e+9),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

################################################################################
## One and two step models can be somewhat sensitive to the starting values
## here we need to do another procedure to ensure the models converge

Index_One_Step_Graph <- lapply(fit_One_Step_Graph, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_One_Step_Graph, length) > 0)){
  for(i in 1:d){
    if(is_empty(Index_One_Step_Graph[[i]])){
      next()
    }
    else{
      start_par_One_Step <- cbind(runif(d-1, 0.1, 0.5),
                                  runif(d-1, 0.1, 0.5),
                                  runif(d-1, -2, 2),
                                  runif(d-1, 0.5, 3),
                                  runif(d-1, 0.5, 3),
                                  runif(d-1, 0.5, 3))
      
      ind <- Index_One_Step_Graph[[i]]
      for(j in 1:length(ind)){
        fit_One_Step_Graph[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVAGG(data = Y_Yi_large[[i]][[ind[j]]],
                                  cond = i,
                                  graph = graph_full,
                                  start = start_par_One_Step,
                                  maxit = 1e+9,
                                  nOptim = 5),
              silent = TRUE)
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_One_Step_Graph <- lapply(fit_One_Step_Graph, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  if(all(sapply(Index_One_Step_Graph, length) == 0)){
    print("Model Fitting Complete")
  }
}

Index_Two_Step_Graph <- lapply(fit_Two_Step_Graph, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_Two_Step_Graph, length) > 0)){
  for(i in 1:d){
    if(is_empty(Index_Two_Step_Graph[[i]])){
      next()
    }
    else{
      start_par_Two_Step <- cbind(runif(d-1, -0.5, 0.5),
                                  runif(d-1, 0.5, 5),
                                  runif(d-1, 0.5, 5),
                                  runif(d-1, 0.5, 3))
      
      ind <- Index_Two_Step_Graph[[i]]
      for(j in 1:length(ind)){
        fit_Two_Step_Graph[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVAGG_Two_Step(data = Y_Yi_large[[i]][[ind[j]]],
                                           cond = i,
                                           graph = graph_full,
                                           v = v,
                                           start_AGG = start_par_Two_Step,
                                           maxit = 1e+9),
              silent = TRUE)
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_Two_Step_Graph <- lapply(fit_Two_Step_Graph, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  if(all(sapply(Index_Two_Step_Graph, length) == 0)){
    print("Model Fitting Complete")
  }
}

## Might have some issues with the three-step model due to the differing
## dependence structures
## Code may need a slight edit to update ranges of sampling for starting parameters
Index_Three_Step_Indep <- lapply(fit_Three_Step_Indep, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_Three_Step_Indep, length) > 0)){
  for(i in 1:d){
    if(is_empty(Index_Three_Step_Indep[[i]])){
      next()
    }
    else{
      start_par_HT <- c(runif(2, min = 0.1, max = 0.3))
      start_par_AGG <- c(runif(1, -2, 2), runif(1, 0.5, 5), runif(1, 0.5, 5), runif(1, 0.5, 3))
      
      ind <- Index_Three_Step_Indep[[i]]
      for(j in 1:length(ind)){
        fit_Three_Step_Indep[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large[[i]][[ind[j]]],
                                             cond = i,
                                             graph = NA,
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
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_Three_Step_Indep <- lapply(fit_Three_Step_Indep, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  if(all(sapply(Index_Three_Step_Indep, length) == 0)){
    print("Model Fitting Complete")
  }
}

Index_Three_Step_Graph <- lapply(fit_Three_Step_Graph, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_Three_Step_Graph, length) > 0)){
  for(i in 1:d){
    if(is_empty(Index_Three_Step_Graph[[i]])){
      next()
    }
    else{
      start_par_HT <- c(runif(2, min = 0.1, max = 0.3))
      start_par_AGG <- c(runif(1, -2, 2), runif(1, 0.5, 5), runif(1, 0.5, 5), runif(1, 0.5, 3))
      
      ind <- Index_Three_Step_Graph[[i]]
      for(j in 1:length(ind)){
        fit_Three_Step_Graph[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large[[i]][[ind[j]]],
                                             cond = i,
                                             graph = g_true,
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
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_Three_Step_Graph <- lapply(fit_Three_Step_Graph, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  if(all(sapply(Index_Three_Step_Graph, length) == 0)){
    print("Model Fitting Complete")
  }
}

Index_Three_Step_Full <- lapply(fit_Three_Step_Full, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
count <- 0
while(any(sapply(Index_Three_Step_Full, length) > 0)){
  for(i in 1:d){
    if(is_empty(Index_Three_Step_Full[[i]])){
      next()
    }
    else{
      start_par_HT <- c(runif(2, min = 0.1, max = 0.3))
      start_par_AGG <- c(runif(1, -2, 2), runif(1, 0.5, 5), runif(1, 0.5, 5), runif(1, 0.5, 3))
      
      ind <- Index_Three_Step_Full[[i]]
      for(j in 1:length(ind)){
        fit_Three_Step_Full[[i]][[ind[j]]] <-
          try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large[[i]][[ind[j]]],
                                             cond = i,
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
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_Three_Step_Full <- lapply(fit_Three_Step_Full, function(x){which(sapply(x, function(y){is.na(y$convergence)}))})
  if(all(sapply(Index_Three_Step_Full, length) == 0)){
    print("Model Fitting Complete")
  }
}

################################################################################
## Extract the parameter estimates
mat <- matrix(NA, d, d)
lower_tri_elements <- lower.tri(mat, diag = TRUE)

## One-step Graphical residuals
a_hat_One_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})
b_hat_One_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})
loc_hat_One_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})
scale_1_hat_One_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})
scale_2_hat_One_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})
shape_hat_One_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})
Gamma_hat_One_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})

## Two-step Graphical residuals
a_hat_Two_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})
b_hat_Two_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})
loc_hat_Two_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})
scale_1_hat_Two_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})
scale_2_hat_Two_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})
shape_hat_Two_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})
Gamma_hat_Two_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})

## Three-Step Independent Residuals
a_hat_Three_Step_Indep <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})
b_hat_Three_Step_Indep <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})
loc_hat_Three_Step_Indep <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})
scale_1_hat_Three_Step_Indep <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})
scale_2_hat_Three_Step_Indep <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})
shape_hat_Three_Step_Indep <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})
Gamma_hat_Three_Step_Indep <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})

## Three-Step Graphical residuals
a_hat_Three_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})
b_hat_Three_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})
loc_hat_Three_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})
scale_1_hat_Three_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})
scale_2_hat_Three_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})
shape_hat_Three_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})
Gamma_hat_Three_Step_Graph <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})

## Three-Step Saturated Residuals
a_hat_Three_Step_Full <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})
b_hat_Three_Step_Full <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})
loc_hat_Three_Step_Full <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})
scale_1_hat_Three_Step_Full <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})
scale_2_hat_Three_Step_Full <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})
shape_hat_Three_Step_Full <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})
Gamma_hat_Three_Step_Full <- lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})

## Engelke and Hitz model
Gamma_hat_EH <- lapply(1:d, function(j){
  t(sapply(fit_EH, function(y){Add_NA_Matrix(Gamma2Theta(y, k = j), j)[lower_tri_elements]}))})

################################################################################
## Assess convergence of the parameters
method_vec <- c("One-step - Graphical", "Two-step - Graphical", "Three-Step")

# Alpha
pdf(file = "Images/Simulation_Study/Mixed_Data/Alpha.pdf", width = 15, height = 10)
boxplot_MLEs(
  data = lapply(1:d, function(i){do.call(cbind, lapply(1:d, function(j){
    cbind(a_hat_One_Step_Graph[[i]][,j],
          a_hat_Two_Step_Graph[[i]][,j],
          a_hat_Three_Step_Indep[[i]][,j])}))}),
  methods = method_vec, y_lab = expression(hat(alpha)[j ~ "|" ~ i]))
dev.off()

# Beta
pdf(file = "Images/Simulation_Study/Mixed_Data/Beta.pdf", width = 15, height = 10)
boxplot_MLEs(
  data = lapply(1:d, function(i){do.call(cbind, lapply(1:d, function(j){
    cbind(b_hat_One_Step_Graph[[i]][,j],
          b_hat_Two_Step_Graph[[i]][,j],
          b_hat_Three_Step_Indep[[i]][,j])}))}),
  methods = method_vec, y_lab = expression(hat(beta)[j ~ "|" ~ i]))
dev.off()

# Location
pdf(file = "Images/Simulation_Study/Mixed_Data/Location.pdf", width = 15, height = 10)
boxplot_MLEs(
  data = lapply(1:d, function(i){do.call(cbind, lapply(1:d, function(j){
    cbind(loc_hat_One_Step_Graph[[i]][,j],
          loc_hat_Two_Step_Graph[[i]][,j],
          loc_hat_Three_Step_Indep[[i]][,j])}))}),
  methods = method_vec, y_lab = expression(hat(nu)[j ~ "|" ~ i]))
dev.off()

# Scale (Left)
pdf(file = "Images/Simulation_Study/Mixed_Data/Scale_Left.pdf", width = 15, height = 10)
boxplot_MLEs(
  data = lapply(1:d, function(i){do.call(cbind, lapply(1:d, function(j){
    cbind(scale_1_hat_One_Step_Graph[[i]][,j],
          scale_1_hat_Two_Step_Graph[[i]][,j],
          scale_1_hat_Three_Step_Indep[[i]][,j])}))}),
  methods = method_vec, y_lab = expression(hat(kappa[1])[j ~ "|" ~ i]))
dev.off()

# Scale (Right)
pdf(file = "Images/Simulation_Study/Mixed_Data/Scale_Right.pdf", width = 15, height = 10)
boxplot_MLEs(
  data = lapply(1:d, function(i){do.call(cbind, lapply(1:d, function(j){
    cbind(scale_2_hat_One_Step_Graph[[i]][,j],
          scale_2_hat_Two_Step_Graph[[i]][,j],
          scale_2_hat_Three_Step_Indep[[i]][,j])}))}),
  methods = method_vec, y_lab = expression(hat(kappa[2])[j ~ "|" ~ i]))
dev.off()

## check that the left and right scale are not the same
scale_plot_data <- data.frame(scale_1 = do.call(c, lapply(1:d, function(i){do.call(c, lapply(1:d, function(j){scale_1_hat_Three_Step_Graph[[i]][,j]}))})),
                              scale_2 = do.call(c, lapply(1:d, function(i){do.call(c, lapply(1:d, function(j){scale_2_hat_Three_Step_Graph[[i]][,j]}))})),
                              Conditioning_Variable = rep(1:d, each = n_sim*d),
                              Dependent_Variable = rep(rep(1:d, each = n_sim), d))

scale_plot_data$Conditioning_Variable <- factor(scale_plot_data$Conditioning_Variable, levels = 1:d, labels = sapply(1:d, function(k){substitute(i ~ "=" ~ j, list(j = k))}))
scale_plot_data$Dependent_Variable <- factor(scale_plot_data$Dependent_Variable, levels = 1:d, labels = sapply(1:d, function(k){substitute(Y[j] ~ "|" ~ Y[i] > u[Y[i]], list(j = k))}))

pdf(file = "Images/Simulation_Study/Mixed_Data/Scale_Comp.pdf", width = 15, height = 15)
ggplot(data = scale_plot_data, aes(x = scale_1, y = scale_2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed", linewidth = 1) +
  labs(x = substitute(hat(kappa[1])), y = substitute(hat(kappa[2]))) + 
  lims(x = c(floor(min(scale_plot_data$scale_1, scale_plot_data$scale_2, na.rm = TRUE)/0.1)*0.1,
             ceiling(max(scale_plot_data$scale_1, scale_plot_data$scale_2, na.rm = TRUE)/0.1)*0.1),
       y = c(floor(min(scale_plot_data$scale_1, scale_plot_data$scale_2, na.rm = TRUE)/0.1)*0.1,
             ceiling(max(scale_plot_data$scale_1, scale_plot_data$scale_2, na.rm = TRUE)/0.1)*0.1)) +
  facet_grid(Conditioning_Variable ~ Dependent_Variable, labeller = label_parsed)
dev.off()

# Shape
pdf(file = "Images/Simulation_Study/Mixed_Data/Shape.pdf", width = 15, height = 10)
boxplot_MLEs(
  data = lapply(1:d, function(i){do.call(cbind, lapply(1:d, function(j){
    cbind(shape_hat_One_Step_Graph[[i]][,j],
          shape_hat_Two_Step_Graph[[i]][,j],
          shape_hat_Three_Step_Indep[[i]][,j])}))}),
  methods = method_vec, y_lab = expression(hat(delta)[j ~ "|" ~ i]))
dev.off()

## Precision matrix
method_vec <- c("One-step - Graphical", "Two-step - Graphical",
                "Three-step - Graphical", "Three-step - Saturated")

pdf(file = "Images/Simulation_Study/Mixed_Data/Gamma.pdf", width = 15, height = 10)
boxplot_MLEs_Cov_Mat(
  data = lapply(1:d, function(i){
    list(Gamma_hat_One_Step_Graph[[i]],
         Gamma_hat_Two_Step_Graph[[i]],
         Gamma_hat_Three_Step_Graph[[i]],
         Gamma_hat_Three_Step_Full[[i]])}),
  methods = method_vec, y_lab = expression(hat(Gamma)[~ "|" ~ i]))
dev.off()
################################################################################
## Prediction

## "true" probabilities from the underlying probability distribution which are calculated
## using a point estimate from simulating a large number of data points from the
## true distribution
X_prob_calc <- mixed_data_generation(n_data = 1e+7,
                                     Gamma_MVP = Gamma_MVP,
                                     Sigma_MVN = Sigma_MVN)

## Get the probabilities from the Engelke and Hitz model
## Do this by simulating a data set of comparable size and transforming this back onto
## the scale of the original data
X_EH_Pareto <- lapply(fit_EH, function(x){rmpareto(n = n_data, model = "HR", par = x)})
X_EH_Uniform <- lapply(X_EH_Pareto, function(x){apply(x, 2, rank)/(1 + n_data)})
X_EH_Original_Margins <- lapply(1:n_sim, function(i){sapply(1:d, function(j){
  texmex:::revTransform(x = X_EH_Uniform[[i]][,j],
                        data = X[[i]][,j],
                        qu = qu_final[[i]][j],
                        th = u_final[[i]][j],
                        sigma = scale_final[[i]][j],
                        xi = shape_final[[i]][j],
                        method = "mixture")})})

## Heffernan and Tawn model
X_HT <- mcmapply(FUN = Sim_Surface_HT,
                 transforms = lapply(1:n_sim, function(i){lapply(1:d, function(j){X_to_Y[[j]][[i]]})}),
                 CMEVM_fits = lapply(1:n_sim, function(i){lapply(1:d, function(j){fit_HT[[j]][[i]]})}),
                 MoreArgs = list(n_sim = 50*n_data, q = dqu),
                 SIMPLIFY = FALSE,
                 mc.cores = detectCores() - 1)

## One-step graphical model
X_One_Step_Graph <- mcmapply(FUN = Sim_Surface_MVAGG,
                             transforms = lapply(1:n_sim, function(i){lapply(1:d, function(j){X_to_Y[[j]][[i]]})}),
                             CMEVM_fits = lapply(1:n_sim, function(i){lapply(1:d, function(j){fit_One_Step_Graph[[j]][[i]]})}),
                             MoreArgs = list(n_sim = 50*n_data, q = dqu),
                             SIMPLIFY = FALSE,
                             mc.cores = detectCores() - 1)

## Two-step graphical model
X_Two_Step_Graph <- mcmapply(FUN = Sim_Surface_MVAGG,
                             transforms = lapply(1:n_sim, function(i){lapply(1:d, function(j){X_to_Y[[j]][[i]]})}),
                             CMEVM_fits = lapply(1:n_sim, function(i){lapply(1:d, function(j){fit_Two_Step_Graph[[j]][[i]]})}),
                             MoreArgs = list(n_sim = 50*n_data, q = dqu),
                             SIMPLIFY = FALSE,
                             mc.cores = detectCores() - 1)

## Three-step independence model
X_Three_Step_Indep <- mcmapply(FUN = Sim_Surface_MVAGG,
                               transforms = lapply(1:n_sim, function(i){lapply(1:d, function(j){X_to_Y[[j]][[i]]})}),
                               CMEVM_fits = lapply(1:n_sim, function(i){lapply(1:d, function(j){fit_Three_Step_Indep[[j]][[i]]})}),
                               MoreArgs = list(n_sim = 50*n_data, q = dqu),
                               SIMPLIFY = FALSE,
                               mc.cores = detectCores() - 1)

## Three-step graphical model
X_Three_Step_Graph <- mcmapply(FUN = Sim_Surface_MVAGG,
                               transforms = lapply(1:n_sim, function(i){lapply(1:d, function(j){X_to_Y[[j]][[i]]})}),
                               CMEVM_fits = lapply(1:n_sim, function(i){lapply(1:d, function(j){fit_Three_Step_Graph[[j]][[i]]})}),
                               MoreArgs = list(n_sim = 50*n_data, q = dqu),
                               SIMPLIFY = FALSE,
                               mc.cores = detectCores() - 1)

## Three-step saturated model
X_Three_Step_Full <- mcmapply(FUN = Sim_Surface_MVAGG,
                              transforms = lapply(1:n_sim, function(i){lapply(1:d, function(j){X_to_Y[[j]][[i]]})}),
                              CMEVM_fits = lapply(1:n_sim, function(i){lapply(1:d, function(j){fit_Three_Step_Full[[j]][[i]]})}),
                              MoreArgs = list(n_sim = 50*n_data, q = dqu),
                              SIMPLIFY = FALSE,
                              mc.cores = detectCores() - 1)

## get the unconditioned random variables for probability estimation
uncon <- lapply(1:d, function(i){sapply(1:(d-1), function(j){
  combinations(n = d-1, r = j, v = (1:d)[-i])})})
uncon <- lapply(uncon, function(x){do.call(c, lapply(x, function(y){
  lapply(apply(y, 1, list), function(z){z[[1]]})
}))})

## threshold above which to calculate the probabilities
q_X <- 0.95
u_X <- apply(sapply(X, function(x){apply(x, 2, quantile, q_X)}), 1, max)

p_true_X <- t(sapply(1:d, function(i){
  sapply(uncon[[i]], function(z){
    p_data_surv_multi(data = X_prob_calc, cond = i, u = u_X, uncon = z)})}))

p_EH <- lapply(1:d, function(i){
  t(sapply(X_EH_Original_Margins, function(y){
    sapply(uncon[[i]], function(z){
      p_data_surv_multi(data = y, cond = i, u = u_X, uncon = z)})}))})

p_HT <- lapply(1:d, function(i){
  t(sapply(lapply(X_HT, function(x){x$Data_Margins}), function(y){
    sapply(uncon[[i]], function(z){
      p_data_surv_multi(data = y, cond = i, u = u_X, uncon = z)})}))})

p_One_Step_Graph <- lapply(1:d, function(i){
  t(sapply(lapply(X_One_Step_Graph, function(x){x$Data_Margins}), function(y){
    sapply(uncon[[i]], function(z){
      p_data_surv_multi(data = y, cond = i, u = u_X, uncon = z)})}))})

p_Two_Step_Graph <- lapply(1:d, function(i){
  t(sapply(lapply(X_Two_Step_Graph, function(x){x$Data_Margins}), function(y){
    sapply(uncon[[i]], function(z){
      p_data_surv_multi(data = y, cond = i, u = u_X, uncon = z)})}))})

p_Three_Step_Indep <- lapply(1:d, function(i){
  t(sapply(lapply(X_Three_Step_Indep, function(x){x$Data_Margins}), function(y){
    sapply(uncon[[i]], function(z){
      p_data_surv_multi(data = y, cond = i, u = u_X, uncon = z)})}))})

p_Three_Step_Graph <- lapply(1:d, function(i){
  t(sapply(lapply(X_Three_Step_Graph, function(x){x$Data_Margins}), function(y){
    sapply(uncon[[i]], function(z){
      p_data_surv_multi(data = y, cond = i, u = u_X, uncon = z)})}))})

p_Three_Step_Full <- lapply(1:d, function(i){
  t(sapply(lapply(X_Three_Step_Full, function(x){x$Data_Margins}), function(y){
    sapply(uncon[[i]], function(z){
      p_data_surv_multi(data = y, cond = i, u = u_X, uncon = z)})}))})


## Compare the probability estimation

## Collate the output
p_comp <- rep(list(vector("list", length(uncon[[1]]))), d)
for(i in 1:d){
  for(j in 1:length(uncon[[1]])){
    p_comp[[i]][[j]] <- cbind(rep(p_true_X[i,j], n_sim),
                              p_EH[[i]][,j],
                              p_HT[[i]][,j],
                              p_One_Step_Graph[[i]][,j],
                              p_Two_Step_Graph[[i]][,j],
                              p_Three_Step_Indep[[i]][,j],
                              p_Three_Step_Graph[[i]][,j],
                              p_Three_Step_Full[[i]][,j])
  }
}

plot_titles <- rep(list(vector("list", length(uncon[[1]]))), d)
for(i in 1:d){
  for(j in 1:length(uncon[[1]])){
    if(length(uncon[[i]][[j]]) == 1){
      plot_titles[[i]][[j]] <- substitute(P(X[j] > u[j] ~ "|" ~ X[i] > u[i]), list(i = i, j = uncon[[i]][[j]]))
    }
    else if(length(uncon[[i]][[j]]) == 2){
      plot_titles[[i]][[j]] <- substitute(P(X[j] > u[j], X[k] > u[k] ~ "|" ~ X[i] > u[i]), list(i = i, j = uncon[[i]][[j]][1], k = uncon[[i]][[j]][2]))
    }
    else if(length(uncon[[i]][[j]]) == 3){
      plot_titles[[i]][[j]] <- substitute(P(X[j] > u[j], X[k] > u[k], X[l] > u[l] ~ "|" ~ X[i] > u[i]), list(i = i, j = uncon[[i]][[j]][1], k = uncon[[i]][[j]][2], l = uncon[[i]][[j]][3]))
    }
    else{
      plot_titles[[i]][[j]] <- substitute(P(X[j] > u[j], X[k] > u[k], X[l] > u[l], X[m] > u[m] ~ "|" ~ X[i] > u[i]), list(i = i, j = uncon[[i]][[j]][1], k = uncon[[i]][[j]][2], l = uncon[[i]][[j]][3], m = uncon[[i]][[j]][4]))
    }
  }
}

# Compare the output above
method_vec <- c("True", "Engelke & Hitz", "Heffernan & Tawn",
                "One-step - Graphical", "Two-step - Graphical",
                "Three-step - Independence", "Three-step - Graphical", "Three-step - Saturated")
method_vec_1 <- method_vec[-1]
for(i in 1:d){
  for(j in 1:length(uncon[[1]])){
    ymin = floor(min(p_comp[[i]][[j]][,-1] - p_comp[[i]][[j]][,1])/0.1)*0.1
    ymax = ceiling(max(p_comp[[i]][[j]][,-1] - p_comp[[i]][[j]][,1])/0.1)*0.1
    
    data_to_plot <- data.frame(Value = c(p_comp[[i]][[j]][,-1] - p_comp[[i]][[j]][,1]))
    data_to_plot$Model = rep(method_vec_1, each = n_sim)
    data_to_plot$Model <- factor(data_to_plot$Model, levels = method_vec_1)
    
    pdf(paste0("Images/Simulation_Study/Mixed_Data/Probabilities/Site_", i, "/Prob_", j, ".pdf"), height = 10, width = 10)
    par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
    p_plot <- ggplot() + geom_boxplot(data = data_to_plot, aes(y = Value, fill = Model)) +
      lims(y = c(ymin, max(0.4, ymax))) +
      labs(y = "Bias", title = plot_titles[[i]][[j]]) +
      theme(legend.position = c(.95, .95),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6),
            legend.key.size = unit(0.75, 'cm'),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 10),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank()) +
      geom_hline(yintercept = 0, col = 2, linetype = "dashed", linewidth = 1)
    print(p_plot)
    dev.off()
  }
}

## Summary output from the model

## RMSE
RMSE_out <- lapply(1:d, function(i){t(sapply(1:length(uncon[[1]]), function(j){
  sapply(2:length(method_vec), function(k){
    RMSE(p_comp[[i]][[j]][,k], p_comp[[i]][[j]][,1])})}))})

RMSE_min <- sapply(RMSE_out, function(x){apply(x, 1, which.min)})

RMSE_Winner <- data.frame(matrix(0, ncol = d, nrow = length(method_vec_1)))
rownames(RMSE_Winner) <- method_vec_1
for(i in 1:d){
  RMSE_Winner[as.numeric(names(table(RMSE_min[,i]))), i] <- table(RMSE_min[,i])
}
RMSE_Winner$Total <- rowSums(RMSE_Winner)
RMSE_Winner

## Bias
Bias_out <- lapply(1:d, function(i){t(sapply(1:length(uncon[[1]]), function(j){
  sapply(2:length(method_vec), function(k){
    Bias(p_comp[[i]][[j]][,k], p_comp[[i]][[j]][,1])})}))})

Bias_min <- sapply(Bias_out, function(x){apply(x, 1, which.min)})

Bias_Winner <- data.frame(matrix(0, ncol = d, nrow = length(method_vec_1)))
rownames(Bias_Winner) <- method_vec_1
for(i in 1:d){
  Bias_Winner[as.numeric(names(table(Bias_min[,i]))), i] <- table(Bias_min[,i])
}
Bias_Winner$Total <- rowSums(Bias_Winner)
Bias_Winner

################################################################################
# out <- list(transforms = X_to_Y,
#             par_true = list(n_sim = n_sim, n_data = n_data,
#                             Gamma_MVP = Gamma_MVP, Gamma_MVN = Gamma_MVN, graph = graph_full))
# saveRDS(out, file = "Data/Mixed_Data_D5.RData")
################################################################################
