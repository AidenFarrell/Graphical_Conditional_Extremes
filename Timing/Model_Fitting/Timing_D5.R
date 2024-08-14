################################################################################
## Load in required packages
rm(list = ls())
required_pckgs <- c("dplyr", "fake", "gtools", "igraph", "kableExtra", "LaplacesDemon", "rlang", "parallel")
# install.packages(required_pckgs, dependencies = TRUE)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Load in Rs script with functions for fitting the conditional extremes model with

## Reading in general functions
source("Miscellaneous_Functions/MVAGG_Functions.R")
source("Miscellaneous_Functions/General_Functions.R")
source("Miscellaneous_Functions/Transformations.R")

## read in threshold selection functions
source("threshold_selection_paper/helper_functions.R")
source("threshold_selection_paper/thresh_qq_metric.R")

## Reading in functions required for model fitting
source("Model_Fitting/Cond_Extremes_MVN_Residuals.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_One_Step.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Two_Step.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")

################################################################################

## Set up the simulation study

## True graph
d <- 5
c1 <- 1:3
c2 <- 3:5
cliques_true <- list(c1, c2)
edges <- do.call(rbind, lapply(cliques_true, function(x){combinations(length(x), r = 2, v = x)}))
edges <- t(edges[which(duplicated(edges) == FALSE),])
g_true <- graph(edges = edges, directed = FALSE)
plot.igraph(g_true)

## Generate a precision matrix with this structure and get Sigma
seed <- 57978638
set.seed(seed)
simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(g_true)), v_sign = -1,
                           v_within = c(0.3, 0.7))
Gamma_true <- simul$omega
Sigma_true <- solve(Gamma_true)
rho_true <- cov2cor(Sigma_true)
rho_true

## Mean parameter
mu_true <- runif(n = d, min = -5, max = 5)

## True number of data points
# n_data <- c(125, 250, 500, 1000, 2500, 5000, 10000, 20000)
n_exceedances <- c(250, 500, 1000, 2000, 4000)
dqu <- 0.8
n_data <- n_exceedances/(1 - dqu)

## Get X data as list
X <- vector("list", length(n_data))
for(i in 1:length(n_data)){
  X[[i]] <- rmvn.sparse(n = n_data[i], mu = mu_true, CH = Cholesky(as(solve(rho_true), "sparseMatrix")))
}

## transform the data onto Laplace margins
X_list_by_data <- lapply(X, function(x){lapply(apply(x, 2, list), function(y){y[[1]]})})
q_poss <- seq(0.55, 0.99, length.out = 100)
X_to_Y <- lapply(X_list_by_data, function(x){
  mcmapply(x = x,
           FUN = X_to_Laplace,
           MoreArgs = list(q = q_poss),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 1)})

## Get the output
Y <- lapply(X_to_Y, function(x){sapply(x, function(y){y$data$Y})})

## Now we want to subset the data so that each component is large in turn
Y_u <- qlaplace(dqu)

## level i corresponds to data set ith component in the data set given that i is large
## level j corresponds to the jth data set
Y_Yi_large <- rep(list(list()), length(n_data))
for(i in 1:d){
  for(j in 1:length(n_data)){
    Y_Yi_large[[j]][[i]] <- Y[[j]][which(Y[[j]][,i] > Y_u),]
  }
}

################################################################################
## Fix the starting parameters

start_par <- matrix(c(rep(0.1, d-1),
                      rep(0.1, d-1),
                      rep(0, d-1),
                      rep(1, d-1),
                      rep(2, d-1),
                      rep(1.5, d-1)),
                    nrow = d-1)

## Fir the models
fit_one_step_graph <- fit_one_step_full <- fit_two_step_graph <- fit_two_step_full <-
  fit_three_step_graph <- fit_three_step_full <- vector("list", length(n_data))
times <- matrix(NA, nrow = 6, ncol = length(n_data))
colnames(times) <- n_exceedances
rownames(times) <- c("One-step - Graph", "Two-step - Graph", "Three-step - Graph", 
                     "One-step - Saturated", "Two-step - Saturated", "Three-step - Saturated")

for(j in 1:length(n_data)){
  ## One-step graphical model
  start_time_one_step_graph <- Sys.time()
  fit_one_step_graph[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG, 
                                      data = Y_Yi_large[[j]], 
                                      cond = 1:d, 
                                      MoreArgs = list(graph = g_true,
                                                      maxit = 1e+9,
                                                      start = start_par), 
                                      SIMPLIFY = FALSE, 
                                      mc.cores = detectCores() - 1)
  end_time_one_step_graph <- Sys.time()
  times[1,j] <- difftime(end_time_one_step_graph, start_time_one_step_graph, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_data)){
  ## One-step saturated model
  start_time_one_step_full <- Sys.time()
  fit_one_step_full[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG, 
                                     data = Y_Yi_large[[j]], 
                                     cond = 1:d, 
                                     MoreArgs = list(graph = make_full_graph(n = d),
                                                     maxit = 1e+9,
                                                     start = start_par), 
                                     SIMPLIFY = FALSE, 
                                     mc.cores = detectCores() - 1)
  end_time_one_step_full <- Sys.time()
  times[4,j] <- difftime(end_time_one_step_full, start_time_one_step_full, units = "secs")
  print(paste(j, "data sets complete"))
}


for(j in 1:length(n_data)){
  ## Two-step graphical model
  start_time_two_step_graph <- Sys.time()
  fit_two_step_graph[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG_Two_Step, 
                                      data = Y_Yi_large[[j]], 
                                      cond = 1:d, 
                                      MoreArgs = list(graph = g_true,
                                                      v = ceiling(max(sapply(Y, max))) + 1,
                                                      maxit = 1e+9), 
                                      SIMPLIFY = FALSE, 
                                      mc.cores = detectCores() - 1)
  end_time_two_step_graph <- Sys.time()
  times[2,j] <- difftime(end_time_two_step_graph, start_time_two_step_graph, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_data)){
  ## Two-step saturated model
  start_time_two_step_full <- Sys.time()
  fit_two_step_full[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG_Two_Step, 
                                     data = Y_Yi_large[[j]], 
                                     cond = 1:d, 
                                     MoreArgs = list(graph = make_full_graph(n = d),
                                                     v = ceiling(max(sapply(Y, max))) + 1,
                                                     maxit = 1e+9), 
                                     SIMPLIFY = FALSE, 
                                     mc.cores = detectCores() - 1)
  end_time_two_step_full <- Sys.time()
  times[5,j] <- difftime(end_time_two_step_full, start_time_two_step_full, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_data)){
  ## Three-step graphical model
  start_time_three_step_graph <- Sys.time()
  fit_three_step_graph[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step, 
                                        data = Y_Yi_large[[j]], 
                                        cond = 1:d, 
                                        MoreArgs = list(graph = g_true,
                                                        v = ceiling(max(sapply(Y, max))) + 1,
                                                        maxit = 1e+9), 
                                        SIMPLIFY = FALSE, 
                                        mc.cores = detectCores() - 1)
  end_time_three_step_graph <- Sys.time()
  times[3,j] <- difftime(end_time_three_step_graph, start_time_three_step_graph, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_data)){
  ## Three-step saturated model
  start_time_three_step_full <- Sys.time()
  fit_three_step_full[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step, 
                                       data = Y_Yi_large[[j]], 
                                       cond = 1:d, 
                                       MoreArgs = list(graph = make_full_graph(n = d),
                                                       v = ceiling(max(sapply(Y, max))) + 1,
                                                       maxit = 1e+9), 
                                       SIMPLIFY = FALSE, 
                                       mc.cores = detectCores() - 1)
  end_time_three_step_full <- Sys.time()
  times[6,j] <- difftime(end_time_three_step_full, start_time_three_step_full, units = "secs")
  print(paste(j, "data sets complete"))
}


## Get the output
round(times, 2)

round(times, 2) %>% 
  kbl(format = "latex",
      col.names = colnames(times), align = "c") %>% 
  kable_classic(full_width = F, html_font = "Source Sans Pro")


################################################################################


