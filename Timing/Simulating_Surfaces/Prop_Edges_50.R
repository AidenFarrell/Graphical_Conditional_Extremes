################################################################################
## Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "gtools", "igraph", "LaplacesDemon", "rlang", "parallel")
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
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")

## Reading in functions for prediction
source("Prediction/Sim_Surfaces.R")

################################################################################
## Set up the simulation study

## True graph
d <- 100
prop_edges <- 0.5
all_edges <- combinations(n = d, r = 2, v = 1:d)
edges <- t(all_edges[sample(x = 1:nrow(all_edges), size = nrow(all_edges)*prop_edges, replace = FALSE),])
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

## Mean parameter
mu_true <- runif(n = d, min = -5, max = 5)

## True number of data points
n_exceedances <- 500
dqu <- 0.8
n_data <- 2500

## Get X data as list
X <- rmvn.sparse(n = n_data, mu = mu_true, CH = Cholesky(as(solve(rho_true), "sparseMatrix")))

## transform the data onto Laplace margins
X_list_by_data <- lapply(apply(X, 2, list), function(y){y[[1]]})
q_poss <- seq(0.55, 0.99, length.out = 100)
X_to_Y <- mcmapply(X_list_by_data,
                   FUN = X_to_Laplace,
                   MoreArgs = list(q = q_poss),
                   SIMPLIFY = FALSE,
                   mc.cores = detectCores() - 1)

## Get the output
Y <- sapply(X_to_Y, function(y){y$data$Y})

## Now we want to subset the data so that each component is large in turn
Y_u <- qlaplace(dqu)

## level i corresponds to data set ith component in the data set given that i is large
## level j corresponds to the jth data set
Y_Yi_large <- vector("list", d)
for(i in 1:d){
  Y_Yi_large[[i]] <- Y[which(Y[,i] > Y_u),]
}

################################################################################
## Fit the three-step graphical and saturated models
fit_three_step_graph <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
                                 data = Y_Yi_large,
                                 cond = 1:d,
                                 MoreArgs = list(graph = g_true,
                                                 v = ceiling(max(Y)) + 1,
                                                 maxit = 1e+9),
                                 SIMPLIFY = FALSE,
                                 mc.cores = detectCores() - 1)

fit_three_step_full <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
                                data = Y_Yi_large,
                                cond = 1:d,
                                MoreArgs = list(graph = make_full_graph(n = d),
                                                v = ceiling(max(Y)) + 1,
                                                maxit = 1e+9),
                                SIMPLIFY = FALSE,
                                mc.cores = detectCores() - 1)

## Simulate a single surface from both models
n_sim_surfaces <- 200
times <- matrix(NA, nrow = 2, ncol = 200)
rownames(times) <- c("Three-step - Graph", "Three-step - Saturated")

for(i in 1:n_sim_surfaces){
  start_time_three_step_graph <- Sys.time()
  sim_surface_three_step_graph <- Sim_Surface_MVAGG(n_sim = 5*n_data, q = dqu,
                                                    transforms = X_to_Y,
                                                    CMEVM_fits = fit_three_step_graph)
  end_time_three_step_graph <- Sys.time()
  times[1,i] <- difftime(end_time_three_step_graph, start_time_three_step_graph, units = "secs") 
}

for(i in 1:n_sim_surfaces){
  start_time_three_step_full <- Sys.time()
  sim_surface_three_step_full <- Sim_Surface_MVAGG(n_sim = 5*n_data, q = dqu,
                                                   transforms = X_to_Y,
                                                   CMEVM_fits = fit_three_step_full)
  end_time_three_step_full <- Sys.time()
  times[2,i] <- difftime(end_time_three_step_full, start_time_three_step_full, units = "secs")
}

round(apply(times, 1, mean), 2)

################################################################################