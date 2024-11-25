################################################################################
## Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "gtools", "igraph", "kableExtra", "parallel")
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
seed <- -1776067839
set.seed(seed)

## True graph
d <- 15
prop_edges <- 0.6
all_edges <- combinations(n = d, r = 2, v = 1:d)
edges <- t(all_edges[sample(x = 1:nrow(all_edges), size = nrow(all_edges)*prop_edges, replace = FALSE),])
g_true <- graph(edges = edges, directed = FALSE)
plot.igraph(g_true)

## Generate a precision matrix with this structure and get Sigma
set.seed(seed)
alpha_true <- runif(n = d, 0.1, 0.5)
beta_true <- runif(n = d, 0.1, 0.3)
loc_true <- runif(n = d, -5, 5)
scale_1_true <- runif(n = d, 0.5, 2)
scale_2_true <- runif(n = d, 1.5, 3)
shape_true <- runif(n = d, 0.8, 2.5)

simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(g_true)), 
                           v_sign = -1, v_within = c(0.1, 0.6))
Gamma_true <- simul$omega
Sigma_true <- solve(Gamma_true)
rho_true <- cov2cor(Sigma_true)

Sigma_true_i <- lapply(1:d, function(i){Cond_Sigma(Sigma_true, i)})
rho_true_i <- lapply(Sigma_true_i, function(x){cov2cor(x)})

## True number of data points
# n_data <- c(125, 250, 500, 1000, 2500, 5000, 10000, 20000)
n_excesses <- c(250, 500, 1000, 2000, 4000)
dqu <- 0.8

## Obtain Yi | Yi > u
Yi_large <- lapply(1:length(n_excesses), function(i){
  replicate(n = d, expr = qlaplace(runif(n_excesses[i], dqu, 1)), simplify = TRUE)})

## Obtain Y_i
a_yi <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){sapply(alpha_true[-j], function(x){Yi_large[[i]][,j]*x})})})
b_yi <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){sapply(beta_true[-j], function(x){Yi_large[[i]][,j]^x})})})

Z_not_i <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    rmvagg(n = n_excesses[i], 
           loc = loc_true[-j],
           scale_1 = scale_1_true[-j],
           scale_2 = scale_2_true[-j],
           shape = shape_true[-j],
           Gamma = as(solve(rho_true_i[[j]]), "sparseMatrix"))   
  })
})

Y_i <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    a_yi[[i]][[j]] + b_yi[[i]][[j]]*Z_not_i[[i]][[j]]
  })})

## Now get Y (i.e. add Yi_large into Y_i)
Y_Yi_large <- rep(list(vector("list", d)), length(n_excesses))
for(i in 1:length(n_excesses)){
  for(j in 1:d){
    if(j == 1){
      Y_Yi_large[[i]][[j]] <- cbind(Yi_large[[i]][,j], Y_i[[i]][[j]])
    }
    else if(j == d){
      Y_Yi_large[[i]][[j]] <- cbind(Y_i[[i]][[j]], Yi_large[[i]][,j])
    }
    else{
      Y_Yi_large[[i]][[j]] <- cbind(Y_i[[i]][[j]][,1:(j-1)], Yi_large[[i]][,j], Y_i[[i]][[j]][,j:(d-1)])
    }
  }
}

################################################################################
## Fix the starting parameters

start_par <- cbind(rep(0.1, d),
                   rep(0.1, d),
                   loc_true + runif(d, -0.25, 0.25),
                   pmax(0.5, scale_1_true + runif(d, -0.25, 0.25)),
                   pmin(2, pmax(0.5, scale_2_true + runif(d, -0.25, 0.25))),
                   pmin(2, pmax(0.5, shape_true + runif(d, -0.25, 0.25))))

start_par <- lapply(1:d, function(i){start_par[-i,]})
start_HT <- lapply(start_par, function(x){x[,c(1:2)]})
start_AGG <- lapply(start_par, function(x){x[,-c(1:2)]})

## Fix v
v <- ceiling(max(sapply(Y_Yi_large, function(x){sapply(x, max)}))) + 1

## Fit the models
fit_one_step_graph <- fit_one_step_full <- fit_two_step_graph <- fit_two_step_full <-
  fit_three_step_graph <- fit_three_step_full <- vector("list", length(n_excesses))
times <- matrix(NA, nrow = 6, ncol = length(n_excesses))
colnames(times) <- n_excesses
rownames(times) <- c("One-step - Graph", "Two-step - Graph", "Three-step - Graph", 
                     "One-step - Saturated", "Two-step - Saturated", "Three-step - Saturated")

for(j in 1:length(n_excesses)){
  ## One-step graphical model
  start_time_one_step_graph <- Sys.time()
  fit_one_step_graph[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG, 
                                      data = Y_Yi_large[[j]], 
                                      cond = 1:d, 
                                      start = start_par,
                                      MoreArgs = list(graph = g_true,
                                                      maxit = 1e+9), 
                                      SIMPLIFY = FALSE, 
                                      mc.cores = detectCores() - 1)
  end_time_one_step_graph <- Sys.time()
  times[1,j] <- difftime(end_time_one_step_graph, start_time_one_step_graph, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_excesses)){
  ## One-step saturated model
  start_time_one_step_full <- Sys.time()
  fit_one_step_full[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG, 
                                     data = Y_Yi_large[[j]], 
                                     cond = 1:d, 
                                     start = start_par,
                                     MoreArgs = list(graph = make_full_graph(n = d),
                                                     maxit = 1e+9), 
                                     SIMPLIFY = FALSE, 
                                     mc.cores = detectCores() - 1)
  end_time_one_step_full <- Sys.time()
  times[4,j] <- difftime(end_time_one_step_full, start_time_one_step_full, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_excesses)){
  ## Two-step graphical model
  start_time_two_step_graph <- Sys.time()
  fit_two_step_graph[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG_Two_Step, 
                                      data = Y_Yi_large[[j]], 
                                      cond = 1:d, 
                                      start_HT = start_HT,
                                      start_AGG = start_AGG,
                                      MoreArgs = list(graph = g_true,
                                                      v = v,
                                                      maxit = 1e+9), 
                                      SIMPLIFY = FALSE, 
                                      mc.cores = detectCores() - 1)
  end_time_two_step_graph <- Sys.time()
  times[2,j] <- difftime(end_time_two_step_graph, start_time_two_step_graph, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_excesses)){
  ## Two-step saturated model
  start_time_two_step_full <- Sys.time()
  fit_two_step_full[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG_Two_Step, 
                                     data = Y_Yi_large[[j]], 
                                     cond = 1:d, 
                                     start_HT = start_HT,
                                     start_AGG = start_AGG,
                                     MoreArgs = list(graph = make_full_graph(n = d),
                                                     v = v,
                                                     maxit = 1e+9), 
                                     SIMPLIFY = FALSE, 
                                     mc.cores = detectCores() - 1)
  end_time_two_step_full <- Sys.time()
  times[5,j] <- difftime(end_time_two_step_full, start_time_two_step_full, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_excesses)){
  ## Three-step graphical model
  start_time_three_step_graph <- Sys.time()
  fit_three_step_graph[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step, 
                                        data = Y_Yi_large[[j]], 
                                        cond = 1:d,
                                        start_HT = start_HT,
                                        start_AGG = start_AGG,
                                        MoreArgs = list(graph = g_true,
                                                        v = v,
                                                        maxit = 1e+9), 
                                        SIMPLIFY = FALSE, 
                                        mc.cores = detectCores() - 1)
  end_time_three_step_graph <- Sys.time()
  times[3,j] <- difftime(end_time_three_step_graph, start_time_three_step_graph, units = "secs")
  print(paste(j, "data sets complete"))
}

for(j in 1:length(n_excesses)){
  ## Three-step saturated model
  start_time_three_step_full <- Sys.time()
  fit_three_step_full[[j]] <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step, 
                                       data = Y_Yi_large[[j]], 
                                       cond = 1:d,
                                       start_HT = start_HT,
                                       start_AGG = start_AGG,
                                       MoreArgs = list(graph = make_full_graph(n = d),
                                                       v = v,
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


any(sapply(fit_one_step_graph, function(x){sapply(x, function(y){is.na(y$convergence)})}))
any(sapply(fit_one_step_full, function(x){sapply(x, function(y){is.na(y$convergence)})}))
any(sapply(fit_two_step_graph, function(x){sapply(x, function(y){is.na(y$convergence)})}))
any(sapply(fit_two_step_full, function(x){sapply(x, function(y){is.na(y$convergence)})}))
any(sapply(fit_three_step_graph, function(x){sapply(x, function(y){is.na(y$convergence)})}))
any(sapply(fit_three_step_full, function(x){sapply(x, function(y){is.na(y$convergence)})}))

################################################################################