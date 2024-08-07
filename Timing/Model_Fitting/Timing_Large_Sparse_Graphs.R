################################################################################
## Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "gtools", "glasso", "igraph", "LaplacesDemon", "rlang", "parallel")
# install.packages(required_pckgs, dependencies = TRUE)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Load in Rs script with functions for fitting the conditional extremes model with

## Reading in general functions
source("Miscellaneous_Functions/MVAGG_Functions.R")
source("Miscellaneous_Functions/General_Functions.R")
source("Miscellaneous_Functions/Transformations.R")

## Read in R Script to fit the Heffernan and Tawn Model
source("Model_Fitting/Cond_Extremes_MVN_Residuals.R")

################################################################################
fit_comp_three_step <- function(alpha, beta, loc, scale_1, scale_2, shape, rho, cond, n_excesses, dqu, graph = NA,
                                constrain = TRUE, q = c(0,1), v = 20, aLow = -1, maxit = 1e+6, nOptim = 1, start_HT, start_AGG){
  
  ## Get information required
  d <- length(alpha)
  if(length(cond) > 1){
    stop("cond must be of length 1")
  }
  else if(cond <= 0 | cond > d | cond%%1 != 0){
    stop("cond must be an integer in the interval (0,d]")
  }
  
  ## Obtain Yi | Yi > u
  Yi_large <- qlaplace(runif(n_excesses, dqu, 1))
  
  ## Obtain Y_i
  a_yi <- sapply(alpha[-cond], function(x){Yi_large*x})
  b_yi <- sapply(beta[-cond], function(x){Yi_large^x})
  Z_not_i <- rmvagg(n = n_excesses, 
                    loc = loc_true[-cond],
                    scale_1 = scale_1_true[-cond],
                    scale_2 = scale_2_true[-cond],
                    shape = shape_true[-cond],
                    Gamma = as(solve(rho_true_i[[cond]]), "sparseMatrix"))
  
  Y_i <- a_yi + b_yi*Z_not_i
  
  ## Now get Y (i.e. add Yi_large into Y_i)
  if(cond == 1){
    Y_Yi_large <- cbind(Yi_large, Y_i)
  }
  else if(cond == d){
    Y_Yi_large <- cbind(Y_i, Yi_large)
  }
  else{
    Y_Yi_large <- cbind(Y_i[,1:(cond-1)], Yi_large, Y_i[,cond:(d-1)])
  }
  
  ## Now fit the models
  ## Obtain the starting parameters
  if(missing(start_HT)){
    start_HT <- rep(0.1, 2)
  }
  if(missing(start_AGG)){
    start_AGG <- matrix(rep(c(0, 1.5, 2, 1.5), d-1), nrow = d-1, ncol = 4, byrow = TRUE)
  }
  else if(!is.numeric(start_HT) | !is.numeric(start_AGG)){
    stop("start_HT and start_AGG must be vectors")
  }
  if(length(start_AGG) == 4){
    start_AGG <- matrix(rep(start_AGG, d-1), ncol = 4, byrow = TRUE)
  }
  if(length(start_HT) != 2 | length(start_AGG) != 4*(d-1)){
    stop("start_HT and start_AGG must be vectors of length 2 or 2(d-1) and 4 or 4(d-1), respectively")
  }
  else if(any(abs(start_HT) > 1)){
    stop("Initial starting values are outside the parameter space for Heffernan and Tawn model")
  }
  else if(any(start_AGG[,-1] <= 0)){
    stop("Initial starting values are outside the parameter space for MVAGG")
  }
  
  ## Step 1
  ## Fit the original conditional multivaraite extreme value model to the data to
  ## obtain the fitted residuals and dependence parameters
  yex <- as.matrix(Y_Yi_large[,cond])
  ydep <- lapply(apply(Y_Yi_large[,-cond], 2, list), function(x){x[[1]]})
  
  start_time_HT <- Sys.time()
  res_HT <- mcmapply(FUN = qfun_MVN_indep,
                     ydep = ydep,
                     MoreArgs = list(yex = yex,
                                     constrain = constrain, aLow = aLow, q = q, v = ceiling(max(Y_Yi_large)) + 1,
                                     maxit = maxit, start = start_HT, nOptim = nOptim),
                     SIMPLIFY = FALSE,
                     mc.cores = detectCores() - 1)
  end_time_HT <- Sys.time()
  
  ## Get the necessary output
  z <- lapply(res_HT, function(x){x$Z})
  
  start_AGG[,1] <- sapply(z, mean)
  start_AGG <- lapply(apply(start_AGG, 1, list), function(x){x[[1]]})
  
  start_time_AGG <- Sys.time()
  res_AGG <- mcmapply(FUN = fit_agg,
                      data = z,
                      par = start_AGG,
                      SIMPLIFY = FALSE,
                      mc.cores = detectCores() - 1)
  index <- which(sapply(res_AGG, function(x){is.na(x$convergence)}))
  while(!is_empty(index)){
    for(j in 1:length(index)){
      res_AGG[[index[j]]] <-  fit_agg(data = z[,index[j]],
                                      par = c(start_AGG[index[j],1],
                                              runif(1, 0.5, 2),
                                              runif(1, 1, 3),
                                              runif(1, 0.5, 2.5)))
    }
    index <- which(sapply(res_AGG, function(x){is.na(x$convergence)}))
  }
  end_time_AGG <- Sys.time()
  
  ## Get the necessary output
  loc_hat <- sapply(res_AGG, function(x){x$par[1]})
  scale_1_hat <- sapply(res_AGG, function(x){x$par[2]})
  scale_2_hat <- sapply(res_AGG, function(x){x$par[3]})
  shape_hat <- sapply(res_AGG, function(x){x$par[4]})
  
  ## First transform the data onto standard Gaussian margins
  start_time_transformation <- Sys.time()
  Z_Uniform <- mcmapply(FUN = pagg,
                        q = z,
                        loc = loc_hat,
                        scale_1 = scale_1_hat,
                        scale_2 = scale_2_hat,
                        shape = shape_hat,
                        SIMPLIFY = TRUE,
                        mc.cores = detectCores() - 1)
  Z_Gaussian <- qnorm(Z_Uniform)
  end_time_transformation <- Sys.time()
  
  ## Graph is provided so we need to determine its structure
  start_time_graph <- Sys.time()
  graph_cond <- delete_vertices(graph, cond)
  ## Determine the missing edges in the graph
  all_edges <- as.data.frame(combinations(n = d-1, r = 2, v = 1:(d-1)))
  g_edges <- as.data.frame(as_edgelist(graph_cond))
  all_edges$exists <- do.call(paste0, all_edges) %in% do.call(paste0, g_edges)
  non_edges <- as.matrix(all_edges[which(all_edges$exists == FALSE), 1:2])
  
  ## Fit the graphical model
  fit_glasso <- try(suppressWarnings(glasso(s = cor(Z_Gaussian), rho = 0, nobs = n_excesses,
                                            zero = non_edges, thr = 1e-8, maxit = 1e+6)), silent = TRUE)
  end_time_graph <- Sys.time()
  
  ## Fit the saturated model
  start_time_full <- Sys.time()
  Gamma_full <- as(solve(cor(Z_Gaussian)), "sparseMatrix")
  end_time_full <- Sys.time()
  
  times <- c(difftime(end_time_HT, start_time_HT, units = "secs"),
             difftime(end_time_AGG, start_time_AGG, units = "secs"),
             difftime(end_time_transformation, start_time_transformation, units = "secs"),
             difftime(end_time_graph, start_time_graph, units = "secs"),
             difftime(end_time_full, start_time_full, units = "secs"))
  return(times)
}

################################################################################
## Set up the simulation study
seed <- -1776067839
set.seed(seed)

## True graph
d <- 200
prop_edges <- 0.1
all_edges <- combinations(n = d, r = 2, v = 1:d)
edges <- t(all_edges[sample(x = 1:nrow(all_edges), size = nrow(all_edges)*prop_edges, replace = FALSE),])
g_true <- graph(edges = edges, directed = FALSE)

## Get the dependence and AGG parameters
set.seed(seed)
alpha_true <- runif(n = d, 0.1, 0.5)
beta_true <- runif(n = d, 0.1, 0.3)
loc_true <- runif(n = d, -5, 5)
scale_1_true <- runif(n = d, 0.5, 2)
scale_2_true <- runif(n = d, 1.5, 3)
shape_true <- runif(n = d, 0.8, 2.5)

## Generate a precision matrix with this structure and get Sigma
simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(g_true)), v_sign = -1,
                           v_within = c(0.3, 0.7))
Gamma_true <- simul$omega
Sigma_true <- solve(Gamma_true)
rho_true <- cov2cor(Sigma_true)

Sigma_true_i <- mcmapply(FUN = Cond_Sigma,
                         j = 1:d,
                         MoreArgs = list(x = Sigma_true),
                         SIMPLIFY = FALSE,
                         mc.cores = detectCores() - 1)
rho_true_i <- lapply(Sigma_true_i, function(x){cov2cor(x)})

## Generate the data
dqu <- 0.8
n_excesses <- 4000

################################################################################

## Perform the timing experiment
times <- matrix(NA, nrow = d, ncol = 5)
for(i in 1:d){
  times[i,] <- fit_comp_three_step(alpha = alpha_true,
                                     beta = beta_true,
                                     loc = loc_true,
                                     scale_1 = scale_1_true,
                                     scale_2 = scale_2_true,
                                     shape = shape_true,
                                     rho = rho_true_i,
                                     n_excesses = n_excesses,
                                     dqu = dqu,
                                     graph = g_true,
                                     cond = i) 
  print(paste(i, "of", d, "fits complete"))
}

colnames(times) <- c("HT", "AGG", "Transform", "Graph", "Full")
round(apply(times, 2, mean), 2)