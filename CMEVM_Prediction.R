################################################################################
## Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "ggplot2", "gtools", "igraph", "kableExtra", "LaplacesDemon", "rlang", "parallel", "purrr", "pryr")
# install.packages(required_pckgs, dependencies = TRUE)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Load in Rs script with functions for fitting the conditional extremes model with
# setwd("/home/farrel11/Graphical_Models/")

## General functions
source("Miscellaneous_Functions/General_Functions.R")

## Functions for MVAGG
source("Miscellaneous_Functions/MVAGG_Functions.R")

## Model fitting scripts
source("Model_Fitting/Cond_Extremes_MVN_Residuals.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")

################################################################################

## Function to generate data from the model
Y_from_MVAGG <- function(n_data, p_y, con, alpha, beta, loc, scale_1, scale_2, shape, Gamma){
  ## checks on the inputs
  if(length(n_data) > 1 | n_data <= 0 | n_data%%1 != 0){
    stop("n_data must be a single positive integer")
  }
  if(p_y < 0 | p_y > 1){
    stop("p_y must be a probability")
  }
  
  #checks on alpha and beta
  if(!is.vector(alpha) | any(abs(alpha) > 1)){
    stop("alpha must be a vector in the region [-1,1]^(d-1)")
  }
  if(!is.vector(beta) | any(beta > 1)){
    stop("beta must be a vector in the region [-Infty,1]^(d-1)")
  }
  if(!is.vector(loc) | !is.vector(loc) | !is.vector(loc) | !is.vector(loc)){
    stop("loc, scale_1, scale_2, and shape must all be vectors")
  }
  d <- length(alpha) + 1
  if(length(con) > 1 | con <= 0 | con > d | con%%1 != 0){
    stop("con must be a single positive integer in the region [1,d]")
  }
  
  ## Simulate Yi | Yi > ui
  Yi <- qlaplace(runif(n_data, dqu, 1))
  
  ## Simulate Z from the MVAGG
  Z_not_i <- rmvagg(n = n_data,
                    loc = loc, scale_1 = scale_1, scale_2 = scale_2, shape = shape,
                    Gamma = Gamma)
  
  ## Transform Z onto standard laplace scale
  a_yi <- sapply(alpha, function(a){Yi*a})
  b_yi <- sapply(beta, function(b){Yi^b})
  Y_i <- a_yi + b_yi*Z_not_i
  
  if(con == 1){
    Y <- cbind(Yi, Y_i)
  }
  else if(con == d){
    Y <- cbind(Y_i, Yi)
  }
  else{
    Y <- cbind(Y_i[,1:(con-1)], Yi, Y_i[,con:(d-1)])
  }
  return(Y)
}

## Function to generate data from the model
Y_from_resid <- function(n_data, p_y, con, alpha, beta, Z){
  ## checks on the inputs
  if(length(n_data) > 1 | n_data <= 0 | n_data%%1 != 0){
    stop("n_data must be a single positive integer")
  }
  if(p_y < 0 | p_y > 1){
    stop("p_y must be a probability")
  }
  
  #checks on alpha and beta
  if(!is.vector(alpha) | any(abs(alpha) > 1)){
    stop("alpha must be a vector in the region [-1,1]^(d-1)")
  }
  if(!is.vector(beta) | any(beta > 1)){
    stop("beta must be a vector in the region [-Infty,1]^(d-1)")
  }
  
  d <- length(alpha) + 1
  if(length(con) > 1 | con <= 0 | con > d | con%%1 != 0){
    stop("con must be a single positive integer in the region [1,d]")
  }
  if(!is.matrix(Z)){
    stop("Z must be a matrix of dimension n = x (d-1)")
  }
  
  ## Simulate Yi | Yi > ui
  Yi <- qlaplace(runif(n_data, dqu, 1))
  
  ## Simulate Z from the MVAGG
  n_Z <- nrow(Z)
  Z_not_i <- Z[sample(1:n_Z, n_data, replace = TRUE),]
  
  ## Transform Z onto standard laplace scale
  a_yi <- sapply(alpha, function(a){Yi*a})
  b_yi <- sapply(beta, function(b){Yi^b})
  Y_i <- a_yi + b_yi*Z_not_i
  
  if(con == 1){
    Y <- cbind(Yi, Y_i)
  }
  else if(con == d){
    Y <- cbind(Y_i, Yi)
  }
  else{
    Y <- cbind(Y_i[,1:(con-1)], Yi, Y_i[,con:(d-1)])
  }
  return(Y)
}

## Function to calcualate joint survivor probabilites
surv_prob <- function(data, v, con, uncon) {
  
  ## Determine when the conditioning component is greater than v
  con_surv <- data[, con] > v
  
  ## Get the minimum across all the columns
  data_min <- do.call(pmin, lapply(c(con, uncon), function(x) data[, x]))
  
  ## Compute the conditional probability
  p_min <- sum(data_min > v)
  p_surv_con <- sum(con_surv)
  
  # Avoid division by zero
  if(p_surv_con == 0){
    p_out <- 0
  }
  else{
    p_out = p_min/p_surv_con
  }
  return(p_out)
}

## Function to:
## 1) Generate data from the true model
## 2) Fit the CMEVM and SCMEVMs
## 3) Simulate data from the fitted model
## 4) Calculate P[Y_{A} > u | Y_{i} > u]
Prob_Comp_Function <- function(par, n_excesses, dqu, cond_site, graph, n_model_sim, uncon_sites, u){
  ## Checks
  if(!is.list(par)){
    stop("par must be a list with each element of the list corresponding to the SCMEVM parameters")
  }
  d <- length(par$alpha)
  if(n_excesses <= 0 | n_excesses%%1 != 0){
    stop("n_excesses must be a positive intger")
  }
  if(dqu <= 0 | dqu >= 1){
    stop("dqu, the dependence threshold, is a probability in the region [0,1]")
  }
  if(cond_site <= 0 | cond_site > d | cond_site%%1 != 0){
    stop("cond_site must be a positive intger in the region [1,d]")
  }
  if(!is_igraph(graph) | length(V(graph)) != d){
    stop("graph is not an igraph object or does not have the correct number of vertices")
  }
  if(n_model_sim <= 0 | n_model_sim%%1 != 0){
    stop("n_model_sim must be a positive intger")
  }
  if(!is.list(uncon_sites)){
    stop("uncon_sites should be a list of the unconditioned sites such for calculating \n
         P[X_{A} > u | X_{i} > u")
  }
  
  ## simulate data from the model
  Y_Yi_large <- Y_from_MVAGG(n_data = n_excesses,
                             p_y = dqu,
                             con = cond_site,
                             alpha = par$alpha[-cond_site],
                             beta = par$beta[-cond_site],
                             loc = par$loc[-cond_site],
                             scale_1 = par$scale_1[-cond_site],
                             scale_2 = par$scale_2[-cond_site],
                             shape = par$shape[-cond_site],
                             Gamma = as(solve(par$rho), "sparseMatrix"))
  
  v <- ceiling(max(Y_Yi_large)) + 1
  
  ## Fit the various models
  fit_CMEVM <- Cond_Extremes_MVN(data = Y_Yi_large,
                                 graph = NA,
                                 cond = cond_site,
                                 v = v,
                                 maxit = 1e+9)
  
  fit_SCMEVM_Graph <-  try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large,
                                                          graph = graph,
                                                          cond = cond_site,
                                                          v = v,
                                                          maxit = 1e+9),
                           silent = TRUE)
  
  fit_SCMEVM_Full <-  try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large,
                                                         graph = make_full_graph(n = d),
                                                         cond = cond_site,
                                                         v = v,
                                                         maxit = 1e+9),
                          silent = TRUE)
  
  ## capture if we need different starting parameters for the SCMEVMs
  if(is.character(fit_SCMEVM_Graph)){
    count <- 0
    while(is.character(fit_SCMEVM_Graph)){
      start_par <- cbind(runif(d-1, -0.5, 0.5),
                         runif(d-1, 0.5, 5),
                         runif(d-1, 0.5, 5),
                         runif(d-1, 0.5, 3))
      fit_SCMEVM_Graph <-  try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large,
                                                              graph = graph,
                                                              cond = cond_site,
                                                              v = v,
                                                              start_AGG = start_par,
                                                              maxit = 1e+9),
                               silent = TRUE)
      count <- count + 1
      if(count >= 50){
        stop("50 starting values attempted for SCMEVM_Graph")
      }
    }
  }
  
  if(is.character(fit_SCMEVM_Full)){
    count <- 0
    while(is.character(fit_SCMEVM_Full)){
      start_par <- cbind(runif(d-1, -0.5, 0.5),
                         runif(d-1, 0.5, 5),
                         runif(d-1, 0.5, 5),
                         runif(d-1, 0.5, 3))
      fit_SCMEVM_Full <-  try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large,
                                                             graph = make_full_graph(n = d),
                                                             cond = cond_site,
                                                             v = v,
                                                             start_AGG = start_par,
                                                             maxit = 1e+9),
                              silent = TRUE)
      count <- count + 1
      if(count >= 50){
        stop("50 starting values attempted for SCMEVM_Graph")
      }
    }
  }
  
  ## Get simualted data sets from the fitted models
  Y_Yi_large_star_CMEVM <- Y_from_resid(n_data = n_model_sim,
                                        p_y = dqu,
                                        con = cond_site,
                                        alpha = fit_CMEVM$par$main[1,],
                                        beta = fit_CMEVM$par$main[2,],
                                        Z = fit_CMEVM$Z)
  
  Y_Yi_large_star_SCMEVM_Graph <- Y_from_MVAGG(n_data = n_model_sim,
                                               p_y = dqu,
                                               con = cond_site,
                                               alpha = fit_SCMEVM_Graph$par$main[1,],
                                               beta = fit_SCMEVM_Graph$par$main[2,],
                                               loc = fit_SCMEVM_Graph$par$main[3,],
                                               scale_1 = fit_SCMEVM_Graph$par$main[4,],
                                               scale_2 = fit_SCMEVM_Graph$par$main[5,],
                                               shape = fit_SCMEVM_Graph$par$main[6,],
                                               Gamma = fit_SCMEVM_Graph$par$Gamma)
  
  Y_Yi_large_star_SCMEVM_Full <- Y_from_MVAGG(n_data = n_model_sim,
                                              p_y = dqu,
                                              con = cond_site,
                                              alpha = fit_SCMEVM_Full$par$main[1,],
                                              beta = fit_SCMEVM_Full$par$main[2,],
                                              loc = fit_SCMEVM_Full$par$main[3,],
                                              scale_1 = fit_SCMEVM_Full$par$main[4,],
                                              scale_2 = fit_SCMEVM_Full$par$main[5,],
                                              shape = fit_SCMEVM_Full$par$main[6,],
                                              Gamma = fit_SCMEVM_Full$par$Gamma)
  
  ## Now get the predictions
  p_CMEVM <- mcmapply(FUN = surv_prob,
                      uncon = uncon_sites,
                      MoreArgs = list(data = Y_Yi_large_star_CMEVM,
                                      v = u,
                                      con = cond_site),
                      SIMPLIFY = TRUE,
                      mc.cores = detectCores() - 1)
  
  p_SCMEVM_Graph <- mcmapply(FUN = surv_prob,
                             uncon = uncon_sites,
                             MoreArgs = list(data = Y_Yi_large_star_SCMEVM_Graph,
                                             v = u,
                                             con = cond_site),
                             SIMPLIFY = TRUE,
                             mc.cores = detectCores() - 1)
  
  p_SCMEVM_Full <- mcmapply(FUN = surv_prob,
                            uncon = uncon_sites,
                            MoreArgs = list(data = Y_Yi_large_star_SCMEVM_Full,
                                            v = u,
                                            con = cond_site),
                            SIMPLIFY = TRUE,
                            mc.cores = detectCores() - 1)
  
  ## Save the output
  p_out <- data.frame(p_CMEVM = p_CMEVM,
                      p_SCMEVM_Graph = p_SCMEVM_Graph,
                      p_SCMEVM_Full = p_SCMEVM_Full)
  
  return(p_out)
}

## Function to:
## 1) Generate data from the true model
## 2) Fit the CMEVM and SCMEVMs
## 3) Simulate data from the fitted model
## 4) Create scatter plots of the data to show their different uses
Scatter_Plot_Comp_Function <- function(par, n_excesses, dqu, cond_site, graph, n_model_sim){
  ## Checks
  if(!is.list(par)){
    stop("par must be a list with each element of the list corresponding to the SCMEVM parameters")
  }
  d <- length(par$alpha)
  if(n_excesses <= 0 | n_excesses%%1 != 0){
    stop("n_excesses must be a positive intger")
  }
  if(dqu <= 0 | dqu >= 1){
    stop("dqu, the dependence threshold, is a probability in the region [0,1]")
  }
  if(cond_site <= 0 | cond_site > d | cond_site%%1 != 0){
    stop("cond_site must be a positive intger in the region [1,d]")
  }
  if(!is_igraph(graph) | length(V(graph)) != d){
    stop("graph is not an igraph object or does not have the correct number of vertices")
  }
  if(n_model_sim <= 0 | n_model_sim%%1 != 0){
    stop("n_model_sim must be a positive intger")
  }
  
  ## simulate data from the model
  Y_Yi_large <- Y_from_MVAGG(n_data = n_excesses,
                             p_y = dqu,
                             con = cond_site,
                             alpha = par$alpha[-cond_site],
                             beta = par$beta[-cond_site],
                             loc = par$loc[-cond_site],
                             scale_1 = par$scale_1[-cond_site],
                             scale_2 = par$scale_2[-cond_site],
                             shape = par$shape[-cond_site],
                             Gamma = as(solve(par$rho), "sparseMatrix"))
  
  v <- ceiling(max(Y_Yi_large)) + 1
  
  ## Fit the various models
  fit_CMEVM <- Cond_Extremes_MVN(data = Y_Yi_large,
                                 graph = NA,
                                 cond = cond_site,
                                 v = v,
                                 maxit = 1e+9)
  
  fit_SCMEVM_Graph <-  try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large,
                                                          graph = graph,
                                                          cond = cond_site,
                                                          v = v,
                                                          maxit = 1e+9),
                           silent = TRUE)
  
  ## capture if we need different starting parameters for the SCMEVMs
  if(is.character(fit_SCMEVM_Graph)){
    count <- 0
    while(is.character(fit_SCMEVM_Graph)){
      start_par <- cbind(runif(d-1, -0.5, 0.5),
                         runif(d-1, 0.5, 5),
                         runif(d-1, 0.5, 5),
                         runif(d-1, 0.5, 3))
      fit_SCMEVM_Graph <-  try(Cond_Extremes_MVAGG_Three_Step(data = Y_Yi_large,
                                                              graph = graph,
                                                              cond = cond_site,
                                                              v = v,
                                                              start_AGG = start_par,
                                                              maxit = 1e+9),
                               silent = TRUE)
      count <- count + 1
      if(count >= 50){
        stop("50 starting values attempted for SCMEVM_Graph")
      }
    }
  }
  
  ## Get simualted data sets from the fitted models
  Y_Yi_large_star_CMEVM <- Y_from_resid(n_data = n_model_sim,
                                        p_y = dqu,
                                        con = cond_site,
                                        alpha = fit_CMEVM$par$main[1,],
                                        beta = fit_CMEVM$par$main[2,],
                                        Z = fit_CMEVM$Z)
  
  Y_Yi_large_star_SCMEVM_Graph <- Y_from_MVAGG(n_data = n_model_sim,
                                               p_y = dqu,
                                               con = cond_site,
                                               alpha = fit_SCMEVM_Graph$par$main[1,],
                                               beta = fit_SCMEVM_Graph$par$main[2,],
                                               loc = fit_SCMEVM_Graph$par$main[3,],
                                               scale_1 = fit_SCMEVM_Graph$par$main[4,],
                                               scale_2 = fit_SCMEVM_Graph$par$main[5,],
                                               shape = fit_SCMEVM_Graph$par$main[6,],
                                               Gamma = fit_SCMEVM_Graph$par$Gamma)
  
  ## Scatter plot of the data
  
  ## Look at bivariate plots with the conditioning variable
  i <- sample((1:d)[-cond_site], 1)
  
  ## Only plot 2000 of the simulated data points as otherwise it will take too long to plot
  ## and distort the message
  n_points <- 2000
  index_plot <- sort(sample(1:n_model_sim, size = min(n_points, n_model_sim)))
  
  methods <- c("True", "CMEVM", "SCMEVM - Graphical")
  n_methods <- length(methods)
  
  dat_to_plot <- data.frame(Sample_Number = c(1:n_excesses, rep(index_plot, times = n_methods - 1)),
                            Method = c(rep(x = methods[1], times = n_excesses),
                                       rep(x = methods[-1], each = n_points)),
                            x = c(Y_Yi_large[,cond_site],
                                  Y_Yi_large_star_CMEVM[index_plot,cond_site],
                                  Y_Yi_large_star_SCMEVM_Graph[index_plot,cond_site]),
                            y = c(Y_Yi_large[,i],
                                  Y_Yi_large_star_CMEVM[index_plot,i],
                                  Y_Yi_large_star_SCMEVM_Graph[index_plot,i]))
  dat_to_plot$Method = factor(dat_to_plot$Method)
  
  ## Plot the sampel cloud
  plot_1 <- ggplot(data = dat_to_plot, aes(x = x, y = y)) +
    geom_point() +
    labs(x = substitute(Y[k] ~ "|" ~ Y[k] > u, list(k = cond_site)),
         y = substitute(Y[i] ~ "|" ~ Y[k] > u, list(i = i, k = cond_site))) +
    theme(strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          title = element_text(size = 16)) +
    facet_wrap(~ Method, ncol = n_methods)
  
  ## Now look when we don't include the conditioning site
  j <- sample((1:d)[-c(i, cond_site)], 1)
  
  dat_to_plot <- data.frame(Sample_Number = c(1:n_excesses, rep(index_plot, times = n_methods - 1)),
                            Method = c(rep(x = methods[1], times = n_excesses),
                                       rep(x = methods[-1], each = n_points)),
                            x = c(Y_Yi_large[,i],
                                  Y_Yi_large_star_CMEVM[index_plot,i],
                                  Y_Yi_large_star_SCMEVM_Graph[index_plot,i]),
                            y = c(Y_Yi_large[,j],
                                  Y_Yi_large_star_CMEVM[index_plot,j],
                                  Y_Yi_large_star_SCMEVM_Graph[index_plot,j]))
  dat_to_plot$Method = factor(dat_to_plot$Method)
  
  ## Plot the sampel cloud
  plot_2 <- ggplot(data = dat_to_plot, aes(x = x, y = y)) +
    geom_point() +
    labs(x = substitute(Y[i] ~ "|" ~ Y[k] > u, list(i = i, k = cond_site)),
         y = substitute(Y[j] ~ "|" ~ Y[k] > u, list(j = j, k = cond_site))) +
    theme(strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          title = element_text(size = 16)) +
    facet_wrap(~ Method, ncol = n_methods)
  
  ## Get the output
  return(list(plot_1, plot_2))
}

################################################################################
## Set the seed
seed <- 895266713
set.seed(seed)

## True graph
d <- 20
prop_edges <- 0.2
all_edges <- combinations(n = d, r = 2, v = 1:d)
g_true <- make_empty_graph()
while(length(V(g_true)) != d && is_connected(g_true) == FALSE){
  edges <- t(all_edges[sample(x = 1:nrow(all_edges), size = nrow(all_edges)*prop_edges, replace = FALSE),])
  g_true <- graph(edges = edges, directed = FALSE) 
}

## Get the dependence and AGG parameters
alpha_true <- runif(n = d, 0.1, 0.3)
beta_true <- runif(n = d, 0.1, 0.2)
loc_true <- runif(n = d, -1, 1)
scale_1_true <- runif(n = d, 0.5, 1)
scale_2_true <- runif(n = d, 1.5, 2)
shape_true <- runif(n = d, 0.8, 2.5)

## True correlation matrix
simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(g_true)), 
                           v_sign = -1, v_within = c(0, 0.1))
Gamma_true <- simul$omega
Sigma_true <- solve(Gamma_true)
rho_true <- cov2cor(Sigma_true)

## True sample sizes and threshold
n_boot <- 200
dqu <- 0.8
n_excesses <- 250
n_model_sim <- 5e+6

## Choose the site we wish to condition on
cond_site <- sample(1:d, size = 1)

## Get the conditional matrix for the conditioning site
rho_true_cond_site <- cov2cor(Cond_Sigma(Sigma_true, cond_site))

## Get a list for the parameters
par_true <- list(alpha = alpha_true,
                 beta = beta_true,
                 loc = loc_true,
                 scale_1 = scale_1_true,
                 scale_2 = scale_2_true,
                 shape = shape_true,
                 rho = rho_true_cond_site)

## We want to assess P[X_{A} > u | X_{i} > u]

## Let u be the 0.999 quantile of the standard Laplace distribution
u <- qlaplace(0.999)

## Consider A is a set of size 3
## Only look at 500 of the possible probabilities as there are too many otherwise
n_combos <- 500
tri_combos <- combinations(n = d-1, r = 3, v = (1:d)[-cond_site])
tri_combos <- tri_combos[sample(1:nrow(tri_combos), n_combos),]
tri_combos <- lapply(seq_len(n_combos), function(i){tri_combos[i,]})

## Obtain the probabilities from the model
p_model <- vector("list", n_boot)
for(i in 1:n_boot){
  p_model[[i]] <- Prob_Comp_Function(par = par_true,
                                               n_excesses = n_excesses,
                                               dqu = dqu,
                                               cond_site = cond_site,
                                               graph = g_true,
                                               n_model_sim = n_model_sim,
                                               uncon_sites = tri_combos,
                                               u = u)
  print(paste(i, "samples of", n_boot, "complete"))
}

## Obtain probability from the true distribution
## Simulate one very large sample from the true distribution for "true" probabilities
Y_Yi_large_true <- 
  Y_from_MVAGG(n_data = 1e+7,
               p_y = dqu,
               con = cond_site,
               alpha = alpha_true[-cond_site],
               beta = beta_true[-cond_site],
               loc = loc_true[-cond_site],
               scale_1 = scale_1_true[-cond_site],
               scale_2 = scale_2_true[-cond_site],
               shape = shape_true[-cond_site],
               Gamma = as(solve(rho_true_cond_site), "sparseMatrix"))

## True probability
p_true <- mcmapply(FUN = surv_prob,
                   uncon = tri_combos,
                   MoreArgs = list(data = Y_Yi_large_true,
                                   v = u,
                                   con = cond_site),
                   SIMPLIFY = TRUE,
                   mc.cores = detectCores() - 1)

## Obtain MAE and RMSE metrics
## Look at summary statistics and plots
p_CMEVM <- sapply(p_model, function(x){x$p_CMEVM})
p_SCMEVM_Graph <- sapply(p_model, function(x){x$p_SCMEVM_Graph})
p_SCMEVM_Full <- sapply(p_model, function(x){x$p_SCMEVM_Full})

Bias_out <- sapply(1:n_combos, function(j){
  c(Bias(x = p_CMEVM[j,], xhat = p_true[j]),
    Bias(x = p_SCMEVM_Graph[j,], xhat = p_true[j]),
    Bias(x = p_SCMEVM_Full[j,], xhat = p_true[j]))
})
Bias_tab <- as.numeric(table(apply(Bias_out, 2, which.min)))

RMSE_out <- sapply(1:n_combos, function(j){
  c(RMSE(x = p_CMEVM[j,], xhat = p_true[j]),
    RMSE(x = p_SCMEVM_Graph[j,], xhat = p_true[j]),
    RMSE(x = p_SCMEVM_Full[j,], xhat = p_true[j]))
})
RMSE_tab <- as.numeric(table(apply(RMSE_out, 2, which.min)))

## Print the output
methods <- c("CMEVM", "SCMEVM - Graphical", "SCMEVM - Saturated")
n_methods <- length(methods)
MAE_RMSE_out <- data.frame(Bias = Bias_tab,
                           RMSE = RMSE_tab)
rownames(MAE_RMSE_out) <- methods
print(MAE_RMSE_out)

## Plot the output
## Plot the output
p_comp <- data.frame(Prob_Index = rep(1:n_combos, times = n_methods),
                     Method = rep(methods, each = n_combos),
                     True_Prob = rep(p_true, n_methods),
                     Model_Prob = c(apply(p_CMEVM, 1, quantile, probs = 0.5),
                                    apply(p_SCMEVM_Graph, 1, quantile, probs = 0.5),
                                    apply(p_SCMEVM_Full, 1, quantile, probs = 0.5)),
                     Model_SE = c(apply(p_CMEVM, 1, sd),
                                  apply(p_SCMEVM_Graph, 1, sd),
                                  apply(p_SCMEVM_Full, 1, sd)))

pdf("/home/farrel11/Graphical_Models/Images/Prob_Comp.pdf", height = 10, width = 10)
ggplot(data = p_comp, aes(x = exp(True_Prob), y = exp(Model_Prob), col = Model_SE)) +
  geom_point() +
  scale_colour_gradient(low = "blue", high = "red") +
  labs(x = "True Probability (Exponential Scale)",
       y = "Model Probability (Exponential Scale)",
       col = "Model SE",
       title = substitute(P(Y[A] > v ~ "|" ~ Y[i] > v))) +
  theme(legend.position = "right",
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        title = element_text(size = 16)) +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Method, ncol = n_methods)
dev.off()

################################################################################
## Get some scatter plots of the data
plot_out <- Scatter_Plot_Comp_Function(par = par_true,
                                       n_excesses = n_excesses,
                                       dqu = dqu,
                                       cond_site = cond_site,
                                       graph = g_true,
                                       n_model_sim = n_model_sim)

pdf("/home/farrel11/Graphical_Models/Images/Sample_Clouds_Cond_Site.pdf", height = 10, width = 10)
print(plot_out[[1]])
dev.off()

pdf("/home/farrel11/Graphical_Models/Images/Sample_Clouds_Not_Cond_Site.pdf", height = 10, width = 10)
print(plot_out[[2]])
dev.off()
## Shows the rays in the simualted data from the fitted model using the Heffernan
## and Tawn approach
