################################################################################
## Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "gtools", "igraph", "kableExtra", "LaplacesDemon", "rlang", "parallel", "purrr")
# install.packages(required_pckgs, dependencies = TRUE)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Load in Rs script with functions for fitting the conditional extremes model with

## General functions
source("Miscellaneous_Functions/General_Functions.R")

## Functions for MVAGG
source("Miscellaneous_Functions/MVAGG_Functions.R")

## Model fitting scripts
source("Model_Fitting/Cond_Extremes_MVN_Residuals.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_One_Step.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Two_Step.R")
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")

## Plotting functions
source("Simulation_Study/True_Distribution/Plotting_Functions.R")

################################################################################
## True graph
d <- 5
c1 <- 1:3
c2 <- 3:5
cliques_true <- list(c1, c2)
edges <- do.call(rbind, lapply(cliques_true, function(x){combinations(length(x), r = 2, v = x)}))
edges <- t(edges[which(duplicated(edges) == FALSE),])
g_true <- graph(edges = edges, directed = FALSE)
plot.igraph(g_true)

## True parameters in Z
# seed <- sample(.Random.seed, 1)
seed <- -1776067839
set.seed(seed)
alpha_true <- runif(n = d, 0.1, 0.5)
beta_true <- runif(n = d, 0.1, 0.3)
loc_true <- runif(n = d, -5, 5)
scale_1_true <- runif(n = d, 0.5, 2)
scale_2_true <- runif(n = d, 1.5, 3)
shape_true <- runif(n = d, 0.8, 2.5)

simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(g_true)), 
                           v_sign = 1, v_within = c(0.5, 0.9))
Gamma_true <- simul$omega

## Print for inclusion in paper
round(Gamma_true, 3) %>% 
  kbl(format = "latex", col.names = 1:d, align = "c") %>% 
  kable_classic(full_width = F, html_font = "Source Sans Pro")

Sigma_true <- solve(Gamma_true)
rho_true <- cov2cor(Sigma_true)

## Print for inclusion in paper
round(rho_true, 3) %>% 
  kbl(format = "latex", col.names = 1:d, align = "c") %>% 
  kable_classic(full_width = F, html_font = "Source Sans Pro")

Sigma_true_i <- lapply(1:d, function(i){Cond_Sigma(Sigma_true, i)})
rho_true_i <- lapply(Sigma_true_i, function(x){cov2cor(x)})

## Simulate large Yi
n_sim <- 200
dqu <- 0.8
n_excesses <- c(250, 500)

Yi_large <- vector("list", length = length(n_excesses))
for(i in 1:length(n_excesses)){
  Yi_large[[i]] <- replicate(n = d,
                             expr = replicate(n = n_sim, expr = qlaplace(runif(n_excesses[i], dqu, 1)), simplify = FALSE),
                             simplify = FALSE)
  
}


## Obtain Y_i
a_yi <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(Yi_large[[i]][[j]], function(x){sapply(alpha_true[-j], function(y){x*y})})})})
b_yi <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(Yi_large[[i]][[j]], function(x){sapply(beta_true[-j], function(y){x^y})})})})

Z_not_i <- rep(list(vector("list", d)), length(n_excesses))
for(i in 1:length(n_excesses)){
  for(j in 1:d){
    Z_not_i[[i]][[j]] <- replicate(n = n_sim, expr = rmvagg(n = n_excesses[i], 
                                                            loc = loc_true[-j],
                                                            scale_1 = scale_1_true[-j],
                                                            scale_2 = scale_2_true[-j],
                                                            shape = shape_true[-j],
                                                            Gamma = as(solve(rho_true_i[[j]]), "sparseMatrix")),
                                   simplify = FALSE)
  } 
}

Y_i <- rep(list(rep(list(vector("list", n_sim)), d)), length(n_excesses))
for(i in 1:length(n_excesses)){
  for(j in 1:d){
    for(k in 1:n_sim){
      Y_mat <- matrix(NA, nrow = n_excesses[i], ncol = d-1)
      for(l in 1:(d-1)){
        Y_mat[,l] <- a_yi[[i]][[j]][[k]][,l] + b_yi[[i]][[j]][[k]][,l]*Z_not_i[[i]][[j]][[k]][,l]
      }
      Y_i[[i]][[j]][[k]] <- Y_mat
    }
  } 
}

## Now get Y (i.e. add Yi_large into Y_i)
Y_Yi_large <- rep(list(rep(list(vector("list", n_sim)), d)), length(n_excesses))
for(i in 1:length(n_excesses)){
  for(j in 1:d){
    for(k in 1:n_sim){
      if(j == 1){
        Y_Yi_large[[i]][[j]][[k]] <- cbind(Yi_large[[i]][[j]][[k]], Y_i[[i]][[j]][[k]])
      }
      else if(j == d){
        Y_Yi_large[[i]][[j]][[k]] <- cbind(Y_i[[i]][[j]][[k]], Yi_large[[i]][[j]][[k]])
      }
      else{
        Y_Yi_large[[i]][[j]][[k]] <- cbind(Y_i[[i]][[j]][[k]][,1:(j-1)], Yi_large[[i]][[j]][[k]], Y_i[[i]][[j]][[k]][,j:(d-1)])
      }
    }
  }
}

################################################################################
## Fit models to the data

v <- ceiling(max(sapply(1:length(n_excesses), function(i){sapply(1:d, function(j){sapply(Y_Yi_large[[i]][[j]], max)})}))) + 1

## Heffernan and Tawn model
fit_HT <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVN,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = NA,
                             cond = j,
                             v = v,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

Z_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){lapply(y, function(z){z$Z})})})
a_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){lapply(y, function(z){z$par$main[1,]})})})
b_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){lapply(y, function(z){z$par$main[2,]})})})

## One-step model fits
## Independence
## Set nOptim = 5 so that the optimiser runs several times
## This is necessary as convergence can be difficult
fit_One_Step_Indep <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = NA,
                             cond = j,
                             maxit = 1e+9,
                             nOptim = 5),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

start_par_One_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    lapply(1:n_sim, function(k){
      cbind(rep(0.1, d-1), rep(0.1, d-1),
            fit_One_Step_Indep[[i]][[j]][[k]]$par$main[3,],
            rep(1.5, d-1), rep(2, d-1), rep(1.5, d-1))
    })
  })
})

## Graphical
fit_One_Step_Graph <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG,
             data = Y_Yi_large[[i]][[j]],
             start = start_par_One_Step[[i]][[j]],
             MoreArgs = list(graph = g_true,
                             cond = j,
                             maxit = 1e+9,
                             nOptim = 5),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Saturated
fit_One_Step_Full <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG,
             data = Y_Yi_large[[i]][[j]],
             start = start_par_One_Step[[i]][[j]],
             MoreArgs = list(graph = make_full_graph(n = d),
                             cond = j,
                             maxit = 1e+9,
                             nOptim = 5),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})


## Two-step model fits
## Give an informed start to the location parameter
loc_start_Two_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(1:n_sim, function(k){
    apply(Z_hat_HT[[i]][[j]][[k]], 2, mean)})})})
start_par_Two_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(1:n_sim, function(k){
    cbind(loc_start_Two_Step[[i]][[j]][[k]], rep(1.5, d-1), rep(1.5, d-1), rep(1.5, d-1))})})})

## Independence
fit_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG_Two_Step,
             data = Y_Yi_large[[i]][[j]],
             start_AGG = start_par_Two_Step[[i]][[j]],
             MoreArgs = list(graph = NA,
                             cond = j,
                             v = v,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Give some informed starting place to the graphical and saturated models
start_par_Two_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    lapply(1:n_sim, function(k){
      cbind(fit_Two_Step_Indep[[i]][[j]][[k]]$par$main[3,],
            rep(1.5, d-1), rep(2, d-1), rep(1.5, d-1))
    })
  })
})

## Graphical
fit_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG_Two_Step,
             data = Y_Yi_large[[i]][[j]],
             start_AGG = start_par_Two_Step[[i]][[j]],
             MoreArgs = list(graph = g_true,
                             cond = j,
                             v = v,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Saturated
fit_Two_Step_Full <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG_Two_Step,
             data = Y_Yi_large[[i]][[j]],
             start_AGG = start_par_Two_Step[[i]][[j]],
             MoreArgs = list(graph = make_full_graph(n = d),
                             cond = j,
                             v = v,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Three-step model fits
## Independence
fit_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = NA,
                             cond = j,
                             v = v,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Graphical
fit_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = g_true,
                             cond = j,
                             v = v,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Saturated
fit_Three_Step_Full <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = make_full_graph(n = d),
                             cond = j,
                             v = v,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

################################################################################
## One and two step models can be somewhat sensitive to the starting values
## here we need to do another procedure to ensure the models converge

Index_One_Step_Indep <- lapply(fit_One_Step_Indep, function(x){
  lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
count <- 0
while(any(sapply(Index_One_Step_Indep, function(x){sapply(x, length)}) > 0)){
  
  for(i in 1:length(n_excesses)){
    
    start_par_One_Step <- cbind(runif(d-1, 0.1, 0.5),
                                runif(d-1, 0.1, 0.5),
                                runif(d-1, -2, 2),
                                runif(d-1, 0.5, 3),
                                runif(d-1, 0.5, 3),
                                runif(d-1, 0.5, 3))
    
    for(j in 1:d){
      if(is_empty(Index_One_Step_Indep[[i]][[j]])){
        next()
      }
      else{
        ind <- Index_One_Step_Indep[[i]][[j]]
        for(k in 1:length(ind)){
          fit_One_Step_Indep[[i]][[j]][[ind[k]]] <-
            try(Cond_Extremes_MVAGG(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                    cond = j,
                                    graph = NA,
                                    start = start_par_One_Step,
                                    maxit = 1e+9,
                                    nOptim = 5),
                silent = TRUE)
        }
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_One_Step_Indep <- lapply(fit_One_Step_Indep, function(x){
    lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
  if(any(sapply(Index_One_Step_Indep, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_One_Step_Graph <- lapply(fit_One_Step_Graph, function(x){
  lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
count <- 0
while(any(sapply(Index_One_Step_Graph, function(x){sapply(x, length)}) > 0)){
  
  for(i in 1:length(n_excesses)){
    
    start_par_One_Step <- cbind(runif(d-1, 0.1, 0.5),
                                runif(d-1, 0.1, 0.5),
                                runif(d-1, -2, 2),
                                runif(d-1, 0.5, 3),
                                runif(d-1, 0.5, 3),
                                runif(d-1, 0.5, 3))
    
    for(j in 1:d){
      if(is_empty(Index_One_Step_Graph[[i]][[j]])){
        next()
      }
      else{
        ind <- Index_One_Step_Graph[[i]][[j]]
        for(k in 1:length(ind)){
          fit_One_Step_Graph[[i]][[j]][[ind[k]]] <-
            try(Cond_Extremes_MVAGG(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                    cond = j,
                                    graph = g_true,
                                    start = start_par_One_Step,
                                    maxit = 1e+9,
                                    nOptim = 5),
                silent = TRUE)
        }
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_One_Step_Graph <- lapply(fit_One_Step_Graph, function(x){
    lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
  if(any(sapply(Index_One_Step_Graph, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_One_Step_Full <- lapply(fit_One_Step_Full, function(x){
  lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
count <- 0
while(any(sapply(Index_One_Step_Full, function(x){sapply(x, length)}) > 0)){
  
  for(i in 1:length(n_excesses)){
    
    start_par_One_Step <- cbind(runif(d-1, 0.1, 0.5),
                                runif(d-1, 0.1, 0.5),
                                runif(d-1, -2, 2),
                                runif(d-1, 0.5, 3),
                                runif(d-1, 0.5, 3),
                                runif(d-1, 0.5, 3))
    
    for(j in 1:d){
      if(is_empty(Index_One_Step_Full[[i]][[j]])){
        next()
      }
      else{
        ind <- Index_One_Step_Full[[i]][[j]]
        for(k in 1:length(ind)){
          fit_One_Step_Full[[i]][[j]][[ind[k]]] <-
            try(Cond_Extremes_MVAGG(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                    cond = j,
                                    graph = make_full_graph(n = d),
                                    start = start_par_One_Step,
                                    maxit = 1e+9,
                                    nOptim = 5),
                silent = TRUE)
        }
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_One_Step_Full <- lapply(fit_One_Step_Full, function(x){
    lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
  if(any(sapply(Index_One_Step_Full, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_Two_Step_Indep <- lapply(fit_Two_Step_Indep, function(x){
  lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
count <- 0
while(any(sapply(Index_Two_Step_Indep, function(x){sapply(x, length)}) > 0)){
  
  for(i in 1:length(n_excesses)){
    
    start_par_Two_Step <- cbind(runif(d-1, -0.5, 0.5),
                                runif(d-1, 0.5, 5),
                                runif(d-1, 0.5, 5),
                                runif(d-1, 0.5, 3))
    
    for(j in 1:d){
      if(is_empty(Index_Two_Step_Indep[[i]][[j]])){
        next()
      }
      else{
        ind <- Index_Two_Step_Indep[[i]][[j]]
        for(k in 1:length(ind)){
          fit_Two_Step_Indep[[i]][[j]][[ind[k]]] <-
            try(Cond_Extremes_MVAGG_Two_Step(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                             cond = j,
                                             graph = NA,
                                             start_AGG = start_par_Two_Step,
                                             maxit = 1e+9),
                silent = TRUE)
        }
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_Two_Step_Indep <- lapply(fit_Two_Step_Indep, function(x){
    lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
  if(any(sapply(Index_Two_Step_Indep, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_Two_Step_Graph <- lapply(fit_Two_Step_Graph, function(x){
  lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
count <- 0
while(any(sapply(Index_Two_Step_Graph, function(x){sapply(x, length)}) > 0)){
  
  for(i in 1:length(n_excesses)){
    
    start_par_Two_Step <- cbind(runif(d-1, -0.5, 0.5),
                                runif(d-1, 0.5, 5),
                                runif(d-1, 0.5, 5),
                                runif(d-1, 0.5, 3))
    
    for(j in 1:d){
      if(is_empty(Index_Two_Step_Graph[[i]][[j]])){
        next()
      }
      else{
        ind <- Index_Two_Step_Graph[[i]][[j]]
        for(k in 1:length(ind)){
          fit_Two_Step_Graph[[i]][[j]][[ind[k]]] <-
            try(Cond_Extremes_MVAGG_Two_Step(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                             cond = j,
                                             graph = g_true,
                                             start_AGG = start_par_Two_Step,
                                             maxit = 1e+9),
                silent = TRUE)
        }
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_Two_Step_Graph <- lapply(fit_Two_Step_Graph, function(x){
    lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
  if(any(sapply(Index_Two_Step_Graph, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_Two_Step_Full <- lapply(fit_Two_Step_Full, function(x){
  lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
count <- 0
while(any(sapply(Index_Two_Step_Full, function(x){sapply(x, length)}) > 0)){
  
  for(i in 1:length(n_excesses)){
    
    start_par_Two_Step <- cbind(runif(d-1, -0.5, 0.5),
                                runif(d-1, 0.5, 5),
                                runif(d-1, 0.5, 5),
                                runif(d-1, 0.5, 3))
    
    for(j in 1:d){
      if(is_empty(Index_Two_Step_Full[[i]][[j]])){
        next()
      }
      else{
        ind <- Index_Two_Step_Full[[i]][[j]]
        for(k in 1:length(ind)){
          fit_Two_Step_Full[[i]][[j]][[ind[k]]] <-
            try(Cond_Extremes_MVAGG_Two_Step(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                             cond = j,
                                             graph = make_full_graph(n = d),
                                             start_AGG = start_par_Two_Step,
                                             maxit = 1e+9),
                silent = TRUE)
        }
      }
    }
  }
  count <- count + 1
  print(paste0(count, " set of starting parameters tested"))
  if(count >= 50){
    stop("50 starting values attempted")
  }
  
  Index_Two_Step_Full <- lapply(fit_Two_Step_Full, function(x){
    lapply(x, function(y){which(sapply(y, function(z){is.na(z$convergence)}))})})
  if(any(sapply(Index_Two_Step_Full, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

################################################################################
## Extract the parameter estimates
lower_tri_elements <- lower.tri(Gamma_true, diag = TRUE)

## One-step Independent Residuals
a_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

## One-step Graphical residuals
a_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

## One-step Saturated Residuals
a_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

## Two-step Independent Residuals
a_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

## Two-step Graphical residuals
a_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

## Two-step Saturated Residuals
a_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

## Three-Step Independent Residuals
a_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

## Three-Step Graphical residuals
a_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

## Three-Step Saturated Residuals
a_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[1,]), j)}))})})
b_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[2,]), j)}))})})
loc_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[3,]), j)}))})})
scale_1_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[4,]), j)}))})})
scale_2_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[5,]), j)}))})})
shape_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(y$par$main[6,]), j)}))})})
Gamma_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(as.matrix(y$par$Gamma), j)[lower_tri_elements]}))})})

################################################################################
## Assess convergence to true parameters
method_vec <- c("One-step - Independence", "One-step - Graphical", "One-step - Saturated", "Two/Three-step")

## Alpha
Alpha_plot <- comp_plots(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(a_hat_One_Step_Indep[[i]][[j]], 
         a_hat_One_Step_Graph[[i]][[j]], 
         a_hat_One_Step_Full[[i]][[j]],
         a_hat_Two_Step_Indep[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(alpha)[j ~ "|" ~ i]), 
  ylims = c(-1.25, 1), par_true = alpha_true)
pdf("Images/Simulation_Study/True_Distribution/Negative_Dependence/Alpha.pdf", width = 15, height = 10)
print(Alpha_plot)
dev.off()

## Beta
Beta_plot <- comp_plots(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(b_hat_One_Step_Indep[[i]][[j]], 
         b_hat_One_Step_Graph[[i]][[j]], 
         b_hat_One_Step_Full[[i]][[j]],
         b_hat_Two_Step_Indep[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(beta)[j ~ "|" ~ i]), 
  ylims = c(-0.5, 0.5), par_true = beta_true)
pdf("Images/Simulation_Study/True_Distribution/Negative_Dependence/Beta.pdf", width = 15, height = 10)
print(Beta_plot)
dev.off()

method_vec <- c("One-step - Independence", "One-step - Graphical", "One-step - Saturated", 
                "Two-step - Independence", "Two-step - Graphical", "Two-step - Saturated",
                "Three-step")

## Location
Loc_plot <- comp_plots(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(loc_hat_One_Step_Indep[[i]][[j]], 
         loc_hat_One_Step_Graph[[i]][[j]], 
         loc_hat_One_Step_Full[[i]][[j]],
         loc_hat_Two_Step_Indep[[i]][[j]], 
         loc_hat_Two_Step_Graph[[i]][[j]], 
         loc_hat_Two_Step_Full[[i]][[j]],
         loc_hat_Three_Step_Indep[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(nu)[j ~ "|" ~ i]), 
  ylims = c(-1, 1), par_true = loc_true)
pdf("Images/Simulation_Study/True_Distribution/Negative_Dependence/Location.pdf", width = 20, height = 15)
print(Loc_plot)
dev.off()

## Scale (Left)
Scale_1_plot <- comp_plots(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(scale_1_hat_One_Step_Indep[[i]][[j]], 
         scale_1_hat_One_Step_Graph[[i]][[j]], 
         scale_1_hat_One_Step_Full[[i]][[j]],
         scale_1_hat_Two_Step_Indep[[i]][[j]], 
         scale_1_hat_Two_Step_Graph[[i]][[j]], 
         scale_1_hat_Two_Step_Full[[i]][[j]],
         scale_1_hat_Three_Step_Indep[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(kappa[1])[j ~ "|" ~ i]), 
  ylims = c(-0.5, 1), par_true = scale_1_true)
pdf("Images/Simulation_Study/True_Distribution/Negative_Dependence/Scale_Left.pdf", width = 20, height = 15)
print(Scale_1_plot)
dev.off()

## Scale (Right)
Scale_2_plot <- comp_plots(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(scale_2_hat_One_Step_Indep[[i]][[j]], 
         scale_2_hat_One_Step_Graph[[i]][[j]], 
         scale_2_hat_One_Step_Full[[i]][[j]],
         scale_2_hat_Two_Step_Indep[[i]][[j]], 
         scale_2_hat_Two_Step_Graph[[i]][[j]], 
         scale_2_hat_Two_Step_Full[[i]][[j]],
         scale_2_hat_Three_Step_Indep[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(kappa[2])[j ~ "|" ~ i]), 
  ylims = c(-1, 1.5), par_true = scale_2_true)
pdf("Images/Simulation_Study/True_Distribution/Negative_Dependence/Scale_Right.pdf", width = 20, height = 15)
print(Scale_2_plot)
dev.off()

## Shape
Shape_plot <- comp_plots(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(shape_hat_One_Step_Indep[[i]][[j]], 
         shape_hat_One_Step_Graph[[i]][[j]], 
         shape_hat_One_Step_Full[[i]][[j]],
         shape_hat_Two_Step_Indep[[i]][[j]], 
         shape_hat_Two_Step_Graph[[i]][[j]], 
         shape_hat_Two_Step_Full[[i]][[j]],
         shape_hat_Three_Step_Indep[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(delta)[j ~ "|" ~ i]), 
  ylims = c(-1, 1.5), par_true = shape_true)
pdf("Images/Simulation_Study/True_Distribution/Negative_Dependence/Shape.pdf", width = 20, height = 15)
print(Shape_plot)
dev.off()

method_vec <- c("One-step - Graphical", "One-step - Saturated", 
                "Two-step - Graphical", "Two-step - Saturated",
                "Three-step - Graphical", "Three-step - Saturated")

## Precision Matrix
Gamma_plot <- comp_plots_matrix(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(Gamma_hat_One_Step_Graph[[i]][[j]], Gamma_hat_One_Step_Full[[i]][[j]],
         Gamma_hat_Two_Step_Graph[[i]][[j]], Gamma_hat_Two_Step_Full[[i]][[j]],
         Gamma_hat_Three_Step_Graph[[i]][[j]], Gamma_hat_Three_Step_Full[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(Gamma)[~ "|" ~ i]), 
  ylims = c(-0.3, 0.3), cov_mat_true = Sigma_true, precision = TRUE)

pdf("Images/Simulation_Study/True_Distribution/Negative_Dependence/Gamma.pdf", width = 20, height = 15)
print(Gamma_plot)
dev.off()

################################################################################

## Compare the log-likelihoods from the various models to affirm that we don't lose
## information by doing the iterative approach

## Extract the log-likelihood for the various models
log_like_true <- 
  lapply(1:length(n_excesses), function(i){
    sapply(1:d, function(j){
      sapply(1:n_sim, function(k){
        sum(dmvagg(data = Z_not_i[[i]][[j]][[k]],
                   loc = loc_true[-j],
                   scale_1 = scale_1_true[-j],
                   scale_2 = scale_2_true[-j],
                   shape = shape_true[-j],
                   Gamma = as(solve(rho_true_i[[j]]), "sparseMatrix"), log = TRUE)) -
          sum(sapply(beta_true[-j], function(x){x*log(Yi_large[[i]][[j]][[k]])}))
      })})})

log_like_One_Step_Indep <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_One_Step_Indep[[i]][[j]][[k]]$loglike})})})

log_like_One_Step_Graph <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_One_Step_Graph[[i]][[j]][[k]]$loglike})})})

log_like_One_Step_Full <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_One_Step_Full[[i]][[j]][[k]]$loglike})})})

log_like_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_Two_Step_Indep[[i]][[j]][[k]]$loglike})})})

log_like_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_Two_Step_Graph[[i]][[j]][[k]]$loglike})})})

log_like_Two_Step_Full <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_Two_Step_Full[[i]][[j]][[k]]$loglike})})})

log_like_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_Three_Step_Indep[[i]][[j]][[k]]$loglike})})})

log_like_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_Three_Step_Graph[[i]][[j]][[k]]$loglike})})})

log_like_Three_Step_Full <- lapply(1:length(n_excesses), function(i){
  sapply(1:d, function(j){
    sapply(1:n_sim, function(k){fit_Three_Step_Full[[i]][[j]][[k]]$loglike})})})

## Compare median and 95% confidence interval of bias between the model-based and
## true log-likelihood
probs = c(0.5, 0.025, 0.975)

CIs_One_Step_Indep <- pmap(.l = list(x = log_like_true, y = log_like_One_Step_Indep),
                           .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_Two_Step_Indep <- pmap(.l = list(x = log_like_true, y = log_like_Two_Step_Indep),
                           .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_Three_Step_Indep <- pmap(.l = list(x = log_like_true, y = log_like_Three_Step_Indep),
                             .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_One_Step_Graph <- pmap(.l = list(x = log_like_true, y = log_like_One_Step_Graph),
                           .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_Two_Step_Graph <- pmap(.l = list(x = log_like_true, y = log_like_Two_Step_Graph),
                           .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_Three_Step_Graph <- pmap(.l = list(x = log_like_true, y = log_like_Three_Step_Graph),
                             .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_One_Step_Full <- pmap(.l = list(x = log_like_true, y = log_like_One_Step_Full),
                          .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_Two_Step_Full <- pmap(.l = list(x = log_like_true, y = log_like_Two_Step_Full),
                          .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_Three_Step_Full <- pmap(.l = list(x = log_like_true, y = log_like_Three_Step_Full),
                            .f = function(x, y){round(t(apply(y - x, 2, quantile, p = probs)), 1)})

CIs_out <- rep(list(matrix(NA, nrow = d, ncol = 10)), length(n_excesses))
for(i in 1:length(n_excesses)){
  for(j in 1:d){
    CIs_out[[i]][j,] <- c(j,
                          paste0(CIs_One_Step_Indep[[i]][j,1], " (", CIs_One_Step_Indep[[i]][j,2], ", ", CIs_One_Step_Indep[[i]][j,3], ")"),
                          paste0(CIs_Two_Step_Indep[[i]][j,1], " (", CIs_Two_Step_Indep[[i]][j,2], ", ", CIs_Two_Step_Indep[[i]][j,3], ")"),
                          paste0(CIs_Three_Step_Indep[[i]][j,1], " (", CIs_Three_Step_Indep[[i]][j,2], ", ", CIs_Three_Step_Indep[[i]][j,3], ")"),
                          paste0(CIs_One_Step_Graph[[i]][j,1], " (", CIs_One_Step_Graph[[i]][j,2], ", ", CIs_One_Step_Graph[[i]][j,3], ")"),
                          paste0(CIs_Two_Step_Graph[[i]][j,1], " (", CIs_Two_Step_Graph[[i]][j,2], ", ", CIs_Two_Step_Graph[[i]][j,3], ")"),
                          paste0(CIs_Three_Step_Graph[[i]][j,1], " (", CIs_Three_Step_Graph[[i]][j,2], ", ", CIs_Three_Step_Graph[[i]][j,3], ")"),
                          paste0(CIs_One_Step_Full[[i]][j,1], " (", CIs_One_Step_Full[[i]][j,2], ", ", CIs_One_Step_Full[[i]][j,3], ")"),
                          paste0(CIs_Two_Step_Full[[i]][j,1], " (", CIs_Two_Step_Full[[i]][j,2], ", ", CIs_Two_Step_Full[[i]][j,3], ")"),
                          paste0(CIs_Three_Step_Full[[i]][j,1], " (", CIs_Three_Step_Full[[i]][j,2], ", ", CIs_Three_Step_Full[[i]][j,3], ")"))
  }
}

col_names <- c("Conditioning Variable", rep(c("One-Step", "Two-Step", "Three-Step"), 3))
CIs_out[[1]] %>% 
  kbl(format = "latex",
      col.names = col_names, align = "c") %>% 
  kable_classic(full_width = F, html_font = "Source Sans Pro")

CIs_out[[2]] %>% 
  kbl(format = "latex",
      col.names = col_names, align = "c") %>% 
  kable_classic(full_width = F, html_font = "Source Sans Pro")
