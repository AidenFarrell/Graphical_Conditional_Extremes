################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("fake", "ggpattern", "ggplot2", "glasso", "gridExtra", "gtools", "igraph", "jmuOutlier", "LaplacesDemon", "matrixcalc", "mev", "moments", "rlang", "parallel", "pracma", "purrr", "texmex")
# install.packages(required_pckgs, dependencies = TRUE)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
setwd("/home/farrel11/Documents/Section_4_1/Par_Ests/Low_Dependence/")

################################################################################
## Load in Rs script with functions for fitting the conditional extremes model with

## a graphical structure for the residuals
# source("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/Conditional_Extremes/231114_Cond_Extremes_Graphical_Residuals.R")
source("/home/farrel11/PhD/Project_2/Conditional_Extremes/231114_Cond_Extremes_Graphical_Residuals.R")

## Files for the Murphy et al. (2023) threshold selection method
# source("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/threshold_selection_paper/helper_functions.R")
# source("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/threshold_selection_paper/thresh_qq_metric.R")
source("/home/farrel11/PhD/Project_2/threshold_selection_paper/helper_functions.R")
source("/home/farrel11/PhD/Project_2/threshold_selection_paper/thresh_qq_metric.R")

# source("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/Conditional_Extremes/General_Functions.R")
source("/home/farrel11/PhD/Project_2/Conditional_Extremes/General_Functions.R")

# source("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/Conditional_Extremes/MVAGG_Functions.R")
source("/home/farrel11/PhD/Project_2/Conditional_Extremes/MVAGG_Functions.R")

################################################################################
## Colour blind friendly pallete:
cbbPalette <- c(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                palette.colors(palette = "Dark 2"))

################################################################################
## Comparison plots for low and high dependence structures
comp_plots <- function(data, methods, y_lab, par_true){
  
  ## Obtain some infromation from the data
  if(!is.list(data)){
    stop("data must be a list of parameter estimates to plot")
  }
  n_methods <- length(methods)
  d_data <- length(data)
  n <- nrow(data[[1]][[1]][[1]])
  p <- ncol(data[[1]][[1]][[1]])
  
  ## Get the true parameters in the correct form
  par_true_NA <- lapply(1:d, function(i){Add_NA_Vector(par_true[-i], i)})
  
  ## get the bias
  data_bias <- lapply(1:d_data, function(i){lapply(1:2, function(j){lapply(1:n_methods, function(k){
    data[[i]][[j]][[k]] - matrix(rep(par_true_NA[[i]], n), nrow = n, byrow = TRUE)
  })})})
  
  ## Construct the plotting data
  plot_data <- data.frame(y = do.call(c, lapply(1:d_data, function(i){
    do.call(c, lapply(1:p, function(j){
      do.call(c, lapply(1:2, function(k){
        do.call(c, lapply(1:n_methods, function(l){
          data_bias[[i]][[k]][[l]][,j]}))}))}))})))
  
  plot_data$Method = rep(rep(rep(rep(methods, each = n), 2), p), d_data)
  plot_data$Method <- factor(plot_data$Method, levels = methods)
  plot_data$Conditioning_Varaible <- rep(1:d_data, each = 2*p*n*n_methods)
  plot_data$Conditioning_Varaible <- factor(plot_data$Conditioning_Varaible, levels = 1:d_data)
  plot_data$Dependent_Varaible <- rep(rep(1:p, each = 2*n*n_methods), d_data)
  plot_data$Dependent_Varaible <- factor(plot_data$Dependent_Varaible, levels = 1:p)
  plot_data$Type <- rep(rep(rep(c("250", "500"), each = n*n_methods), p), d_data)
  plot_data$Type <- factor(plot_data$Type, levels = c("250", "500"))
  
  ## Order the data by Pair and Type
  plot_data <- plot_data[order(plot_data$Dependent_Varaible, plot_data$Type), ]
  
  ## Make the plot
  plot_out <- ggplot(data = plot_data, aes(x = Dependent_Varaible, y = y, fill = Method, pattern = Type)) + 
    geom_boxplot_pattern(position = position_dodge(preserve = "single"),
                         color = "black", pattern_fill = "black",
                         pattern_angle = 45, pattern_density = 0.1,
                         pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    scale_pattern_manual(values = c("250" = "none", "500" = "stripe")) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme(legend.position = "top") +
    labs(x = "Dependent Variable", y = y_lab, pattern = "Number of Excesses") +
    facet_grid(rows = vars(Conditioning_Varaible)) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.5)
  
  return(plot_out)
}

## Comparison plots for low and high dependence structures
comp_plots_matrix <- function(data, methods, y_lab, cov_mat_true, precision = FALSE){
  
  ## Obtain some infromation from the data
  if(!is.list(data)){
    stop("data must be a list of parameter estimates to plot")
  }
  n_methods <- length(methods)
  d_data <- length(data)
  n <- nrow(data[[1]][[1]][[1]])
  p <- ncol(data[[1]][[1]][[1]])
  
  ## get some plotting parameters
  x_labels <- apply(combinations(n = d_data, r = 2, v = 1:d_data, repeats.allowed = TRUE), 1, paste0, collapse = "")
  
  ## Get the true parameters in the correct form
  lower_tri_elements <- lower.tri(cov_mat_true, diag = TRUE)
  cov_mat_true_i <- lapply(1:d_data, function(i){Cond_Sigma(cov_mat_true, i)})
  cor_mat_true_i <- lapply(cov_mat_true_i, Sigma2rho)
  if(precision){
    cor_mat_true_i <- lapply(cor_mat_true_i, solve)
  }
  cor_mat_true_i_NA <- lapply(1:d_data, function(i){Add_NA_Matrix(cor_mat_true_i[[i]], i)[lower_tri_elements]})
  
  ## get the bias
  data_bias <- lapply(1:d_data, function(i){lapply(1:2, function(j){lapply(1:n_methods, function(k){
    data[[i]][[j]][[k]] - matrix(rep(cor_mat_true_i_NA[[i]], n), nrow = n, byrow = TRUE)
  })})})
  
  ## Construct the plotting data
  plot_data <- data.frame(y = do.call(c, lapply(1:d_data, function(i){
    do.call(c, lapply(1:p, function(j){
      do.call(c, lapply(1:2, function(k){
        do.call(c, lapply(1:n_methods, function(l){
          data_bias[[i]][[k]][[l]][,j]}))}))}))})))
  
  plot_data$Method = rep(rep(rep(rep(methods, each = n), 2), p), d_data)
  plot_data$Method <- factor(plot_data$Method, levels = methods)
  plot_data$Conditioning_Varaible <- rep(1:d_data, each = 2*p*n*n_methods)
  plot_data$Conditioning_Varaible <- factor(plot_data$Conditioning_Varaible, levels = 1:d_data)
  plot_data$Pair <- rep(rep(x_labels, each = 2*n*n_methods), d_data)
  plot_data$Pair <- factor(plot_data$Pair, levels = x_labels)
  plot_data$Type <- rep(rep(rep(c("250", "500"), each = n*n_methods), p), d_data)
  plot_data$Type <- factor(plot_data$Type, levels = c("250", "500"))
  
  ## Make the plot
  plot_out <- ggplot(data = plot_data, aes(x = Pair, y = y, fill = Method, pattern = Type)) + 
    geom_boxplot_pattern(position = position_dodge(preserve = "single"),
                         color = "black", pattern_fill = "black",
                         pattern_angle = 45, pattern_density = 0.1,
                         pattern_spacing = 0.01, pattern_key_scale_factor = 0.6) +
    scale_pattern_manual(values = c("250" = "none", "500" = "stripe")) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme(legend.position = "top") +
    labs(x = "Pair", y = y_lab, pattern = "Number of Excesses") +
    facet_grid(rows = vars(Conditioning_Varaible)) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.5)
  
  return(plot_out)
}

################################################################################
## DO NOT RUN

## True graph
d <- 5
sep <- d1 <- ceiling(d/2)
edges <- cbind(t(combinations(n = d1, r = 2, v = 1:d1)), t(combinations(n = d1, r = 2, v = d1:d)))
# edges <- t(rbind(c(1,2), c(2,3), c(2,4), c(2,5), c(3,4), c(3,5), c(4,5), c(4,6), c(5,6)))
g_true <- graph(edges = edges, directed = FALSE)
plot.igraph(g_true)

## True parameters in Z
n_par <- d*(d-1)/2
alpha_true <- runif(n = d, 0.1, 0.5)
beta_true <- runif(n = d, 0.1, 0.3)
loc_true <- runif(n = d, -5, 5)
scale_1_true <- runif(n = d, 0.5, 2)
scale_2_true <- runif(n = d, 1.5, 3)
shape_true <- runif(n = d, 0.8, 2.5)

# alpha_true <- c(0.142255902010947, 0.231802942370996, 0.241576673928648, 0.284543012687936, 0.251147609343752)
# beta_true <- c(0.164718892658129, 0.29614326544106, 0.147650257870555, 0.147218536119908, 0.1729731567204)
# loc_true <- c(-2.28005731012672, 4.38922758214176, 0.479259195271879, -4.08727945527062, -0.409800345078111)
# scale_1_true <- c(1.86106580833439, 1.99781919436064, 1.18922822321765, 0.933195473253727, 0.561478922399692)
# scale_2_true <- c(0.562924829428084, 1.58315827895422, 1.59529121881351, 1.41660505274776, 1.68265633923002)
# shape_true <- c(1.3071317654103, 2.36861560784746, 0.829040728230029, 1.38464096179232, 1.85569009417668)

simul <- SimulatePrecision(theta = as.matrix(as_adjacency_matrix(g_true)), 
                           v_sign = -1, v_within = c(0.1, 0.6))
Gamma_true <- simul$omega
Sigma_true <- solve(Gamma_true)
rho_true <- Sigma2rho(Sigma_true)
rho_true

Sigma_true_i <- lapply(1:d, function(i){Cond_Sigma(Sigma_true, i)})
rho_true_i <- lapply(Sigma_true_i, function(x){Sigma2rho(x)})

## Simulate large Yi
n_sim <- 200
dqu <- 0.8
n_excesses <- c(250, 500)

Yi_large <- vector("list", length = length(n_excesses))
for(i in 1:length(n_excesses)){
  Yi_large[[i]] <- replicate(n = d,
                             expr = replicate(n = n_sim, expr = unif_to_laplace(runif(n_excesses[i], dqu, 1)), simplify = FALSE),
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
    Z_not_i[[i]][[j]] <- replicate(n = n_sim, expr = rmvagg(n = n_excesses[i], dim = d-1, 
                                                            loc = loc_true[-j],
                                                            scale_1 = scale_1_true[-j],
                                                            scale_2 = scale_2_true[-j],
                                                            shape = shape_true[-j],
                                                            Sigma = rho_true_i[[j]]),
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
## DO NOT RUN

## Heffernan and Tawn model
fit_HT <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_graph_new,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = NA,
                             cond = j,
                             v = ceiling(max(sapply(Y_Yi_large, function(x){sapply(x, function(y){sapply(y, max)})}))) + 1,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

Z_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){lapply(y, function(z){z$Z})})})
a_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){lapply(y, function(z){z$par$main[1,]})})})
b_hat_HT <- lapply(fit_HT, function(x){lapply(x, function(y){lapply(y, function(z){z$par$main[2,]})})})

## One-step model fits
## One-step model is sensitive to the starting parameter of the the location parameter in the AGG
## This is heavily correlated with the alpha dependence parameter
## Give an "informed" start for the location parameter in the AGG

loc_start_One_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(1:n_sim, function(k){
    apply((Y_Yi_large[[i]][[j]][[k]][,-j] - matrix(0.1*Y_Yi_large[[i]][[j]][[k]][,j], nrow = n_excesses[i], ncol = d-1))/
            matrix(Y_Yi_large[[i]][[j]][[k]][,j]^0.1, nrow = n_excesses[i], ncol = d-1), 2, mean)
  })})})
start_par_One_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(1:n_sim, function(k){
    cbind(rep(0.1, d-1), rep(0.1, d-1), loc_start_One_Step[[i]][[j]][[k]],
          rep(1.5, d-1), rep(1.5, d-1), rep(1.5, d-1))})})})

## Independence
fit_One_Step_Indep <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Gamma_new,
             data = Y_Yi_large[[i]][[j]],
             start = start_par_One_Step[[i]][[j]],
             MoreArgs = list(graph = NA,
                             cond = j,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## give some informed starting place to the graphical and saturated models
start_par_One_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(1:n_sim, function(k){
    cbind(pmin(pmax(unname(unlist(fit_One_Step_Indep[[i]][[j]][[k]]$par$main[1,])), 0.1), 0.5),
          pmin(pmax(unname(unlist(fit_One_Step_Indep[[i]][[j]][[k]]$par$main[2,])), 0.1), 0.5),
          unname(unlist(fit_One_Step_Indep[[i]][[j]][[k]]$par$main[3,])),
          pmin(pmax(unname(unlist(fit_One_Step_Indep[[i]][[j]][[k]]$par$main[4,])), 0.5), 2),
          pmin(pmax(unname(unlist(fit_One_Step_Indep[[i]][[j]][[k]]$par$main[5,])), 0.5), 2),
          pmin(pmax(unname(unlist(fit_One_Step_Indep[[i]][[j]][[k]]$par$main[6,])), 0.5), 2))})})})

## Graphical
fit_One_Step_Graph <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Gamma_new,
             data = Y_Yi_large[[i]][[j]],
             start = start_par_One_Step[[i]][[j]],
             MoreArgs = list(graph = g_true,
                             cond = j,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Saturated
fit_One_Step_Full <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Gamma_new,
             data = Y_Yi_large[[i]][[j]],
             start = start_par_One_Step[[i]][[j]],
             MoreArgs = list(graph = make_full_graph(n = d),
                             cond = j,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})


## Two-step model fits
## Again, give an informed start to the location parameter
loc_start_Two_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(1:n_sim, function(k){
    apply(Z_hat_HT[[i]][[j]][[k]], 2, mean)})})})
start_par_Two_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(1:n_sim, function(k){
    cbind(loc_start_Two_Step[[i]][[j]][[k]], rep(1.5, d-1), rep(1.5, d-1), rep(1.5, d-1))})})})

## Independence
fit_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Gamma_TS_new_Z,
             z = Z_hat_HT[[i]][[j]],
             a = a_hat_HT[[i]][[j]],
             b = b_hat_HT[[i]][[j]],
             start = start_par_Two_Step[[i]][[j]],
             MoreArgs = list(graph = NA,
                             cond = j,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Give more informed starting positions for the graphical and saturated models
start_par_Two_Step <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){lapply(1:n_sim, function(k){
    cbind(unname(unlist(fit_Two_Step_Indep[[i]][[j]][[k]]$par$main[3,])),
          pmin(pmax(unname(unlist(fit_Two_Step_Indep[[i]][[j]][[k]]$par$main[4,])), 0.5), 2),
          pmin(pmax(unname(unlist(fit_Two_Step_Indep[[i]][[j]][[k]]$par$main[5,])), 0.5), 2),
          pmin(pmax(unname(unlist(fit_Two_Step_Indep[[i]][[j]][[k]]$par$main[6,])), 0.5), 2))})})})

## Graphical
fit_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Gamma_TS_new_Z,
             z = Z_hat_HT[[i]][[j]],
             a = a_hat_HT[[i]][[j]],
             b = b_hat_HT[[i]][[j]],
             start = start_par_Two_Step[[i]][[j]],
             MoreArgs = list(graph = g_true,
                             cond = j,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Saturated
fit_Two_Step_Full <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Gamma_TS_new_Z,
             z = Z_hat_HT[[i]][[j]],
             a = a_hat_HT[[i]][[j]],
             b = b_hat_HT[[i]][[j]],
             start = start_par_Two_Step[[i]][[j]],
             MoreArgs = list(graph = make_full_graph(n = d),
                             cond = j,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Update the function to take starting parameter arguments
## Then set location parameter to mean of fitted Z values

## Three-step model fits
## Independence
fit_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Three_step,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = NA,
                             cond = j,
                             v = ceiling(max(sapply(Y_Yi_large, function(x){sapply(x, function(y){sapply(y, max)})}))) + 1,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Graphical
fit_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Three_step,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = g_true,
                             cond = j,
                             v = ceiling(max(sapply(Y_Yi_large, function(x){sapply(x, function(y){sapply(y, max)})}))) + 1,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

## Saturated
fit_Three_Step_Full <- lapply(1:length(n_excesses), function(i){
  lapply(1:d, function(j){
    mcmapply(FUN = Cond_extremes_MVAGGD_Three_step,
             data = Y_Yi_large[[i]][[j]],
             MoreArgs = list(graph = make_full_graph(n = d),
                             cond = j,
                             v = ceiling(max(sapply(Y_Yi_large, function(x){sapply(x, function(y){sapply(y, max)})}))) + 1,
                             maxit = 1e+9),
             SIMPLIFY = FALSE,
             mc.cores = detectCores() - 1)})})

################################################################################
## DO NOT RUN

## One and two step models can be somewhat sensitive to the starting values
## here we need to do another procedure to ensure the models converge

Index_One_Step_Indep <- lapply(fit_One_Step_Indep, function(x){
  lapply(x, function(y){which(sapply(y, is.character))})})
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
            tryCatch(Cond_extremes_MVAGGD_Gamma_new(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                                    cond = j,
                                                    graph = NA,
                                                    start = start_par_One_Step,
                                                    maxit = 1e+9),
                     error = function(e){"Error in model fitting"})
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
    lapply(x, function(y){which(sapply(y, is.character))})})
  if(any(sapply(Index_One_Step_Indep, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}


Index_One_Step_Graph <- lapply(fit_One_Step_Graph, function(x){
  lapply(x, function(y){which(sapply(y, is.character))})})
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
            tryCatch(Cond_extremes_MVAGGD_Gamma_new(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                                    cond = j,
                                                    graph = g_true,
                                                    start = start_par_One_Step,
                                                    maxit = 1e+9),
                     error = function(e){"Error in model fitting"})
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
    lapply(x, function(y){which(sapply(y, is.character))})})
  if(any(sapply(Index_One_Step_Graph, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_One_Step_Full <- lapply(fit_One_Step_Full, function(x){
  lapply(x, function(y){which(sapply(y, is.character))})})
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
            tryCatch(Cond_extremes_MVAGGD_Gamma_new(data = Y_Yi_large[[i]][[j]][[ind[k]]],
                                                    cond = j,
                                                    graph = make_full_graph(n = d),
                                                    start = start_par_One_Step,
                                                    maxit = 1e+9),
                     error = function(e){"Error in model fitting"})
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
    lapply(x, function(y){which(sapply(y, is.character))})})
  if(any(sapply(Index_One_Step_Full, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_Two_Step_Indep <- lapply(fit_Two_Step_Indep, function(x){
  lapply(x, function(y){which(sapply(y, is.character))})})
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
            tryCatch(Cond_extremes_MVAGGD_Gamma_TS_new_Z(z = Z_hat_HT[[i]][[j]][[ind[k]]],
                                                         a = a_hat_HT[[i]][[j]][[ind[k]]],
                                                         b = b_hat_HT[[i]][[j]][[ind[k]]],
                                                         cond = j,
                                                         graph = NA,
                                                         start = start_par_Two_Step,
                                                         maxit = 1e+9),
                     error = function(e){"Error in model fitting"})
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
    lapply(x, function(y){which(sapply(y, is.character))})})
  if(any(sapply(Index_Two_Step_Indep, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_Two_Step_Graph <- lapply(fit_Two_Step_Graph, function(x){
  lapply(x, function(y){which(sapply(y, is.character))})})
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
            tryCatch(Cond_extremes_MVAGGD_Gamma_TS_new_Z(z = Z_hat_HT[[i]][[j]][[ind[k]]],
                                                         a = a_hat_HT[[i]][[j]][[ind[k]]],
                                                         b = b_hat_HT[[i]][[j]][[ind[k]]],
                                                         cond = j,
                                                         graph = g_true,
                                                         start = start_par_Two_Step,
                                                         maxit = 1e+9),
                     error = function(e){"Error in model fitting"})
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
    lapply(x, function(y){which(sapply(y, is.character))})})
  if(any(sapply(Index_Two_Step_Graph, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}

Index_Two_Step_Full <- lapply(fit_Two_Step_Full, function(x){
  lapply(x, function(y){which(sapply(y, is.character))})})
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
            tryCatch(Cond_extremes_MVAGGD_Gamma_TS_new_Z(z = Z_hat_HT[[i]][[j]][[ind[k]]],
                                                         a = a_hat_HT[[i]][[j]][[ind[k]]],
                                                         b = b_hat_HT[[i]][[j]][[ind[k]]],
                                                         cond = j,
                                                         graph = make_full_graph(n = d),
                                                         start = start_par_Two_Step,
                                                         maxit = 1e+9),
                     error = function(e){"Error in model fitting"})
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
    lapply(x, function(y){which(sapply(y, is.character))})})
  if(any(sapply(Index_Two_Step_Full, function(x){sapply(x, any)})) == FALSE){
    print("Model Fitting Complete")
  }
}


################################################################################
## set up the simulation study
out <- readRDS("/home/farrel11/Documents/Data/Section_4_1/True_Dist_Low_Dep.RData")
# out <- readRDS("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/Code_for_Paper/Data/Section_4_1/True_Dist_Low_Dep.RData")
g_true <- out$par_true$graph
d <- length(V(g_true))

alpha_true <- out$par_true$alpha
beta_true <- out$par_true$beta
loc_true <- out$par_true$location
scale_1_true <- out$par_true$scale_left
scale_2_true <- out$par_true$scale_right
shape_true <- out$par_true$shape

Gamma_true <- out$par_true$Gamma
Sigma_true <- solve(Gamma_true)
rho_true <- Sigma2rho(Sigma_true)
rho_true

Sigma_true_i <- lapply(1:d, function(i){Cond_Sigma(Sigma_true, i)})
rho_true_i <- lapply(Sigma_true_i, function(x){Sigma2rho(x)})

## Simulate large Yi
n_sim <- out$par_true$n_sim
dqu <- out$par_true$dqu
n_excesses <- out$par_true$n_excesses
Y_Yi_large <- out$data$Y_Yi_large

fit_HT <- out$Model_out$HT
fit_One_Step_Indep <- out$Model_out$One_Step_Indep
fit_One_Step_Graph <- out$Model_out$One_Step_Graph
fit_One_Step_Full <- out$Model_out$One_Step_Full
fit_Two_Step_Indep <- out$Model_out$Two_Step_Indep
fit_Two_Step_Graph <- out$Model_out$Two_Step_Graph
fit_Two_Step_Full <- out$Model_out$Two_Step_Full
fit_Three_Step_Indep <- out$Model_out$Three_Step_Indep
fit_Three_Step_Graph <- out$Model_out$Three_Step_Graph
fit_Three_Step_Full <- out$Model_out$Three_Step_Full

################################################################################
## Extract the parameter estimates
lower_tri_elements <- lower.tri(Gamma_true, diag = TRUE)

## One-step Independent Residuals
a_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(solve(y$par$Gamma), j)[lower_tri_elements]}))})})
Gamma_hat_One_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

## One-step Graphical residuals
a_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(solve(y$par$Gamma), j)[lower_tri_elements]}))})})
Gamma_hat_One_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

## One-step Saturated Residuals
a_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(solve(y$par$Gamma), j)[lower_tri_elements]}))})})
Gamma_hat_One_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_One_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

## Two-step Independent Residuals
a_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(solve(y$par$Gamma), j)[lower_tri_elements]}))})})
Gamma_hat_Two_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

## Two-step Graphical residuals
a_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(solve(y$par$Gamma), j)[lower_tri_elements]}))})})
Gamma_hat_Two_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

## Two-step Saturated Residuals
a_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(solve(y$par$Gamma), j)[lower_tri_elements]}))})})
Gamma_hat_Two_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Two_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

## Three-Step Independent Residuals
a_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(y$par$Sigma, j)[lower_tri_elements]}))})})
Gamma_hat_Three_Step_Indep <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Indep[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

## Three-Step Graphical residuals
a_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(y$par$Sigma, j)[lower_tri_elements]}))})})
Gamma_hat_Three_Step_Graph <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Graph[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

## Three-Step Saturated Residuals
a_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[1,])), j)}))})})
b_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[2,])), j)}))})})
loc_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[3,])), j)}))})})
scale_1_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[4,])), j)}))})})
scale_2_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[5,])), j)}))})})
shape_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Vector(unname(unlist(y$par$main[6,])), j)}))})})
Sigma_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(y$par$Sigma, j)[lower_tri_elements]}))})})
Gamma_hat_Three_Step_Full <- lapply(1:length(n_excesses), function(i){lapply(1:d, function(j){
  t(sapply(fit_Three_Step_Full[[i]][[j]], function(y){Add_NA_Matrix(y$par$Gamma, j)[lower_tri_elements]}))})})

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
  methods = method_vec, y_lab = expression("Bias in" ~ hat(alpha)), par_true = alpha_true)
pdf("Alpha.pdf", width = 15, height = 10)
print(Alpha_plot)
dev.off()

## Beta
Beta_plot <- comp_plots(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(b_hat_One_Step_Indep[[i]][[j]], 
         b_hat_One_Step_Graph[[i]][[j]], 
         b_hat_One_Step_Full[[i]][[j]],
         b_hat_Two_Step_Indep[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(beta)), par_true = beta_true)
pdf("Beta.pdf", width = 15, height = 10)
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
  methods = method_vec, y_lab = expression("Bias in" ~ hat(nu)), par_true = loc_true)
pdf("Location.pdf", width = 15, height = 10)
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
  methods = method_vec, y_lab = expression("Bias in" ~ hat(kappa[1])), par_true = scale_1_true)
pdf("Scale_Left.pdf", width = 15, height = 10)
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
  methods = method_vec, y_lab = expression("Bias in" ~ hat(kappa[2])), par_true = scale_2_true)
pdf("Scale_Right.pdf", width = 15, height = 10)
print(Scale_2_plot)
dev.off()

## Shape
shape_plot <- comp_plots(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(shape_hat_One_Step_Indep[[i]][[j]], 
         shape_hat_One_Step_Graph[[i]][[j]], 
         shape_hat_One_Step_Full[[i]][[j]],
         shape_hat_Two_Step_Indep[[i]][[j]], 
         shape_hat_Two_Step_Graph[[i]][[j]], 
         shape_hat_Two_Step_Full[[i]][[j]],
         shape_hat_Three_Step_Indep[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(delta)), par_true = shape_true)
pdf("Shape.pdf", width = 15, height = 10)
print(shape_plot)
dev.off()

method_vec <- c("One-step - Independence", "One-step - Graphical", "One-step - Saturated", 
                "Two-step - Independence", "Two-step - Graphical", "Two-step - Saturated",
                "Three-step - Independence", "Three-step - Graphical", "Three-step - Saturated")

## Correlation Matrix
Sigma_plot <- comp_plots_matrix(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(Sigma_hat_One_Step_Indep[[i]][[j]], Sigma_hat_One_Step_Graph[[i]][[j]], Sigma_hat_One_Step_Full[[i]][[j]],
         Sigma_hat_Two_Step_Indep[[i]][[j]], Sigma_hat_Two_Step_Graph[[i]][[j]], Sigma_hat_Two_Step_Full[[i]][[j]],
         Sigma_hat_Three_Step_Indep[[i]][[j]], Sigma_hat_Three_Step_Graph[[i]][[j]], Sigma_hat_Three_Step_Full[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(Sigma)), cov_mat_true = Sigma_true, precision = FALSE)

pdf("Sigma.pdf", width = 20, height = 15)
print(Sigma_plot)
dev.off()

## Precision Matrix
Gamma_plot <- comp_plots_matrix(data = lapply(1:d, function(j){
  lapply(1:length(n_excesses), function(i){
    list(Gamma_hat_One_Step_Indep[[i]][[j]], Gamma_hat_One_Step_Graph[[i]][[j]], Gamma_hat_One_Step_Full[[i]][[j]],
         Gamma_hat_Two_Step_Indep[[i]][[j]], Gamma_hat_Two_Step_Graph[[i]][[j]], Gamma_hat_Two_Step_Full[[i]][[j]],
         Gamma_hat_Three_Step_Indep[[i]][[j]], Gamma_hat_Three_Step_Graph[[i]][[j]], Gamma_hat_Three_Step_Full[[i]][[j]])})}),
  methods = method_vec, y_lab = expression("Bias in" ~ hat(Gamma)), cov_mat_true = Sigma_true, precision = TRUE)

pdf("Gamma.pdf", width = 20, height = 15)
print(Gamma_plot)
dev.off()

################################################################################
p <- c(0.5, 0.025, 0.975)

CI_Latex <- function(x, p = c(0.5, 0.025, 0.975), dp = 2){
  if(length(p) != 3 | any(abs(p) > 1)){
    stop("p must be a vector of probabilities of length 3")
  }
  if(dp <= 0 | length(dp) > 1 | dp%%1 != 0){
    stop("dp must be a single positive integer represneting the number of decimial places to print to")
  }
  if(all(is.na(x))){
    paste0("")
  }
  else{
    paste0(format(round(quantile(x, probs = p[1], na.rm = TRUE), dp), nsmall = dp), 
           " (",
           format(round(quantile(x, probs = p[2], na.rm = TRUE), dp), nsmall = dp),
           ", ",
           format(round(quantile(x, probs = p[3], na.rm = TRUE), dp), nsmall = dp),
           ")")  
  }
} 

lapply(1:d, function(i){
  paste(c(paste(c(paste(c("", method_vec[1], paste(sapply(1:d, function(j){CI_Latex(a_hat_One_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[2], paste(sapply(1:d, function(j){CI_Latex(a_hat_One_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[3], paste(sapply(1:d, function(j){CI_Latex(a_hat_One_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[4], paste(sapply(1:d, function(j){CI_Latex(a_hat_Two_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " ")), collapse = "")
})

lapply(1:d, function(i){
  paste(c(paste(c(paste(c("", method_vec[1], paste(sapply(1:d, function(j){CI_Latex(b_hat_One_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[2], paste(sapply(1:d, function(j){CI_Latex(b_hat_One_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[3], paste(sapply(1:d, function(j){CI_Latex(b_hat_One_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[4], paste(sapply(1:d, function(j){CI_Latex(b_hat_Two_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " ")), collapse = "")
})

lapply(1:d, function(i){
  paste(c(paste(c(paste(c("", method_vec[1], paste(sapply(1:d, function(j){CI_Latex(loc_hat_One_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[2], paste(sapply(1:d, function(j){CI_Latex(loc_hat_One_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[3], paste(sapply(1:d, function(j){CI_Latex(loc_hat_One_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[4], paste(sapply(1:d, function(j){CI_Latex(loc_hat_Two_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[5], paste(sapply(1:d, function(j){CI_Latex(loc_hat_Two_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[6], paste(sapply(1:d, function(j){CI_Latex(loc_hat_Two_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[7], paste(sapply(1:d, function(j){CI_Latex(loc_hat_Three_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " ")), collapse = "")
})

lapply(1:d, function(i){
  paste(c(paste(c(paste(c("", method_vec[1], paste(sapply(1:d, function(j){CI_Latex(scale_1_hat_One_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[2], paste(sapply(1:d, function(j){CI_Latex(scale_1_hat_One_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[3], paste(sapply(1:d, function(j){CI_Latex(scale_1_hat_One_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[4], paste(sapply(1:d, function(j){CI_Latex(scale_1_hat_Two_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[5], paste(sapply(1:d, function(j){CI_Latex(scale_1_hat_Two_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[6], paste(sapply(1:d, function(j){CI_Latex(scale_1_hat_Two_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[7], paste(sapply(1:d, function(j){CI_Latex(scale_1_hat_Three_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " ")), collapse = "")
})

lapply(1:d, function(i){
  paste(c(paste(c(paste(c("", method_vec[1], paste(sapply(1:d, function(j){CI_Latex(scale_2_hat_One_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[2], paste(sapply(1:d, function(j){CI_Latex(scale_2_hat_One_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[3], paste(sapply(1:d, function(j){CI_Latex(scale_2_hat_One_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[4], paste(sapply(1:d, function(j){CI_Latex(scale_2_hat_Two_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[5], paste(sapply(1:d, function(j){CI_Latex(scale_2_hat_Two_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[6], paste(sapply(1:d, function(j){CI_Latex(scale_2_hat_Two_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[7], paste(sapply(1:d, function(j){CI_Latex(scale_2_hat_Three_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " ")), collapse = "")
})


lapply(1:d, function(i){
  paste(c(paste(c(paste(c("", method_vec[1], paste(sapply(1:d, function(j){CI_Latex(shape_hat_One_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[2], paste(sapply(1:d, function(j){CI_Latex(shape_hat_One_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[3], paste(sapply(1:d, function(j){CI_Latex(shape_hat_One_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[4], paste(sapply(1:d, function(j){CI_Latex(shape_hat_Two_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[5], paste(sapply(1:d, function(j){CI_Latex(shape_hat_Two_Step_Graph[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[6], paste(sapply(1:d, function(j){CI_Latex(shape_hat_Two_Step_Full[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " "),
          paste(c(paste(c("", method_vec[7], paste(sapply(1:d, function(j){CI_Latex(shape_hat_Three_Step_Indep[[1]][[i]][,j])}), collapse = " & ")), collapse = " & "), "\\"), collapse = " ")), collapse = "")
})


################################################################################
## Save the output
# out <- list(data = list(Y_Yi_large = Y_Yi_large),
#             par_true = list(n_sim = n_sim, n_excesses = n_excesses, dqu = dqu,
#                             graph = g_true,
#                             alpha = alpha_true,
#                             beta = beta_true,
#                             location = loc_true,
#                             scale_left = scale_1_true,
#                             scale_right = scale_2_true,
#                             shape = shape_true,
#                             Gamma = Gamma_true),
#             Model_out = list(HT = fit_HT,
#                              One_Step_Indep = fit_One_Step_Indep,
#                              One_Step_Graph = fit_One_Step_Graph,
#                              One_Step_Full = fit_One_Step_Full,
#                              Two_Step_Indep = fit_Two_Step_Indep,
#                              Two_Step_Graph = fit_Two_Step_Graph,
#                              Two_Step_Full = fit_Two_Step_Full,
#                              Three_Step_Indep = fit_Three_Step_Indep,
#                              Three_Step_Graph = fit_Three_Step_Graph,
#                              Three_Step_Full = fit_Three_Step_Full))
# 
# saveRDS(out, file = "/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/Code_for_Paper/Data/Section_4_1/True_Dist_Low_Dep.RData")
# out <- readRDS("/Users/aidenfarrell/Library/CloudStorage/OneDrive-LancasterUniversity/PhD/Project_2/Code_for_Paper/Data/Section_4_1/True_Dist_Low_Dep.RData")
################################################################################