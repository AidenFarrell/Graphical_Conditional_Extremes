################################################################################
## Reading in required packages
required_pckgs <- c("ggpattern", "ggplot2", "gtools")
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################

## Comparison plots parameters
comp_plots <- function(data, methods, y_lab, ylims, par_true){
  
  ## Obtain some information from the data
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
  
  ## Define custom labels for facets
  facet_labels <- setNames(paste0("i = ", 1:d_data), 1:d_data)
  
  ## Make the plot
  plot_out <- ggplot(data = plot_data, aes(x = Dependent_Varaible, y = y, fill = Method, pattern = Type)) + 
    geom_boxplot_pattern(position = position_dodge(preserve = "single"),
                         color = "black", pattern_fill = "black",
                         pattern_angle = 45, pattern_density = 0.1,
                         pattern_spacing = 0.025, pattern_key_scale_factor = 0.6,
                         outlier.shape = NA) +
    scale_pattern_manual(values = c("250" = "none", "500" = "stripe")) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme(legend.position = "top",
          legend.box = "horizontal",
          legend.box.just = "center",
          legend.spacing.x = unit(0.5, 'cm'),
          legend.spacing.y = unit(0.5, 'cm'),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0)) +
    lims(y = ylims) +
    labs(x = "Dependent Variable (j)", y = y_lab, pattern = "Number of Excesses") +
    facet_grid(rows = vars(Conditioning_Varaible), labeller = labeller(Conditioning_Varaible = facet_labels)) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.5)
  
  return(plot_out)
}

## Comparison plots for covaraince/precision matrices
comp_plots_matrix <- function(data, methods, y_lab, ylims, cov_mat_true, precision = FALSE){
  
  ## Obtain some information from the data
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
  cor_mat_true_i <- lapply(cov_mat_true_i, cov2cor)
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
  
  ## Define custom labels for facets
  facet_labels <- setNames(paste0("i = ", 1:d_data), 1:d_data)
  
  ## Make the plot
  plot_out <- ggplot(data = plot_data, aes(x = Pair, y = y, fill = Method, pattern = Type)) + 
    geom_boxplot_pattern(position = position_dodge(preserve = "single"),
                         color = "black", pattern_fill = "black",
                         pattern_angle = 45, pattern_density = 0.1,
                         pattern_spacing = 0.01, pattern_key_scale_factor = 0.6,
                         outlier.shape = NA) +
    scale_pattern_manual(values = c("250" = "none", "500" = "stripe")) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme(legend.position = "top",
          legend.box = "horizontal",
          legend.box.just = "center",
          legend.spacing.x = unit(0.5, 'cm'),
          legend.spacing.y = unit(0.5, 'cm'),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0)) +
    lims(y = ylims) +
    labs(x = "Pair", y = y_lab, pattern = "Number of Excesses") +
    facet_grid(rows = vars(Conditioning_Varaible), labeller = labeller(Conditioning_Varaible = facet_labels)) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.5)
  
  return(plot_out)
}
