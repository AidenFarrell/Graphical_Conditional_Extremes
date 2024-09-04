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
  facet_labels_rows <- setNames(paste0("i = ", 1:d_data), 1:d_data)
  facet_labels_cols <- setNames(c("n = 250", "n = 500"), levels(plot_data$Type))
  
  ## Make the plot
  plot_out <- ggplot(data = plot_data, aes(x = Dependent_Varaible, y = y, fill = Method)) + 
    geom_boxplot(outlier.shape = NA) + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top",
          legend.box = "horizontal",
          legend.box.just = "center",
          legend.spacing.x = unit(0.5, 'cm'),
          legend.spacing.y = unit(0.5, 'cm'),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 16),
          strip.text = element_text(size = 16),
          panel.grid.minor.y = element_blank()) +
    lims(y = ylims) +
    labs(x = "Dependent Variable (j)", y = y_lab) +
    facet_grid(rows = vars(Conditioning_Varaible), cols = vars(Type), 
               labeller = labeller(Conditioning_Varaible = facet_labels_rows, Type = facet_labels_cols)) +
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
  facet_labels_rows <- setNames(paste0("i = ", 1:d_data), 1:d_data)
  facet_labels_cols <- setNames(c("n = 250", "n = 500"), levels(plot_data$Type))
  
  ## Make the plot
  plot_out <- ggplot(data = plot_data, aes(x = Pair, y = y, fill = Method)) + 
    geom_boxplot(outlier.shape = NA) + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top",
          legend.box = "horizontal",
          legend.box.just = "center",
          legend.spacing.x = unit(0.5, 'cm'),
          legend.spacing.y = unit(0.5, 'cm'),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 16),
          strip.text = element_text(size = 16),
          panel.grid.minor.y = element_blank()) +
    lims(y = ylims) +
    labs(x = "Pair", y = y_lab) +
    facet_grid(rows = vars(Conditioning_Varaible), cols = vars(Type), 
               labeller = labeller(Conditioning_Varaible = facet_labels_rows, Type = facet_labels_cols)) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.5)
  
  return(plot_out)
}
