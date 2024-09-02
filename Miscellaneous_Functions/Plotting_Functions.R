## boxplots for parameter estiamtes
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

## Boxplot of bias in the precision matrix
boxplot_MLEs_Cov_Mat_Bias <- function(data, methods, y_lab, cov_mat_true, precision = FALSE){
  
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
  
  ## Get the true parameters in the correct form
  lower_tri_elements <- lower.tri(cov_mat_true, diag = TRUE)
  cov_mat_true_i <- lapply(1:d_data, function(i){Cond_Sigma(cov_mat_true, i)})
  cor_mat_true_i <- lapply(cov_mat_true_i, cov2cor)
  if(precision){
    cor_mat_true_i <- lapply(cor_mat_true_i, solve)
  }
  cor_mat_true_i_NA <- lapply(1:d_data, function(i){Add_NA_Matrix(cor_mat_true_i[[i]], i)[lower_tri_elements]})
  
  ## get the bias
  data_bias <- lapply(1:d_data, function(i){lapply(1:n_methods, function(j){
    data[[i]][[j]] - matrix(rep(cor_mat_true_i_NA[[i]], n), nrow = n, byrow = TRUE)
  })})
  
  ## Construct the plotting data
  plot_data <- data.frame(y = do.call(c, lapply(1:d_data, function(i){
    do.call(c, lapply(1:p, function(j){
      do.call(c, lapply(1:n_methods, function(k){
        data_bias[[i]][[k]][,j]}))}))})))
  
  plot_data$Method = rep(rep(rep(methods, each = n), p), d_data)
  plot_data$Method <- factor(plot_data$Method, levels = methods)
  plot_data$Conditioning_Varaible <- rep(1:d_data, each = p*n*n_methods)
  plot_data$Conditioning_Varaible <- factor(plot_data$Conditioning_Varaible, levels = 1:d_data)
  plot_data$Pair <- rep(rep(x_labels, each = n*n_methods), d_data)
  plot_data$Pair <- factor(plot_data$Pair, levels = x_labels)
  
  ## Define custom labels for facets
  facet_labels <- setNames(paste0("i = ", 1:d_data), 1:d_data)
  
  plot_out <- ggplot(data = plot_data, aes(x = Pair, y = y, fill = Method)) + 
    geom_boxplot() +
    theme(legend.position = "top") +
    guides(colour = guide_legend(nrow = 1)) +
    labs(x = "Pair", y = y_lab) +
    facet_grid(rows = vars(Conditioning_Varaible), labeller = labeller(Conditioning_Varaible = facet_labels)) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.5)
  
  return(plot_out)
}

## Boxplot of the precision matrix
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
    guides(colour = guide_legend(nrow = 1)) +
    labs(x = "Pair", y = y_lab) +
    facet_grid(rows = vars(Conditioning_Varaible), labeller = labeller(Conditioning_Varaible = facet_labels),
               scales = "free_y") +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed", linewidth = 0.25)
  
  return(plot_out)
}