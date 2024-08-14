## Probability Calculations

## P[X_{A} > u_{A} | X_{i} > u_{i}] some set A
p_data_surv_multi <- function(data, cond, uncon, u){
  
  ## Checks
  if(!is.matrix(data)){
    stop("data must be a matrix")
  }
  
  d <- dim(data)[2]
  if(length(cond) != 1 | cond <= 0 | cond > d | cond%%1 != 0){
    stop("cond must be a single positive integer no greater than the dimension of the data")
  }
  if(!is.numeric(uncon) | length(uncon) >= d | any(uncon <= 0) | any(uncon > d) | any(uncon%%1 != 0)){
    stop("uncon must be a vector of positive integer no greater than the dimension of the data")
  }
  if(!is.numeric(uncon) | length(u) != d){
    stop("u must be a vector equal to the dimension of the data ")
  }
  
  ## Obtain the number of times Xi > ui
  data_cond_large <- data[which(data[,cond] > u[cond]),]
  nu_x <- nrow(data_cond_large)
  ## Now get the conditional probability
  p <- length(which(apply(t(data_cond_large[,uncon]) > u[uncon], 2, all)))/nu_x
  return(p)
}

## P[X_{A} <= u_{A} | X_{i} > u_{i}] some set A
p_data_cond_cdf <- function(data, cond, uncon, u){
  
  ## Checks
  if(!is.matrix(data)){
    stop("data must be a matrix")
  }
  
  d <- dim(data)[2]
  if(length(cond) != 1 | cond <= 0 | cond > d | cond%%1 != 0){
    stop("cond must be a single positive integer no greater than the dimension of the data")
  }
  if(!is.numeric(uncon) | length(uncon) >= d | any(uncon <= 0) | any(uncon > d) | any(uncon%%1 != 0)){
    stop("uncon must be a vector of positive integer no greater than the dimension of the data")
  }
  if(!is.numeric(uncon) | length(u) != d){
    stop("u must be a vector equal to the dimension of the data ")
  }
  
  ## Obtain the number of times Xi > ui
  data_cond_large <- data[which(data[,cond] > u[cond]),]
  nu_x <- nrow(data_cond_large)
  ## Now get the conditional probability
  p <- length(which(apply(t(data_cond_large[,uncon]) <= u[uncon], 2, all)))/nu_x
  return(p)
}