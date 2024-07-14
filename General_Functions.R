## Function to get the covaraince matrix when we condition on one random variable
Cond_Sigma <- function(x, j){
  if(!is.square.matrix(x)){
    stop("x must be a square matrix")
  }
  else if(any(x - t(x) > 1e-6)){
    stop("x is not a symmetric matrix")
  }
  else if(any(eigen(x)$values <= 0) | is.complex(eigen(x)$values)){
    stop("The matrix has negative/Complex eigenvalues. This is not a valid variance-covariance matrix")
  }
  
  d <- dim(x)[1] + 1
  if(j < 0 | j > d | j%%1 != 0){
    stop("j must be a positive integer between 1 and d")
  }
  
  Sigma_11 <- x[-j,-j]
  Sigma_12 <- x[j,-j]
  Sigma_22 <- x[j,j]
  Sigma_bar <- Sigma_11 - Sigma_12%*%solve(Sigma_22)%*%t(Sigma_12)
  return(Sigma_bar)
}

## Function to convert a covariance matrix to a correlation matrix
Sigma2rho <- function(Sigma){
  d <- dim(Sigma)[1]
  rho <- diag(1, d)
  for(i in 1:d){
    for(j in 1:d){
      if(i >= j){
        next()
      }
      else{
        rho[i,j] <- rho[j,i] <- Sigma[i,j]/sqrt(Sigma[i,i]*Sigma[j,j])
      }
    }
  }
  return(rho)
}

## add NA into a vector position
Add_NA_Vector <- function(x, j){
  d <- length(x) + 1
  if(j < 0 | j > d | j%%1 != 0){
    stop("j must be a positive integer between 1 and d")
  }
  if(j == 1){
    x <- c(NA, x)
  }
  else if(j == d){
    x <- c(x, NA)
  }
  else{
    x <- c(x[1:(j-1)], NA, x[j:(d-1)])
  }
  return(x)
}

## add row and column of NAs into a matrix
Add_NA_Matrix <- function(x, j){
  if(!is.square.matrix(x)){
    stop("x must be a square matrix")
  }
  
  d <- dim(x)[1] + 1
  if(j < 0 | j > d | j%%1 != 0){
    stop("j must be a positive integer between 1 and d")
  }
  
  y <- matrix(NA, nrow = d, ncol = d)
  if(j == 1 | j == d){
    y[-j,-j] <- x
  }
  else{
    y[1:(j-1), 1:(j-1)] <- x[1:(j-1), 1:(j-1)]
    y[(j+1):d, 1:(j-1)] <- x[j:(d-1), 1:(j-1)]
    y[1:(j-1), (j+1):d] <- x[1:(j-1), j:(d-1)]
    y[(j+1):d, (j+1):d] <- x[j:(d-1), j:(d-1)]  
  }
  return(y)
}