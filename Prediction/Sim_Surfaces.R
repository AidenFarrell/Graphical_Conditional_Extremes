################################################################################
## Reading in required packages
required_pckgs <- c("Matrix", "LaplacesDemon")
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in required scripts
source("Miscellaneous_Functions/MVAGG_Functions.R")
source("Miscellaneous_Functions/DL_funcs.R")

################################################################################
## Simulation of the spatial process
## See Richards et al. 2022 Algorithm 1
Sim_Surface_HT <- function(n_sim, q, transforms, CMEVM_fits){
  ## Checks
  if(!is.list(transforms) | !is.list(CMEVM_fits)){
    stop("transforms and CMEVM_fits must be lists")
  }
  if(!all(sapply(transforms, function(x){class(x) == "Laplace_Transform"}))){
    stop("transforms must be lists where each element in the list is of class 'Laplace_Transform'")
  }
  if(!all(sapply(CMEVM_fits, function(x){class(x) == "Cond_Extremes_MVN"}))){
    stop("transforms must be lists where each element in the list is of class 'Cond_Extremes_MVN'")
  }
  if(length(transforms) != length(CMEVM_fits)){
    stop("lengths of transforms and CMEVM_fits do not match")
  }
  
  ## Get some information on the data now
  d <- length(transforms)
  if(d == 1){
    stop("The dimension (length of transfoms) must be at least 2")
  }
  n <- length(transforms[[1]]$data$X)
  
  ## More checks
  if(!is.numeric(n_sim) | length(n_sim) != 1 | n_sim%%1 != 0 | n_sim <= n){
    stop("n_sim must be a single positive integer greater than n (the number of data points)")
  }
  if(!is.numeric(q) | length(q) != 1){
    stop("q must be a single positive real number in the range(0, 1)") 
  }
  else if(any(q <= 0) | any(q >= 1)){
    stop("q must be a single positive real number in the range(0, 1)")
  }
  
  ## Get the conditioning sites
  cond_sites <- sample(x = 1:d, size = n_sim, replace = TRUE)
  ## Convert the threshold onto Laplace margins and simulate from the conditioning site
  v <- qlaplace(q)
  large_Laplace <- v +  rexp(n = n_sim, rate = 1)
  ## Sample residuals from the fitted model
  Resid <- lapply(CMEVM_fits, function(x){x$Z})
  Z_sample <- t(sapply(cond_sites, function(i){
    unlist(unname(Resid[[i]][sample(1:nrow(Resid[[i]]), 1),]))}))
  ## Convert the simulated data onto Laplace margins
  a_hats <- lapply(CMEVM_fits, function(x){unname(x$par$main[1,])})
  b_hats <- lapply(CMEVM_fits, function(x){unname(x$par$main[2,])})
  a_cond_sites <- t(sapply(cond_sites, function(i){a_hats[[i]]}))
  b_cond_sites <- t(sapply(cond_sites, function(i){b_hats[[i]]}))
  Laplace_Samples <- matrix(NA, nrow = n_sim, ncol = d-1)
  for(i in 1:(d-1)){
    Laplace_Samples[,i] <- a_cond_sites[,i]*large_Laplace + (large_Laplace^(b_cond_sites[,i]))*Z_sample[,i]
  }
  ## Now combine Y{-i}|Y_{i} > u_{Y_{i}} and Y{-i}|Y_{i} > u_{Y_{i}}
  Laplace_Samples_Full <- matrix(NA, nrow = n_sim, ncol = d)
  for(i in 1:n_sim){
    r <- cond_sites[i]
    Laplace_Samples_Full[i,r] <- large_Laplace[i]
    if(r == 1){
      Laplace_Samples_Full[i,-1] <- Laplace_Samples[i,]
    }
    else if(r == d){
      Laplace_Samples_Full[i,-d] <- Laplace_Samples[i,]
    }
    else{
      Laplace_Samples_Full[i,(1:(r-1))] <- Laplace_Samples[i,(1:(r-1))]
      Laplace_Samples_Full[i,((r+1):d)] <- Laplace_Samples[i,(r:(d-1))]
    }
  }
  ## Now use importance sampling to resample the simulated data
  ## Idea is to up-weight those close to the boundary and down-weight those in the
  ## centre of the spatial process
  weights <- d/apply(Laplace_Samples_Full, 1, function(x){length(which(x > v))})
  Index <- sample(x = 1:n_sim, size = n, prob = weights/sum(weights))
  ## Now we convert the data onto the original scale
  Final_Laplace_Samples <- Laplace_Samples_Full[Index,]
  ## convert onto original scale
  Final_Uniform_Samples <- apply(Final_Laplace_Samples, 2, plaplace, simplify = FALSE)
  Final_Original_Samples <- mcmapply(FUN = texmex:::revTransform,
                                     x = Final_Uniform_Samples,
                                     data = lapply(transforms, function(x){x$data$X}),
                                     qu = sapply(transforms, function(x){x$par$qu}),
                                     th = sapply(transforms, function(x){x$par$u}),
                                     sigma = sapply(transforms, function(x){x$par$scale}),
                                     xi = sapply(transforms, function(x){x$par$shape}),
                                     MoreArgs = list(method = "mixture"),
                                     SIMPLIFY = TRUE,
                                     mc.cores = detectCores() - 1)
  ## Above simulated X(s) | max(X(s)) > v
  ## We want to get the unconditioned process
  ## First get the data points were we have no extremes
  data <- sapply(transforms, function(x){x$data$X})
  data_Laplace <- sapply(transforms, function(x){x$data$Y})
  Index_No_Extremes <- which(apply(data_Laplace, 1, max) < v)
  Data_Body <- data[Index_No_Extremes,]
  Data_Body_Laplace <- data_Laplace[Index_No_Extremes,]
  ## then get the probability we will draw from the joint tail
  Index_At_Least_One_Extreme <- which(apply(data_Laplace, 1, max) > v)
  p_tail <- length(Index_At_Least_One_Extreme)/(n + 1)
  p_accecpt <- runif(n)
  ## Get the index of body/tail
  Index_Body <- sample(x = 1:length(Index_No_Extremes), size = length(which(p_accecpt >= p_tail)), replace = TRUE)
  Index_Tail <- which(p_accecpt < p_tail)
  ## Get the final data sets
  Data_Final_Laplace_Margins <- rbind(Data_Body_Laplace[Index_Body,],
                                      Final_Laplace_Samples[Index_Tail,])
  Data_Final_Original_Margins <- rbind(Data_Body[Index_Body,],
                                       Final_Original_Samples[Index_Tail,])
  return(list(Data_Margins = Data_Final_Original_Margins,
              Laplace_Margins = Data_Final_Laplace_Margins))
}

## As above but we simulate from the fitted model for the residuals rather than
## sampling from the fitted residuals
## The fitted model is the asymmetric generalised Gaussian distribution
Sim_Surface_MVAGG <- function(n_sim, q, transforms, CMEVM_fits){
  ## Checks
  if(!is.list(transforms) | !is.list(CMEVM_fits)){
    stop("transforms and CMEVM_fits must be lists")
  }
  if(!all(sapply(transforms, function(x){class(x) == "Laplace_Transform"}))){
    stop("transforms must be lists where each element in the list is of class 'Laplace_Transform'")
  }
  if(!all(sapply(CMEVM_fits, function(x){class(x) == "Cond_Extremes_MVAGG"}))){
    stop("transforms must be lists where each element in the list is of class 'Cond_Extremes_MVAGG'")
  }
  if(length(transforms) != length(CMEVM_fits)){
    stop("lengths of transforms and CMEVM_fits do not match")
  }
  
  ## Get some infromation on the data now
  d <- length(transforms)
  if(d == 1){
    stop("The dimension (length of transfoms) must be at least 2")
  }
  n <- length(transforms[[1]]$data$X)
  
  ## More checks
  if(!is.numeric(n_sim) | length(n_sim) != 1 | n_sim%%1 != 0 | n_sim <= n){
    stop("n_sim must be a single positive integer greater than n (the number of data points)")
  }
  if(!is.numeric(q) | length(q) != 1){
    stop("q must be a single positive real number in the range(0, 1)") 
  }
  else if(any(q <= 0) | any(q >= 1)){
    stop("q must be a single positive real number in the range(0, 1)")
  }
  
  ## Get the conditioning sites
  cond_sites <- sample(x = 1:d, size = n_sim, replace = TRUE)
  cond_sites_unique <- sort(unique(cond_sites))
  n_cond_sites <- as.numeric(table(cond_sites))
  ## Convert the threshold onto Laplace margins and simulate from the conditioning site
  v <- qlaplace(q)
  large_Laplace <- v +  rexp(n = n_sim, rate = 1)
  ## Simulate residuals from the fitted model
  Z <- mcmapply(FUN = rmvagg,
                loc = lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[3,]}),
                scale_1 = lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[4,]}),
                scale_2 = lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[5,]}),
                shape = lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[6,]}),
                Gamma = lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$Gamma}),
                n = n_cond_sites,
                SIMPLIFY = FALSE,
                mc.cores = detectCores() - 1)
  ## Convert the simulated data onto Laplace margins
  a_cond_sites <- lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[1,]})
  b_cond_sites <- lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[2,]})
  
  Laplace_Samples <- vector("list", d)
  for(i in 1:d){
    large_Laplace_cond <- large_Laplace[cond_sites == i]
    Y_Yi_large <- sapply(a_cond_sites[[i]], function(a){a*large_Laplace_cond}) + 
      sapply(b_cond_sites[[i]], function(b){large_Laplace_cond^b})*Z[[i]]
    
    ## Now combine Y{-i}|Y_{i} > u_{Y_{i}} and Y{-i}|Y_{i} > u_{Y_{i}}
    if(i == 1){
      Laplace_Samples[[i]] <- cbind(large_Laplace_cond, Y_Yi_large)
    }
    else if(i == d){
      Laplace_Samples[[i]] <- cbind(Y_Yi_large, large_Laplace_cond)
    }
    else{
      Laplace_Samples[[i]] <- cbind(Y_Yi_large[,1:(i-1)], large_Laplace_cond, Y_Yi_large[,i:(d-1)])
    }
  }
  Laplace_Samples_Full <- do.call(rbind, Laplace_Samples)
  
  ## Now use importance sampling to resample the simulated data
  ## Idea is to up-weight those close to the boundary and down-weight those in the
  ## centre of the spatial process
  weights <- d/apply(Laplace_Samples_Full, 1, function(x){length(which(x > v))})
  Index <- sample(x = 1:n_sim, size = n, prob = weights/sum(weights))
  ## Now we convert the data onto the original scale
  Final_Laplace_Samples <- Laplace_Samples_Full[Index,]
  ## convert onto original scale
  Final_Uniform_Samples <- apply(Final_Laplace_Samples, 2, plaplace, simplify = FALSE)
  Final_Original_Samples <- mcmapply(FUN = texmex:::revTransform,
                                     x = Final_Uniform_Samples,
                                     data = lapply(transforms, function(x){x$data$X}),
                                     qu = sapply(transforms, function(x){x$par$qu}),
                                     th = sapply(transforms, function(x){x$par$u}),
                                     sigma = sapply(transforms, function(x){x$par$scale}),
                                     xi = sapply(transforms, function(x){x$par$shape}),
                                     MoreArgs = list(method = "mixture"),
                                     SIMPLIFY = TRUE,
                                     mc.cores = detectCores() - 1)
  ## Above simulated X(s) | max(X(s)) > v
  ## We want to get the unconditioned process
  ## First get the data points were we have no extremes
  data <- sapply(transforms, function(x){x$data$X})
  data_Laplace <- sapply(transforms, function(x){x$data$Y})
  Index_No_Extremes <- which(apply(data_Laplace, 1, max) < v)
  Data_Body <- data[Index_No_Extremes,]
  Data_Body_Laplace <- data_Laplace[Index_No_Extremes,]
  ## then get the probability we will draw from the joint tail
  Index_At_Least_One_Extreme <- which(apply(data_Laplace, 1, max) > v)
  p_tail <- length(Index_At_Least_One_Extreme)/(n + 1)
  p_accecpt <- runif(n = n)
  ## Get the index of body/tail
  Index_Body <- sample(x = 1:length(Index_No_Extremes), size = length(which(p_accecpt >= p_tail)), replace = TRUE)
  Index_Tail <- which(p_accecpt < p_tail)
  ## Get the final data sets
  Data_Final_Laplace_Margins <- rbind(Data_Body_Laplace[Index_Body,],
                                      Final_Laplace_Samples[Index_Tail,])
  Data_Final_Original_Margins <- rbind(Data_Body[Index_Body,],
                                       Final_Original_Samples[Index_Tail,])
  return(list(Data_Margins = Data_Final_Original_Margins,
              Laplace_Margins = Data_Final_Laplace_Margins))
}

## As above but we simulate from the fitted model for the residuals rather than
## sampling from the fitted residuals
## The fitted model is the generalised Gaussian (delta-Laplace) distribution
Sim_Surface_MVGG <- function(n_sim, q, transforms, CMEVM_fits){
  ## Checks
  if(!is.list(transforms) | !is.list(CMEVM_fits)){
    stop("transforms and CMEVM_fits must be lists")
  }
  if(!all(sapply(transforms, function(x){class(x) == "Laplace_Transform"}))){
    stop("transforms must be lists where each element in the list is of class 'Laplace_Transform'")
  }
  if(!all(sapply(CMEVM_fits, function(x){class(x) == "Cond_Extremes_MVAGG"}))){
    stop("transforms must be lists where each element in the list is of class 'Cond_Extremes_MVAGG'")
  }
  if(length(transforms) != length(CMEVM_fits)){
    stop("lengths of transforms and CMEVM_fits do not match")
  }
  
  ## Get some infromation on the data now
  d <- length(transforms)
  if(d == 1){
    stop("The dimension (length of transfoms) must be at least 2")
  }
  n <- length(transforms[[1]]$data$X)
  
  ## More checks
  if(!is.numeric(n_sim) | length(n_sim) != 1 | n_sim%%1 != 0 | n_sim <= n){
    stop("n_sim must be a single positive integer greater than n (the number of data points)")
  }
  if(!is.numeric(q) | length(q) != 1){
    stop("q must be a single positive real number in the range(0, 1)") 
  }
  else if(any(q <= 0) | any(q >= 1)){
    stop("q must be a single positive real number in the range(0, 1)")
  }
  
  ## Get the conditioning sites
  cond_sites <- sample(x = 1:d, size = n_sim, replace = TRUE)
  cond_sites_unique <- sort(unique(cond_sites))
  n_cond_sites <- as.numeric(table(cond_sites))
  ## Convert the threshold onto Laplace margins and simulate from the conditioning site
  v <- qlaplace(q)
  large_Laplace <- v +  rexp(n = n_sim, rate = 1)
  ## Simulate residuals from the fitted model
  ## Need to round the matrix otherwise it claims it is not symmetric when it is
  Z <- mcmapply(FUN = rmvdlaplace,
                mu = lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[3,]}),
                sigmad = lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[4,]}),
                delta = lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[5,]}),
                Sigma = lapply(cond_sites_unique, function(i){round(solve(as.matrix(CMEVM_fits[[i]]$par$Gamma)), 8)}),
                n = n_cond_sites,
                MoreArgs = list(dim = d - 1),
                SIMPLIFY = FALSE,
                mc.cores = detectCores() - 1)
  ## Convert the simulated data onto Laplace margins
  a_cond_sites <- lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[1,]})
  b_cond_sites <- lapply(cond_sites_unique, function(i){CMEVM_fits[[i]]$par$main[2,]})
  
  Laplace_Samples <- vector("list", d)
  for(i in 1:d){
    large_Laplace_cond <- large_Laplace[cond_sites == i]
    Y_Yi_large <- sapply(a_cond_sites[[i]], function(a){a*large_Laplace_cond}) + 
      sapply(b_cond_sites[[i]], function(b){large_Laplace_cond^b})*Z[[i]]
    
    ## Now combine Y{-i}|Y_{i} > u_{Y_{i}} and Y{-i}|Y_{i} > u_{Y_{i}}
    if(i == 1){
      Laplace_Samples[[i]] <- cbind(large_Laplace_cond, Y_Yi_large)
    }
    else if(i == d){
      Laplace_Samples[[i]] <- cbind(Y_Yi_large, large_Laplace_cond)
    }
    else{
      Laplace_Samples[[i]] <- cbind(Y_Yi_large[,1:(i-1)], large_Laplace_cond, Y_Yi_large[,i:(d-1)])
    }
  }
  Laplace_Samples_Full <- do.call(rbind, Laplace_Samples)
  
  ## Now use importance sampling to resample the simulated data
  ## Idea is to up-weight those close to the boundary and down-weight those in the
  ## centre of the spatial process
  weights <- d/apply(Laplace_Samples_Full, 1, function(x){length(which(x > v))})
  Index <- sample(x = 1:n_sim, size = n, prob = weights/sum(weights))
  ## Now we convert the data onto the original scale
  Final_Laplace_Samples <- Laplace_Samples_Full[Index,]
  ## convert onto original scale
  Final_Uniform_Samples <- apply(Final_Laplace_Samples, 2, plaplace, simplify = FALSE)
  Final_Original_Samples <- mcmapply(FUN = texmex:::revTransform,
                                     x = Final_Uniform_Samples,
                                     data = lapply(transforms, function(x){x$data$X}),
                                     qu = sapply(transforms, function(x){x$par$qu}),
                                     th = sapply(transforms, function(x){x$par$u}),
                                     sigma = sapply(transforms, function(x){x$par$scale}),
                                     xi = sapply(transforms, function(x){x$par$shape}),
                                     MoreArgs = list(method = "mixture"),
                                     SIMPLIFY = TRUE,
                                     mc.cores = detectCores() - 1)
  ## Above simulated X(s) | max(X(s)) > v
  ## We want to get the unconditioned process
  ## First get the data points were we have no extremes
  data <- sapply(transforms, function(x){x$data$X})
  data_Laplace <- sapply(transforms, function(x){x$data$Y})
  Index_No_Extremes <- which(apply(data_Laplace, 1, max) < v)
  Data_Body <- data[Index_No_Extremes,]
  Data_Body_Laplace <- data_Laplace[Index_No_Extremes,]
  ## then get the probability we will draw from the joint tail
  Index_At_Least_One_Extreme <- which(apply(data_Laplace, 1, max) > v)
  p_tail <- length(Index_At_Least_One_Extreme)/(n + 1)
  p_accecpt <- runif(n = n)
  ## Get the index of body/tail
  Index_Body <- sample(x = 1:length(Index_No_Extremes), size = length(which(p_accecpt >= p_tail)), replace = TRUE)
  Index_Tail <- which(p_accecpt < p_tail)
  ## Get the final data sets
  Data_Final_Laplace_Margins <- rbind(Data_Body_Laplace[Index_Body,],
                                      Final_Laplace_Samples[Index_Tail,])
  Data_Final_Original_Margins <- rbind(Data_Body[Index_Body,],
                                       Final_Original_Samples[Index_Tail,])
  return(list(Data_Margins = Data_Final_Original_Margins,
              Laplace_Margins = Data_Final_Laplace_Margins))
}
