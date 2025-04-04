################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("ggplot2", "graphicalExtremes", "igraph", "mev", "parallel")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Reading in general functions
source("Miscellaneous_Functions/General_Functions.R")

## Reading in functions to transform the data
source("Miscellaneous_Functions/Transformations.R")

## Reading in functions required for model fitting  - Three-step AGG
source("Model_Fitting/Cond_Extremes_MVAGG_Residuals_Three_Step.R")

## Reading in functions required for model fitting - Three-step GG
source("Model_Fitting/Cond_Extremes_MVGG_Residuals_Three_Step.R")

## Reading in functions required for simulating data from the model
source("Prediction/Sim_Surfaces.R")

################################################################################
## Plot the upper Danube River basin
danube_graph <- graph(t(danube$flow_edges), directed = FALSE)
layout <- cbind(danube$info$PlotCoordX, danube$info$PlotCoordY)
V(danube_graph)$size <- 15
V(danube_graph)$label.cex <- 2
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot(danube_graph, layout = layout)

## Get the flow connections for later
d <- length(V(danube_graph))
flow_connections <- matrix(0, nrow = d, ncol = d)
flow_connections[1,-1] <- 1
flow_connections[2,-c(1:2,13,28:31)] <- 1
flow_connections[3,-c(1:3,13,14:19,28:31)] <- 1
flow_connections[4,-c(1:4,13,14:19,28:31)] <- 1
flow_connections[5,c(6:12,20:22)] <- 1
flow_connections[6,c(7:12,20:22)] <- 1
flow_connections[7,c(8:12,20:22)] <- 1
flow_connections[8,c(9,10,11,12)] <- 1
flow_connections[9,c(10,11,12)] <- 1
flow_connections[10,c(11,12)] <- 1
flow_connections[11,12] <- 1
flow_connections[13,28:31] <- 1
flow_connections[14,15:19] <- 1
flow_connections[15,16:19] <- 1
flow_connections[16,17:19] <- 1
flow_connections[17,18:19] <- 1
flow_connections[18,19] <- 1
flow_connections[20,21:22] <- 1
flow_connections[21,22] <- 1
flow_connections[23,24] <- 1
flow_connections[25,26:27] <- 1
flow_connections[26,27] <- 1
flow_connections[28,29:31] <- 1
flow_connections[29,30:31] <- 1
flow_connections[30,31] <- 1

lower_tri_flow <- lower.tri(flow_connections)
flow_connections[lower_tri_flow] <- t(flow_connections)[lower_tri_flow]

## Euclidean distances between the sites
danube_distances <- as.matrix(dist(x = cbind(danube$info$Long, danube$info$Lat),
                                   diag = TRUE, upper = TRUE)*60*1.852)

## Sites that are connected (lowest site first)
flow_connected <- which(flow_connections == 1, arr.ind = TRUE)
flow_connected <- flow_connected[which(flow_connected[,1] < flow_connected[,2]),]
flow_connected <- flow_connected[order(flow_connected[,1]),]
flow_connected <- as.data.frame(flow_connected)

################################################################################
## DO NOT RUN

# ## Bootstrap the data
# n_boot <- 200
# n_data <- nrow(danube$data_clustered)
# set.seed(seed)
# danube_boot <- replicate(n = n_boot, 
#                          expr = danube$data_clustered[sample(x = 1:n_data, n_data, replace = TRUE),],
#                          simplify = FALSE)
# 
# ## Transform data onto standard Laplace margins using the Coles and Tawn Method
# danube_boot_list <- lapply(danube_boot, function(x){
#   lapply(apply(x, 2, list), function(x){x[[1]]})})
# 
# X_to_Y <- lapply(1:n_boot, function(i){
#   mcmapply(FUN = X_to_Laplace,
#            x = danube_boot_list[[i]],
#            MoreArgs = list(q = seq(0.55, 0.925, by = 0.01)),
#            SIMPLIFY = FALSE,
#            mc.cores = detectCores() - 1)})
# 
# saveRDS(X_to_Y, file = "Data/Danube_Bootstrapped.RData")
################################################################################
# ## Read in data
X_to_Y <- readRDS("Data/Danube_Bootstrapped.RData")

## Get the output
X <- lapply(X_to_Y, function(x){sapply(x, function(y){y$data$X})})
Y <- lapply(X_to_Y, function(x){sapply(x, function(y){y$data$Y})})

n_boot <- length(X)
n_data <- nrow(X[[1]])

################################################################################
## Now we want to subset the data so that each component is large in turn
dqu <- 0.8

Y_u <- qlaplace(dqu)
Y_Yi_large <- rep(list(vector("list", d)), n_boot)
for(i in 1:n_boot){
  for(j in 1:d){
    Y_Yi_large[[i]][[j]] <- Y[[i]][which(Y[[i]][,j] > Y_u),]
  } 
}

################################################################################
## Fit the Three-Step Saturated models using the MVAGG and MVGG for the residual distribution
v <- ceiling(max(sapply(Y, max))) + 1

## MVAGG model
fits_MVAGG <- lapply(1:n_boot, function(i){
  mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
           data = Y_Yi_large[[i]],
           cond = 1:d,
           MoreArgs = list(graph = make_full_graph(n = d),
                           v = v,
                           start_HT = c(0.5, 0.1)),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

## MVGG model
fits_MVGG <- lapply(1:n_boot, function(i){
  mcmapply(FUN = Cond_Extremes_MVGG_Three_Step,
           data = Y_Yi_large[[i]],
           cond = 1:d,
           MoreArgs = list(graph = make_full_graph(n = d),
                           v = v,
                           start_HT = c(0.5, 0.1)),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

################################################################################
## Simulate surfaces from the model

## MVAGG
X_MVAGG <- mcmapply(Sim_Surface_MVAGG,
                    transforms = X_to_Y,
                    CMEVM_fits = fits_MVAGG,
                    MoreArgs = list(n_sim = 20*n_data, q = dqu),
                    SIMPLIFY = FALSE,
                    mc.cores = detectCores() - 1)

## MVGG
X_MVGG <- mcmapply(Sim_Surface_MVGG,
                   transforms = X_to_Y,
                   CMEVM_fits = fits_MVGG,
                   MoreArgs = list(n_sim = 20*n_data, q = dqu),
                   SIMPLIFY = FALSE,
                   mc.cores = detectCores() - 1)

################################################################################
## Obtain summary statistics from the data
## Look at eta, chi and chibar

u <- c(0.8, 0.85, 0.9)
n_u <- length(u)

grid <- combinations(d, r = 2, v = 1:d)
grid_list <- lapply(apply(grid, 1, list), function(x){x[[1]]})
n_grid <- length(grid_list)

eta_boot <- eta_MVAGG_boot <- eta_MVGG_boot <- array(NA, dim = c(n_grid, n_u, n_boot))
eta_boot_med <- eta_MVAGG_med <- eta_MVGG_med <- array(NA, dim = c(n_grid, n_u))
eta_boot_se <- eta_MVAGG_se <- eta_MVGG_se <- array(NA, dim = c(n_grid, n_u))

eta_data <- lapply(X, function(x){
  mcmapply(FUN = taildep,
           data = lapply(grid_list, function(j){x[,j]}),
           MoreArgs = list(u = u, 
                           plot = FALSE),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

eta_MVAGG <- lapply(X_MVAGG, function(x){
  mcmapply(FUN = taildep,
           data = lapply(grid_list, function(j){x$Data_Margins[,j]}),
           MoreArgs = list(u = u, 
                           plot = FALSE),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

eta_MVGG <- lapply(X_MVGG, function(x){
  mcmapply(FUN = taildep,
           data = lapply(grid_list, function(j){x$Data_Margins[,j]}),
           MoreArgs = list(u = u, 
                           plot = FALSE),
           SIMPLIFY = FALSE,
           mc.cores = detectCores() - 2)
})

## Extract the eta values
for(k in 1:n_boot){
  eta_boot[,,k] <- t(sapply(eta_data[[k]], function(x){x$eta[,1]}))
  eta_MVAGG_boot[,,k] <- t(sapply(eta_MVAGG[[k]], function(x){x$eta[,1]}))
  eta_MVGG_boot[,,k] <- t(sapply(eta_MVGG[[k]], function(x){x$eta[,1]}))
}

## Obtain the mean and sd over the samples
for(i in 1:n_grid){
  for(j in 1:n_u){
    
    ## Data
    eta_boot_med[i,j] <- mean(eta_boot[i,j,])
    eta_boot_se[i,j] <- sd(eta_boot[i,j,])
    
    ## MVAGG
    eta_MVAGG_med[i,j] <- mean(eta_MVAGG_boot[i,j,])
    eta_MVAGG_se[i,j] <- sd(eta_MVAGG_boot[i,j,])
    
    ## MVGG
    eta_MVGG_med[i,j] <- mean(eta_MVGG_boot[i,j,])
    eta_MVGG_se[i,j] <- sd(eta_MVGG_boot[i,j,])
  }
}

## plot the output
methods <- c("MVAGG", "MVGG")
n_methods <- length(methods)

## eta plot
eta_out_comp <- data.frame(Site_1 = rep(grid[,1], n_methods*n_u),
                           Site_2 = rep(grid[,2], n_methods*n_u))
eta_out_comp$Connected <- do.call(paste0, eta_out_comp[,1:2]) %in% do.call(paste0, flow_connected)
eta_out_comp$u <- rep(rep(u, each = n_grid), n_methods)
eta_out_comp$Method <- rep(methods, each = n_grid*n_u)
eta_out_comp$x <- rep(c(eta_boot_med), n_methods)
eta_out_comp$y <- c(c(eta_MVAGG_med), c(eta_MVGG_med))
eta_out_comp$se <- c(c(eta_MVAGG_se), c(eta_MVGG_se))

eta_out_comp$u <- factor(eta_out_comp$u, levels = u)
eta_out_comp$Method <- factor(eta_out_comp$Method, levels = methods)

label_y <- function(labels) {
  sapply(labels, function(label) {
    substitute(eta(label), list(label = label))
  })
}

# pdf("Images/Danube/Bootstrapped_Ouput/ETA_Comp.pdf", height = 15, width = 15)
ggplot(data = eta_out_comp) + geom_point(aes(x = x, y = y, col = se, shape = Connected, alpha = Connected)) +
  lims(x = c(0.6, 1), y = c(0.6, 1)) +
  labs(x = "Empirical", y = "Theoretical", shape = "Flow-connected", color = "Standard error") +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed", linewidth = 0.5) +
  scale_shape_manual(values = c(1, 2)) +
  scale_colour_gradient(low = "blue", high = "red",
                        breaks = c(0.02, 0.03, 0.04, 0.05)) +
  scale_alpha_manual(values = c(0.75, 0.3)) +  # Set transparency: 0.3 for FALSE, 1 for TRUE
  theme(legend.position = "top",
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  guides(col = guide_colorbar(title.position = "left", title.vjust = 0.9,
                              barwidth = 15),  
         shape = guide_legend(title.position = "left"),
         alpha = "none") +
  facet_grid(cols = vars(Method), rows = vars(u),
             labeller = labeller(u = as_labeller(label_y, default = label_parsed)))
# dev.off()

################################################################################
## Transform the data onto standard Laplace margins
d <- ncol(danube$data_clustered)
n_data <- nrow(danube$data_clustered)

danube_data_list <- lapply(apply(danube$data_clustered, 2, list), function(x){x[[1]]})
X_to_Y <- mcmapply(FUN = X_to_Laplace,
                   x = danube_data_list,
                   MoreArgs = list(q = seq(0.55, 0.925, by = 0.01)),
                   SIMPLIFY = FALSE,
                   mc.cores = detectCores() - 1)

## Get the output
u_final <- unname(sapply(X_to_Y, function(x){unname(x$par$u)}))
qu_final <- unname(sapply(X_to_Y, function(x){unname(x$par$qu)}))
scale_final <- unname(sapply(X_to_Y, function(x){unname(x$par$scale)}))
shape_final <- unname(sapply(X_to_Y, function(x){unname(x$par$shape)}))
Y <- sapply(X_to_Y, function(x){x$data$Y})

## Now we want to subset the data so that each component is large in turn
dqu <- 0.8
Y_u <- qlaplace(dqu)

Y_Yi_large <- rep(list(list()), d)
for(i in 1:d){
  Y_Yi_large[[i]] <- Y[which(Y[,i] > Y_u),]
}

################################################################################
## Fit the H+T model
v <- ceiling(max(Y)) + 1
fit_HT <- mcmapply(FUN = Cond_Extremes_MVN,
                   data = Y_Yi_large, 
                   cond = 1:d,
                   MoreArgs = list(graph = NA,
                                   v = v,
                                   start = c(0.5, 0.1)),
                   SIMPLIFY = FALSE,
                   mc.cores = detectCores() - 1)

## Extract the residuals
Z_hat <- lapply(fit_HT, function(x){x$Z})

fits_MVAGG <- mcmapply(FUN = Cond_Extremes_MVAGG_Three_Step,
                       data = Y_Yi_large,
                       cond = 1:d,
                       MoreArgs = list(graph = make_full_graph(n = d),
                                       v = v,
                                       start_HT = c(0.5, 0.1)),
                       SIMPLIFY = FALSE,
                       mc.cores = detectCores() - 2)

fits_MVGG <- mcmapply(FUN = Cond_Extremes_MVGG_Three_Step,
                      data = Y_Yi_large,
                      cond = 1:d,
                      MoreArgs = list(graph = make_full_graph(n = d),
                                      v = v,
                                      start_HT = c(0.5, 0.1)),
                      SIMPLIFY = FALSE,
                      mc.cores = detectCores() - 2)

scale_left_AGG <- sapply(fits_MVAGG, function(x){unname(x$par$main[4,])})
scale_right_AGG <- sapply(fits_MVAGG, function(x){unname(x$par$main[5,])})

scale_GG <- sapply(fits_MVGG, function(x){unname(x$par$main[4,])})

for(i in 1:d){
  plot(scale_left_AGG[,i],
       scale_right_AGG[,i],
       main = paste("Conditioning Station - ", i))
  points(scale_GG[,i],
         scale_GG[,i], 
         col = "blue")
  abline(a = 0, b = 1, col = 2, lty = 2)
}

################################################################################
## simualte surfaces from the fitted models
sim_mvagg <- Sim_Surface_MVAGG(transforms = X_to_Y,
                               CMEVM_fits = fits_MVAGG,
                               n_sim = 20*n_data, 
                               q = dqu)

sim_mvgg <- Sim_Surface_MVGG(transforms = X_to_Y,
                             CMEVM_fits = fits_MVGG,
                             n_sim = 20*n_data, 
                             q = dqu)

## Compare the eta values
grid <- expand.grid(1:d, 1:d)
grid_list <- lapply(apply(grid, 1, list), function(x){x[[1]]})

u = c(0.8, 0.85, 0.9)

eta_data <- mcmapply(FUN = taildep,
                     data = lapply(grid_list, function(i){danube$data_clustered[,i]}),
                     MoreArgs = list(u = u,
                                     depmeas = "eta"),
                     SIMPLIFY = FALSE,
                     mc.cores = detectCores() - 2)

eta_mvagg <- mcmapply(FUN = taildep,
                      data = lapply(grid_list, function(i){sim_mvagg$Data_Margins[,i]}),
                      MoreArgs = list(u = u,
                                      depmeas = "eta"),
                      SIMPLIFY = FALSE,
                      mc.cores = detectCores() - 2)

eta_mvgg <- mcmapply(FUN = taildep,
                      data = lapply(grid_list, function(i){sim_mvgg$Data_Margins[,i]}),
                      MoreArgs = list(u = u,
                                      depmeas = "eta"),
                      SIMPLIFY = FALSE,
                      mc.cores = detectCores() - 2)

n_grid <- nrow(grid)
n_u <- length(u)

methods <- c("MVAGG", "MVGG")
n_methods <- length(methods)
df <- data.frame(Station_1 = rep(rep(grid[,1], each = n_u), n_methods),
                 Station_2 = rep(rep(grid[,2], each = n_u), n_methods),
                 u = rep(u, times = n_grid*n_methods),
                 Method = rep(methods, each = n_u*n_grid),
                 Empirical = rep(do.call(c, lapply(eta_data, function(x){x$eta[,1]})), times = n_methods),
                 Theoretical = c(do.call(c, lapply(eta_mvagg, function(x){x$eta[,1]})),
                                 do.call(c, lapply(eta_mvgg, function(x){x$eta[,1]}))))

df$Station_1 <- factor(df$Station_1, levels = 1:d, labels = 1:d)
df$Station_2 <- factor(df$Station_2, levels = 1:d, labels = 1:d)
df$u <- factor(df$u, levels = u, labels = u)
df$Method <- factor(df$Method, levels = methods, labels = methods)

ggplot(data = df) +
  geom_point(aes(x = Empirical, y = Theoretical)) +
  geom_abline(intercept = 0, slope = 1, col = 2, linetype = "dashed", linewidth = 1) + 
  facet_grid(rows = vars(u),
             cols = vars(Method))

sapply(eta_data, function(x){x$eta[,1]})
eta_data[[1]]$eta

i <- 19
j <- 29

plot(Y[,i], Y[,j])
points(sim_mvagg$Laplace_Margins[,i], sim_mvagg$Laplace_Margins[,j], col = "blue")
points(sim_mvgg$Laplace_Margins[,i], sim_mvgg$Laplace_Margins[,j], col = "red")


plot(danube$data_clustered[,i], danube$data_clustered[,j])
points(sim_mvagg$Data_Margins[,i], sim_mvagg$Data_Margins[,j], col = "blue")
points(sim_mvgg$Data_Margins[,i], sim_mvgg$Data_Margins[,j], col = "red")



sim_mvgg <- rmvdlaplace(n = 85, 
                        dim = d - 1,
                        mu = fits_MVGG[[1]]$par$main[3,],
                        sigmad = fits_MVGG[[1]]$par$main[4,],
                        delta = fits_MVGG[[1]]$par$main[5,],
                        Sigma = round(solve(as.matrix(fits_MVGG[[1]]$par$Gamma)), 10))



################################################################################
## Sites 19 and 16
dat_19_16 <- Z_hat[[19]][,16]

## DL and AGG fits
fit_ddlaplce_19_16 <- optim(par = c(0, 1, 1),
                            fn = fit_ddlaplce,
                            z = dat_19_16,
                            negative = TRUE)

fit_agg_19_16 <- fit_agg(par = c(0, 1, 1, 1),
                         data = dat_19_16)

## Empirical and Theoretical densities
dens_z_19_16 <- density(dat_19_16)
dens_z_19_16_dlaplace <- ddlaplace(z = dens_z_19_16$x,
                                   mu = fit_ddlaplce_19_16$par[1],
                                   sigma = fit_ddlaplce_19_16$par[2],
                                   delta = fit_ddlaplce_19_16$par[3])
dens_z_19_16_agg <- dagg(x = dens_z_19_16$x,
                         loc = fit_agg_19_16$par[1],
                         scale_1 = fit_agg_19_16$par[2],
                         scale_2 = fit_agg_19_16$par[3],
                         shape = fit_agg_19_16$par[4])

## Sites 19 and 29
dat_19_29 <- Z_hat[[19]][,28]

## DL and AGG fits
fit_ddlaplce_19_29 <- optim(par = c(0, 1, 1),
                            fn = fit_ddlaplce,
                            z = dat_19_29,
                            negative = TRUE)
fit_agg_19_29 <- fit_agg(par = c(0, 1, 1, 1),
                         data = dat_19_29)

## Empirical and Theoretical densities
dens_z_19_29 <- density(dat_19_29)
dens_z_19_29_dlaplace <- ddlaplace(z = dens_z_19_29$x,
                                   mu = fit_ddlaplce_19_29$par[1],
                                   sigma = fit_ddlaplce_19_29$par[2],
                                   delta = fit_ddlaplce_19_29$par[3])
dens_z_19_29_agg <- dagg(x = dens_z_19_29$x,
                         loc = fit_agg_19_29$par[1],
                         scale_1 = fit_agg_19_29$par[2],
                         scale_2 = fit_agg_19_29$par[3],
                         shape = fit_agg_19_29$par[4])


methods <- c("KDE", "Generlaised Gaussian", "AGG")
n_methods <- length(methods)

density_name <- c(substitute(z[16 ~ "|" ~ 19]), substitute(z[29 ~ "|" ~ 19]))
n_densities <- length(density_name)

n_data <- c(length(dens_z_19_16$x), length(dens_z_19_29$x))

Resid_dens <- data.frame(z = c(rep(dens_z_19_16$x, n_methods),
                               rep(dens_z_19_29$x, n_methods)),
                         y = c(c(dens_z_19_16$y, dens_z_19_16_dlaplace, dens_z_19_16_agg),
                               c(dens_z_19_29$y, dens_z_19_29_dlaplace, dens_z_19_29_agg)),
                         method = rep(c(sapply(n_data, function(x){rep(methods, each = x)})), times = n_densities),
                         density = c(sapply(1:n_densities, function(i){rep(i, times = n_data[i]*n_methods)})))
Resid_dens$method <- factor(Resid_dens$method, levels = methods, labels = methods)
Resid_dens$density <- factor(Resid_dens$density, levels = 1:n_densities, labels = density_name)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(Resid_dens) +
  geom_line(aes(x = z, y = y, col = method)) + 
  labs(x =  "z", y = "Density", col = "Method", 
       title = "Comparison of marginal residual distributions") +
  scale_color_manual(values = cbbPalette[1:n_methods],
                     breaks = methods) +
  theme(legend.position = "top") +
  facet_wrap(~density, ncol = 2,
             scales = "free",
             labeller = label_parsed)

plot(Z_hat[[19]][,28])
points(rdlaplace(n = n_data[1],
                 mu = fit_ddlaplce_19_29$par[1],
                 sigma = fit_ddlaplce_19_29$par[2],
                 delta = fit_ddlaplce_19_29$par[3]),
       col = cbbPalette[2])
points(ragg(n = n_data[1],
            loc = fit_agg_19_29$par[1],
            scale_1 = fit_agg_19_29$par[2],
            scale_2 = fit_agg_19_29$par[3],
            shape = fit_agg_19_29$par[4]),
       col = cbbPalette[3])

plot(x = dens_z_19_16$x,
     y = dens_z_19_16$y,
     type = "l",
     ylim = c(0, max(dens_z_19_16$y, dens_z_19_16_dlaplace, dens_z_19_16_agg)))
lines(x = dens_z_19_16$x,
      y = dens_z_19_16_dlaplace,
      col = "gold")
lines(x = dens_z_19_16$x,
      y = dens_z_19_16_agg,
      col = "blue")

sqrt(mean(((dens_z_19_16$y - dens_z_19_16_agg)^2)))
sqrt(mean(((dens_z_19_16$y - dens_z_19_16_dlaplace)^2)))

n_points <- length(dat_19_16)
p <- (1:n_points)/(n_points + 1)

q_emp <- sort(dat_19_16)
q_theo_19_16_dlaplce <- qdlaplace(p,
                                  mu = fit_ddlaplce_19_16$par[1],
                                  sigma = fit_ddlaplce_19_16$par[2],
                                  delta = fit_ddlaplce_19_16$par[3])
q_theo_19_16_agg <- qagg(p = p,
                         loc = fit_agg_19_16$par[1],
                         scale_1 = fit_agg_19_16$par[2],
                         scale_2 = fit_agg_19_16$par[3],
                         shape = fit_agg_19_16$par[4])

plot(q_emp,
     q_theo_19_16_dlaplce,
     xlab = "Empirical",
     ylab = "Theoretical",
     ylim = c(min(q_emp, q_theo_19_16_dlaplce, q_theo_19_16_agg),
              max(q_emp, q_theo_19_16_dlaplce, q_theo_19_16_agg)))
points(q_emp,
       q_theo_19_16_agg,
       col = "blue")
abline(a = 0, b = 1, col = 2)

mean(abs(q_emp - q_theo_19_16_agg))
mean(abs(q_emp - q_theo_19_16_dlaplce))

sqrt(mean(((q_emp - q_theo_19_16_agg)^2)))
sqrt(mean(((q_emp - q_theo_19_16_dlaplce)^2)))

## PP plot
p_theo_19_16_dlaplce <- pdlaplace(z = sort(dat_19_16),
                                  mu = fit_ddlaplce_19_16$par[1],
                                  sigma = fit_ddlaplce_19_16$par[2],
                                  delta = fit_ddlaplce_19_16$par[3])
p_theo_19_16_agg <- pagg(q = sort(dat_19_16),
                         loc = fit_agg_19_16$par[1],
                         scale_1 = fit_agg_19_16$par[2],
                         scale_2 = fit_agg_19_16$par[3],
                         shape = fit_agg_19_16$par[4])

plot(p,
     p_theo_19_16_dlaplce,
     xlab = "Empirical",
     ylab = "Theoretical",
     xlim = c(0, 1),
     ylim = c(0, 1))
points(p,
       p_theo_19_16_agg,
       col = "blue")
abline(a = 0, b = 1, col = 2)

mean(abs(p - p_theo_19_16_agg))
mean(abs(p - p_theo_19_16_dlaplce))

sqrt(mean(((p - p_theo_19_16_agg)^2)))
sqrt(mean(((p - p_theo_19_16_dlaplce)^2)))

################################################################################
## Sites 19 and 29
dat_19_29 <- Z_hat[[19]][,28]
fit_ddlaplce_19_29 <- optim(par = c(0, 1, 1),
                            fn = fit_ddlaplce,
                            z = dat_19_29,
                            negative = TRUE)

fit_agg_19_29 <- fit_agg(par = c(0, 1, 1, 1),
                         data = dat_19_29)

fit_HT[[19]]$par$main[,28]
plot(Y_Yi_large[[19]][,19],
     Y_Yi_large[[19]][,29])
plot(density(Y_Yi_large[[19]][,29]))

dens_z_19_29 <- density(dat_19_29)
dens_z_19_29_dlaplace <- ddlaplace(z = dens_z_19_29$x,
                                   mu = fit_ddlaplce_19_29$par[1],
                                   sigma = fit_ddlaplce_19_29$par[2],
                                   delta = fit_ddlaplce_19_29$par[3])
dens_z_19_29_agg <- dagg(x = dens_z_19_29$x,
                         loc = fit_agg_19_29$par[1],
                         scale_1 = fit_agg_19_29$par[2],
                         scale_2 = fit_agg_19_29$par[3],
                         shape = fit_agg_19_29$par[4])

plot(x = dens_z_19_29$x,
     y = dens_z_19_29$y,
     type = "l",
     ylim = c(0, max(dens_z_19_29$y, dens_z_19_29_dlaplace, dens_z_19_29_agg)))
lines(x = dens_z_19_29$x,
      y = dens_z_19_29_dlaplace,
      col = "gold")
lines(x = dens_z_19_29$x,
      y = dens_z_19_29_agg,
      col = "blue")

sqrt(mean(((dens_z_19_29$y - dens_z_19_29_agg)^2)))
sqrt(mean(((dens_z_19_29$y - dens_z_19_29_dlaplace)^2)))

n_points <- length(dat_19_29)
p <- (1:n_points)/(n_points + 1)

q_emp <- sort(dat_19_29)
q_theo_19_29_dlaplce <- qdlaplace(p,
                                  mu = fit_ddlaplce_19_29$par[1],
                                  sigma = fit_ddlaplce_19_29$par[2],
                                  delta = fit_ddlaplce_19_29$par[3])
q_theo_19_29_agg <- qagg(p = p,
                         loc = fit_agg_19_29$par[1],
                         scale_1 = fit_agg_19_29$par[2],
                         scale_2 = fit_agg_19_29$par[3],
                         shape = fit_agg_19_29$par[4])

plot(q_emp,
     q_theo_19_29_dlaplce,
     xlab = "Empirical",
     ylab = "Theoretical",
     ylim = c(min(q_emp, q_theo_19_29_dlaplce, q_theo_19_29_agg),
              max(q_emp, q_theo_19_29_dlaplce, q_theo_19_29_agg)))
points(q_emp,
       q_theo_19_29_agg,
       col = "blue")
abline(a = 0, b = 1, col = 2)

mean(abs(q_emp - q_theo_19_29_agg))
mean(abs(q_emp - q_theo_19_29_dlaplce))

sqrt(mean(((q_emp - q_theo_19_29_agg)^2)))
sqrt(mean(((q_emp - q_theo_19_29_dlaplce)^2)))

## PP plot
p_theo_19_29_dlaplce <- pdlaplace(z = sort(dat_19_29),
                                  mu = fit_ddlaplce_19_29$par[1],
                                  sigma = fit_ddlaplce_19_29$par[2],
                                  delta = fit_ddlaplce_19_29$par[3])
p_theo_19_29_agg <- pagg(q = sort(dat_19_29),
                         loc = fit_agg_19_29$par[1],
                         scale_1 = fit_agg_19_29$par[2],
                         scale_2 = fit_agg_19_29$par[3],
                         shape = fit_agg_19_29$par[4])

plot(p,
     p_theo_19_29_dlaplce,
     xlab = "Empirical",
     ylab = "Theoretical",
     xlim = c(0, 1),
     ylim = c(0, 1))
points(p,
       p_theo_19_29_agg,
       col = "blue")
abline(a = 0, b = 1, col = 2)

mean(abs(p - p_theo_19_29_agg))
mean(abs(p - p_theo_19_29_dlaplce))

sqrt(mean(((p - p_theo_19_29_agg)^2)))
sqrt(mean(((p - p_theo_19_29_dlaplce)^2)))

