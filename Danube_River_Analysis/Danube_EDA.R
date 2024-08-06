################################################################################
#Load in required packages
rm(list = ls())
required_pckgs <- c("evd", "ggplot2", "ggpubr", "graphicalExtremes", "igraph", "mev")
# install.packages(required_pckgs, dependencies = TRUE, Ncpus = detectCores() - 1)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################
## Plot the upper Danube River basin
danube_graph <- graph(t(danube$flow_edges), directed = FALSE)
layout <- cbind(danube$info$PlotCoordX, danube$info$PlotCoordY)
V(danube_graph)$size <- 15
V(danube_graph)$label.cex <- 2

pdf(file = "Images/Danube/Inference_On_Data/Danube_River.pdf", width = 10, height = 10)
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot(danube_graph, layout = layout, edge.width = 10)
dev.off()
################################################################################
## Investigate the dependence in the river network
data_frechet <- apply(apply(danube$data_clustered, 2, rank)/(nrow(danube$data_clustered) + 1), 2, qfrechet)

ai_pairs <- rbind(c(12, 23),
                  c(13, 23),
                  c(14, 28), c(14, 29), c(14, 30), c(14, 31),
                  c(15, 28), c(15, 29), c(15, 30), c(15, 31),
                  c(16, 31), 
                  c(17, 27), c(17, 31),
                  c(18, 31),
                  c(19, 31),
                  c(20, 23), c(20, 24), c(20, 25), c(20, 26), c(20, 27), c(20, 31),
                  c(23, 28), c(23, 29), c(23, 30), c(23, 31))

for(i in 1:nrow(ai_pairs)){
  chi_test <- texmex::chi(data = danube$data_clustered[,ai_pairs[i,]])
  par(mfrow = c(1, 4), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
  
  plot(x = data_frechet[,ai_pairs[i,1]], y = data_frechet[,ai_pairs[i,2]], pch = 19,
       xlab = paste("Site", ai_pairs[i,1]), ylab = paste("Site", ai_pairs[i,2]))
  
  plot(x = chi_test$quantile, y = chi_test$chibar[,2], 
       xlim = c(0, 1), ylim = c(-1, 1), type = "l", 
       xlab = "Quantile", ylab = substitute(bar(chi)(u)))
  lines(x = chi_test$quantile, y = chi_test$chibar[,1], lty = "dashed")
  lines(x = chi_test$quantile, y = chi_test$chibar[,3], lty = "dashed")
  
  plot(x = chi_test$quantile, y = chi_test$chi[,2], 
       xlim = c(0, 1), ylim = c(-1, 1), type = "l", 
       xlab = "Quantile", ylab = substitute(chi(u)))
  lines(x = chi_test$quantile, y = chi_test$chi[,1], lty = "dashed")
  lines(x = chi_test$quantile, y = chi_test$chi[,3], lty = "dashed")
  
  eta_test <- taildep(data = danube$data_clustered[,ai_pairs[i,]], u = chi_test$quantile,
                      plot = FALSE)
  plot(x = chi_test$quantile, y = eta_test$eta[,1], 
       xlim = c(0, 1), ylim = c(-1, 1), type = "l", 
       xlab = "Quantile", ylab = substitute(eta(u)))
  lines(x = chi_test$quantile, y = eta_test$eta[,2], lty = "dashed")
  lines(x = chi_test$quantile, y = eta_test$eta[,3], lty = "dashed")
}

## display sites 19 and and 29 and 19 and 31

## 19 and 29
dep_plots <- list()
dep_plots[[1]] <- ggplot(data = data.frame(x = data_frechet[,19], y = data_frechet[,29]), aes(x = x, y = y)) +
  geom_point() +
  labs(x = "Station 19", y = "Station 29")

chi_test <- texmex::chi(data = danube$data_clustered[,c(19,29)])
chi_test$chiulb[is.infinite(chi_test$chiulb)] <- NA
dep_plots[[2]] <- ggplot() + geom_polygon(data = data.frame(x = c(chi_test$quantile, rev(chi_test$quantile)),
                                                            y = c(chi_test$chibar[,1], rev(chi_test$chibar[,3]))),
                                          aes(x = x, y = y), fill = alpha(colour = "orange", alpha = 0.5)) +
  geom_line(data = data.frame(x = chi_test$quantile, y = chi_test$chibar[,2]), aes(x = x, y = y), col = "blue") +
  geom_line(data = data.frame(x = chi_test$quantile, y = -1), aes(x = x, y = y), col = "black", linetype = "dashed") +
  geom_line(data = data.frame(x = chi_test$quantile, y = 0), aes(x = x, y = y), col = "black", linetype = "dashed") +
  geom_line(data = data.frame(x = chi_test$quantile, y = 1), aes(x = x, y = y), col = "black", linetype = "dashed") +
  lims(x = c(0, 1), y = c(-1, 1)) +
  labs(x = "Quantile", y = expression(bar(chi)(u)))

dep_plots[[3]] <- ggplot() + geom_polygon(data = data.frame(x = c(chi_test$quantile, rev(chi_test$quantile)),
                                                            y = c(chi_test$chi[,1], rev(chi_test$chi[,3]))),
                                          aes(x = x, y = y), fill = alpha(colour = "orange", alpha = 0.5)) +
  geom_line(data = data.frame(x = chi_test$quantile, y = chi_test$chi[,2]), aes(x = x, y = y), col = "blue") +
  geom_line(data = data.frame(x = chi_test$quantile, y = chi_test$chiulb), aes(x = x, y = y), col = "black", linetype = "dashed") +
  geom_line(data = data.frame(x = chi_test$quantile, y = 0), aes(x = x, y = y), col = "black", linetype = "dashed") +
  geom_line(data = data.frame(x = chi_test$quantile, y = 1), aes(x = x, y = y), col = "black", linetype = "dashed") +
  lims(x = c(0, 1), y = c(-1, 1)) +
  labs(x = "Quantile", y = expression(chi(u)))

## 19 and 31
dep_plots[[4]] <- ggplot(data = data.frame(x = data_frechet[,19], y = data_frechet[,31]), aes(x = x, y = y)) +
  geom_point() +
  labs(x = "Station 19", y = "Station 31")

chi_test <- texmex::chi(data = danube$data_clustered[,c(19,31)])
chi_test$chiulb[is.infinite(chi_test$chiulb)] <- NA
dep_plots[[5]] <- ggplot() + geom_polygon(data = data.frame(x = c(chi_test$quantile, rev(chi_test$quantile)),
                                                            y = c(chi_test$chibar[,1], rev(chi_test$chibar[,3]))),
                                          aes(x = x, y = y), fill = alpha(colour = "orange", alpha = 0.5)) +
  geom_line(data = data.frame(x = chi_test$quantile, y = chi_test$chibar[,2]), aes(x = x, y = y), col = "blue") +
  geom_line(data = data.frame(x = chi_test$quantile, y = -1), aes(x = x, y = y), col = "black", linetype = "dashed") +
  geom_line(data = data.frame(x = chi_test$quantile, y = 0), aes(x = x, y = y), col = "black", linetype = "dashed") +
  geom_line(data = data.frame(x = chi_test$quantile, y = 1), aes(x = x, y = y), col = "black", linetype = "dashed") +
  lims(x = c(0, 1), y = c(-1, 1)) +
  labs(x = "Quantile", y = expression(bar(chi)(u)))

dep_plots[[6]] <- ggplot() + geom_polygon(data = data.frame(x = c(chi_test$quantile, rev(chi_test$quantile)),
                                                            y = c(chi_test$chi[,1], rev(chi_test$chi[,3]))),
                                          aes(x = x, y = y), fill = alpha(colour = "grey", alpha = 0.5)) +
  geom_line(data = data.frame(x = chi_test$quantile, y = chi_test$chi[,2]), aes(x = x, y = y), col = "grey") +
  geom_line(data = data.frame(x = chi_test$quantile, y = chi_test$chiulb), aes(x = x, y = y), col = "black", linetype = "dashed") +
  geom_line(data = data.frame(x = chi_test$quantile, y = 0), aes(x = x, y = y), col = "black", linetype = "dashed") +
  geom_line(data = data.frame(x = chi_test$quantile, y = 1), aes(x = x, y = y), col = "black", linetype = "dashed") +
  lims(x = c(0, 1), y = c(-1, 1)) +
  labs(x = "Quantile", y = expression(chi(u)))

par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
pdf(file = "Images/Danube/Inference_On_Data/EDA.pdf", width = 10, height = 10)
ggarrange(plotlist = dep_plots, nrow = 2, ncol = 3)
dev.off()
################################################################################

