################################################################################
rm(list = ls())
required_pckgs <- c("gtools", "igraph")
# install.packages(required_pckgs, dependencies = TRUE)
t(t(sapply(required_pckgs, require, character.only = TRUE)))

################################################################################

## True graph
d <- 5
c1 <- 1:3
c2 <- 3:5
cliques_true <- list(c1, c2)
edges <- do.call(rbind, lapply(cliques_true, function(x){combinations(length(x), r = 2, v = x)}))
edges <- t(edges[which(duplicated(edges) == FALSE),])
g_true <- graph(edges = edges, directed = FALSE)

## construct the layout
layout <- rbind(c(-4, -4),
                c(-4, 4),
                c(0, 0),
                c(4, -4),
                c(4, 4))
## plotting parameters
V(g_true)$size <- 40

## Obtain the subgraph
i <- 1
g_cond <- delete_vertices(g_true, i)
V(g_cond)$label <- (1:d)[-i]

## Plot the output
pdf(file = "Images/Methods/Graph.pdf", width = 10, height = 10)
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot.igraph(g_true, layout = layout, edge.width = 25, vertex.label.cex = 5)
dev.off()

pdf(file = "Images/Methods/Graph_Cond.pdf", width = 10, height = 10)
par(mfrow = c(1, 1), mgp = c(2.3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
plot.igraph(g_cond, layout = layout[-i,], edge.width = 25, vertex.label.cex = 5)
dev.off()
