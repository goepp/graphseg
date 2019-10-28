setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(graphseg)
library(Matrix)
library(sf)
library(spdep)
library(spatialreg)
library(igraph)
library(microbenchmark)
library(tidyverse)
library(profvis)
library(profmem)

## Plot options
fine <- 500
pal <- colorRampPalette(c("gray80", "dark red"))
# opar <- par(mar = rep(0, 4))

lambda <- 1e-4
delta <- 1e-5
tol <- 1e-10
thresh <- 0.01
graph <- make_lattice(c(10, 10))
p <- length(V(graph))
layout_graph <- layout_on_grid(graph)
# edgelist <- as_edgelist(graph, names = FALSE)
edgelist_tmp <- as_edgelist(graph, names = FALSE)
edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
# edgelist <- edgelist[order(edgelist[, 2]), ]
adj <- as(as_adjacency_matrix(graph), "symmetricMatrix")
sel <- sel_old <- adj

set.seed(0)
# gamma <- sample(1:p)
gamma_pc <- c(rep(0, floor(p / 2)), rep(1, p - floor(p / 2)))
filter <- expm((-0.5) * (Diagonal(x = colSums(adj)) - adj))
gamma <- filter %*% gamma_pc + rnorm(p, 0, 0.01)
graph_col <- heat.colors(fine)[as.numeric(cut(as.vector(gamma),
                                              breaks = quantile(gamma, seq(0, 1, length = fine + 1)),
                                              include.lowest = TRUE))]
plot(graph, edge.width = 1,
     vertex.size = 5, vertex.label.cex = 0.0001,
     vertex.color = graph_col, layout = layout_graph)
weights <- rep(1, p)
# weights <- c(rep(1, p / 2), rep(0, p / 2))
weighted_laplacian_init <- lambda * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
chol_init <- Cholesky(weighted_laplacian)

weighted_laplacian <- lambda * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
chol <- update(chol_init, weighted_laplacian)
theta <- solve(chol, weights * gamma)
adj@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
  ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
all(abs(sel_old@x - sel@x) < tol)
max(abs(sel_old@x - sel@x))
sel_old <- sel
plot(sel@x)

sel_graph <- graph_from_adjacency_matrix(sel, mode = "undirected", weighted = TRUE)
graph_col <- heat.colors(fine)[as.numeric(cut(as.vector(theta),
                                              breaks = quantile(gamma, seq(0, 1, length = fine + 1)),
                                              include.lowest = TRUE))]
plot(sel_graph, edge.width = (1 - E(sel_graph)$weight) * 10, edge.arrow.size = 0.5,
     vertex.size = 5, vertex.label.cex = 0.0001,
     vertex.color = graph_col, layout = layout_graph)

graph_del <- graph_from_adjacency_matrix((sel > 1 - thresh), mode = "undirected")
graph_del <- graph %>% delete_edges(which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
theta <- ave(as.vector(gamma), components(graph_del)$membership)

res <- sapply(10 ^ seq(-5, 4, 0.1), agraph_one_lambda, gamma = gamma, graph = graph,
              weights = weights) %>% t()
matplot(res, type = "l")
colSums(res) %>% plot()
points(gamma, col = "red")

res2 <- agraph(gamma, graph)
res2$result %>% t() %>% matplot(type = "l")
plot(res2$bic)
