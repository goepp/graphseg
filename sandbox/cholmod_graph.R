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
opar <- par(mar = rep(0, 4))

data <- graphseg::departement

p <- length(data$geometry)
kappa <- 1e0
delta <- 1e-10
centroids <- st_transform(data, 29101) %>%
  st_centroid() %>%
  # this is the crs from d, which has no EPSG code:
  st_transform(., '+proj=longlat +ellps=GRS80 +no_defs') %>%
  # since you want the centroids in a second geometry col:
  st_geometry()
layout_graph <- do.call(rbind, st_geometry(centroids))

neighbour <- poly2nb(data)
listw <- nb2listw(neighbour, style = "B", zero.policy = TRUE)
adj <- as(as(listw, "CsparseMatrix"), "symmetricMatrix")
sel <- adj
graph <- graph_from_adjacency_matrix(adj, mode = "undirected")



set.seed(0)
gamma <- sample(1:p)
weights <- rep(1, p)

L <- Diagonal(x = colSums(adj)) - adj
L_chol <- Cholesky(L, Imult = 1 / kappa)
class(L_chol)
sparsechol <- solve(L_chol, gamma, system = "A") / kappa
sparse <- solve(kappa * L + diag(p), gamma)
all.equal(sparsechol, sparse)

edgelist <- as_edgelist(graph, names = FALSE)

adj@x <- 1 / ((sparsechol[edgelist[, 1]] - sparsechol[edgelist[, 2]]) ^ 2 + delta)
sel@x <- (sparsechol[edgelist[, 1]] - sparsechol[edgelist[, 2]]) ^ 2 /
  ((sparsechol[edgelist[, 1]] - sparsechol[edgelist[, 2]]) ^ 2 + delta)
L <- Diagonal(x = colSums(adj)) - adj
L_chol <- Cholesky(L, Imult = 1 / kappa)
sparsechol <- solve(L_chol, gamma, system = "A") / kappa

sel_graph <- graph_from_adjacency_matrix(sel, mode = "undirected", weighted = TRUE)
graph_col <- heat.colors(fine)[as.numeric(cut(as.vector(sparsechol),
                                              breaks = quantile(gamma, seq(0, 1, length = fine + 1)),
                                              include.lowest = TRUE))]
plot(sel_graph, edge.width = E(sel_graph)$weight + 0.5, edge.arrow.size = 0.5,
     vertex.size = 10, vertex.label.cex = 0.001,
     vertex.color = graph_col, layout = layout_graph)
