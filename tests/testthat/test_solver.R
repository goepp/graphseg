context("loglik.R")

library(tidyverse)
library(gdata)
library(RColorBrewer)
library(microbenchmark)
library(rgeos)
library(rgdal)
library(lattice)
library(glmnet)
library(igraph)
library(mgcv)
library(Matrix)
library(limSolve)
library(graphseg)
library(sf)

adj <- st_intersects(departement) %>%
  as.matrix() %>%
  '[<-'(cbind(1:nrow(.), 1:ncol(.)), 0)
adj_graph <- graph_from_adjacency_matrix(adj, mode = "undirected")

set.seed(0)
departement$gamma <- rnorm(length(departement$code))

weight <- rep(1, length(E(adj_graph)))
sigma_sq <- rep(1, length(departement$gamma))
pen <- 1
K <- laplacian_matrix(adj_graph, weights = weight, sparse = FALSE)

test_that("solver and solver_nr are equal",
          expect_equal(solver(departement$gamma, sigma_sq, K, pen),
                       solver_nr(departement$gamma, sigma_sq, K, pen)))

## Remark: The explicit solver is way faster
microbenchmark(solver(departement$gamma, sigma_sq, K, pen),
               solver_nr(departement$gamma, sigma_sq, K, pen))
