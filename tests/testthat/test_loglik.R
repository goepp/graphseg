context("loglik.R")

library(tidyverse)
library(gdata)
library(RColorBrewer)
library(microbenchmark)
library(car)
library(rgeos)
library(rgdal)
library(lattice)
library(glmnet)
library(igraph)
library(sf)
library(profvis)
library(mgcv)
library(Matrix)
library(pryr)
library(limSolve)
library(graphseg)

adj <- st_intersects(departement) %>%
  as.matrix() %>%
  '[<-'(cbind(1:nrow(.), 1:ncol(.)), 0)
adj_graph <- graph_from_adjacency_matrix(adj, mode = "undirected")

set.seed(0)
departement$gamma <- rnorm(length(departement$code))

weight <- rep(1, length(E(adj_graph)))
sigma_sq <- 1
pen <- 1
K <- laplacian_matrix(adj_graph, weights = weight, sparse = FALSE)

test_that("loglik and loglik_old are equal",
          expect_equal(loglik(departement$gamma, departement$gamma, sigma_sq, K, pen),
                     loglik_old(departement$gamma, departement$gamma, sigma_sq, K, pen)))

## Remark: It runs faster
microbenchmark(loglik(departement$gamma, departement$gamma, sigma_sq, K, pen),
               loglik_old(departement$gamma, departement$gamma, sigma_sq, K, pen))
