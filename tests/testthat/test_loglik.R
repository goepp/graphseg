context("loglik.R")

library(gdata)
library(glmnet)
library(Matrix)
library(limSolve)

library(microbenchmark)
library(car)
library(tidyverse)
library(igraph)
library(graphseg)
library(sf)

# data(departement)
adj <- sf::st_intersects(departement) %>%
  as.matrix() %>%
  '[<-'(cbind(1:nrow(.), 1:ncol(.)), 0)
adj_graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")

set.seed(0)
departement$gamma <- rnorm(length(departement$code))

weight <- rep(1, length(igraph::E(adj_graph)))
sigma_sq <- rep(1, length(departement$gamma))
pen <- 1
K <- igraph::laplacian_matrix(adj_graph, weights = weight, sparse = FALSE)

test_that("loglik and loglik_old are equal",
          expect_equal(loglik(departement$gamma, departement$gamma, sigma_sq, K, pen),
                       loglik_old(departement$gamma, departement$gamma, sigma_sq, K, pen)))

## Remark: It runs faster
microbenchmark(loglik(departement$gamma, departement$gamma, sigma_sq, K, pen),
               loglik_old(departement$gamma, departement$gamma, sigma_sq, K, pen))
