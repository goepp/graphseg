# Note: This code is deprecated and should not be used. The main
# function of the package is "agraph".
#' Vectorized \code{\link{str_detect}} from package stringr
#'
#' @param string string to check for presence of the pattern
#' @param pattern_vect pattern to detec
#' @return A vector of boolean of the same length as \code{string}
#' @author Vivien Goepp \email{vivien.goepp@@gmail.com}
#' @seealso \code{\link{str_detect}} from package \code{stringr}
#' @importFrom stringr str_detect
#' @importFrom mgcv bam
#' @importFrom stats setNames gaussian
#' @import graphics
#' @importFrom igraph decompose graph_from_adjacency_matrix clusters V E spectrum
#' @import magrittr
#' @import methods
str_detect_2 <- function(string, pattern_vect) {
  sapply(pattern_vect, str_detect, string = string) %>%
    apply(1, any)
}
sel_columns <- function(X, t) {
  Xd <- as.data.frame(X)
  X_new <- c()
  for (ind in 1:max(t)) {
    # colnum <- which(names(Xd) %in% names(t)[which(t == ind)])
    colnum <- str_detect_2(names(Xd), names(t)[which(t == ind)]) %>% which()
    temp <- Xd[colnum]
    X_new <- cbind(X_new, Reduce(`+`, temp))
  }
  unname(X_new)
}
# solver_aridge_lambda <- function(y, X, touches, lambda = 10 ^ -0.5, maxiter = 100, epsilon = 1e-5,
#                                  thresh = 1e-8) {
#   adjacency <- touches
#   theta <- rep(0, ncol(X))
#   iter <- 1
#   sel <- adjacency
#   while (iter < maxiter) {
#     old_sel <- sel
#     K <- diag(apply(adjacency, 1, sum)) - adjacency
#     p <- 1 / (1 + exp(-X %*% theta))
#     Gamma <- 1 / p + 1 / (1 - p)
#     z <- X %*% theta + Gamma * (y - p)
#     temp <- as.vector(exp(-X %*% theta))
#     w <- temp / (1 + temp) ^ 2
#     reg <- gam(z ~ X, family = gaussian(), data = data.frame(z = z, X = X),
#                drop.intercept = TRUE, weights = w, H = lambda * K)
#     theta <- reg$coefficients
#     gamma <- theta
#     adjacency <- touches *
#       ((replicate(length(gamma), gamma) - t(replicate(length(gamma), gamma))) ^ 2 + epsilon ^ 2) ^ (-1)
#     sel <- adjacency *
#       ((replicate(length(gamma), gamma) - t(replicate(length(gamma), gamma))) ^ 2 + epsilon ^ 2)
#     converge <- max(abs(old_sel - sel)) <= thresh
#     converge
#     if (converge) {
#       components <- graph_from_adjacency_matrix((sel < 0.99) * touches)
#       t <- clusters(components)$membership
#       X_sel <- sel_columns(X, t)
#       fit <- gam(y ~ X_sel, family = "binomial", drop.intercept = TRUE)
#       # temp <- function(ind) rep(fit$coefficients[ind], sum(t == ind))
#       temp2 <- function(ind) fit$coefficients[t[ind]]
#       par <- sapply(1:length(t), temp2) %>%
#         setNames(colnames(X))
#       break
#     }
#     iter <- iter + 1
#   }
#   par
# }
# solver_aridge <- function(
#   y, X, adj, lambda = 10 ^ seq(-2, 2, length = 50),
#   maxiter = 100, epsilon = 1e-5, thresh = 1e-8
# ) {
#   par_ls <- vector("list", length = length(lambda))
#   weight_adj <- adj
#   theta <- rep(0, ncol(X))
#   iter <- 1
#   ind <- 1
#   sel <- adjacency
#   while (iter < maxiter) {
#     old_sel <- sel
#     K <- diag(apply(weight_adj, 1, sum)) - adja
#     p <- 1 / (1 + exp(-X %*% theta))
#     Gamma <- 1 / p + 1 / (1 - p)
#     z <- X %*% theta + Gamma * (y - p)
#     temp <- as.vector(exp(-X %*% theta))
#     w <- temp / (1 + temp) ^ 2
#     reg <- bam(z ~ X, family = gaussian(), data = data.frame(z = z, X = X),
#                drop.intercept = TRUE, weights = w, H = lambda * K)
#     theta <- reg$coefficients
#     weight_adj <- adj *
#       ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2) ^ (-1)
#     sel <- weight_adj *
#       ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2)
#     converge <- max(abs(old_sel - sel)) <= thresh
#     if (converge) {
#       components <- graph_from_adjacency_matrix((sel < 0.99) * adj)
#       t <- clusters(components)$membership
#       X_sel <- sel_columns(X, t)
#       fit <- bam(y ~ X_sel, family = "binomial", drop.intercept = TRUE)
#       par_ls[[ind]] <- sapply(seq_along(t), function(ind) fit$coefficients[t[ind]]) %>%
#         setNames(colnames(X))
#       ind <-  ind + 1
#       break
#     }
#     iter <- iter + 1
#     if (ind > length(par_ls)) break
#   }
#   par
# }
# solver_graph_aridge_old <- function(gamma, sigma_sq, adj,
#                                     lambda = 10 ^ seq(-4, 4, length = 100),
#                                     maxiter = 1000,
#                                     epsilon = 1e-5, thresh = 1e-8) {
#   X <-  diag(length(gamma)) %>%
#     # Matrix(sparse = TRUE) %>%
#     "colnames<-"(colnames(adj))
#   weight_adj <- adj
#   theta <- gamma * 0
#   par_ls <- vector("list", length(lambda))
#   iter <- 1
#   ind <- 1
#   sel <- weight_adj
#   while (iter < maxiter) {
#     old_sel <- sel
#     K <- diag(apply(weight_adj, 1, sum)) - weight_adj
#     reg <- gam(gamma ~ X - 1, family = gaussian(),
#                data = data.frame(gamma = gamma, X = X),
#                drop.intercept = TRUE, weights = 1 / sigma_sq, H = (lambda[ind] * K) %>% as.matrix())
#     theta <- reg$coefficients
#     weight_adj <- adj *
#       ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2) ^ (-1)
#     sel <- weight_adj * ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2)
#     converge <- max(abs(old_sel - sel)) <= thresh
#     converge
#     if (converge) {
#       components <- graph_from_adjacency_matrix((sel < 0.99) * adj, mode = "undirected")
#       t <- clusters(components)$membership
#       X_sel <- sel_columns(X, t)
#       fit <- mgcv::gam(gamma ~ X_sel, family = gaussian(), drop.intercept = TRUE)
#       temp <- function(ind) fit$coefficients[t[ind]]
#       par_ls[[ind]] <- sapply(1:length(t), temp) %>%
#         setNames(colnames(X))
#       # weight_adj <- adj
#       # sel <- weight_adj
#       ind <- ind + 1
#     }
#     iter <- iter + 1
#     if (ind > length(par_ls)) break
#   }
#   par_ls
# }
loglik_old <- function(par, gamma, sigma_sq, K, pen) {
  pen_mat <- (t(replicate(length(par), par)) - replicate(length(par), par)) ^ 2
  sum((par - gamma) ^ 2 / (2 * sigma_sq)) - pen / 2 * sum(pen_mat * K)
}
graph_aridge_old <- function(gamma, sigma_sq, adj,
                             pen = 10 ^ seq(-4, 4, length = 100),
                             maxiter = 10000,
                             epsilon = 1e-5,
                             thresh = 1e-8) {
  par_ls <- vector("list", length(pen))
  bic <- bic_1 <- bic_2 <- laplace_bic <- pen * 0
  X <- "colnames<-"(diag(length(gamma)), colnames(adj))
  weight_adj <- adj
  theta <- gamma * 0
  iter <- 1
  ind <- 1
  sel <- weight_adj
  while (iter < maxiter) {
    old_sel <- sel
    K <- diag(apply(weight_adj, 1, sum)) - weight_adj
    theta <- solver_nr(gamma, sigma_sq, K, pen = pen[ind])
    weight_adj <- adj *
      ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2) ^ (-1)
    sel <- weight_adj * ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2)
    converge <- max(abs(old_sel - sel)) <= thresh
    converge
    if (converge) {
      components <- graph_from_adjacency_matrix((sel < 0.99) * adj, mode = "undirected")
      t <- clusters(components)$membership
      X_sel <- sel_columns(X, t)
      fit <- mgcv::gam(gamma ~ X_sel, family = stats::gaussian(), drop.intercept = TRUE)
      temp <- function(ind) fit$coefficients[t[ind]]
      par_ls[[ind]] <- setNames(sapply(1:length(t), temp), colnames(X))
      bic[ind] <- log(length(gamma)) * ncol(X_sel) +
        2 * loglik(par_ls[[ind]], gamma, sigma_sq, K, pen = 0)
      bic_1[ind] <- log(length(gamma)) * ncol(X_sel)
      bic_2[ind] <- 2 * loglik(par_ls[[ind]], gamma, sigma_sq, K, pen = 0)
      laplace_bic[ind] <- bic[ind] - ncol(X_sel) * log(2 * pi) - sum(log(sigma_sq))
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(par_ls)) break
  }
  list(par = par_ls, bic = bic, laplace_bic = laplace_bic,
       bic_1 = bic_1, bic_2 = bic_2)
}