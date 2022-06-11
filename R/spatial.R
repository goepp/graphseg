#' Segmentation using graph structure
#' @param gamma entry vector to regularize
#' @param graph \url{igraph} giving the regularization structure
#' @param lambda regularizing constant
#' @param weights weights for gamma. Default value is one.
#' @param delta Computational constant in the adaptive ridge reweighting formula.
#' @param tol Tolerance to test for convergence of the adaptive ridge
#' @param thresh Thresholding constant used to fuse two adjacent regions with close value of \code{gamma}.
agraph_one_lambda <- function(gamma, graph, lambda = 1, weights = NULL,
                              delta = 1e-10, tol = 1e-8,
                              thresh = 0.01) {
  gamma <- as.vector(gamma)
  lambda <- as.vector(lambda)
  p <- length(gamma)
  if (is.null(weights)) {
    weights <- rep(1, p)
  }
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be in igraph format")
  }
  edgelist_tmp <- igraph::as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adj <- Matrix::forceSymmetric(igraph::as_adjacency_matrix(graph))
  sel <- sel_old <- adj
  converge <- FALSE
  weighted_laplacian_init <- lambda * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
  chol_init <- Cholesky(weighted_laplacian_init)
  while (!converge) {
    weighted_laplacian <- lambda * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
    chol <- update(chol_init, weighted_laplacian)
    theta <- solve(chol, weights * gamma)
    adjacency@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
      ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    converge <- all(abs(sel@x - sel_old@x) < tol)
    sel_old <- sel
  }
  graph_del <- igraph::delete_edges(graph, which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
  theta <- stats::ave(as.vector(gamma), igraph::components(graph_del)$membership)
  return(theta)
}
#' Segmentation using graph structure
#' @param gamma entry vector to regularize
#' @param graph graph (an \url{igraph} object) giving the regularization structure
#' @param lambda regularizing constant
#' @param shrinkage Boolean, defaults TRUE. Whether to return the adaptive ridge estimate as output. If FALSE, the adaptive ridge is used to define a segmentation into zones, and the signal is estimated on each zone using non-penalized estimation.
#' @param weights weights for gamma. Default value is one.
#' @param delta Computational constant in the adaptive ridge reweighting formula.
#' @param tol Tolerance to test for convergence of the adaptive ridge
#' @param thresh Thresholding constant used to fuse two adjacent regions with close value of \code{gamma}.
#' @param itermax Total number of iterations. Default value is 10000. Setting a low value can make the procedure return NULL entries for some values of \code{lambda}.
#' @importFrom Matrix Cholesky solve update Diagonal
#' @export
agraph <- function(gamma, graph, lambda = 10 ^ seq(-4, 4, length.out = 50),
                   weights = NULL, shrinkage = TRUE,
                   delta = 1e-10, tol = 1e-8,
                   thresh = 0.01, itermax = 50000) {
  p <- length(gamma)
  if (is.null(weights)) {
    weights <- rep(1, p)
  }
  prec <- Matrix::Diagonal(x = weights)
  nll <- model_dim <- rep(0, length(lambda))
  result <- matrix(NA, length(lambda), p)
  ind <- 1
  iter <- 1
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be an igraph object.")
  }
  if (length(gamma) != igraph::vcount(graph)) {
    stop("gamma must be a vector of length the number of vertices in the graph.")
  }
  edgelist_tmp <- igraph::as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adj <- Matrix::forceSymmetric(igraph::as_adjacency_matrix(graph))
  sel <- adj
  converge <- FALSE
  weighted_laplacian_init <- lambda[ind] * (Matrix::Diagonal(x = Matrix::colSums(adj)) - adj) + prec
  chol_init <- Matrix::Cholesky(weighted_laplacian_init)
  while (iter < itermax) {
    sel_old <- sel
    weighted_laplacian <- lambda[ind] * (Diagonal(x = Matrix::colSums(adj)) - adj) + prec
    chol <- Matrix::update(chol_init, weighted_laplacian)
    theta <- Matrix::solve(chol, weights * gamma)
    adj@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
      ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    converge <- all(abs(sel@x - sel_old@x) < tol)
    if (converge) {
      graph_del <- igraph::delete_edges(graph, which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
      segmentation <- igraph::components(graph_del)$membership
      if (shrinkage) {
        result[ind, ] <- stats::ave(as.vector(theta), segmentation)
      } else {
        result[ind, ] <- stats::ave(as.vector(gamma), segmentation)
      }
      nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% prec %*% (result[ind, ] - gamma)
      model_dim[ind] <- sum(diag(Matrix::solve(weighted_laplacian, prec)))
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(lambda)) break
    # if (model_dim == 1) break
  }
  bic <- 2 * nll + log(p) * model_dim
  aic <- 2 * nll + 2 * model_dim
  gcv <- 2 * nll / (p * (1 - model_dim / p) ^ 2)
  return(list(result = result, bic = bic, gcv = gcv, model_dim = model_dim,
              nll = nll, aic = aic))
}
rgraph <- function(gamma, graph,
                   lambda = 10 ^ seq(-4, 4, length.out = 50),
                   weights = NULL) {
  gamma <- as.vector(gamma)
  lambda <- as.vector(lambda)
  p <- length(gamma)
  if (is.null(weights)) {
    weights <- rep(1, p)
  }
  prec <- Diagonal(x = weights)
  nll <- model_dim <- rep(0, length(lambda))
  result <- matrix(NA, length(lambda), p)
  ind <- 1
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be in igraph format")
  }
  adj <- Matrix::forceSymmetric(igraph::as_adjacency_matrix(graph))
  weighted_laplacian_init <- lambda[ind] * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
  chol_init <- Cholesky(weighted_laplacian_init)
  for (ind in seq_along(lambda)) {
    weighted_laplacian <- lambda[ind] * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
    chol <- update(chol_init, weighted_laplacian)
    result[ind, ] <- as.vector(solve(chol, weights * gamma))
    nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% prec %*% (result[ind, ] - gamma)
    model_dim[ind] <- sum(diag(Matrix::solve(weighted_laplacian, prec)))
  }
  bic <- 2 * nll + log(p) * model_dim
  aic <- 2 * nll + 2 * model_dim
  gcv <- 2 * nll / (p * (1 - model_dim / p) ^ 2)
  return(list(result = result, bic = bic, gcv = gcv, model_dim = model_dim,
              nll = nll, aic = aic))
}
#' @rdname agraph
#' @param prec precision matrix (inverse of the variance-covariance matrix). Has to be a sparse matrix for efficiency.
#' @export
agraph_prec <- function(gamma, graph, prec,
                        lambda = 10 ^ seq(-4, 4, length.out = 50),
                        weights = NULL, shrinkage = TRUE,
                        delta = 1e-10, tol = 1e-8,
                        thresh = 0.01, itermax = 10000) {
  gamma <- as.vector(gamma)
  lambda <- as.vector(lambda)
  if (!"sparseMatrix" %in% methods::is(prec)) {
    prec <- Matrix::Matrix(prec, sparse = TRUE)
    warning("Precision matrix was not sparse. It was coerced to sparse format.")
  }
  if (!isSymmetric(prec)) {
    prec <- Matrix::forceSymmetric(prec)
    warning("Precision matrix was not symmetric. It was coerced to symmetric.")
  }
  p <- length(gamma)
  prec_gamma <- prec %*% gamma
  nll <- model_dim <- rep(0, length(lambda))
  result <- matrix(NA, length(lambda), p)
  ind <- 1
  iter <- 1
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be in igraph format")
  }
  edgelist_tmp <- igraph::as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adj <- Matrix::forceSymmetric(igraph::as_adjacency_matrix(graph))
  sel <- adj
  converge <- FALSE
  weighted_laplacian_init <- lambda[ind] * (Matrix::Diagonal(x = colSums(adj)) - adj) + prec
  chol_init <- Matrix::Cholesky(weighted_laplacian_init)
  while (iter < itermax) {
    sel_old <- sel
    weighted_laplacian <- lambda[ind] * (Matrix::Diagonal(x = colSums(adj)) - adj) + prec
    chol <- Matrix::update(chol_init, weighted_laplacian)
    theta <- Matrix::solve(chol, prec_gamma)
    adj@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
      ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    converge <- all(abs(sel@x - sel_old@x) < tol)
    if (converge) {
      graph_del <- igraph::delete_edges(graph, which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
      segmentation <- igraph::components(graph_del)$membership
      if (shrinkage) {
        result[ind, ] <- stats::ave(as.vector(theta), segmentation)
      } else {
        result[ind, ] <- stats::ave(as.vector(gamma), segmentation)
      }
      nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% prec %*% (result[ind, ] - gamma)
      model_dim[ind] <- sum(diag(Matrix::solve(weighted_laplacian, prec)))
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(lambda)) break
  }
  bic <- 2 * nll + log(p) * model_dim
  aic <- 2 * nll + 2 * model_dim
  gcv <- 2 * nll / ((p - model_dim) ^ 2)
  return(list(result = result, aic = aic, bic = bic,
              gcv = gcv, model_dim = model_dim, nll = nll))
}
#' @importFrom magrittr %>%
graph2connlist <- function(graph) {
  connlist <- igraph::as_adj_list(graph) %>%
    lapply(as.character) %>%
    lapply(as.numeric) %>%
    lapply(function(a) a - 1) %>%
    lapply(function(a) {
      attributes(a) <- NULL
      return(a)
    }) %>%
    lapply(as.integer) %>%
    'class<-'("connListObj")
  return(connlist)
}
#' Segmentation using graph structure and the fused lasso estimate
#' @description Wrapper around the function \code{\link[flsa:flsa]{flsa::flsa}}.
#' @param gamma entry vector to regularize
#' @param graph graph (an \url{igraph} object) giving the regularization structure
#' @param lambda regularizing constant
#' @importFrom Matrix Cholesky solve update Diagonal
#' @references Hoefling, H., A Path Algorithm for the Fused Lasso Signal Approximator,
#' Journal of Computational and Graphical Statistics (2010)
#' \doi{10.1198/jcgs.2010.09208}
#' @export
flsa_graph <- function(gamma, graph, lambda) {
  flsa_fit <- flsa::flsa(gamma, connListObj = graph2connlist(graph),
                         lambda2 = lambda)
  model_dim <- apply(flsa_fit, 1, function(a) length(unique(a))) %>% unname
  nll <- 1 / 2 * apply(flsa_fit - t(replicate(length(lambda), gamma)),
                       1,
                       function(a) sum(a ^ 2))
  list("result" = flsa_fit, "aic" = 2 * nll + 2 * model_dim,
       "bic" = 2 * nll + log(length(gamma)) * model_dim,
       "gcv" = 2 * nll / ((length(gamma) - model_dim) ^ 2),
       "model_dim" = model_dim, "nll" = nll)
}

