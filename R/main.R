# Suppress R CMD check note
# This just imports a function from the package sf, which is needed.
#' @importFrom sf st_centroid
NULL


#' Segmentation using graph structure
#'
#' @param gamma entry vector to regularize
#' @param graph an igraph object (from package \code{igraph}) giving the regularization structure
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
#' @description
#' These functions provide a clustering of a signal on graph into a piecewise constant
#' signal on graph. Given a graph and a signal `gamma` assigning a value to each node,
#' it returns another signal which is constant over subgraphs where `gamma` has
#' close to equal value. See references.
#'
#' Only parameters `gamma` and `graph` need be provided. The other parameters concern
#' the internals of the estimating procedure and usually do not need to be changed.
#' `agraph` is the general-purpose function. `agraph_prec` does the same thing as
#' `agraph` in the case where `gamma` as a covariance structure. It is provided as
#' the precision matrix `prec`, which has to be a sparse matrix (`Matrix::sparseMatrix`)
#' for fast computation. See Goepp and van de Kassteele (2021).
#'
#' @param gamma input vector to regularize
#' @param graph an igraph object (from package \code{igraph}) giving the regularization structure
#' @param lambda regularizing constant
#' @param shrinkage Boolean, defaults TRUE. Whether to return the adaptive ridge estimate as output. If FALSE, the adaptive ridge is used to define a segmentation into zones, and the signal is estimated on each zone using non-penalized estimation.
#' @param weights weights for gamma. Default value is one.
#' @param delta Computational constant in the adaptive ridge reweighting formula.
#' @param tol Tolerance to test for convergence of the adaptive ridge
#' @param thresh Thresholding constant used to fuse two adjacent regions with close value of \code{gamma}.
#' @param itermax Total number of iterations. Default value is 10000. Setting a low value can make the procedure return NULL entries for some values of \code{lambda}.
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{result}: matrix whose rows are the segmented output of input signal \code{gamma}, for each value of \code{lambda}}
#' \item{\code{bic}, \code{gcv}, and \code{aic}: vectors of length \code{length(lambda)}, giving the BIC, GCV, and AIC criteria for each value of lambda. See references below.}
#' \item{\code{model_dim}, \code{nll}: vectors of length \code{length(lambda)}, giving the model dimension and negative log-likelihood for each value of lambda. See reference below for the definition of these terms.}
#' }
#' @references
#' Schwarz G. (1978)
#' Estimating the Dimension of a Model.
#' Ann. Statist. 6 (2) 461 - 464, March, 1978.
#' \doi{10.1214/aos/1176344136}
#'
#' Akaike H. (1974)
#' A new look at the statistical model identification,
#' in IEEE Transactions on Automatic Control, vol. 19, no. 6, pp. 716-723, December 1974
#' \doi{10.1109/TAC.1974.1100705}
#'
#' Hastie T., Friedman J., and Tibshirani R. (2009)
#' The elements of statistical learning: data mining, inference, and prediction (Vol. 2, pp. 1-758).
#' New York: Springer
#' \doi{10.1007/978-0-387-21606-5}
#'
#' Goepp V. and van de Kassteele J. (2021)
#' Graph-Based Spatial Segmentation of Health-Related Areal Data,
#' arxiv preprint.
#' \doi{10.48550/arXiv.2206.06752}
#' @importFrom Matrix Cholesky solve update Diagonal
#' @export
#' @seealso [graphseg::flsa_graph()]
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
    # lapply(function(a) a - 1) %>%
    lapply(function(a) {
      attributes(a) <- NULL
      return(a)
    }) %>%
    lapply(as.integer) %>%
    'class<-'("connListObj")
  return(connlist)
}
#' Segmentation using graph structure and the fused lasso estimate
#'
#' @description Wrapper around the function \code{flsa::flsa}, which computes the
#' fused lasso signal approximator (see reference). Like `agraph`, this function
#' takes a signal on graph and returns a clustering thereof into a piecewise-constant
#' signal. The difference with `agraph` is the estimation method: `agraph` works well when the
#' true signal is sparse and its computation time scales well to large graphs.
#'
#' @param gamma entry vector to regularize
#' @param graph graph (an \link[igraph]{igraph} object) giving the regularization structure
#' @param lambda regularizing constant
#' @importFrom Matrix Cholesky solve update Diagonal
#' @seealso [graphseg::agraph()]
#' @references
#' Hoefling, H., A Path Algorithm for the Fused Lasso Signal Approximator,
#' Journal of Computational and Graphical Statistics (2010)
#' \doi{10.1198/jcgs.2010.09208}
#' @export
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{result}: matrix whose rows are the segmented output of input signal \code{gamma}, for each value of \code{lambda}}
#' \item{\code{bic}, \code{gcv}, and \code{aic}: vectors of length \code{length(lambda)}, giving the BIC, GCV, and AIC criteria for each value of lambda. See references below.}
#' \item{\code{model_dim}, \code{nll}: vectors of length \code{length(lambda)}, giving the model dimension and negative log-likelihood for each value of lambda. See reference below for the definition of these terms.}
#' }
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

