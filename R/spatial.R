#' Penalized log-likelihood, score, and hessian of the graph-based signal.
#'
#'
#'
#' @param par current parameter.
#' @param gamma observed signal on a graph.
#' @param sigma_sq variance of \code{gamma} (functions as weights for \code{gamma}).
#' @param K weighted laplacian matrix of the graph.
#' @param pen penalty constant (must be a positive number). Higher values of \code{pen} yields a more regularized fit.
#' @return The penalized negative log-likelihood of the parameter (function \code{loglik}), or its first (function \code{score}) and second (function \code{hessian}) derivatives.
#' @author Vivien Goepp \email{vivien.goepp@@gmail.com}
#' @rdname loglik
#' @title this is a title
#' @export
#' @importFrom stats setNames
#' @import graphics
#' @import sf
#' @importFrom igraph decompose graph_from_adjacency_matrix clusters
#' @importFrom magrittr %>%
loglik <- function(par, gamma, sigma_sq, K, pen) {
  if (length(pen) != 1 | !is.numeric(pen)) {
    stop("Error: pen must be a numeric vector of length 1")
  }
  par <- as.vector(par)
  gamma <- as.vector(gamma)
  sigma_sq <- as.vector(sigma_sq)
  if (length(par) != length(gamma) | length(par) != length(sigma_sq)) {
    stop("Error: length of par, gamma, and sigma_sq must agree")
  }
  if (is.matrix(K)) {
    pen_mat <- car::wcrossprod(par, par, K)[1]
  } else {# when K is a sparse Matrix
    pen_mat <- car::wcrossprod(par, par, as.matrix(K))[1]
  }
  sum((par - gamma) ^ 2 / (2 * sigma_sq)) + pen * pen_mat
}
#' @rdname loglik
#' @export
hessian <- function(par, gamma, sigma_sq, K, pen) {
  if (length(pen) != 1 | !is.numeric(pen)) {
    stop("Error: pen must be a numeric vector of length 1")
  }
  par <- as.vector(par)
  gamma <- as.vector(gamma)
  sigma_sq <- as.vector(sigma_sq)
  if (length(par) != length(gamma) | length(par) != length(sigma_sq)) {
    stop("Error: length of par, gamma, and sigma_sq must agree")
  }
  diag(1 / sigma_sq) +
    2 * pen * K -
    2 * pen * diag(rowSums(K))
}
#' @rdname loglik
#' @export
score <- function(par, gamma, sigma_sq, K, pen) {
  if (length(pen) != 1 | !is.numeric(pen)) {
    stop("Error: pen must be a numeric vector of length 1")
  }
  par <- as.vector(par)
  gamma <- as.vector(gamma)
  sigma_sq <- as.vector(sigma_sq)
  if (length(par) != length(gamma) | length(par) != length(sigma_sq)) {
    stop("Error: length of par, gamma, and sigma_sq must agree")
  }
  pen_matrix <- (replicate(length(par), par) - t(replicate(length(par), par))) *
    K
  as.vector(
    (par - gamma) / sigma_sq - 2 * pen * Matrix::rowSums(pen_matrix)
  )
}
#' Compute the weighted ridge penalized estimate
#'
#' @param gamma observed signal on a graph
#' @param sigma_sq variance of \code{gamma} (functions as weights for \code{gamma})
#' @param K weighted laplacian matrix of the graph
#' @param pen penalty constant (must be a positive number). Higher values of \code{pen} yields a more regularized fit
#' @param ginv Boolean, defaults FALSE. Whether to use the generalized inverse solver \code{ginv} from
#' package \code{limSolve}
#' @seealso \code{\link{loglik}}, \code{\link{score}}, \code{\link{hessian}}
#'
#' @importFrom limSolve Solve
#' @export
solver <- function(gamma, sigma_sq, K, pen, ginv = FALSE) {
  if (ginv) {
    Solve(diag(length(gamma)) + pen * sweep(K, 1, 2 * sigma_sq, `*`),
          gamma)
  } else {
    solve(diag(length(gamma)) + pen * sweep(K, 1, 2 * sigma_sq, `*`),
          gamma)
  }
}
solver_nr <- function(gamma, sigma_sq, K,
                      pen, thresh = 1e-8, maxiter = 100) {
  par <- rep(1, length(gamma))
  iter <- 1
  while (iter < maxiter) {
    old_par <- par
    if (any(is(K) %in% "sparseMatrix")) {
      par <- old_par - solve_spd(hessian(old_par, gamma, sigma_sq, K, pen),
                                 score(old_par, gamma, sigma_sq, K, pen))
    } else {
      par <- old_par - limSolve::Solve(hessian(old_par, gamma, sigma_sq, K, pen),
                                       score(old_par, gamma, sigma_sq, K, pen))
    }
    converge <-  max(abs(old_par - par)) < thresh
    if (converge) break
    iter <- iter + 1
  }
  if (iter == maxiter) warning("Warning: Newton-Raphson did not converge")
  par
}
#' Regularization of spatial data using adjacency graph
#'
#' @param gamma vector of estimated values
#' @param sigma_sq vector of estimated variances
#' @param adj adjacency matrix of the regions corresponding to each values
#' @param pen penalty value parameter
#' @param maxiter maximal number of iterations in the adaptive ridge algorithm, see details.
#' @param epsilon numerical constant used in the adaptive ridge algorithm. Should be small compared to \code{gamma} and large compared to machine precision.
#' @param thresh relative tolerance for testing the convergence of the adaptive ridge iterations.
#' @seealso agraph
#' @export
graph_aridge <- function(gamma, sigma_sq, adj,
                         pen = 10 ^ seq(-4, 4, length = 100),
                         maxiter = 10000,
                         epsilon = 1e-5,
                         thresh = 1e-8) {
  .Deprecated("agraph")
  par_ls <- vector("list", length(pen))
  bic <- laplace_bic <- pen * 0
  X <- "colnames<-"(diag(length(gamma)), colnames(adj))
  weight_adj <- adj
  theta <- gamma * 0
  iter <- 1
  ind <- 1
  sel <- weight_adj
  while (iter < maxiter) {
    old_sel <- sel
    K <- diag(apply(weight_adj, 1, sum)) - weight_adj
    theta <- solver(gamma, sigma_sq, K, pen = pen[ind], ginv = FALSE)[, 1]
    weight_adj <- adj /
      ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2)
    sel <- weight_adj * ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2)
    converge <- max(abs(old_sel - sel)) <= thresh
    converge
    if (converge) {
      components <- igraph::graph_from_adjacency_matrix((sel < 0.99) * adj, mode = "undirected")
      t <- igraph::clusters(components)$membership
      X_sel <- sel_columns(X, t)
      fit <- mgcv::gam(gamma ~ X_sel, family = stats::gaussian(), drop.intercept = TRUE)
      temp <- function(ind) fit$coefficients[t[ind]]
      par_ls[[ind]] <- setNames(sapply(1:length(t), temp), colnames(X))
      bic[ind] <- log(length(gamma)) * ncol(X_sel) +
        2 * loglik(par_ls[[ind]], gamma, sigma_sq, K, pen = 0)
      laplace_bic[ind] <- bic[ind] - ncol(X_sel) * log(2 * pi) - sum(log(sigma_sq))
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(par_ls)) break
  }
  list(par = par_ls, bic = bic, laplace_bic = laplace_bic)
}
#' Linear solver when the matrix is symmetric positive definite
#'
#' @param mat the matrix
#' @param vect the vector
#' @return \code{mat}\eqn{^{-1}} \code{vec}
solve_spd <- function(mat, vect) {
  R <- chol(mat)
  temp <- backsolve(R, vect, transpose = TRUE)
  result <- backsolve(R, temp, transpose = FALSE)
  result
}
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
  adjacency <- Matrix::as(igraph::as_adjacency_matrix(graph), "symmetricMatrix")
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
  precision <- Matrix::Diagonal(x = weights)
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
  adj <- as(igraph::as_adjacency_matrix(graph), "symmetricMatrix")
  sel <- adj
  converge <- FALSE
  weighted_laplacian_init <- lambda[ind] * (Matrix::Diagonal(x = Matrix::colSums(adj)) - adj) +
    Matrix::Diagonal(x = weights)
  chol_init <- Matrix::Cholesky(weighted_laplacian_init)
  while (iter < itermax) {
    sel_old <- sel
    weighted_laplacian <- lambda[ind] * (Diagonal(x = Matrix::colSums(adj)) - adj) +
      Matrix::Diagonal(x = weights)
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
      nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% precision %*% (result[ind, ] - gamma)
      model_dim[ind] <- sum(diag(Matrix::solve(weighted_laplacian, precision)))
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
  adj <- Matrix::as(igraph::as_adjacency_matrix(graph), "symmetricMatrix")
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
agraph_prec <- function(gamma, graph, prec,
                        lambda = 10 ^ seq(-4, 4, length.out = 50),
                        weights = NULL, shrinkage = TRUE,
                        delta = 1e-10, tol = 1e-8,
                        thresh = 0.01, itermax = 10000) {
  gamma <- as.vector(gamma)
  lambda <- as.vector(lambda)
  if (!isSymmetric(prec)) {
    precision <- Matrix::forceSymmetric(prec)
    warning("Precision matrix was not symmetric. It was coerced to symmetric.")
  }
  precision <- precision %>% as("dsCMatrix")
  p <- length(gamma)
  prec_gamma <- precision %*% gamma
  nll <- model_dim <- rep(0, length(lambda))
  result <- matrix(NA, length(lambda), p)
  ind <- 1
  iter <- 1
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be in igraph format")
  }
  edgelist_tmp <- igraph::as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adj <- Matrix::as(igraph::as_adjacency_matrix(graph), "symmetricMatrix")
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
      nll[ind] <- 1 / 2 * t(result[ind, ] - gamma) %*% precision %*% (result[ind, ] - gamma)
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

