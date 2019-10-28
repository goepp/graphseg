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
#' @importFrom stringr str_detect
#' @importFrom mgcv bam
#' @importFrom stats setNames
#' @importFrom car wcrossprod
#' @importFrom Matrix rowSums
#' @import graphics
#' @import sf
#' @importFrom igraph decompose graph_from_adjacency_matrix clusters V E spectrum
#' @import magrittr
#' @import methods
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
#' @param thresh relative tolerance for the convergence of the adaptive ridge iterations.
#' @export
graph_aridge <- function(gamma, sigma_sq, adj,
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
#' Linear solver when the matrix is symmetric positive definite
#'
#' @param mat the matrix
#' @param vect the vector
#' @return \code{mat}\eqn{^{-1}} \code{vec}
#' @export
solve_spd <- function(mat, vect) {
  R <- chol(mat)
  temp <- backsolve(R, vect, transpose = TRUE)
  result <- backsolve(R, temp, transpose = FALSE)
  result
}
#' Build junction tree of a graph
#'
#'
#' @param g graph
#' @param verbose whether to print progress while running. Defaults to \code{FALSE}.
#'
#' @importFrom igraph vcount neighborhood graph.neighborhood graph V "V<-"
#' @export
build_jt <- function(g, verbose = FALSE) {
  # clique and separator sets
  N <- 0
  order <- NULL
  C <- gC <- S <- list()
  active <- rep(TRUE, vcount(g))
  # main loop
  while (sum(active) > 0) {
    # compute neighborhood
    vois <- neighborhood(g, order = 1)
    gvois <- graph.neighborhood(g, order = 1)
    # fillin
    fillin <- unlist(lapply(gvois, countfillin))
    # elim min fillin (among active vertices)
    elim <- which(active)[sort(fillin[active], index.return = TRUE)$ix[1]]
    # compute remove set from elim
    remove <- NULL
    for (v in vois[[elim]]) {
      # if vois[[v]] subset vois[[elim]]
      if (sum(is.element(vois[[v]], vois[[elim]])) == length(vois[[v]]))
        remove <- c(remove, v)
    }
    # add fillin
    for (v in vois[[elim]]) g[v, setdiff(vois[[elim]], v)] <- TRUE
    # disconnect and inactive remove set
    order <- c(order, remove)
    g[remove, ] <- FALSE
    active[remove] <- FALSE
    # update clique and separator set
    N <- N + 1
    C[[N]] <- vois[[elim]]
    gC[[N]] <- gvois[[elim]]
    S[[N]] <- setdiff(C[[N]], remove)
    if (verbose)
      print(paste0(N, " cliques, ", sum(active), " nodes left"))
  }
  # connect separators to cliques
  connect <- list()
  edges <- NULL
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      if (sum(is.element(S[[i]], C[[j]])) == length(S[[i]])) {
        connect[[i]] <- j
        edges <- c(edges, c(i, j))
        break
      }
    }
    connect[[N]] <- numeric(0)
  }
  # build the JT
  jt <- graph(edges, directed = FALSE)
  if (!is.null(V(g)$name)) {
    V(jt)$name <- paste0("C",
                                 1:length(C),
                                 "={",
                                 sapply(C, function(z) paste0(V(g)$name[as.numeric(z)], collapse = ",")),
                                 "}")
  } else {
    V(jt)$name <- paste0("C",
                                 1:length(C),
                                 "={",
                                 sapply(C, function(z) paste0(as.numeric(z), collapse = ",")),
                                 "}")
  }
  list(C = C, S = S, elim.order = order,
       connect = connect, graph = jt, N = N,
       treewidth = max(unlist(lapply(C, length))))
}
#' Utility function for \code{buid_jt}
#'
#' @importFrom igraph vcount ecount
#'
#' @param g graph
countfillin <- function(g) {
  n <- vcount(g)
  m <- ecount(g)
  return(n * (n - 1)/2 - m)
}
#' Segmentation using graph structure
#' @param gamma entry vector to regularize
#' @param graph \url{igraph} giving the regularization structure
#' @param lambda regularizing constant
#' @param weights weights for gamma. Default value is one.
#' @param delta Computational constant in the adaptive ridge reweighting formula.
#' @param tol Tolerance to test for convergence of the adaptive ridge
agraph_one_lambda <- function(gamma, graph, lambda = 1e0, weights = NULL,
                   delta = 1e-6, tol = 1e-10,
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
  edgelist_tmp <- as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adjacency <- as(as_adjacency_matrix(graph), "symmetricMatrix")
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
  graph_del <- graph %>% delete_edges(which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
  theta <- ave(as.vector(gamma), components(graph_del)$membership)
  return(theta)
}
#' Segmentation using graph structure
#' @param gamma entry vector to regularize
#' @param graph \url{igraph} giving the regularization structure
#' @param lambda regularizing constant
#' @param weights weights for gamma. Default value is one.
#' @param delta Computational constant in the adaptive ridge reweighting formula.
#' @param tol Tolerance to test for convergence of the adaptive ridge
#' @export
agraph <- function(gamma, graph, lambda = 10 ^ seq(-2, 2, length.out = 100),
                   weights = NULL,
                   delta = 1e-6, tol = 1e-10,
                   thresh = 0.01, itermax = 5000) {
  gamma <- as.vector(gamma)
  lambda <- as.vector(lambda)
  p <- length(gamma)
  if (is.null(weights)) {
    weights <- rep(1, p)
  }
  bic <- rep(0, length(lambda))
  result <- matrix(NA, p, length(lambda))
  ind <- 1
  iter <- 1
  if (!(class(graph) %in% "igraph")) {
    stop("Graph must be in igraph format")
  }
  edgelist_tmp <- as_edgelist(graph, names = FALSE)
  edgelist <- edgelist_tmp[order(edgelist_tmp[, 2]), c(2, 1)]
  adjacency <- as(as_adjacency_matrix(graph), "symmetricMatrix")
  sel <- adj
  converge <- FALSE
  weighted_laplacian_init <- lambda[1] * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
  chol_init <- Cholesky(weighted_laplacian_init)
  while (iter < itermax) {
    sel_old <- sel
    weighted_laplacian <- lambda[ind] * (Diagonal(x = colSums(adj)) - adj) + Diagonal(x = weights)
    chol <- update(chol_init, weighted_laplacian)
    theta <- solve(chol, weights * gamma)
    adjacency@x <- 1 / ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    sel@x <- (theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 /
      ((theta[edgelist[, 1]] - theta[edgelist[, 2]]) ^ 2 + delta)
    converge <- all(abs(sel@x - sel_old@x) < tol)
    if (converge) {
      graph_del <- graph %>% delete_edges(which((sel@x > 1 - thresh)[order(edgelist[, 2])]))
      segmentation <- components(graph_del)$membership
      result[, ind] <- ave(as.vector(gamma), segmentation)
      bic[ind] <- sum(weights * (result[, ind] - gamma) ^ 2) + log(p) * max(segmentation)
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(lambda)) break
  }
  return(list(result = result, bic = bic))
}
