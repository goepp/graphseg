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
#' @import graphics
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
  pen_mat <- car::wcrossprod(par, par, K)[1]
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
    2 * pen * diag(Matrix::rowSums(K))
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
#' @export
solver <- function(gamma, sigma_sq, K, pen) {
  limSolve::Solve(diag(length(gamma)) + pen * sweep(K, 1, 2 * sigma_sq, `*`),
        gamma)
}
solver_nr <- function(gamma, sigma_sq, K,
                      pen, thresh = 1e-8, maxiter = 100) {
  par <- rep(1, length(gamma))
  iter <- 1
  while (iter < maxiter) {
    old_par <- par
    if (any(is(K) %in% "sparseMatrix")) {
      par <- old_par - solve_ssdp(hessian(old_par, gamma, sigma_sq, K, pen),
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
    theta <- solver(gamma, sigma_sq, K, pen = pen[ind])
    weight_adj <- adj *
      ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2) ^ (-1)
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
#' @export
solve_ssdp <- function(mat, vect) {
  R <- chol(mat)
  temp <- backsolve(R, vect, transpose = TRUE)
  result <- backsolve(R, temp, transpose = FALSE)
  result
}
#'
#' @export
build_jt <- function(g, verbose = FALSE) {
  # clique and separator sets
  N <- 0
  order <- NULL
  C <- gC <- S <- list()
  active <- rep(TRUE, igraph::vcount(g))
  # main loop
  while (sum(active) > 0) {
    # compute neighborhood
    vois <- igraph::neighborhood(g, order = 1)
    gvois <- igraph::graph.neighborhood(g, order = 1)
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
  jt <- igraph::graph(edges, directed = FALSE)
  if (!is.null(igraph::V(g)$name)) {
    igraph::V(jt)$name <- paste0("C",
                                 1:length(C),
                                 "={",
                                 sapply(C, function(z) paste0(igraph::V(g)$name[as.numeric(z)], collapse = ",")),
                                 "}")
  } else {
    igraph::V(jt)$name <- paste0("C",
                                 1:length(C),
                                 "={",
                                 sapply(C, function(z) paste0(as.numeric(z), collapse = ",")),
                                 "}")
  }
  list(C = C, S = S, elim.order = order,
       connect = connect, graph = jt, N = N,
       treewidth = max(unlist(lapply(C, length))))
}
countfillin <- function(g) {
  n <- igraph::vcount(g)
  m <- igraph::ecount(g)
  return(n * (n - 1)/2 - m)
}
