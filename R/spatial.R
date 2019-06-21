#' @export
str_detect_2 <- function(string, pattern_vect) {
  sapply(pattern_vect, str_detect, string = string) %>%
    apply(1, any)
}
#' @export
sel_columns <- function(X, t) {
  Xd <- as.data.frame(X)
  X_new <- c()
  for (ind in 1:max(t)) {
    # colnum <- which(names(Xd) %in% names(t)[which(t == ind)])
    colnum <- str_detect_2(names(Xd), names(t)[which(t == ind)]) %>% which()
    temp <- Xd[colnum]
    X_new <- cbind(X_new, Reduce(`+`, temp))
  }
  X_new %>% unname()
}
#' @export
solver_aridge_lambda <- function(y, X, touches, lambda = 10 ^ -0.5, maxiter = 100, epsilon = 1e-5,
                                 thresh = 1e-8) {
  adjacency <- touches
  theta <- rep(0, ncol(X))
  iter <- 1
  sel <- adjacency
  while (iter < maxiter) {
    old_sel <- sel
    K <- diag(apply(adjacency, 1, sum)) - adjacency
    p <- 1 / (1 + exp(-X %*% theta))
    Gamma <- 1 / p + 1 / (1 - p)
    z <- X %*% theta + Gamma * (y - p)
    temp <- as.vector(exp(-X %*% theta))
    w <- temp / (1 + temp) ^ 2
    reg <- gam(z ~ X, family = gaussian(), data = data.frame(z = z, X = X),
               drop.intercept = TRUE, weights = w, H = lambda * K)
    theta <- reg$coefficients
    gamma <- theta
    adjacency <- touches *
      ((replicate(length(gamma), gamma) - t(replicate(length(gamma), gamma))) ^ 2 + epsilon ^ 2) ^ (-1)
    sel <- adjacency *
      ((replicate(length(gamma), gamma) - t(replicate(length(gamma), gamma))) ^ 2 + epsilon ^ 2)
    converge <- max(abs(old_sel - sel)) <= thresh
    converge
    if (converge) {
      components <- graph_from_adjacency_matrix((sel < 0.99) * touches)
      t <- clusters(components)$membership
      X_sel <- sel_columns(X, t)
      fit <- gam(y ~ X_sel, family = "binomial", drop.intercept = TRUE)
      # temp <- function(ind) rep(fit$coefficients[ind], sum(t == ind))
      temp2 <- function(ind) fit$coefficients[t[ind]]
      par <- sapply(1:length(t), temp2) %>%
        setNames(colnames(X))
      break
    }
    iter <- iter + 1
  }
  par
}
#' @export
solver_aridge <- function(
  y, X, adj, lambda = 10 ^ seq(-2, 2, length = 50),
  maxiter = 100, epsilon = 1e-5, thresh = 1e-8
) {
  par_ls <- vector("list", length = length(lambda))
  weight_adj <- adj
  theta <- rep(0, ncol(X))
  iter <- 1
  ind <- 1
  sel <- adjacency
  while (iter < maxiter) {
    old_sel <- sel
    K <- diag(apply(weight_adj, 1, sum)) - adja
    p <- 1 / (1 + exp(-X %*% theta))
    Gamma <- 1 / p + 1 / (1 - p)
    z <- X %*% theta + Gamma * (y - p)
    temp <- as.vector(exp(-X %*% theta))
    w <- temp / (1 + temp) ^ 2
    reg <- bam(z ~ X, family = gaussian(), data = data.frame(z = z, X = X),
               drop.intercept = TRUE, weights = w, H = lambda * K)
    theta <- reg$coefficients
    weight_adj <- adj *
      ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2) ^ (-1)
    sel <- weight_adj *
      ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2)
    converge <- max(abs(old_sel - sel)) <= thresh
    if (converge) {
      components <- graph_from_adjacency_matrix((sel < 0.99) * adj)
      t <- clusters(components)$membership
      X_sel <- sel_columns(X, t)
      fit <- bam(y ~ X_sel, family = "binomial", drop.intercept = TRUE)
      par_ls[[ind]] <- sapply(seq_along(t), function(ind) fit$coefficients[t[ind]]) %>%
        setNames(colnames(X))
      ind <-  ind + 1
      break
    }
    iter <- iter + 1
    if (ind > length(par_ls)) break
  }
  par
}
#' @export
solver_graph_aridge_old <- function(gamma, sigma_sq, adj,
                                    lambda = 10 ^ seq(-4, 4, length = 100),
                                    maxiter = 1000,
                                    epsilon = 1e-5, thresh = 1e-8) {
  X <-  diag(length(gamma)) %>%
    # Matrix(sparse = TRUE) %>%
    "colnames<-"(colnames(adj))
  weight_adj <- adj
  theta <- gamma * 0
  par_ls <- vector("list", length(lambda))
  iter <- 1
  ind <- 1
  sel <- weight_adj
  while (iter < maxiter) {
    old_sel <- sel
    K <- diag(apply(weight_adj, 1, sum)) - weight_adj
    reg <- gam(gamma ~ X - 1, family = gaussian(),
               data = data.frame(gamma = gamma, X = X),
               drop.intercept = TRUE, weights = 1 / sigma_sq, H = (lambda[ind] * K) %>% as.matrix())
    theta <- reg$coefficients
    weight_adj <- adj *
      ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2 + epsilon ^ 2) ^ (-1)
    sel <- weight_adj * ((replicate(length(theta), theta) - t(replicate(length(theta), theta))) ^ 2)
    converge <- max(abs(old_sel - sel)) <= thresh
    converge
    if (converge) {
      components <- graph_from_adjacency_matrix((sel < 0.99) * adj, mode = "undirected")
      t <- clusters(components)$membership
      X_sel <- sel_columns(X, t)
      fit <- gam(gamma ~ X_sel, family = gaussian(), drop.intercept = TRUE)
      temp <- function(ind) fit$coefficients[t[ind]]
      par_ls[[ind]] <- sapply(1:length(t), temp) %>%
        setNames(colnames(X))
      # weight_adj <- adj
      # sel <- weight_adj
      ind <- ind + 1
    }
    iter <- iter + 1
    if (ind > length(par_ls)) break
  }
  par_ls
}
#' @export
loglik <- function(par, gamma, sigma_sq, K, pen) {
  pen_mat <- wcrossprod(par, par, K)[1]
  sum((par - gamma) ^ 2 / (2 * sigma_sq)) + pen * pen_mat
}
#' @export
loglik_old <- function(par, gamma, sigma_sq, K, pen) {
  pen_mat <- (t(replicate(length(par), par)) - replicate(length(par), par)) ^ 2
  sum((par - gamma) ^ 2 / (2 * sigma_sq)) - pen / 2 * sum(pen_mat * K)
}
#' @export
hessian <- function(par, gamma, sigma_sq, K, pen) {
  diag(1 / sigma_sq) +
    2 * pen * K -
    2 * pen * diag(Matrix::rowSums(K))
}
#' @export
score <- function(par, gamma, sigma_sq, K, pen) {
  pen_matrix <- (replicate(length(par), par) - t(replicate(length(par), par))) *
    K
  as.vector(
    (par - gamma) / sigma_sq - 2 * pen * Matrix::rowSums(pen_matrix)
  )
}
#' @export
solver <- function(gamma, sigma_sq, K, pen) {
  Solve(diag(length(gamma)) + pen * sweep(K, 1, 2 * sigma_sq, `*`),
        gamma)
}
#' @export
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
      par <- old_par - Solve(hessian(old_par, gamma, sigma_sq, K, pen),
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
  X <-  diag(length(gamma)) %>%
    "colnames<-"(colnames(adj))
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
      fit <- gam(gamma ~ X_sel, family = gaussian(), drop.intercept = TRUE)
      temp <- function(ind) fit$coefficients[t[ind]]
      par_ls[[ind]] <- sapply(1:length(t), temp) %>%
        setNames(colnames(X))
      bic[ind] <- log(length(gamma)) * ncol(X_sel) +
        2 * loglik(par_ls[[ind]], gamma, sigma_sq, K, pen = 0)
      bic_1[ind] <- log(length(gamma)) * ncol(X_sel)
      bic_2[ind] <- 2 * loglik(par_ls[[ind]], gamma, sigma_sq, K, pen = 0)
      laplace_bic[ind] <- bic - ncol(X_sel) * log(2 * pi) - sum(log(sigma_sq))
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
    V(jt)$name <- paste0("C", 1:length(C), "={", sapply(C, function(z) paste0(V(g)$name[as.numeric(z)], collapse = ",")),
                         "}")
  } else {
    V(jt)$name <- paste0("C", 1:length(C), "={", sapply(C, function(z) paste0(as.numeric(z), collapse = ",")), "}")
  }
  list(C = C, S = S, elim.order = order,
       connect = connect, graph = jt, N = N,
       treewidth = max(unlist(lapply(C, length))))
}
countfillin <- function(g) {
  n <- vcount(g)
  m <- ecount(g)
  return(n * (n - 1)/2 - m)
}
