#' Build junction tree of a graph
#'
#'
#' @param g graph
#' @param verbose whether to print progress while running. Defaults to \code{FALSE}.
#'
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
#'
#' @param g graph
countfillin <- function(g) {
  n <- vcount(g)
  m <- ecount(g)
  return(n * (n - 1)/2 - m)
}