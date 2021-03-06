
library(tidyverse)
library(Matrix)


## help(Cholesky)
data(KNex)
mtm <- with(KNex, crossprod(mm))
str(mtm@factors) # empty list()
(C1 <- Cholesky(mtm))             # uses show(<MatrixFactorization>)
str(mtm@factors) # 'sPDCholesky' (simpl)
(Cm <- Cholesky(mtm, super = TRUE))
c(C1 = isLDL(C1), Cm = isLDL(Cm))
str(mtm@factors) # 'sPDCholesky'  *and* 'SPdCholesky'
str(cm1  <- as(C1, "sparseMatrix"))
str(cmat <- as(Cm, "sparseMatrix"))# hmm: super is *less* sparse here
cm1[1:20, 1:20]

aaa <- expand(C1)
bbb <- expand(Cm)
all.equal(as(t(aaa$P) %*% tcrossprod(cm1) %*% aaa$P, "symmetricMatrix"), mtm)
all.equal(as(t(bbb$P) %*% tcrossprod(cm1) %*% bbb$P, "symmetricMatrix"), mtm)

b <- matrix(c(rep(0, 711), 1), ncol = 1)
## solve(Cm, b) by default solves  Ax = b, where A = Cm'Cm (= mtm)!
## hence, the identical() check *should* work, but fails on some GOTOblas:
x <- solve(Cm, b)
stopifnot(identical(x, solve(Cm, b, system = "A")),
          all.equal(x, solve(mtm, b)))
all.equal(solve(Cm, b), solve(Cm, b, system = "A"))

Cn <- Cholesky(mtm, perm = FALSE)# no permutation -- much worse:
sizes <- c(simple = object.size(C1),
           super  = object.size(Cm),
           noPerm = object.size(Cn))
## simple is 100, super= 137, noPerm= 812 :
noquote(cbind(format(100 * sizes / sizes[1], digits = 4)))


## Visualize the sparseness:
dq <- function(ch) paste('"', ch, '"', sep = "") ## dQuote(<UTF-8>) gives bad plots
image(mtm, main = paste("crossprod(mm) : Sparse", dq(class(mtm))))
image(cm1, main = paste("as(Cholesky(crossprod(mm)),\"sparseMatrix\"):",
                       dq(class(cm1))))


## Smaller example, with same matrix as in  help(chol) :
(mm <- Matrix(toeplitz(c(10, 0, 1, 0, 3)), sparse = TRUE)) # 5 x 5
(opts <- expand.grid(perm = c(TRUE, FALSE), LDL = c(TRUE, FALSE), super = c(FALSE, TRUE)))
rr <- lapply(seq_len(nrow(opts)), function(i)
  do.call(Cholesky, c(list(A = mm), opts[i,])))
nn <- do.call(expand.grid, c(attr(opts, "out.attr")$dimnames,
                             stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE))
names(rr) <- apply(nn, 1, function(r)
  paste(sub("(=.).*","\\1", r), collapse = ","))
str(rr, max = 1)

str(re <- lapply(rr, expand), max = 2)
## each has a 'P' and a 'L' matrix

R0 <- chol(mm, pivot = FALSE)
R1 <- chol(mm, pivot = TRUE )
stopifnot(all.equal(t(R1), re[[1]]$L),
          all.equal(t(R0), re[[2]]$L),
          identical(as(1:5, "pMatrix"), re[[2]]$P), # no pivoting
          TRUE)

# Version of the underlying SuiteSparse library by Tim Davis :
.SuiteSparse_version()


## help(CHMfactor)
showMethods("Cholesky")

