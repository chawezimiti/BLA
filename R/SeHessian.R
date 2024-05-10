#' Hessian matrix
#'
#' This is a set of functions that support the censored bivariate normal model.
#'
#' @param hessian If `True`, hessian matrix is used.
#' @param silent Condition of matrix.
#' @param a The hessian matrix.
#' @keywords internal
#'
seHessian<-function(a, hessian = FALSE, silent = FALSE){
  namesp <- colnames(a)
  mathessian <- a
  mathessian <- ifelse(mathessian == -Inf, -1e+09, mathessian)
  mathessian <- ifelse(mathessian == +Inf, 1e+09, mathessian)
  sigma <- try(solve(mathessian), silent = TRUE)
  if (inherits(sigma, "try-error")) {
    if (!silent)
      warning("Error in Hessian matrix inversion")
    mathessianx <- try(as.matrix(getFromNamespace("nearPD",
                                                  ns = "Matrix")(mathessian)$mat), silent = TRUE)
    if (inherits(mathessianx, "try-error")) {
      if (!silent)
        warning("Error in estimation of the Nearest Positive Definite Matrix. Calculates the Moore-Penrose generalized inverse. Use result with caution.")
      sigma <- try(ginv(mathessian), silent = TRUE)
      if (is.null(colnames(sigma)) | is.null(rownames(sigma))) {
        colnames(sigma) <- rownames(sigma) <- colnames(mathessian)
      }
    }
    else {
      if (!silent)
        warning("Calculates the Nearest Positive Definite Matrix. Use result with caution.")
      sigma <- try(solve(mathessianx), silent = TRUE)
    }
  }
  if (!inherits(sigma, "try-error")) {
    if (all(diag(sigma) >= 0)) {
      res <- sqrt(diag(sigma))
    }
    else {
      s. <- svd(sigma)
      R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
      res <- structure((matrix(rep(1, nrow(R)), nrow = 1,
                               byrow = TRUE) %*% R)[1, ], .Names = colnames(mathessian))
      if (any(res < 0)) {
        d <- diag(as.matrix(getFromNamespace("nearPD",
                                             ns = "Matrix")(sigma)$mat))
        names(d) <- colnames(mathessian)
        res <- ifelse(d < 0, NA, sqrt(d))
      }
      if (any(is.na(res))) {
        a <- sigma
        n = dim(a)[1]
        root = matrix(0, n, n)
        for (i in 1:n) {
          sum = 0
          if (i > 1) {
            sum = sum(root[i, 1:(i - 1)]^2)
          }
          x = a[i, i] - sum
          if (x < 0) {
            x = 0
          }
          root[i, i] = sqrt(x)
          if (i < n) {
            for (j in (i + 1):n) {
              if (root[i, i] == 0) {
                x = 0
              }
              else {
                sum = 0
                if (i > 1) {
                  sum = root[i, 1:(i - 1)] %*% t(t(root[j,
                                                        1:(i - 1)]))
                }
                x = (a[i, j] - sum)/root[i, i]
              }
              root[j, i] = x
            }
          }
        }
        colnames(root) <- rownames(root) <- colnames(mathessian)
        pseudoV <- root %*% t(root)
        d <- diag(pseudoV)
        if (any(d != 0) & all(d >= 0)) {
          res <- sqrt(d)
          if (!silent)
            warning("Estimates using pseudo-variance based on Gill & King (2004)")
        }
        else {
          if (!silent)
            warning("Approximation of Cholesky matrix based on Rebonato and Jackel (2000)")
          if (!silent)
            warning("Estimates using pseudo-variance based on Gill & King (2004)")
          newMat <- a
          cholError <- TRUE
          iter <- 0
          while (cholError) {
            iter <- iter + 1
            newEig <- eigen(newMat)
            newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
            newMat <- newEig$vectors %*% diag(newEig2) %*%
              t(newEig$vectors)
            newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
            cholStatus <- try(u <- chol(newMat), silent = TRUE)
            cholError <- ifelse(inherits(cholStatus,
                                         "try-error"), TRUE, FALSE)
          }
          root <- cholStatus
          colnames(root) <- rownames(root) <- colnames(mathessian)
          pseudoV <- root %*% t(root)
          res <- sqrt(diag(pseudoV))
        }
      }
    }
  }
  SEInf <- namesp[!namesp %in% names(res)]
  res <- c(res, structure(rep(+Inf, length(SEInf)), .Names = SEInf))
  if (hessian) {
    return(list(SE = res, hessian = mathessian))
  }
  else {
    return(res)
  }
}


