mglasso <- function(X, weights, pendiag = FALSE, nrho = 50L, rho.min.ratio, rho, maxR2, maxit_em = 1.0E+3, thr_em = 1.0E-4, maxit_bcd = 1.0E+4, thr_bcd = 1.0E-4, trace = 0L){
    this.call <- match.call()
    # checking 'X'
    R <- is.na(X)
    if(!any(R)) {
        warning("the dataset is full observed; 'glasso' is called to fit the l1-penalized Gaussian graphical model")
        out <- glasso(X, weights = weights, pendiag = pendiag, nrho = nrho, rho.min.ratio, rho, maxR2, maxit = maxit_bcd, thr = thr_bcd)
        return(out)
    }
    if(missing(X)) stop(sQuote("X"), " is missing")
    if(!is.matrix(X)) stop(sQuote("X"), " is not a matrix")
    if(any(!is.finite(X[!R]))) stop("some element in 'X' is not finite")
    # for internal purposes, NA are treated as right-censored values
    X[R] <- .Machine$double.xmax
    X <- datacggm(X, up = .Machine$double.xmax)
    X$R <- X$R[-1L, -1L]
    n <- dim(X)[1L]
    p <- dim(X)[2L]
    vnames <- colnames(X$X)
    # checking 'weights'
    if(missing(weights)) weights <- matrix(rep(-1, p * p), nrow = p, ncol = p, dimnames = list(vnames, vnames))
    else {
        if(!is.matrix(weights)) stop(sQuote("weights"), " is not a matrix")
        if(dim(weights)[1] != dim(weights)[2]) stop(sQuote("weights"), " is not a square matrix")
        if(!isSymmetric(weights)) stop(sQuote("weights"), " is not a symmetric matrix")
        if(any(weights[weights != 0] < 0)) stop("negative weights are not allowed")
        weights[weights == +Inf] <- .Machine$double.xmax
    }
    # checking 'pendiag'
    if(!is.logical(pendiag)) stop(sQuote("pendiag"), " is not an object of type ", dQuote("logical"))
    if(!is.vector(pendiag)) stop(sQuote("pendiag"), " is not a vector of length ", sQuote(1))
    if(length(pendiag) != 1) stop(sQuote("pendiag"), " is not a vector of length ", sQuote(1))
    # checking 'nrho'
    if(!is.vector(nrho)) stop(sQuote("nrho"), " is not a vector of length ", sQuote(1))
    if(length(nrho) != 1) stop(sQuote("nrho"), " is not a vector of length ", sQuote(1))
    if(abs(as.integer(nrho) - nrho) > 0) stop(sQuote("nrho"), " is not an object of type ", dQuote("integer"))
    if(nrho <= 0) stop(sQuote("nrho"), " is not a positive integer")
    # checking 'rho.min.ratio'
    if(missing(rho.min.ratio)) rho.min.ratio <- ifelse(p < n, 1.0e-4, 1.0e-2)
    else {
        if(!is.vector(rho.min.ratio)) stop(sQuote("rho.min.ratio"), " is not a vector of length ", sQuote(1))
        if(length(rho.min.ratio) != 1) stop(sQuote("rho.min.ratio"), " is not a vector of length ", sQuote(1))
        if(rho.min.ratio < 0 | rho.min.ratio > 1) stop(sQuote("rho.min.ratio"), " does not belong to the interval [0, 1]")
    }
    # checking 'rho'
    if(missing(rho)) rho <- rep(-1, nrho)
    else {
        if(!is.vector(rho)) stop(sQuote("rho"), " is not a vector")
        if(any(rho < 0)) stop(sQuote("rho"), " is not a vector of positive values")
        rho <- sort(rho, decreasing = TRUE)
        nrho <- length(rho)
        rho.min.ratio <- rho[nrho] / rho[1L]
    }
    # checking 'maxR2'
    if(missing(maxR2)) maxR2 <- ifelse(p < n, 1, 0.9)
    else {
        if(!is.vector(maxR2)) stop(sQuote("maxR2"), " is not a vector of length ", sQuote(1))
        if(length(maxR2) != 1) stop(sQuote("maxR2"), " is not a vector of length ", sQuote(1))
        if(maxR2 < 0 | maxR2 > 1) stop(sQuote("maxR2"), " does not belong to the interval [0, 1]")
    }
    # checking 'maxi_em'
    if(!is.vector(maxit_em)) stop(sQuote("maxit_em"), " is not a vector of length ", sQuote(1))
    if(length(maxit_em) != 1) stop(sQuote("maxit_em"), " is not a vector of length ", sQuote(1))
    if(abs(as.integer(maxit_em) - maxit_em) > 0) stop(sQuote("maxit_em"), " is not an object of type ", dQuote("integer"))
    if(maxit_em <= 0) stop(sQuote("maxit_em"), " is not a positive integer")
    # checking 'thr_em'
    if(!is.vector(thr_em)) stop(sQuote("thr_em"), " is not a vector of length ", sQuote(1))
    if(length(thr_em) != 1) stop(sQuote("thr_em"), " is not a vector of length ", sQuote(1))
    if(thr_em < 0 ) stop(sQuote("thr_em"), " is not a vector of positive values")
    # checking 'maxit_bcd'
    if(!is.vector(maxit_bcd)) stop(sQuote("maxit_bcd"), " is not a vector of length ", sQuote(1))
    if(length(maxit_bcd) != 1) stop(sQuote("maxit_bcd"), " is not a vector of length ", sQuote(1))
    if(abs(as.integer(maxit_bcd) - maxit_bcd) > 0) stop(sQuote("maxit_bcd"), " is not an object of type ", dQuote("integer"))
    if(maxit_bcd <= 0) stop(sQuote("maxit_bcd"), " is not a positive integer")
    # checking 'thr_bcd'
    if(!is.vector(thr_bcd)) stop(sQuote("thr_bcd"), " is not a vector of length ", sQuote(1))
    if(length(thr_bcd) != 1) stop(sQuote("thr_bcd"), " is not a vector of length ", sQuote(1))
    if(thr_bcd < 0 ) stop(sQuote("thr_bcd"), " is not a vector of positive values")
    # checking 'trace'
    if(!is.vector(trace)) stop(sQuote("trace"), " is not a vector of length ", sQuote(1))
    if(length(trace) != 1) stop(sQuote("trace"), " is not a vector of length ", sQuote(1))
    if(abs(as.integer(trace) - trace) > 0) stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
    if(!is.element(trace, c(0L, 1L, 2L))) stop(sQuote("trace"), " can be equal to '0', '1' or '2'")
    # calling mglasso.fit
    out <- mglasso.fit(object = X, w = weights, pendiag = pendiag, nrho = nrho,
                        rhoratio = rho.min.ratio, rho = rho, maxR2 = maxR2,
                        maxit_em = maxit_em, thr_em = thr_em, maxit_bcd = maxit_bcd,
                        thr_bcd = thr_bcd, trace = trace)
    # checking convergence
    if(out$conv != 0) {
        nrho <- out$nrho
        out$conv <- switch(as.character(out$conv),
                    "-1" = "memory allocation error",
                     "1" = "maximum number of iterations has been exceeded",
                     "2" = "error in computing the coditional moments",
                     "3" = "matrix inversion failed")
        out$subrout <- switch(as.character(out$subrout),
                    "0" = "Ok",
                    "1" = "impute",
                    "2" = "glasso",
                    "3" = "mglasso")
        if(nrho == 0) stop("mglasso does not congerce: ", sQuote(out$conv), "\n  Error in subroutine ", sQuote(out$subrout))
        msg <- paste("mglasso does not congerce at nrho =", nrho + 1, " with error code", sQuote(out$conv), "\n  Error in subroutine ", sQuote(out$subrout))
        msg <- paste(msg, ifelse(nrho == 1, "\n  The first solution is reported", paste("\n  The first", nrho, "solutions are reported")))
        warning(msg)
        seq_id <- seq_len(nrho)
        out$rho <- out$rho[1:nrho]
        out$Xipt <- out$Xipt[, , seq_id, drop = FALSE]
        out$S <- out$S[, , seq_id, drop = FALSE]
        out$mu <- out$mu[, seq_id, drop = FALSE]
        out$Sgm <- out$Sgm[, , seq_id, drop = FALSE]
        out$Tht <- out$Tht[, , seq_id, drop = FALSE]
        out$Adj <- out$Adj[, , seq_id, drop = FALSE]
        out$df <- out$df[seq_id]
        out$R2 <- out$R2[seq_id]
        out$ncomp <- out$ncomp[seq_id]
        out$Ck <- out$Ck[, seq_id, drop = FALSE]
        out$pk <- out$pk[, seq_id, drop = FALSE]
        out$nit <- out$nit[seq_id, , drop = FALSE]
    }
    # checking if R2 >= maxR2
    if(out$nrho != nrho) {
        nrho <- out$nrho
        seq_id <- seq_len(nrho)
        out$rho <- out$rho[seq_id]
        out$Xipt <- out$Xipt[, , seq_id, drop = FALSE]
        out$S <- out$S[, , seq_id, drop = FALSE]
        out$mu <- out$mu[, seq_id, drop = FALSE]
        out$Sgm <- out$Sgm[, , seq_id, drop = FALSE]
        out$Tht <- out$Tht[, , seq_id, drop = FALSE]
        out$Adj <- out$Adj[, , seq_id, drop = FALSE]
        out$df <- out$df[seq_id]
        out$R2 <- out$R2[seq_id]
        out$ncomp <- out$ncomp[seq_id]
        out$Ck <- out$Ck[, seq_id, drop = FALSE]
        out$pk <- out$pk[, seq_id, drop = FALSE]
        out$nit <- out$nit[seq_id, , drop = FALSE]
    }
    # preparing object with S3 class 'mglasso'
    out <- c(call = this.call, out)
    class(out) <- c("mglasso", "glasso")
    out
}
