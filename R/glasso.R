glasso <- function(X, weights, pendiag = FALSE, nrho = 50L, rho.min.ratio, rho, maxR2, maxit = 1.0E+4, thr = 1.0e-04, trace = 0L){
    this.call <- this.call <- match.call()
    # checking 'X'
    if(missing(X)) stop(sQuote("X"), " is missing")
    if(!is.matrix(X)) stop(sQuote("X"), " is not a matrix")
    if(is.null(colnames(X))) {
        vnames <- paste("X", 1:p, sep = "")
        colnames(X) <- vnames
    } else vnames <- colnames(X)
    n <- dim(X)[1]
    p <- dim(X)[2]
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
        rho.min.ratio <- rho[1] / rho[nrho]
    }
    # checking 'maxR2'
    if(missing(maxR2)) maxR2 <- ifelse(p < n, 1, 0.9)
    else {
        if(!is.vector(maxR2)) stop(sQuote("maxR2"), " is not a vector of length ", sQuote(1))
        if(length(maxR2) != 1) stop(sQuote("maxR2"), " is not a vector of length ", sQuote(1))
        if(maxR2 < 0 | maxR2 > 1) stop(sQuote("maxR2"), " does not belong to the interval [0, 1]")
    }
    # checking 'maxit'
    if(!is.vector(maxit)) stop(sQuote("maxit"), " is not a vector of length ", sQuote(1))
    if(length(maxit) != 1) stop(sQuote("maxit"), " is not a vector of length ", sQuote(1))
    if(abs(as.integer(maxit) - maxit) > 0) stop(sQuote("maxit"), " is not an object of type ", dQuote("integer"))
    if(maxit <= 0) stop(sQuote("maxit"), " is not a positive integer")
    # checking 'thr'
    if(!is.vector(thr)) stop(sQuote("thr"), " is not a vector of length ", sQuote(1))
    if(length(thr) != 1) stop(sQuote("thr"), " is not a vector of length ", sQuote(1))
    if(thr < 0 ) stop(sQuote("thr"), " is not a vector of positive values")
    # checking 'trace'
    if(!is.vector(trace)) stop(sQuote("trace"), " is not a vector of length ", sQuote(1))
    if(length(trace) != 1) stop(sQuote("trace"), " is not a vector of length ", sQuote(1))
    if(abs(as.integer(trace) - trace) > 0) stop(sQuote("trace"), " is not an object of type ", dQuote("integer"))
    if(!is.element(trace, c(0L, 1L, 2L))) stop(sQuote("trace"), " can be equal to '0', '1' or '2'")
    # calling glasso.fit
    mx <- colMeans(X)
    S <- crossprod(X) / n - outer(mx, mx)
    out <- glasso.fit(S = S, w = weights, pendiag = pendiag, nrho = nrho,
                        rhoratio = rho.min.ratio, rho = rho, maxR2 = maxR2,
                        maxit = maxit, thr = thr, trace = trace)
    # checking convergence
    if(out$conv != 0) {
        nrho <- out$nrho
        if(nrho == 0) stop("glasso does not congerce with error code ", sQuote(out$conv))
        msg <- paste("glasso does not congerce at nrho =", nrho + 1, " with error code", sQuote(out$conv))
        msg <- paste(msg, ifelse(nrho == 1, "\n  The first solution is reported", paste("\n  The first", nrho, "solutions are reported")))
        warning(msg)
        out$rho <- out$rho[1:nrho]
        out$Sgm <- out$Sgm[, , 1:nrho, drop = FALSE]
        out$Tht <- out$Tht[, , 1:nrho, drop = FALSE]
        out$Adj <- out$Adj[, , 1:nrho, drop = FALSE]
        out$df <- out$df[1:nrho]
        out$R2 <- out$R2[1:nrho]
        out$ncomp <- out$ncomp[1:nrho]
        out$Ck <- out$Ck[, 1:nrho, drop = FALSE]
        out$pk <- out$pk[, 1:nrho, drop = FALSE]
        out$nit <- out$nit[1:nrho]
    }
    if(out$nrho != nrho) {
        nrho <- out$nrho
        out$rho <- out$rho[1:nrho]
        out$Sgm <- out$Sgm[, , 1:nrho, drop = FALSE]
        out$Tht <- out$Tht[, , 1:nrho, drop = FALSE]
        out$Adj <- out$Adj[, , 1:nrho, drop = FALSE]
        out$df <- out$df[1:nrho]
        out$R2 <- out$R2[1:nrho]
        out$ncomp <- out$ncomp[1:nrho]
        out$Ck <- out$Ck[, 1:nrho, drop = FALSE]
        out$pk <- out$pk[, 1:nrho, drop = FALSE]
        out$nit <- out$nit[1:nrho]
    }
    # preparing object with S3 class 'glasso'
    out[[1]] <- X
    out <- c(call = this.call, out)
    out$w[out$w == -1] <- 1
    out$pendiag <- as.logical(out$pendiag)
    attributes(out)$names[c(2, 4, 7)] <- c("X", "weights", "rho.min.ratio")
    class(out) <- "glasso"
    out
}
