cglasso <- function(X, lo, up, weights, nrho = 50, rho.min.ratio, rho, maxR2, maxit_em = 1.0E+3, thr_em = 1.0e-4, maxit_bcd = 10000, thr_bcd = 1.0e-4, trace = 0L){
    this.call <- this.call <- match.call()
    # checking 'X'
    if(class(X) != "datacggm") X <- datacggm(X, lo = lo, up = up)
    if(X$startmis == 0) {
        warning("the dataset is full observed; 'glasso' is called to fit the l1-penalized Gaussian graphical model")
        out <- glasso(X$X, weights = weights, pendiag = FALSE, nrho = nrho, rho.min.ratio, rho, maxR2, maxit = maxit_bcd, thr = thr_bcd)
        return(out)
    }
    n <- dim(X$X)[1]
    p <- dim(X$X)[2]
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
    # checking 'maxi_emt'
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
    # computing marginal estimates
    par.ini <- parini(X)
    xm <- par.ini$xm
    vm <- par.ini$vm
    # calling cglasso.fit
    out <- cglasso.fit(object = X, w = weights, xm = xm, vm = vm,
                        nrho = nrho, rhoratio = rho.min.ratio, rho = rho, maxR2 = maxR2,
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
                    "1" = "Fit_MarginalDistributions",
                    "2" = "update",
                    "3" = "glasso",
                    "4" = "cglasso")
        if(nrho == 0) stop("cglasso does not congerce: ", sQuote(out$conv), "\n  Error in subroutine ", sQuote(out$subrout))
        msg <- paste("cglasso does not congerce at nrho =", nrho + 1, " with error code", sQuote(out$conv), "\n  Error in subroutine ", sQuote(out$subrout))
        msg <- paste(msg, ifelse(nrho == 1, "\n  The first solution is reported", paste("\n  The first", nrho, "solutions are reported")))
        warning(msg)
        out$rho <- out$rho[1:nrho]
        out$Xipt <- out$Xipt[, , 1:nrho, drop = FALSE]
        out$S <- out$S[, , 1:nrho, drop = FALSE]
        out$mu <- out$mu[, 1:nrho, drop = FALSE]
        out$Sgm <- out$Sgm[, , 1:nrho, drop = FALSE]
        out$Tht <- out$Tht[, , 1:nrho, drop = FALSE]
        out$Adj <- out$Adj[, , 1:nrho, drop = FALSE]
        out$df <- out$df[1:nrho]
        out$R2 <- out$R2[1:nrho]
        out$ncomp <- out$ncomp[1:nrho]
        out$Ck <- out$Ck[, 1:nrho, drop = FALSE]
        out$pk <- out$pk[, 1:nrho, drop = FALSE]
        out$nit <- out$nit[1:nrho, , drop = FALSE]
    }
    # checking if R2 >= maxR2
    if(out$nrho != nrho) {
        nrho <- out$nrho
        out$rho <- out$rho[1:nrho]
        out$Xipt <- out$Xipt[, , 1:nrho, drop = FALSE]
        out$S <- out$S[, , 1:nrho, drop = FALSE]
        out$mu <- out$mu[, 1:nrho, drop = FALSE]
        out$Sgm <- out$Sgm[, , 1:nrho, drop = FALSE]
        out$Tht <- out$Tht[, , 1:nrho, drop = FALSE]
        out$Adj <- out$Adj[, , 1:nrho, drop = FALSE]
        out$df <- out$df[1:nrho]
        out$R2 <- out$R2[1:nrho]
        out$ncomp <- out$ncomp[1:nrho]
        out$Ck <- out$Ck[, 1:nrho, drop = FALSE]
        out$pk <- out$pk[, 1:nrho, drop = FALSE]
        out$nit <- out$nit[1:nrho, , drop = FALSE]
    }
    # preparing object with S3 class 'cglasso'
    out <- c(call = this.call, out)
    class(out) <- c("cglasso", "glasso")
    out
}
