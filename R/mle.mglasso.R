mle.mglasso <- function(object, ..., maxit_em = 1.0e+03, thr_em = 1.0e-4, maxit_bcd = 1.0e+4, thr_bcd = 1.0e-4, trace = 0L) {
    this.call <- this.call <- match.call()
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
    # extracting elements from 'object'
    X <- object$X
    R <- is.na(X)
    X[R] <- .Machine$double.xmax
    X <- datacggm(X, up = .Machine$double.xmax)
    R <- X$R[-1L, -1L]
    startmis <- X$startmis
    X <- X$X
    n <- dim(X)[1L]
    p <- dim(X)[2L]
    nrho <- object$nrho
    Xipt <- object$Xipt
    S <- object$S
    mu <- object$mu
    Sgm <- object$Sgm
    Tht <- object$Tht
    R2 <- object$R2
    nit <- matrix(0L, nrho, 2, dimnames = list(NULL, c("EM", "nit")))
    conv <- object$conv
    subrout <- object$subrout
    # storage.mode setting
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(X) <- "double"
    storage.mode(R) <- "integer"
    storage.mode(startmis) <- "integer"
    storage.mode(nrho) <- "integer"
    storage.mode(maxit_em) <- "integer"
    storage.mode(thr_em) <- "double"
    storage.mode(maxit_bcd) <- "integer"
    storage.mode(thr_bcd) <- "double"
    storage.mode(Xipt) <- "double"
    storage.mode(S) <- "double"
    storage.mode(mu) <- "double"
    storage.mode(Sgm) <- "double"
    storage.mode(Tht) <- "double"
    R2 <- double(nrho)
    storage.mode(nit) <- "integer"
    conv <- integer(1L)
    subrout <- integer(1L)
    storage.mode(trace) <- "integer"
    out <- .Fortran(C_mlemglasso, n = n, p = p, X = X, R = R, startmis = startmis,
                    nrho = nrho, maxit_em = maxit_em, thr_em = thr_em,
                    maxit_bcd = maxit_bcd, thr_bcd = thr_bcd, Xipt = Xipt, S = S, mu = mu,
                    Sgm = Sgm, Tht = Tht, R2 = R2, nit = nit, conv = conv, subrout = subrout,
                    trace = trace)
    # removing unused elements
    out <- out[-c(1L, 2L, 4L, 5L)]
    # adding elements form 'object'
    out$X[R[, -1L] == 1L] <- NA
    out$Adj <- object$Adj
    out$rho <- object$rho
    out$Ck <- object$Ck
    out$pk <- object$pk
    out$ncomp <- object$ncomp
    out$df <- object$df
    # checking convergence
    if(out$conv != 0L) {
        nrho <- out$nrho
        out$subrout <- switch(as.character(out$subrout),
                        "0" = "Ok",
                        "1" = "impute",
                        "2" = "glasso",
                        "3" = "mlemglasso")
        out$conv <- switch(as.character(out$conv),
                        "-1" = "memory allocation error",
                        "1" = "maximum number of iterations has been exceeded",
                        "2" = "error in computing the coditional moments",
                        "3" = "matrix inversion failed")
        if(nrho == 0) stop("mle.mglasso does not congerce: ", sQuote(out$conv), "\n  Error in subroutine ", sQuote(out$subrout))
        msg <- paste("mle.mglasso does not congerce at nrho =", nrho + 1L, " with error code", sQuote(out$conv), "\n  Error in subroutine ", sQuote(out$subrout))
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
        out$nit <- out$nit[1:nrho, , drop = TRUE]
    }
    out <- c(call = this.call, out)
    class(out) <- c("mggm", "mglasso", "glasso")
    out
}
