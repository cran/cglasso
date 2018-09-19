mle.glasso <- function(object, ..., maxit = 1.0e+04, thr = 1.0e-04, trace = 0L) {
    this.call <- this.call <- match.call()
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
    # extracting elements from 'object'
    S <- object$S
    p <- dim(S)[1]
    nrho <- object$nrho
    Sgm <- object$Sgm
    Tht <- object$Tht
    R2 <- object$R2
    # storage.mode setting
    storage.mode(p) <- "integer"
    storage.mode(S) <- "double"
    storage.mode(nrho) <- "integer"
    storage.mode(maxit) <- "integer"
    storage.mode(thr) <- "double"
    storage.mode(Sgm) <- "double"
    storage.mode(Tht) <- "double"
    storage.mode(R2) <- "double"
    nit <- integer(nrho)
    conv <- integer(1)
    storage.mode(trace) <- "integer"
    out <- .Fortran(C_mleglasso, p = p, S = S, nrho = nrho, maxit = maxit,
                    thr = thr, Sgm = Sgm, Tht = Tht, R2 = R2, nit = nit,
                    conv = conv, trace = trace)
    # removing unused elements
    out <- out[-1]
    # adding elements form 'object'
    out$X <- object$X
    out$Adj <- object$Adj
    out$rho <- object$rho
    out$Ck <- object$Ck
    out$pk <- object$pk
    out$ncomp <- object$ncomp
    out$df <- object$df
    # checking convergence
    if(out$conv != 0) {
        nrho <- out$nrho
        if(nrho == 0) stop("mle.glasso does not congerce with error code ", sQuote(out$conv))
        msg <- paste("mle.glasso does not congerce at nrho =", nrho + 1, " with error code", sQuote(out$conv))
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
    out <- c(call = this.call, out)
    class(out) <- c("ggm", "glasso")
    out
}
