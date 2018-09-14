cglasso.fit <- function(object, w, xm, vm, nrho, rhoratio, rho, maxR2, maxit_em, thr_em, maxit_bcd, thr_bcd, trace){
    X <- object$X
    up <- object$up
    lo <- object$lo
    R <- object$R
    startmis <- object$startmis
    vnames <- colnames(X)
    n <- dim(X)[1]
    p <- dim(X)[2]
    mu <- matrix(0, p, nrho, dimnames = list(vnames, NULL))
    Xipt <- array(0, dim = c(n, p, nrho), dimnames = list(NULL, vnames, NULL))
    S <- array(0, dim = c(p, p, nrho), dimnames = list(vnames, vnames, NULL))
    Sgm <- array(0, dim = c(p, p, nrho), dimnames = list(vnames, vnames, NULL))
    Tht <- array(0, dim = c(p, p, nrho), dimnames = list(vnames, vnames, NULL))
    Adj <- array(0, dim = c(p, p, nrho), dimnames = list(vnames, vnames, NULL))
    Ck <- matrix(0, p, nrho)
    pk <- matrix(0, p, nrho)
    nit <- matrix(0, nrho, 2, dimnames = list(NULL, c("EM", "nit")))
    # storage.mode setting
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(X) <- "double"
    storage.mode(R) <- "integer"
    storage.mode(startmis) <- "integer"
    storage.mode(lo) <- "double"
    storage.mode(up) <- "double"
    storage.mode(w) <- "double"
    storage.mode(xm) <- "double"
    storage.mode(vm) <- "double"
    storage.mode(nrho) <- "integer"
    storage.mode(rhoratio) <- "double"
    storage.mode(rho) <- "double"
    storage.mode(maxR2) <- "double"
    storage.mode(maxit_em) <- "integer"
    storage.mode(thr_em) <- "double"
    storage.mode(maxit_bcd) <- "integer"
    storage.mode(thr_bcd) <- "double"
    storage.mode(Xipt) <- "double"
    storage.mode(S) <- "double"
    storage.mode(mu) <- "double"
    storage.mode(Sgm) <- "double"
    storage.mode(Tht) <- "double"
    storage.mode(Adj) <- "integer"
    df <- integer(nrho)
    R2 <- double(nrho)
    ncomp <- integer(nrho)
    storage.mode(Ck) <- "integer"
    storage.mode(pk) <- "integer"
    storage.mode(nit) <- "integer"
    conv <- integer(1)
    subrout <- integer(1)
    storage.mode(trace) <- "integer"
    out <- .Fortran(C_cglasso, n = n, p = p, X = X, R = R, startmis = startmis,
                    lo = lo, up = up, w = w, xm = xm, vm = vm,
                    nrho = nrho, rhoratio = rhoratio, rho = rho, maxR2 = maxR2,
                    maxit_em = maxit_em, thr_em = thr_em, maxit_bcd = maxit_bcd,
                    thr_bcd = thr_bcd, Xipt = Xipt, S = S, mu = mu, Sgm = Sgm,
                    Tht = Tht, Adj = Adj, df = df, R2 = R2, ncomp = ncomp, Ck = Ck,
                    pk = pk, nit = nit, conv = conv, subrout = subrout, trace = trace)
    out$X <- object
    out$w[out$w == -1] <- 1
    # removing unused elements
    out <- out[-c(1, 2, 4, 5, 6, 7)]
    attributes(out)$names[c(2, 6)] <- c("weights", "rho.min.ratio")
    out
}
