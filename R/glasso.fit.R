glasso.fit <- function(S, w, pendiag, nrho, rhoratio, rho, maxR2, maxit, thr, trace){
    vnames <- colnames(S)
    p <- dim(S)[1]
    Sgm <- array(0, dim = c(p, p, nrho), dimnames = list(vnames, vnames, NULL))
    Tht <- array(0, dim = c(p, p, nrho), dimnames = list(vnames, vnames, NULL))
    Adj <- array(0, dim = c(p, p, nrho), dimnames = list(vnames, vnames, NULL))
    Ck <- matrix(0, p, nrho)
    pk <- matrix(0, p, nrho)
    # storage.mode setting
    storage.mode(p) <- "integer"
    storage.mode(S) <- "double"
    storage.mode(w) <- "double"
    storage.mode(pendiag) <- "integer"
    storage.mode(nrho) <- "integer"
    storage.mode(rhoratio) <- "double"
    storage.mode(rho) <- "double"
    storage.mode(maxR2) <- "double"
    storage.mode(maxit) <- "integer"
    storage.mode(thr) <- "double"
    storage.mode(Sgm) <- "double"
    storage.mode(Tht) <- "double"
    storage.mode(Adj) <- "integer"
    df <- integer(nrho)
    R2 <- double(nrho)
    ncomp <- integer(nrho)
    storage.mode(Ck) <- "integer"
    storage.mode(pk) <- "integer"
    nit <- integer(nrho)
    conv <- integer(1)
    storage.mode(trace) <- "integer"
    out <- .Fortran(C_glasso, p = p, S = S, w = w, pendiag = pendiag, nrho = nrho,
                    rhoratio = rhoratio, rho = rho, maxR2 = maxR2, maxit = maxit,
                    thr = thr, Sgm = Sgm, Tht = Tht, Adj = Adj, df = df, R2 = R2,
                    ncomp = ncomp, Ck = Ck, pk = pk, nit = nit, conv = conv,
                    trace = trace)
    out
}
