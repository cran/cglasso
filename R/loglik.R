loglik <- function(object) {
    if(is.element(class(object)[1], c("glasso", "ggm"))) {
        X <- object$X
        fun <- "log-Lik"
    }
    if(is.element(class(object)[1], c("cglasso", "cggm"))) {
        X <- object$X$X
        fun <- "Q-Fun"
    }
    n <- dim(X)[1]
    p <- dim(X)[2]
    nrho <- object$nrho
    S <- object$S
    Tht <- object$Tht
    Ck <- object$Ck
    pk <- object$pk
    ncomp <- object$ncomp
    df <- object$df + 2 * p
    if(length(dim(S)) == 2) {
        dim(S) <- c(p, p, 1)
        ii <- rep.int(1, nrho)
    } else ii <- 1:nrho
    c1 <- 0.5 * n
    c2 <- log(2 * pi)
    value <- vector(mode = "numeric", length = nrho)
    for(i in 1:nrho) {
        if(ncomp[i] == 1) {
            p_k1 <- pk[1, i]
            S_j <- S[, , ii[i]]
            Tht_j <- Tht[, , i]
            value[i] <- c1 * (determinant(Tht_j)$modulus - sum(S_j * Tht_j) - p_k1 * c2)
        } else {
            for(j in 1:ncomp[i]) {
                p_k1 <- pk[1, i]
                p_k2 <- p - p_k1
                id1 <- Ck[1:p_k1, i]
                id2 <- Ck[(p_k1 + 1):p, i]
                S_j <- S[id1, id1, ii[i]]
                Tht_j <- Tht[id1, id1, i]
                if(p_k1 == 1)
                    value[i] <- value[i] + c1 * (log(Tht_j) - S_j * Tht_j - c2)
                else
                    value[i] <- value[i] + c1 * (determinant(Tht_j)$modulus - sum(S_j * Tht_j) - p_k1 * c2)
                Ck[1:p_k2, i] <- id2
                Ck[(p_k2 + 1):p, i] <- id1
                pk[1:(ncomp[i] - 1), i] = pk[2:ncomp[i], i]
                pk[ncomp[i], i] = p_k1
            }
        }
    }
    out <- list(value = value, df = df, n = n, p = p, rho = object$rho, model = class(object)[1], fun = fun)
    class(out) <- "loglik"
    out
}
