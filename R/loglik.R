loglik <- function(object) {
    fun <- ifelse(is.element(class(object)[1], c("glasso", "ggm")), "log-Lik", "Q-Fun")
    X <- as.matrix(object$X)
    n <- dim(X)[1L]
    p <- dim(X)[2L]
    nrho <- object$nrho
    S <- object$S
    Tht <- object$Tht
    Ck <- object$Ck
    pk <- object$pk
    ncomp <- object$ncomp
    df <- object$df + 2 * p
    if(length(dim(S)) == 2L) {
        dim(S) <- c(p, p, 1L)
        ii <- rep.int(1L, nrho)
    } else ii <- seq_len(nrho)
    c1 <- 0.5 * n
    c2 <- log(2 * pi)
    value <- double(nrho)
    for(i in seq_len(nrho)) {
        if(ncomp[i] == 1L) {
            p_k1 <- pk[1L, i]
            S_j <- S[, , ii[i]]
            Tht_j <- Tht[, , i]
            value[i] <- c1 * (determinant(Tht_j)$modulus - sum(S_j * Tht_j) - p_k1 * c2)
        } else {
            for(j in seq_len(ncomp[i])) {
                p_k1 <- pk[1L, i]
                p_k2 <- p - p_k1
                id1 <- Ck[1L:p_k1, i]
                id2 <- Ck[(p_k1 + 1L):p, i]
                S_j <- S[id1, id1, ii[i]]
                Tht_j <- Tht[id1, id1, i]
                if(p_k1 == 1L)
                    value[i] <- value[i] + c1 * (log(Tht_j) - S_j * Tht_j - c2)
                else
                    value[i] <- value[i] + c1 * (determinant(Tht_j)$modulus - sum(S_j * Tht_j) - p_k1 * c2)
                Ck[1L:p_k2, i] <- id2
                Ck[(p_k2 + 1L):p, i] <- id1
                pk[1L:(ncomp[i] - 1L), i] = pk[2L:ncomp[i], i]
                pk[ncomp[i], i] = p_k1
            }
        }
    }
    out <- list(value = value, df = df, n = n, p = p, rho = object$rho, model = class(object)[1L], fun = fun)
    class(out) <- "loglik"
    out
}
