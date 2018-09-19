coef.glasso <- function(object, ..., nrho = 1L, type = c("theta", "sigma"), print.info = FALSE, digits = 3L){
    # checkin 'nrho'
    if(!is.vector(nrho)) stop(sQuote("nrho"), " is not a vector of length ", sQuote(1))
    if(length(nrho) != 1) stop(sQuote("nrho"), " is not a vector of length ", sQuote(1))
    if(abs(as.integer(nrho) - nrho) > 0) stop(sQuote("nrho"), " is not an object of type ", dQuote("integer"))
    if(nrho <= 0) stop(sQuote("nrho"), " is not a positive integer")
    if(nrho > object$nrho) stop("'nrho' can not be larger than ", sQuote(object$nrho))
    #stop("model is fitted using nrho = ", object$nrho)
    type <- match.arg(type)
    if(!is.logical(print.info)) stop(sQuote("print.info"), " is not an object of type ", dQuote("logical"))
    out <- switch(type,
        "theta" = list(object$Tht[, , nrho], "Precision matrix"),
        "sigma" = list(object$Sgm[, , nrho], "Covariance matrix"))
    coef.mat <- out[[1]]
    msg <- out[[2]]
    if(print.info) {
        # extracting elements from 'object'
        X <- object$X
        p <- dim(X)[2]
        Ck <- object$Ck[, nrho]
        pk <- object$pk[, nrho]
        ncomp <- object$ncomp[nrho]
        rho <- object$rho[nrho]
        R2 <- round(object$R2[nrho], digits = digits)
        # printing first part
        cat("\nResults of the estimated", sQuote(class(object)[1]), "model")
        cat("\n\n",
        apply(cbind(format(c("rho:", "R2:", "Number of connected components:", "Number of vertices per component:"), justify = "right"),
        c(rho, R2, ncomp, paste(pk[1:ncomp], collapse = " ")), rep("\n", 4L)), 1L, paste, collapse = " "), "\n")
        # printing second part
        if(ncomp == 1) {
            cat(msg, ":\n\n", sep = "")
            print.table(coef.mat, digits = digits, zero.print = ".", ...)
            cat("\n")
        } else {
            msg <- paste(msg, "(connected component: ")
            for(k in 1:ncomp) {
                p_k1 <- pk[1]
                p_k2 <- p - p_k1
                id1 <- Ck[1:p_k1]
                id2 <- Ck[(p_k1 + 1):p]
                cat(msg, k, ")\n\n", sep = "")
                print.table(coef.mat[sort(id1), sort(id1), drop = FALSE], digits = digits, zero.print = ".", ...)
                cat("\n")
                Ck[1:p_k2] <- id2
                Ck[(p_k2 + 1):p] <- id1
                pk[1:(ncomp - 1)] = pk[2:ncomp]
                pk[ncomp] = p_k1
            }
        }
        invisible(coef.mat)
    } else coef.mat
}
