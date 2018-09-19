make_summary_table <- function(object, gof, par.gof, digits){
    X <- as.matrix(object$X)
    n <- dim(X)[1L]
    p <- dim(X)[2L]
    rho <- object$rho
    R2 <- object$R2
    df <- object$df + 2 * p
    maxdf <- p + p * (p + 1) / 2
    df.per <- round(df / maxdf * 100, digits = digits)
    df.per <- paste("(", df.per, "%)", sep = "")
    out_gof <- switch(gof,
                    "AIC" = aic(object, k = par.gof),
                    "BIC" = bic(object),
                    "eBIC" = ebic(object, g = par.gof))
    GoF <- out_gof$value_gof
    rnk <- rank(GoF)
    rnk_min <- which.min(rnk)
    rnk <- as.character(rnk)
    rnk[-rnk_min] <- paste(rnk[-rnk_min], "  ")
    rnk[rnk_min] <- paste(rnk[rnk_min], "<-")
    type <- out_gof$type
    tbl <- data.frame(1:object$nrho, rho, df, df.per, R2, GoF, rnk)
    names(tbl) <- c("", "rho", "df", "", "R2", type, "Rank")
    ncomp <- object$ncomp[rnk_min]
    pk <- object$pk[1:ncomp, rnk_min]
    out <- list(tbl = tbl, rnk_min = rnk_min, ncomp = ncomp, pk = pk)
    out
}
