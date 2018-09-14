print.glasso <- function (x, digits = max(3, getOption("digits") - 3), ...){
    p <- dim(x$X)[2]
    rho <- x$rho
    R2 <- x$R2
    df <- x$df + 2 * p
    maxdf <- p + p * (p + 1) / 2
    df.per <- round(df / maxdf * 100, digits = digits)
    df.per <- paste("(", df.per, "%)", sep = "")
    ncomp <- x$ncomp
    tbl <- data.frame(1:x$nrho, rho, R2, df, df.per, ncomp)
	names(tbl) <- c("", "rho", "R2", "df", "", "N. Comp.")
	cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	print.data.frame(tbl, digits = digits, print.gap = 2, quote = FALSE, row.names = FALSE, ...)
	invisible(tbl)
}
