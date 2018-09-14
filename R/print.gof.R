print.gof <- function (x, digits = max(3, getOption("digits") - 3), ...){
    type <- x$type
    rho <- x$rho
    df <- x$df
    val <- x$value_gof
    rnk <- rank(val)
    rnk_min <- which.min(rnk)
    rnk <- as.character(rnk)
    rnk[-rnk_min] <- paste(rnk[-rnk_min], "  ")
    rnk[rnk_min] <- paste(rnk[rnk_min], "<-")
    tbl <- data.frame(rho = rho, df = df, val = val, Rank = rnk)
    names(tbl)[3] <- type
	cat("\nSequence of", sQuote(type), "values for the estimated", sQuote(x$model),"models")
	cat("\n\nDetails:\n")
    print.data.frame(tbl, digits = digits, quote = FALSE, row.names = FALSE, ...)
	invisible(tbl)
}
