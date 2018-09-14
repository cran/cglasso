print.loglik <- function (x, digits = max(3, getOption("digits") - 3), ...){
    tbl <- data.frame(rho = x$rho, df = x$df, value = x$value)
    names(tbl)[3] <- ifelse(x$fun == "log-Lik", "log-Lik.", "Q-values")
    cat("\nSequence of the", ifelse(x$fun == "log-Lik", "log-likelihood values", "Q-values") , "for the estimated", sQuote(x$model), "models")
    cat("\n\nDetails:\n")
    print.data.frame(tbl, digits = digits, quote = FALSE, row.names = FALSE, ...)
    invisible(tbl)
}
