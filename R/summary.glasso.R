summary.glasso <- function(object, ..., gof = c("BIC", "AIC", "eBIC"), par.gof, digits = 4L){
    gof <- match.arg(gof)
    if(gof == "AIC" & missing(par.gof)) par.gof <- 2
    if(gof == "eBIC" & missing(par.gof)) par.gof <- 0.5
    out_tbl <- make_summary_table(object, gof = gof, par.gof = par.gof, digits = digits)
    rnk_min <- out_tbl$rnk_min
    cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    print.data.frame(out_tbl$tbl, print.gap = 2, quote = FALSE, row.names = FALSE, digits = digits, ...)
    # printing second part
    lbl <- c("Model", "rho", "which.min", names(out_tbl$tbl)[6L], "R2", "df", "Number of Components", "Sizes")
    lbl <- paste("\n", format(lbl, justify = "right"), ":", sep = "")
    cat("\n\n===============================================================")
    cat("\n\nSummary of the Selected Model\n")
    cat(lbl[1], sQuote(class(object)[1L]))
    cat(lbl[2], out_tbl$tbl[rnk_min, "rho"])
    cat(lbl[3], rnk_min)
    cat(lbl[4], out_tbl$tbl[rnk_min, 6L])
    cat(lbl[5], out_tbl$tbl[rnk_min, "R2"])
    cat(lbl[6], paste(out_tbl$tbl[rnk_min, "df"], as.character(out_tbl$tbl[rnk_min, 4L]), sep = " "))
    cat(lbl[7], out_tbl$ncomp)
    cat(lbl[8], out_tbl$pk)
    cat("\n\n===============================================================\n\n")
    out <- out_tbl[1L:2L]
    names(out) <- c("table", "which.min")
    invisible(out)
}
