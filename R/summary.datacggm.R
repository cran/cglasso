summary.datacggm <- function(object, digits = max(3L, getOption("digits") - 3L), ...) {
    R <- event(X)
    X <- object$X
    nm <- colnames(X)
    n <- dim(X)[1]
    p <- dim(X)[2]
    lo <- object$lo
    lo[lo == -.Machine$double.xmax] <- -Inf
    up <- object$up
    up[up == .Machine$double.xmax] <- Inf
    out <- vector(mode = "list", length = p)
    lbs1 <- format(c("Lower", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Upper"))
    lbs2 <- format(c("N. Left C.", "N. Obs.", "N. Right C."))
    nr <- length(lbs1) + length(lbs2) + 1
    for(m in seq_len(p)) {
        xobs <- X[R[, m] == 0, m]
        nobs <- length(xobs)
        nrc <- sum(R[, m] == 1)
        nlc <- n - nobs - nrc
        qq <- stats::quantile(xobs)
        qq <- c(lo[m], qq[1L:3L], mean(xobs), qq[4L:5L], up[m])
        qq <- format(qq, digits = digits, ...)
        ninfo <- c(nlc, nobs, nrc)
        out[[m]] <- c(paste0(lbs1, " : ", qq, "  "), "---", paste0(lbs2, " : ", ninfo, "  "))
    }
    out <- unlist(out)
    nchar1 <- max(nchar(out))
    nchar2 <- max(nchar(colnames(X)))
    if(nchar2 >= nchar1) nm <- format(colnames(X), justify = "centre")
    else {
        nspace <- floor((nchar1 - nchar2) / 2 + nchar2 / 2)
        space <- paste(rep(" ", nspace), collapse = "")
        nm <- format(c(space, colnames(X)), justify = "right")[-1]
    }
    dim(out) <- c(nr, p)
    dimnames(out) <- list(rep.int("", nr), nm)
    attr(out, "class") <- c("table")
    out
}
