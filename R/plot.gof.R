plot.gof <- function(x, ...) {
    rho <- x$rho
    val <- x$value_gof
    minval <- which.min(val)
    dots <- list(...)
    if(is.null(dots$xlab)) dots$xlab <- expression(rho)
    if(is.null(dots$ylab)) {
        ylab <- switch(x$type,
            "AIC" = "Akaike Information Criterion",
            "BIC" = "Bayesian Information Criterion",
            "GoF" = "Measure of Goodness of Fit",
            "eBIC" = "extended Bayesian Information Criterion")
        dots$ylab <- ylab
    }
    if(is.null(dots$type)) dots$type <- "b"
    if(is.null(dots$main)) dots$main <- "Tuning Parameter Selection"
    do.call(function(...) plot(rho, val, ...), dots)
    abline(v = rho[minval], lty = 2, lwd = 2)
}
