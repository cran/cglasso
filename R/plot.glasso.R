plot.glasso <- function(x, typeplot = c("path", "graph"), gof, diag = FALSE, nrho, weighted = FALSE, ...) {
    typeplot <- match.arg(typeplot)
    # checking 'gof'
    if(!missing(gof)) {
        if(class(gof) != "gof") stop("argument 'gof' is not an object with class 'gof'")
        val <- gof$value_gof
        nrho_opt <- which.min(val)
    } else gof <- NULL
    # checking 'diag'
    if(!is.logical(diag)) stop(sQuote("diag"), " is not an object of type ", dQuote("logical"))
    if(!is.vector(diag)) stop(sQuote("diag"), " is not a vector of length ", sQuote(1))
    if(length(diag) != 1) stop(sQuote("diag"), " is not a vector of length ", sQuote(1))
    # checking 'nrho'
    if(!missing(nrho)) {
        if(!is.vector(nrho)) stop(sQuote("nrho"), " is not a vector of length ", sQuote(1))
        if(length(nrho) != 1) stop(sQuote("nrho"), " is not a vector of length ", sQuote(1))
        if(abs(as.integer(nrho) - nrho) > 0) stop(sQuote("nrho"), " is not an object of type ", dQuote("integer"))
        if(nrho <= 0) stop(sQuote("nrho"), " is not a positive integer")
        if(nrho > x$nrho) stop(sQuote("nrho"), " is equal to ", nrho, ", but it must be less than or equal to ", x$nrho)
    }
    if(x$nrho == 1) nrho <- x$nrho
    # if 'x$nrho == 1' we can not plot the coefficients path
    if((typeplot == "path") & (x$nrho == 1)) typeplot <- "graph"
    # plotting coefficient path
    if(typeplot == "path") {
        # setting graphical parameters
        linetype <- "l"
        col <- 1
        lty <- 1
        rho <- x$rho
        U <- upper.tri(x$Tht[, , 1], diag = diag)
        Tht <- t(apply(x$Tht, 3, function(M) M[U]))
        matplot(rho, Tht, type = linetype, col = col, lty = lty, xlab = expression(rho), ylab = expression(hat(Theta)), ...)
        if(!is.null(gof)) abline(v = rho[nrho_opt], lty = 2, lwd = 2)
    }
    # plotting graph
    if(typeplot == "graph") {
        if(!is.null(gof)) {
            if(!missing(nrho)) warning("argument 'nrho' is overwritten using 'gof'")
            nrho <- nrho_opt
        }
        out_graph <- to_graph(x, nrho = nrho, weighted = weighted)
        plot(out_graph, ...)
        invisible(out_graph)
    }
}
