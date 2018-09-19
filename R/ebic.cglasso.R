ebic.cglasso <- function(object, g = 0.5){
    # checking 'g'
    if(!is.vector(g)) stop(sQuote("g"), " is not a vector of length ", sQuote(1))
    if(length(g) != 1) stop(sQuote("g"), " is not a vector of length ", sQuote(1))
    if((g < 0) | (g > 1)) stop(sQuote("g"), " does not belong to the interval [0, 1]")
    # extracting elements from 'object'
    n <- dim(object$X$X)[1L]
    p <- dim(object$X$X)[2L]
    # computing mle estimates
    if(class(object)[1] == "cglasso") {
        out_mle <- mle(object)
        cl <- "cglasso"
    }
    if(class(object)[1] == "cggm") {
        out_mle <- object
        cl <- "cggm"
    }
    out_loglik <- loglik(out_mle)
    df <- out_loglik$df
    value <- out_loglik$value
    rho <- out_loglik$rho
    val <-  -2 * value + df * (log(n) + 4 * g * log(p))
    out <- list(value_gof = val, rho = rho, value = value, df = df, n = n, p = p, model = cl, type = "eBIC")
    class(out) <- "gof"
    out
}
