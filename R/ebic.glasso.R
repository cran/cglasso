ebic.glasso <- function(object, g = 0.5){
    # checking 'g'
    if(!is.vector(g)) stop(sQuote("g"), " is not a vector of length ", sQuote(1))
    if(length(g) != 1) stop(sQuote("g"), " is not a vector of length ", sQuote(1))
    if((g < 0) | (g > 1)) stop(sQuote("g"), " does not belong to the interval [0, 1]")
    # extracting elements from 'object'
    n <- dim(object$X)[1]
    p <- dim(object$X)[2]
    # computing mle estimates
    if(class(object)[1] == "glasso") {
        out_mle <- mle(object)
        cl <- "glasso"
    }
    if(class(object)[1] == "ggm") {
        out_mle <- object
        cl <- "ggm"
    }
    out_loglik <- loglik(object)
    df <- out_loglik$df
    value <- out_loglik$value
    rho <- out_loglik$rho
    val <-  -2 * value + df * (log(n) + 4 * g * log(p))
    out <- list(value_gof = val, rho = rho, value = value, df = df, n = n, p = p, model = out_loglik$model, type = "eBIC")
    class(out) <- "gof"
    out
}
