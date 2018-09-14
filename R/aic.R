aic <- function(object, k = 2){
    if(k <= 0) stop(sQuote("k"), " is not a positive value")
    type <- ifelse(k == 2, "AIC", "GoF")
    out_loglik <- loglik(object)
    value <- out_loglik$value
    df <- out_loglik$df
    p <- out_loglik$p
    n <- out_loglik$n
    rho <- out_loglik$rho
    val <- - 2 * value + k * df
    out <- list(value_gof = val, rho = rho, value = value, df = df, n = n, p = p, model = out_loglik$model, type = type)
    class(out) <- "gof"
    out
}

