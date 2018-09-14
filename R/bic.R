bic <- function(object){
    n <- dim(object$X)[1L]
    out <- aic(object, k = log(n))
    out$type <- "BIC"
    out
}

