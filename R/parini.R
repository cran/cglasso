parini <- function(object){
    fn <- function(par, x, lo, up){
        n <- length(x)
        mu <- par[1]
        sgm <- par[2]
        val <- mean(dnorm(x, mean = mu, sd = sgm, log = TRUE))
        prob <- 0
        if(up <  big) prob <- prob + pnorm(up, mean = mu, sd = sgm, log.p = FALSE)
        else prob <- prob + 1
        if(lo > -big) prob <- prob - pnorm(lo, mean = mu, sd = sgm, log.p = FALSE)
        val - log(prob)
    }
    grd <- function(par, x, lo, up){
        dfn <- double(2)
        mu <- par[1]
        sgm <- par[2]
        n <- length(x)
        z <- (x - mu) / sgm
        dfn[1] <- mean(z)
        dfn[2] <- mean(z^2) - 1
        prob <- 0
        if(up <  big) prob <- prob + pnorm(up, mean = mu, sd = sgm, log.p = FALSE)
        else prob <- prob + 1
        if(lo > -big) prob <- prob - pnorm(lo, mean = mu, sd = sgm, log.p = FALSE)
        if(up < big) {
            beta <- (up - mu) / sgm
            dfn[1] <- dfn[1] + dnorm(up, mean = mu, sd = sgm) / prob
            dfn[2] <- dfn[2] + beta * dnorm(up, mean = mu, sd = sgm) / prob
        }
        if(lo > -big) {
            alpha <- (lo - mu) / sgm
            dfn[1] <- dfn[1] - dnorm(lo, mean = mu, sd = sgm) / prob
            dfn[2] <- dfn[2] - alpha * dnorm(lo, mean = mu, sd = sgm) / prob
        }
        dfn
    }
    big <- .Machine$double.xmax / 2
    X <- object$X
    n <- dim(X)[1]
    p <- dim(X)[2]
    lo <- object$lo
    up <- object$up
    R <- event(object)
    xm <- double(p)
    vm <- double(p)
    for(j in 1:p){
        x <- X[R[, j] == 0 , j]
        if(length(x) != n){
            par.ini <- c(mean(x), sd(x))
            out <- optim(par.ini, fn = fn, gr = grd, x = x, lo = lo[j], up = up[j], control = list(fnscale = -1),
                    method = "L-BFGS-B", lower = c(-Inf, 0), upper = c(+Inf, +Inf))
            if(out$conv == 0) {
                xm[j] <- out$par[1]
                vm[j] <- out$par[2]^2
            } else {
                xm[j] <- mean(x)
                vm[j] <- mean(x^2) - xm[j]^2
            }
        } else {
            xm[j] <- mean(x)
            vm[j] <- mean(x^2) - xm[j]^2
        }
    }
    out <- list(xm = xm, vm = vm)
    out
}
