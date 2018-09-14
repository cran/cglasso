rdatacggm <- function(n, mu, Sigma, probl, probr, lo, up, ...) {
    if(missing(mu) & missing(Sigma)) stop("Please, specify at least one of the following arguments: 'mu' or 'Sigma'")
    if(missing(probl) & missing(lo) & missing(probr) & missing(up)) stop("Please, specify at least one of the following arguments: 'probl', 'probr', 'lo' or 'up'")
    # setting default values for 'mu' or 'Sigma'
    if(missing(mu)) {
        p <- dim(Sigma)[1L]
        mu <- rep.int(0L, p)
    }
    if(missing(Sigma)) {
        p <- length(mu)
        Sigma <- diag(p)
    }
    # checking dimensions for 'mu' and 'Sigma'
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p))) stop("'mu' and 'Sigma' have incompatible dimensions")
    xsd <- sqrt(diag(Sigma))
    # handle 'probl' and 'lo'
    if(!missing(probl)) {
        if(!missing(lo)) warning("'lo' is overwritten using 'probl'")
        if(!is.vector(probl)) stop("'probl' is not a vector")
        if(length(probl) == 1) probl <- rep(probl, p)
        if(!all(0 <= probl & probl <= 1)) stop("some element in 'probl' does not belong to the interval [0, 1]")
        lo <- qnorm(p = probl, mean = mu, sd = xsd, lower.tail = TRUE)
    } else {
        if(missing(lo)) probl <- rep(0L, p)
        else probl <- pnorm(lo, mean = mu, sd = xsd, lower.tail = TRUE)
    }
    # handling 'probr' and 'up'
    if(!missing(probr)) {
        if(!missing(up)) warning("'up' is overwritten using 'probr'")
        if(!is.vector(probr)) stop("'probr' is not a vector")
        if(length(probr) == 1) probr <- rep(probr, p)
        if(!all(0 <= probr & probr <= 1)) stop("some element in 'probr' does not belong to the interval [0, 1]")
        up <- qnorm(p = probr, mean = mu, sd = xsd, lower.tail = FALSE)
    } else {
        if(missing(up)) probr <- rep(0L, p)
        else probr <- pnorm(up, mean = mu, sd = xsd, lower.tail = FALSE)
    }
    # checking 'probl' and 'probr'
    if(any(probl + probr >= 1)) stop("'probl' plus 'probr' can not be greater than or equal to 1")
    # sampling from a multivariate Gaussian distribution
    if(all(Sigma[upper.tri(Sigma)] == 0)) X <- apply(cbind(mu, sqrt(diag(Sigma))), 1, function(par) rnorm(n = n, mean = par[1], sd = par[2]))
    else X <- mvrnorm(n = n, mu = mu, Sigma = Sigma, ...)
    if(n == 1) X <- t(X)
    # handling censored vaules
    if(!missing(lo)) {
        for(m in seq_len(p)) X[X[, m] <= lo[m], m] <- lo[m]
    }
    if(!missing(up)) {
        for(m in seq_len(p)) X[X[, m] >= up[m], m] <- up[m]
    }
    # generating an object with S3 class 'datacggm'
    X <- datacggm(X = X, lo = lo, up = up)
    X
}
