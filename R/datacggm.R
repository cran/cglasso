datacggm <- function(X, lo, up) {
    big <- .Machine$double.xmax
    thr <- big / 2
    # checking 'X'
    if(missing(X)) stop(sQuote("X"), " is missing")
    if(!is.matrix(X)) stop(sQuote("X"), " is not a matrix")
    n <- dim(X)[1]
    p <- dim(X)[2]
    if(is.null(colnames(X))) {
        vnames <- paste("X", 1:p, sep = "")
        colnames(X) <- vnames
    } else vnames <- colnames(X)
    # checking 'lo' and 'up'
    if(missing(lo) & missing(up)) stop("arguments ", sQuote("lo"), " and ", sQuote("up"), " are missing")
    if(missing(lo)) lo <- rep(-big, p)
    else {
        if(!is.vector(lo)) stop(sQuote("lo"), " is not a vector")
        if(length(lo) == 1) lo <- rep(lo, p)
        id <- -big <= lo & lo <= -thr
        if(any(id)) {
            lo[id] <- -big
            warning("some elements in 'lo' are below the tolerance. These values are treated as -Inf")
        }
        id <- lo == -Inf
        if(any(id)) lo[id] <- -big
    }
    names(lo) <- vnames
    if(missing(up)) up <- rep(big, p)
    else {
        if(!is.vector(up)) stop(sQuote("up"), " is not a vector")
        if(length(up) == 1) up <- rep(up, p)
        id <- thr <= up & up <= big
        if(any(id)) {
            up[id] <- big
            warning("some elements in 'up' are over the tolerance. These values are treated as +Inf")
        }
        id <- up == Inf
        if(any(id)) up[id] <- big
    }
    names(up) <- vnames
    if(!all(lo < up)) stop(sQuote("lo"), " is not less than ", sQuote("up"))
    R <- matrix(0, nrow = (n + 1), ncol = p + 2)
    colnames(R) <- c("nlc", "nrc", vnames)
    # storage.mode
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(X) <- "double"
    storage.mode(up) <- "double"
    storage.mode(R) <- "integer"
    startmis <- integer(1)
    out <- .Fortran(C_setup, n = n, p = p, X = X, lo = lo, up = up, R = R, startmis = startmis)
    if(out$startmis == 0) warning("the dataset is full observed. Use 'glasso' to fit the l1-penalized Gaussian graphical model")
    out <- out[-c(1, 2)]
    class(out) <- "datacggm"
    out
}

"[.datacggm" <- function(x, i, j, drop = FALSE) {
    if(!missing(j)) {
        x$X <- x$X[, j, drop = drop]
        x$R <- cbind(x$R[, c(1, 2), drop = FALSE], x$R[, -c(1, 2)][, j, drop = drop])
        x$lo <- x$lo[j]
        x$up <- x$up[j]
    }
    if(!missing(i)) {
        X <- x$X[i, , drop = drop]
        lo <- x$lo
        lo[lo == -.Machine$double.xmax] <- -Inf
        up <- x$up
        up[up == .Machine$double.xmax] <- +Inf
        x <- datacggm(X = X, lo = lo, up = up)
    }
    x
}

is.datacggm <- function(x) inherits(x, 'datacggm')

dim.datacggm <- function(x) dim(x$X)

dimnames.datacggm <- function(x) dimnames(x$X)

`dimnames<-.datacggm` <- function(x, value) {
    x$X <- `dimnames<-`(x$X, value)
    x
}

Math.datacggm <- function(...)  stop("Invalid operation on a 'datacggm' object")

Ops.datacggm  <- function(...)  stop("Invalid operation on a 'datacggm' object")

as.matrix.datacggm <- function(x, ...) x$X






