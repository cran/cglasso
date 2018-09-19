scale.datacggm <- function(x, center = TRUE, scale = TRUE) {
    big <- .Machine$double.xmax
    p <- dim(x)[2L]
    X <- x$X
    R <- event(x)
    up <- x$up
    lo <- x$lo
    if(is.logical(center)) {
        if(center) {
            center <- double(p)
            for(j in seq_len(p)) center[j] <- mean(X[R[, j] == 0, j])
        }
    }
    if (is.logical(scale)) {
        if (scale) {
            f <- function(v) {
                v <- v[!is.na(v)]
                sqrt(sum(v^2)/max(1, length(v) - 1L))
            }
            scale <- double(p)
            for(j in seq_len(p)) scale[j] <- f(X[R[, j] == 0, j])
        }
    }
    X <- scale.default(x = X, center = center, scale = scale)
    center <- attr(X, "scaled:center")
    if(is.numeric(center)) {
        id <- up < big
        if(any(id)) up[id] <- up[id] - center[id]
        id <- -big < lo
        if(any(id)) lo[id] <- lo[id] - center[id]
        attr(x, "scaled:center") <- center
    }
    scale <- attr(X, "scaled:scale")
    if(is.numeric(scale)) {
        id <- up < big
        if(any(id)) up[id] <- up[id] / scale[id]
        id <- -big < lo
        if(any(id)) lo[id] <- lo[id] / scale[id]
        attr(x, "scaled:scale") <- scale
    }
    x$X <- X
    x$up <- up
    x$lo <- lo
    x
}
