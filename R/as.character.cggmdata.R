as.character.datacggm <- function(x, ...) {
    X <- format(x$X, ...)
    marker <- event(x)
    X[marker == 1] <- paste0(X[marker == 1], "+", sep = "")
    X[marker == -1] <- paste0(X[marker == -1], "-", sep = "")
    X
}
