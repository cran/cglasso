event <- function(x) {
    if(!is.datacggm(x)) stop("'x' is not an object of class 'datacggm'")
    R <- x$R[-1, -c(1, 2), drop = FALSE]
    R
}
