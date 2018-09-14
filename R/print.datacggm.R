print.datacggm <- function(x, quote = FALSE, ...) {
    invisible(print(as.character.datacggm(x, ...), quote = quote))
}
