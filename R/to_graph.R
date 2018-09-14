to_graph <- function(object, nrho = 1L, weighted = FALSE) {
    # checking 'nrho'
    if(!is.vector(nrho)) stop(sQuote("nrho"), " is not a vector of length ", sQuote(1))
    if(length(nrho) != 1) stop(sQuote("nrho"), " is not a vector of length ", sQuote(1))
    if(abs(as.integer(nrho) - nrho) > 0) stop(sQuote("nrho"), " is not an object of type ", dQuote("integer"))
    if(nrho <= 0) stop(sQuote("nrho"), " is not a positive integer")
    if(nrho > object$nrho) stop("'nrho' can not be larger than ", sQuote(object$nrho))
    # checking 'weighted'
    if(!is.logical(weighted)) stop(sQuote("weighted"), " is not an object of type ", dQuote("logical"))
    if(!weighted) {
        Adj <- object$Adj[, , nrho]
        weighted <- NULL
    }
    else Adj <- object$Tht[, , nrho]
    out <- graph_from_adjacency_matrix(adjmatrix = Adj, mode = "undirected", diag = FALSE, weighted = weighted)
    if(!is.null(weighted)) E(out)$lty <- ifelse(E(out)$weight > 0, "solid", "dashed")
    out
}
