\name{to_graph}
\alias{to_graph}

\title{Create Undirected Graphs}

\description{
\sQuote{\code{to_graph}} function is used to create an undirected graph from objects returned by the function \sQuote{\code{glasso}} or \sQuote{\code{cglasso}}.
}

\usage{
to_graph(object, nrho = 1L, weighted = FALSE)
}

\arguments{
\item{object}{an object of class \sQuote{\code{glasso}}, \sQuote{\code{ggm}}, \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}.}
\item{nrho}{integer specifying the model used to create the undireced graph. Default is \code{1L}.}
\item{weighted}{flag specifying whether to create a weighted graph from the adjacency matrix. Default is \code{FALSE}.}
}

\details{
The adjacency matrix \sQuote{\code{object$Adj[, , nrho]}} is passed to \link[igraph]{graph_from_adjacency_matrix}, available in the \R package \pkg{igraph}, to create the undirected graph estimated by the glasso or cglasso model. If \sQuote{\code{weighted = TRUE}} then the precision matrix \sQuote{\code{object$Tht[, , nrho]}} is used to create a weighted undirected graph. In this case, an edge associated to a negative partial correlation coefficient is plotted using a dashed line (see examples below).
}

\value{
\sQuote{\code{to_graph}} returns an \link[igraph]{igraph} object.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Abbruzzo, A. and Vinciotti, V. (2018).
\eqn{\ell_1}{l1}-Penalized gaussian graphical model.
\emph{Biostatistics} (to appear).
}

\seealso{
\link{glasso} and \link{cglasso}. For more details about the object of class \sQuote{\code{igraph}}, the interested reader is referred to the package \pkg{\link[igraph]{igraph}}.
}

\examples{
library("cglasso")
set.seed(123)

#################
# cglasso model #
#################
n <- 100L
p <- 10L
mu <- rep.int(0L, times = p)
X <- rdatacggm(n = n, mu = mu, probr = 0.05)
out <- cglasso(X)
out

# creating the undirected graph associated to the fitted
# model number 3
graph_cglasso <- to_graph(out, nrho = 3L)
graph_cglasso
V(graph_cglasso)
E(graph_cglasso)
plot(graph_cglasso, layout = layout_in_circle(graph_cglasso))

# creating the weighted graph associated to the fitted
# model number 3
graph_cglasso <- to_graph(out, nrho = 3L, weighted = TRUE)
graph_cglasso
E(graph_cglasso)
E(graph_cglasso)$weight
plot(graph_cglasso, layout = layout_in_circle(graph_cglasso))

#################
# glasso model #
#################
n <- 100L
p <- 10L
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
out <- glasso(X)
out

# creating the undirected graph associated to the fitted
# model number 3
graph_cglasso <- to_graph(out, nrho = 3L)
graph_cglasso
V(graph_cglasso)
E(graph_cglasso)
plot(graph_cglasso, layout = layout_in_circle(graph_cglasso))

# creating the weighted graph associated to the fitted
# model number 3
graph_cglasso <- to_graph(out, nrho = 3L, weighted = TRUE)
graph_cglasso
E(graph_cglasso)
E(graph_cglasso)$weight
plot(graph_cglasso, layout = layout_in_circle(graph_cglasso))
}

\keyword{multivariate}
\keyword{graphs}
