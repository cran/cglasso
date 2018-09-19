\name{plot}

\alias{plot}
\alias{plot.glasso}
\alias{plot.mglasso}

\title{Plot for \sQuote{\code{glasso}} and \sQuote{\code{cglasso}} Object}

\description{
The method functions \sQuote{\code{plot.glasso}} and \sQuote{\code{plot.mglasso}} produce plots to study the sequence of models estimated by \code{\link{glasso}}, \code{\link{mglasso}} or \code{\link{cglasso}}.
}

\usage{
\method{plot}{glasso}(x, typeplot = c("path", "graph"),
     gof, diag = FALSE, nrho, weighted = FALSE,
     isolated = FALSE, \dots)

\method{plot}{mglasso}(x, typeplot = c("path", "graph"),
     gof, diag = FALSE, nrho, weighted = FALSE,
     isolated = FALSE, \dots)
}

\arguments{
\item{x}{a fitted model object.}
\item{typeplot}{a string specifying the produced plot.}
\item{gof}{an object with class \sQuote{\code{gof}}. See examples below.}
\item{diag}{flag specifying whether the diagonal elements of the concentration matrix are plotted.}
\item{nrho}{availables only if \sQuote{\code{typeplot = graph}}: integer specifying the model used to procude the \sQuote{\code{igraph}} object. This argument is overwritten if \sQuote{\code{gof}} is available. See examples below.}
\item{weighted}{argument availables only if \sQuote{\code{typeplot = graph}}. Flag specifying whether to create a weighted graph.}
\item{isolated}{flag specifying whether the isolated vertices are removed from the graph. Default is \code{FALSE}, i.e. the isolated vertices are removed.}
\item{\dots}{additional arguments passed to the method function \sQuote{\code{plot.default}}.}
}

\details{
The plot produced by the method functions \sQuote{\code{plot.glasso}} and \sQuote{\code{plot.mglasso}} depends on the argument \sQuote{\code{typeplot}}.

If \sQuote{\code{typeplot = path}}, the regularization paths are produced; in this case, if an object with class \sQuote{\code{gof}} is passed by the argument \sQuote{\code{gof}}, then a vertical dashed line is added to identify the optimal \eqn{\rho}{rho}-value.

If \sQuote{\code{typeplot = graph}}, the method functions \sQuote{\code{plot.glasso}} and \sQuote{\code{plot.mglasso}} produce the undirected graph associated to the model specified by the argument \sQuote{\code{nrho}}. If an object with class \sQuote{\code{gof}} is passed by \sQuote{\code{gof}}, the undirected graph of the model selected by the function \sQuote{\link{aic}}, \sQuote{\link{bic}} or \sQuote{\link{ebic}} is produced.
}

\value{
If \code{typeplot = "graph"} then the igraph object is returned (see example below).
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{glasso}}, \code{\link{mglasso}} and \code{\link{cglasso}}.
}

\examples{
library("cglasso")
set.seed(123)

#################
# cglasso model #
#################
n <- 100L
p <- 5L
mu <- rep.int(0L, times = p)
X <- rdatacggm(n = n, mu = mu, probr = 0.05)
out <- cglasso(X = X, nrho = 100L)
out_aic <- aic(out)

# plotting the regularization paths + 'gof' object
plot(out, typeplot = "path")
plot(out, typeplot = "path", gof = out_aic)

# plotting the graph associated to the fitted model
# specified by 'nrho'
out_graph <- plot(out, typeplot = "graph", nrho = 10L)
out_graph
# plotting the graph associated to the fitted model
# specified by 'gof'
out_graph <- plot(out, typeplot = "graph", gof = out_aic)
out_graph

#################
# mglasso model #
#################
R <- event(X)
X <- as.matrix(X)
X[R == 1L] <- NA
out <- mglasso(X = X, nrho = 100L)
out_aic <- aic(out)

# plotting the regularization paths + 'gof' object
plot(out, typeplot = "path")
plot(out, typeplot = "path", gof = out_aic)

# plotting the graph associated to the fitted model
# specified by 'nrho'
out_graph <- plot(out, typeplot = "graph", nrho = 10L)
out_graph
# plotting the graph associated to the fitted model
# specified by 'gof'
out_graph <- plot(out, typeplot = "graph", gof = out_aic)
out_graph

#################
# glasso model #
#################
X <- MASS::mvrnorm(n = n, mu = mu, Sigma = diag(p))
out <- glasso(X, nrho = 100L)
out_aic <- aic(out)

# plotting the regularization paths + 'gof' object
plot(out, typeplot = "path")
plot(out, typeplot = "path", gof = out_aic)

# plotting the graph associated to the fitted model
# specified by 'nrho'
out_graph <- plot(out, typeplot = "graph", nrho = 10L)
out_graph
# plotting the graph associated to the fitted model
# specified by 'gof'
out_graph <- plot(out, typeplot = "graph", gof = out_aic)
out_graph
}

\keyword{hplot}
\keyword{graphs}
