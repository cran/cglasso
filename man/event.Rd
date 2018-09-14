\name{event}
\alias{event}

\title{Return the Indicator Matrix from an Object with class \sQuote{\code{datacggm}}}

\description{
The \sQuote{\code{event}} function is used to create a status indicator matrix from an object with class \sQuote{\code{datacggm}}. The elements of the matrix, denoted by \code{R}, are used to specify the status of an observation:
\itemize{
\item \sQuote{\code{R[i, j] =  0}} means that the \eqn{i}{i}th observation of the \eqn{j}{j}th random variable is observed;
\item \sQuote{\code{R[i, j] = -1}} means that the \eqn{i}{i}th observation of the \eqn{j}{j}th random variable is left-censored;
\item \sQuote{\code{R[i, j] = +1}} means that the \eqn{i}{i}th observation of the \eqn{j}{j}th random variable is right-censored.
}
See examples below.
}

\usage{
event(x)
}

\arguments{
\item{x}{an object with class \sQuote{\code{datacggm}}.}
}

\value{
\code{event} returns a \eqn{(n\times p)}{(n x p)}-dimensional matrix.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Abbruzzo, A. and Vinciotti, V. (2018).
\eqn{\ell_1}{l1}-Penalized gaussian graphical model.
\emph{Biostatistics} (to appear).
}

\seealso{
\code{\link{datacggm}}, \code{\link{rdatacggm}} and the method function \code{\link{summary.datacggm}}.
}

\examples{
set.seed(123)
library("cglasso")

# dataset from a left-censored Gaussian graphical model
n <- 100L
p <- 5L
X <- matrix(rnorm(n * p), n, p)
lo <- -1
X[X <= lo] <- lo
X <- datacggm(X, lo = lo)
event(X)

# dataset from a right-censored Gaussian graphical model
n <- 100L
p <- 5L
X <- matrix(rnorm(n * p), n, p)
up <- 1
X[X >= up] <- up
X <- datacggm(X, up = up)
event(X)

# dataset from a censored Gaussian graphical model
n <- 100L
p <- 5L
X <- matrix(rnorm(n * p), n, p)
up <- 1
lo <- -1
X[X >= up] <- up
X[X <= lo] <- lo
X <- datacggm(X, lo = lo, up = up)
event(X)
}

\keyword{datagen}
\keyword{multivariate}
