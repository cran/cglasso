\name{datacggm}
\alias{datacggm}
\alias{[.datacggm}
\alias{dimnames.datacggm}
\alias{dimnames<-.datacggm}
\alias{is.datacggm}
\alias{print.datacggm}
\alias{dim.datacggm}
\alias{as.matrix.datacggm}

\title{Create a Dataset from a Censored Gaussian Graphical Model}

\description{
\sQuote{\code{datacggm}} function is used to create a dataset from a censored Gaussian graphical model.
}

\usage{
datacggm(X, lo, up)
}

\arguments{
\item{X}{a \eqn{(n\times p)}{(n x p)}-dimensional matrix; each row is an observation from a censored Gaussian graphical model with censoring vectors \code{lo} and \code{up}.}
\item{lo}{the lower censoring vector; \code{lo[j]} is used to specify the lower censoring value for the random variable \eqn{X_j}{X_j}.}
\item{up}{the upper censoring vector; \code{up[j]} is used to specify the upper censoring value for the random variable \eqn{X_j}{X_j}.}
}

\details{
The function \sQuote{\code{datacggm}} returns a named list with class \sQuote{\code{datacggm}} containing the elements needed to fit a censored graphical lasso (cglasso) model. In output, the matrix \code{X} is ordered according to the patter of censoring values.

There are specific method functions developed to help the user to deal with the censored values. The \sQuote{\code{print.datacggm}} method function print out the left and right-censored values using the following rules: a right-censored value is labeld adding the symbol \sQuote{\code{+}} at the end of the value, whereas the symbol \sQuote{\code{-}} is used for the left-censored values (see examples bellow). The summary statistics about the censored values can be obtained using the method function \sQuote{\code{\link{summary.datacggm}}}.

Finally, the status indicator matrix, denoted by \code{R}, can be obtained by the function \code{\link{event}}. The elements of this matrix specify the status of an observation as follows:
\itemize{
\item \sQuote{\code{R[i, j] =  0}} means that the \eqn{i}{i}th observation of the \eqn{j}{j}th random variable is observed;
\item \sQuote{\code{R[i, j] = -1}} means that the \eqn{i}{i}th observation of the \eqn{j}{j}th random variable is left-censored;
\item \sQuote{\code{R[i, j] = +1}} means that the \eqn{i}{i}th observation of the \eqn{j}{j}th random variable is right-censored.
}
}

\value{
\sQuote{\code{datacggm}} returns an object with S3 class \dQuote{\code{datacggm}}, i.e. a list containing the
following components:
\item{X}{the \eqn{(n\times p)}{(n x p)}-dimensional matrix \code{X} ordered according to the patterns of censored values.}
\item{lo}{the lower censoring vector.}
\item{up}{the upper censoring vector.}
\item{R}{the augmented status indicator matrix encoding the patterns of censored values (for internal purposes only); the status indicator matrix is returned by function \code{\link{event}}.}
\item{startmis}{the row of the matrix \code{X} where are starting the patterns of censored values (for internal purposes only).}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Abbruzzo, A. and Vinciotti, V. (2018).
\eqn{\ell_1}{l1}-Penalized gaussian graphical model.
\emph{Biostatistics} (to appear).
}

\seealso{
\code{\link{event}}, \code{\link{rdatacggm}}, \code{\link{cglasso}} and the method function \code{\link{summary.datacggm}}.
}

\examples{
set.seed(123)
library("cglasso")

# a dataset from a left-censored Gaussian graphical model
n <- 100L
p <- 5L
X <- matrix(rnorm(n * p), n, p)
lo <- -1
X[X <= lo] <- lo
X <- datacggm(X, lo = lo)
X

# a dataset from a right-censored Gaussian graphical model
n <- 100L
p <- 5L
X <- matrix(rnorm(n * p), n, p)
up <- 1
X[X >= up] <- up
X <- datacggm(X, up = up)
X

# a dataset from a censored Gaussian graphical model
n <- 100L
p <- 5L
X <- matrix(rnorm(n * p), n, p)
up <- 1
lo <- -1
X[X >= up] <- up
X[X <= lo] <- lo
X <- datacggm(X, lo = lo, up = up)
X
}

\keyword{datagen}
\keyword{multivariate}
