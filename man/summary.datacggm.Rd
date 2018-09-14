\name{summary.datacggm}
\alias{summary.datacggm}

\title{Summarizing Objects of Class \sQuote{\code{cggmdata}}}

\description{
\sQuote{\code{summary.datacggm}} function is used to produce summaries of an object of class \sQuote{\code{datacggm}}.
}

\usage{
\method{summary}{datacggm}(object, digits = max(3L, getOption("digits") - 3L), \dots)
}

\arguments{
\item{object}{an object of class \sQuote{\code{datacggm}}.}
\item{digits}{integer used for number formatting with \sQuote{\code{format()}}.}
\item{\dots}{further arguments passed to \sQuote{\code{format()}}.}
}

\details{
The method function \sQuote{\code{summary.datacggm}} extends the results given by \sQuote{\code{summary.matrix()}} taking into account the censoring values. For each variable, the mean and the quartiles are computing using only the observed values; the lower and upper censoring values (denoted by \sQuote{\code{Lower}} and \sQuote{\code{Upper}}) are also reported. The number of observed and censured values are computed and showed in the second part of the output (see example below).
}

\value{
\sQuote{\code{summary.dataccgm}} returns a matrix with class \sQuote{\code{table}}, obtained by computing the summary measures to each variable and collating the results.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Abbruzzo, A. and Vinciotti, V. (2018).
\eqn{\ell_1}{l1}-Penalized gaussian graphical model.
\emph{Biostatistics} (to appear).
}

\seealso{
\code{\link{datacggm}} and \code{\link{rdatacggm}}.
}

\examples{
set.seed(123)
library("cglasso")

n <- 1000L
p <- 3L
mu <- rep(1L, p)
rho <- 0.3
Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
X <- rdatacggm(n = n, mu = mu, Sigma = Sigma, probr = 0.05, probl = 0.05)
summary(X)
}

\keyword{datagen}
\keyword{multivariate}

