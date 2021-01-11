\name{nObs + nResp + nPred}

\alias{nObs}
\alias{nResp}
\alias{nPred}

\title{Extract the Number of Observations/Responses/Predictors from a datacggm Object}

\usage{
nObs(x)
nResp(x)
nPred(x)
}

\description{
Extract the number of observations, response variables and predictors from an object of class \code{datacggm}.
}

\arguments{
  \item{x}{an object of class \code{datacggm}.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{datacggm}}, \code{\link{rcggm}} and \code{\link{dim.datacggm}}.
}

\examples{
set.seed(123)
n <- 100
p <- 3
q <- 2
b0 <- rep(1, p)
X <- matrix(rnorm(n * q), n, q)
B <- matrix(rnorm(q * p), q, p)
Sigma <- outer(1:p, 1:p, function(i, j) 0.3^abs(i - j))
probl <- 0.05
probr <- 0.05
probna <- 0.05

Z <- rcggm(n = n, b0 = b0, Sigma = Sigma, probl = probl, probr = probr, 
           probna = probna)
nObs(Z)
nResp(Z)
nPred(Z)

Z <- rcggm(b0 = b0, X = X, B = B, Sigma = Sigma, probl = probl, probr = probr, 
           probna = probna)
nObs(Z)
nResp(Z)
nPred(Z)
}

\keyword{methods}
\keyword{classes}
