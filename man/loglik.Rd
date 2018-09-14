\name{loglik}
\alias{loglik}
\alias{print.loglik}

\title{Extract Log-Likelihood or Q-Function}

\description{
\sQuote{\code{loglik}} function extracts the values of the log-liklihood function from an object with class \sQuote{\code{glasso}} or \sQuote{\code{ggm}}, whereas the values of the Q-function are returned if the object has class \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}.
}

\usage{
loglik(object)
}

\arguments{
\item{object}{an object of class \sQuote{\code{glasso}}, \sQuote{\code{ggm}}, \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}.}
}

\details{
If \sQuote{\code{object}} has class \sQuote{\code{glasso}} or \sQuote{\code{ggm}}, the function \sQuote{\code{loglik()}} returns the value of the log-likelihood function: \deqn{\frac{n}{2}\left\{\log det\Theta - tr(S\Theta) - p\log(2\pi)\right\},}{n/2 \{log det Tht - tr(S Tht) - p log(2 pi)\},} where \eqn{\Theta}{Tht} is estimated using the function \link{glasso} or \link{mle.glasso}.

If \sQuote{\code{object}} has class \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}, the function \sQuote{\code{loglik()}} returns the value of the Q-function, i.e. the function maximaxed in the M-step of the EM-like algorithm (see details in \link{cglasso} or \link{mle.cglasso}). The Q-function is defined as follows: \deqn{\frac{n}{2}\left\{\log det\Theta - tr(S'\Theta) - p\log(2\pi)\right\},}{n/2 \{log det Tht - tr(S' Tht) - p log(2 pi)\},} where \eqn{S'}{S'} is computed in the E-step.

The method function \sQuote{\code{print.loglik}} is used the improve the readability of the results.
}

\value{
\sQuote{\code{loglik}} returns an object with S3 class \dQuote{\code{loglik}}, i.e. a list containing the following components:
\item{value}{the values of the log-likelihood function or the values of the Q-function.}
\item{df}{the number of the estimated non-zero parameters, i.e. the number of non-zero partial correlations plus \eqn{2p}{2p}.}
\item{n}{the sample size.}
\item{p}{the number of variables.}
\item{rho}{the values of the tuning parameter used to fit the model.}
\item{model}{the name of the fitted model.}
\item{fun}{the name of the used function, i.e. the log-likelihood or the Q-function.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{cglasso}, \code{\link{glasso}}, \code{\link{mle}} and the method functions, \code{\link{plot}}, \code{\link{aic}}, \code{\link{bic}} and \code{\link{ebic}}.
}

\examples{
library("cglasso")
set.seed(123)

#################
# cglasso model #
#################
p <- 5L
n <- 100L
mu <- rep(0L, p)
Tht <- diag(p)
diag(Tht[-1L, -p]) <- diag(Tht[-p, -1L]) <- 0.3
Sgm <- solve(Tht)
X <- rdatacggm(n = n, mu = mu, Sigma = Sgm, probr = 0.05)
out <- cglasso(X = X)
out_loglik <- loglik(out)
out_loglik

##############
# cggm model #
##############
out_mle <- mle(out)
out_loglik <- loglik(out_mle)
out_loglik

#################
# glasso model #
#################
library(MASS)
X <- mvrnorm(n = n, mu = mu, Sigma = Sgm)
out <- glasso(X)
out_loglik <- loglik(out)
out_loglik

#############
# ggm model #
#############
out_mle <- mle(out)
out_loglik <- loglik(out_mle)
out_loglik
}

\keyword{models}
