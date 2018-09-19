\name{ebic}
\alias{ebic}
\alias{ebic.glasso}
\alias{ebic.mglasso}
\alias{ebic.cglasso}

\title{Extended Bayesian Information Criterion}

\description{
\sQuote{\code{ebic}} function computes the extended Bayesian Information Criterion.
}

\usage{
ebic(object, g)

\method{ebic}{glasso}(object, g = 0.5)

\method{ebic}{mglasso}(object, g = 0.5)

\method{ebic}{cglasso}(object, g = 0.5)
}

\arguments{
\item{object}{a fitted model object.}
\item{g}{the parameter indexing the extended BIC: a value belonging to the interval \eqn{[0, 1]}{[0, 1]}. Default is \code{0.5}}
}

\details{
The measure of goodness-of-fit (gof) returned by the function \sQuote{\code{ebic}} depends on the class of the fitted model.

If \sQuote{\code{object}} has class \sQuote{\code{glasso}} or \sQuote{\code{ggm}}, then \sQuote{\code{ebic}} computes the extended Bayesian Information Criterion (eBIC) proposed in Foygel \emph{and others} (2010): \deqn{\mbox{eBIC} = -2\,\mbox{log-likelihood} + a(\rho)(\log n + 4 \gamma log p),}{eBIC = -2*log-likelihood + a(rho)*(log(n) + 4*g*log(p)),} where \eqn{a(\rho)}{a(rho)} denotes the number of non-zero off-diagonal elements in \eqn{\hat{\Theta}^\rho}{hat{Tht}^rho} and \eqn{\gamma}{g} is a value belonging to the interval \eqn{[0, 1]}{[0, 1]} indexing the measure of goodness-of-fit. As explained in Foygel \emph{and others} (2010), the log-likelihood function is evaluated using the maximum likelihood estimates of the model select by glasso. For this reason, \sQuote{\code{ebic}} calls the generic function \code{\link{mle}} to fit the Gaussian graphical model (GGM) selected by \code{\link{glasso}}.

For the remaining models, eBIC is defined as: \deqn{\mbox{eBIC} = -2\,Q\mbox{-function} + a(\rho)(\log n + 4 \gamma log p),}{eBIC = -2*Q-function + a(rho)*(log(n) + 4*g*log(p)),} where the \eqn{Q}{Q}-function is evaluated at the M-step of the EM-like algorithm describted in \code{\link{mle}}.

\sQuote{\code{ebic}} returns an object with S3 class \sQuote{\code{gof}} for which are available the method functions \sQuote{\code{print.gof}} and \sQuote{\code{plot.gof}}. These method functions are developed with the aim of helping the user in finding the optimal value of the tuning parameter, defined as the \eqn{\rho}{rho}-value minimizing the eBIC meaure. For this reason, \sQuote{\code{print.gof}} shows also the ranking of the fitted models (the best model is pointed out with an arrow) whereas \sQuote{\code{plot.gof}} points out the optimal \eqn{\rho}{rho}-value by a vertical dashed line (see below for some examples).
}

\value{
\sQuote{\code{ebic}} returns an object with S3 class \dQuote{\code{gof}}, i.e. a list containing the following components:
\item{value_gof}{the values of the measure of goodness-of-fit used to evaluate the fitted models.}
\item{rho}{the values of the tuning parameter used to fit the models.}
\item{value}{the values of the log-likelihood or Q-function.}
\item{df}{the number of the estimated non-zero parameters, i.e. the number of non-zero partial correlations plus \eqn{2p}{2p}.}
\item{n}{the sample size.}
\item{p}{the number of variables.}
\item{model}{the name of the fitted model.}
\item{type}{the measure of goodness-of-fit used to evaluate the fitted models.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Foygel, R. and Drton, M. (2010).
\href{https://papers.nips.cc/paper/4087-extended-bayesian-information-criteria-for-gaussian-graphical-models}{Extended Bayesian Information Criteria for Gaussian Graphical Models}.
In: Lafferty, J., Williams, C., Shawe-taylor, J., Zemel, R.s. and Culott, A. (editors),
\emph{Advances in Neural Information Processing Systems 23}. pp. 604--612.
}

\seealso{
\code{\link{loglik}}, \code{\link{cglasso}}, \code{\link{mglasso}}, \code{\link{glasso}}, \code{\link{mle}}, \code{\link{aic}}, \code{\link{bic}} and the method funtions \code{\link{plot}} and \code{\link{summary}}.
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
out <- cglasso(X = X)
out_ebic <- ebic(out)
out_ebic
plot(out_ebic)

##############
# cggm model #
##############
out_mle <- mle(out)
out_ebic <- ebic(out_mle)
out_ebic
plot(out_ebic)

#################
# mglasso model #
#################
X <- rnorm(n * p)
id.na <- sample.int(n = n * p, size = n * p * 0.05)
X[id.na] <- NA
dim(X) <- c(n, p)
out <- mglasso(X = X)
out_ebic <- ebic(out)
out_ebic
plot(out_ebic)

##############
# mggm model #
##############
out_mle <- mle(out)
out_ebic <- ebic(out_mle)
out_ebic
plot(out_ebic)

#################
# glasso model #
#################
X <- rnorm(n * p)
dim(X) <- c(n, p)
out <- glasso(X)
out_ebic <- ebic(out)
out_ebic
plot(out_ebic)

#############
# ggm model #
#############
out_mle <- mle(out)
out_ebic <- ebic(out_mle)
out_ebic
plot(out_ebic)
}

\keyword{models}
