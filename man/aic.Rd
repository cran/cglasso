\name{aic}
\alias{aic}
\alias{bic}
\alias{print.gof}
\alias{plot.gof}

\title{Akaike's An Information Criterion}

\description{
\sQuote{\code{aic}} computes the \sQuote{Akaike Information Criterion} whereas \sQuote{\code{bic}} computes the \sQuote{Bayesian Information Criterion}.
}

\usage{
aic(object, k = 2)

bic(object)
}

\arguments{
\item{object}{an object with class \sQuote{\code{glasso}}, \sQuote{\code{ggm}}, \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}.}
\item{k}{the \emph{penalty} per parameter to be used; the default \code{k = 2} is the classical AIC.}
}

\details{
The measure of goodness-of-fit (gof) returned by the functions \sQuote{\code{aic}} and \sQuote{\code{bic}} depends on the class of the fitted model.

If \sQuote{\code{object}} has class \sQuote{\code{glasso}} or \sQuote{\code{ggm}}, then \sQuote{\code{aic}} computes the following measure of goodness-of-fit: \deqn{-2\,\mbox{log-likelihood} + k\,\mbox{df},}{-2*log-likelihood + k*df,} where \eqn{k}{k} is the \emph{penalty} per parameter and \eqn{\mbox{df}}{df} represents the number of parameters in the fitted model. The values of the log-likelihood function are computed using the function \code{\link{loglik}}. The usual Akaike Information Criterion (AIC) is computed letting \eqn{k = 2} (default value of the function \sQuote{\code{aic}}) whereas the \sQuote{Bayesian Information Criterion} (BIC) is computed letting \eqn{k = \log(n)}{k = log(n)}, where \eqn{n} is the sample size.

If \sQuote{\code{object}} has class \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}, then \sQuote{\code{aic}} computes the following measure of goodness-of-fit: \deqn{-2\,Q\mbox{-function} + k\,df,}{-2*Q-function + k df,} in other words the log-likelihood is replaced with the \eqn{Q}{Q}-function maximized in the M-step of the EM-like algorithm describted in \code{\link{cglasso}} and \code{\link{mle}}. This measure of goodness-of-fit was proposed in Ibrahim \emph{and others} (2008) for statistical model with missing-data.

\sQuote{\code{aic}} and \sQuote{\code{bic}} return an object with S3 class \sQuote{\code{gof}} for which  are available the method functions \sQuote{\code{print.gof}} and \sQuote{\code{plot.gof}}. These method functions are developed with the aim of helping the user in finding the optimal value of the tuning parameter,  defined as the \eqn{\rho}{rho}-value minimizing the chosen measure of goodness-of-fit. For this reason, \sQuote{\code{print.gof}} shows also the ranking of the fitted models (the best model is pointed out with an arrow) whereas \sQuote{\code{plot.gof}} pointed out the optimal \eqn{\rho}{rho}-value by a vertical dashed line (see below for some examples).
}

\value{
\sQuote{\code{aic}} and \sQuote{\code{bic}} return an object with S3 class \dQuote{\code{gof}}, i.e. a list containing the following components:
\item{value_gof}{the values of the measure of goodness-of-fit used to evaluate the fitted models.}
\item{rho}{the values of the tuning parameter used to fit the model.}
\item{value}{the values of the log-likelihood function or the Q-function.}
\item{df}{the number of the estimated non-zero parameters, i.e. the number of non-zero partial correlations plus \eqn{2p}{2p}.}
\item{n}{the sample size.}
\item{p}{the number of variables.}
\item{model}{the name of the fitted models.}
\item{type}{the measure of goodness-of-fit used to evaluate the fitted models.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Ibrahim, J.G., Zhu, H. and Tang, N. (2008).
\href{https://amstat.tandfonline.com/doi/abs/10.1198/016214508000001057#.W4eujC1aYmI}{Model selection criteria for missing-data problems using the EM algorithm}.
\emph{Journal of the American Statistical Association} \bold{103}, 1648--1658.

Sakamoto, Y., Ishiguro, M., and Kitagawa, G. (1986).
\href{https://www.springer.com/gb/book/9789027722539}{\emph{Akaike Information Criterion Statistics}}.
D. Reidel Publishing Company.
}

\seealso{
\code{\link{loglik}}, \code{\link{cglasso}}, \code{\link{glasso}}, \code{\link{mle}}, \code{\link{ebic}} and the method funtions \sQuote{\link{plot}} and \code{\link{summary}}.
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
out_aic <- aic(out)
out_aic
plot(out_aic)

out_bic <- bic(out)
out_bic
plot(out_bic)

##############
# cggm model #
##############
out_mle <- mle(out)
out_aic <- aic(out_mle)
out_aic
plot(out_aic)

out_bic <- bic(out_mle)
out_bic
plot(out_bic)

#################
# glasso model #
#################
X <- MASS::mvrnorm(n = n, mu = mu, Sigma = diag(p))
out <- glasso(X)
out_aic <- aic(out)
out_aic
plot(out_aic)

out_bic <- bic(out)
out_bic
plot(out_bic)

#############
# ggm model #
#############
out_mle <- mle(out)
out_aic <- aic(out_mle)
out_aic
plot(out_aic)

out_bic <- bic(out_mle)
out_bic
plot(out_bic)
}

\keyword{models}
