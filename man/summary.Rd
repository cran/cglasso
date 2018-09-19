\name{summary}
\alias{summary}
\alias{summary.glasso}

\title{Summary Method}

\description{
\sQuote{\code{summary}} produces a summary of the sequence of fitted models.
}

\usage{
\method{summary}{glasso}(object, \dots, gof = c("BIC", "AIC", "eBIC"), par.gof, digits = 4L)
}

\arguments{
\item{object}{an object of class \sQuote{\code{glasso}}, \sQuote{\code{ggm}}, \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}.}
\item{gof}{string specifying the measure of goodness-of-fit used to evaluate the fitted models. Default is \sQuote{\code{BIC}}.}
\item{par.gof}{the parameter of the measure of goodness-of-fit used to evaluate the fitted models.}
\item{digits}{the minimum number of significant digits to be used: see \sQuote{\link[base]{print.default}}.}
\item{\dots}{further arguments passed to the method function \sQuote{\link[base]{print.data.frame}}.}
}

\details{
The method function \code{summary.glasso} gives information about the sequence of fitted models. The output is divided in two sections.

First section shows the call that produced \code{object} followed by a \code{data.frame} reporting the values of the tuning parameter used to fit the model (\code{rho}), the number of non-zero estimates (\code{df}), the values of the pseudo R-squared (\code{R2}) describted in \code{\link{glasso}}, \code{\link{mglasso}} and \code{\link{cglasso}}, the values of the measure of goodness-of-fit used to evaluate the fitted models and the ranking of the fitted models (\code{Rank}). The model with the lowest measure of goodness-of-fit is pointed out by an arrow.

Second section shows the details of the selected model plus the number of connected components and the number of vertices per component.
}

\value{
The function \sQuote{\code{summary.glasso}} computes and returns a list of summary statistics with the following elements:
\item{table}{a data.frame containing the summary statistics used to evaluate the sequence of fitted models.}
\item{which.min}{the number of the model with the lowest measure of goodness-of-fit.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{link{glasso}}, \code{\link{mglasso}}, \code{\link{cglasso}}, \code{\link{mle}}, \code{\link{aic}} \code{\link{bic}} and \code{\link{ebic}}.
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
summary(out, gof = "AIC")
summary(out, gof = "BIC")
summary(out, gof = "eBIC")

##############
# cggm model #
##############
out_mle <- mle(out)
summary(out_mle, gof = "AIC")
summary(out_mle, gof = "BIC")
summary(out_mle, gof = "eBIC")

#################
# cglasso model #
#################
R <- event(X)
X <- as.matrix(X)
X[R == 1L] <- NA
out <- mglasso(X = X)
summary(out, gof = "AIC")
summary(out, gof = "BIC")
summary(out, gof = "eBIC")

##############
# mggm model #
##############
out_mle <- mle(out)
summary(out_mle, gof = "AIC")
summary(out_mle, gof = "BIC")
summary(out_mle, gof = "eBIC")

#################
# glasso model #
#################
X <- MASS::mvrnorm(n = n, mu = mu, Sigma = diag(p))
out <- glasso(X)
summary(out, gof = "AIC")
summary(out, gof = "BIC")
summary(out, gof = "eBIC")

#############
# ggm model #
#############
out_mle <- mle(out)
summary(out_mle, gof = "AIC")
summary(out_mle, gof = "BIC")
summary(out_mle, gof = "eBIC")
}

\keyword{multivariate}
\keyword{models}

