\name{coef}
\alias{coef}
\alias{coef.glasso}
\alias{coef.cglasso}

\title{Extract Model Coefficients}

\description{
\sQuote{\code{coef}} extracts model coefficients from a fitted model.
}

\usage{
\method{coef}{glasso}(object, nrho = 1L, type = c("theta", "sigma"), print.info = FALSE,
     digits = 3L, \dots)
\method{coef}{cglasso}(object, nrho = 1L, type = c("theta", "sigma", "mu"), print.info = FALSE,
     digits = 3L, \dots)
}

\arguments{
\item{object}{an object of class \sQuote{\code{glasso}}, \sQuote{\code{ggm}}, \sQuote{\code{cglasso}} or \sQuote{\code{cggm}}.}
\item{nrho}{integer used to specify the model from which to extract the coefficients. Default is \code{nrho = 1}.}
\item{type}{a string specifying the returned parameters. If \sQuote{\code{object}} has class \sQuote{\code{glasso}}, the user can choice between the precision matrix (\sQuote{\code{type = "theta"}}) and the covariance matrix (\sQuote{\code{type = "sigma"}}). For an \sQuote{\code{object}} with class \sQuote{\code{cglasso}}, the user can also extract the estimates of the expected values (\sQuote{\code{type = "mu"}}). Default is \code{"theta"}}
\item{print.info}{flag specifying if information about the model is printed out. Default is \code{FALSE}.}
\item{digits}{the minimum number of significant digits to be used. Default is \code{3L}.}
\item{\dots}{additional argument added for backward compatibility with the generic function \link[stats]{coef}.}
}

\details{
By default, the method functions \sQuote{\code{coef.glasso}} and \sQuote{\code{coef.cglasso}} return the parameters specified by the argument \sQuote{\code{type}}.

If \sQuote{\code{print.info = TRUE}} then the estimated parameters are silently returned and information about the chosen model is printed out, i.e. the value of the tuning parameter, the value of the pseudo R-squared, the number of connected components and the number of vertices per connected component. Furthermore, to improve the readability of the results the estimates are printed out taken into account the connected components (see the examples below).
}

\value{
Coefficients extracted from \sQuote{\code{object}} are returned.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\seealso{
\code{\link{glasso}}, \code{\link{cglasso}} and \code{\link{mle}}.
}

\examples{
library("cglasso")

# cglasso model
set.seed(123)
p <- 5L
n <- 100L
mu <- rep(0L, p)
Tht <- diag(p)
diag(Tht[-1L, -p]) <- diag(Tht[-p, -1L]) <- 0.3
Sgm <- solve(Tht)
X <- rdatacggm(n = n, mu = mu, Sigma = Sgm, probr = 0.05)
out <- cglasso(X = X)

coef(out, nrho = 3L, type = "theta", print.info = TRUE)
Tht_hat <- coef(out, nrho = 3L, type = "theta")
Tht_hat

coef(out, nrho = 3L, type = "sigma", print.info  = TRUE)
Sgm_hat <- coef(out, nrho = 3L, type = "sigma")
Sgm_hat

coef(out, nrho = 3L, type = "mu", print.info = TRUE)
mu_hat <- coef(out, nrho = 3L, type = "mu")
mu_hat

# glasso model
X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sgm)
out <- glasso(X = X)

coef(out, nrho = 3L, type = "theta", print.info = TRUE)
Tht_hat <- coef(out, nrho = 3L, type = "theta")
Tht_hat

coef(out, nrho = 3L, type = "sigma", print.info  = TRUE)
Sgm_hat <- coef(out, nrho = 3L, type = "sigma")
Sgm_hat
}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}

