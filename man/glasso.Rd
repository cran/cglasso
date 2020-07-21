\name{glasso}
\alias{glasso}
\alias{glasso.fit}
\alias{print.glasso}

\title{Lasso Estimator for Gaussian Graphical Models}

\description{
\sQuote{\code{glasso}} fits the l1-penalized Gaussian graphical model.
}

\usage{
glasso(X, weights, pendiag = FALSE, nrho = 50L, rho.min.ratio, rho, maxR2,
       maxit = 1.0e+4, thr = 1.0e-04, trace = 0L)
}

\arguments{
\item{X}{the \eqn{(n\times p)}{(n x p)}-dimensional matrix used to compute the covariance matrix.}
\item{weights}{an optional symmetric matrix of non-negative weights. This matrix can be used to specify the unpenalized partial correlation coefficients (\sQuote{\code{weights[i, j] = 0}}) or the structural zeros in the precision matrix (\sQuote{\code{weights[i, j] = +Inf}}). See below for an example. By default, glasso model is fitted without weights.}
\item{pendiag}{flag used to specify if the diagonal elements of the concentration matrix are penalized (\sQuote{\code{pendiag = TRUE}}) or unpenalized (\sQuote{\code{pendiag = FALSE}}).}
\item{nrho}{the integer specifying the number of tuning parameters used to fit the glasso model. Default is \sQuote{\code{nrho = 50L}}.}
\item{rho.min.ratio}{the smallest value for the tuning parameter \eqn{\rho}{rho}, as a fraction of the smallest tuning parameter for which all the estimated partial correlation coefficients are zero. The default depends on the sample size \sQuote{\eqn{n}{n}} relative to the number of variables \sQuote{\eqn{p}{p}}. If \sQuote{\eqn{p < n}{p < n}}, the default is \sQuote{1.0E-4} otherwise the value \sQuote{1.0E-2} is used as default. A very small value of \sQuote{\code{rho.min.ratio}} will lead to a saturated fitted model in the \sQuote{\eqn{p < n}{p < n}} case.}
\item{rho}{optional argument. A user supplied \code{rho} sequence. WARNING: avoid supplying a single value for the tuning parameter; supply instead a decreasing sequence of \eqn{\rho}{rho}-values.}
\item{maxR2}{a value belonging to the interval \eqn{[0, 1]}{[0, 1]} specifying the largest value of the pseudo R-squared measure (see Section \bold{Details}). The regularization path is stopped when \eqn{R^2}{R2} exceeds \sQuote{\code{maxR2}}. The default depends on the sample size \sQuote{\eqn{n}{n}} relative to the number of variables \sQuote{\eqn{p}{p}}. If \sQuote{\eqn{p < n}{p < n}}, the default is \sQuote{1} otherwise the value \sQuote{0.9} is used as default.}
\item{maxit}{maximum number of iterations of the glasso algorithm. Default is \code{1.0E+4}.}
\item{thr}{threshold for the convergence of the glasso algorithm. Default is \code{1.0E-4}.}
\item{trace}{integer for printing out information as iterations proceed: \code{trace = 0} no information is printed out on video; \code{trace = 1} basic information is printed out on video; \code{trace = 2} detailed information is printed out on video.}
}

\details{
For a fixed value of the tuning parameter, \code{glasso} solves the following maximization problem: \deqn{\max_{\Theta}\log det\Theta - tr\{S\Theta\}-\rho \sum_{h,k}w_{hk}|\theta_{hk}|,}{max_{Tht} log det Tht - tr\{S Tht\} - rho sum_{hk} w_{hk} |tht_{hk}|,} where \eqn{w_{hk}}{w_{hk}} is the non-negative weight for \eqn{\theta_{hk}}{tht_{hk}}. The previous maximization problem is solved effeciently combining the block-coordinate descent algorithm (Friedman \emph{and others}, 2008) with the screening rule proposed in Witten \emph{and others} (2011).

In order to avoid the overfitting of the model, we use the following pseudo R-squared measure: \deqn{R^2 = 1 - \frac{\|S - \hat{\Sigma}^\rho\|_F}{\|S - \hat{\Sigma}^{\rho_{\max}}\|_F},}{R2 = 1 - ||S - Sgm^rho||_F / ||S - Sgm^{rho_max}||_F,} where \eqn{\|\cdot\|_F}{|| . ||_F} denotes the Frobenius norm and \eqn{\rho_{\max}}{rho_max} denotes the smallest value of the tuning parameter for which all the estimated partial correlation coefficients are zero. By straightforward algebra, it is easy to show that the proposed pseudo R-squared belongs to the closed interval \eqn{[0, 1]}{[0, 1]}: \eqn{R^2 = 0}{R2 = 0} when the tuning parameter is equal to \eqn{\rho_{\max}}{rho_max} and \eqn{R^2 = 1}{R2 = 1} when \eqn{\rho = 0}{rho = 0}. The regularization path is stopped when \eqn{R^2}{R2} exceeds \sQuote{\code{maxR2}}.
}

\value{
\sQuote{\code{glasso}} returns an object with S3 class \dQuote{\code{glasso}}, i.e. a list containing the following components:
\item{call}{the call that produced this object.}
\item{X}{the matrix used to compute the covariance matrix.}
\item{S}{the covariance matrix used to fit the glasso model.}
\item{weights}{the used weights.}
\item{pendiag}{the flag specifying if the diagonal elements of the precisione matrix are penalized.}
\item{nrho}{the number of fitted glasso model.}
\item{rho.min.ratio}{the scale factor used to compute the smallest value of the tuning parameter.}
\item{rho}{the \eqn{p}{p}-dimensional vector reporting the values of the tuning parameter used to fit the glasso model.}
\item{maxR2}{the threshold value used for the pseudo R-squared measure.}
\item{maxit}{the maximum number of iterations of the glasso algorithm.}
\item{thr}{the threshold for the convergence of the glasso algorithm.}
\item{Sgm}{an array of dimension \eqn{(p\times p\times\texttt{nrho})}{p x p x nrho)}. \code{Sgm[, , k]} is the estimate of the covariance matrix of the glasso model fitted using \code{rho[k]}.}
\item{Tht}{an array of dimension \eqn{(p\times p\times\texttt{nrho})}{(p x p x nrho)}. \code{Tht[, , k]} is the estimate of the precision matrix of the glasso model fitted using \code{rho[k]}.}
\item{Adj}{an array of dimension \eqn{(p\times p\times\texttt{nrho})}{(p x p x nrho)}. \code{Adj[, , k]} is the adjacency matrix associated to \code{Tht[, , k]}, i.e. \code{Adj[i, j, k] }\eqn{= 1}{= 1} iff \code{Tht[i, j, k] }\eqn{\neq 0}{!= 0} and \code{0} otherwise.}
\item{df}{the \eqn{nrho}{nrho}-dimensional vector reporting the number of non-zero partial correlation coefficients.}
\item{R2}{the \eqn{nrho}{nrho}-dimensional vector reporting the values of the measure \eqn{R^2}{R2} described in the section \bold{Details}.}
\item{ncomp}{the \eqn{nrho}{nrho}-dimensional vector reporting the number of connected components (for internal purposes only).}
\item{Ck}{the \eqn{(p\times\texttt{nrho})}{(p x nrho)}-dimensional matrix encoding the connected components (for internal purposes only).}
\item{pk}{the \eqn{(p\times\texttt{nrho})}{(p x nrho)}-dimensional matrix reporting the number of vertices per connected component (for internal purposes only).}
\item{nit}{the \eqn{p}{p}-dimensional vector reporting the number of iterations.}
\item{conv}{a description of the error that has occurred.}
\item{trace}{the integer used for printing out information.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Friedman, J.H., Hastie, T., and Tibshirani, R. (2008).
Sparse inverse covariance estimation with the graphical lasso.
\emph{Biostatistics} \bold{9}, 432--441.

Witten, D.M., Friedman, J.H, and Simon, N. (2011).
New insights and faster computations for the graphical lasso.
\emph{Journal of Computational and Graphical Statistics} \bold{20}, 892--900.
}

\seealso{
\code{\link{mle}}, \code{\link{to_graph}} and the method functions \code{\link{summary}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{aic}}, \code{\link{bic}} and \code{\link{ebic}}.
}

\examples{
library("cglasso")
set.seed(123)

p <- 5L
n <- 100L
mu <- rep(0L, p)
Tht <- diag(p)
diag(Tht[-1L, -p]) <- diag(Tht[-p, -1L]) <- 0.3
Sgm <- solve(Tht)
X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sgm)
out <- glasso(X)
out

# in this example we use the argument 'weights' to specify
# the unpenalized partial correlation coefficients and the
# structural zeros in the precision matrix

w <- rep(1, p * p)
dim(w) <- c(p, p)

# specifing the unpenalized partial correlation coefficients
diag(w) <- diag(w[-1L, -p]) <- diag(w[-p, -1L]) <- 0

# specifing the structural zero
w[1L, 4L:5L] <- w[4L:5L, 1L] <- +Inf
w[2L, 5L] <- w[5L, 2L] <- +Inf
w

out <- glasso(X = X, weights = w)

# checking structural zeros
out$Tht[, , out$nrho][w == +Inf]

# checking stationarity conditions of the MLE estimators
# (the unpenalized partial correlation coefficients)
(out$Sgm[, , out$nrho] - out$S)[w == 0L]
}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}
