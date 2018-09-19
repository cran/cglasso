\name{mglasso}
\alias{mglasso}
\alias{mglasso.fit}

\title{Graphical Lasso Estimator with Missing-at-Random Data}

\description{
\sQuote{\code{mglasso}} function is used to fit an l1-penalized Gaussian graphical model with missing-at-random data.
}

\usage{
mglasso(X, weights, pendiag = FALSE, nrho = 50L, rho.min.ratio, rho,
        maxR2, maxit_em = 1.0e+3, thr_em = 1.0e-4, maxit_bcd = 1.0e+4,
        thr_bcd = 1.0e-4, trace = 0L)
}

\arguments{
\item{X}{the \eqn{(n\times p)}{(n x p)}-dimensional matrix used to fit the model.}
\item{weights}{an optional symmetric matrix of non-negative weights. This matrix can be used to specify the unpenalized partial correlation coefficients (\sQuote{\code{weights[i, j] = 0}}) or the structural zeros in the precision matrix (\sQuote{\code{weights[i, j] = +Inf}}). See below for an example. By default, mglasso model is fitted without weights.}
\item{pendiag}{flag used to specify if the diagonal elements of the concentration matrix are penalized (\sQuote{\code{pendiag = FALSE}}) or unpenalized (\sQuote{\code{pendiag = TRUE}}).}
\item{nrho}{the integer specifying the number of tuning parameters used to fit the mglasso model. Default is \sQuote{\code{nrho = 50}}.}
\item{rho.min.ratio}{the smallest value for the tuning parameter \eqn{\rho}{rho}, as a fraction of the smallest tuning parameter for which all the estimated partial correlation coefficients are zero. The default depends on the sample size \sQuote{\eqn{n}{n}} relative to the number of variables \sQuote{\eqn{p}{p}}. If \sQuote{\eqn{p < n}{p < n}}, the default is \sQuote{1.0E-4} otherwise the value \sQuote{1.0E-2} is used as default. A very small value of \sQuote{\code{rho.min.ratio}} will lead to a saturated fitted model in the \sQuote{\eqn{p < n}{p < n}} case.}
\item{rho}{optional argument. A user supplied \code{rho} sequence. WARNING: avoid supplying a single value for the tuning parameter; supply instead a decreasing sequence of \eqn{\rho}{rho}-values.}
\item{maxR2}{a value belonging to the interval \eqn{[0, 1]}{[0, 1]} specifying the largest value of the pseudo R-squared measure (see Section \bold{Details}). The regularization path is stopped when \eqn{R^2}{R2} exceeds \sQuote{\code{maxR2}}. Default depends on the sample size \sQuote{\eqn{n}{n}} relative to the number of variables \sQuote{\eqn{p}{p}}. If \sQuote{\eqn{p < n}{p < n}}, the default is \sQuote{1} otherwise the value \sQuote{0.9} is used as default.}
\item{maxit_em}{maximum number of iterations of the EM algorithm. Default is \code{1.0E+3}.}
\item{thr_em}{threshold for the convergence of the EM algorithm. Default value is \code{1.0E-4}.}
\item{maxit_bcd}{maximum number of iterations of the glasso algorithm. Default is \code{1.0E+4}.}
\item{thr_bcd}{threshold for the convergence of the glasso algorithm. Default is \code{1.0E-4}.}
\item{trace}{integer for printing out information as iterations proceed: \code{trace = 0} no information is printed out on video; \code{trace = 1} basic information is printed out on video; \code{trace = 2} detailed information is printed out on video.}
}

\details{
The missglasso estimator (Stadler \emph{and other}, 2012) is an extension of the classical graphical lasso (glasso) estimator (Yuan \emph{and other}, 2007) developed to fit a sparse  Gaussian graphical model under the assumption that data are missing-at-random.

\code{mglasso} function fits the model using the following EM algorithm:
\tabular{cl}{
Step \tab Description \cr
1.  \tab Let \eqn{\{\hat{\mu}^\rho_{ini};\hat{\Theta}^\rho_{ini}\}}{\{hat{mu}^{rho}_{ini}; hat{Tht}^{rho}_{ini}\}} be initial estimates;\cr
2.  \tab \bold{E-step}\cr
    \tab use the expected values of the conditional normal distribution to impute the missing data\cr
    \tab let \eqn{X^{\rho}}{X^{rho}} the completed data and \eqn{S^{\rho}}{S^{rho}} the corresponding empirical covariance matrix\cr
3.  \tab \bold{M-step}\cr
    \tab let \eqn{\hat{\mu}_h^\rho = \sum_{i=1}^n x_{ih}^\rho}{hat{mu}_h^{rho} = sum_i x^{rho}_{ih} / n}; \cr
    \tab compute \eqn{\hat{\Theta}^\rho}{hat{Tht}^{rho}} using \eqn{S^{\rho}}{S^{rho}} and the glasso algorithm (Friedman \emph{and other}, 2008);\cr
4.  \tab repeat steps 2. and 3. until a convergence criterion is met.
}

In order to avoid the overfitting of the model, we use the following pseudo R-squared measure: \deqn{R^2 = 1 - \frac{\|S^\rho - \hat{\Sigma}^\rho\|_F}{\|S^{\rho_{\max}} - \hat{\Sigma}^{\rho_{\max}}\|_F},}{R2 = 1 - ||S^rho - Sgm^rho||_F / ||S^{rho_max} - Sgm^{rho_max}||_F,} where \eqn{\|\cdot\|_F}{|| . ||_F} denotes the Frobenius norm and \eqn{\rho_{\max}}{rho_max} denotes the smallest value of the tuning parameter for which all the estimated partial correlation coefficients are zero. By straightforward algebra, it is easy to show that the proposed pseudo R-squared belongs to the closed interval \eqn{[0, 1]}{[0, 1]}: \eqn{R^2 = 0}{R2 = 0} when the tuning parameter is equal to \eqn{\rho_{\max}}{rho_max} and \eqn{R^2 = 1}{R2 = 1} when \eqn{\rho = 0}{rho = 0}. The regularization path is stopped when \eqn{R^2}{R2} exceeds the threshold specify by  \sQuote{\code{maxR2}}.
}

\value{
\code{mglasso} returns an object with S3 class \dQuote{\code{mglasso}}, i.e., a list containing the
following components:
\item{call}{the call that produced this object.}
\item{X}{the original matrix used to fit the missglasso model.}
\item{weights}{the weights used to fit the missglasso model.}
\item{pendiag}{the flag specifying if the diagonal elements of the precisione matrix are penalized.}
\item{nrho}{the number of fitted missglasso model.}
\item{rho.min.ratio}{the scale factor used to compute the smallest value of the tuning parameter.}
\item{rho}{the \eqn{p}{p}-dimensional vector reporting the values of the tuning parameter used to fit the missglasso model.}
\item{maxR2}{the threshold value used to stop the regularization path.}
\item{maxit_em}{the maximum number of iterations of the EM algorithm.}
\item{thr_em}{the threshold for the convergence of the EM algorithm.}
\item{maxit_bcd}{the maximum number of iterations of the glasso algorithm.}
\item{thr_bcd}{the threshold for the convergence of the glasso algorithm.}
\item{Xipt}{an array of dimension \eqn{n\times p\times\texttt{nrho}}{n x p x nrho}. \code{Xipt[, , k]} is the matrix where the missing values are replaced with the conditional expected values computed in the E-step of the algorithm describted in section \bold{Details}.}
\item{S}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{S[, , k]} is the matrix \eqn{S^\rho} used to fit the glasso model in the M-step of the algorithm describted in section \bold{Details}.}
\item{mu}{a matrix of dimension \eqn{p\times\texttt{nrho}}{p x nrho}. The \eqn{k}{k}th column is the estimate of the expected values of the missglasso model fitted using \code{rho[k]}.}
\item{Sgm}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Sgm[, , k]} is the estimate of the covariance matrix of the missglasso model fitted using \code{rho[k]}.}
\item{Tht}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Tht[, , k]} is the estimate of the precision matrix of the missglasso model fitted using \code{rho[k]}.}
\item{Adj}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Adj[, , k]} is the adjacency matrix associated to \code{Tht[, , k]}, i.e. \code{Adj[i, j, k] }\eqn{= 1}{= 1} iff \code{Tht[i, j, k] }\eqn{\neq 0}{!= 0} and \code{0} otherwise.}
\item{df}{the \eqn{p}{p}-dimensional vector reporting the number of non-zero partial correlation coefficients.}
\item{R2}{the \eqn{p}{p}-dimensional vector reporting the values of the measure \eqn{R^2}{R2} described in section \bold{Details}.}
\item{ncomp}{the \eqn{p}{p}-dimensional vector reporting the number of connected components (for internal purposes only).}
\item{Ck}{the \eqn{(p\times\texttt{nrho})}{(p x nrho)}-dimensional matrix encoding the connected components (for internal purposes only).}
\item{pk}{the \eqn{(p\times\texttt{nrho})}{(p x nrho)}-dimensional matrix reporting the number of vertices per connected component (for internal purposes only).}
\item{nit}{the \eqn{(\texttt{nrho}\times 2)}{(nrho x 2)}-dimensional matrix reporting the number of iterations.}
\item{conv}{a description of the error that has occurred.}
\item{subrout}{the name of the Fortran subroutine where the error has occurred (for internal debug only).}
\item{trace}{the integer used for printing out information.}
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Friedman, J.H., Hastie, T., and Tibshirani, R. (2008) <DOI:10.1093/biostatistics/kxm045>.
\href{https://academic.oup.com/biostatistics/article/9/3/432/224260}{Sparse inverse covariance estimation with the graphical lasso}.
\emph{Biostatistics} \bold{9}, 432--441.

Stadler N., and Buhlmann P. (2012) <DOI:10.1007/s11222-010-9219-7>.
\href{https://link.springer.com/article/10.1007/s11222-010-9219-7}{Missing values: sparse inverse covariance estimation and an extension to sparse regression}.
\emph{Statistics and Computing} \bold{22}, 219--235.

Yuan, M., and Lin, Y. (2007) <DOI:10.1093/biomet/asm018>.
\href{https://academic.oup.com/biomet/article-abstract/94/1/19/228821?redirectedFrom=fulltext}{Model selection and estimation in the Gaussian graphical model}.
\emph{Biometrika} \bold{94}, 19--35.
}

\seealso{
\code{\link{glasso}}, \code{\link{to_graph}}, \code{\link{mle}} and the method functions \code{\link{summary}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{aic}}, \code{\link{bic}}, \code{\link{ebic}}.
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
X <- as.vector(X)
id.na <- sample.int(n = n * p, size = n * p * 0.05)
X[id.na] <- NA
dim(X) <- c(n, p)
out <- mglasso(X = X)
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

out <- mglasso(X = X, weights = w)

# checking structural zeros
out$Tht[, , out$nrho][w == +Inf]

# checking stationarity conditions of the MLE estimators
# (the unpenalized partial correlation coefficients)
(out$Sgm[, , out$nrho] - out$S[, , out$nrho])[w == 0]
}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}

