\name{mle}
\alias{mle}
\alias{mle.glasso}
\alias{mle.mglasso}
\alias{mle.cglasso}

\title{Maximum Likelihood Estimation}

\description{
The generic function \sQuote{\code{mle}} fits the graphical model selected by \code{\link{glasso}}, \code{\link{mglasso}} or \code{\link{cglasso}}.
}

\usage{
mle(object, \dots)

\method{mle}{glasso}(object, \dots, maxit = 1.0e+04, thr = 1.0e-04, trace = 0L)

\method{mle}{mglasso}(object, \dots, maxit_em = 1.0e+03, thr_em = 1.0e-4, maxit_bcd = 1.0e+4,
    thr_bcd = 1.0e-4, trace = 0L)

\method{mle}{cglasso}(object, \dots, maxit_em = 1.0e+03, thr_em = 1.0e-4, maxit_bcd = 1.0e+4,
    thr_bcd = 1.0e-4, trace = 0L)
}

\arguments{
\item{object}{an object of class \sQuote{\code{glasso}}, \sQuote{\code{mglasso}} or \sQuote{\code{cglasso}}.}
\item{maxit}{maximum number of iterations of the glasso algorithm. Default is \code{1.0E+4}.}
\item{thr}{threshold for the convergence of the glasso algorithm. Default is \code{1.0E-4}.}
\item{maxit_em}{maximum number of iterations of the EM algorithm. Default is \code{1.0E+03}.}
\item{thr_em}{threshold for the convergence of the EM algorithm. Default is \code{1.0E-4}.}
\item{maxit_bcd}{maximum number of iterations of the glasso algorithm. Default is \code{1.0E+4}.}
\item{thr_bcd}{threshold for the convergence of the glasso algorithm. Default is \code{1.0E-4}.}
\item{trace}{integer for printing out information as iterations proceed: \code{trace = 0} no information is printed out on video; \code{trace = 1} basic information is printed out on video; \code{trace = 2} detailed information is printed out on video.}
\item{\dots}{additional argument added for backward compatibility with the generic function \sQuote{\code{mle}}.}
}

\details{
The generic function \sQuote{\code{mle}} computes the maximum likelihood estimates of the graphical model selected by the function \code{\link{glasso}}, \code{\link{mglasso}} or \code{\link{cglasso}}.

If \sQuote{\code{object}} has class \sQuote{\code{glasso}}, the method function \sQuote{\code{mle.glasso}} computes the maximum likelihood estimates of the parameters of the Gaussian graphical models (ggm) associated to the sequence of glasso estimates. Formally, for a given value of the tuning parameter let \eqn{\hat{\Theta}^\rho}{\hat\{|Tht\}^rho} be the glasso estimate of the precision matrix, then \sQuote{\code{mle.glasso}} solves the following maximization problem: \deqn{\max_{\bar{\Theta}}\log det\bar{\Theta} - tr\{S\bar{\Theta}\},}{max_{bar{Tht}} log det bar{Tht} - tr\{S bar{Tht}\},} where \eqn{\bar{\Theta} = \{\bar{\theta}_{hk}\}}{bar{Tht} = {tht_{hk}}} and \eqn{\bar{\theta}_{hk} = 0}{bar\{tht\}_{hk} = 0} if \eqn{\hat{\theta}^{\rho}_{hk} = 0}{hat\{tht\}^rho_{hk} = 0} otherwise it is estimated.

If \sQuote{\code{object}} has class \sQuote{\code{mglasso}}, the method function \sQuote{\code{mle.mglasso}} computes the maximum likelihood estimates of the parameters of the Gaussian graphical models with missing-at-random data (mggm) associated to the sequence of mglasso estimates. Formally, for a given value of the tuning parameter let \eqn{\hat{\Theta}^\rho}{hat{Tht}^rho} be the mglasso estimate of the precision matrix, then \sQuote{\code{mle.mglasso}} computes the maximum likelihood estimate by the following EM-like algorithm:
    \tabular{cl}{
    Step \tab Description \cr
    1.  \tab let \eqn{\hat{\Theta}^\rho}{hat{Tht}^rho} be the mglasso estimate of the precision matrix; \cr
    2.  \tab \bold{E-step}\cr
        \tab use the moments of the conditional normal distribution to impute the missing values; \cr
    3.  \tab \bold{M-step}\cr
        \tab let \eqn{X'}{X'} the completed matrix and \eqn{S'}{S'} the corresponding empirical variance matrix, then:\cr
        \tab let \eqn{\hat{\mu}_h = \sum_{i = 1}^n x'_{ih} / n}{hat{mu}_h = sum_i x'_{ih} / n}\cr
        \tab estimate \eqn{\bar{\Theta}}{bar{Tht}} maximixing the function: \eqn{\log det\bar{\Theta} - tr\{S'\bar{\Theta}\},}{log det bar{Tht} - Tr{S'bar{Tht}}} where \eqn{\bar{\theta}_{hk} = 0}{bar{tht}_{hk} = 0} if \eqn{\hat{\theta}^{\rho}_{hk} = 0}{hat{tht}^rho_{hk} = 0} \cr
        \tab otherwise it is estimated; \cr
    4.  \tab repeat steps 2. and 3. until a convergence criterion is met.
    }

If \sQuote{\code{object}} has class \sQuote{\code{cglasso}}, the method function \sQuote{\code{mle.cglasso}} computes the maximum likelihood estimates of the parameters of the censored Gaussian graphical models (cggm) associated to the sequence of cglasso estimates. Formally, for a given value of the tuning parameter let \eqn{\hat{\Theta}^\rho}{hat{Tht}^rho} be the cglasso estimate of the precision matrix, then \sQuote{\code{mle.cglasso}} computes the maximum likelihood estimate by the following EM-like algorithm:
    \tabular{cl}{
    Step \tab Description \cr
    1.  \tab let \eqn{\hat{\Theta}^\rho}{hat{Tht}^rho} be the cglasso estimate of the precision matrix; \cr
    2.  \tab \bold{E-step}\cr
        \tab use the moments of the truncated normal distribution to compute the current estimates of the \cr
        \tab marginal means, denoted by \eqn{\bar{x}'}{bar{x}'}, and to complete the empirical covariance matrix \eqn{S'}{S'};\cr
    3.  \tab \bold{M-step}\cr
        \tab let \eqn{\hat{\mu} = \bar{x}'}{hat{mu} = bar{x}'}; \cr
        \tab estimate \eqn{\bar{\Theta}}{bar{Tht}} maximixing the \eqn{Q}{Q}-function: \eqn{\log det\bar{\Theta} - tr\{S'\bar{\Theta}\},}{log det bar{Tht} - Tr{S'bar{Tht}}} where \eqn{\bar{\theta}_{hk} = 0}{bar{tht}_{hk} = 0} if \eqn{\hat{\theta}^{\rho}_{hk} = 0}{hat{tht}^rho_{hk} = 0} \cr
        \tab otherwise it is estimated; \cr
    4.  \tab repeat steps 2. and 3. until a convergence criterion is met.
    }
}

\value{
If \sQuote{\code{object}} has class \sQuote{\code{glasso}}, then \sQuote{\code{mle}} returns and object with S3 class \sQuote{\code{ggm}}, which inherits the class \sQuote{\code{glasso}}. An object with class \sQuote{\code{ggm}} is a list containing the following components:

\item{call}{the call that produced this object.}
\item{X}{the matrix used to compute the covariance matrix.}
\item{S}{the covariance matrix used to fit the ggm model.}
\item{nrho}{the number of fitted ggm model.}
\item{rho}{the \eqn{p}{p}-dimensional vector reporting the values of the tuning parameter used to fit the glasso model.}
\item{maxit}{the maximum number of iterations of the glasso algorithm.}
\item{thr}{the threshold for the convergence of the glasso algorithm.}
\item{Sgm}{an array of dimension \eqn{(p\times p\times\texttt{nrho})}{p x p x nrho)}. \code{Sgm[, , k]} is the estimate of the covariance matrix of the \eqn{k}{k}th ggm model.}
\item{Tht}{an array of dimension \eqn{(p\times p\times\texttt{nrho})}{(p x p x nrho)}. \code{Tht[, , k]} is the estimate of the precision matrix of the \eqn{k}{k}th ggm model.}
\item{Adj}{an array of dimension \eqn{(p\times p\times\texttt{nrho})}{(p x p x nrho)}. \code{Adj[, , k]} is the adjacency matrix associated to \code{Tht[, , k]}, i.e. \code{Adj[i, j, k] }\eqn{= 1}{= 1} iff \code{Tht[i, j, k] }\eqn{\neq 0}{!= 0} and \code{0} otherwise.}
\item{df}{the \eqn{p}{p}-dimensional vector reporting the number of non-zero partial correlation coefficients.}
\item{R2}{the \eqn{p}{p}-dimensional vector reporting the values of the measure \eqn{R^2}{R2}. See section \bold{Details}, in \code{\link{glasso}}.}
\item{ncomp}{the \eqn{p}{p}-dimensional vector reporting the number of connected components (for internal purposes only).}
\item{Ck}{the \eqn{(p\times\texttt{nrho})}{(p x nrho)}-dimensional matrix encoding the connected components (for internal purposes only).}
\item{pk}{the \eqn{(p\times\texttt{nrho})}{(p x nrho)}-dimensional matrix reporting the number of vertices per connected component (for internal purposes only).}
\item{nit}{the \eqn{p}{p}-dimensional vector reporting the number of iterations.}
\item{conv}{a description of the error that has occurred.}
\item{trace}{the integer used for printing out information.}

If \sQuote{\code{object}} has class \sQuote{\code{mglasso}}, then \sQuote{\code{mle}} returns and object with S3 class \sQuote{\code{mggm}}, which inherits the class \sQuote{\code{mglasso}}. An object with class \sQuote{\code{mggm}} is a list containing the following components:
\item{call}{the call that produced this object.}
\item{X}{the object with S3 class \sQuote{\code{datacggm}} used to fit the cggm model.}
\item{nrho}{the number of fitted cggm model.}
\item{rho}{the \eqn{p}{p}-dimensional vector reporting the values of the tuning parameter used to fit the cglasso model.}
\item{maxit_em}{maximum number of iterations of the EM algorithm.}
\item{thr_em}{threshold for the convergence of the EM algorithm.}
\item{maxit_bcd}{maximum number of iterations of the glasso algorithm.}
\item{thr_bcd}{threshold for the convergence of the glasso algorithm.}
\item{Xipt}{an array of dimension \eqn{n\times p\times\texttt{nrho}}{n x p x nrho}. \code{Xipt[, , k]} is the matrix where the censored vaules are replaced using the conditional expected vaules computed in the E-step of the algorithm describted in section \bold{Details}.}
\item{S}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{S[, , k]} is the matrix \eqn{S'}{S'} used to fit the cggm model (see the section \bold{Details}).}
\item{mu}{a matrix of dimension \eqn{p\times\texttt{nrho}}{p x nrho}. The \eqn{k}{k}th column is the estimate of the expected values of the \eqn{k}{k}th cggm model.}
\item{Sgm}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Sgm[, , k]} is the estimate of the covariance matrix of the \eqn{k}{k}th cggm model.}
\item{Tht}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Tht[, , k]} is the estimate of the precision matrix of the \eqn{k}{k}th cggm model.}
\item{Adj}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Adj[, , k]} is the adjacency matrix associated to \code{Tht[, , k]}, i.e. \code{Adj[i, j, k] }\eqn{= 1}{= 1} iff \code{Tht[i, j, k] }\eqn{\neq 0}{!= 0} and \code{0} otherwise.}
\item{df}{the \eqn{p}{p}-dimensional vector reporting the number of non-zero partial correlation coefficients.}
\item{R2}{the \eqn{p}{p}-dimensional vector reporting the values of the measure \eqn{R^2}{R2}. See section \bold{Details}, in \code{\link{cglasso}}.}
\item{ncomp}{the \eqn{p}{p}-dimensional vector reporting the number of connected components (for internal purposes only).}
\item{Ck}{the \eqn{(p\times\texttt{nrho})}{(p x nrho)}-dimensional matrix encoding the connected components (for internal purposes only).}
\item{pk}{the \eqn{(p\times\texttt{nrho})}{(p x nrho)}-dimensional matrix reporting the number of vertices per connected component (for internal purposes only).}
\item{nit}{the \eqn{(\texttt{nrho}\times 2)}{(nrho x 2)}-dimensional matrix reporting the number of iterations.}
\item{conv}{a description of the error that has occurred.}
\item{subrout}{the name of the Fortran subroutine where the error has occurred (for internal debug only).}
\item{trace}{the integer used for printing out information.}

If \sQuote{\code{object}} has class \sQuote{\code{cglasso}}, then \sQuote{\code{mle}} returns and object with S3 class \sQuote{\code{cggm}}, which inherits the class \sQuote{\code{cglasso}}. An object with class \sQuote{\code{cggm}} is a list containing the following components:
\item{call}{the call that produced this object.}
\item{X}{the object with S3 class \sQuote{\code{datacggm}} used to fit the cggm model.}
\item{nrho}{the number of fitted cggm model.}
\item{rho}{the \eqn{p}{p}-dimensional vector reporting the values of the tuning parameter used to fit the cglasso model.}
\item{maxit_em}{maximum number of iterations of the EM algorithm.}
\item{thr_em}{threshold for the convergence of the EM algorithm.}
\item{maxit_bcd}{maximum number of iterations of the glasso algorithm.}
\item{thr_bcd}{threshold for the convergence of the glasso algorithm.}
\item{Xipt}{an array of dimension \eqn{n\times p\times\texttt{nrho}}{n x p x nrho}. \code{Xipt[, , k]} is the matrix where the censored vaules are replaced using the conditional expected vaules computed in the E-step of the algorithm describted in section \bold{Details}.}
\item{S}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{S[, , k]} is the matrix \eqn{S'}{S'} used to fit the cggm model (see the section \bold{Details}).}
\item{mu}{a matrix of dimension \eqn{p\times\texttt{nrho}}{p x nrho}. The \eqn{k}{k}th column is the estimate of the expected values of the \eqn{k}{k}th cggm model.}
\item{Sgm}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Sgm[, , k]} is the estimate of the covariance matrix of the \eqn{k}{k}th cggm model.}
\item{Tht}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Tht[, , k]} is the estimate of the precision matrix of the \eqn{k}{k}th cggm model.}
\item{Adj}{an array of dimension \eqn{p\times p\times\texttt{nrho}}{p x p x nrho}. \code{Adj[, , k]} is the adjacency matrix associated to \code{Tht[, , k]}, i.e. \code{Adj[i, j, k] }\eqn{= 1}{= 1} iff \code{Tht[i, j, k] }\eqn{\neq 0}{!= 0} and \code{0} otherwise.}
\item{df}{the \eqn{p}{p}-dimensional vector reporting the number of non-zero partial correlation coefficients.}
\item{R2}{the \eqn{p}{p}-dimensional vector reporting the values of the measure \eqn{R^2}{R2}. See section \bold{Details}, in \code{\link{cglasso}}.}
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
Augugliaro, L., Abbruzzo, A., and Vinciotti, V. (2018).
\eqn{\ell_1}{l1}-Penalized gaussian graphical model.
\emph{Biostatistics} (to appear).

Stadler N., and Buhlmann P. (2012) <DOI:10.1007/s11222-010-9219-7>.
\href{https://link.springer.com/article/10.1007/s11222-010-9219-7}{Missing values: sparse inverse covariance estimation and an extension to sparse regression}.
\emph{Statistics and Computing} \bold{22}, 219--235.
}

\seealso{
\code{\link{glasso}}, \code{\link{mglasso}}, \code{\link{cglasso}}, \code{\link{to_graph}}, and the method functions \code{\link{summary}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{aic}}, \code{\link{bic}} and \code{\link{ebic}}.
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
out_mle <- mle(out)
out_mle

inherits(out_mle, "cglasso")
class(out_mle)

# inheriting method functions from 'cglasso': some examples
coef(out_mle, nrho = 10L, print.info = TRUE)
to_graph(out_mle, nrho = 10L, weighted = TRUE)
out_aic <- aic(out_mle)
out_aic
plot(out_mle, typeplot = "graph", gof = out_aic)

#################
# mglasso model #
#################
R <- event(X)
X <- as.matrix(X)
X[R == 1L] <- NA
out <- mglasso(X = X)
out_mle <- mle(out)
out_mle

inherits(out_mle, "mglasso")
class(out_mle)

# inheriting method functions from 'mglasso': some examples
coef(out_mle, nrho = 10L, print.info = TRUE)
to_graph(out_mle, nrho = 10L, weighted = TRUE)
out_aic <- aic(out_mle)
out_aic
plot(out_mle, typeplot = "graph", gof = out_aic)

################
# glasso model #
################
X <- MASS::mvrnorm(n = n, mu = mu, Sigma = diag(p))
out <- glasso(X)
out_mle <- mle(out)
out_mle

inherits(out_mle, "glasso")
class(out_mle)

# inheriting method functions from 'glasso': some examples
coef(out_mle, nrho = 10L, print.info = TRUE)
to_graph(out_mle, nrho = 10L, weighted = TRUE)
out_aic <- aic(out_mle)
out_aic
plot(out_mle, typeplot = "graph", gof = out_aic)
}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}

