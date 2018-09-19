\name{rdatacggm}
\alias{rdatacggm}

\title{Simulate from a Censored Gaussian Graphical Model}

\description{
\sQuote{\code{rdatacggm}} function is used to produce one or more samples from the specified censored Gaussian graphical model.
}

\usage{
rdatacggm(n, mu, Sigma, probl, probr, lo, up, \dots)
}

\arguments{
\item{n}{the number of samples required.}
\item{mu}{a vector giving the means of the variables. By default all the expected vaules are equal to zero.}
\item{Sigma}{a positive-definite symmetric matrix specifying the covariance matrix of the variables. By default the identity matrix is used as covariance matrix.}
\item{probl}{a vector giving the probabilities that the random variables are left-censored.}
\item{probr}{a vector giving the probabilities that the random variables are right-censored.}
\item{lo}{a vector giving the left-censoring values.}
\item{up}{a vector giving the right-censoring values.}
\item{\dots}{futher arguments passed to the function \sQuote{\code{mvrnorm}}.}
}

\details{
\sQuote{\code{rdatacggm}} function simulates a dataset from a censored Gaussian graphical model and returns an object with class \sQuote{\code{datacggm}}.

The dataset is simulated in two steps:
\enumerate{
\item in the first step the arguments \code{n}, \code{mu}, \code{Sigma} and \code{\dots} are passed to the function \code{\link[MASS]{mvrnorm}} to simulate one or more samples from the specified multivariate Gaussian distribution.
\item in the second step, the values that are below or upper the censoring values are replaced.
}

The user can specify the censoring values in two equivalent ways. The simplest way is to use the arguments \code{lo} and \code{up}; a warning is produced if a full-observed dataset is simulated (see the last example). Alternatively, the censoring values can be implicitly specified using the arguments \code{probl} and \code{probr}. The \eqn{j}{j}th lower censoring value, denoted by \eqn{l_j}{l_j}, is such that: \deqn{\texttt{probl[j]} = \Pr\{X_j \le l_j\}.}{probl[j] = Pr\{X_j <= l_j\}.} In the same way, the \eqn{j}{j}th upper censoring value, denoted by \eqn{u_j}{u_j}, is such that: \deqn{\texttt{probr[j]} = \Pr\{X_j \ge u_j\}.}{probr[j] = Pr\{X_j >= u_j\}.}
}

\value{
\code{rdatacggm} returns an object with class \sQuote{\code{datacggm}}. See \code{\link{datacggm}} for more details.
}

\author{Luigi Augugliaro (\email{luigi.augugliaro@unipa.it})}

\references{
Augugliaro, L., Abbruzzo, A. and Vinciotti, V. (2018).
\eqn{\ell_1}{l1}-Penalized gaussian graphical model.
\emph{Biostatistics} (to appear).
}

\seealso{
\code{\link{datacggm}}, \code{\link{event}}, \code{\link{cglasso}} and the method function \code{\link{summary.datacggm}}.
}

\examples{
set.seed(123)

n <- 1000L
p <- 3L
mu <- rep(1L, p)
rho <- 0.3
Sigma <- outer(1L:p, 1L:p, function(i, j) rho^abs(i - j))

# simulating a dataset from a left-censored Gaussian graphical model
X <- rdatacggm(n = n, mu = mu, Sigma = diag(p), probl = 0.05)
# the same: X <- rdatacggm(n = n, mu = mu, Sigma = Sigma, probl = 0.05, probr = 0)
# the same: X <- rdatacggm(n = n, mu = mu, Sigma = Sigma, probl = 0.05, up = +Inf)
summary(X)

# simulating a dataset from a right-censored Gaussian graphical model
X <- rdatacggm(n = n, mu = mu, Sigma = diag(p), probr = 0.05)
# the same: X <- rdatacggm(n = n, mu = mu, Sigma = Sigma, probr = 0.05, probl = 0)
# the same: X <- rdatacggm(n = n, mu = mu, Sigma = Sigma, probr = 0.05, lo = -Inf)
summary(X)

# simulating a dataset from a censored Gaussian graphical model
X <- rdatacggm(n = n, mu = mu, Sigma = Sigma, probr = 0.05, probl = 0.05)
summary(X)

# simulating a full observed dataset: a warning is produced
X <- rdatacggm(n = n, mu = mu, Sigma = Sigma, lo = -10, up = 10)
summary(X)
}

\keyword{distribution}
\keyword{datagen}
\keyword{multivariate}

