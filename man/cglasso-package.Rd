\name{cglasso-package}
\alias{cglasso-package}

\docType{package}

\title{L1-Penalized Censored Gaussian Graphical Model}

\description{The \eqn{\ell_1}{l1}-penalized censored Gaussian graphical model (Augugliaro \emph{and other}, 2018) is an extension of the graphical lasso estimator (Yuan \emph{and other}, 2007) developed to handle datasets from a censored Gaussian graphical model. An EM-like algorithm is implemented to fit the model. The graphical lasso algorithm (Friedman \emph{and other}, 2008) is used to solve the maximization problem in the M-step.}

\details{
\tabular{ll}{
Package: \tab cglasso\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2018-09-09\cr
License: \tab GPL (>=2)
}
}

\author{
Luigi Augugliaro\cr
Maintainer: Luigi Augugliaro  \email{luigi.augugliaro@unipa.it}
}

\references{
Augugliaro, L., Abbruzzo, A. and Vinciotti, V. (2018).
\eqn{\ell_1}{l1}-Penalized gaussian graphical model.
\emph{Biostatistics} (to appear).

Friedman, J.H., Hastie, T., and Tibshirani, R. (2008).
\href{https://academic.oup.com/biostatistics/article/9/3/432/224260}{Sparse inverse covariance estimation with the graphical lasso}.
\emph{Biostatistics} \bold{9}, 432--441.

Yuan, M. and Lin, Y. (2007).
\href{https://academic.oup.com/biomet/article-abstract/94/1/19/228821?redirectedFrom=fulltext}{Model selection and estimation in the Gaussian graphical model}.
\emph{Biometrika} \bold{94}, 19--35.

}

\keyword{multivariate}
\keyword{package}
\keyword{models}
