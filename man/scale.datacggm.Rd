\name{scale.datacggm}
\alias{scale.datacggm}

\title{Scaling and Centering of \dQuote{datacggm} Objects}

\description{
The method function \code{scale.datacggm} centers and/or scales the columns of a numeric matrix storaged in a \sQuote{\code{datacggm}} object.
}

\usage{
\method{scale}{datacggm}(x, center = TRUE, scale = TRUE)
}

\arguments{
\item{x}{an object of class \sQuote{\code{datacggm}}.}
\item{center}{either a logical value or numeric-alike vector of length equal to the number of columns of \code{x}, where \sQuote{numeric-alike} means that \code{\link[base]{as.numeric}(.)} will be applied successfully if \code{\link[base]{is.numeric}(.)} is not true.}
\item{scale}{either a logical value or a numeric-alike vector of length equal to the number of columns of \code{x}.}
}

\value{
The method function \sQuote{\code{scale.datacggm}} returns an object of class \code{\link{datacggm}}. The numeric centering and scalings used (if any) are returned as attributes \code{"scaled:center"} and \code{"scaled:scale"}.
}

\details{
  The value of \code{center} determines how column centering is performed. If \code{center} is a numeric-alike vector with length equal to the number of columns of \code{x}, then each column of \code{x} has the corresponding value from \code{center} subtracted from it.  If \code{center} is \code{TRUE} then centering is done by subtracting the column means (omitting censoring values) of \code{x$X} from their corresponding columns, and if \code{center} is \code{FALSE}, no centering is done. The same is done for \code{x$lo} and \code{x$up}.

  The value of \code{scale} determines how column scaling is performed (after centering).  If \code{scale} is a numeric-alike vector with length equal to the number of columns of \code{x}, then each column of \code{x$X} is divided by the corresponding value from \code{scale}. If \code{scale} is \code{TRUE} then scaling is done by dividing the (centered) columns of \code{x$X} by their standard deviations if \code{center} is \code{TRUE}, and the root mean square otherwise. If \code{scale} is \code{FALSE}, no scaling is done. The same is done for \code{x$lo} and \code{x$up}.

The root-mean-square for a (possibly centered) column is defined as \eqn{\sqrt{\sum(x^2)/(n-1)}}{sqrt(sum(x^2)/(n-1))}, where \eqn{x} is a vector of observed values and \eqn{n} is the number of observed values.  In the case \code{center = TRUE}, this is the same as the standard deviation, but in general it is not.  (To scale by the standard deviations without centering, use \code{scale(x, center = FALSE, scale = apply(x, 2, sd, na.rm = TRUE))}.)
}

\seealso{
  \code{\link{datacggm}}.
}

\examples{
set.seed(123)

n <- 100L
p <- 3L
mu <- rep(1L, p)
X <- rdatacggm(n = n, mu = mu, probr = 0.05, probl = 0.05)
centered.X <- scale(X)
centered.X
attr(centered.X, "scaled:center")
attr(centered.X, "scaled:scale")
}

\keyword{array}
