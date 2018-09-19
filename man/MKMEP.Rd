\name{MKMEP}
\docType{data}
\alias{MKMEP}

\title{
Megakaryocyte-Erythroid Progenitors
}

\description{
In a study about the formation of blood cells, Psaila \emph{and others} (2016) have recently identified three distinct subpopulations of cells, which are all derived from hematopoietic stem cells through cell differentiation. One of these sub-populations, denoted by MK-MEP, is a previously unknown, rare population of cells that are bipotent but primarily generate megakaryocytic progeny.

\sQuote{\code{MKMEP}} in an object of class \sQuote{\link{datacggm}} containing a subset of the data available from Psaila \emph{and others} (2016).

Multiplex RT-qPCR is used to profile 63 genes and 48 single human MK-MEP cells. RT-qPCR data are typically right-censored with a limit of detection fixed by the manufacturer to 40. Raw data have been mean normalized using the method proposed in Pipelers \emph{and others} (2017). See Section 5 in Augugliaro \emph{and others} (2018) for more details.
}

\usage{data("MKMEP")}

\references{
Augugliaro, L., Abbruzzo, A. and Vinciotti, V. (2018).
\eqn{\ell_1}{l1}-Penalized gaussian graphical model.
\emph{Biostatistics} (to appear).

Pipelers, P., Clement, L., Vynck, M., Hellemans, J., Vandesompele, J. and Thas, O. (2017).
\href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182832}{A unified censored normal regression model for qPCR differential gene expression analysis}.
\emph{PLoS One} \bold{12}, e0182832.

Psaila, B., Barkas, N., Iskander, D., Roy, A., Anderson, S., Ashley, N., Caputo, V. S., Lichtenberg, J., Loaiza, S., Bodine, D. M. \emph{and others}. (2016).
\href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0939-7}{Single-cell profiling of human megakaryocyte-erythroid progenitors identifies distinct megakaryocyte and erythroid differentiation pathways}.
\emph{Genome Biology} \bold{17}, 83--102.
}

\seealso{
\code{\link{cglasso}}, \code{\link{to_graph}}, and the method functions \code{\link{summary}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{aic}}, \code{\link{bic}} and \code{\link{ebic}}.
}

\examples{
data("MKMEP")
out_cglasso <- cglasso(MKMEP, pendiag = TRUE, nrho = 200L, rho.min.ratio = 0.35)
out_ebic <- ebic(out_cglasso)
plot(out_ebic, type = "l")

out_graph <- to_graph(out_cglasso, nrho = which.min(out_ebic$value_gof))
V(out_graph)$color <- "white"
V(out_graph)$frame.color <- NA
V(out_graph)$label.color <- "black"
E(out_graph)$color <- "gray75"

plot(out_graph, layout = layout_with_lgl(out_graph))
mtext(text = "Megakaryocytic MEP population", cex = 1.5, line = 1)
}

\keyword{datasets}
