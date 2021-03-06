% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module.compare.R
\docType{methods}
\name{module.compare}
\alias{module.compare}
\title{Compare and plot the overlap among predicted modules}
\usage{
module.compare(res.module1, res.module2, used.mods1 = NULL,
  used.mods2 = NULL, type = c("model", "max.patients", "max.genes")[1],
  max.n1 = 30, max.n2 = max.n1, show.overlap = TRUE, cex = 10)
}
\arguments{
\item{res.module1}{a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}}

\item{res.module2}{a 'seed.module' or 'cluster.module' object returned by \code{\link{seed.module}} or \code{\link{cluster.module}}}

\item{used.mods1}{the modules to display}

\item{used.mods2}{the modules to display}

\item{type}{the module type to display}

\item{max.n1}{the maximum number of modules to display. If "used.mods1" is set, this option will be ignored.}

\item{max.n2}{the maximum number of modules to display. If "used.mods2" is set, this option will be ignored.}

\item{show.overlap}{boolean, display the overlap number}

\item{cex}{the font cex to display the overlap number}
}
\value{
The heatmap plot for gene overlaps.
}
\description{
Plot the overlap among predicted DEG modules
}
\details{
This function is to compare the modules from different studies, e.g. the different diseases or the different data for the same disease.
}
\examples{
\dontrun{
module.compare(res.mod1,res.mod2, type="model", max.n=20)
}
}
\author{
Guofeng Meng
}
\references{
Gu Z, Eils R and Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
}

