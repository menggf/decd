% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot.module.R
\name{Plot.cluster.module}
\alias{Plot.cluster.module}
\title{Plot the DEGs modules}
\usage{
\method{Plot}{cluster.module}(res.module, ann = NULL, deg = NULL,
  col.order = NULL, show.mods = NULL, overlap = NULL,
  dissimilarity = NULL, max.n = min(length(res.module), 30),
  type = c("model", "max.patients", "max.genes")[1], label.col = "#008000")
}
\arguments{
\item{res.module}{a 'cluster.module' object returned by  \code{\link{cluster.module}}}

\item{ann}{a data.frame for the patient annotation}

\item{deg}{a "deg" to display. It is returned by \code{\link{bi.deg}}}

\item{col.order}{the order of column in heatmap}

\item{show.mods}{a vector, the modules to display}

\item{overlap}{the similarity cutoff to display as carrying the module}

\item{dissimilarity}{the similarity cutoff to display as not carrying the module}

\item{max.n}{the maximum number of modules to display}

\item{type}{the module type to display}

\item{label.col}{the color to label}
}
\value{
A heatmap plot
}
\description{
Plot the DEGs modules
}
\details{
This function is to dispaly the relationship of the predicted DEG modules and the patients.

'deg' can be set to display the modules from different datasets, e.g. the modules predicted from disease A and display them in the binary DEG matrix of disease B.

The output is a heatmap Plot where the modules with maximum observations are showed.
}
\examples{
\dontrun{
Plot(cluster.mods1, er.ann, max.n=15)
Plot(cluster.mods2, er.ann, deg=deg, max.n=15)
# to display the ER-positive module in all the patients
}
}
\author{
Guofeng Meng
}
\references{
Gu Z, Eils R and Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
}

