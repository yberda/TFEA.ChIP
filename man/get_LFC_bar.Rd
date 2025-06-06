% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/All-functions.R
\name{get_LFC_bar}
\alias{get_LFC_bar}
\title{Plots a color bar from log2(Fold Change) values.}
\usage{
get_LFC_bar(LFC)
}
\arguments{
\item{LFC}{Vector of log2(fold change) values arranged from higher
to lower. Use only the values of genes that have an Entrez ID.}
}
\value{
Plotly heatmap plot -log2(fold change) bar-.
}
\description{
Function to plot a color bar from log2(Fold Change)
values from an expression experiment.
}
\examples{
get_LFC_bar(arranged.log2FC.array)
}
