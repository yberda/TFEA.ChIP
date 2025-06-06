% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/All-functions.R
\name{plot_CM}
\alias{plot_CM}
\title{Interactive HTML Plot for Transcription Factor Enrichment}
\usage{
plot_CM(CM.statMatrix, plot_title = NULL, specialTF = NULL, TF_colors = NULL)
}
\arguments{
\item{CM.statMatrix}{Output of the function 'getCMstats', a data frame 
containing Accession ID, Transcription Factor, Odds Ratio, p-value, and 
adjusted p-value.}

\item{plot_title}{The title for the plot (default: "Transcription Factor 
Enrichment").}

\item{specialTF}{(Optional) Named vector of TF symbols to be highlighted 
in the plot, allowing for grouped color representation.}

\item{TF_colors}{(Optional) Colors to highlight TFs specified in 
specialTF.}
}
\value{
A plotly scatter plot.
}
\description{
Generates an interactive HTML plot from a transcription 
factor enrichment table, output of the function 'getCMstats'.
}
