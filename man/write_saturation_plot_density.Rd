% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_saturation_plot_density.R
\name{write_saturation_plot_density}
\alias{write_saturation_plot_density}
\title{Write a density plot of the transitions and transversions mean values from each saturation test random sampling to disk}
\usage{
write_saturation_plot_density(data, dir = NULL, format = "png",
  show = TRUE, ...)
}
\arguments{
\item{data}{a \dQuote{oligodiag} object}

\item{dir}{path where the plot should be saved}

\item{format}{image format for the saved plot. Should be \dQuote{png}, \dQuote{jpg} or \dQuote{pdf}}

\item{show}{a logical whether to show the plot}
}
\description{
Write a density plot of the transitions and transversions mean values from each saturation test random sampling to disk
}
\examples{
write_saturation_plot_density(data)
}

