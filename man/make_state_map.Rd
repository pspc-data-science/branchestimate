% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch-estimating-functions.R
\name{make_state_map}
\alias{make_state_map}
\title{Renders a state choropleth for R effective.}
\usage{
make_state_map(choro_data)
}
\arguments{
\item{choro_data}{The tibble provided by make_state_choro_data().}
}
\value{
A state level choropleth for R effective in those states seeing an epidemic outbreak (R_eff>1).
}
\description{
This function creates a state choropleth depending on the output of make_state_choro_data().
}
