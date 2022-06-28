% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakid.R
\name{integrate_aup_trapz}
\alias{integrate_aup_trapz}
\title{Integrate the Area Under a Plume (aup) using a trapezoidal method.}
\usage{
integrate_aup_trapz(conc, time, plumes, dx = 1)
}
\arguments{
\item{conc}{The concentration time-series as a vector with the background
removed.}

\item{plumes}{A Data Frame with 'start' and 'end' columns
containing plume boundaries, as returned by \code{detect_plumes}.}

\item{dx}{Sampling period, passed onto the dz argument of
np.trapz. I.e. the time between consecutive measurements.}
}
\value{
A Data Frame with one row per plume and 3 columns \code{start}, \code{end}, and
\code{area}. The first 2 are the same as in the input \code{plumes}, while \code{area}
contains the integrated area.
}
\description{
Integrate the Area Under a Plume (aup) using a trapezoidal method.
}