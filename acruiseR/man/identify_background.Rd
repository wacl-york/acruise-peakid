% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakid.R
\name{identify_background}
\alias{identify_background}
\title{Identifies background from a concentration time-series.}
\usage{
identify_background(
  concentration,
  method = c("gam", "rolling"),
  k = 10,
  bg_sd_window = 180,
  bg_sd_threshold = 0.5,
  bg_mean_window = 660
)
}
\arguments{
\item{concentration}{Concentration time-series as a vector.}

\item{method}{The background identification method to be employed.
\itemize{
\item{\code{gam}: Fits a Generalised Additive Model (gam), which is a non-linear
regression technique. Requires parameter \code{k} to be tuned.}
\item{\code{rolling}: The original 2-step rolling average method, which
firstly fits a rolling standard deviation and then a rolling mean to extract
the background. Requires parameters \code{bg_sd_window}, \code{bg_sd_threshold},
and \code{bg_mean_window} to be tuned.}
}}

\item{k}{The number of spline points in the GAM, a higher number
should provide a better fit but at the cost of overfitting to
noise. Used only
when \code{method='gam'}.}

\item{bg_sd_window}{Window size for the rolling standard deviation
smooth to identify the background, as an integer. Used only
when \code{method='rolling'}.}

\item{bg_sd_threshold}{Background measurements are considered as
those whose rolling sd is within this threshold
when \code{method='rolling'}. Used only}

\item{bg_mean_window}{The rolling mean to smooth the interpolated
background, as an integer. Used only
when \code{method='rolling'}.}
}
\value{
A list containing the background (slot \code{bg}) which is
a vector with the same length as \code{concentration}, and the arguments
used in this call (slot \code{call}).
}
\description{
This process takes 3 steps:
\enumerate{
\item Obtaining a smoothed rolling standard deviation of the background
\item Identifying background measurements as those that lie within a set
threshold
\item Interpolating values outside of this background threshold in a
linear manner and then using a rolling mean so that 'background'
measurements are available for the entire time-series.
}
}
