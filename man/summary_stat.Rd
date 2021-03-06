% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_stat.R
\name{summary_stat}
\alias{summary_stat}
\title{Summarize fields in \code{bfastSpatial()} output}
\usage{
summary_stat(bfast_in, stat, fun1, fun2, st_date = -Inf, ed_date = Inf)
}
\arguments{
\item{bfast_in}{A data frame generated by \code{bfastSpatial()}.}

\item{stat}{A character containing the name of a field in \code{bfast_in}.}

\item{fun1}{A function that returns a single value and accepts an
\code{na.rm} argument. This argument summarizes the values in \code{stat}
by cell number.}

\item{fun2}{A function that returns a single value and accepts an
\code{na.rm} argument. This argument summarizes the results of the
application of \code{fun1}.}

\item{st_date}{A Date indicating the earliest timing of observations to
include in the summary.}

\item{ed_date}{A Date indicating the latest timing of observations to inclide
in the summary.}
}
\value{
A number that summarizes \code{bfast_in} with \code{fun1} and
  \code{fun2}.
}
\description{
This function manipulates the data frame produced by \code{bfastSpatial} to
provide summary statistics for a specified field.
}
\details{
This function may be used for zonal statistics if combined with
\code{subset_bfast()}.
}
\examples{
\dontrun{
summary_stat(bf_df, "brk", sum, mean)
}
}
