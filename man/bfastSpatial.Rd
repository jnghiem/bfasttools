% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bfastSpatial.R
\name{bfastSpatial}
\alias{bfastSpatial}
\title{Spatial Breaks for Additive Season and Trend (BFAST)}
\usage{
bfastSpatial(data_layers, dates, obs_per_year, processingGeometry = NULL,
  subset_values = NULL, nodata_threshold = c(0.05, 4),
  track_progress = NULL, h = 0.15, season = "harmonic", max.iter = 5,
  breaks = NULL, level = 0.05, type = "OLS-MOSUM", impute = FALSE,
  mc.cores = 1)
}
\arguments{
\item{data_layers}{A RasterBrick or RasterStack object. A RasterBrick is
preferable to reduce processing time.}

\item{dates}{A vector of class Date whose entries correspond to the date of
each raster observation.}

\item{obs_per_year}{A number indicating the number of observations in the
time series per year. For example, this argument should be set to 12 for
monthly data. For MODIS 16-day aggregates, this should be set to 23
(365/16, rounded to nearest integer.)}

\item{processingGeometry}{An optional RasterLayer, SpatialPolygons, or
SpatialPolygonsDataFrame object that demarcates the boundary in which the
function will operate. This argument will also take a RasterStack or
RasterBrick object and compute the boundary based on the intersected
region.}

\item{subset_values}{A vector or list of vectors that specifies which pixel
values to select in \code{processingGeometry}. If
\code{class(processingGeometry)} is \code{"RasterStack"} or
\code{"RasterBrick"}, this argument must be a list of vectors and must be
matched in order of layer. Ignored if \code{processingGeometry} is a vector
dataset.}

\item{nodata_threshold}{A numeric vector of length 2. The first entry is the
proportion of missing data points in each time series above which the
function will not execute and return \code{NA} for all fields. The second
entry is the number of consecutive missing data points above which the
function will not execute and return \code{NA} for all fields.}

\item{track_progress}{An optional character file path to a .txt file. If this
argument is supplied, the function will write to the file the iteration
number and time every 20,000 iterations through the cells.}

\item{h}{See documentation for \code{strucchange::efp()} for more details.}

\item{season}{See documentation for \code{bfast::bfast()} for more details.}

\item{breaks}{An integer indicating the maximum of breaks that should be
estimated per time series. The default is the maximum number allowed by
\code{h}.}

\item{level}{A number for the significance level of \code{sctest()}.
Alternatively, this can be a vector of cells whose structural change
p-values have been deemed significant (for running with \code{sctest_p}).}

\item{type}{See documentation for \code{strucchange::efp()} for more details.}

\item{impute}{A logical indicating whether the function should linearly
interpolate through time series with missing data if the time series does
not exceed the thresholds defined by \code{nodata_threshold}.}

\item{mc.cores}{A numeric indicating the number of cores to be used in
parallel computing.}
}
\value{
A data frame with with columns specified in Details.
}
\description{
BFAST detects breakpoints in the linear trend component of a time series.
This function extends the methodology to spatial datasets. This function is
based on the code for \code{bfast()} in the \code{bfast} package.
Specifically, the code from line 123 to 207 are derived from the code for
\code{bfast()}, and this package's author takes no credit for it.
}
\details{
This function is quite computationally intensive, and may take on the order
of days to complete. One benchmark is that it takes approximately 20 hours to
process 180,000 cells in a raster dataset on a newer Windows 10 machine with
\code{mc.cores=6}. It is recommended that the user first executes
\code{sctest_p} on the dataset to perform a p-value adjustment and then use
this function's \code{level} argument to specify exactly which cells should
have breakpoints estimated.

The returned value will be a data frame with the following format. Each cell
will have at least two rows in the data frame unless it contained too much
missing data (as defined by \code{nodata_threshold}). Given a cell number,
the first row contains information for the start date of the time series. Any
breakpoints detected will be listed in the rows subsequent. The final row for
a cell contains information on the ending date of the time series.

The exact definitions for each field are: \enumerate{ \item \code{no_cell}:
The cell number. This number indexes all cells based on the raster grid of
\code{data_layers} starting from 1 in the upper left corner and increasing to
the right. \item \code{start_date}: The starting date for the time period up
to the next breakpoint or endpoint (except for the endpoint itself since
there is no data afterwards). \item \code{brk}: Either 0 or 1. 0 indicates
the row does not correspond to a breakpoint (e.g. it is a start or endpoint).
1 indicates that the row is a breakpoint. \item \code{no_brk}: Total number
of breakpoints for a given cell. Only valid for a \code{start_date} for the
start point of the time series. \item \code{slope}: The slope of the trend
after the start or breakpoint in units of value/days. Not valid for the
endpoint. \item \code{length}: The length of time in days until the next
breakpoint or endpoint. \item \code{resid.lower}: The 25\% quantile of
residuals for the length of the time series from \code{start_date} to the
next breakpoint or endpoint. \item \code{resid.median}: The median of
residuals for the length of the time series from \code{start_date} to the
next breakpoint or endpoint. \item \code{resid.upper}: The 75\% quantile of
residuals for the length of the time series from \code{start_date} to the
next breakpoint or endpoint. \item \code{resid.max, resid.min}: The
maximum/minimum residual for the length of the time series from
\code{start_date} to the next breakpoint or endpoint. \item
\code{resid.max_date, resid.min_date}: The date of \code{resid.max} or
\code{resid.min}. \item \code{resid.rsq}: The coefficient of determination
(\emph{R^2}) for the length of the time series from \code{start_date} to the
next breakpoint or endpoint. \item \code{resid.sd}: The standard deviation of
the residuals for the length of the time series from \code{start_date} to the
next breakpoint or endpoint. \item \code{refit_slope_p}: The p-value of the
t-statistics for the slope coefficient for an OLS model of time vs.
residuals. If this value is large, we do not reject the null hypothesis that
this "re-fitted" slope is 0. Note that this is the p-value for a single test,
so multiple testing corrections are necessary to draw conclusions for more
than one test. \item \code{avg.trend}: The mean value in the trend component
of the length of the time series from \code{start_date} to the next
breakpoint or endpoint. \item \code{shift}: The positive or negative change
occurring at a breakpoint. Not valid for start or endpoints.}

All dates in the returned data frame are in the format "YYYYMMDD." -9999 is a
placeholder value for all values that are invalid for a given field.

Since this function tends to have long processing times, it is recommended
that the output be saved to a tabular file format to avoid repeated function
calls. However, writing this file and reading it back into R may be
time-consuming because file size tends to be large as well. As a result,
\code{write.csv()} and \code{read.csv()} should not be used. \code{fwrite()}
and \code{fread()} in the \code{data.table} package are much faster
alternatives and should be used instead.
}
\examples{
\dontrun{
bfastSpatial(mod.brick, dates=monthly_date(2, 2000, 5, 2018, 15), monthly=TRUE, processingGeometry=pg, subset_values=1, track_progress="C:/Desktop/progress_log.txt", h=0.05, breaks=5, level=sig.cells, impute=TRUE, mc.cores=6)
}
}
