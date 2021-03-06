% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_raster.R
\name{create_raster}
\alias{create_raster}
\title{Create RasterLayer object from \code{bfastSpatial()} output}
\usage{
create_raster(cell, variable, data, template, no_data = NA)
}
\arguments{
\item{cell}{A character with the name of the field in \code{data} containing
the cell numbers.}

\item{variable}{A character with the name of the field in \code{data}
containing the variable to be set as the pixel values in the raster.}

\item{data}{The data frame containing \code{cell} and \code{variable}.}

\item{template}{A Raster object with the same extent, resolution, and
projection as the raster dataset used to generate \code{data}.
Alternatively, the path to an XML file generated by
\code{create_raster_metadata}.}

\item{no_data}{A numeric value to fill in the rest of the returned
RasterLayer object where no value is specified. The default is \code{NA}.}
}
\value{
A RasterLayer object containing the specified pixel values.
}
\description{
This function formats tabular data into a RasterLayer.
}
\examples{
\dontrun{
create_raster(no_cell, total_breaks, bfast_data_frame, template)
}
}
