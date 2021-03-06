#' Normalized difference vegetation index for a region near Salton Sea,
#' California
#'
#' A spatial dataset containing the normalized difference vegetation index
#' (NDVI) derived from Terra MODIS collection 6 16-day composites (MOD13Q1). The
#' dataset begins on February 18, 2000 and ends on May 25, 2018, with one
#' observation every 16 days. The NDVI values are scaled by 10000.
#'
#' This data was downloaded using the \code{MODIStsp} package. The layer names
#' reflect the package's naming convention of \code{"MODIS_(product)_(year)_(Julian
#' date)."}
#'
#' @format A RasterBrick object with 421 layers. Each layer has 5040 cells (70
#'   rows, 72 columns).
#'
#' @source
#' \url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1}
#'
"CA_ndvi"
