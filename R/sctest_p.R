#' Compute the structural change p-values for a spatiotemporal dataset
#'
#' This function executes \code{strucchange::sctest} for a spatiotemporal
#' dataset as a pre-processing step before BFAST analysis.
#'
#' The minimum p-value that \code{sctest} can return is 0.01. It is unclear if
#' this is a representative value for all p-values less than or equal to 0.01 or
#' if all p-values at 0.01 are indeed exactly equal to 0.01. This subtle detail
#' is important for performing p-value adjustments. For example, if many
#' p-values are actually much smaller than 0.01, the adjusted p-value may flag
#' them as not significant as a result of having them set at the nominal 0.01
#' level.
#'
#' @param data_layers A RasterBrick or RasterStack object. A RasterBrick is
#'   preferable to reduce processing time.
#' @param obs_per_year A number indicating the number of observations in the
#'   time series per year. For example, this argument should be set to 12 for
#'   monthly data. For MODIS 16-day aggregates, this should be set to 23
#'   (365/16, rounded to nearest integer.)
#' @param processingGeometry An optional RasterLayer, SpatialPolygons, or
#'   SpatialPolygonsDataFrame object that demarcates the boundary in which the
#'   function will operate. This argument will also take a RasterStack or
#'   RasterBrick object and compute the boundary based on the intersected
#'   region.
#' @param subset_values A vector or list of vectors that specifies which pixel
#'   values to select in \code{processingGeometry}. If
#'   \code{class(processingGeometry)} is \code{"RasterStack"} or
#'   \code{"RasterBrick"}, this argument must be a list of vectors and must be
#'   matched in order of layer. Ignored if \code{processingGeometry} is a vector
#'   dataset.
#' @param nodata_threshold A numeric vector of length 2. The first entry is the
#'   proportion of missing data points in each time series above which the
#'   function will not execute and return \code{NA} for all fields. The second
#'   entry is the number of consecutive missing data points above which the
#'   function will not execute and return \code{NA} for all fields.
#' @param h See documentation for \code{strucchange::sctest} for more details.
#' @param season See documentation for \code{strucchange::sctest} for more
#'   details.
#' @param type See documentation for \code{strucchange::sctest} for more
#'   details.
#' @param impute A logical indicating whether the function should linearly
#'   interpolate through time series with missing data if the time series does
#'   not exceed the thresholds defined by \code{nodata_threshold}.
#' @param mc.cores A numeric indicating the number of cores to be used in
#'   parallel computing.
#' @return A data frame with two columns: \code{no_cell} and \code{sc.p}.
#'   \code{no_cell} is the cell number (starting from 1 at the top-left of the
#'   raster grid and increasing by row). \code{sc.p} is the structural change
#'   test p-value.
#' @examples
#' \dontrun{
#' sctest_p(mod.brick, monthly=TRUE, processingGeometry=raster("C:/Desktop/boundary.tif"), subset_values=1, mc.cores=6)
#' }
#' @importFrom strucchange sctest efp
#' @importFrom imputeTS na.interpolation
#' @export
sctest_p <- function(data_layers, obs_per_year, processingGeometry=NULL, subset_values=NULL, nodata_threshold=c(0.05, 4), h=0.15, season="harmonic", type="OLS-MOSUM", impute=FALSE, mc.cores=1) {

  class.ras <- class(data_layers)
  if (class.ras=="RasterStack") {
    data_layers <- brick(data_layers)
  } else if (class.ras!="RasterBrick") {
    stop("Argument data_layers must have class RasterStack or RasterBrick.")
  }

  if (is.null(processingGeometry)) {
    cells <- 1:ncell(data_layers)
  } else {
    cells <- processingGeometry_cells(processingGeometry, subset_values, data_layers)
  }

  if (mc.cores>1) {
    cl <- snow::makeCluster(mc.cores, type="SOCK")
    doSNOW::registerDoSNOW(cl)
    res <- foreach::foreach(i=1:length(cells), .combine=rbind, .inorder=FALSE, .packages=c("imputeTS", "strucchange")) %dopar% {
      cell <- cells[i]
      ts <- ts(as.vector(data_layers[cell]), frequency=obs_per_year, start=0)
      ts.na <- is.na(ts)
      check.any <- any(ts.na)
      check.all <- all(ts.na)
      prop.na <- sum(ts.na)/length(ts)
      rle.na <- rle(ts.na)
      if (check.all | (!check.all & check.any & !impute) | impute & (prop.na>nodata_threshold[1] | max(rle.na$lengths[rle.na$values==1])>=nodata_threshold[2])) {
        return(data.frame(no_cell=cell, sc.p=NA))
      } else {
        if (check.any) {
          ts <- na.interpolation(ts, option="linear")
        }
        return(data.frame(no_cell=cell, sc.p=bf.sc(ts, h=h, season=season, type=type)))
      }
    }
    snow::stopCluster(cl)
  } else if (mc.cores==1) {
    res <- data.frame()
    for (i in 1:length(cells)) {
      cell <- cells[i]
      ts <- ts(as.vector(data_layers[cell]), frequency=obs_per_year, start=0)
      ts.na <- is.na(ts)
      check.any <- any(ts.na)
      check.all <- all(ts.na)
      prop.na <- sum(ts.na)/length(ts)
      rle.na <- rle(ts.na)
      if (check.all | (!check.all & check.any & !impute) | impute & (prop.na>nodata_threshold[1] | max(rle.na$lengths[rle.na$values==1])>=nodata_threshold[2])) {
        res <- rbind(res, data.frame(no_cell=cell, sc.p=NA))
      } else {
        if (check.any) {
          ts <- na.interpolation(ts, option="linear")
        }
        res <- rbind(res, data.frame(no_cell=cell, sc.p=bf.sc(ts, h=h, season=season, type=type)))
      }
    }
  } else {
    stop("Invalid mc.cores argument.")
  }
  return(res)
}
