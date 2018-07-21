#' Compute BFAST interactively
#'
#' This function allows the selection of a cell by clicking on a raster plot,
#' and then executes \code{bfast::bfast()} for the cell's time series. The
#' function then runs and plots \code{bfast()} on the time series.
#'
#' @inheritParams bfastSpatial
#' @param display A number or a RasterLayer object. If \code{display} is a
#'   number, the corresponding layer from \code{data_layers} will be plotted for
#'   the cell selection. If \code{display} is a RasterLayer object, it will be
#'   plotted for cell selection.
#' @param add_layers An optional list of Spatial objects to be plotted on top of
#'   \code{display}.
#' @param cell_number A number indicating the cell number through which
#'   \code{bfast()} should be run from \code{data_layers}. Ignored if
#'   \code{display} is supplied.
#' @param return_bfast A logical indicating if the BFAST results should be
#'   returned.
#' @param ... Arguments passed on to \code{bfast}.
#' @return A bfast object (if \code{return_bfast=TRUE}) or nothing otherwise.
#' @examples
#' \dontrun{
#' bfast_interactve(mod.brick, monthly=TRUE, display=raster("C:/Desktop/Mojave.tif", add_layers=list(boundary), impute=TRUE, nodata_threshold=c(0.05, 3))
#' }
#' @importFrom imputeTS na.interpolation
#' @importFrom bfast bfast
#' @importFrom raster subset
#' @export
bfast_interactive <- function(data_layers, obs_per_year, display=NULL, add_layers=NULL, cell_number=NULL, impute, nodata_threshold=c(0.05, 4), return_bfast=FALSE, ...) {

  if (!is.null(cell_number)) {
    ts <- ts(as.vector(data_layers[cell_number]), start=0, frequency=obs_per_year)
    if (impute) {
      if (check_ts(ts, nodata_threshold)) {
        ts <- na.interpolation(ts, "linear")
      } else {
        stop("The time series exceeds missing data thresholds.")
      }
    }
    bf <- bfast(ts, ...)
    plot(bf)
    if (return_bfast) {
      return(bf)
    }
  } else if (!is.null(display)) {
    class_display <- class(display)
    if (class_display=="numeric") {
      display <- subset(data_layers, display)
    } else if (class_display!="RasterLayer") {
      stop("Argument display must be a numeric or list.")
    }
    plot(display)
    if (!is.null(add_layers)) {
      for (i in 1:length(add_layers)) {
        plot(add_layers[[i]], add=TRUE)
      }
    }
    cell_info <- click(ras, n=1, xy=TRUE, cell=TRUE)
    ts <- ts(as.vector(data_layers[cell_info[1,"cell"]]), start=0, frequency=obs_per_year)
    if (impute) {
      if (check_ts(ts, nodata_threshold)) {
        ts <- na.interpolation(ts, "linear")
      } else {
        stop("The time series exceeds missing data thresholds.")
      }
    }
    bf <- bfast(ts, ...)
    plot(bf)
    if (return_bfast) {
      return(bf)
    }
  } else {
    stop("One of display and cell_number must be supplied.")
  }
}
