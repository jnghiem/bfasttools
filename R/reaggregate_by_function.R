#' Re-aggregate raster layers according to function
#'
#' This function re-aggregates a multilayered raster dataset according to a
#' function (e.g. maximum, minimum, mean) and a time interval. All raster
#' outputs will be written to GeoTIFF files.
#'
#' Note that writing to a single multilayered raster means that layer names will
#' not be preserved. That is to say, all layers will simply be named according
#' to the file name with a numerical index appended when the multilayered
#' GeoTIFF is read back into R. The drawback is that any date information
#' contained in the layer name will be lost. Thus, it is advisable to keep the
#' individual raster files so that one can reference back to the date
#' information.
#'
#' @inheritParams process_PR
#' @param fun A function that returns a single value and has an \code{na.rm}
#'   argument, which will be set to \code{TRUE}. Examples are \code{mean, sum,
#'   max, median}, and \code{min}.
#' @param dates A vector of class Date or a number. If \code{dates} is a date
#'   vector, the data will be re-aggregated to monthly data by \code{fun}. If
#'   \code{dates} is a number, the data will be re-aggregated in groups
#'   specified by the number. For example, if \code{dates=4} the data will be
#'   reaggregated by grouping layers in groups of 4 and executing \code{fun}.
#' @return Nothing. This function has a side-effect of writing files to a
#'   directory.
#' @examples
#' \dontrun{
#' reaggregate_by_function(data_layers, max, "MODIS_max", "C:/Desktop", mc.cores=6)
#' }
#' @importFrom lubridate year month
#' @export
reaggregate_by_function <- function(data_layers, fun, dates, basename, output_directory, bylayer=FALSE, clean=FALSE, mc.cores=1) {
  class_dates <- class(dates)
  if (class_dates=="Date") {
    years <- year(dates); yr.range <- min(years):max(years)
    months <- month(dates); mn.range <- min(months):max(months)
  } else if (class_dates=="numeric" & length(dates)==1) {
    group_sequence <- seq(0, nlayers(data_layers), by=dates)[-1]
    no_groups_vec <- 1:length(group_sequence)
  } else {
    stop("Invalid format for argument dates.")
  }

  if ((mc.cores>1) & (class_dates=="Date")) {
    cl <- snow::makeCluster(mc.cores, type="SOCK")
    doSNOW::registerDoSNOW(cl)
    foreach::foreach(i=yr.range, .inorder=FALSE, .packages=c("raster", "stringr", "rgdal")) %dopar% {
      for (j in mn.range) {
        layers <- which(years==i & months==j)
        no_layers <- length(layers)
        if (no_layers==0) next
        sub <- subset(data_layers, layers)
        if (no_layers>1) {
          sub <- calc(sub, fun=fun, na.rm=TRUE)
        }
        writeRaster(sub, filename=paste0(output_directory, "\\", basename, "_reaggregated_", i, "_", str_pad(j, width=2, side="left", pad="0"), ".tif"), format="GTiff", overwrite=TRUE)
      }
    }
    snow::stopCluster(cl)
  } else if ((mc.cores==1) & (class_dates=="Date")) {
    for (i in yr.range) {
      for (j in mn.range) {
        layers <- which(years==i & months==j)
        no_layers <- length(layers)
        if (no_layers==0) next
        sub <- subset(data_layers, layers)
        if (no_layers>1) {
          sub <- calc(sub, fun=fun, na.rm=TRUE)
        }
        writeRaster(sub, filename=paste0(output_directory, "\\", basename, "_reaggregated_", i, "_", str_pad(j, width=2, side="left", pad="0"), ".tif"), format="GTiff", overwrite=TRUE)
      }
    }
  } else if ((mc.cores>1) & (class_dates=="numeric")) {
    cl <- snow::makeCluster(mc.cores, type="SOCK")
    doSNOW::registerDoSNOW(cl)
    foreach::foreach(i=no_groups_vec, .inorder=FALSE, .packages=c("raster", "stringr", "rgdal")) %dopar% {
      top_bnd <- group_sequence[i]
      sub <- subset(data_layers, (top_bnd-dates+1):top_bnd)
      no_layers <- nlayers(sub)
      if (no_layers>1) {
        sub <- calc(sub, fun=fun, na.rm=TRUE)
      }
      writeRaster(sub, filename=paste0(output_directory, "\\", basename, "_reaggregated_group_", str_pad(i, width=nchar(max(no_groups_vec)), side="left", pad="0")), format="GTiff", overwrite=TRUE)
    }
    snow::stopCluster(cl)
  } else {
    for (i in no_groups_vec) {
      top_bnd <- group_sequence[i]
      sub <- subset(data_layers, (top_bnd-dates+1):top_bnd)
      no_layers <- nlayers(sub)
      if (no_layers>1) {
        sub <- calc(sub, fun=fun, na.rm=TRUE)
      }
      writeRaster(sub, filename=paste0(output_directory, "\\", basename, "_reaggregated_group_", str_pad(i, width=nchar(max(no_groups_vec)), side="left", pad="0")), format="GTiff", overwrite=TRUE)
    }
  }
  if (!bylayer) {
    output_file_head <- paste0(output_directory, "\\", basename, "_reaggregated_")
    if (class_dates=="Date") {
      output_files <- paste0(output_file_head, as.vector(sapply(yr.range, paste0, "_", str_pad(mn.range, width=2, side="left", pad="0"), ".tif")))
    } else {
      output_files <- paste0(output_file_head, "group_", str_pad(no_groups_vec, width=nchar(max(no_groups_vec)), side="left", pad="0"), ".tif")
    }
    output_files <- output_files[file.exists(output_files)]
    print("Writing individual raster files to a multilayered raster file:")
    output_files %>%
      stack() %>%
      writeRaster(filename=paste0(output_directory, "\\", basename, "_reaggregated.tif"), format="GTiff", bylayer=FALSE, overwrite=TRUE, progress="text")
    if (clean) {
      unlink(output_files)
    }
  }
}
