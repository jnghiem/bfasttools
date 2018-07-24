#' Process MODIS raster layers for pixel reliability
#'
#' This function sets values in a MODIS raster dataset to \code{NA} according to
#' pixel reliability values. All raster outputs will be written to GeoTIFF
#' files.
#'
#' NB: If you are using \code{MODIStsp} to derive pixel reliability, please note
#' that it seems that pixel reliability 0 corresponds to \code{NA}.
#'
#' Since this function matches the layers in \code{data_layers} and
#' \code{pixel_reliability} by position, please ensure that all raster layers in
#' both datasets are sorted in the corrected order. For example, the first layer
#' in both datasets should correspond to the same time point.
#'
#' Note that writing to a single multilayered raster means that layer names will
#' not be preserved. That is to say, all layers will simply be named according
#' to the file name with a numerical index appended when the multilayered
#' GeoTIFF is read back into R. The drawback is that any date information
#' contained in the layer name will be lost. Thus, it is advisable to keep the
#' individual raster files so that one can reference back to the date
#' information.
#'
#' The \code{selected_PR} argument requires list containing a single named
#' numerical vector. The vector must be named either \code{"include"} or
#' \code{"exclude"} to determine the pixel reliability selection. If the vector
#' is named \code{"include"}, the numbers in the vector will be
#' \strong{included} in the pixel reliability values whose corresponding cells
#' will be filtered out (e.g. set to \code{NA}). If the vector is named
#' \code{"exclude"}, the numbers in the vector will be \strong{excluded} from
#' the pixel reliability values that will be used to filter out cells. Any pixel
#' reliability value that is not in the exclude vector will be used to filter
#' out cells.
#'
#' @param data_layers A RasterStack or RasterBrick object containing the MODIS
#'   data. A RasterBrick object is preferable to reduce processing time.
#' @param pixel_reliability A RasterStack or RasterBrick object containing the
#'   MODIS pixel reliability layers. A RasterBrick object is preferable to
#'   reduce processing time.
#' @param selected_PR A list of length 1 containing a named vector. This name
#'   must be either \code{"include"} or \code{"exclude"}. See Details for more
#'   information.
#' @param basename A character that will be prepended to all file names to be
#'   written.
#' @param output_directory A character file path to a directory where all files
#'   will be written.
#' @param bylayer A logical indicating if the individual raster outputs should
#'   remain separate (\code{TRUE}) or combined into a single multilayered raster
#'   file (\code{FALSE}).
#' @param clean A logical indicating if the individual raster outputs should be
#'   deleted once the multilayered raster has been created. Ignored if
#'   \code{bylayer=TRUE}.
#' @param mc.cores A numeric indicating the number of cores to be used in
#'   parallel computing.
#' @return Nothing. This function has the side-effect of writing to files.
#' @examples
#' \dontrun{
#' process_PR(mod.brick, pr, list(exclude=c(NA, 1)), "MODIS_EVI", "C:/Desktop/PR_processed", mc.cores=6)
#' }
#' @importFrom raster subset
#' @importFrom stringr str_pad
#' @export
process_PR <- function(data_layers, pixel_reliability, selected_PR, basename, output_directory, bylayer=FALSE, clean=FALSE, mc.cores=1) {
  num.width <- nchar(nlayers(data_layers))
  spr_names <- names(selected_PR)
  if (length(spr_names)!=1) {
    stop("Incorrect length for argument selected_PR.")
  } else if (!(spr_names %in% c("include", "exclude"))) {
    stop("The vector in selected_PR is incorrectly named.")
  }
  process_fun <- include_exclude_PR(spr_names)
  if (mc.cores>1) {
    cl <- snow::makeCluster(mc.cores, type="SOCK")
    doSNOW::registerDoSNOW(cl)
    foreach(i=1:nlayers(data_layers), .inorder=FALSE, .packages=c("raster", "rgdal", "stringr")) %dopar% {
      layer <- subset(data_layers, i)
      faulty_cells <- process_fun(subset(pixel_reliability, i), selected_PR[[1]])
      layer[faulty_cells] <- NA
      writeRaster(layer, filename=paste0(output_directory, "\\", basename, "_PR_processed_", str_pad(i, width=num.width, side="left", pad="0"), ".tif"), format="GTiff", overwrite=TRUE)
    }
    snow::stopCluster(cl)
  } else if (mc.cores==1) {
    for (i in 1:nlayers(data_layers)) {
      layer <- subset(data_layers, i)
      faulty_cells <- process_fun(subset(pixel_reliability, i), selected_PR[[1]])
      layer[faulty_cells] <- NA
      writeRaster(layer, filename=paste0(output_directory, "\\", basename, "_PR_processed_", str_pad(i, width=num.width, side="left", pad="0"), ".tif"), format="GTiff", overwrite=TRUE)
    }
  } else {
    stop("Incorrect mc.cores argument.")
  }
  if (!bylayer) {
    output_files <- paste0(output_directory, "\\", basename, "_PR_processed_", str_pad(1:nlayers(data_layers), width=num.width, side="left", pad="0"), ".tif")
    output_files <- output_files[file.exists(output_files)]
    print("Writing individual raster files to a multilayered raster file:")
    output_files %>%
      stack() %>%
      writeRaster(filename=paste0(output_directory, "\\", basename, "_PR_processed.tif"), format="GTiff", bylayer=FALSE, overwrite=TRUE, progress="text")
    if (clean) {
      unlink(output_files)
    }
  }
}
