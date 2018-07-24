#' Select MODIS files by date
#'
#' This function copies or combines to a multilayer raster MODIS files that
#' match the specified date criteria. Dates may be specified using Julian date
#' or traditional date.
#'
#' @param files A character vector of absolute file paths to the files from
#'   which to select.
#' @param output_directory A character file path to a directory where files will
#'   be written.
#' @param output_filename A character for the name of the output multilayered
#'   raster. Ignored if \code{bylayer=TRUE}.
#' @param selection A numeric vector of length 2 for the selection range
#'   (inclusive). The first value should be the date or Julian date of the start
#'   of the selection range. The second value should be the date or Julian date
#'   of the end of the selection range.
#' @param inverse A logical indicating if the selection range should be inverted
#'   (e.g. select files that are outside the selection range.) If \code{TRUE},
#'   the inverted selection range is exclusive.
#' @param convention A character. Set this argument to "standard" if the file
#'   names are in the standard MODIS naming convention. Set this argument to
#'   "MODIStsp" if the files have been processed using \code{MODIStsp} (and thus
#'   follow its naming convention).
#' @param bylayer A logical indicating if the individual raster outputs should
#'   remain separate (\code{TRUE}) or combined into a single multilayered raster
#'   file (\code{FALSE}).
#' @return Nothing. This function has the side-effect of writing to a file.
#' @examples
#' \dontrun{
#' select_MODIS_by_date(list.files("C:/Desktop/MODIS_files", pattern="\\.tif$", full.names=TRUE), "C:/Desktop/selected_MODIS", c(20, 120))
#' }
#' @importFrom stringr str_match
#' @importFrom dplyr mutate filter
#' @export
select_MODIS_by_date <- function(files, output_directory, output_filename=NULL, selection, inverse=FALSE, convention="standard", bylayer=FALSE) {
  if (convention=="standard") {
    regexpr <- "^MOD[[:alnum:]]+\\.A([[:digit:]]{4})([[:digit:]]{3})\\."
  } else if (convention=="MODIStsp") {
    regexpr <- "^MOD[[:alnum:]]+_[[:alnum:]]+_([[:digit:]]+)_([[:digit:]]+)"
  } else {
    stop('Argument convention must be "standard" or "MODIStsp."')
  }
  file_df <- str_match(basename(files), regexpr) %>%
    as.data.frame() %>%
    mutate(V2=as.numeric(as.character(V2)), V3=as.numeric(as.character(V3))) %>%
    cbind(data.frame(path=files))
  class_selection <- class(selection)
  len_selection <- length(selection)
  if (class_selection=="numeric" & len_selection==2) {
    if (inverse) {
      sel <- filter(file_df, V3<selection[1], V3>selection[2])
    } else {
      sel <- filter(file_df, V3>=selection[1], V3<=selection[2])
    }
  } else if (class_selection=="Date" & len_selection==2) {
    file_df <- mutate(file_df, date=base::as.Date(V3, origin=base::as.Date(paste(V2-1, "12", "31", sep="-"))))
    if (inverse) {
      sel <- filter(file_df, date<selection[1], date>selection[2])
    } else {
      sel <- filter(file_df, date>=selection[1], date<=selection[2])
    }
  } else {
    stop("Incorrect class or length for argument selection.")
  }
  sel <- as.character(sel[,"path"])
  if (bylayer) {
    file.copy(from=sel, to=output_directory, overwrite=TRUE)
    return(invisible(NULL))
  } else {
    if (is.null(output_filename)) {
      stop("output_filename must be supplied for writing to a multilayered raster.")
    }
    print("Writing individual raster files to a multilayered raster file:")
    sel %>%
      stack() %>%
      writeRaster(filename=paste0(output_directory, "\\", output_filename, ".tif"), format="GTiff", bylayer=FALSE, progress="text")
  }
}
