#' Create raster XML metadata file
#'
#' This function writes an XML file containing the information on the extent,
#' projection, and dimensions of the input dataset. This is designed to be an
#' alternative for the \code{template} argument for other functions.
#'
#' The information in the XML file does not create an exactly identical raster
#' in terms of extent and resolution compared to the original input raster, most
#' likely because of infinite precision issues. In other words, testing the
#' extent and resolution of a raster created from the XML file and the original
#' raster with \code{identical()} yields \code{FALSE}. However, testing them
#' with \code{all.equal()}, which allows for small numerical differences, yields
#' \code{TRUE}. It is likely that this is not a major problem because the
#' \code{raster} package still performs map algebra without issue. However, this
#' has not yet been tested in ArcGIS (and may apprently be remedied with the
#' Snap Raster environment setting if it is a problem). For all display
#' purposes, the rasters may be considered as exactly identical.
#'
#' @param data_layers A Raster object for which metadata will be created.
#' @param output_file A character file path to the XML file to which the
#'   metadata will be written.
#' @return Nothing. This function has a side-effect of writing to a file.
#' @examples
#' \dontrun{
#' create_raster_metadata(mod.brick, "C:/Desktop/metadata.xml")
#' }
#' @importFrom XML newXMLDoc newXMLNode saveXML xpathSApply xmlValue xmlParse
#' @importFrom sp proj4string
#' @export
create_raster_metadata <- function(data_layers, output_file) {
  ext <- extent(data_layers)
  doc <- newXMLDoc()
  root <- newXMLNode("raster", parent=doc)
  ext_node <- newXMLNode("extent", parent=root)
  newXMLNode("coord", attrs=list(type="xmin"), ext[1], parent=ext_node); newXMLNode("coord", attrs=list(type="xmax"), ext[2], parent=ext_node)
  newXMLNode("coord", attrs=list(type="ymin"), ext[3], parent=ext_node); newXMLNode("coord", attrs=list(type="ymax"), ext[4], parent=ext_node)
  newXMLNode("projection", proj4string(data_layers), parent=root)
  dim <- newXMLNode("dimensions", parent=root)
  newXMLNode("cells", attrs=list(dim="nrow"), nrow(data_layers), parent=dim)
  newXMLNode("cells", attrs=list(dim="ncol"), ncol(data_layers), parent=dim)
  saveXML(doc, file=output_file)
}
