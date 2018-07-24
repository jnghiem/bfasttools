rsq <- function(observed, modeled) {
  y_bar <- mean(observed)
  tss <- sum((observed-y_bar)^2)
  rss <- sum((observed-modeled)^2)
  return(1-(rss/tss))
}
bf.sc <- function(ts, h, season, type) {
  ti <- time(ts)
  if (season=="harmonic") {
    St <- stl(ts, "periodic")$time.series[, "seasonal"]
  } else if (season=="dummy") {
    # Start the iterative procedure and for first iteration St=decompose result
    St <- stl(ts, "periodic")$time.series[, "seasonal"]
    D <- forecast::seasonaldummy(Yt)
    D[rowSums(D) == 0,] <- -1
    smod <- Wt ~ -1 + D
  } else if (season=="none") {
    print("No seasonal model will be fitted!")
    St <- 0
  } else {
    stop("Not a correct seasonal model is selected ('harmonic' or 'dummy') ")
  }
  Vt <- ts-St
  p.Vt <- sctest(efp(Vt ~ ti, h=h, type=type))
  return(p.Vt$p.value)
}
processingGeometry_cells <- function(processingGeometry, subset_values, data_layers) {
  class.pg <- class(processingGeometry)
  if (class.pg %in% c("RasterStack", "RasterBrick")) {
    no_layers <- nlayers(processingGeometry)
    for (i in 1:no_layers) {
      processingGeometry[[i]] <- processingGeometry[[i]] %in% subset_values[[i]]
    }
    check.ras <- calc(processingGeometry, fun=sum, na.rm=TRUE)
    cells <- which(as.vector(check.ras)==no_layers)
  }  else if (grepl("SpatialPolygon", class.pg)) {
    check.ras <- rasterize(processingGeometry, data_layers, field=1)
    cells <- which(as.vector(check.ras)==1)
  } else if (class.pg=="RasterLayer") {
    check.ras <- processingGeometry %in% subset_values
    cells <- which(as.vector(check.ras)==1)
  } else if (class.pg!="NULL"){
    stop("The processingGeometry input must be a RasterLayer or SpatialPolygons.")
  }
  return(cells)
}
create_annotation <- function(vector, position="topleft") {
  stat <- signif(c(mean(vector), median(vector), sd(vector)), 6)
  stat <- paste0(c("Mean ", "Median ", "Standard deviation "), "= ", stat)
  if (position=="topleft") {
    xpos <- -Inf; ypos <- Inf
  } else if (position=="topright") {
    xpos <- Inf; ypos <- Inf
  } else if (position=="bottomleft") {
    xpos <- -Inf; ypos <- -Inf
  } else if (position=="bottomright") {
    xpos <- Inf; ypos <- -Inf
  } else {
    stop("Position argument is not recognized.")
  }
  set1 <- c(1, 2.5, 4); set2 <- sort(-set1)
  hjustvar <- ifelse(sign(xpos)==1, 1, 0)
  if (sign(ypos)==1) {
    vjustvar <- set1
  } else {
    vjustvar <- set2
  }
  annotation <- data.frame(label=stat, xpos=xpos, ypos=ypos, hjustvar=hjustvar, vjustvar=vjustvar)
  return(geom_text(data=annotation, aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar, label=label)))
}
create_template_raster <- function(x) {
  class_x <- class(x)
  if (class_x=="RasterLayer") {
    return(x)
  } else if (any(class_x %in% paste0("Raster", c("Stack", "Brick")))) {
    return(raster::subset(x, 1))
  } else if (class_x=="character" & tools::file_ext(x)=="xml") {
    xml_doc <- xmlParse(x)
    ext <- extent(as.numeric(xpathSApply(xml_doc, "/raster/extent/coord", xmlValue)))
    proj <- xpathSApply(xml_doc, "//projection", xmlValue)
    dim <- as.numeric(xpathSApply(xml_doc, "//cells", xmlValue))
    return(raster(ext, nrows=dim[1], ncols=dim[2], crs=proj))
  } else {
    stop("Incorrect template argument.")
  }
}
effective_range <- function(model) {
  nugget <- model$psill[1]
  psill <- model$psill[2]
  target_gamma <- (nugget+psill)*0.95
  range_coef <- model$range[2]
  return(-range_coef*log(1-((target_gamma-nugget)/psill)))
}
clean_bfast <- function(bfast_in) {
  bfast_in %>%
    mutate(start_date=base::as.Date(as.character(start_date), format="%Y%m%d")) %>%
    mutate(year=year(start_date), month=month(start_date)) %>%
    return()
}
include_exclude_PR <- function(inex) {
  if (inex=="include") {
    return(function(x, vec) which(as.vector(x) %in% vec))
  } else if (inex=="exclude") {
    return(function(x, vec) which(!(as.vector(x) %in% vec)))
  } else {
    stop("Incorrect selected_PR argument.")
  }
}
check_ts <- function(ts, vec) {
  ts.na <- is.na(ts)
  prop.na <- sum(ts.na)/length(ts)
  rle.na <- rle(ts.na)
  test2 <- suppressWarnings(max(rle.na$lengths[rle.na$values==1]))
  return((prop.na<vec[1]) & (test2<vec[2]))
}
