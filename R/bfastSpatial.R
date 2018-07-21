#' Spatial Breaks for Additive Season and Trend (BFAST)
#'
#' BFAST detects breakpoints in the linear trend component of a time series.
#' This function extends the methodology to spatial datasets. This function is
#' based on the code for \code{bfast()} in the \code{bfast} package.
#'
#' This function is quite computationally intensive, and may take on the order
#' of days to complete. One benchmark is that it takes approximately 20 hours to
#' process 180,000 cells in a raster dataset on a newer Windows 10 machine with
#' \code{mc.cores=6}. It is recommended that the user first executes
#' \code{sctest_p} on the dataset to perform a p-value adjustment and then use
#' this function's \code{level} argument to specify exactly which cells should
#' have breakpoints estimated.
#'
#' The returned value will be a data frame with the following format. Each cell
#' will have at least two rows in the data frame unless it contained too much
#' missing data (as defined by \code{nodata_threshold}). Given a cell number,
#' the first row contains information for the start date of the time series. Any
#' breakpoints detected will be listed in the rows subsequent. The final row for
#' a cell contains information on the ending date of the time series.
#'
#' The exact definitions for each field are: \enumerate{ \item \code{no_cell}:
#' The cell number. This number indexes all cells based on the raster grid of
#' \code{data_layers} starting from 1 in the upper left corner and increasing to
#' the right. \item \code{start_date}: The starting date for the time period up
#' to the next breakpoint or endpoint (except for the endpoint itself since
#' there is no data afterwards). \item \code{brk}: Either 0 or 1. 0 indicates
#' the row does not correspond to a breakpoint (e.g. it is a start or endpoint).
#' 1 indicates that the row is a breakpoint. \item \code{no_brk}: Total number
#' of breakpoints for a given cell. Only valid for a \code{start_date} for the
#' start point of the time series. \item \code{slope}: The slope of the trend
#' after the start or breakpoint in units of value/days. Not valid for the
#' endpoint. \item \code{length}: The length of time in days until the next
#' breakpoint or endpoint. \item \code{resid.upper}: The 75\% quantile of
#' residuals for the length of the time series from \code{start_date} to the
#' next breakpoint or endpoint. \item \code{resid.lower}: The 25\% quantile of
#' residuals for the length of the time series from \code{start_date} to the
#' next breakpoint or endpoint. \item \code{resid.sd}: The standard deviation of
#' the residuals for the length of the times series from \code{start_date} to
#' the next breakpoint or endpoint. \item \code{avg.trend}: The mean value in
#' the trend component of the length of the time series from \code{start_date}
#' to the next breakpoint or endpoint. \item \code{shift}: The positive or
#' negative change occurring at a breakpoint. Not valid for start or endpoints.
#' \item \code{resid.max, resid.min}: The maximum/minimum residual for the
#' length of the time series from \code{start_date} to the next breakpoint or
#' endpoint. \item \code{resid.max_date, resid.min_date}: The date of
#' \code{resid.max} or \code{resid.min}. }
#'
#' All dates in the returned data frame are in the format "YYYYMMDD." -9999 is a
#' placeholder value for all values that are invalid for a given field.
#'
#' @inheritParams sctest_p
#' @param dates A vector of class Date whose entries correspond to the date of
#'   each raster observation.
#' @param track_progress An optional character file path to a .txt file. If this
#'   argument is supplied, the function will write to the file the iteration
#'   number and time every 20,000 iterations through the cells.
#' @param breaks An integer indicating the maximum of breaks that should be
#'   estimated per time series. The default is the maximum number allowed by
#'   \code{h}.
#' @param level A number for the significance level of \code{sctest()}.
#'   Alternatively, this can be a vector of cells whose structural change
#'   p-values have been deemed significant (for running with \code{sctest_p}).
#' @return A data frame with with columns specified in Details.
#' @examples
#' \dontrun{
#' bfastSpatial(mod.brick, dates=monthly_date(2, 2000, 5, 2018, 15), monthly=TRUE, processingGeometry=pg, subset_values=1, track_progress="C:/Desktop/progress_log.txt", h=0.05, breaks=5, level=sig.cells, impute=TRUE, mc.cores=6)
#' }
#' @import gstat
#' @import raster
#' @importFrom foreach %dopar% foreach
#' @importFrom imputeTS na.interpolation
#' @importFrom strucchange sctest efp
#' @importFrom dplyr arrange mutate
#' @importFrom raster %in%
#' @export
bfastSpatial <- function(data_layers, dates, obs_per_year, processingGeometry=NULL, subset_values=NULL, nodata_threshold=c(0.05, 4),
                         track_progress=NULL, h=0.15, season="harmonic", max.iter=5, breaks=NULL, level=0.05, type="OLS-MOSUM", impute=FALSE, mc.cores=1) {

  class.ras <- class(data_layers)
  if (class.ras=="RasterStack") {
    data_layers <- brick(data_layers)
  } else if (class.ras!="RasterBrick") {
    stop("Argument data_layers must have class RasterStack or RasterBrick.")
  }

  if (class(dates)!="Date") {
    stop("dates must have class Date.")
  }

  if (length(dates)!=nlayers(data_layers)) {
    stop("The length of date vector must match the number of layers in the raster file.")
  }

  min.date <- min(dates) #start date
  max.date <- max(dates) #end date

  if (is.null(processingGeometry)) {
    cells <- 1:ncell(data_layers)
  } else {
    cells <- processingGeometry_cells(processingGeometry, subset_values, data_layers)
  }

  bfastC <- function(Yt, h, season =c("dummy", "harmonic", "none"), max.iter=NULL, breaks=NULL, level, type= "OLS-MOSUM") {
    season <- match.arg(season)
    ti <- time(Yt)
    f <- frequency(Yt)      # on cycle every f time points (seasonal cycle)
    if(class(Yt)!="ts")
      stop ("Not a time series object")
    ## return value
    Tt <- 0

    # seasonal model setup
    if (season=="harmonic") {
      w <- 1/f # f = 23 when freq=23 :-)
      tl <- 1:length(Yt)
      co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
      co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2)
      co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3)
      smod <- Wt ~ co+si+co2+si2+co3+si3
      # Start the iterative procedure and for first iteration St=decompose result
      St <- stl(Yt, "periodic")$time.series[, "seasonal"]
    } else if (season=="dummy") {
      # Start the iterative procedure and for first iteration St=decompose result
      St <- stl(Yt, "periodic")$time.series[, "seasonal"]
      D <- seasonaldummy(Yt)
      D[rowSums(D) == 0,] <- -1
      smod <- Wt ~ -1 + D
    } else if (season == "none") {
      print("No seasonal model will be fitted!")
      St <- 0
    } else {
      stop("Not a correct seasonal model is selected ('harmonic' or 'dummy') ")
    }
    # number/timing of structural breaks in the trend/seasonal component
    Vt.bp <- 0
    Wt.bp <- 0
    CheckTimeTt <- 1
    CheckTimeSt <- 1
    i <- 0
    while ( (!identical(CheckTimeTt,Vt.bp) | !identical(CheckTimeSt,Wt.bp)) & i < max.iter) {
      CheckTimeTt <- Vt.bp
      CheckTimeSt <- Wt.bp
      # TREND
      Vt <- Yt-St
      if (length(level)==1) {
        p.Vt <- sctest(efp(Vt ~ ti, h=h, type=type))
        bp.est <- p.Vt$p.value <= level
      } else {
        bp.est <- cell %in% level
      }

      if (bp.est) {
        bp.Vt <- breakpoints(Vt~ti, h=h, breaks=breaks)
        nobp.Vt <- is.na(breakpoints (bp.Vt)[1])
      } else {
        nobp.Vt <- TRUE
        bp.Vt <- NA
      }

      if (nobp.Vt) {
        fm0 <- lm(Vt ~  ti)
        Vt.bp <- 0      # no breaks times
        Tt <- ts(fitted(fm0))     # Data minus trend
        tsp(Tt) <- tsp(Yt)
        ci.Vt <- NA
      } else {
        fm1 <- lm(Vt ~ breakfactor(bp.Vt)/ti)
        ci.Vt <- confint(bp.Vt, het.err = FALSE)
        Vt.bp <- ci.Vt$confint[,2]
        Tt <- ts(fitted(fm1))     # Data minus trend
        tsp(Tt) <- tsp(Yt)
      }

      # SEASONAL COMPONENT
      Wt <- Yt-Tt
      sm0 <- lm(smod)
      St <- ts(fitted(sm0))  #  The fitted seasonal component
      tsp(St) <- tsp(Yt)
      Wt.bp <- 0             # no seasonal breaks

      i <- i+1
      output <- list(Tt=Tt,St=St,Nt=Yt-Tt-St, Vt=Vt, bp.Vt=bp.Vt, Wt=Wt)
    }
    return(output)
  }
  bfastCell <- function(cell) {
    ts <- ts(as.vector(data_layers[cell]), frequency=obs_per_year, start=0)
    ts.na <- is.na(ts)
    check.any <- any(ts.na)
    check.all <- all(ts.na)
    prop.na <- sum(ts.na)/length(ts)
    rle.na <- rle(ts.na)
    if (check.all | (!check.all & check.any & !impute) | impute & (prop.na>nodata_threshold[1] | max(rle.na$lengths[rle.na$values==1])>=nodata_threshold[2])) {
      add.df <- data.frame(no_cell=cell, start_date=NA, brk=NA, no_brk=NA, slope=NA, length=NA, resid.upper=NA, resid.lower=NA, resid.sd=NA,
                           avg.trend=NA, shift=NA, resid.max=NA, resid.max_date=NA, resid.min=NA, resid.min_date=NA)
    } else {
      if (check.any) {
        ts <- na.interpolation(ts, option="linear")
      }
      #Running bfast
      bf <- bfastC(ts, h=h, season=season, max.iter=max.iter, breaks=breaks, level=level, type=type)
      #Creating index for breakpoints and endpoints of time series
      bk_index <- bf$bp.Vt[[1]]
      bk_index_cleaned <- bk_index[!is.na(bk_index)]
      bk_index_full <- c(1, bk_index_cleaned, length(ts))
      #Creating output for diagnostic table
      brk <- c(0, rep(1, times=length(bk_index_cleaned)), 0) #indicator field to signal whether a breakpoint occurs at corresponding date
      date <- c(min.date, dates[bk_index_cleaned], max.date) #vector of dates for endpoints and breakpoints (if found)
      no_rows <- length(date)
      slope <- double(no_rows-1); resid.upper <- slope; resid.lower <- slope; resid.sd <- slope; avg.trend <- slope
      resid.max <- slope; resid.min <- slope; resid.max_date <- slope; resid.min_date <- slope
      shift <- c(-9999, rep(0, times=no_rows-2))
      trend <- bf$Tt
      overall_resids <- bf$Nt
      max.resid <- max(overall_resids); max.resid.date <- as.numeric(gsub("-", "", as.character(dates[which.max(overall_resids)])))
      min.resid <- min(overall_resids); min.resid.date <- as.numeric(gsub("-", "", as.character(dates[which.min(overall_resids)])))
      for (k in 2:no_rows) {
        rng <- bk_index_full[k-1]:bk_index_full[k]
        dates_subset <- dates[rng]
        #Calculating the slope
        values.lm <- trend[rng]
        days.lm <- seq(from=0, to=as.numeric(dates[max(rng)]-dates[min(rng)]), along.with=values.lm)
        slope[k-1] <- lm(values.lm~days.lm)[[1]][[2]]
        #Calculating average of trend segment
        avg.trend[k-1] <- mean(values.lm)
        #Calculating residual statistics
        resids <- bf$Nt[rng]
        resid.quant <- quantile(resids)
        resid.upper[k-1] <- resid.quant[[3]]; resid.lower[k-1] <- resid.quant[[2]]
        resid.sd[k-1] <- sd(resids)
        resid.max[k-1] <- max(resids); resid.min[k-1] <- min(resids)
        resid.max_date[k-1] <- as.numeric(gsub("-", "", as.character(dates_subset[which.max(resids)])))
        resid.min_date[k-1] <- as.numeric(gsub("-", "", as.character(dates_subset[which.min(resids)])))
        #Calculating shift in trend
        if (k!=no_rows) {
          shift[k] <- trend[bk_index_full[k]+1]-trend[bk_index_full[k]]
        }
      }
      no_cell <- rep(cell, times=no_rows) #cell number
      index <- 1:no_rows
      no_brk <- c(sum(brk), rep(-9999, times=no_rows-1))
      length <- as.numeric(diff(date))
      add.stats <- rbind(data.frame(slope, length, resid.upper, resid.lower, resid.sd, avg.trend, shift, resid.max, resid.max_date, resid.min, resid.min_date), -9999)
      start_date <- as.numeric(gsub("-", "", as.character(date))) #convert data type of date to be able to write to file
      add.df <- data.frame(no_cell, start_date, brk, no_brk, add.stats)
    }
    return(add.df)
  }
  if (mc.cores>1) {
    cl <- snow::makeCluster(mc.cores, type="SOCK")
    doSNOW::registerDoSNOW(cl)
    res <- foreach(i=1:length(cells), .combine=rbind, .inorder=FALSE, .packages=c("imputeTS", "strucchange")) %dopar% {
      if ((i %% 20000==0) & !is.null(track_progress)) {
        sink(track_progress, append=TRUE)
        cat(paste("Starting iteration", i, "on", date(), "\n"))
        sink()
      }
      cell <- cells[i]
      return(bfastCell(cell))
    }
    snow::stopCluster(cl)
  } else {
    res <- plyr::rbind.fill(lapply(cells, bfastCell))
  }
  return(arrange(res, no_cell, start_date))
}
