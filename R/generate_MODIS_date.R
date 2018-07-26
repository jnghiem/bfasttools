#' Generate a date vector for MODIS aggregate products
#'
#' This function generates a date vector that correponds to the observation
#' dates for MODIS aggregated products (e.g. 8-day and 16-day aggregagtes).
#'
#' This function assumes that the dates for a given year begin on January 1 and
#' increment by \code{interval}. This pattern matches that of 8-day and 16-day
#' aggregate products, but tests against other products have not been performed.
#'
#' @param start_date A Date or character of format YYYY-MM-DD for the starting
#'   date of the vector.
#' @param end_date A Date or character of format YYYY-MM-DD for the ending date
#'   of the vector.
#' @param interval The time in days between each entry in the vector. For
#'   example, \code{interval} should be 16 for 16-day data.
#' @param subset A numeric vector of length 2 or 4. If the length is 2,
#'   \code{subset} defines a yearly subset of the date vector by Julian date,
#'   inclusive. For example, \code{subset=c(50, 150)} would define a subset to
#'   keep only the dates every year whose Julian dates fall within that range.
#'   If the length is 4, \code{subset} is defined in terms of a starting month
#'   and date and an ending month and date (e.g. \code{c(2, 4, 6, 27)} refers to
#'   a subset starting February 4 and ending June 27 for each year).
#' @param inverse A logical indicating if the subset defined by \code{subset}
#'   should be inverted to find only dates \strong{outside} the subset.
#' @return A vector of class Date containing dates corresponding to MODIS data
#'   dates.
#' @examples
#' \dontrun{
#' generate_MODIS_date("2000-02-28", "2009-07-02", 8)
#' }
#' @importFrom lubridate year yday
#' @export
generate_MODIS_date <- function(start_date, end_date, interval, subset=NULL, inverse=FALSE) {
  start_date <- base::as.Date(start_date); end_date <- base::as.Date(end_date)
  min_length <- ceiling(365/interval)
  date_vec <- integer(0)
  class(date_vec) <- "Date"
  for (i in year(start_date):year(end_date)) {
    add_vec <- seq.Date(from=base::as.Date(paste0(i, "-01-01")), by=interval, length.out=min_length)
    add_vec <- add_vec[year(add_vec)==i]
    date_vec <- c(date_vec, add_vec)
  }
  date_vec <- date_vec[date_vec>=start_date & date_vec<=end_date]
  if (!is.null(subset)) {
    sub_length <- length(subset)
    if (sub_length==2) {
      jul <- yday(date_vec)
      if (inverse) {
        date_vec <- date_vec[jul<min(subset) | jul>max(subset)]
      } else {
        date_vec <- date_vec[jul>=min(subset) & jul<=max(subset)]
      }
    } else if (sub_length==4) {
      yr <- range(year(date_vec))
      sub_dates <- integer(0)
      class(sub_dates) <- "Date"
      for (i in yr[1]:yr[2]) {
        add_sub <- seq.Date(from=base::as.Date(paste(i, subset[1], subset[2], sep="-")), to=base::as.Date(paste(i, subset[3], subset[4], sep="-")), by=1)
        sub_dates <- c(sub_dates, add_sub)
      }
      if (inverse) {
        sub_dates <- base::setdiff(seq.Date(from=base::as.Date(paste0(yr[1], "-01-01")), to=base::as.Date(paste0(yr[2], "-12-31")), by=1), sub_dates)
      }
      date_vec <- date_vec[date_vec %in% sub_dates]
    } else {
      stop("The length for subset must be 2 or 4.")
    }
  }
  return(date_vec)
}
