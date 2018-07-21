#' Create a monthly date vector
#'
#' This function produces a vector of dates that increment by month.
#'
#' The output of this function is most useful for passing as an argument to
#' \code{bfastSpatial()}, which requires a date vector.
#'
#' @param st_month A number indicating the month of the starting date.
#' @param st_year A number indicating the year of the starting date.
#' @param ed_month A number indicating the month of the ending date.
#' @param ed_year A number indicating the year of the ending date.
#' @param day_numer A number indicating the dummy date number that should be set
#'   for each date.
#' @return A vector of class Date spanning the start and end dates by month.
#' @examples
#' \dontrun{
#' monthly_date(2, 2000, 5, 2018)
#' }
#' @importFrom magrittr %>%
#' @export
monthly_date <- function(st_month, st_year, ed_month, ed_year, date_number=1) {
  dates <- sapply(st_year:ed_year, paste, 1:12, date_number, sep="-") %>%
    as.vector() %>%
    base::as.Date()
  start <- base::as.Date(paste(st_year, st_month, "1", sep="-"))
  if (ed_month==12) {
    ed_year <- ed_year+1
    ed_month <- 1
  } else {
    ed_month <- ed_month+1
  }
  end <- base::as.Date(paste(ed_year, ed_month, "1", sep="-"))
  return(dates[dates>=start & dates<end])
}
