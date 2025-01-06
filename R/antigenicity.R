#' Antigenicity Accelerated Stability Data
#'
#' An example dataset containing antigenicity concentration data at different
#' temperatures over a period of up to 407 days. Two points over 180 days are
#' to be used for validation instead of fitting.
#'
#' @docType data
#'
#' @usage data(antigenicity)
#'
#' @format An object of class \code{"data.frame"} with 56 rows and 6 variables
#' \describe{
#'  \item{time}{Number of days in years for which the datapoints are gathered.}
#'  \item{Celsius}{The temperature in celsius.}
#'  \item{K}{The temperature in Kelvin.}
#'  \item{conc}{The concentration at a time.}
#'  \item{N.days}{Number of days for which the datapoints are gathered.}
#'  \item{validA}{Whether the data point is to be used for validation or fitting.}
#'
#' }
#'
#' @keywords dataset
#'
"antigenicity"
