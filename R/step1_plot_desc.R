#' @title Plot Stability Data
#'
#' @description Plot raw accelerated stability data.
#'
#' @details Plot the raw accelerated stability data by selecting the columns -
#' response, time and temperature.
#'
#' @param data Dataframe containing accelerated stability data.
#' @param y Name of decreasing variable (e.g. concentration) contained within data
#' @param .time Time variable contained within data.
#' @param K Kelvin variable (numeric or column name) (optional).
#' @param C Celsius variable (numeric or column name) (optional).
#' @param xname Label for the x-axis (optional).
#' @param yname Label for the y-axis (optional).
#' @param xlim x-axis limits (optional).
#' @param ylim y-axis limits (optional).
#'
#' @return Plot of raw accelerated stability data.
#'
#' @examples
#' #load example datasets
#' data(antigenicity)
#' data(potency)
#'
#' plot1 <- step1_plot_desc(data=antigenicity, y="conc", .time="time", C = "Celsius")
#'
#' plot2 <- step1_plot_desc(data=potency, y="Potency", .time="Time", C = "Celsius")
#'
#' @import ggplot2
#'
#' @export step1_plot_desc

step1_plot_desc <- function (data, y, .time, K = NULL, C = NULL,
                             xname = NULL, yname = NULL, xlim = NULL, ylim = NULL){

  if (is.null(K) & is.null(C))
    stop("Select the temperature variable in Kelvin or Celsius")
  if (is.null(xname))
    xname = "Time"
  if (is.null(yname))
    yname = "Response Variable"

  dat = data

  if (is.null(C)){
    dat$C = dat[, K] - 273.15
    dat$Celsius = as.factor(dat$C)
  }else{
    dat$Celsius = as.factor(dat[, C])
  }

  dat$time = dat[, .time]
  dat$y = dat[, y]

  mytheme <- ggplot2::theme(legend.position = "bottom", strip.background = element_rect(fill = "white"),
                            legend.key = element_rect(fill = "white"), legend.key.width = unit(2,"cm"),
                            axis.text = element_text(size = 13), axis.title = element_text(size = 13),
                            strip.text = element_text(size = 13),
                            legend.text = element_text(size = 13),
                            legend.title = element_text(size = 13))

  plot <- ggplot2::ggplot(dat, aes(time, y, colour = Celsius)) + ggplot2::geom_point() +
   ggplot2::stat_summary(fun = mean, geom = "line") +
   labs( x = xname, y = yname) +
   {if(!is.null(xlim))scale_x_continuous(limits = xlim)} +
   {if(!is.null(ylim))scale_y_continuous(limits = ylim)} +
   mytheme

  return(plot)

}










