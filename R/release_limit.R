#' @title release_limit
#'
#' @description Calculate the Release Limit for a product.
#'
#' @details Use the output from step1.down to run a release limit calculation.
#'
#' @param step1_down_object The fit object from the step1.down function (required).
#' @param shelf_temp The temperature at which the product will be stored.
#' @param shelf_time The time for which the product will be stored.
#' @param LSL A numeric value for the lower specification.
#' @param interval Choose between confidence intervals or prediction intervals (default is CI).
#' @param confidence_interval Confidence level for the confidence and prediction intervals
#'  around the predictions (default 0.95).
#' @param xname Label for the x-axis (optional).
#' @param yname Label for the y-axis (optional).
#' @param xlim x-axis limits (optional).
#' @param ylim y-axis limits (optional).
#'
#' @return Returns release limit calculations and plot.
#'
#' @examples
#' #data randomly generated for the purpose of the examples.
#' data(example_data)
#'
#' #run step1.down fit
#' fit1 <- step1_down(data = example_data, y = "conc", .time = "time",
#'  C = "Celsius",max_time_pred = 3)
#'
#' release_lim <- release_limit(step1_down_object = fit1, shelf_temp = 5,shelf_time = 3,
#'                              LSL = 45, interval = "CI", confidence_interval = 0.95)
#' release_lim$RL_plot
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export release_limit


release_limit <- function(step1_down_object, shelf_temp, shelf_time, LSL,
                          interval = "CI", confidence_interval = 0.95,
                          xname = NULL, yname = NULL, xlim = NULL, ylim = NULL){

  fit_object <- step1_down_object$fit

  dat = step1_down_object$data

  Kref = mean(dat$K)

  preds <- data.frame( # Making empty prediction frame
    temps = rep(shelf_temp,301),
    K = rep(shelf_temp + 273.15,301),
    conc = rep(NA,301),
    degrad = rep(NA,301),
    total_time = seq(0,shelf_time,length.out = 301))

  SIG = vcov(fit_object)
  sigma = summary(fit_object)$sigma

  if(step1_down_object$user_parameters$reparameterisation == F &&
     step1_down_object$user_parameters$zero_order == F){

  k1 = fit_object$par$k1
  k2 = fit_object$par$k2
  k3 = fit_object$par$k3
  c0 = fit_object$par$c0

  preds$degrad <- (1 - ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2 / (preds$K))))^(1/(1-k3)))
  preds$conc <- c0 - c0 * preds$degrad


    preds$derivk1 = c0 * preds$total_time * (-exp(k1 - k2/preds$K)) * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K)))^(1/(1 - k3) - 1)
    preds$derivk2 = (c0 * preds$total_time * exp(k1 - k2/preds$K) * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K)))^(1/(1 - k3) - 1)) / preds$K
    preds$derivk3 = c0 * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K)))^(1/(1 - k3)) * ((preds$total_time * exp(k1 - k2/preds$K)) / ((1 - k3)^2 * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K))) + log((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K)))/(1 - k3)^2)
    preds$derivc0 = 1 - preds$degrad

    preds$varY = (preds$derivk1)^2 * SIG[1,1] + (preds$derivk2)^2 * SIG[2,2] + (preds$derivk3)^2 * SIG[3,3] + (preds$derivc0)^2 * SIG[4,4] +
      2*preds$derivk1*preds$derivk2 * SIG[1,2] + 2*preds$derivk1*preds$derivk3 * SIG[1,3] + 2*preds$derivk1*preds$derivc0 * SIG[1,4] +
      2*preds$derivk2*preds$derivk3 * SIG[2,3] + 2*preds$derivk2*preds$derivc0 * SIG[2,4] + 2*preds$derivk3*preds$derivc0 * SIG[3,4]

    preds$derivk1 = preds$derivk2 = preds$derivk3 = NULL}else if(step1_down_object$user_parameters$reparameterisation == T &&
                                                                 step1_down_object$user_parameters$zero_order == T){
      k1 = fit_object$par$k1
      k2 = fit_object$par$k2
      c0 = fit_object$par$c0

      preds$degrad = preds$total_time * exp(k1 - k2/preds$K + k2/Kref)
      preds$conc = c0 - c0 * preds$degrad


      preds$derivk1 = -c0 * preds$degrad
      preds$derivk2 = -c0 * (1/Kref - 1/preds$K) * preds$degrad
      preds$derivc0 = 1 - preds$degrad
      preds$varY = (preds$derivk1)^2 * SIG[1, 1] + (preds$derivk2)^2 *
        SIG[2, 2] + (preds$derivc0)^2 * SIG[3, 3] + 2 * preds$derivk1 *
        preds$derivk2 * SIG[1, 2] + 2 * preds$derivk1 * preds$derivc0 *
        SIG[1, 3] + 2 * preds$derivk2 * preds$derivc0 * SIG[2,
                                                          3]
      preds$derivk1 = preds$derivk2 = preds$derivc0 = NULL}else if(step1_down_object$user_parameters$reparameterisation == F &&
                                                                   step1_down_object$user_parameters$zero_order == T){
        k1 = fit_object$par$k1
        k2 = fit_object$par$k2
        c0 = fit_object$par$c0

        preds$degrad = preds$total_time * exp(k1 - k2 / preds$K)
        preds$conc = c0 - c0 * preds$degrad

        preds$derivk1 = -c0 * preds$degrad
        preds$derivk2 = c0 / preds$K * preds$degrad
        preds$derivc0 = 1 - preds$degrad
        preds$varY = (preds$derivk1)^2 * SIG[1,1] + (preds$derivk2)^2 * SIG[2,2] + (preds$derivc0)^2 * SIG[3,3] +
          2*preds$derivk1*preds$derivk2 * SIG[1,2] + 2*preds$derivk1*preds$derivc0 * SIG[1,3] + 2*preds$derivk2*preds$derivc0 * SIG[2,3]

        preds$derivk1 = preds$derivk2 = preds$derivc0 = NULL}else if(step1_down_object$user_parameters$reparameterisation == T &&
                                                                     step1_down_object$user_parameters$zero_order == F){
          k1 = fit_object$par$k1
          k2 = fit_object$par$k2
          k3 = fit_object$par$k3
          c0 = fit_object$par$c0

          preds$degrad = 1 - ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2 / preds$K + k2 / Kref)))^(1/(1-k3))
          preds$conc = c0 - c0*preds$degrad

          preds$derivk1 = c0 * preds$total_time * (-exp(k1 - k2/preds$K + k2/Kref)) * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K + k2/Kref)))^(1/(1 - k3) - 1)
          preds$derivk2 = c0 * preds$total_time * (1/Kref - 1/preds$K) * (-exp(k1 - k2/preds$K + k2/Kref)) * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K + k2/Kref)))^(1/(1 - k3) - 1)
          preds$derivk3 = c0 * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K + k2/Kref)))^(1/(1 - k3)) * ((preds$total_time * exp(k1 - k2/preds$K + k2/Kref)) / ((1 - k3)^2 * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K + k2/Kref))) + log((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K + k2/Kref)))/(1 - k3)^2)
          preds$derivc0 = 1 - preds$degrad

          preds$varY = (preds$derivk1)^2 * SIG[1,1] + (preds$derivk2)^2 * SIG[2,2] + (preds$derivk3)^2 * SIG[3,3] + (preds$derivc0)^2 * SIG[4,4] +
              2*preds$derivk1*preds$derivk2 * SIG[1,2] + 2*preds$derivk1*preds$derivk3 * SIG[1,3] + 2*preds$derivk1*preds$derivc0 * SIG[1,4] +
              2*preds$derivk2*preds$derivk3 * SIG[2,3] + 2*preds$derivk2*preds$derivc0 * SIG[2,4] + 2*preds$derivk3*preds$derivc0 * SIG[3,4]
          preds$derivk1 = preds$derivk2 = preds$derivk3 = preds$derivc0 = NULL
        }

    if(interval == "CI"){
      preds$CI1b = preds$conc - qt((1+confidence_interval)/2, summary(fit_object)$df[2]) * sqrt(preds$varY)
      preds$CI2b = preds$conc + qt((1+confidence_interval)/2, summary(fit_object)$df[2]) * sqrt(preds$varY)
      # Add the degradation intervals
      preds$degrad_CI_up <- (c0 - preds$CI1b) / c0
      preds$degrad_CI_low <- (c0 - preds$CI2b) / c0

    }else{
      preds$PI1b = preds$conc - qt((1+confidence_interval)/2, summary(fit_object)$df[2]) * sqrt(preds$varY + sigma^2)
      preds$PI2b = preds$conc + qt((1+confidence_interval)/2, summary(fit_object)$df[2]) * sqrt(preds$varY + sigma^2)
      # Add the degradation intervals
      preds$degrad_PI_up <- (c0 - preds$PI1b) / c0
      preds$degrad_PI_low <- (c0 - preds$PI2b) / c0

    }

  # Calculate the RL
  if(interval == "CI"){degrad = as.numeric(preds %>% filter(total_time == shelf_time) %>% select(degrad_CI_up))}else{
    degrad = as.numeric(preds %>% filter(total_time == shelf_time) %>% select(degrad_PI_up))
  }

  RL = LSL / (1 - degrad)

  print(paste0("Release Limit: ",signif(RL,6))) # Print it

  print(paste0("Lower Specification Limit: ",LSL))

  # Now to add the new predictions with the release limit as c0 instead
  c0 = RL
  if(step1_down_object$user_parameters$reparameterisation == F &&
     step1_down_object$user_parameters$zero_order == F){
  preds$degradRL <- (1 - ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2 / (preds$K))))^(1/(1-k3)))
  }else if(step1_down_object$user_parameters$reparameterisation == T &&
           step1_down_object$user_parameters$zero_order == T){
    preds$degradRL <- preds$total_time * exp(k1 - k2/preds$K + k2/Kref)
  }else if(step1_down_object$user_parameters$reparameterisation == F &&
           step1_down_object$user_parameters$zero_order == T){
    preds$degradRL <- preds$total_time * exp(k1 - k2 / preds$K)
  }else if(step1_down_object$user_parameters$reparameterisation == T &&
           step1_down_object$user_parameters$zero_order == F){
    preds$degradRL <- 1 - ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2 / preds$K + k2 / Kref)))^(1/(1-k3))
    }

  preds$concRL <- c0 - c0 * preds$degradRL

  if (interval == "CI"){
    preds$concRL_up <- c0 - c0 * preds$degrad_CI_low
    preds$concRL_low <- c0 - c0 * preds$degrad_CI_up
  }else{
    preds$concRL_up <- c0 - c0 * preds$degrad_PI_low
    preds$concRL_low <- c0 - c0 * preds$degrad_PI_up
  }

  dat = step1_down_object$data
  pred = step1_down_object$prediction %>% filter(Celsius != shelf_temp)

  preds <- mutate(preds, Celsius = as.factor(temps))

  confidence_i <- paste0(confidence_interval * 100," % CI")
  prediction_i <- paste0(confidence_interval * 100," % PI")

  lines_t <- c("solid","dotted","longdash")
  names(lines_t) <- c("Prediction",confidence_i,prediction_i)

  if (is.null(xname))
    xname = "Time"
  if (is.null(yname))
    yname = "Response Variable"

  mytheme <- ggplot2::theme(legend.position = "bottom", strip.background = element_rect(fill = "white"),
                            legend.key = element_rect(fill = "white"), legend.key.width = unit(2,"cm"),
                            axis.text = element_text(size = 13), axis.title = element_text(size = 13),
                            strip.text = element_text(size = 13),
                            legend.text = element_text(size = 13),
                            legend.title = element_text(size = 13))

  if(interval == "CI"){
    plot1 = ggplot() + geom_point(data=dat, mapping=aes(x= time, y = y, colour = Celsius)) +
      ggtitle(paste0("Release Limit Estimation based on a ",confidence_interval*100,"% Confidence Interval and a Lower Specification Limit of ",LSL))+
      labs( x = xname, y = yname) +
      {if(!is.null(xlim))scale_x_continuous(limits = xlim)} +
      {if(!is.null(ylim))scale_y_continuous(limits = ylim)} +
      mytheme +
      geom_line(data=pred, mapping=aes(x= time, y = Response, colour = Celsius, linetype = "Prediction")) +
      geom_hline(yintercept = LSL, linetype = "dashed") +
      geom_line(data = preds, aes(x = total_time, y=conc, colour = Celsius, linetype = "Prediction")) +
      geom_line(data = preds, aes(x = total_time, y=concRL, colour = Celsius, linetype = "Prediction")) +
      geom_line(data = preds, aes(x = total_time, y=CI1b, colour = Celsius, linetype = confidence_i)) +
      geom_line(data = preds, aes(x = total_time, y=CI2b, colour = Celsius, linetype = confidence_i)) +
      geom_ribbon(data = preds, aes(x = total_time,ymin = concRL_low,ymax = concRL_up,fill = Celsius), alpha = 0.2) +
      scale_linetype_manual(name = NULL, values=lines_t) +
      scale_fill_discrete(name = "", labels = paste0("At Release ", round(RL, 1))) +
      theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))
  }else{
    plot1 = ggplot() + geom_point(data=dat, mapping=aes(x= time, y = y, colour = Celsius)) +
      ggtitle(paste0("Release Limit Estimation based on a ",confidence_interval*100,"% Prediction Interval and a Lower Specification Limit of ",LSL))+
      labs( x = xname, y = yname) +
      {if(!is.null(xlim))scale_x_continuous(limits = xlim)} +
      {if(!is.null(ylim))scale_y_continuous(limits = ylim)} +
      mytheme +
      geom_line(data=pred, mapping=aes(x= time, y = Response, colour = Celsius, linetype = "Prediction")) +
      geom_hline(yintercept = LSL, linetype = "dashed") +
      geom_line(data = preds, aes(x = total_time, y=conc, colour = Celsius, linetype = "Prediction")) +
      geom_line(data = preds, aes(x = total_time, y=concRL, colour = Celsius, linetype = "Prediction")) +
      geom_line(data = preds, aes(x = total_time, y=PI1b, colour = Celsius, linetype = prediction_i)) +
      geom_line(data = preds, aes(x = total_time, y=PI2b, colour = Celsius, linetype = prediction_i)) +
      geom_ribbon(data = preds, aes(x = total_time,ymin = concRL_low,ymax = concRL_up,fill = Celsius), alpha = 0.2) +
      scale_linetype_manual(name = NULL, values=lines_t) +
      scale_fill_discrete(name = "", labels = paste0("At Release ", round(RL, 1))) +
      theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))
  }

  if (is.null(step1_down_object$user_parameters$max_time_pred)){
    max_time_pred = max(dat$time)}else{
      max_time_pred = step1_down_object$user_parameters$max_time_pred
    }

  if(max_time_pred != shelf_time){
  print("To ensure predictions are the same length please align the values of max_time_pred in your step1.down call and shelf_time in your RL call")}

  results <- list(preds,plot1)
  names(results) <- c("Predictions","RL_plot")

  return(results)


}

globalVariables(c('total_time','degrad_CI_up','degrad_PI_up','temps','conc',
                  'concRL','CI1b','CI2b','concRL_low','concRL_up','PI1b','PI2b'))
