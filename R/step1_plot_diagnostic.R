#' @title  Create Diagnostic Plots
#'
#' @description Generate residual diagnostic plots from a step1_down fit.
#'
#' @details Use the fit object obtained from the step1_down function to plot the
#' residual diagnostic plots, assess the quality of fit and search for anomalies.
#' Change the type of Residuals assessed.
#' Plots created are: Residuals Histogram, QQ Plot of Residuals, Observed Vs Predicted results, Residuals
#'  Vs Predicted results and Residuals By Time.
#'
#' @param step1_down_object The fit object from the step1_down function (required).
#' @param bins The number of bins in the Histogram plot (default 7).
#' @param residuals The type of residuals to plot classic/studentized/standardized (default classic).
#'
#' @return A list containing the five ggplot2 plots.
#'
#' @examples
#' #load antigenicity data
#' data(antigenicity)
#'
#' #run step1_down fit
#' fit1 <- step1_down(data = antigenicity, y = "conc", .time = "time",
#'  C = "Celsius", max_time_pred = 3)
#'
#' #plot diagnostic plots to asses the fit
#' step1_plot_diagnostic(fit1)
#'
#' @import ggplot2
#' @importFrom stats dnorm sd qqnorm
#'
#' @export step1_plot_diagnostic

step1_plot_diagnostic <- function(step1_down_object, bins = 7, residuals = "classic")
{
  if (is.null(step1_down_object))
    stop("First, run the model")

  dat = step1_down_object$data

  mytheme <- ggplot2::theme(legend.position = "bottom", strip.background = element_rect(fill = "white"),
                            legend.key = element_rect(fill = "white"), legend.key.width = unit(2,"cm"),
                            axis.text = element_text(size = 13), axis.title = element_text(size = 13),
                            strip.text = element_text(size = 13),
                            legend.text = element_text(size = 13),
                            legend.title = element_text(size = 13))

  validation = step1_down_object$user_parameters$validation

  if(!is.null(validation)){
 val_dat <- dat[dat$validation == "Validation",] ## Selects only validation rows
 dat <- dat[dat$validation == "Fit",] ## Selects only non-validation rows

  }

  dat$residuals <- summary(step1_down_object$fit)$residuals
  #dat$y = dat[, step1_down_object$user_parameters$y]
  dat$predicted <- dat$y - dat$residuals

# Generate predictions and residuals for validation data if present and append to data used in model fitting
if(!is.null(validation)){
     Kref = mean(dat$K)
     c0 = step1_down_object$fit$par$c0
     k1 = step1_down_object$fit$par$k1
     k2 = step1_down_object$fit$par$k2
     if(step1_down_object$user_parameters$zero_order == FALSE) {
     k3 = step1_down_object$fit$par$k3 }

     if(step1_down_object$user_parameters$reparameterisation == TRUE & step1_down_object$user_parameters$zero_order == TRUE) {     ## Model type 1: reparameterisation and k3 = 0
     val_degrad = val_dat$time * exp(k1 - k2/val_dat$K + k2/Kref)
     val_predicted = c0 - c0 * val_degrad
     } else if(step1_down_object$user_parameters$reparameterisation == FALSE & step1_down_object$user_parameters$zero_order == TRUE){  ## Model type 2: no reparameterisation and k3 = 0
     val_degrad = val_dat$time * exp(k1 - k2 / val_dat$K)
     val_predicted = c0 - c0 * val_degrad
     } else if(step1_down_object$user_parameters$reparameterisation == TRUE & step1_down_object$user_parameters$zero_order == FALSE) {  ## Model type 3: reparameterisation and k3 is not zero
     val_degrad = 1 - ((1 - k3) * (1/(1 - k3) - val_dat$time * exp(k1 - k2 / val_dat$K + k2 / Kref)))^(1/(1-k3))
     val_predicted = c0 - c0 * val_degrad
     } else if(step1_down_object$user_parameters$reparameterisation == FALSE & step1_down_object$user_parameters$zero_order == FALSE) {  ## Model type 4: no reparameterisation and k3 is not 0
     val_degrad = 1 - ((1 - k3) * (1/(1 - k3) - val_dat$time * exp(k1 - k2 / val_dat$K)))^(1/(1-k3))
     val_predicted = c0 - c0 * val_degrad                    }

     val_residuals = val_dat$y - val_predicted

     val_dat = data.frame(val_dat, residuals = val_residuals, predicted = val_predicted)
     dat = rbind(dat, val_dat)

    shape_types <- c(16,1)
    names(shape_types) <- c("Fit", "Validation")    }

# Other types of residuals for fitted data
  title_addition = ""

  if (residuals == "studentized") {
    # Initialize a vector to store studentized residuals
    studentized_residuals <- numeric(nrow(dat))

    # Loop through each data point
    for (i in seq_len(nrow(dat))) {
      # Exclude the i-th data point
      dat_excluded <- dat[-i, ]

      # Re-fit the model without the i-th data point
      sink(tempfile())
      fit_excluded <- step1_down_basic(
        data = dat_excluded,
        y = "y",
        .time = "time",
        #C = "Celsius",
        K = "K",
        parms = step1_down_object$user_parameters$parms,
        reparameterisation = step1_down_object$user_parameters$reparameterisation,
        zero_order = step1_down_object$user_parameters$zero_order
      )
      sink()

      # find the predicted value when fitted without the point
      k1 = coef(fit_excluded)["k1"]
      k2 = coef(fit_excluded)["k2"]
      k3 = coef(fit_excluded)["k3"]
      c0 = coef(fit_excluded)["c0"]
      if (step1_down_object$user_parameters$reparameterisation == TRUE){
        Kref = mean(dat_excluded$K)
        predicted_i = c0 - c0 * (1 - ((1 - k3) * (1/(1 - k3) - dat$time[i] * exp(k1 - k2/dat$K[i] + k2/Kref)))^(1/(1 - k3)))

        }else{
          predicted_i = c0 - c0 * (1 - ((1 - k3) * (1/(1 - k3) - dat$time[i] * exp(k1 - k2 / dat$K[i])))^(1/(1-k3)))
      }

      # Calculate the residual for the excluded data point
      residual_i <- dat$y[i] - predicted_i

      # Calculate the standard error of the residual
      se_residual <- sqrt(mean((summary(fit_excluded)$residuals)^2))

      # Calculate the studentized residual
      studentized_residuals[i] <- residual_i / se_residual
      }

    # Store the studentized residuals in the data frame
    dat$residuals <- studentized_residuals

    title_addition <- "Studentized "
  }else if(residuals == "standardized"){
    # Extract the fitted values and residuals from the model
    ordinary_residuals <- dat$residuals

    # Calculate the standard error of the residuals
    sigma_squared <- mean(ordinary_residuals^2)
    se_residuals <- sqrt(sigma_squared)

    # Calculate the standardized residuals
    standardized_residuals <- ordinary_residuals / se_residuals

    # Store the standardized residuals in the data frame
    dat$residuals <- standardized_residuals

    title_addition <- "Standardized "

  }


  # Histogram plot (excludes validation if present)
   if(!is.null(validation)) {
   histo_dat = dat[dat$validation == "Fit",]
   } else { histo_dat = dat }

   res_histo = ggplot(histo_dat, aes(x = residuals)) +
    geom_histogram(aes(y = after_stat(density)),
                   breaks = seq(min(histo_dat$residuals), max(histo_dat$residuals), by = (max(histo_dat$residuals) - min(histo_dat$residuals))/bins),
                   colour = "black",
                   fill = "white") +
    stat_function(fun = dnorm, args = list(mean = mean(histo_dat$residuals), sd = sd(histo_dat$residuals)),
                  xlim = c(min(histo_dat$residuals), max(histo_dat$residuals)),
                  col = "turquoise",
                  linewidth = 1,
                  alpha = 0.6) + ggtitle (paste0(title_addition,"Residuals Histogram")) + xlab(paste0(title_addition,"Residuals")) + ylab("Density") +
    mytheme


   # Observed vs predicted (includes validation points if present)
  obs_pred = ggplot() + geom_point(data = dat, mapping = aes(x = predicted, y = y, colour = Celsius, shape = validation)) +
    geom_smooth(data = dat[dat$validation == "Fit",], method ="lm", formula = y ~ x, mapping = aes(x = predicted, y = y)) +
    labs( x = "Predicted response", y = "Observed data")+
    ggtitle ("Observed Vs Predicted") +
    mytheme +
   {if(!is.null(validation))scale_shape_manual(values = shape_types, name = NULL)} +
   theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))


    # Residuals vs predicted (includes validation points if present)
  res_pred = ggplot() + geom_point(data = dat, mapping = aes(x = predicted, y = residuals, colour = Celsius, shape = validation)) +
    labs( x = "Predicted response", y = paste0(title_addition,"Residuals")) +
    geom_hline(yintercept=0, linetype="solid", color = "black")+
    ggtitle (paste0(title_addition,"Residuals Vs Predicted")) +
    mytheme +
   {if(!is.null(validation))scale_shape_manual(values = shape_types, name = NULL)} +
   theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))


   # QQ plot (excludes validation points if present)
  if(!is.null(validation)){
  qqplot <- ggplot() +
  stat_qq_line(data = dat[dat$validation == "Fit",], aes(sample = residuals)) +
  geom_point(data = cbind(validation = dat[dat$validation == "Fit",]$validation, Celsius = dat[dat$validation == "Fit",]$Celsius, as.data.frame(stats::qqnorm(dat[dat$validation == "Fit",]$residuals, plot.it = FALSE))), aes(x = x, y = y, colour = Celsius)) +
  ggtitle ("Q-Q Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")+
  mytheme
  } else {
  qqplot <- ggplot() +
  stat_qq_line(data = dat, aes(sample = residuals)) +
  geom_point(data = cbind(Celsius = dat$Celsius, as.data.frame(stats::qqnorm(dat$residuals, plot.it = FALSE))), aes(x = x, y = y, colour = Celsius)) +
  ggtitle ("Q-Q Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")+
  mytheme  }


   # Plot showing residuals by time (includes validation points if present)
  res_time = ggplot() + geom_point(data = dat, mapping = aes(x = time, y = residuals, colour = Celsius, shape = validation)) +
    labs( x = "Time", y = paste0(title_addition,"Residuals")) +
    geom_hline(yintercept=0, linetype="solid", color = "black")+
    ggtitle (paste0(title_addition,"Residuals By Time")) +
    mytheme +
   {if(!is.null(validation))scale_shape_manual(values = shape_types, name = NULL)} +
   theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))


  results = list(res_histo, qqplot, obs_pred, res_pred, res_time)
  names(results) = c("Residuals_Histogram","Q_Q_Plot", "Observed_V_Predicted","Residuals_V_Predicted","Residuals_By_Time")
  return(results)
}

globalVariables(c('residuals','density','predicted','x'))
