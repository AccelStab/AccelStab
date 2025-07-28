#' @title  Create Diagnostic Plots
#'
#' @description Generate residual diagnostic plots from a step1_down fit.
#'
#' @details Use the fit object obtained from the step1_down or step1_down_batch functions to plot the
#' residual diagnostic plots, assess the quality of fit and search for anomalies. The type of residuals to be assessed
#' can also be selected.
#' Plots created are: Residuals Histogram, QQ Plot of Residuals, Observed Vs Predicted results, Residuals
#'  Vs Predicted results and Residuals By Time.
#' If the step1_down object specified validation data, residuals of these data points are displayed only in the
#' latter three plots. If present in the step1_down object, batches are detected and shown in separate panels per plot.
#'
#' @param step1_down_object The fit object from the step1_down or step1_down_batch functions (required).
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

  dat$batcheff = 0
  
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

     # For validation containing batches, also predict expected batch effects
     if(!is.null(step1_down_object$user_parameters$batch)) {
       if(step1_down_object$user_parameters$zero_order == TRUE){
               n.from = 4}else{n.from=5}
       n.fixeffs =  length(levels(step1_down_object$data$batch)) - 1
       batchcoding = data.frame(
       batch = c(levels(step1_down_object$data$batch)),
       coding = contrasts(step1_down_object$data$batch),
       t(coef(step1_down_object $fit)[n.from:(n.from+n.fixeffs-1)]))
       for( j in 1:n.fixeffs){
                eff <- c(batchcoding[, (j + 1)]) * c(batchcoding[, j + (1 + n.fixeffs) ])
                if (j == 1){   effs = eff
 	              }else{effs = data.frame(cbind(effs, eff)) }  }
       batchcoding$totaleff = rowSums(effs[,1:n.fixeffs])
       val_dat$batcheff = batchcoding$totaleff[match(val_dat$batch, batchcoding$batch)]
                                                              } else {val_dat$batcheff = 0 }
  
     if(step1_down_object$user_parameters$reparameterisation == TRUE & step1_down_object$user_parameters$zero_order == TRUE) {     ## Model type 1: reparameterisation and k3 = 0
     val_degrad = (val_dat$time * exp(k1 - k2/val_dat$K + k2/Kref))
     val_predicted = (c0 + val_dat$batcheff) - (c0 + val_dat$batcheff) * val_degrad
     } else if(step1_down_object$user_parameters$reparameterisation == FALSE & step1_down_object$user_parameters$zero_order == TRUE){  ## Model type 2: no reparameterisation and k3 = 0
     val_degrad = (val_dat$time * exp(k1 - k2 / val_dat$K))
     val_predicted = (c0 + val_dat$batcheff) - (c0 + val_dat$batcheff) * val_degrad
     } else if(step1_down_object$user_parameters$reparameterisation == TRUE & step1_down_object$user_parameters$zero_order == FALSE) {  ## Model type 3: reparameterisation and k3 is not zero
     val_degrad = (1 - ((1 - k3) * (1/(1 - k3) - val_dat$time * exp(k1 - k2 / val_dat$K + k2 / Kref)))^(1/(1-k3)))
     val_predicted = (c0 + val_dat$batcheff) - (c0 + val_dat$batcheff) * val_degrad
     } else if(step1_down_object$user_parameters$reparameterisation == FALSE & step1_down_object$user_parameters$zero_order == FALSE) {  ## Model type 4: no reparameterisation and k3 is not 0
     val_degrad = (1 - ((1 - k3) * (1/(1 - k3) - val_dat$time * exp(k1 - k2 / val_dat$K)))^(1/(1-k3)))
     val_predicted = (c0 + val_dat$batcheff) - (c0 + val_dat$batcheff) * val_degrad                   }

     val_residuals = val_dat$y - val_predicted

     val_dat = data.frame(val_dat, residuals = val_residuals, predicted = val_predicted)
     dat = rbind(dat, val_dat)
     contrasts(dat$batch) = contrasts(step1_down_object$data$batch) # Restore any custom factor coding to batch variable

    shape_types <- c(16,1)
    names(shape_types) <- c("Fit", "Validation")    }

# Other types of residuals for fitted data
  title_addition = ""

# Studentized residuals

 # Studentized residuals when the original model lacked batch effects
  if (residuals == "studentized" & is.null(step1_down_object$user_parameters$batch)) {   
                                    
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
      c0 = coef(fit_excluded)["c0"]
      
   if(step1_down_object$user_parameters$reparameterisation == TRUE & step1_down_object$user_parameters$zero_order == TRUE) {     ## Model type 1: reparameterisation and k3 = 0
        Kref = mean(dat_excluded$K)
 	degrad_i = (dat$time[i] * exp(k1 - k2/dat$K[i] + k2/Kref))
        predicted_i = c0 - c0 * degrad_i
   } else if(step1_down_object$user_parameters$reparameterisation == FALSE & step1_down_object$user_parameters$zero_order == TRUE){  ## Model type 2: no reparameterisation and k3 = 0
       degrad_i = (dat$time[i] * exp(k1 - k2 / dat$K[i]))
       predicted_i = c0 - c0 * degrad_i
    } else if(step1_down_object$user_parameters$reparameterisation == TRUE & step1_down_object$user_parameters$zero_order == FALSE) {  ## Model type 3: reparameterisation and k3 is not zero
     Kref = mean(dat_excluded$K)
     k3 = coef(fit_excluded)["k3"] 
     degrad_i = (1 - ((1 - k3) * (1/(1 - k3) - dat$time[i] * exp(k1 - k2 / dat$K[i] + k2 / Kref)))^(1/(1-k3)))
     predicted_i = c0 - c0 * degrad_i
  } else if(step1_down_object$user_parameters$reparameterisation == FALSE & step1_down_object$user_parameters$zero_order == FALSE) {  ## Model type 4: no reparameterisation and k3 is not 0
     k3 = coef(fit_excluded)["k3"]
     degrad_i = (1 - ((1 - k3) * (1/(1 - k3) - dat$time[i] * exp(k1 - k2 / dat$K[i])))^(1/(1-k3)))
     predicted_i = c0 - c0 * degrad_i                   }

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
 
# Studentized residuals when the original model included batch effects
  }else if(residuals == "studentized" & !is.null(step1_down_object$user_parameters$batch)) {   

 # Initialize a vector to store studentized residuals
    studentized_residuals <- numeric(nrow(dat))

    # Loop through each data point
    for (i in seq_len(nrow(dat))) {
 
     # Exclude the i-th data point
      dat_excluded <- dat[-i, ]

      # Re-fit the model without the i-th data point
      sink(tempfile())
      fit_excluded <- step1_down_basic_batch(
        data = dat_excluded,
        y = "y",
        .time = "time",
        #C = "Celsius",
        K = "K",
        batch = "batch",
        parms = step1_down_object$user_parameters$parms,
        reparameterisation = step1_down_object$user_parameters$reparameterisation,
        zero_order = step1_down_object$user_parameters$zero_order
      )
      sink()

      # find the predicted value when fitted without the point
      k1 = coef(fit_excluded)["k1"]
      k2 = coef(fit_excluded)["k2"]
      c0 = coef(fit_excluded)["c0"]

  if(step1_down_object$user_parameters$zero_order == TRUE){
     n.from = 4 } else {
     n.from = 5
     k3 = coef(fit_excluded)["k3"] }    
  n.fixeffs =  length(levels(step1_down_object$data$batch)) - 1
       batchcoding = data.frame(
       batch = c(levels(step1_down_object$data$batch)),
       coding = contrasts(step1_down_object$data$batch),
       t(coef(fit_excluded)[n.from:(n.from+n.fixeffs-1)]))
       for( j in 1:n.fixeffs){
                eff <- c(batchcoding[, (j + 1)]) * c(batchcoding[, j + (1 + n.fixeffs) ])
                if (j == 1){   effs = eff
 	              }else{effs = data.frame(cbind(effs, eff)) }  }
       batchcoding$totaleff = rowSums(effs[,1:n.fixeffs])
       dat$batcheff = batchcoding$totaleff[match(dat$batch, batchcoding$batch)]
	
   if(step1_down_object$user_parameters$reparameterisation == TRUE & step1_down_object$user_parameters$zero_order == TRUE) {     ## Model type 1: reparameterisation and k3 = 0
        Kref = mean(dat_excluded$K)
 	degrad_i = (dat$time[i] * exp(k1 - k2/dat$K[i] + k2/Kref))
        predicted_i = (c0 + dat$batcheff[i]) - (c0 + dat$batcheff[i]) * degrad_i
   } else if(step1_down_object$user_parameters$reparameterisation == FALSE & step1_down_object$user_parameters$zero_order == TRUE){  ## Model type 2: no reparameterisation and k3 = 0
       degrad_i = (dat$time[i] * exp(k1 - k2 / dat$K[i]))
        predicted_i = (c0 + dat$batcheff[i]) - (c0 + dat$batcheff[i]) * degrad_i
    } else if(step1_down_object$user_parameters$reparameterisation == TRUE & step1_down_object$user_parameters$zero_order == FALSE) {  ## Model type 3: reparameterisation and k3 is not zero
     Kref = mean(dat_excluded$K)
     k3 = coef(fit_excluded)["k3"] 
     degrad_i = (1 - ((1 - k3) * (1/(1 - k3) - dat$time[i] * exp(k1 - k2 / dat$K[i] + k2 / Kref)))^(1/(1-k3)))
     predicted_i = (c0 + dat$batcheff[i]) - (c0 + dat$batcheff[i]) * degrad_i
  } else if(step1_down_object$user_parameters$reparameterisation == FALSE & step1_down_object$user_parameters$zero_order == FALSE) {  ## Model type 4: no reparameterisation and k3 is not 0
     k3 = coef(fit_excluded)["k3"]
     degrad_i = (1 - ((1 - k3) * (1/(1 - k3) - dat$time[i] * exp(k1 - k2 / dat$K[i])))^(1/(1-k3)))
     predicted_i = (c0 + dat$batcheff[i]) - (c0 + dat$batcheff[i]) * degrad_i     }

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

    } else if(residuals == "standardized"){
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
  if(!is.null(step1_down_object$user_parameters$batch)){       # Panels for individual batches if present in the model
      res_histo =  res_histo + facet_wrap(~batch) }

   # Observed vs predicted (includes validation points if present)
  obs_pred = ggplot() + geom_point(data = dat, mapping = aes(x = predicted, y = y, colour = Celsius, shape = validation)) +
    geom_smooth(data = dat[dat$validation == "Fit",], method ="lm", formula = y ~ x, mapping = aes(x = predicted, y = y)) +
    labs( x = "Predicted response", y = "Observed data")+
    ggtitle ("Observed Vs Predicted") +
    mytheme +
   {if(!is.null(validation))scale_shape_manual(values = shape_types, name = NULL)} +
   theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))
  if(!is.null(step1_down_object$user_parameters$batch)){       # Panels for individual batches if present in the model
      obs_pred =  obs_pred + facet_wrap(~batch) }

    # Residuals vs predicted (includes validation points if present)
  res_pred = ggplot() + geom_point(data = dat, mapping = aes(x = predicted, y = residuals, colour = Celsius, shape = validation)) +
    labs( x = "Predicted response", y = paste0(title_addition,"Residuals")) +
    geom_hline(yintercept=0, linetype="solid", color = "black")+
    ggtitle (paste0(title_addition,"Residuals Vs Predicted")) +
    mytheme +
   {if(!is.null(validation))scale_shape_manual(values = shape_types, name = NULL)} +
   theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))
  if(!is.null(step1_down_object$user_parameters$batch)){       # Panels for individual batches if present in the model
      res_pred =  res_pred + facet_wrap(~batch) }

   # QQ plot (excludes validation points if present)
   if(!is.null(validation)) {
   qq_dat = dat[dat$validation == "Fit",]
   } else { qq_dat = dat }  

  qqs = as.data.frame(stats::qqnorm(qq_dat$residuals, plot.it = FALSE))
  qq_dat = data.frame(qq_dat, x.n = qqs$x, y.n = qqs$y)

  qqplot <- ggplot(data = qq_dat) +
  stat_qq_line(aes(sample = residuals)) +
  geom_point(aes(x = x.n, y = y.n, colour = Celsius)) +
  ggtitle ("Q-Q Plot") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")+
  mytheme
      if(!is.null(step1_down_object$user_parameters$batch)){       # Panels for individual batches if present in the model
      qqplot =  qqplot + facet_wrap(~batch) }
  
   # Plot showing residuals by time (includes validation points if present)
  res_time = ggplot() + geom_point(data = dat, mapping = aes(x = time, y = residuals, colour = Celsius, shape = validation)) +
    labs( x = "Time", y = paste0(title_addition,"Residuals")) +
    geom_hline(yintercept=0, linetype="solid", color = "black")+
    ggtitle (paste0(title_addition,"Residuals By Time")) +
    mytheme +
   {if(!is.null(validation))scale_shape_manual(values = shape_types, name = NULL)} +
   theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))
     if(!is.null(step1_down_object$user_parameters$batch)){       # Panels for individual batches if present in the model
      res_time =  res_time + facet_wrap(~batch) }

  results = list(res_histo, qqplot, obs_pred, res_pred, res_time)
  names(results) = c("Residuals_Histogram","Q_Q_Plot", "Observed_V_Predicted","Residuals_V_Predicted","Residuals_By_Time")
  return(results)
closeAllConnections()
}

globalVariables(c('residuals','density','predicted','x'))
