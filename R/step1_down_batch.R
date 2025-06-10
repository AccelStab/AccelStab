#' @title Step1 Down Model With Batch Effects
#'
#' @description Fit the one-step Šesták–Berggren kinetic model including batch effects.
#'
#' @details Fit the one-step Šesták–Berggren kinetic (non-linear) model using
#' accelerated stability and batch data that has been stored in an R data frame. Additionally,
#' predictions of the mean at each tested temperature are returned for each batch, including associated
#' confidence and prediction intervals, which can be subsequently visualised with
#' step1_plot_pred(), step1_plot_CI(), step1_plot_PI() and step1_plot_T(). Kinetic
#' parameters (k1, k2 and, if used, k3) are retained in the model even if one or more of
#' these parameters turn out to be non-significant. Further arguments relating to
#' model fitting, such as setting lower bounds for one or more model parameters,
#' may be passed.
#'
#' @param data Dataframe containing accelerated stability data and batches (required).
#' @param y Name of decreasing variable (e.g. concentration) contained within data
#'  (required).
#' @param .time Time variable contained within data (required).
#' @param batch Batch variable denoting individual batch or lot IDs; factor coding may be used and will be recognised (required).
#' @param K Kelvin variable (numeric or column name) (optional).
#' @param C Celsius variable (numeric or column name) (optional).
#' @param validation Validation dummy variable, the column must contain only
#'  1s and 0s, 1 for validation data and 0 for fit data. (column name) (optional).
#' @param draw Number of simulations used to estimate confidence intervals.
#'  When set to NULL the calculus method is used, however this is not recommended.
#' @param parms Starting values for the parameters as a list - k1, k2, k3, and c0.
#' @param temp_pred_C Integer or numeric value to predict the response for a
#'  given temperature (in Celsius).
#' @param max_time_pred Maximum time to predict the response variable.
#' @param confidence_interval Confidence level for the confidence and prediction intervals
#'  around the predictions (default 0.95).
#' @param by Number of points (on the time scale) to smooth the statistical
#'  intervals around the predictions.
#' @param reparameterisation Use alternative parameterisation of the one-step
#'  model which aims to reduce correlation between k1 and k2.
#' @param zero_order Set kinetic order, k3, to zero (straight lines).
#' @param ... Further arguments to passed to minpack.lm.
#'
#' @return An SB class object, a list including the following elements:
#' \itemize{
#'  \item *fit* - The non-linear fit.
#'  \item *data* - The data set.
#'  \item *prediction* - A data frame containing the predictions with the confidence and prediction intervals.
#'  \item *user_parameters* - List of users input parameters which is utilised by other
#'    functions in the package.
#'  \item *sample_coefficients* - A matrix containing the coefficients sampled during bootstrapping.
#'    }
#'
#' @examples #Create a df with 3 batches
#' [Further example text]
#' [Further example text]
#' [Further example text]
#'
#' @importFrom stats vcov coef runif confint rnorm rchisq quantile qt complete.cases
#' @importFrom minpack.lm nls.lm
#'
#' @export step1_down

step1_down_batch <- function (data, y, .time, batch, K = NULL, C = NULL, validation = NULL,
                        draw = 10000, parms = NULL, temp_pred_C = NULL,
                        max_time_pred = NULL, confidence_interval = 0.95, by = 101,
                        reparameterisation = FALSE, zero_order = FALSE, ...){

 if (is.null(batch))
    stop("Select the batch or lot variable")

  if (is.null(K) & is.null(C))
    stop("Select the temperature variable in Kelvin or Celsius")

  if (!is.null(parms) & !is.list(parms))
    stop("The starting values for parameters must be a list, or keep as NULL")

  if (!is.null(validation))
    if (!all(data[,validation] %in% c(0,1)))
      stop("Validation column must contain 1s and 0s only")

  user_parameters <- list(
    data = data, y = y, .time = .time, batch = batch, K = K, C = C, validation = validation,draw = draw,
    parms = parms, temp_pred_C = temp_pred_C, max_time_pred = max_time_pred,
    confidence_interval = confidence_interval, by = by,
    reparameterisation = reparameterisation, zero_order = zero_order)

  ## Additional arguments in the call will be passed to model fitting with minpack.lm
  minpack_args = list(...)

  ## Temperature: both C and K are provided
  if(!is.null(C) & !is.null(K)) {
    data[, C] <- ifelse(is.na(data[, C]) & !is.na(data[, K]),
                        data[, K] - 273.15,
                        data[, C])
    data[, K] <- ifelse(is.na(data[, K]) & !is.na(data[, C]),
                        data[, C] + 273.15,
                        data[, K])   }

  ## Temperature: only C or only K is provided
  if (!is.null(C) & is.null(K)) {
   K = 'K'
   data[, K] = data[, C] + 273.15  }
  if (!is.null(K) & is.null(C)) {
    C = 'C'
    data[, C] = data[, K] - 273.15 }

  data <- data[complete.cases(data[, c(C,K,y,.time, batch)]), ]

  dat = data
  dat$K = dat[, K]
  Kref = mean(dat$K)
  dat$Celsius = as.factor(dat[, C])
  dat$time = dat[, .time]
  dat$y = dat[, y]
  if (is.factor(dat[, batch])) {
			dat$batch = dat[, batch] } else {
			dat$batch = factor(dat[, batch])}
  if(!is.null(validation)){
    dat$validation = ifelse(dat[,validation] == 0, "Fit", "Validation")
    if(validation != "validation"){
      dat <- dat[, !names(dat) %in% c(validation)]
    }
  }
  if(.time != "time"){
    dat <- dat[, !names(dat) %in% c(.time)]
  }
  if(y != "y"){
    dat <- dat[, !names(dat) %in% c(y)]
  }
  if (batch != "batch") {
    dat <- dat[, !names(dat) %in% c(batch)]
                         }

  Temps = sort(unique(dat$K))
  batches = unique(dat$batch)
  no_batches <- length(batches)
  no_fixed_effect <- no_batches - 1
  batch_names = as.data.frame(levels(dat$batch))
  batch_coding = as.data.frame(contrasts(dat$batch))
  
  if (!is.null(temp_pred_C))
    Temps = unique(sort(c(Temps, temp_pred_C + 273.15)))
  if (is.null(max_time_pred))
    max_time_pred = max(dat$time, na.rm = TRUE)
  times.pred = seq(0, max_time_pred, length.out = by)

  dat_full <- dat
  if(!is.null(validation)){
    dat <- dat[dat$validation == "Fit",]
  }

  if(is.null(parms)){
    sorted_data <- dat[order(dat$time), ]

    min_time <- min(sorted_data$time)

    if (sum(sorted_data$time == min_time) > 3) {
      selected_rows <- sorted_data$time == min_time
    } else {
      selected_rows <- seq_len(min(3, nrow(sorted_data)))
    }
    c0_initial <- mean(sorted_data$y[selected_rows])
  }

## Model type 1 - reparameterisation and k3 = 0
  if(reparameterisation & zero_order){

## Print a message informing lower bounds = 0 may not be suitable with the reparamerised version
cat("The alternative parameterisation of the one-step model was used. Note that the lower bounds for all parameters are set to 0 unless other lower bounds are specified in step1_down() or step1_down_basic().\n\n")

  MyFctNL = function(parms) { # Make function
      k1 = parms$k1
      k2 = parms$k2
      c0 = parms$c0
       
      for (i in 1:no_fixed_effect) {
        var_name <- paste0("b", i)
        assign(var_name, as.numeric(parms[i + 3]))
                                    }
      
      degrad = dat$time * exp(k1 - k2/dat$K + k2/Kref) 
      effs = paste0("(c0 + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + ") , ") - (c0 + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + "), ") * degrad")
      test = rep("", no_batches)
      
      for (i in 1: no_batches-1) {
        test[i] = paste0("ifelse(dat$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
        for (j in no_batches: no_batches) {
        test[j] = gsub("bn", j, effs)
                                          }
                                 }

      test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
    
      residual = dat$y - eval(parse(text = test))
      return(residual)
                                  }

  if (!"fn" %in% names(minpack_args)) {
    minpack_args$fn =  MyFctNL    }

    # Fit model:
    if (!is.null(parms)) {
      minpack_args$par =  parms
    if (!"lower" %in% names(minpack_args)) 	{
    minpack_args$lower =  c(rep(0, 3), rep(-Inf, no_fixed_effect))    }
	if(length(minpack_args$par) != length(minpack_args$lower))
	stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")

  fit = do.call(minpack.lm::nls.lm, minpack_args)
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 40), k2 = stats::runif(1, 1000, 20000), c0 = c0_initial)

        for (i in 1:no_fixed_effect) {
          var_name <- paste0("b", i)
          parms[[var_name]] <- rnorm(1, 0, 2)
        }

minpack_args$par = parms

if (!"lower" %in% names(minpack_args)) 	{	##
	minpack_args$lower =  c(rep(0, 3), rep(-Inf, no_fixed_effect))    }

	if(length(minpack_args$par) != length(minpack_args$lower))
	stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")

  fit = suppressWarnings(do.call(minpack.lm::nls.lm, minpack_args))

   fit <- tryCatch({
          suppressWarnings(do.call(minpack.lm::nls.lm, minpack_args))
        },
        error = function(e){"error"},
        warning = function(w){"warning"})

        vcov_test <- tryCatch({
          stats::vcov(fit)
        },
        error = function(e){"error"},
        warning = function(w){"warning"})

        if(all(!(fit %in% c("error","warning"))) && all(!(vcov_test %in% c("error","warning", NaN)))){
          break
        }
      }
      fit = do.call(minpack.lm::nls.lm, minpack_args)

 }

    # Generate the predictions
    k1 = stats::coef(fit)[1]
    k2 = stats::coef(fit)[2]
    c0 = stats::coef(fit)[3]

   for (i in 1:no_fixed_effect) {
      var_name <- paste0("b", i)
      assign(var_name, fit$par[[var_name]])##
    }

    SIG = stats::vcov(fit)
    sigma = summary(fit)$sigma
    DF = summary(fit)$df[2]
    n.params = summary(fit)$df[1]

    pred = expand.grid(time = times.pred, K = Temps, batch = batches)
    pred$Degradation = pred$time * exp(k1 - k2/pred$K + k2/Kref)
 
    effs = paste0("(c0 + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + ") , ") - (c0 + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + "), ") * pred$Degradation"       )
    test = rep("", no_batches)
    for (i in 1: no_batches-1) {
    test[i] = paste0("ifelse(pred$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
    for (j in no_batches: no_batches) {
    test[j] = gsub("bn", j, effs)
      }
                                                   }
    test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
    pred$Response = eval(parse(text = test))

# Function to be used in drawing from multi-t
pred_fct = function(coef.fit) {
  degrad = pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K + coef.fit[2] / Kref)
  effs = paste0("(coef.fit[3] + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * coef.fit[", 4:(no_fixed_effect+3),"]", collapse = " + ") , ") - (coef.fit[3] + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * coef.fit[", 4:(no_fixed_effect+3), "]", collapse = " + "), ") * degrad"       )
  test = rep("", no_batches)
  for (i in 1: no_batches-1) {
    test[i] = paste0("ifelse(pred$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
    for (j in no_batches: no_batches) {
      test[j] = gsub("bn", j, effs)
                                       } }
  test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
  conc =  eval(parse(text = test))
  return(conc)
                                }

      # Multi T samples
      rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
      rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)
      res.draw = matrix(nrow = draw, ncol = nrow(pred), byrow = TRUE, apply(rand.coef, 1, pred_fct))

      CI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      CI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

      sigma.dist = sqrt((DF * sigma^2) / rchisq(n = draw*length(pred$time), df = DF))
      res.draw = res.draw + rnorm(n = draw*length(pred$time), mean = 0, sd = sigma.dist)

      PI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      PI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)


## Model type 2 - no reparameterisation and k3 = 0
  }else if(!reparameterisation & zero_order){
  MyFctNL = function(parms) { # Make function
      k1 = parms$k1
      k2 = parms$k2
      c0 = parms$c0
       
      for (i in 1:no_fixed_effect) {
        var_name <- paste0("b", i)
        assign(var_name, as.numeric(parms[i + 3]))
                                    }
      
      degrad = dat$time * exp(k1 - k2/dat$K) 

      effs = paste0("(c0 + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + ") , ") - (c0 + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + "), ") * degrad")
      test = rep("", no_batches)
      
      for (i in 1: no_batches-1) {
        test[i] = paste0("ifelse(dat$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
        for (j in no_batches: no_batches) {
        test[j] = gsub("bn", j, effs)
                                          }
                                 }

      test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
    
      residual = dat$y - eval(parse(text = test))
      return(residual)
                                  }

  if (!"fn" %in% names(minpack_args)) 	{
    minpack_args$fn =  MyFctNL    }

    # Fit model
    if (!is.null(parms)) {
    minpack_args$par =  parms
    if (!"lower" %in% names(minpack_args)) 	{
    minpack_args$lower =  c(rep(0, 3), rep(-Inf, no_fixed_effect))   }
	if(length(minpack_args$par) != length(minpack_args$lower))
	stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")

      fit = do.call(minpack.lm::nls.lm, minpack_args)
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 40), k2 = stats::runif(1, 1000, 20000), c0 = c0_initial)

        for (i in 1:no_fixed_effect) {
          var_name <- paste0("b", i)
          parms[[var_name]] <- rnorm(1, 0, 2)##
        }

  minpack_args$par = parms ##

if (!"lower" %in% names(minpack_args)) 	{
	  minpack_args$lower =  c(rep(0, 3), rep(-Inf, no_fixed_effect))    }

	if(length(minpack_args$par) != length(minpack_args$lower))
  stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")

 fit <- tryCatch({
          suppressWarnings(do.call(minpack.lm::nls.lm, minpack_args))
        },
        error = function(e){"error"},
        warning = function(w){"warning"})

        vcov_test <- tryCatch({
          stats::vcov(fit)
        },
        error = function(e){"error"},
        warning = function(w){"warning"})

        if(all(!(fit %in% c("error","warning"))) && all(!(vcov_test %in% c("error","warning", NaN)))){
          break
        }
      }
      fit = do.call(minpack.lm::nls.lm, minpack_args)
    }

   # Calculate the predictions
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    c0 = coef(fit)[3]
   
    for (i in 1:no_fixed_effect) {
      var_name <- paste0("b", i)
      assign(var_name, fit$par[[var_name]])
    }

    SIG = vcov(fit)
    sigma = summary(fit)$sigma
    DF = summary(fit)$df[2]
    n.params = summary(fit)$df[1]  

    pred = expand.grid(time = times.pred, K = Temps, batch = batches)
    pred$Degradation = pred$time * exp(k1 - k2 / pred$K)
    
    effs = paste0("(c0 + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + ") , ") - (c0 + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + "), ") * pred$Degradation"       )
    test = rep("", no_batches)
    for (i in 1: no_batches-1) {
    test[i] = paste0("ifelse(pred$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
    for (j in no_batches: no_batches) {
    test[j] = gsub("bn", j, effs)
      }
                                                   }
    test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
    pred$Response = eval(parse(text = test))

# Function to be used in drawing from multi-T
pred_fct = function(coef.fit) {
  degrad = pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K)
  effs = paste0("(coef.fit[3] + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * coef.fit[", 4:(no_fixed_effect+3),"]", collapse = " + ") , ") - (coef.fit[3] + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * coef.fit[", 4:(no_fixed_effect+3), "]", collapse = " + "), ") * degrad"       )
  test = rep("", no_batches)
  for (i in 1: no_batches-1) {
    test[i] = paste0("ifelse(pred$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
    for (j in no_batches: no_batches) {
      test[j] = gsub("bn", j, effs)
                                       } }
  test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
  conc =  eval(parse(text = test))
  return(conc)
                                }

      # Multi T samples
      rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
      rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)
      res.draw = matrix(nrow = draw, ncol = nrow(pred), byrow = TRUE, apply(rand.coef, 1, pred_fct))

      CI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      CI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

      sigma.dist = sqrt((DF * sigma^2) / rchisq(n = draw*length(pred$time), df = DF))
      res.draw = res.draw + rnorm(n = draw*length(pred$time), mean = 0, sd = sigma.dist)
      PI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      PI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)


## Model type 3 - reparameterisation and k3 is not zero
  }else if(reparameterisation & !zero_order){

  ## Print a message informing lower bounds = 0 may not be suitable with the reparamerised version
  cat("The alternative parameterisation of the one-step model was used. Note that the lower bounds for all parameters are set to 0 unless other lower bounds are specified in step1_down() or step1_down_basic().\n\n")

   MyFctNL = function(parms) {
      k1 = parms$k1
      k2 = parms$k2
      k3 = parms$k3
      c0 = parms$c0

      for (i in 1:no_fixed_effect) {
        var_name <- paste0("b", i)
        assign(var_name, as.numeric(parms[i + 4]))
                                    }
      
      degrad = (1 - ((1 - k3) * (1/(1 - k3) - dat$time * exp(k1 - k2/dat$K + k2/Kref)))^(1/(1 - k3)))
      
      effs = paste0("(c0 + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + ") , ") - (c0 + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + "), ") * degrad"       )

      test = rep("", no_batches)
            for (i in 1: no_batches-1) {
              test[i] = paste0("ifelse(dat$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
                 for (j in no_batches: no_batches) {
                     test[j] = gsub("bn", j, effs)
                                                   }
                                        }

      test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
      residual = dat$y - eval(parse(text = test))##
      return(residual)
    }

  if (!"fn" %in% names(minpack_args)) 	{
    minpack_args$fn =  MyFctNL    }

  # Fit model:
  if (!is.null(parms)) { 
        minpack_args$par =  parms
  if (!"lower" %in% names(minpack_args)) 	{
    minpack_args$lower =  c(rep(0, 4), rep(-Inf, no_fixed_effect))   }

	if(length(minpack_args$par) != length(minpack_args$lower))
        stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")

      fit = do.call(minpack.lm::nls.lm, minpack_args)
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 60), k2 = stats::runif(1, 1000, 20000), k3 = stats::runif(1, 0, 11), c0 = c0_initial)

      for (i in 1:no_fixed_effect) {
          var_name <- paste0("b", i)
          parms[[var_name]] <- rnorm(1, 0, 2)
                                    }

  minpack_args$par = parms

  if (!"lower" %in% names(minpack_args)) 	{
	    minpack_args$lower =  c(rep(0, 4), rep(-Inf, no_fixed_effect))  }

	if(length(minpack_args$par) != length(minpack_args$lower))
        stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")

   fit <- tryCatch({
          suppressWarnings(do.call(minpack.lm::nls.lm, minpack_args))
        },
        error = function(e){"error"},
        warning = function(w){"warning"})

        vcov_test <- tryCatch({
          stats::vcov(fit)
        },
        error = function(e){"error"},
        warning = function(w){"warning"})

        if(all(!(fit %in% c("error","warning"))) && all(!(vcov_test %in% c("error","warning", NaN)))){
          break
        }
      }
      fit = do.call(minpack.lm::nls.lm, minpack_args)
 }

    # Calculate the predictions
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    k3 = coef(fit)[3]
    if (k3 == 0){cat(paste("k3 is fitted to be exactly 0, we strongly suggest using option zero_order = TRUE","The model will continue with k3 = 0, so degradation is linear over time"," "," ", sep = "\n"))
    }else if(confint(fit,'k3')[1] < 0 && confint(fit,'k3')[2] > 0){print(paste0("The 95% Wald Confidence Interval for k3 includes 0, k3 is estimated as ",signif(k3,4),". We suggest considering option zero_order = TRUE"))}
    c0 = coef(fit)[4]

    for (i in 1:no_fixed_effect) {
      var_name <- paste0("b", i)
      assign(var_name, fit$par[[var_name]])##
    }

    SIG = vcov(fit)
    sigma = summary(fit)$sigma
    DF = summary(fit)$df[2]
    n.params = summary(fit)$df[1]

    pred = expand.grid(time = times.pred, K = Temps, batch = batches)
    pred$Degradation = 1 - ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2 / pred$K + k2 / Kref)))^(1/(1-k3))
    
## Predict response by batch  
effs = paste0("(c0 + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + ") , ") - (c0 + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + "), ") * pred$Degradation"       )

test = rep("", no_batches)
for (i in 1: no_batches-1) {
test[i] = paste0("ifelse(pred$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
for (j in no_batches: no_batches) {
test[j] = gsub("bn", j, effs)
                                  }
                           }

test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")

pred$Response = eval(parse(text = test))

# Function to be used in drawing from multi-T
      pred_fct = function(coef.fit){
          degrad = 1 - ((1 - coef.fit[3]) * (1/(1 - coef.fit[3]) - pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K + coef.fit[2] / Kref)))^(1/(1-coef.fit[3]))
          effs = paste0("(coef.fit[4] + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * coef.fit[", 5:(no_fixed_effect+4),"]", collapse = " + ") , ") - (coef.fit[4] + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * coef.fit[", 5:(no_fixed_effect+4), "]", collapse = " + "), ") * degrad"       )
          test = rep("", no_batches)
          for (i in 1: no_batches-1) {
            test[i] = paste0("ifelse(pred$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
            for (j in no_batches: no_batches) {
                 test[j] = gsub("bn", j, effs)
                                              }
                                      }
         test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
         conc =  eval(parse(text = test))
         return(conc)
                                   }

      # Multi T samples
      rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
      rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)
      res.draw = matrix(nrow = draw, ncol = nrow(pred), byrow = TRUE, apply(rand.coef, 1, pred_fct))

      no_k3_below0 <- sum(rand.coef[,3] < 0)
      if(no_k3_below0 > 0.5){
        cat(paste(paste0(no_k3_below0*100/draw, "% of the draws for k3 are below zero, this might have an adverse effect on the confidence interval, particularly if this value exceeds the confidence %."),"We suggest considering option zero_order = TRUE", sep = "\n"))
                             }

      CI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      CI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

      sigma.dist = sqrt((DF * sigma^2) / rchisq(n = draw*length(pred$time), df = DF))
      res.draw = res.draw + rnorm(n = draw*length(pred$time), mean = 0, sd = sigma.dist)
      PI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      PI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)


## Model type 4 - no reparameterisation and k3 is not 0
  }else if(!reparameterisation & !zero_order){
    MyFctNL = function(parms) {
      k1 = parms$k1
      k2 = parms$k2
      k3 = parms$k3
      c0 = parms$c0

      for (i in 1:no_fixed_effect) {
        var_name <- paste0("b", i)
        assign(var_name, as.numeric(parms[i + 4]))
                                    }
     degrad = (1 - ((1 - k3) * (1/(1 - k3) - dat$time * exp(k1 - k2 / dat$K)))^(1/(1-k3)))##
     effs = paste0("(c0 + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + ") , ") - (c0 + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + "), ") * degrad")
     test = rep("", no_batches)
     for (i in 1: no_batches-1) {
        test[i] = paste0("ifelse(dat$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
        for (j in no_batches: no_batches) {
            test[j] = gsub("bn", j, effs)
                                           }
                                 }
     test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
     residual = dat$y - eval(parse(text = test))
    return(residual)
                            }

  if (!"fn" %in% names(minpack_args)) 	{
    minpack_args$fn =  MyFctNL    }

# Fitting the model:
    if (!is.null(parms)) { 
    minpack_args$par =  parms
  if (!"lower" %in% names(minpack_args)) 	{
    minpack_args$lower =  c(rep(0, 4), rep(-Inf, no_fixed_effect))    }
	if(length(minpack_args$par) != length(minpack_args$lower))
        stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")

      fit = do.call(minpack.lm::nls.lm, minpack_args)
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 60), k2 = stats::runif(1, 1000, 20000), k3 = stats::runif(1, 0, 11), c0 = c0_initial)

       for (i in 1:no_fixed_effect) {
          var_name <- paste0("b", i)
          parms[[var_name]] <- rnorm(1, 0, 2)
        }

  minpack_args$par = parms

  if (!"lower" %in% names(minpack_args)) 	{
          minpack_args$lower =  c(rep(0, 4), rep(-Inf, no_fixed_effect))    }

  if(length(minpack_args$par) != length(minpack_args$lower))
  stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")

 fit <- tryCatch({
          suppressWarnings(do.call(minpack.lm::nls.lm, minpack_args))
        },
        error = function(e){"error"},
        warning = function(w){"warning"})

        vcov_test <- tryCatch({
          stats::vcov(fit)
        },
        error = function(e){"error"},
        warning = function(w){"warning"})

        if(all(!(fit %in% c("error","warning"))) && all(!(vcov_test %in% c("error","warning", NaN)))){
          break
        }
      }
      fit = do.call(minpack.lm::nls.lm, minpack_args)
    }

    # Calculate the predictions
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    k3 = coef(fit)[3]
    if (k3 == 0){cat(paste("k3 is fitted to be exactly 0, we strongly suggest using option zero_order = TRUE","The model will continue with k3 = 0, so degradation is linear over time"," ", " ", sep = "\n"))
    }else if(confint(fit,'k3')[1] < 0 && confint(fit,'k3')[2] > 0){print(paste0("The 95% Wald Confidence Interval for k3 includes 0, k3 is estimated as ",signif(k3,4),". We suggest considering option zero_order = TRUE"))}
    c0 = coef(fit)[4]
    
   for (i in 1:no_fixed_effect) {
      var_name <- paste0("b", i)
      assign(var_name, fit$par[[var_name]])
    }

     SIG = stats::vcov(fit)
    sigma = summary(fit)$sigma
    DF = summary(fit)$df[2]
    n.params = summary(fit)$df[1]

    pred = expand.grid(time = times.pred, K = Temps, batch = batches)
    pred$Degradation = 1 - ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2 / pred$K)))^(1/(1-k3))
   
    effs = paste0("(c0 + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + ") , ") - (c0 + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * b", 1:no_fixed_effect, collapse = " + "), ") * pred$Degradation"       )
    test = rep("", no_batches)
    for (i in 1: no_batches-1) {
       test[i] = paste0("ifelse(pred$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
          for (j in no_batches: no_batches) {
              test[j] = gsub("bn", j, effs)
                                             }
                                }
    test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
    pred$Response = eval(parse(text = test))##

# Function to be used in drawing from multi-T
      pred_fct = function(coef.fit){
        degrad = 1 - ((1 - coef.fit[3]) * (1/(1 - coef.fit[3]) - pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K)))^(1/(1-coef.fit[3]))
        effs = paste0("(coef.fit[4] + ", paste0("batch_coding[bn,", 1:no_fixed_effect,"] * coef.fit[", 5:(no_fixed_effect+4),"]", collapse = " + ") , ") - (coef.fit[4] + ",  paste0("batch_coding[bn,", 1:no_fixed_effect,"] * coef.fit[", 5:(no_fixed_effect+4), "]", collapse = " + "), ") * degrad"       )
        test = rep("", no_batches)
        for (i in 1: no_batches-1) {
             test[i] = paste0("ifelse(pred$batch == batch_names[", gsub("bn", i, paste0("bn, 1],", effs)),",")
             for (j in no_batches: no_batches) {
                   test[j] = gsub("bn", j, effs)
                                                }
                                    }
       test = paste(c(test, paste(rep(")",no_fixed_effect))), collapse = "")
       conc =  eval(parse(text = test))
       return(conc)
                                  }

      # Multi T samples
      rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
      rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)
      res.draw = matrix(nrow = draw, ncol = nrow(pred), byrow = TRUE, apply(rand.coef, 1, pred_fct))

      no_k3_below0 <- sum(rand.coef[,3] < 0)
      if(no_k3_below0 > 0.5){
        cat(paste(paste0(no_k3_below0*100/draw, "% of the draws for k3 are below zero, this might have an adverse effect on the confidence interval, particularly if this value exceeds the confidence %."),"We suggest considering option zero_order = TRUE", sep = "\n"))
      }

      CI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      CI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

      sigma.dist = sqrt((DF * sigma^2) / rchisq(n = draw*length(pred$time), df = DF))
      res.draw = res.draw + rnorm(n = draw*length(pred$time), mean = 0, sd = sigma.dist)
      PI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      PI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)
  }

  pred$Celsius = as.factor(pred$K - 273.15)
  pred$K = as.factor(pred$K)
  pred$fit = "Prediction"
  pred$CI = paste(100*confidence_interval, "% CI")
  pred$PI = paste(100*confidence_interval, "% PI")

   batch_effs = rep(NA, no_fixed_effect)
     for (i in 1:no_fixed_effect) {
      batch_effs[i] <- paste0("b", i)  }

    pred$CI1 = CI1b
    pred$CI2 = CI2b
    pred$PI1 = PI1b
    pred$PI2 = PI2b

    if(zero_order){
      colnames(rand.coef) <- c("k1","k2","c0", batch_effs)
                  }else{
      colnames(rand.coef) <- c("k1","k2","k3","c0", batch_effs)
                        }

  results = list(fit, dat_full, pred,user_parameters, rand.coef)
  names(results) = c("fit", "data", "prediction","user_parameters","sample_coefficients")
  class(results) = "SB"
  return(results)
  }
