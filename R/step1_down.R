#' @title Step1 Down Model
#'
#' @description Fit the one-step Šesták–Berggren kinetic model.
#'
#' @details Fit the one-step Šesták–Berggren kinetic (non-linear) model using
#'  accelerated stability data from an R dataframe format. Parameters are kept
#'  in even when not significant. Predictions, including confidence and prediction
#'  intervals, are also computed. Further arguments relating to model fitting,
#'  such as parameter lower bounds, may be passed.
#'
#' @param data Dataframe containing accelerated stability data (required).
#' @param y Name of decreasing variable (e.g. concentration) contained within data
#'  (required).
#' @param .time Time variable contained within data (required).
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
#' @examples #load antigenicity and potency data.
#' data(antigenicity)
#' data(potency)
#'
#' #Basic use of the step1_down function with C column defined.
#' fit1 <- step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", draw = 5000)
#'
#' #Basic use of the step1_down function with K column defined & Validation data segmented out.
#' fit2 <- step1_down(data = antigenicity, y = "conc", .time = "time", K = "K",
#' validation = "validA", draw = 5000)
#'
#' #When zero_order = FALSE, the output suggests using zero_order = TRUE for Potency dataset.
#' fit3 <- step1_down(data = potency, y = "Potency", .time = "Time",C = "Celsius",
#'   reparameterisation = FALSE, zero_order = TRUE, draw = 5000)
#'
#' #reparameterisation is TRUE.
#' fit4 <- step1_down(data = antigenicity, y = "conc", .time = "time",C = "Celsius",
#'   reparameterisation = TRUE, draw = 5000)
#'
#' #Use a custom lower bound for k1 (default is 0).
#' fit5 <- step1_down(data = potency, y = "Potency", .time = "Time",C = "Celsius",
#'   reparameterisation = TRUE, zero_order = TRUE, draw = 5000, lower = c(-Inf, 0, 0))
#'
#' @importFrom stats vcov coef runif confint rnorm quantile qt complete.cases
#' @importFrom minpack.lm nls.lm
#' @importFrom mvtnorm rmvt
#'
#' @export step1_down

step1_down <- function (data, y, .time, K = NULL, C = NULL, validation = NULL,
                        draw = 10000, parms = NULL, temp_pred_C = NULL,
                        max_time_pred = NULL, confidence_interval = 0.95, by = 101,
                        reparameterisation = FALSE, zero_order = FALSE, ...){

  if (is.null(K) & is.null(C))
    stop("Select the temperature variable in Kelvin or Celsius")

  if (!is.null(parms) & !is.list(parms))
    stop("The starting values for parameters must be a list, or keep as NULL")

  if (!is.null(validation))
    if (!all(data[,validation] %in% c(0,1)))
      stop("Validation column must contain 1s and 0s only")

  user_parameters <- list(
    data = data, y = y, .time = .time, K = K, C = C, validation = validation,draw = draw,
    parms = parms, temp_pred_C = temp_pred_C, max_time_pred = max_time_pred,
    confidence_interval = confidence_interval, by = by,
    reparameterisation = reparameterisation, zero_order = zero_order)

  ## Additional arguments in the call will be passed to model fitting with minipak.lm
  minpack_args = list(...)                    ##

  ## Temperature: both C and K are provided
  if(!is.null(C) & !is.null(K)) {
    data[, C] <- ifelse(is.na(data[, C]) & !is.na(data[, K]),
                        data[, K] - 273.15,
                        data[, C])
    data[, K] <- ifelse(is.na(data[, K]) & !is.na(data[, C]),
                        data[, C] + 273.15,
                        data[, K])   }

  ## Temperature: only C or only K is provided
  if (!is.null(C) & is.null(K)) {             ##
   K = 'K'                                  ##
   data[, K] = data[, C] + 273.15  }        ##
  if (!is.null(K) & is.null(C)) {             ##
    C = 'C'                                 ##
    data[, C] = data[, K] - 273.15 }        ##

  data <- data[complete.cases(data[, c(C,K,y,.time)]), ]

  dat = data
  dat$K = dat[, K]                         ##
  Kref = mean(dat$K)
  dat$Celsius = as.factor(dat[, C])
  dat$time = dat[, .time]
  dat$y = dat[, y]
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

  Temps = sort(unique(dat$K))
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
      Model = c0 - c0 * dat$time * exp(k1 - k2/dat$K +
                                         k2/Kref)
      residual = dat$y - Model
      return(residual)
    }

  if (!"fn" %in% names(minpack_args)) {	##
    minpack_args$fn =  MyFctNL    }     ##

    # Fit model :
    if (!is.null(parms)) {
     minpack_args$par =  parms                          ##
    if (!"lower" %in% names(minpack_args)) 	{	##
    minpack_args$lower =  rep(0, length(parms))   }     ##

	if(length(minpack_args$par) != length(minpack_args$lower))                             ##
	stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")  ##

  fit = do.call(minpack.lm::nls.lm, minpack_args)
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 40), k2 = stats::runif(1, 1000, 20000), c0 = c0_initial)

minpack_args$par = parms

if (!"lower" %in% names(minpack_args)) 	{	##
	minpack_args$lower =  rep(0, length( parms ))   }     ##

	if(length(minpack_args$par) != length(minpack_args$lower))                             ##
	stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")  ##

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

    # Calculate the predictions
    k1 = stats::coef(fit)[1]
    k2 = stats::coef(fit)[2]
    c0 = stats::coef(fit)[3]
    SIG = stats::vcov(fit)
    sigma = summary(fit)$sigma
    DF = summary(fit)$df[2]
    pred = expand.grid(time = times.pred, K = Temps)
    pred$Degradation = pred$time * exp(k1 - k2/pred$K + k2/Kref)
    pred$Response = c0 - c0 * pred$Degradation

    if(is.null(draw)){
    pred$derivk1 = -c0 * pred$Degradation
    pred$derivk2 = -c0 * (1/Kref - 1/pred$K) * pred$Degradation
    pred$derivc0 = 1 - pred$Degradation
    pred$varY = (pred$derivk1)^2 * SIG[1, 1] + (pred$derivk2)^2 *
      SIG[2, 2] + (pred$derivc0)^2 * SIG[3, 3] + 2 * pred$derivk1 *
      pred$derivk2 * SIG[1, 2] + 2 * pred$derivk1 * pred$derivc0 *
      SIG[1, 3] + 2 * pred$derivk2 * pred$derivc0 * SIG[2,
                                                        3]
    pred$derivk1 = pred$derivk2 = pred$derivc0 = NULL}else{ # Bootstrap
      pred_fct = function(coef.fit)
      {
        degrad = pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K + coef.fit[2] / Kref)
        conc = coef.fit[3] - coef.fit[3]*degrad
        return(conc)
      }
      # Multi T bootstrap
      rand.coef = mvtnorm::rmvt(draw, sigma = SIG, df = nrow(dat) - 3) + matrix(nrow = draw, ncol = 3, byrow = TRUE, coef(fit))
      res.boot = matrix(nrow = draw, ncol = nrow(pred), byrow = TRUE, apply(rand.coef, 1, pred_fct))

      CI1b = apply(res.boot, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      CI2b = apply(res.boot, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

      res.boot = res.boot + rnorm(draw*length(pred$time), 0, sigma)
      PI1b = apply(res.boot, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      PI2b = apply(res.boot, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)
    }

## Model type 2 - no reparameterisation and k3 = 0
  }else if(!reparameterisation & zero_order){
    MyFctNL = function(parms) { # make function
      k1 = parms$k1
      k2 = parms$k2
      c0 = parms$c0

      Model = c0 - c0 * dat$time * exp(k1 - k2 / dat$K)


      residual = dat$y - Model
      return(residual)
    }

  if (!"fn" %in% names(minpack_args)) 	{	##
    minpack_args$fn =  MyFctNL    }             ##

    if (!is.null(parms)) { # fit model
    minpack_args$par =  parms                       ##
    if (!"lower" %in% names(minpack_args)) 	{   ##
    minpack_args$lower =  rep(0, length(parms))   } ##

	if(length(minpack_args$par) != length(minpack_args$lower))                             ##
	stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")  ##

      fit = do.call(minpack.lm::nls.lm, minpack_args)
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 40), k2 = stats::runif(1, 1000, 20000), c0 = c0_initial)

  minpack_args$par = parms ##

if (!"lower" %in% names(minpack_args)) 	{	##
	    minpack_args$lower =  rep(0, length( parms ))   }   ##

	if(length(minpack_args$par) != length(minpack_args$lower))                             ##
  stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")  ##

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

# Predict
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    c0 = coef(fit)[3]
    SIG = vcov(fit)
    sigma = summary(fit)$sigma
    DF = summary(fit)$df[2]

    pred = expand.grid("time" = times.pred, K = Temps)
    pred$Degradation = pred$time * exp(k1 - k2 / pred$K)
    pred$Response = c0 - c0*pred$Degradation

    if(is.null(draw)){
    pred$derivk1 = -c0 * pred$Degradation
    pred$derivk2 = c0 / pred$K * pred$Degradation
    pred$derivc0 = 1 - pred$Degradation
    pred$varY = (pred$derivk1)^2 * SIG[1,1] + (pred$derivk2)^2 * SIG[2,2] + (pred$derivc0)^2 * SIG[3,3] +
      2*pred$derivk1*pred$derivk2 * SIG[1,2] + 2*pred$derivk1*pred$derivc0 * SIG[1,3] + 2*pred$derivk2*pred$derivc0 * SIG[2,3]

    pred$derivk1 = pred$derivk2 = pred$derivc0 = NULL}else{ # Bootstrap
      pred_fct = function(coef.fit)
      {

        degrad = pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K)
        conc = coef.fit[3] - coef.fit[3]*degrad
        return(conc)
      }
      # Multi T bootstrap
      rand.coef = mvtnorm::rmvt(draw, sigma = SIG, df = nrow(dat) - 3) + matrix(nrow = draw, ncol = 3, byrow = TRUE, coef(fit))
      res.boot = matrix(nrow = draw, ncol = nrow(pred), byrow = TRUE, apply(rand.coef, 1, pred_fct))

      CI1b = apply(res.boot, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      CI2b = apply(res.boot, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

      res.boot = res.boot + rnorm(draw*length(pred$time), 0, sigma)
      PI1b = apply(res.boot, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      PI2b = apply(res.boot, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)


    }

## Model type 3 - reparameterisation and k3 is not zero
  }else if(reparameterisation & !zero_order){

  ## Print a message informing lower bounds = 0 may not be suitable with the reparamerised version
  cat("The alternative parameterisation of the one-step model was used. Note that the lower bounds for all parameters are set to 0 unless other lower bounds are specified in step1_down() or step1_down_basic().\n\n")

   MyFctNL = function(parms) {
      k1 = parms$k1
      k2 = parms$k2
      k3 = parms$k3
      c0 = parms$c0
      Model = c0 - c0 * (1 - ((1 - k3) * (1/(1 - k3) - dat$time *
                                            exp(k1 - k2/dat$K + k2/Kref)))^(1/(1 - k3)))
      residual = dat$y - Model
      return(residual)
    }

  if (!"fn" %in% names(minpack_args)) 	{	##
    minpack_args$fn =  MyFctNL    }             ##

  if (!is.null(parms)) { # Fit the model
        minpack_args$par =  parms               ##
  if (!"lower" %in% names(minpack_args)) 	{	##
    minpack_args$lower =  rep(0, length(parms))   }     ##

	if(length(minpack_args$par) != length(minpack_args$lower))                             ##
  stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")  ##

      fit = do.call(minpack.lm::nls.lm, minpack_args)
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 60), k2 = stats::runif(1, 1000, 20000), k3 = stats::runif(1, 0, 11), c0 = c0_initial)

  minpack_args$par = parms

  if (!"lower" %in% names(minpack_args)) 	{	##
	    minpack_args$lower =  rep(0, length( parms ))   }     ##

	if(length(minpack_args$par) != length(minpack_args$lower))                             ##
  stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")  ##

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

    # Predict
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    k3 = coef(fit)[3]
    if (k3 == 0){cat(paste("k3 is fitted to be exactly 0, we strongly suggest using option zero_order = TRUE","The model will continue with k3 = 0, so degradation is linear over time"," "," ", sep = "\n"))
    }else if(confint(fit,'k3')[1] < 0 && confint(fit,'k3')[2] > 0){print(paste0("The 95% Wald Confidence Interval for k3 includes 0, k3 is estimated as ",signif(k3,4),". We suggest considering option zero_order = TRUE"))}
    c0 = coef(fit)[4]
    SIG = vcov(fit)
    sigma = summary(fit)$sigma
    DF = summary(fit)$df[2]

    pred = expand.grid("time" = times.pred, K = Temps)
    pred$Degradation = 1 - ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2 / pred$K + k2 / Kref)))^(1/(1-k3))
    pred$Response = c0 - c0*pred$Degradation

    if(is.null(draw)){
    pred$derivk1 = c0 * pred$time * (-exp(k1 - k2/pred$K + k2/Kref)) * ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K + k2/Kref)))^(1/(1 - k3) - 1)
    pred$derivk2 = c0 * pred$time * (1/Kref - 1/pred$K) * (-exp(k1 - k2/pred$K + k2/Kref)) * ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K + k2/Kref)))^(1/(1 - k3) - 1)
    pred$derivk3 = c0 * ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K + k2/Kref)))^(1/(1 - k3)) * ((pred$time * exp(k1 - k2/pred$K + k2/Kref)) / ((1 - k3)^2 * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K + k2/Kref))) + log((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K + k2/Kref)))/(1 - k3)^2)
    pred$derivc0 = 1 - pred$Degradation

    pred$varY = (pred$derivk1)^2 * SIG[1,1] + (pred$derivk2)^2 * SIG[2,2] + (pred$derivk3)^2 * SIG[3,3] + (pred$derivc0)^2 * SIG[4,4] +
      2*pred$derivk1*pred$derivk2 * SIG[1,2] + 2*pred$derivk1*pred$derivk3 * SIG[1,3] + 2*pred$derivk1*pred$derivc0 * SIG[1,4] +
      2*pred$derivk2*pred$derivk3 * SIG[2,3] + 2*pred$derivk2*pred$derivc0 * SIG[2,4] + 2*pred$derivk3*pred$derivc0 * SIG[3,4]
    pred$derivk1 = pred$derivk2 = pred$derivk3 = pred$derivc0 = NULL}else{ # Bootstrap
      pred_fct = function(coef.fit)
      {
        degrad = 1 - ((1 - coef.fit[3]) * (1/(1 - coef.fit[3]) - pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K + coef.fit[2] / Kref)))^(1/(1-coef.fit[3]))
        conc = coef.fit[4] - coef.fit[4]*degrad
        return(conc)
      }
      # Multi T bootstrap
      rand.coef = mvtnorm::rmvt(draw, sigma = SIG, df = nrow(dat) - 4) + matrix(nrow = draw, ncol = 4, byrow = TRUE, coef(fit))
      res.boot = matrix(nrow = draw, ncol = nrow(pred), byrow = TRUE, apply(rand.coef, 1, pred_fct))

      no_k3_below0 <- sum(rand.coef[,3] < 0)
      if(no_k3_below0 > 0.5){
        cat(paste(paste0(no_k3_below0*100/draw, "% of the bootstraps for k3 are below zero, this might have an adverse effect on the confidence interval, particularly if this value exceeds the confidence %."),"We suggest considering option zero_order = TRUE", sep = "\n"))
      }

      CI1b = apply(res.boot, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      CI2b = apply(res.boot, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

      res.boot = res.boot + rnorm(draw*length(pred$time), 0, sigma)
      PI1b = apply(res.boot, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      PI2b = apply(res.boot, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)


    }

## Model type 4 - no reparameterisation and k3 is not 0
  }else if(!reparameterisation & !zero_order){
    MyFctNL = function(parms) {
      k1 = parms$k1
      k2 = parms$k2
      k3 = parms$k3
      c0 = parms$c0

      test = c0 - c0 * (1 - ((1 - k3) * (1/(1 - k3) - dat$time * exp(k1 - k2 / dat$K)))^(1/(1-k3)))

      residual = dat$y - test
      return(residual)

    }

  if (!"fn" %in% names(minpack_args)) 	{	##
    minpack_args$fn =  MyFctNL    }             ##

    if (!is.null(parms)) { # Fitting the model
    minpack_args$par =  parms                   ##
  if (!"lower" %in% names(minpack_args)) 	{	##
    minpack_args$lower =  rep(0, length(parms))   }   ##

	if(length(minpack_args$par) != length(minpack_args$lower))                             ##
  stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")  ##

      fit = do.call(minpack.lm::nls.lm, minpack_args)
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 60), k2 = stats::runif(1, 1000, 20000), k3 = stats::runif(1, 0, 11), c0 = c0_initial)

  minpack_args$par = parms

  if (!"lower" %in% names(minpack_args)) 	{	##
	    minpack_args$lower =  rep(0, length( parms ))   } ##

	if(length(minpack_args$par) != length(minpack_args$lower))                             ##
  stop("The number of parameters (",length(minpack_args$par),") does not match the number of specified lower bounds (",length(minpack_args$lower),").")  ##

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

    # Predict
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    k3 = coef(fit)[3]
    if (k3 == 0){cat(paste("k3 is fitted to be exactly 0, we strongly suggest using option zero_order = TRUE","The model will continue with k3 = 0, so degradation is linear over time"," ", " ", sep = "\n"))
    }else if(confint(fit,'k3')[1] < 0 && confint(fit,'k3')[2] > 0){print(paste0("The 95% Wald Confidence Interval for k3 includes 0, k3 is estimated as ",signif(k3,4),". We suggest considering option zero_order = TRUE"))}
    c0 = coef(fit)[4]
    SIG = vcov(fit)
    sigma = summary(fit)$sigma
    DF = summary(fit)$df[2]

    pred = expand.grid("time" = times.pred, K = Temps)
    pred$Degradation = 1 - ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2 / pred$K)))^(1/(1-k3))
    pred$Response = c0 - c0*pred$Degradation

    if(is.null(draw)){ # Derivatives
    pred$derivk1 = c0 * pred$time * (-exp(k1 - k2/pred$K)) * ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K)))^(1/(1 - k3) - 1)
    pred$derivk2 = (c0 * pred$time * exp(k1 - k2/pred$K) * ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K)))^(1/(1 - k3) - 1)) / pred$K
    pred$derivk3 = c0 * ((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K)))^(1/(1 - k3)) * ((pred$time * exp(k1 - k2/pred$K)) / ((1 - k3)^2 * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K))) + log((1 - k3) * (1/(1 - k3) - pred$time * exp(k1 - k2/pred$K)))/(1 - k3)^2)
    pred$derivc0 = 1 - pred$Degradation
    pred$varY = (pred$derivk1)^2 * SIG[1,1] + (pred$derivk2)^2 * SIG[2,2] + (pred$derivk3)^2 * SIG[3,3] + (pred$derivc0)^2 * SIG[4,4] +
      2*pred$derivk1*pred$derivk2 * SIG[1,2] + 2*pred$derivk1*pred$derivk3 * SIG[1,3] + 2*pred$derivk1*pred$derivc0 * SIG[1,4] +
      2*pred$derivk2*pred$derivk3 * SIG[2,3] + 2*pred$derivk2*pred$derivc0 * SIG[2,4] + 2*pred$derivk3*pred$derivc0 * SIG[3,4]
    pred$derivk1 = pred$derivk2 = pred$derivk3 = pred$derivc0 = NULL }else{ # Bootstrap
      pred_fct = function(coef.fit)
      {
        degrad = 1 - ((1 - coef.fit[3]) * (1/(1 - coef.fit[3]) - pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K)))^(1/(1-coef.fit[3]))
        conc = coef.fit[4] - coef.fit[4]*degrad
        return(conc)
      }
      # Multi T bootstrap
      rand.coef = mvtnorm::rmvt(draw, sigma = SIG, df = nrow(dat) - 4) + matrix(nrow = draw, ncol = 4, byrow = TRUE, coef(fit))
      res.boot = matrix(nrow = draw, ncol = nrow(pred), byrow = TRUE, apply(rand.coef, 1, pred_fct))

      no_k3_below0 <- sum(rand.coef[,3] < 0)
      if(no_k3_below0 > 0.5){
        cat(paste(paste0(no_k3_below0*100/draw, "% of the bootstraps for k3 are below zero, this might have an adverse effect on the confidence interval, particularly if this value exceeds the confidence %."),"We suggest considering option zero_order = TRUE", sep = "\n"))
      }

      CI1b = apply(res.boot, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      CI2b = apply(res.boot, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

      res.boot = res.boot + rnorm(draw*length(pred$time), 0, sigma)
      PI1b = apply(res.boot, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
      PI2b = apply(res.boot, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)
    }


  }

  pred$Celsius = as.factor(pred$K - 273.15)
  pred$K = as.factor(pred$K)
  pred$fit = "Prediction"
  pred$CI = paste(100*confidence_interval, "% CI")
  pred$PI = paste(100*confidence_interval, "% PI")

  if(is.null(draw)){
  pred$CI1 = pred$Response - qt(0.5 + confidence_interval/2, summary(fit)$df[2]) * sqrt(pred$varY)
  pred$CI2 = pred$Response + qt(0.5 + confidence_interval/2, summary(fit)$df[2]) * sqrt(pred$varY)
  pred$PI1 = pred$Response - qt(0.5 + confidence_interval/2, summary(fit)$df[2]) * sqrt(pred$varY + sigma^2)
  pred$PI2 = pred$Response + qt(0.5 + confidence_interval/2, summary(fit)$df[2]) * sqrt(pred$varY + sigma^2)
  }else{
    pred$CI1 = CI1b
    pred$CI2 = CI2b
    pred$PI1 = PI1b
    pred$PI2 = PI2b}

  if(is.null(draw)){
    rand.coef = "Calculus method used"
  }else{
    if(zero_order){
      colnames(rand.coef) <- c("k1","k2","c0")
    }else{
      colnames(rand.coef) <- c("k1","k2","k3","c0")
    }
  }

  results = list(fit, dat_full, pred,user_parameters, rand.coef)
  names(results) = c("fit", "data", "prediction","user_parameters","sample_coefficients")
  class(results) = "SB"
  return(results)

  }
