#' @title Sample the Multivariate t Distribution
#'
#' @description Take a selected number of samples from the multivariate t distribution (mvt).
#'
#' @details Using the provided data the function creates a fit of the
#'  Šesták–Berggren kinetic model and then draws a selected number of
#'   samples from the mvt of the model parameters.
#'
#' @param data Dataframe containing accelerated stability data (required).
#' @param y Name of decreasing variable (e.g. concentration) contained within data
#'  (required).
#' @param .time Time variable contained within data (required).
#' @param K Kelvin variable (numeric or column name) (optional).
#' @param C Celsius variable (numeric or column name) (optional).
#' @param validation Validation dummy variable (column name) (optional).
#' @param draw Number of samples to draw from mvt (required).
#' @param parms Starting values for the parameters as a list - k1, k2, k3, and c0 (optional).
#' @param reparameterisation Use alternative parameterisation of the one-step
#'  model which aims to reduce correlation between k1 and k2.
#' @param zero_order Set kinetic order, k3, to zero (straight lines).
#' @param ... Further arguments to passed to minpack.lm.
#'
#' @return A matrix containing parameter draws from the mvt distribution.
#'
#' @examples #load antigenicity data.
#' data(antigenicity)
#'
#' #Basic use of the step1_sample_mvt function with C column defined and 1000 draws.
#' sample1 <- step1_sample_mvt(data = antigenicity, y = "conc", .time = "time",
#'  C = "Celsius", draw = 1000)
#'
#' #Basic use of the step1_sample_mvt function with K column defined and 50000 draws
#' sample2 <- step1_sample_mvt(data = antigenicity, y = "conc", .time = "time",
#'  K = "K", draw = 50000)
#'
#' #reparameterisation is TRUE and 10000 draws.
#' sample3 <- step1_sample_mvt(data = antigenicity, y = "conc", .time = "time",
#' C = "Celsius", reparameterisation = TRUE, draw = 10000)
#'
#' @importFrom stats vcov coef runif confint rnorm quantile qt complete.cases
#' @importFrom minpack.lm nls.lm
#'
#' @export step1_sample_mvt

step1_sample_mvt <- function (data, y, .time, K = NULL, C = NULL, validation = NULL,
                        draw, parms = NULL, reparameterisation = FALSE, zero_order = FALSE, ...){

  if (is.null(K) & is.null(C))
    stop("Select the temperature variable in Kelvin or Celsius")

  if (!is.null(parms) & !is.list(parms))
    stop("The starting values for parameters must be a list, or keep as NULL")

    if (!is.null(validation))
    if (!all(data[,validation] %in% c(0,1)))
      stop("Validation column must contain 1s and 0s only")

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

  data <- data[complete.cases(data[, c(C,K,y,.time)]), ]

  dat = data
  dat$K = dat[, K]
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
  if(reparameterisation & zero_order){ # reparameterisation and k3 is 0

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

    SIG = stats::vcov(fit)
    DF = summary(fit)$df[2]
    n.params = summary(fit)$df[1]

        pred_fct = function(coef.fit)
        {
          degrad = pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K + coef.fit[2] / Kref)
          conc = coef.fit[3] - coef.fit[3]*degrad
          return(conc)
        }
        # Multi T samples
        rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
        rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)

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

    SIG = vcov(fit)
    DF = summary(fit)$df[2]
    n.params = summary(fit)$df[1]

        pred_fct = function(coef.fit)
        {

          degrad = pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K)
          conc = coef.fit[3] - coef.fit[3]*degrad
          return(conc)
        }
        # Multi T samples
        rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
        rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)

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

    k3 = coef(fit)[3]
     if (k3 == 0){cat(paste("k3 is fitted to be exactly 0, we strongly suggest using option zero_order = TRUE","The model will continue with k3 = 0, so degradation is linear over time"," "," ", sep = "\n"))
    }else if(confint(fit,'k3')[1] < 0 && confint(fit,'k3')[2] > 0){print(paste0("The 95% Wald Confidence Interval for k3 includes 0, k3 is estimated as ",signif(k3,4),". We suggest considering option zero_order = TRUE"))}
    SIG = vcov(fit)
    DF = summary(fit)$df[2]
    n.params = summary(fit)$df[1]

        pred_fct = function(coef.fit)
        {
          degrad = 1 - ((1 - coef.fit[3]) * (1/(1 - coef.fit[3]) - pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K + coef.fit[2] / Kref)))^(1/(1-coef.fit[3]))
          conc = coef.fit[4] - coef.fit[4]*degrad
          return(conc)
        }
        # Multi T samples
        rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
        rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)

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

    k3 = coef(fit)[3]
    if (k3 == 0){cat(paste("k3 is fitted to be exactly 0, we strongly suggest using option zero_order = TRUE","The model will continue with k3 = 0, so degradation is linear over time"," ", " ", sep = "\n"))
    }else if(confint(fit,'k3')[1] < 0 && confint(fit,'k3')[2] > 0){print(paste0("The 95% Wald Confidence Interval for k3 includes 0, k3 is estimated as ",signif(k3,4),". We suggest considering option zero_order = TRUE"))}
    SIG = vcov(fit)
    DF = summary(fit)$df[2]
    n.params = summary(fit)$df[1]

        pred_fct = function(coef.fit)
        {
          degrad = 1 - ((1 - coef.fit[3]) * (1/(1 - coef.fit[3]) - pred$time * exp(coef.fit[1] - coef.fit[2] / pred$K)))^(1/(1-coef.fit[3]))
          conc = coef.fit[4] - coef.fit[4]*degrad
          return(conc)
        }
        # Multi T samples
        rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
        rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)

  }

  if(zero_order){
    colnames(rand.coef) <- c("k1","k2","c0")
  }else{
    colnames(rand.coef) <- c("k1","k2","k3","c0")
  }

  return(rand.coef)

}
utils::globalVariables(c("pred"))
