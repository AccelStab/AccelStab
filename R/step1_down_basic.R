#' @title Basic version Step1 Down Model
#'
#' @description Quickly fit the one-step Šesták–Berggren kinetic model.
#'
#' @details Fit the one-step Šesták–Berggren kinetic (non-linear) model using
#'  accelerated stability data from an R dataframe format. Fit summary is printed
#'  to the console and the fit itself is returned. Fewer arguments and
#'  outputs than the main function for rapid testing. Parameters are kept in
#'  even when not significant.
#'
#' @param data Dataframe containing accelerated stability data (required).
#' @param y Name of decreasing variable (e.g. concentration) contained within data
#'  (required).
#' @param .time Time variable contained within data (required).
#' @param K Kelvin variable (numeric or column name) (optional).
#' @param C Celsius variable (numeric or column name) (optional).
#' @param validation Validation dummy variable, the column must contain only 1s and 0s, 1 for validation data and 0 for fit data. (column name) (optional).
#' @param parms Starting values for the parameters as a list - k1, k2, k3, and c0.
#' @param reparameterisation Use alternative parameterisation of the one-step
#'  model which aims to reduce correlation between k1 and k2.
#' @param zero_order Set kinetic order, k3, to zero (straight lines).
#'
#' @return The fit object
#'
#' @examples #load antigenicity and potency data.
#' data(antigenicity)
#' data(potency)
#'
#' #Use of the step1_down_basic function with C column defined.
#' fit1 <- step1_down_basic(data = antigenicity, y = "conc", .time = "time", C = "Celsius")
#'
#' #Basic use of the step1_down_basic function with K column defined & Validation data segmented out.
#' fit2 <- step1_down_basic(data = antigenicity, y = "conc", .time = "time", K = "K",
#' validation = "validA")
#'
#' #When zero_order = FALSE, the output suggests using zero_order = TRUE for Potency dataset.
#' fit3 <- step1_down_basic(data = potency, y = "Potency", .time = "Time",C = "Celsius",
#'   reparameterisation = FALSE, zero_order = TRUE)
#'
#' #reparameterisation is TRUE.
#' fit4 <- step1_down_basic(data = antigenicity, y = "conc", .time = "time",C = "Celsius",
#'   reparameterisation = TRUE)
#'
#' @importFrom stats coef vcov runif complete.cases
#' @importFrom minpack.lm nls.lm
#'
#' @export step1_down_basic

step1_down_basic <- function (data, y, .time, K = NULL, C = NULL, validation = NULL,
                        parms = NULL, reparameterisation = FALSE, zero_order = FALSE){

  if (is.null(K) & is.null(C))
    stop("Select the temperature variable in Kelvin or Celsius")
  if (!is.null(parms) & !is.list(parms))
    stop("The starting values for parameters must be a list, or keep as NULL")

  user_parameters <- list(
    data = data, y = y, .time = .time, K = K, C = C, validation = validation,
    parms = parms, reparameterisation = reparameterisation, zero_order = zero_order
  )

  if(!is.null(C) & !is.null(K)) {

    data[, C] <- ifelse(is.na(data[, C]) & !is.na(data[, K]),
                        data$K - 273.15,
                        data[, C])

    data[, K] <- ifelse(is.na(data[, K]) & !is.na(data[, C]),
                        data$C + 273.15,
                        data[, K])
  }

  data <- data[complete.cases(data[, c(C,K,y,.time)]), ]

  dat = data

  if (!is.null(validation))
    if (!all(dat[,validation] %in% c(0,1)))
      stop("Validation column must contain 1s and 0s only")

  if (is.null(K))
    dat$K = dat[, C] + 273.15
  if (is.null(C)) {
    dat$C = dat[, K] - 273.15
    C = "C"}

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

  if(reparameterisation & zero_order){ # reparameterisation and k3 is 0
    MyFctNL = function(parms) { # Make function
      k1 = parms$k1
      k2 = parms$k2
      c0 = parms$c0
      Model = c0 - c0 * dat$time * exp(k1 - k2/dat$K +
                                         k2/Kref)
      residual = dat$y - Model
      return(residual)
    }

    # Fit model :
    if (!is.null(parms)) {
      fit = minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0,
                                                                      length(parms)))
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 40), k2 = stats::runif(1,
                                                                    1000, 20000), c0 = c0_initial)
        fit = suppressWarnings(minpack.lm::nls.lm(par = parms,
                                                  fn = MyFctNL, lower = rep(0, length(parms))))
        fit <- tryCatch({
          suppressWarnings(minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0, length(parms))))
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
      fit = minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0,
                                                                      length(parms)))
    }

  }else if(!reparameterisation & zero_order){ # no reparameterisation and k3 is 0
    MyFctNL = function(parms) { # make function
      k1 = parms$k1
      k2 = parms$k2
      c0 = parms$c0

      Model = c0 - c0 * dat$time * exp(k1 - k2 / dat$K)


      residual = dat$y - Model
      return(residual)
    }
    if (!is.null(parms)) { # fit model
      fit = minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0,
                                                                      length(parms)))
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 40), k2 = stats::runif(1,
                                                                    1000, 20000), c0 = c0_initial)
        fit <- tryCatch({
          suppressWarnings(minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0, length(parms))))
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
      fit = minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0,
                                                                      length(parms)))
    }

  }else if(reparameterisation & !zero_order){ #reparameterisation and k3 is not zero
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
    if (!is.null(parms)) { # Fit the model
      fit = minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0,
                                                                      length(parms)))
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 60), k2 = stats::runif(1,
                                                                    1000, 20000), k3 = stats::runif(1, 0, 11), c0 = c0_initial)
        fit <- tryCatch({
          suppressWarnings(minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0, length(parms))))
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
      fit = minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0,
                                                                      length(parms)))
    }

   if (coef(fit)[3] == 0){cat(paste("k3 is fitted to be exactly 0, we strongly suggest using option zero_order = TRUE","The model will continue with k3 = 0, so degradation is linear over time"," "," ", sep = "\n"))
    }else if(confint(fit,'k3')[1] < 0 && confint(fit,'k3')[2] > 0){print(paste0("The 95% Wald Confidence Interval for k3 includes 0, k3 is estimated as ",signif(coef(fit)[3],4),". We suggest considering option zero_order = TRUE"))}

  }else if(!reparameterisation & !zero_order){ # No re-parameterisation and k3 not zero
    MyFctNL = function(parms) {
      k1 = parms$k1
      k2 = parms$k2
      k3 = parms$k3
      c0 = parms$c0

      test = c0 - c0 * (1 - ((1 - k3) * (1/(1 - k3) - dat$time * exp(k1 - k2 / dat$K)))^(1/(1-k3)))

      residual = dat$y - test
      return(residual)

    }
    if (!is.null(parms)) { # Fitting the model
      fit = minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0,
                                                                      length(parms)))
    }
    else {
      repeat {
        suppressWarnings(rm(fit))

        parms = list(k1 = stats::runif(1, 0, 60), k2 = stats::runif(1,
                                                                    1000, 20000), k3 = stats::runif(1, 0, 11), c0 = c0_initial)
        fit <- tryCatch({
          suppressWarnings(minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0, length(parms))))
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
      fit = minpack.lm::nls.lm(par = parms, fn = MyFctNL, lower = rep(0,
                                                                      length(parms)))
    }
    if (coef(fit)[3] == 0){cat(paste("k3 is fitted to be exactly 0, we strongly suggest using option zero_order = TRUE","The model will continue with k3 = 0, so degradation is linear over time"," ", " ", sep = "\n"))
    }else if(confint(fit,'k3')[1] < 0 && confint(fit,'k3')[2] > 0){print(paste0("The 95% Wald Confidence Interval for k3 includes 0, k3 is estimated as ",signif(coef(fit)[3],4),". We suggest considering option zero_order = TRUE"))}
  }

  print(summary(fit))

  return(fit)

}




