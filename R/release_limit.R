#' @title Calculate and Plot Release Limit Information
#'
#' @description Calculate the release limit for a product at a given shelf temperature, shelf time
#' and lower specification limit.
#'
#' @details Release limit information is presented both in a table and in a plot displaying the lower specification
#' limit (LSL) along with the predicted degradation and statistical intervals. The method for computing the intervals
#' is inherited from the given step1_down object (using either draws from the multivariate t-distribution or
#' analytical formulae). In addition to the release limit predictions, the plot displays the original fit for the
#' given shelf temperature as well as any data points at that temperature within the given shelf time.
#' Fits using step1_down_batch() are partially supported; at present, the release limit is calculated
#' only for the reference level of the batch variable (default is the first level).
#'
#' @param step1_down_object The fit object from the step1_down() or step1_down_batch() function (required).
#' @param shelf_temp The temperature (in Celsius) at which the product will be stored (required).
#' @param shelf_time The length of time for which the product will be stored (required).
#' @param LSL A numeric value for the lower specification limit (required).
#' @param interval Choose between confidence intervals or prediction intervals to be displayed in the plot (default is CI).
#' @param confidence_interval Confidence level for the confidence and prediction intervals
#'  around the predictions (default 0.95).
#' @param xname Label for the x-axis (optional).
#' @param yname Label for the y-axis (optional).
#' @param xlim x-axis limits (optional).
#' @param ylim y-axis limits (optional).
#'
#' @return Returns a list containing (1) the release limit calculations, (2) a plot of predicted degradation and
#' (3) predicted values using the release limit.
#'
#' @examples
#'
#' # Run step1_down() fit
#' fit1 <- step1_down(data = antigenicity, y = "conc", .time = "N.days",
#'                     C = "Celsius", validation = "validA")
#' 
#' myRL <- release_limit(step1_down_object = fit1, shelf_temp = 5, shelf_time = 420,
#'                      LSL = 75, interval = "PI", confidence_interval = 0.95,
#'                      yname = "Concentration", xname = "Days")
#'
#' @import ggplot2
#'
#' @export release_limit

# Function options
release_limit <- function(step1_down_object = NULL, shelf_temp = NULL, shelf_time = NULL, LSL = NULL,
                          interval = "CI", confidence_interval = 0.95,
                          xname = NULL, yname = NULL, xlim = NULL, ylim = NULL){

# Examine whether the inputted arguments are appropriate before continuing with the function
if(is.null(step1_down_object)){
stop("A step1_down_object is required for calculating the release limit.") }

if(is.null(shelf_temp) | !is.numeric(shelf_temp)){
stop("Enter shelf_temp (in Celsius) as a numeric variable.") }

if(is.null(shelf_time) | !is.numeric(shelf_time)){
stop("Enter shelf_time as a numeric variable.") }

if(is.null(LSL) | !is.numeric(LSL)){
stop("Enter LSL (lower specification limit) as a numeric variable.") }
  
if(!(interval %in% c("CI", "PI"))){
stop("Ensure interval is entered as either 'CI' or 'PI'.") }

if(!is.null(step1_down_object$user_parameters$batch)){
cat("Models with batch effects are not yet fully supported. Release limit information will be returned only for the reference level of the batch variable (default is the first level).\n\n") }

# Preparations for the release limit (RL) calculation and predictions
  fit <- step1_down_object$fit
 
  Kref = mean(step1_down_object$data$K)

  preds <- data.frame( # Making empty prediction frame
    temps = rep(shelf_temp,301),
    K = rep(shelf_temp + 273.15,301),
    conc = rep(NA,301),
    degrad = rep(NA,301),
    total_time = seq(0,shelf_time,length.out = 301))

  SIG = vcov(fit)
  sigma = summary(fit)$sigma
  DF = summary(fit)$df[2]
  n.params = summary(fit)$df[1]

## STATISTICAL INTERVALS COMPUTED USING ANALYTICAL METHOD
if(is.null(step1_down_object$user_parameters$draw)){

# Model Type 4 ("standard" model) predictions
if(step1_down_object$user_parameters$reparameterisation == F && step1_down_object$user_parameters$zero_order == F){

  k1 = coef(fit)[1]
  k2 = coef(fit)[2]
  k3 = coef(fit)[3]
  c0 = coef(fit)[4]
  preds$degrad <- (1 - ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2 / (preds$K))))^(1/(1-k3)))
  preds$conc <- c0 - c0 * preds$degrad

	# Model Type 4 derivatives for predictions
  preds$derivk1 = c0 * preds$total_time * (-exp(k1 - k2/preds$K)) * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K)))^(1/(1 - k3) - 1)
  preds$derivk2 = (c0 * preds$total_time * exp(k1 - k2/preds$K) * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K)))^(1/(1 - k3) - 1)) / preds$K
  preds$derivk3 = c0 * ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K)))^(1/(1 - k3)) * ((preds$total_time * exp(k1 - k2/preds$K)) / ((1 - k3)^2 * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K))) + log((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2/preds$K)))/(1 - k3)^2)
  preds$derivc0 = 1 - preds$degrad
  preds$varY = (preds$derivk1)^2 * SIG[1,1] + (preds$derivk2)^2 * SIG[2,2] + (preds$derivk3)^2 * SIG[3,3] + (preds$derivc0)^2 * SIG[4,4] +
  2*preds$derivk1*preds$derivk2 * SIG[1,2] + 2*preds$derivk1*preds$derivk3 * SIG[1,3] + 2*preds$derivk1*preds$derivc0 * SIG[1,4] +
  2*preds$derivk2*preds$derivk3 * SIG[2,3] + 2*preds$derivk2*preds$derivc0 * SIG[2,4] + 2*preds$derivk3*preds$derivc0 * SIG[3,4]
  preds$derivk1 = preds$derivk2 = preds$derivk3 = NULL

# Model Type 1 (reparameterisation and k3 = 0)
 } else if(step1_down_object$user_parameters$reparameterisation == T && step1_down_object$user_parameters$zero_order == T){
  k1 = coef(fit)[1]
  k2 = coef(fit)[2]
  c0 = coef(fit)[3]
  preds$degrad = preds$total_time * exp(k1 - k2/preds$K + k2/Kref)
  preds$conc = c0 - c0 * preds$degrad

	# Model Type 1 derivatives for predictions
  preds$derivk1 = -c0 * preds$degrad
  preds$derivk2 = -c0 * (1/Kref - 1/preds$K) * preds$degrad
  preds$derivc0 = 1 - preds$degrad
  preds$varY = (preds$derivk1)^2 * SIG[1, 1] + (preds$derivk2)^2 *  SIG[2, 2] + (preds$derivc0)^2 * SIG[3, 3] + 2 * preds$derivk1 *
  preds$derivk2 * SIG[1, 2] + 2 * preds$derivk1 * preds$derivc0 *
  SIG[1, 3] + 2 * preds$derivk2 * preds$derivc0 * SIG[2, 3]
  preds$derivk1 = preds$derivk2 = preds$derivc0 = NULL

# Model Type 2 (no reparameterisation and k3 = 0)
} else if(step1_down_object$user_parameters$reparameterisation == F && step1_down_object$user_parameters$zero_order == T){
  k1 = coef(fit)[1]
  k2 = coef(fit)[2]
  c0 = coef(fit)[3]
  preds$degrad = preds$total_time * exp(k1 - k2 / preds$K)
  preds$conc = c0 - c0 * preds$degrad

  # Model Type 2 derivatives for predictions
  preds$derivk1 = -c0 * preds$degrad
  preds$derivk2 = c0 / preds$K * preds$degrad
  preds$derivc0 = 1 - preds$degrad
  preds$varY = (preds$derivk1)^2 * SIG[1,1] + (preds$derivk2)^2 * SIG[2,2] + (preds$derivc0)^2 * SIG[3,3] +
  2*preds$derivk1*preds$derivk2 * SIG[1,2] + 2*preds$derivk1*preds$derivc0 * SIG[1,3] + 2*preds$derivk2*preds$derivc0 * SIG[2,3]
  preds$derivk1 = preds$derivk2 = preds$derivc0 = NULL

# Model Type 3 (reparameterisation and k3 is not zero)
}else if(step1_down_object$user_parameters$reparameterisation == T && step1_down_object$user_parameters$zero_order == F){
  k1 = coef(fit)[1]
  k2 = coef(fit)[2]
  k3 = coef(fit)[3]
  c0 = coef(fit)[4]
  preds$degrad = 1 - ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2 / preds$K + k2 / Kref)))^(1/(1-k3))
  preds$conc = c0 - c0*preds$degrad

	# Model Type 3 derivatives for predictions
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
# Obtain CIs using computed variance
  preds$CI1b = preds$conc - qt((1+confidence_interval)/2, summary(fit)$df[2]) * sqrt(preds$varY)
  preds$CI2b = preds$conc + qt((1+confidence_interval)/2, summary(fit)$df[2]) * sqrt(preds$varY)

# Add the degradation intervals
  preds$degrad_CI_up <- (c0 - preds$CI1b) / c0
  preds$degrad_CI_low <- (c0 - preds$CI2b) / c0

  } else {
# Obtain PIs using computed variance
  preds$PI1b = preds$conc - qt((1+confidence_interval)/2, summary(fit)$df[2]) * sqrt(preds$varY + sigma^2)
  preds$PI2b = preds$conc + qt((1+confidence_interval)/2, summary(fit)$df[2]) * sqrt(preds$varY + sigma^2)
  
# Add the degradation intervals
  preds$degrad_PI_up <- (c0 - preds$PI1b) / c0
  preds$degrad_PI_low <- (c0 - preds$PI2b) / c0
 }

## ALTERNATIVELY, OBTAIN STATISTICAL INTERVALS BY USING MULTI-T DRAWS
   } else {

# Model Type 4 ("standard" model) predictions
if(step1_down_object$user_parameters$reparameterisation == F && step1_down_object$user_parameters$zero_order == F){

	# Predictition prep
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    k3 = coef(fit)[3]
    c0 = coef(fit)[4]
    draw = step1_down_object$user_parameters$draw

# Predict response
 preds$degrad <- (1 - ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2 / (preds$K))))^(1/(1-k3)))
 preds$conc <- c0 - c0 * preds$degrad

	# Prediction function
  pred_fct = function(coef.fit) {
  degrad = 1 - ((1 - coef.fit[3]) * (1/(1 - coef.fit[3]) - preds$total_time * exp(coef.fit[1] - coef.fit[2] / preds$K)))^(1/(1-coef.fit[3]))
  conc = coef.fit[4] - coef.fit[4]*degrad
  return(conc)              }
      
	# Multi T samples
  rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
  rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)
  res.draw = matrix(nrow = draw, ncol = nrow(preds), byrow = TRUE, apply(rand.coef, 1, pred_fct))

  if(interval == "CI") {
  preds$CI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
  preds$CI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

  preds$degrad_CI_up <- (c0 - preds$CI1b) / c0
  preds$degrad_CI_low <- (c0 - preds$CI2b) / c0
} else {
  sigma.dist = sqrt((DF * sigma^2) / rchisq(n = draw*nrow(preds), df = DF))
  res.draw = res.draw + rnorm(n = draw*nrow(preds), mean = 0, sd = sigma.dist)
  preds$PI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
  preds$PI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

  preds$degrad_PI_up <- (c0 - preds$PI1b) / c0
  preds$degrad_PI_low <- (c0 - preds$PI2b) / c0
}
  
# Model Type 1 (reparameterisation and k3 = 0)
   } else if(step1_down_object$user_parameters$reparameterisation == T && step1_down_object$user_parameters$zero_order == T){
	
    # Predictition prep
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    c0 = coef(fit)[3]
    draw = step1_down_object$user_parameters$draw

  # Predict response
   preds$degrad <- preds$total_time * exp(k1 - k2 / preds$K + k2 / Kref)
   preds$conc <- c0 - c0 * preds$degrad

	# Prediction function
   pred_fct = function(coef.fit) {
   degrad = preds$total_time * exp(coef.fit[1] - coef.fit[2] / preds$K + coef.fit[2] / Kref)
   conc = coef.fit[3] - coef.fit[3]*degrad
   return(conc)              }
      
	# Multi T samples
    rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
    rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)
    res.draw = matrix(nrow = draw, ncol = nrow(preds), byrow = TRUE, apply(rand.coef, 1, pred_fct))

    if(interval == "CI"){
    preds$CI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
    preds$CI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

    preds$degrad_CI_up <- (c0 - preds$CI1b) / c0
    preds$degrad_CI_low <- (c0 - preds$CI2b) / c0
} else {
    sigma.dist = sqrt((DF * sigma^2) / rchisq(n = draw*nrow(preds), df = DF))
    res.draw = res.draw + rnorm(n = draw*nrow(preds), mean = 0, sd = sigma.dist)
    preds$PI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
    preds$PI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

    preds$degrad_PI_up <- (c0 - preds$PI1b) / c0
    preds$degrad_PI_low <- (c0 - preds$PI2b) / c0
}
    
# Model Type 2 (no reparameterisation and k3 = 0)
   }else if(step1_down_object$user_parameters$reparameterisation == F && step1_down_object$user_parameters$zero_order == T){
	
    # Predictition prep
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    c0 = coef(fit)[3]
    draw = step1_down_object$user_parameters$draw

# Predict response
  preds$degrad <- preds$total_time * exp(k1 - k2 / preds$K)
  preds$conc <- c0 - c0 * preds$degrad

	# Prediction function
    pred_fct = function(coef.fit) {
    degrad = preds$total_time * exp(coef.fit[1] - coef.fit[2] / preds$K)
    conc = coef.fit[3] - coef.fit[3] * degrad
    return(conc)              }
      
	# Multi T samples
    rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
    rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)
    res.draw = matrix(nrow = draw, ncol = nrow(preds), byrow = TRUE, apply(rand.coef, 1, pred_fct))

    if(interval == "CI"){
    preds$CI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
    preds$CI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

    preds$degrad_CI_up <- (c0 - preds$CI1b) / c0
    preds$degrad_CI_low <- (c0 - preds$CI2b) / c0
} else {
    sigma.dist = sqrt((DF * sigma^2) / rchisq(n = draw*nrow(preds), df = DF))
    res.draw = res.draw + rnorm(n = draw*nrow(preds), mean = 0, sd = sigma.dist)
    preds$PI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
    preds$PI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

    preds$degrad_PI_up <- (c0 - preds$PI1b) / c0
    preds$degrad_PI_low <- (c0 - preds$PI2b) / c0
}
    
# Model Type 3 (reparameterisation and k3 is not zero)
   }else if(step1_down_object$user_parameters$reparameterisation == T && step1_down_object$user_parameters$zero_order == F){

    # Predictition prep
    k1 = coef(fit)[1]
    k2 = coef(fit)[2]
    k3 = coef(fit)[3]
    c0 = coef(fit)[4]
    draw = step1_down_object$user_parameters$draw

   # Predict response
   preds$degrad <- 1 - ((1 - k3) * (1/(1 - k3) - preds$total_time * exp(k1 - k2 / preds$K + k2 / Kref)))^(1/(1-k3))
   preds$conc <- c0 - c0 * preds$degrad

	# Prediction function
   pred_fct = function(coef.fit) {
   degrad = degrad = 1 - ((1 - coef.fit[3]) * (1/(1 - coef.fit[3]) - preds$total_time * exp(coef.fit[1] - coef.fit[2] / preds$K + coef.fit[2] / Kref)))^(1/(1-coef.fit[3]))
   conc = coef.fit[4] - coef.fit[4] * degrad
   return(conc)              }
      
	# Multi T samples
   rand.coef = matrix(nrow = n.params, ncol = draw, rnorm(n = n.params * draw, mean = 0, sd = 1))
   rand.coef = t(coef(fit) + t(chol(SIG * DF / (DF - 2))) %*% rand.coef)
   res.draw = matrix(nrow = draw, ncol = nrow(preds), byrow = TRUE, apply(rand.coef, 1, pred_fct))

   if(interval == "CI"){
   preds$CI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
   preds$CI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

   preds$degrad_CI_up <- (c0 - preds$CI1b) / c0
   preds$degrad_CI_low <- (c0 - preds$CI2b) / c0
} else {
   sigma.dist = sqrt((DF * sigma^2) / rchisq(n = draw*nrow(preds), df = DF))
   res.draw = res.draw + rnorm(n = draw*nrow(preds), mean = 0, sd = sigma.dist)
   preds$PI1b = apply(res.draw, 2, quantile, ((1-confidence_interval)/2), na.rm = TRUE)
   preds$PI2b = apply(res.draw, 2, quantile, ((1+confidence_interval)/2), na.rm = TRUE)

   preds$degrad_PI_up <- (c0 - preds$PI1b) / c0
   preds$degrad_PI_low <- (c0 - preds$PI2b) / c0 }
     }      }

# Calculate the RL with CIs or PIs, as specified by the user
if(interval == "CI"){
  degrad = as.numeric(preds$degrad_CI_up[preds$total_time == shelf_time]) 
} else {
  degrad = as.numeric(preds$degrad_PI_up[preds$total_time == shelf_time])
}
  
RL = LSL / (1 - degrad)

# Save RL info as a df
RL_info = data.frame(
"Information" = c("step1_down_object", "shelf_temp", "shelf_time", "LSL", "interval", "release_limit"),
"Detail" = c(deparse(substitute(step1_down_object)), shelf_temp, shelf_time, LSL,
paste0(confidence_interval*100,"% ", interval), round(RL,2)))

# Display the RL printout in console window
cat("Release Limit Information Table:\n\n")
print(RL_info)
cat("\n\n")

## Plotting RL info
# Select which type of RL to make predictions for [using CI or PIs] and
# plot predictions in figure with the selected kind of intervals

# Add predictions using RL as c0, with CIs and PIs
preds$concRL <- RL - RL * preds$degrad

if(interval == "CI"){
preds$concRL_CI_up <- RL - RL * preds$degrad_CI_low
preds$concRL_CI_low <- RL - RL * preds$degrad_CI_up
} else {
preds$concRL_PI_up <- RL - RL * preds$degrad_PI_low
preds$concRL_PI_low <- RL - RL * preds$degrad_PI_up
}

# Reorganise preds df for figure and users
preds$Celsius <- as.factor(preds$temps)

if(interval == "CI"){
predsFit = data.frame(Fit = "Fit",
			time = preds$total_time,
			shelf_temp = preds$Celsius,
			K = preds$K,
			Degradation = preds$degrad,
			Response = preds$conc,
			CI1 = preds$CI1b,
			CI2 = preds$CI2b)
RL_text = paste0("Release Limit = ", round(RL, 2))
predsRL = data.frame(Fit = RL_text,
			time = preds$total_time,
			shelf_temp = preds$Celsius,
			K = preds$K,
			Degradation = preds$degrad,
			Response = preds$concRL,
			CI1 = preds$concRL_CI_low,
			CI2 = preds$concRL_CI_up)
} else {
  
predsFit = data.frame(Fit = "Fit",
    time = preds$total_time,
    shelf_temp = preds$Celsius,
    K = preds$K,
    Degradation = preds$degrad,
    Response = preds$conc,
    PI1 = preds$PI1b,
    PI2 = preds$PI2b)
RL_text = paste0("Release Limit = ", round(RL, 2))
predsRL = data.frame(Fit = RL_text,
 time = preds$total_time,
 shelf_temp = preds$Celsius,
 K = preds$K,
 Degradation = preds$degrad,
 Response = preds$concRL,
 PI1 = preds$concRL_PI_low,
 PI2 = preds$concRL_PI_up)  
}

preds = rbind(predsFit, predsRL)

# Obtain raw data with exactly the same temp as shelf_temp (if present)
# and exclude any of these datapoints which go beyond the shelf_time predictions
if(shelf_temp %in% step1_down_object$data$Celsius){
dat = subset(step1_down_object$data, Celsius == shelf_temp)
    if(max(dat$time) > shelf_time){
      dat = subset(dat, time <= shelf_time) } }                  

# Recognise whether any remaining data points are validation data; make sure these have different shapes in plot
if(exists("dat")){
    if(!is.null(dat$validation)){
      if("Validation" %in% dat$validation){
      shape_types <- c(16, 1)
      names(shape_types) <- c("Fit", "Validation") }} else { validation = NULL } } 

 # Prepatation for the figure and aesthetics
  lines_c <- c("black","blue")
  names(lines_c) <- c("Fit", RL_text)

if(interval == "CI"){
  lines_t <- c("dotted", NA)
  names(lines_t) <- c("Fit", RL_text)
}else{
  lines_t <- c("longdash", NA)
  names(lines_t) <- c("Fit", RL_text)  
}
  
  lines_w = c(0.5, 1)
  names(lines_w) <- c("Fit", RL_text)
  
  ribbons_f = c(NA, "blue")
  names(ribbons_f) = c("Fit", RL_text)

  ribbons_a = c(0, 0.2)
  names(ribbons_a) = c("Fit", RL_text)

if (is.null(xname)){
    xname = "Time"}
if (is.null(yname)){
    yname = "Response Variable"}

mytheme <- ggplot2::theme(legend.position = "bottom", strip.background = element_rect(fill = "white"),
                            legend.key = element_rect(fill = "white"), legend.key.width = unit(2,"cm"),
                            axis.text = element_text(size = 13), axis.title = element_text(size = 13),
                            strip.text = element_text(size = 13),
                            legend.text = element_text(size = 13),
                            legend.title = element_text(size = 13))

# Figure showing predictions with RL and CIs for the given shelf temperature
  if(interval == "CI"){
    plot1 = ggplot2::ggplot() + 
      {if(exists("dat")){geom_point(data = dat, mapping = aes(x= time, y = y, shape = validation), colour = lines_c["Fit"])}} +
      ggtitle(paste0("Release Limit Estimation Using ",confidence_interval*100,"% Confidence Interval,\nShelf Temperature = ", shelf_temp,"°C And Lower Specification Limit = ",LSL))+
      labs( x = xname, y = yname) +
      {if(!is.null(xlim)){coord_cartesian(xlim = xlim)} else { coord_cartesian(xlim = c(0, shelf_time))} } +
      {if(!is.null(ylim)){coord_cartesian(ylim = ylim)}} +
      mytheme +
      geom_hline(yintercept = LSL, linetype = "solid", colour = "red", linewidth = 0.5) +
      geom_line(data = preds, aes(x = time, y = Response, colour = Fit, linewidth = Fit), linetype = "solid") +
      geom_ribbon(data = preds, aes(x = time, ymin = CI1, ymax = CI2, colour = Fit, fill = Fit, alpha = Fit, linetype = Fit)) +
      scale_linetype_manual(name = NULL, values = lines_t) +
      scale_colour_manual(name = NULL, values = lines_c) +
      scale_fill_manual(name = NULL, values = ribbons_f) +
      scale_linewidth_manual(name = NULL, values = lines_w) +
      scale_alpha_manual(name = NULL, values = ribbons_a) +
      {if(exists("dat")){if(!is.null(dat$validation)){if("Validation" %in% dat$validation){scale_shape_manual(values = shape_types, name = NULL) 
        } else {scale_shape(guide = "none")}}}} +
      theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))
   
# Figure showing predictions with PIs and RL with PIs for the given shelf temperature
     }else{
    plot1 = ggplot2::ggplot() + 
    {if(exists("dat")){geom_point(data = dat, mapping = aes(x= time, y = y, shape = validation), colour = lines_c["Fit"])}} +
    ggtitle(paste0("Release Limit Estimation Using ",confidence_interval*100,"% Prediction Interval,\nShelf Temperature = ", shelf_temp,"°C And Lower Specification Limit = ",LSL))+
    labs( x = xname, y = yname) +
    {if(!is.null(xlim)){coord_cartesian(xlim = xlim)} else { coord_cartesian(xlim = c(0, shelf_time))} } +
    {if(!is.null(ylim)){coord_cartesian(ylim = ylim)}} +
    mytheme +
    geom_hline(yintercept = LSL, linetype = "solid", colour = "red", linewidth = 0.5) +
    geom_line(data = preds, aes(x = time, y = Response, colour = Fit, linewidth = Fit), linetype = "solid") +
    geom_ribbon(data = preds, aes(x = time, ymin = PI1, ymax = PI2, colour = Fit, fill = Fit, alpha = Fit, linetype = Fit)) +
    scale_linetype_manual(name = NULL, values = lines_t) +
    scale_colour_manual(name = NULL, values = lines_c) +
    scale_fill_manual(name = NULL, values = ribbons_f) +
    scale_linewidth_manual(name = NULL, values = lines_w) +
    scale_alpha_manual(name = NULL, values = ribbons_a) +
    {if(exists("dat")){if(!is.null(dat$validation)){if("Validation" %in% dat$validation){scale_shape_manual(values = shape_types, name = NULL) 
    } else {scale_shape(guide = "none")}}}} +
    theme(legend.box = "vertical", legend.spacing = unit(-0.4,"line"))
  }

# Display plot in console
print(plot1)

# Save RL info, predictions and plot in a list object
  results <- list(RL_info, plot1, preds)
  names(results) <- c("RL_table", "RL_plot", "RL_predictions")

  return(results)

}

globalVariables(c('total_time','degrad_CI_up','degrad_PI_up','temps','conc',
                  'concRL','CI1b','CI2b','concRL_low','concRL_up','PI1b','PI2b'))
