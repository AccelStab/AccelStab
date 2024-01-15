
step1_down_rmse <- function (data, y, .time, K = NULL, C = NULL,
                        parms, reparameterisation = FALSE){

  if (is.null(K) & is.null(C))
    stop("Select the temperature variable in Kelvin or Celsius")
  if (!is.list(parms))
    stop("The starting values for parameters must be a list")

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

  if (is.null(K))
    dat$K = dat[, C] + 273.15
  if (is.null(C)) {
    dat$C = dat[, K] - 273.15
    C = "C"}

  Kref = mean(dat$K)
  dat$Celsius = as.factor(dat[, C])
  dat$time = dat[, .time]
  dat$y = dat[, y]
  if(.time != "time"){
    dat <- dat[, !names(dat) %in% c(.time)]
  }
  if(y != "y"){
    dat <- dat[, !names(dat) %in% c(y)]
  }

  # I think the way to do it is make an expand grid of all the combos then add a rmse column on the end of that!

  result_grid <- expand.grid(parms) %>% mutate(rmse = NA)

  if(reparameterisation){

  }else{
    for (i in 1:dim(result_grid)[1]){
      c0 <- result_grid[i,]$c0
      k1 <- result_grid[i,]$k1
      k2 <- result_grid[i,]$k2
      k3 <- result_grid[i,]$k3

      dat$Degradation = 1 - ((1 - k3) * (1/(1 - k3) - dat$time * exp(k1 - k2 / dat$K)))^(1/(1-k3))
      dat$Response = c0 - c0*dat$Degradation
      dat$sqrResidual = (dat$Response - dat$y)^2

      result_grid[i,'rmse'] <- sqrt(mean(dat$sqrResidual))
    }
  }

  return(result_grid)

}
