library("mvtnorm")
library(dplyr)

load("~/CMC - Bayesian/AccelStab/data/antigenicity.rda")
load("~/CMC - Bayesian/AccelStab/data/potency.rda")
test1 <- step1_down_rmse(data = antigenicity, y = "conc", .time = "time", C = "Celsius", parms = list(c0 = c(96,98,100), k1 = c(42,45), k2 = c(12000,12500), k3 = c(8,9,10)))

test2 <- step1_down_rmse(data = antigenicity, y = "conc", .time = "time", C = "Celsius", parms = list(c0 = c(100), k1 = c(40), k2 = c(14000), k3 = c(5,6)))

fit1 <- step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", reparameterisation = T)
fit2 <- step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", reparameterisation = F)

test3 <- step1_down_rmse(data = antigenicity, y = "conc", .time = "time", C = "Celsius", parms = list(c0 = c(100,95), k1 = c(2,2.5), k2 = c(12000,13000), k3 = c(9,10)), reparameterisation = T)
test4 <- step1_down_rmse(data = antigenicity, y = "conc", .time = "time", C = "Celsius", parms = list(c0 = c(100,95), k1 = c(2,2.5), k2 = c(12000,13000), k3 = c(9,10)), reparameterisation = F)
