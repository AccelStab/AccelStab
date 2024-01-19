library("mvtnorm")
library(dplyr)
library(ggplot2)
library(scales)

load("~/CMC - Bayesian/AccelStab/data/antigenicity.rda")
load("~/CMC - Bayesian/AccelStab/data/potency.rda")
test1 <- step1_down_rmse(data = antigenicity, y = "conc", .time = "time", C = "Celsius", parms = list(c0 = c(96,98,100), k1 = c(42,45), k2 = c(12000,12500), k3 = c(8,9,10)))

test2 <- step1_down_rmse(data = antigenicity, y = "conc", .time = "time", C = "Celsius", parms = list(c0 = c(100), k1 = c(40), k2 = c(14000), k3 = c(5,6)))

fit1 <- step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", reparameterisation = T)
fit2 <- step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", reparameterisation = F)

test3 <- step1_down_rmse(data = antigenicity, y = "conc", .time = "time", C = "Celsius", parms = list(c0 = c(100,95), k1 = c(2,2.5), k2 = c(12000,13000), k3 = c(9,10)), reparameterisation = T)
test4 <- step1_down_rmse(data = antigenicity, y = "conc", .time = "time", C = "Celsius", parms = list(c0 = c(100,95), k1 = c(2,2.5), k2 = c(12000,13000), k3 = c(9,10)), reparameterisation = F)


# Now to deal with this validation data situation
# First do it for the step1_plot_desc function

antigenicity <- antigenicity %>% mutate(val = seq(1,50,1))
step1_plot_desc(data=antigenicity, y="conc", .time="time", C = "Celsius", validation = "val")

antigenicity <- antigenicity %>% mutate(val = rep(c(0,1),25))
step1_plot_desc(data=antigenicity, y="conc", .time="time", C = "Celsius", validation = "val")

antigenicity <- antigenicity %>% mutate(val = as.factor(rep(c(0,1),25)))
step1_plot_desc(data=antigenicity, y="conc", .time="time", C = "Celsius", validation = "val")
step1_plot_desc(data=antigenicity, y="conc", .time="time", C = "Celsius")

antigenicity <- antigenicity %>% mutate(val = c(rep(c(0,1),24),NA,NA))
step1_plot_desc(data=antigenicity, y="conc", .time="time", C = "Celsius", validation = "val")

# Now it needs to be added to step_1_down and the returned results$data frame before being put
# into every one of the other plotting functions

antigenicity <- antigenicity %>% mutate(val = rep(c(0,1),25))
fit1 <- step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius")
fit2 <- step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", validation = "val")

step1_plot_CI(fit2, ribbon = T)
step1_plot_PI(fit2, ribbon = T)
step1_plot_pred(fit2)
step1_plot_T(fit2, focus_T = 5, ribbon = F)

# Now for excursion and RL
exc1 <- excursion(fit1,temp_changes = c(5,35,5), time_changes = c(0.5,1,3))
exc1$excursion_plot




