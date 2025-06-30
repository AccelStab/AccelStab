## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(AccelStab)

## ----include=FALSE------------------------------------------------------------
library(ggplot2)

## -----------------------------------------------------------------------------
library(AccelStab)
data(antigenicity)
antigenicity$Validate = as.factor(ifelse(antigenicity$time <= 0.5, 0, 1))
head(antigenicity)

## -----------------------------------------------------------------------------
table(antigenicity$Celsius)
unique(antigenicity$time)

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
step1_plot_desc(data = antigenicity, .time = "time", y = "conc", C = "Celsius", validation = "Validate", yname = "Antigenicity")

## -----------------------------------------------------------------------------
args(step1_down)

## -----------------------------------------------------------------------------
res = step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", validation = "Validate")
summary(res$fit)
confint(res$fit)

## -----------------------------------------------------------------------------
res = step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", validation = "Validate", parms = list(k1 = 50, k2 = 10000, k3 = 3, c0 = 100))
summary(res$fit)
confint(res$fit)

## -----------------------------------------------------------------------------
head(res$prediction[,1:5])

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
step1_plot_pred(res, yname = "Antigenicity")

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
res = step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", validation = "Validate", parms = list(k1 = 50, k2 = 10000, k3 = 3, c0 = 100), max_time_pred = 3)
graph = step1_plot_pred(res, yname = "Antigenicity")
graph = graph + geom_hline(aes(yintercept = 65), linetype = "dotted")
graph

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
res = step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", validation = "Validate", parms = list(k1 = 50, k2 = 10000, k3 = 3, c0 = 100), max_time_pred = 3, temp_pred_C = 2)
step1_plot_pred(res, yname = "Antigenicity")


## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
subdat = antigenicity[!(antigenicity$Celsius == "5" & antigenicity$time != 0),]
res = step1_down(data = subdat, y = "conc", .time = "time", C = "Celsius", max_time_pred = 3, temp_pred_C = 5)
step1_plot_pred(res, yname = "Antigenicity")


## -----------------------------------------------------------------------------
res = step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", max_time_pred = 3, validation = "Validate")
head(res$prediction[,-c(3,6:8)])

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
res = step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", max_time_pred = 3, validation = "Validate")
step1_plot_CI(res, yname = "Antigenicity")
step1_plot_PI(res, yname = "Antigenicity", ribbon = TRUE)
step1_plot_T(res, focus_T = 5, yname = "Antigenicity", ribbon = TRUE)

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
res = step1_down(data = antigenicity, y = "conc", .time = "time", C = "Celsius", max_time_pred = 3, validation = "Validate")
exc <- excursion(res, temp_changes = c(5,35,6), time_changes = c(6/12,7/12,24/12), yname = "Antigenicity")
tail(exc$predictions[,-c(4)])
exc$excursion_plot

## -----------------------------------------------------------------------------
draws = step1_sample_mvt(data = antigenicity, y = "conc", .time = "time", C = "Celsius", validation = "Validate", draw = 10^4)
draws = as.data.frame(draws)
head(draws)

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
p1 = draws$c0 - draws$c0 * (1 - ((1 - draws$k3) * (1/(1 - draws$k3) - 1 * exp(draws$k1 - draws$k2 / (5+273.15))))^(1/(1-draws$k3)))
loss_1y = 1 - p1/draws$c0
mean(loss_1y)
quantile(loss_1y, c(0.025, 0.975))
hist(loss_1y, main = "Lost (%) in 1 year at 5C")
abline(v = mean(loss_1y), lwd = 2, col = 3)
abline(v = quantile(loss_1y, c(0.025, 0.975)), lwd = 2, col = 3, lty = 2)

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
data(potency)
potency$Validate = as.factor(ifelse(potency$Time < 8, 0, 1))
head(potency)
step1_plot_desc(data = potency, .time = "Time", y = "Potency", C = "Celsius", validation = "Validate", yname = "Potency")

## -----------------------------------------------------------------------------
res = step1_down(data = potency, y = "Potency", .time = "Time", C = "Celsius", validation = "Validate")
summary(res$fit)
confint(res$fit)

## -----------------------------------------------------------------------------
res = step1_down(data = potency, y = "Potency", .time = "Time", C = "Celsius", validation = "Validate", zero_order = TRUE)
summary(res$fit)
confint(res$fit)

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
step1_plot_diagnostic(res)[1]

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
RMSE = step1_down_rmse(data = potency, y = "Potency", .time = "Time",
 C = "Celsius", parms = list(c0 = 9.5, k1 = seq(38, 42, 0.02),
  k2 = seq(12000, 14000, 5), k3 = 0))
RMSE$logrmse = log(RMSE$rmse)

plot = ggplot() + geom_point(data=RMSE, mapping=aes(x= k1, y = k2, colour = logrmse), size = 1.5, stroke = 0)
plot = plot + labs( x = "k1", y = "k2")
plot = plot + scale_color_gradient2(midpoint = mean(RMSE$logrmse), low = "blue", mid = "yellow", high = "green")
plot

## ----fig.width=8, fig.height=6, out.width='75%',dpi=200, fig.align='center', warning=FALSE----
head(res$prediction[,-(6:8)])
plot = step1_plot_pred(res, yname = "Potency")
plot = plot + scale_y_continuous(limits = c(5,10))
plot
plot = step1_plot_T(res, focus_T = 5, yname = "Potency", ribbon = TRUE)
plot = plot + scale_y_continuous(limits = c(5,10))
plot

## -----------------------------------------------------------------------------
draws = step1_sample_mvt(data = potency, y = "Potency", .time = "Time", C = "Celsius", validation = "Validate", draw = 10^4, zero_order = TRUE)
draws = as.data.frame(draws)
head(draws)
T9 = (1 - 9/draws$c0) / exp(draws$k1 - draws$k2 / (5+273.15))

## ----fig.width=8, fig.height=6, out.width='75%', dpi=200, fig.align='center'----
mean(T9)
quantile(T9, c(0.025, 0.975))
hist(T9, main = "Time to reach the lower specification at 5C")
abline(v = mean(T9), lwd = 2, col = 3)
abline(v = quantile(T9, c(0.025, 0.975)), lwd = 2, col = 3, lty = 2)

