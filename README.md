
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Accelerated Stability Kinetic Modelling

<!-- badges: start -->
<!-- badges: end -->

## Overview

This package utilises the Šesták–Berggren equation alongside the
Arrhenius equation to make a simple and consistent way for a user to
carry out the calculations and predictions required by accelerated
stability studies. Currently the package is only able to work with
decreasing variables, you may choose to transform your increasing
variable into a decreasing one but note that your choice of
transformation can have a large impact on the outcome.

The available functions within the package are as follows:

-   `step1_down()` Fit the one-step Sestak-Berggren kinetic model.
-   `step1_plot_desc()` Plot the stability data.
-   `step1_plot_pred()` Plot the stability data and visualise the
    predictions.
-   `step1_plot_CI()` Plot the stability data and visualise the
    predictions with confidence intervals.
-   `step1_plot_PI()` Plot the stability data and visualise the
    predictions with prediction intervals.
-   `step1_plot_T()` Plot the stability data and visualise the
    predictions with focus on one temperature.
-   `excursion()` Predict a temperature excursion for a product.

## Installation

You can install the development version of AccelStab from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("AccelStab/AccelStab")
```

## Getting Started

A more detailed look at the features and background theory of the
package can be found in the ‘Introduction’ vigentte (ADD LINK HERE),
which gives instructions on how to use the package with the available
functions and their arguments.

## Feedback

Log an issue (add hyoerlink to isseu) or contact Bernard Franq (add
email).
