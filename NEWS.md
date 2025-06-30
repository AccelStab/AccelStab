# AccelStab 2.3.1
This is a minor release that changes the following...

### step1_down() and step1_sample_mvt()
- Replaced sigma with chisq for multi-temperature models
- The residual error standard deviation (sigma) is now sampled from its own distribution rather than being treated as fixed.

### step1_sample_mvt()
- Updated to align with step1_down()
- Added ... (ellipsis) to support optional arguments passed to minpack.lm::nlsLM() during model fitting.

### excursion()
- Updated procedure for excursion() to draw from multi-t dist like used in step1_down() and step1_sample_mvt().

# AccelStab 2.2.1

- Fixed a minor plot issue in the vignette.

# AccelStab 2.2.0

This is a minor release which adds extra functionality to existing functions plus
a vignette has been added that gives more background to the package.

## Enhancements

### `step1_plot_diagnostic()`
- Includes **validation data** in observed vs predicted and residuals vs predicted scatter plots, if provided.
- Adds a new plot: **residuals vs time**, color-coded by temperature, and includes validation data points when present.
- Vignette updated to reflect that the function now generates **five plots** (previously four).

### Plotting Functions (`step1_plot_CI()`, `step1_plot_desc()`, `step1_plot_PI()`, `step1_plot_pred()`, `step1_plot_T()`)
- `xlim` and `ylim` now **truncate the visible plot area only**, without removing data used in the plot.

## Bug Fixes

### `step1_plot_diagnostic()`
- Fixed **QQ plot** to show quantiles for all data, rather than temperature-specific quantiles.

### `step1_down()` and `step1_down_basic()`
- Added `...` to allow passing optional arguments to `minpack.lm`.
- When `reparam = TRUE`, now prints a message: **default lower bounds are 0**, but custom bounds may be specified.
- Fixed issue with **Kelvin temperature** inputs.

## Documentation

- Updated **function descriptions** across all relevant functions for clarity and completeness.

# AccelStab 2.1.1

This is a minor release updating authors and email addresses

# AccelStab 2.1.0

This is a major release adding new functionality, adjusting the data sets and documentation.

## New features

-   `step1_down_basic()` This new function allows for rapid testing of the fit without 
    all of the features.

-   `antigenicity` Changed example data and added a validation column.

-   `step1_down()` Adjusted the help page and added validation options to the examples.
    The function will now also return sampled parameters, reducing the need for the 
    step1_sample_mvt function.

-   `step1_plot_diagnostic()` Added argument to use either classic, standardized or 
    studentized residuals.


# AccelStab 2.0.2

This release has a minor bug fix and additional descriptions added.

## Bug fixes

-   `step1_sample_mvt()` Fixed issue when no values at time point zero

-   `potency.rda` File data adjusted slightly

-   `step1_down()` Added to some arguments' description and new print 
    when any of the k3 draws are below zero.


# AccelStab 2.0.0

This is a major release adding a few new functions and fixing a bug.

## New features

-   `step1_down_rmse()` added which allows the user to calculate the
    root mean squared error for their data and chosen parameters

-   `step1_down_diagnostic()` added which allows the user to plot
    residual diagnostic plots after fitting the model

-   `step1_sample_mvt()` added which allows the user to draw a chosen
    number of sample parameters from the multivariate t distribution for
    their own analyses
    
-   `step1_down()` now accepts an extra argument `validation` which
    sidelines some of the data allowing the user to save it for testing
    purposes

-   `step1_plot_desc()` now accepts an extra argument `validation` which
    sidelines some of the data allowing the user to save it for testing
    purposes

## Bug fixes

-   When selecting a temperature outside the data set for
    `step1_plot_T()` the colours of the prediction line and the ribbon
    are now consistent
-   When using the argument `temp_pred_C` within `step1_down()` no longer
    are predictions duplicated if the temperature is already in the data
-   Issue when no time = 0 rows present in the data and no `parms` provided
    to `step1_down()` now resolved

# AccelStab 1.0.0

-   Added a `NEWS.md` file to track changes to the package.
