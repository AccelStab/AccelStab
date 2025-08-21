## CRAN submission
Resubmission to CRAN that fixes below issue found in CRAN check.

Error(s) in re-building vignettes:
--- re-building ‘AccelStab.Rmd’ using rmarkdown

Quitting from AccelStab.Rmd:214-222 [unnamed-chunk-18]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
<error/rlang_error>
Error in `quantile.default()`:
! missing values and NaN's not allowed if 'na.rm' is FALSE
---
Backtrace:
    ▆
 1. ├─stats::quantile(loss_1y, c(0.025, 0.975))
 2. └─stats:::quantile.default(loss_1y, c(0.025, 0.975))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Error: processing vignette 'AccelStab.Rmd' failed with diagnostics:
missing values and NaN's not allowed if 'na.rm' is FALSE
--- failed re-building ‘AccelStab.Rmd’

SUMMARY: processing the following file failed:
  ‘AccelStab.Rmd’

Error: Vignette re-building failed.

## R CMD check results

0 errors | 0 warnings | 0 note


