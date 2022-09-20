
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WEPPemulator

## Overview

WEPPR an R package that has functionality to construct a statistical
emulator for the Water Erosion Prediction Project

The [water erosion prediction
project](https://www.fs.usda.gov/ccrc/tool/watershed-erosion-prediction-project-wepp)
(WEPP) is a physically-based soil erosion computer model supported by
the United States Department of Agriculture (USDA).

For reading and writing input and output files as well as running WEPP
on Linux systems, check out [WEPPR](https://github.com/jarad/WEPPR)

## Installation

You can install the development version of WEPPemulator from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jarad/WEPPemulator")
```

## Example

### Read the data

``` r
library(WEPPR)
library(WEPPemulator)
library(tidyverse)

fpath_slp <- system.file("extdata", "071000090603_2.slp", package="WEPPR")
fpath_sol <- system.file("extdata", "071000090603_2.sol", package="WEPPR")
fpath_cli <- system.file("extdata", "092.63x040.90.cli", package="WEPPR")

slp <- read_slp(fpath_slp)
sol <- read_sol(fpath_sol)
cli <- read_cli(fpath_cli)

interpolate_slp(slp)

slp_sol <- merge_slp_sol(slp, sol)
interpolate_sol(slp_sol)

interpolate_cli(cli, filename = "092.63x040.90")
```
