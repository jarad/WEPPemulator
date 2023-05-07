library(WEPPR)
library(dplyr)

fpath_slp <- system.file("extdata", "071000090603_2.slp", package = "WEPPR")
slp <- read_slp(fpath_slp)

test_that("correct weighted average calculation", {
    # First method (original)
    method1_df <- get_all_wavg(slp)

    # Second method (exact)
    slp_extra <- tweak_slp(slp)
    exact_df <- data.frame(
        slp_slope_wavg0 = get_wavg_exact(slp_extra, weight = 0),
        slp_slope_wavg1 = get_wavg_exact(slp_extra, weight = 1)
    )

    # Compare approximate with exact method
    expect_equal(method1_df[, 1:2], exact_df, tolerance = 0.05)
})
