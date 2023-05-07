#' Return 3 weighted averages of slope
#'
#' @param slp A slp object
#' @return a data frame with 3 columns with the output of the 3 methods.
#' @export
#' @examples
#' fpath_slp <- system.file("extdata", "071000090603_2.slp", package = "WEPPR")
#' slp <- read_slp(fpath_slp)
#' get_all_wavg(slp)

get_all_wavg <- function(slp) {

    max_dist <- max(WEPPR:::calculate_total_distance(slp))
    b <- 10

    return(data.frame(
        slp_slope_wavg0 = get_wavg(slp, Vectorize(function(t) 1)),
        slp_slope_wavg1 = get_wavg(slp, function(t) 1 - (t / max_dist)),
        slp_slope_wavg2 = get_wavg(slp, function(t) b^t)
    ))
}

#' Calculate the weighted average slope of a slope-length product data frame
#'
#' This function calculates the weighted average slope of a slope-length product
#' (SLP) data frame using a weight function. The weight function determines
#' the contribution of each slope value to the weighted average, and can be
#' specified by the user.
#'
#' @param slp A data frame containing the slope and length values for each
#'   segment of a hillslope. The data frame should have two columns, "x" for
#'   segment length and "slope" for segment slope.
#' @param w A weight function that specifies the contribution of each slope
#'   value to the weighted average. The function should take a single argument,
#'   the distance along the hillslope, and return a weight value for that
#'   distance. The weight function should be non-negative and integrable over
#'   the length of the hillslope.
#'
#' @return The weighted average slope of the SLP data frame.
#'
#' @examples
#' fpath_slp <-system.file("extdata", "071000090603_2.slp", package = "WEPPR")
#' slp <- read_slp(fpath_slp)
#' w_const <- function(t) 1
#' get_wavg(slp, Vectorize(w_const))
#'
get_wavg <- function(slp, w) {
    slp_d <- slp %>%
        WEPPR:::remove_slp_transitions() %>%
        mutate(cuml_dist = WEPPR:::calculate_total_distance(.)) %>%
        mutate(norm_cuml_dist = cuml_dist / max(cuml_dist))

    numerator <- integrate(function(t) {
        y_t <- approx(slp_d$norm_cuml_dist, slp_d$slope, xout = t)$y
        y_t * w(t)
    }, lower = 0, upper = 1)$value

    denominator <- integrate(function(t) {
        w(t)
    }, lower = 0, upper = 1)$value

    weighted_slope <- numerator / denominator

    return(weighted_slope)
}
