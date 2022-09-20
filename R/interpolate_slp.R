#' Interpolate the slope data file
#'
#' Converts the slope file into one WEPP run (data frame with one row)
#'
#' @param slp A \code{slp} object
#' @param n (optional) a positive integer indicating how many slope columns to evaluate
#' @seealso \code{\link{expand_slp}}
#' @return a one-row dataframe containing total distance and slope columns
#' @export
#' @examples
#' slp <- read_slp(system.file("extdata", "071000090603_2.slp", package="WEPPR"))
#' lin_slp <- interpolate_slp(slp)
#'
interpolate_slp <- function(slp, n = 1001) {
    remove_trail_patt <- '^(\\.\\d*?[1-9])0+$'

    total_dist <- tail(WEPPR:::calculate_total_distance(slp), n = 1)

    interpolated_slp <- slp %>%
        WEPPR:::expand_slp(n = n) %>%
        select(slope) %>%
        mutate(ID = str_remove(round(1:n / n, digits = 4), remove_trail_patt)) %>%
        pivot_wider(
            names_from = ID,
            names_prefix = "slp_slope_",
            values_from = slope
        )

    int_slp_dist <- cbind(total_dist, interpolated_slp)

    return(int_slp_dist)
}
