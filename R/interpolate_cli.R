#' Perfoms a pivot wide on ppt inst values rows based on a single date value
#'
#' First it performs an expansion on the cli object and the pivot wides the data
#'
#' @param cli a Climate object
#' @param n (optional) a positive integer indicating how many equally-spaced time to
#'     calculate.
#' @return a data frame with only one row and n columns of cumulative ppt
#' @export
#'
pivot_wide_cli_by_date <- function(cli, n = 1001) {
    remove_trail_patt <- '^(\\.\\d*?[1-9])0+$'

    widened_cli <- cli %>%
        WEPPR:::expand_cli(n = n) %>%
        select(ppt, date) %>%
        mutate(ID = str_remove(round(1:n / n, digits = 4), remove_trail_patt)) %>%
        pivot_wider(
            names_from = ID,
            names_prefix = "cli_ppt_",
            values_from = ppt
        )

    return(widened_cli)
}

#' Interpolates the climate data frame
#'
#' @param cli a Climate object
#' @param n (optional) a positive integer indicating how many equally-spaced time to
#'     calculate.
#' @return a data frame where rows are unique dates and columns are interpolated
#' ppt values.
#' @export
#' @seealso \code{\link{pivot_wide_cli_by_date}}
#' @examples
#' cli <- read_cli(system.file("extdata", "092.63x040.90.cli", package="WEPPR"))
#' interpolate_cli(cli)
#'
interpolate_cli <- function(cli, filename, n = 1001) {
    cli_ppt <- cli$breakpoints
    all_interpolated <- cli_ppt %>%
        group_by(date) %>%
        group_map( ~ pivot_wide_cli_by_date(.x, n), .keep = TRUE) %>%
        bind_rows() %>%
        mutate(date = paste0(filename, date))

    return(all_interpolated)
}
