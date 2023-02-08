#' interpolate the soil data file
#'
#' Converts the soil file into one WEPP run (data frame with one row)
#'
#' @param slp_sol A Soil object combined with slope data
#' @return a one-row dataframe containing soil data
#' @export
#' @examples
#' slp_sol <- merge_slp_sol(slp, sol)
#' int_sol <- interpolate_sol(slp_sol)
#'
interpolate_sol <- function(slp_sol) {

    if (!"distance" %in% colnames(slp_sol))
    {
        cat("Please perform merge_slp_sol before interpolation")

    }

    slp_sol %<>%
        mutate(cum_dist = WEPPR:::calculate_cum_dist(layer, distance)) %>%
        mutate(norm_cum_dist = cum_dist / max(cum_dist))

    # create solthk bin vector
    bins <-
        tibble(solthk = cbind(seq(0, 1800, by = 100))) %>% pull()

    # filter numeric values
    filter_num <- slp_sol %>%
        select_if(is.numeric) %>%
        select(-slope, -distance, -cum_dist)

    # duplicate rows based on count
    dups_sol <- filter_num %>%
        group_by(layer) %>%
        group_map( ~ WEPPemulator:::get_dup_count(.x)) %>%
        unlist() %>%
        cbind(freq = ., filter_num) %>%
        mutate(freq = map(freq, seq_len)) %>%
        unnest(freq) %>%
        select(-freq)

    # interpolate data
    int_sol <- dups_sol %>%
        group_by(layer) %>%
        mutate(solthk_bin = bins, .before = sand) %>%
        ungroup() %>%
        mutate(id = row_number()) %>%
        group_by(layer) %>%
        pivot_wider(
            names_from = c(layer, id, solthk_bin, norm_cum_dist),
            names_glue = "sol_{.value}_{round(norm_cum_dist, 3)}_{solthk_bin}",
            values_from = c(salb:rfg)
        ) %>%
        summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))

    return(int_sol)
}

#' Calculate the duplicate count for each row based on solthk bin for interpolating
#' soil file
#'
#' @param sol A soil object
#' @return A numeric vector of duplicate counts
#' @examples
#' dup_counts <- get_dup_count(sol)
#' dup_counts
#'
get_dup_count <- function(sol) {
    dup <- c()
    j = 0

    # loop through each value
    for (i in 1:nrow(sol)) {
        d = 0

        # while it's greater than bin value, increment
        while (sol[i, ]$solthk > j) {
            d = d + 1
            j = j + 100
        }
        dup <- c(dup, d)
    }

    # calculate duplicates of the final bin (19th)
    left <- 19 - sum(dup)

    # add it to the final duplicate value
    dup[length(dup)] <- dup[length(dup)] + left

    return(dup)
}
