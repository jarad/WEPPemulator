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
    slp_extra <- tweak_slp(slp)

    return(data.frame(
        slp_slope_wavg0 = get_wavg(slp_extra, weight = 0),
        slp_slope_wavg1 = get_wavg(slp_extra, weight = 1),
        slp_slope_wavg2 = get_wavg(slp_extra, weight = 2)
    ))
}

#' Add new columns to slp object for weighted average calculation
#'
#' Creates the following columns:
#' - dist_diff (the difference between total distances)
#' - avg_slp (rolling average with window size of 2)
#' - dist_max (the total distance of each OFE)
#'
#' @param slp A slp object
#' @return a data frame with 3 additional distance and slope columns
tweak_slp <- function(slp) {
    class(slp) <- "data.frame"
    slp$total_dist <- WEPPR:::calculate_total_distance(slp)

    slp_extra <- slp %>%
        mutate(
            dist_diff = total_dist - lag(total_dist, default = 0),
            avg_slp = zoo::rollmean(
                slope,
                k = 2,
                fill = 0,
                align = "right"
            )
        ) %>%
        group_by(n) %>%
        mutate(dist_max = max(distance))

    return(slp_extra)
}

#' Return weighted average of slope
#'
#' @param s vector of slopes
#' @param d vector of distances
#' @return the weighted average of slope by distance
#' @examples
#' s <- c(0, 0.5, 1)
#' d <- c(1, 2, 3)
#' wavg(s, d)
wavg <- function(s, d) {
    avg <- sum(s * d) / sum(d)

    return(avg)
}

#' Return weighted average of slope based on weight parameter
#'
#' @param slp_extra A slp object with additional columns
#' @param weight An integer value (0, 1, 2) specifying the method to use for
#' weighted average calculation
#' @return the weighted average of slope based on the specified weight
#' @examples
#' fpath_slp <- system.file("extdata", "071000090603_2.slp", package = "WEPPR")
#' slp <- read_slp(fpath_slp)
#' slp_extra <- tweak_slp(slp)
#' get_wavg(slp_extra, weight = 0)
#' get_wavg(slp_extra, weight = 1)
#' get_wavg(slp_extra, weight = 2)
get_wavg <- function(slp_extra, weight) {
    if (weight == 0) {
        avg <- slp_extra$slope[nrow(slp_extra)]
    } else if (weight == 1) {
        s <- slp_extra$avg_slp
        d <- slp_extra$dist_diff
        avg <- wavg(s, d)
    } else if (weight == 2) {
        tri_area <- function(adj) {
            opposite <- adj * tan(angle * pi / 180)
            area <- (opposite * adj) / 2
            return(area)
        }

        angle <- 45
        df <- slp_extra %>%
            select(n, total_dist) %>%
            group_by(n) %>%
            summarize(cnt = n(), max_dist = max(total_dist))

        count <- df$cnt
        dist <- df$max_dist

        s <- slp_extra$avg_slp
        a <- c()

        for (i in 1:length(dist)) {
            a <- c(a, tri_area(dist[i]))

            if (i != 1) {
                a[i] <- a[i] - a[i - 1]
            }
        }
        a <- rep(a, count)

        avg <- sum(s * a) / sum(a)
    } else {
        stop("Invalid weight value. Please use 0, 1, or 2.")
    }

    return(avg)
}
