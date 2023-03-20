#' Return weighted averages of slope soil object
#'
#' @param slp_sol A slp_sol object
#' @return a data frame with 3 columns with the output of the 3 methods.
#' @export
#' @examples
#' fpath_sol <- system.file("extdata", "071000090603_2.sol", package="WEPPR")
#' fpath_slp <-system.file("extdata", "071000090603_2.slp", package = "WEPPR")
#' slp <- read_slp(fpath_slp)
#' sol <- read_sol(fpath_sol)
#' get_all_wavg_sol(slp, sol)

get_all_wavg_sol <- function(slp, sol) {

    slp_sol <- WEPPR:::merge_slp_sol(slp, sol)

    # impute missing numerical columns with 0
    slp_sol <- slp_sol %>% replace(is.na(.), 0)

    return(
        data.frame(
            wavg1 = get_wavg_sol(slp_sol, w_dist = 0, w_thk = 0),
            wavg2 = get_wavg_sol(slp_sol, w_dist = 1, w_thk = 1),
            wavg3 = get_wavg_sol(slp_sol, w_dist = 2, w_thk = 1),
            wavg4 = get_wavg_sol(slp_sol, w_dist = 1, w_thk = 2)
        )
    )
}

#' Return weighted average of x
#'
#' @param x vector of slopes
#' @param w vector of distances
#' @return the weighted average of x by w
#' @examples
#' x <- c(0, 0.5, 1)
#' w <- c(1, 2, 3)
#' getRollAvg(x, w)
getRollAvg <- function(x, w) {
    if (length(x) == 1) {
        return(weighted.mean(x, w))
    }
    x_roll <-
        na.omit(data.table::frollmean(x, 2,  hasNA = T, na.rm = T))
    return(weighted.mean(x_roll, diff(w)))
}

#' Return weighted average of slp_sol by solthk
#'
#' @param slp_sol A slp_sol object
#' @return the weighted average of slp_sol
wavg_solthk <- function(slp_sol) {
    solthk_avg <- slp_sol %>%
        group_by(layer) %>%
        summarise(across(
            c(salb:distance, -solthk),
            ~ WEPPemulator:::getRollAvg(.x, w = solthk)
        ))

    return(solthk_avg[-1])
}


#' Return weighted average of slp_sol
#'
#' @param slp_sol A slp_sol object
#' @param w_dist An integer (0, 1, 2) specifying the method to use for
#' weighted average calculation
#' @param w_thk An integer (0, 1, 2) specifying the method to use for
#' weighted average calculation
#' @return the weighted average of slp_sol
get_wavg_sol <- function(slp_sol, w_dist, w_thk) {

    compute_area <- function(cml_x, angle=45) {
        tri_area <- function(adj) {
            opposite <- adj * tan(angle * pi / 180)
            area <- (opposite * adj) / 2
            return(area)
        }

        a <- c()
        for (i in 1:length(cml_x)) {
            a <- c(a, tri_area(cml_x[i]))

            if (i != 1) {
                a[i] <- a[i] - a[i - 1]
            }

        }

        return(a)
    }


    # first row of last layer
    if (w_dist == 0 & w_thk == 0) {
        last_layer <- slp_sol[slp_sol$layer == max(slp_sol$layer), ]
        return(last_layer[1, -c(1, 8, 14, 15)])
    }

    # both 1 (linear)
    if (w_dist == 1 & w_thk == 1) {
        solthk_avg <- wavg_solthk(slp_sol)

        if (nrow(solthk_avg) == 1) {
            return(solthk_avg)
        }

        # weight by distance
        x <- solthk_avg[1:11]
        d <- solthk_avg$distance
        avg <- data.frame(t(colSums(x * d) / sum(d)))
        return(avg)
    }

    # triangle area over distance
    else if (w_dist == 2 & w_thk == 1) {
        solthk_avg <- wavg_solthk(slp_sol)

        solthk_avg$dist_cum <- cumsum(solthk_avg$distance)


        dist <- solthk_avg$dist_cum
        a <- compute_area(dist, angle=45)

        # weight by area
        x <- solthk_avg[1:11]
        avg <- data.frame(t(colSums(x * a) / sum(a)))
        return(avg)
    }

    # triangle over thickness
    else if (w_dist == 1 & w_thk == 2) {
        solthk_avg <- slp_sol %>%
            group_by(layer) %>%
            mutate(area = compute_area(solthk)) %>%
            summarise(across(
                c(salb:distance, -solthk),
                ~ WEPPemulator:::getRollAvg(.x, w = area)
            )) %>%
            select(-layer)

        # weight by distance
        x <- solthk_avg[1:11]
        d <- solthk_avg$distance
        avg <- data.frame(t(colSums(x * d) / sum(d)))

        return(avg)
    }

}





