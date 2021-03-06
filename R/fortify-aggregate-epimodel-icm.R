## March 3, 2020
# Comments for Ben and Shannon
#
# This is the start of a reformat of fortify_epimodel and fortify_pomp into one
# function "fortify_aggregate".  This will take in aggregate level data from
# external sources (currently pomp and epimodel) and output, if possible,
# output a data.frame  with the columns t, X0, .., XK where Xk are the state
# names and the values in Xk are the numbers in state k at time t.





#'  Take external aggregate data and put it in a format used in this package
#'
#' @param data output from external source package.  See details
#' @param states names of states we want aggregate totals of at each time
#' @param package_source optional argument to include the package from which the
#'   output is derived from, which helps with the fortify function when outputs
#'   are of generic classes such as list or data.frame
#' @return a data frame with the following columns
#' \describe{
#' \item{t}{time}
#' \item{Xk}{columns X0, ..., X_K. which are numeric}
#' }
#' @details This function converts external data sources (we currently support
#'   output from the EpiModel and pomp R packages), which is already aggregated
#'   and puts it in a format that can be used by our exploring functions.
#' @export
#' @examples
#' ## For icm
#' out <- fortify_aggregate(EpiModel_icm)
#' head(out)
fortify_aggregate.icm <- function(data,
                              states = NULL,
                              package_source = NULL){

    fortify_aggregate.epimodel_inner(data,
                                     states = dplyr::enquo(states))

}


#' Get the viable names for the states from an icm object
#'
#' @param data of class (icm), output from the EpiModel package
#' @return  vector of names of format {x}.num where {x} is a SINGLE letter
get_epimodel_icm_states <- function(data){
    nms <- names(data$epi)
    out <- grep("^[A-z]{1}.num", x = nms,
                value = TRUE)
    return(out)

}


#' Extract the icm class totals
#'
#' @param nm name of variable to extract
#' @param ii integer between 0 and K - this is the state number
#' @param epi list of matrcies from icm output
#' @return a tibble
extract_icm_cols <- function(nm, ii, epi){

    mat <- epi[[nm]]
    nm <- paste0("X", ii)
    if(tidyr_new_interface()){
         df <- tidyr::pivot_longer(as.data.frame(mat),
                                   cols = tidyr::everything(),
                                   names_to = "sim",
                                   values_to = "X")
         names(df)[names(df) == "X"] <- nm
     } else{
         df <- tidyr::gather(as.data.frame(mat),
                             key = "sim")
         names(df)[names(df) == "value"] <- nm
     }
    return(df)

}

