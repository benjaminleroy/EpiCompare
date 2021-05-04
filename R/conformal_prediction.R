#' Interatively apply the mean-shift algorithm to filamental objects
#'
#' Wrapper of C++ function `psuedo_density_mode_cluster``
#'
#' @param X_list list of filamental objects (data frames or matrices) 
#' that define the psuedo-density estimation
#' @param G_list another list of filament objects (data frames or matrices) 
#' that we will "walk" up to their modes.
#' @param sigma scalar, scale value for the distances between observations
#' @param eps scalar, if difference between steps is less than this - treat as 
#' if the point has converged
#' @param maxT integer, max number of iterations of the algorithm
#' @param verbose boolean, if we should show progress (done with a  C++ 
#' based progressbar)
#' @param usefrac boolean, if the distance should be weighted by the number
#' of points in the filament.
#'
#' @return Final step of points in the \cite{G_list} (hopefully the associated
#' modes).
#' @export
functional_psuedo_density_mode_cluster <- function(X_list, G_list=X_list, 
                                                   sigma,
                                                   eps = 1e-05,  
                                                   maxT = 10,
                                                   verbose = TRUE,
                                                   usefrac = FALSE){
  # wrapper for psuedo_density_mode_cluster C++ function 
  # (some preliminary checks)
  # --
  X_list <- lapply(X_list, as.matrix)
  G_list <- lapply(X_list, as.matrix)
  
  distinct_size_X <- X_list %>% sapply(dim) %>% t() %>% data.frame() %>% 
    dplyr::distinct() 
  assertthat::assert_that(nrow(distinct_size_X) == 1, 
                          msg = "all matrices in X_list must be the same size.")
  distinct_size_G <- G_list %>% sapply(dim) %>% t() %>% data.frame() %>% 
    dplyr::distinct()
  assertthat::assert_that(nrow(distinct_size_G) == 1, 
                          msg = "all matrices in G_list must be the same size.")
  
  assertthat::assert_that(all(distinct_size_X == distinct_size_G),
                          msg = "all matrices in G_list and X_list must be the same size.")
  
  # running code
  out <- psuedo_density_mode_cluster(X_list = X_list,
                                     G_list = G_list,
                                     sigma = sigma,
                                     eps = eps,
                                     maxT = maxT,
                                     verbose = verbose,
                                     usefrac = usefrac)
  
  return(out)
}


#' Find mode clusters of filament objects 
#'
#' Wrapper of C++ function `psuedo_density_mode_cluster`` plus some cluster
#' finding
#'
#' @param g_list list of filamental objects (data frames or matrices) 
#' that we'd like to cluster (will be clustered based on the psuedo-density
#' estimate relative to themselves.)
#' @param g_names grouping information (as a data frame) asssociated with 
#' elements on the \code{g_list}
#' @param position vector of indices of columns that will be used to define the 
#' location of the points in space
#' @param sigma scalar, scale value for the distances between observations
#' @param eps scalar, if difference between steps is less than this - treat as 
#' if the point has converged
#' @param maxT integer, max number of iterations of the algorithm
#' @param diff_eps scalar, if the final step of each of the points is within
#' this distance from each-other they will be grouped together.
#' @param usefrac boolean, if we should use fraction approach to l2 distance
#' or not. Naturally is related to \code{sigma} value.
#' @param verbose  boolean, if we should show progress (done with a  C++ 
#' based progressbar)
#'
#' @return data frame with \code{g_names} and a \code{grouping} column which 
#' indicates which mode group each observation is in.
#' @export 
mode_clustering <- function(g_list, g_names, position, sigma, eps = 1e-05, maxT = 30,
                            diff_eps = 1e-05, 
                            usefrac = F,
                            verbose = T){
  g_list_inner <- g_list %>% 
    lapply(function(df) df[,position])
  
  out <- functional_psuedo_density_mode_cluster(X_list = g_list_inner,
                                                sigma = sigma,
                                                eps = eps,
                                                verbose = verbose,
                                                maxT = maxT,
                                                usefrac = usefrac
  )
  dist_mat <- dist_matrix_innersq_direction(out %>%
                                              lapply(function(x) as.data.frame(x)),
                                            position = 1:ncol(out[[1]]),
                                            usefrac = usefrac,
                                            verbose = verbose) # could updated to include usefrac...
  adjmatrix <- dist_mat <= diff_eps
  ig <- igraph::graph_from_adjacency_matrix(adjmatrix, mode = "undirected")

  groupings <- igraph::components(ig, mode = "strong")$membership
  
  if (inherits(g_names, "data.frame")){
    return(cbind(g_names, groupings))
  } else {
    return(data.frame(g_names, groupings))
  }
  
}



#' Calculates the minimum coverage radius for filamental objects as we increase 
#' the number of observations. 
#'
#' @param data_list list of objects ordered by correct ranking structure.
#' @param e_cols columns in data to compare with (in tidyverse fashion)
#' @param verbose boolean logic to be verbose when calculating
#'
#' @return list of min_cover_vec (minimum at that step) and dist_mat (cumulative 
#' max per path at time step t (column))
#' @export
#'
coverage_down_list <- function(data_list, 
                               e_cols, verbose = T){
  n_obs <- length(data_list)
  n_of_obs <- unique(sapply(data_list, nrow))
  
  if (length(n_of_obs) != 1){
    stop("data.frames in data_list are not the same length.")
  }
  
  
  dist_mat <- matrix(rep(NA, n_obs^2), ncol = n_obs)
  dist_mat[1,1] <- 0
  
  min_cover_vec <- rep(0, n_obs)
  
  if (n_obs > 1){
    if (verbose){
      pb <- progress::progress_bar$new(
        format = "calculating minimum coverage [:bar] :percent eta: :eta",
        total = n_obs, clear = FALSE, width = 60)
    }
    inner_min_dist_mat <- matrix(rep(Inf, n_obs*n_of_obs), ncol = n_of_obs)
    
    if (!pryr::is_promise(e_cols)){
      .e_cols_string <- inherits(e_cols, "character")
    } else {
      .e_cols_string <- F
    }
    
    if (.e_cols_string){
      for (r_idx in 2:n_obs){
        new_point_info <- inner_min_dist_per_point(data_list[[r_idx]], 
                                                   data_list[1:(r_idx-1)], e_cols,
                                                   .e_cols_string = .e_cols_string)
        inner_min_dist_mat[1:(r_idx-1),] <- array(
          c(inner_min_dist_mat[1:(r_idx-1),], new_point_info$nn_info_l),
          dim = c(nrow(new_point_info$nn_info_l),
                  ncol(new_point_info$nn_info_l),
                  2)) %>% 
          apply(1:2, min)
        inner_min_dist_mat[r_idx,] <- new_point_info$nn_info_current
        
        min_cover_vec[r_idx] <- max(inner_min_dist_mat[1:r_idx,]) # same as double max
        
        dist_mat[1:(r_idx-1), r_idx] <- sapply(dist_mat[1:(r_idx-1), r_idx-1],
                                               function(x) max(c(x, min_cover_vec[r_idx])))
        dist_mat[r_idx, r_idx] <- min_cover_vec[r_idx]
        
        
        if (verbose) {
          pb$tick()
        }
      }
    } else {
      for (r_idx in 2:n_obs){
        new_point_info <- inner_min_dist_per_point(data_list[[r_idx]], 
                                                   data_list[1:(r_idx-1)], {{e_cols}},
                                                   .e_cols_string = .e_cols_string)
        inner_min_dist_mat[1:(r_idx-1),] <- array(
          c(inner_min_dist_mat[1:(r_idx-1),], new_point_info$nn_info_l),
          dim = c(nrow(new_point_info$nn_info_l),
                  ncol(new_point_info$nn_info_l),
                  2)) %>% 
          apply(1:2, min)
        inner_min_dist_mat[r_idx,] <- new_point_info$nn_info_current
        
        min_cover_vec[r_idx] <- max(inner_min_dist_mat[1:r_idx,]) # same as double max
        
        dist_mat[1:(r_idx-1), r_idx] <- sapply(dist_mat[1:(r_idx-1), r_idx-1],
                                               function(x) max(c(x, min_cover_vec[r_idx])))
        dist_mat[r_idx, r_idx] <- min_cover_vec[r_idx]
        
        
        if (verbose) {
          pb$tick()
        }
      }
    }
  } else {
    message("A mode cluster has only a single element.")
  }
  
  return(list(min_cover_vec = min_cover_vec, dist_mat = dist_mat))
  
}


#' Calculates the minimum coverage radius for filamental objects as we increase 
#' the number of observations relative to mode groupings.
#'
#' @param data_ll list of objects (filaments)
#' @param e_cols columns in data to compare with (in tidyverse fashion)
#' @param g_order_ll list of indices, each vector is related to a single mode
#' cluster and the ordering defines how we build up the region.
#' @param names_df grouping information (as a data frame) asssociated with 
#' elements on the \code{g_list}
#' @param .td_out boolean, if outcome should a single "distance matrix" 
#' recombining different mode cluster information together. Else a list of 
#' "distance matrices" per mode will be returned
#' @param verbose boolean logic to be verbose when calculating
#'
#' @return matrix (or matrices) with the radius for each filament including at
#' the time step (column). See \code{.td_out} for more information.
#' @export
#'
coverage_down_mlist <- function(data_ll,
                                e_cols,
                                g_order_ll, # has which index of data_ll grouped together and the correct ordering...
                                names_df = NULL,
                                .td_out = F,
                                verbose = T){
  if (!(all(sapply(data_ll, function(x) inherits(x, "data.frame"))))){
    stop("data_ll must be a list of data.frames")
  }
  
  n_lists <- length(g_order_ll)
  n <- length(data_ll)
  
  assertthat::assert_that((sum(sapply(g_order_ll, length)) == n) &
                            (all(sort(unlist(g_order_ll)) == 1:n)),
                          msg = "data_ll length and length of g_order_ll info should match")
  
  
  
  
  if (verbose){
    pb <- progress::progress_bar$new(
      format = "calculating minimum coverage (per mode) [:bar] :percent eta: :eta",
      total = n_lists, clear = FALSE, width = 60)
  }
  coverd_info <- list()
  l_min_cover_vec <- list()
  l_dist_mat <- list()
  t <- 1
  

  
  for (order_vec in g_order_ll){
    if (pryr::is_promise(e_cols)){
      coverd_info[[t]] <- coverage_down_list(data_ll[order_vec], {{e_cols}},
                                             verbose = verbose)
    } else {
      coverd_info[[t]] <- coverage_down_list(data_ll[order_vec], e_cols,
                                             verbose = verbose)
    }
    l_min_cover_vec[[t]] <- coverd_info[[t]]$min_cover_vec
    l_dist_mat[[t]] <- coverd_info[[t]]$dist_mat
    
    if (verbose){
      pb$tick()
    }
    t <- t + 1
  }
  
  if (.td_out){
    # combining
    dist_mat_full <- matrix(NA, nrow = n, ncol = n)
    min_cover_vec_full <- rep(NA, n)
    for (g_id in 1:n_lists){
      ranking_inner <- g_order_ll[[g_id]]
      dist_mat_full[n+1-ranking_inner,n+1-ranking_inner] <- l_dist_mat[[g_id]]
      min_cover_vec_full[n+1-ranking_inner] <- l_min_cover_vec[[g_id]]
    }
    
    for (c_id in 2:n){
      row_loc <- ( (1:n %in% c(1:(c_id-1))) & is.na(dist_mat_full[,c_id]) )
      dist_mat_full[row_loc,c_id] <- dist_mat_full[row_loc,c_id-1]
    }
    
    if (is.null(names_df)){
      names_df <- data.frame(id = 1:n)
    }
    
    td <- tidy_dist_mat(dist_mat_full, 
                        rownames_df = names_df[unlist(g_order_ll),],
                        colnames_df = data.frame(step = 1:n))
    return(td)
  } else {
    return(coverd_info)
  }
  
}


#' calculate the minimum coverage needed for each point in a path relative to 
#' addition of the new path (and the mininum coverage for each point in new 
#' path relative to already seen paths)
#'
#' @param df new path
#' @param df_list list of old paths (ordered)
#' @param e_cols columns associated with location of path points
#' @param .e_cols_string boolean if e_cols is a string (not a promis)
#'
#' @return list of nn_info_l (update relative to single new path) and 
#' nn_info_current 
#' @export
inner_min_dist_per_point <- function(df, df_list, e_cols, .e_cols_string = F){
  nn_info_l <- matrix(NA, ncol = nrow(df), nrow = length(df_list))
  nn_info_r <- matrix(NA, ncol = nrow(df), nrow = length(df_list))
  
  if (!.e_cols_string){
    for (r_idx in 1:length(df_list)){
      nn_info_l[r_idx,] <- RANN::nn2(data = df %>% dplyr::select({{e_cols}}),
                                     query = df_list[[r_idx]] %>% 
                                       dplyr::select({{e_cols}}),
                                     k = 1, treetype = "kd")$nn.dists
      nn_info_r[r_idx,] <- RANN::nn2(query = df %>% dplyr::select({{e_cols}}),
                                     data = df_list[[r_idx]] %>% 
                                       dplyr::select({{e_cols}}),
                                     k = 1, treetype = "kd")$nn.dists
    }
  } else {
    for (r_idx in 1:length(df_list)){
      nn_info_l[r_idx,] <- RANN::nn2(data = df %>% 
                                       dplyr::select(dplyr::one_of(e_cols)),
                                     query = df_list[[r_idx]] %>% 
                                       dplyr::select(dplyr::one_of(e_cols)),
                                     k = 1, treetype = "kd")$nn.dists
      nn_info_r[r_idx,] <- RANN::nn2(query = df %>% 
                                       dplyr::select(dplyr::one_of(e_cols)),
                                     data = df_list[[r_idx]] %>% 
                                       dplyr::select(dplyr::one_of(e_cols)),
                                     k = 1, treetype = "kd")$nn.dists
    }
  }
  nn_info_rr <- apply(nn_info_r, 2, min)
  
  return(list(nn_info_l = nn_info_l, nn_info_current = nn_info_rr))
}


#' Combines psuedo-density estimate data frames and mode groupings and gets 
#' ordering per mode as well as globally
#'
#' Ranking ties are randomly assigned index values 
#'
#' @param psuedo_density_df data frame with object information and 
#' psuedo-density estimates (in column called \code{psuedo_density})
#' @param mode_grouping data frame with object information and 
#' mode clustering index (in column called \code{groupings})
#'
#' @return a list of 2 objects. The first is a data frame defined by the joining
#' of the previous data frames and additional columns for \code{ranking} of 
#' their psuedo_densities and for \code{group_ranking} (ranking conditional on 
#' mode group). The second element of the list is a lists of indices (associated
#' with the rows in the \code{psuedo_density_df}) these indices tell us the 
#' inner group order from highest psuedo-density to lowest.
#' @export
#' 
inner_expanding_info <- function(psuedo_density_df, mode_grouping){
  group_by_names <- names(mode_grouping)[names(mode_grouping) != "groupings"]
  both <- psuedo_density_df %>% 
    dplyr::left_join(mode_grouping, by = group_by_names) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ranking = rank(-.data$psuedo_density, 
                                 ties.method ="random")) %>%
    dplyr::group_by(.data$groupings) %>%
    dplyr::mutate(group_ranking = rank(-.data$psuedo_density, 
                                       ties.method ="random"))
  
  list_split <- both %>% tibble::rownames_to_column() %>% 
    dplyr::group_by(.data$groupings) %>%
    dplyr::group_split() %>% 
    lapply(function(df) as.numeric(df$rowname[order(df$group_ranking)]))
  
  return(list(both, list_split))
}


#' Helper fucntion for simulation-based conformal score calculation based 
#' potentially based on mode clustering and changing radius values.
#'
#' @param df_row_group data frame with filaments information (can be a 
#' grouped_df) to be assessed.
#' @param simulations_group_df grouped data frame with multiple simulated 
#' filaments information
#' @param data_column_names columns of \code{df_row_group} and 
#' \code{simulations_group_df} that correspond to the output space.
#' @param simulation_info_df a dataframe with information of each simulation's
#' psuedo-density estimate, mode clustering, ranking of psuedo-density(within 
#' cluster and overall). See \code{inner_expanding_info} for expected structure.
#' @param list_radius_info list of lists of radius information per each mode.
#' @param list_grouping_id list of vectors of indices (grouped by mode cluster)
#' and ordered by psuedo-density values.
#' @param verbose boolean, if progress should be tracked with a progress bar
#'
#' @return data.frame with a row for each filament in \code{df_row_group} with a 
#' section column \code{containment_val} which is the discrete simulation-based
#' conformal score.
#'
inner_containment_conformal_score_mode_radius <- function(df_row_group,
                                                          simulations_group_df, # need the list structure?
                                                          data_column_names,
                                                          simulation_info_df, # needs grouping info? group ranking and overall ranking? and psuedo_density
                                                          list_radius_info,
                                                          list_grouping_id,
                                                          verbose = F){
  # notes:
  # -
  # internally (until the very end) deals with the conformal score as a 
  # conformal score and  at the very end converts it to a non-conformal score.
  
  options(dplyr.summarise.inform = FALSE)
  all_info_list <- df_row_group %>% dplyr::group_split()
  all_info_df <- do.call(rbind, all_info_list) # not sure why we need this?
  all_info_nrows <- all_info_list %>% sapply(nrow)
  all_info_which <- lapply(1:length(all_info_nrows),
                           function(x) rep(x, all_info_nrows[x])) %>%
    unlist()
  all_info_names <- df_row_group %>% dplyr::group_keys()
  group_names <- names(all_info_names)
  all_info_names <- all_info_names %>%
    dplyr::mutate(single_numeric_idx = 1:nrow(all_info_names))
  
  
  group_names_sim <- simulations_group_df %>% dplyr::group_keys() %>% names()
  
  n_calibrate <-  nrow(all_info_df)
  
  
  all_points <- all_info_df[,data_column_names]
  
  n_draws <- nrow(simulation_info_df)
  n_groups <- length(list_grouping_id)
  
  if (verbose){
    pb <- progress::progress_bar$new(
      format = "Processing [:bar] :percent eta: :eta",
      total = n_draws,
      clear = FALSE, width = 60)
  }
  
  overall_info <- list()
  for (g_idx in 1:n_groups){
    inner_ids <- list_grouping_id[[g_idx]]
    
    nn_contained <- rep(F, n_calibrate) # for all points across calibration data
    nn_containment_num <- rep(n_draws+1, n_calibrate)
    nn_dist_mat <- matrix(0, nrow = length(inner_ids),
                          ncol = n_calibrate)
    
    inner_radius_mat <- list_radius_info[[g_idx]]$dist_mat # to get matrix...
    n_draw_mode <- length(inner_ids)
    current_idx <- 1
    for (sim_idx in inner_ids){
      
      pd_inner <- simulation_info_df[sim_idx,]
      df_col <- simulations_group_df %>% dplyr::inner_join(pd_inner,
                                                           by = group_names_sim)
      if (all(nn_contained)){
        break # then don't need to check if they're contained anymore
      }
      current_not_contained <- !nn_contained # not contained
      
      nn_info <- RANN::nn2(data = df_col[,data_column_names],
                           query = all_points[current_not_contained,],
                           k = 1, treetype = "kd")$nn.dists
      nn_dist_mat[current_idx,current_not_contained] <- nn_info 
      if (current_idx == 1){
        new_contained <- nn_dist_mat[1:current_idx,current_not_contained] < inner_radius_mat[1:current_idx,current_idx]
      } else {
        which_rows_news <- inner_radius_mat[1:(current_idx-1),(current_idx - 1)] < inner_radius_mat[1:(current_idx - 1),current_idx]
        which_rows_news_all <- (1:current_idx)[c(which_rows_news, T)] # if no change - don't look at...
        if (length(which_rows_news_all) == 1){
          new_contained <- nn_dist_mat[which_rows_news_all,current_not_contained] < inner_radius_mat[which_rows_news_all,current_idx]
        } else {
          if (sum(current_not_contained) == 1){
            new_contained <- sweep(x = nn_dist_mat[which_rows_news_all,current_not_contained, drop = F], 
                                   MARGIN = 1,
                                   STATS = inner_radius_mat[which_rows_news_all,current_idx], 
                                   FUN = "<") %>%
              colSums() %>% sapply(function(x) x > 0)
          } else {
            new_contained <- sweep(x = nn_dist_mat[which_rows_news_all,current_not_contained], 
                                   MARGIN = 1,
                                   STATS = inner_radius_mat[which_rows_news_all,current_idx], 
                                   FUN = "<") %>%
              colSums() %>% sapply(function(x) x > 0)
          }
        }
      }
      nn_containment_num[current_not_contained] <- new_contained * (current_idx) + (!new_contained) * (n_draws + 1)
      nn_contained[current_not_contained] <- new_contained
      
      if (verbose){
        pb$tick()
      }
      current_idx <- current_idx + 1
      
    }
    
    all_info_df <- all_info_df %>%
      dplyr::mutate(number_contained = nn_containment_num)
    
    # calculate score per mode
    score_per_mode <- all_info_df %>%
      dplyr::group_by(dplyr::across(tidyselect::one_of(group_names))) %>%
      dplyr::summarize(containment_val = max(.data$number_contained), .groups = "keep") 
    
    simulation_info_inner <- rbind((simulation_info_df %>% dplyr::ungroup() %>% 
      dplyr::select(.data$ranking, .data$group_ranking))[inner_ids,],
      data.frame(ranking = n_draws+1, group_ranking = n_draws+1))
    
    overall_score <- score_per_mode %>% 
      dplyr::left_join(simulation_info_inner,
                       by = c("containment_val" = "group_ranking")) 
    
    overall_info[[g_idx]] <- overall_score
    
    
  }
  # combining across modes
  moverall_info <- do.call(rbind,overall_info) %>%
    dplyr::group_by(dplyr::across(tidyselect::one_of(group_names))) %>% 
    dplyr::summarize(containment_val = n_draws+1 - min(.data$ranking), .groups = "keep") 
  
  return(moverall_info)
}



#' creates an increasing radius list structure for a fixed radius value.
#'
#' @param radius scalar, single radius value
#' @param group_structure list of lists with a vector that captures the number
#' of observations in the group.
#'
#' @return returns list of lists (per each group) with min_cover_vec (vector of
#' radius value) and the dist_mat (radius per each obs included per step (col))
inner_convert_single_radius_to_structure <- function(radius, group_structure){
  n_mode <- length(group_structure)
  inner_mode_size <- sapply(group_structure, length)
  
  out_list <- list()
  for (idx in 1:n_mode){
    inner_size <- inner_mode_size[idx]
    inner_dm <- matrix(NA, ncol = inner_size, nrow = inner_size)
    inner_dm[upper.tri(inner_dm, diag = T)] <- radius
    inner_list <- list(min_cover_vec = rep(radius, inner_size),
                       dist_mat = inner_dm)
    out_list[[idx]] <- inner_list
  }
  return(out_list)
}


#' global wrapper for simulation-based conformal score calculation that allows 
#' for mode clustering and changes of radius.
#'
#' @param truth_grouped_df data frame with filaments information (can be a 
#' grouped_df) to be assessed.
#' @param simulations_grouped_df grouped data frame with multiple simulated 
#' filaments information
#' @param data_column_names columns of \code{df_row_group} and 
#' \code{simulations_group_df} that correspond to the output space.
#' @param number_points  the number of points the filament should be compressed 
#' to (if \code{Inf}) then no compression is applied.
#' @param .change_radius boolean, if we should use the changing radius approach
#' @param .mode_cluster boolean, if we should use mode clustering
#' @param .to_simplex boolean, if we should project points on the unit simplex.
#' @param verbose boolean, be verbose about progression
#'
#' @return data.frame with a row for each filament in \code{truth_grouped_df} 
#' with a second column \code{containment_val} which is the discrete 
#' simulation-based conformal score.
#' @export
simulation_based_conformal3 <- function(truth_grouped_df, simulations_grouped_df,
                                        data_column_names = c("S", "I", "R"),
                                        number_points = 100,
                                        .change_radius = TRUE,
                                        .mode_cluster = TRUE,
                                        .to_simplex = TRUE,
                                        verbose = FALSE){
  
  # "record keeping" (keeping track of keys for sims and new obs)
  assertthat::assert_that(inherits(simulations_grouped_df, "grouped_df"))
  sim_group_names <- names(dplyr::group_keys(simulations_grouped_df))
  simulations_group_df_inner <- simulations_grouped_df[c(sim_group_names,
                                                         data_column_names)]
  
  assertthat::assert_that(inherits(truth_grouped_df, "grouped_df"))
  group_names_containment <- names(dplyr::group_keys(truth_grouped_df))
  truth_df_inner <- truth_grouped_df[c(group_names_containment,
                                       data_column_names)]
  
  sim_list <- simulations_group_df_inner %>% dplyr::group_split()
  sim_names <- simulations_group_df_inner %>% dplyr::group_keys()
  
  all_info_list <- truth_grouped_df %>% dplyr::group_split()
  all_info_df <- do.call(rbind,all_info_list)
  all_info_nrows <- all_info_list %>% sapply(nrow)
  all_info_which <- do.call(c,lapply(1:length(all_info_nrows),
                           function(x) rep(x, all_info_nrows[x])))
  all_info_names <- truth_grouped_df %>% dplyr::group_keys()
  g_names <- names(all_info_names)
  all_info_names <- all_info_names %>%
    dplyr::mutate(single_numeric_idx = 1:nrow(all_info_names))
  
  
  
  
  if (number_points == Inf){
    .do_filament_compression <- FALSE
  } else {
    .do_filament_compression <- TRUE
  }
  
  
  
  
  if (.to_simplex){
    truth_df_inner <- truth_df_inner %>% as.data.frame() %>%
      get_xy_coord(xyz_col = data_column_names) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_names_containment)))
    
    simulations_group_df_inner <- simulations_group_df_inner %>%
      as.data.frame() %>% get_xy_coord(xyz_col = data_column_names) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(g_names)))
    data_column_names <- c("x","y")
  }
  
  #truth_df_inner <- truth_df_inner %>%
  #  mutate_at(vars(one_of(group_names_containment)), as.character)
  
  truth_df_non_group_idx <- which(
    names(truth_df_inner) %in% data_column_names)
  
  #simulations_group_df_inner <- simulations_group_df_inner %>%
  #  mutate_at(vars(one_of(sim_group_names)), as.character)  # why is this done?
  
  simulations_group_df_non_group_idx <- which(
    names(simulations_group_df_inner) %in% data_column_names)
  
  if (.do_filament_compression){
    truth_df_inner <- truth_df_inner %>%
      filament_compression(data_columns = truth_df_non_group_idx,
                           number_points = number_points)
    
    simulations_group_df_inner <- simulations_group_df_inner %>%
      filament_compression(data_columns = simulations_group_df_non_group_idx,
                           number_points = number_points)
  }
  
  tdm_sims <- dist_matrix_innersq_direction(simulations_group_df_inner,
                                            position = simulations_group_df_non_group_idx,
                                            tdm_out = T)
  
  # sigma selection
  sigma_size <- c("20%" = .2, "25%" = .25, "30%" = .3,
                  "35%" = .35, "40%" = .4, "45%" = .45)
  
  percentage_str <- names(sigma_size)[stats::quantile(as.matrix(tdm_sims), sigma_size) > 0][1]
  
  percentage_inner <- check_character_percent(percentage_str, "sigma")
  sigma_val <- stats::quantile(as.matrix(tdm_sims), percentage_inner)
  
  # rank_df
  pseudo_density_df <- distance_psuedo_density_function(
    tdm_sims,
    sigma = sigma_val, df_out = T) %>%
    dplyr::mutate(ranking = rank(-.data$psuedo_density,ties.method = "random")) #spelling error... :(, no ties! need ordering
  
  assertthat::assert_that(all(!is.na(pseudo_density_df$psuedo_density)),
                          msg = paste("internal error in",
                                      "distance_psuedo_density_function",
                                      "function's sigma selection."))
  
  if (!.change_radius){
    top_points <- top_curves_to_points(simulations_group_df_inner,tidy_dm = tdm_sims,
                                       alpha = .2,
                                       quantile_func = distance_psuedo_density_function,
                                       sigma = sigma_val) # 80% curve remain
    
    assertthat::assert_that(nrow(top_points) > 0,
                            msg = paste("the number of simulations are so few",
                                        "thatthere are too few to estimate the",
                                        "shared radius."))
    
    mm_delta <- get_delta_nn(top_points[,data_column_names])
    
    # modes 
    if (.mode_cluster){
      out_groups <- mode_clustering(g_list = sim_list, g_names = sim_names,
                                    position = which(names(sim_list[[1]]) %in% data_column_names),
                                    sigma = sigma_val, 
                                    verbose = verbose)
      
      simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
    } else {
      out_groups <- sim_names %>% dplyr::mutate(groupings = 1)
      simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
    }
    
    tm_radius <- inner_convert_single_radius_to_structure(mm_delta, ordering_list)
    
    conformal_df <- inner_containment_conformal_score_mode_radius(
      df_row_group = truth_df_inner,
      simulations_group_df = simulations_group_df_inner, 
      data_column_names = data_column_names,
      simulation_info_df = simulation_info_df, 
      list_radius_info = tm_radius,
      list_grouping_id = ordering_list,
      verbose = verbose)
    
  } else {
    # varying radius
    mm_delta <- NA
    # get groupings...
    if (.mode_cluster){
      out_groups <- mode_clustering(g_list = sim_list, g_names = sim_names,
                                    position = which(names(sim_list[[1]]) %in% data_column_names),
                                    sigma = sigma_val, maxT = max(length(sim_list), 15),
                                    verbose = verbose)
      
      simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
    } else {
      out_groups <- sim_names %>% dplyr::mutate(groupings = 1)
      simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
    }
    
    
    tm_radius <- coverage_down_mlist(data_ll = sim_list,
                                     e_cols = data_column_names,
                                     g_order_ll =  ordering_list,
                                     names_df = sim_names,
                                     verbose = verbose)
    
    conformal_df <- inner_containment_conformal_score_mode_radius(
      df_row_group = truth_df_inner,
      simulations_group_df = simulations_group_df_inner, 
      data_column_names = data_column_names,
      simulation_info_df = simulation_info_df, 
      list_radius_info = tm_radius,
      list_grouping_id = ordering_list,
      verbose = verbose)
    
    
  }
  
  return(list(conformal_score = conformal_df, containment_df = simulation_info_df,
              mm_delta = mm_delta, tm_radius = tm_radius,
              truth_df_inner = truth_df_inner,
              simulations_group_df_inner = simulations_group_df_inner,
              parameters = c("mm_delta_prop" = .2,
                             "sigma_percentage" = percentage_str,
                             "filament_num_points" = number_points)))
}


#' global wrapper for simulation-based conformal score calculation that allows 
#' for mode clustering and changes of radius. 
#' 
#' Mode clustering also is able to be on more course approximate of the 
#' filaments.
#'
#' @param truth_grouped_df data frame with filaments information (can be a 
#' grouped_df) to be assessed.
#' @param simulations_grouped_df grouped data frame with multiple simulated 
#' filaments information
#' @param data_column_names columns of \code{df_row_group} and 
#' \code{simulations_group_df} that correspond to the output space.
#' @param number_points  the number of points the filament should be compressed 
#' to (if \code{Inf} then no compression is applied.
#' @param .change_radius boolean, if we should use the changing radius approach
#' @param .mode_cluster boolean, if we should use mode clustering
#' @param .to_simplex boolean, if we should project points on the unit simplex.
#' @param .use_frac boolean, if distances should be defined relative to a
#' scaling for the number of points in the path.
#' @param .small_size_mode_cluster int, the number of points the filament should 
#' be compressed to when estimating the mode clustering. We recommend a smaller
#' number to speed up the mode clustering algorithm.  (if Inf, then no change in 
#' size relative to \code{number_points})
#' @param .maxT int, max number of steps for mode clustering mean-shift 
#' algorithm
#' @param .sigma_string string, the quantile of the distance matrix to define
#' the sigma (e.g. \code{"30\%"})
#' @param .diff_eps float, the error allows between mode clustered final steps
#' to still be counted as the same mode.
#' @param verbose boolean, be verbose about progression
#'
#' @return list of information...
#' @export
simulation_based_conformal3.5 <- function(truth_grouped_df, simulations_grouped_df,
                                        data_column_names = c("S", "I", "R"),
                                        number_points = Inf,
                                        .change_radius = TRUE,
                                        .mode_cluster = TRUE,
                                        .to_simplex = TRUE,
                                        .use_frac = FALSE,
                                        .small_size_mode_cluster = Inf,
                                        .maxT = 50,
                                        .sigma_string = "35%",
                                        .diff_eps = 1e-06,
                                        verbose = FALSE){
  
  
  
  # "record keeping" (keeping track of keys for sims and new obs)
  assertthat::assert_that(inherits(simulations_grouped_df, "grouped_df"))
  sim_group_names <- names(dplyr::group_keys(simulations_grouped_df))
  simulations_group_df_inner <- simulations_grouped_df[c(sim_group_names,
                                                         data_column_names)]
  
  assertthat::assert_that(inherits(truth_grouped_df, "grouped_df"))
  group_names_containment <- names(dplyr::group_keys(truth_grouped_df))
  truth_df_inner <- truth_grouped_df[c(group_names_containment,
                                       data_column_names)]
  
  sim_list <- simulations_group_df_inner %>% dplyr::group_split()
  sim_names <- simulations_group_df_inner %>% dplyr::group_keys()
  
  all_info_list <- truth_grouped_df %>% dplyr::group_split()
  all_info_df <- do.call(rbind,all_info_list)
  all_info_nrows <- all_info_list %>% sapply(nrow)
  all_info_which <- do.call(c,lapply(1:length(all_info_nrows),
                                     function(x) rep(x, all_info_nrows[x])))
  all_info_names <- truth_grouped_df %>% dplyr::group_keys()
  g_names <- names(all_info_names)
  all_info_names <- all_info_names %>%
    dplyr::mutate(single_numeric_idx = 1:nrow(all_info_names))
  
  
  
  
  if (number_points == Inf){
    .do_filament_compression <- FALSE
  } else {
    .do_filament_compression <- TRUE
  }
  
  if (.small_size_mode_cluster == Inf & .mode_cluster){
    .do_mode_filament_compression <- FALSE
  } else {
    .do_mode_filament_compression <- TRUE
  }
  
  if(.mode_cluster & .do_mode_filament_compression & !.use_frac){
    message(paste("if you are doing filament compression for mode estimation",
                   "we recommend setting '.use_frac' to TRUE."))
  }
  
  
  if (.to_simplex){
    truth_df_inner <- truth_df_inner %>% as.data.frame() %>%
      get_xy_coord(xyz_col = data_column_names) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_names_containment)))
    
    simulations_group_df_inner <- simulations_group_df_inner %>%
      as.data.frame() %>% get_xy_coord(xyz_col = data_column_names) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(g_names)))
    data_column_names <- c("x","y")
  }
  
  #truth_df_inner <- truth_df_inner %>%
  #  mutate_at(vars(one_of(group_names_containment)), as.character)
  
  truth_df_non_group_idx <- which(
    names(truth_df_inner) %in% data_column_names)
  
  #simulations_group_df_inner <- simulations_group_df_inner %>%
  #  mutate_at(vars(one_of(sim_group_names)), as.character)  # why is this done?
  
  simulations_group_df_non_group_idx <- which(
    names(simulations_group_df_inner) %in% data_column_names)
  
  if (.do_filament_compression){
    truth_df_inner <- truth_df_inner %>%
      filament_compression(data_columns = truth_df_non_group_idx,
                           number_points = number_points)
    
    simulations_group_df_inner <- simulations_group_df_inner %>%
      filament_compression(data_columns = simulations_group_df_non_group_idx,
                           number_points = number_points)
  }
  
  tdm_sims <- dist_matrix_innersq_direction(simulations_group_df_inner,
                                            position = simulations_group_df_non_group_idx,
                                            usefrac = .use_frac,
                                            tdm_out = T,
                                            verbose = verbose)
  
  # sigma selection
  sigma_lower <- check_character_percent(.sigma_string, ".sigma_string")
  sigma_sizes <- sapply(sigma_lower + .05*(0:5), function(v) min(v, 1))
  
  percentage_inner <- sigma_sizes[stats::quantile(as.matrix(tdm_sims), sigma_sizes) > 0][1]
  
  sigma_val <- stats::quantile(as.matrix(tdm_sims), percentage_inner)
  
  # rank_df
  pseudo_density_df <- distance_psuedo_density_function(
    tdm_sims,
    sigma = sigma_val, df_out = T) %>%
    dplyr::mutate(ranking = rank(-.data$psuedo_density,ties.method = "random")) #spelling error... :(, no ties! need ordering
  
  assertthat::assert_that(all(!is.na(pseudo_density_df$psuedo_density)),
                          msg = paste("internal error in",
                                      "distance_psuedo_density_function",
                                      "function's sigma selection."))
  
  if (!.change_radius){
    top_points <- top_curves_to_points(simulations_group_df_inner,tidy_dm = tdm_sims,
                                       alpha = .2,
                                       quantile_func = distance_psuedo_density_function,
                                       sigma = sigma_val) # 80% curve remain
    
    assertthat::assert_that(nrow(top_points) > 0,
                            msg = paste("the number of simulations are so few",
                                        "thatthere are too few to estimate the",
                                        "shared radius."))
    
    mm_delta <- get_delta_nn(top_points[,data_column_names])
    
    # modes 
    if (.mode_cluster){
      if(!.do_mode_filament_compression){
        out_groups <- mode_clustering(g_list = sim_list, g_names = sim_names,
                                      position = which(names(sim_list[[1]]) %in% data_column_names),
                                      sigma = sigma_val, maxT = .maxT,
                                      usefrac = .use_frac, diff_eps = .diff_eps,
                                      verbose = verbose)
        
        simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
        simulation_info_df <- simulation_info_out[[1]]
        ordering_list <- simulation_info_out[[2]]
      } else {
        #browser()
        c_s_out <- compression_and_sigma_estimate(sim_grouped_df = simulations_group_df_inner,
                                                  data_columns = data_column_names,
                                                  usefrac = .use_frac,
                                                  number_points = .small_size_mode_cluster,
                                                  .sigma_string = .sigma_string,
                                                  verbose = verbose)
        c_g_list <- c_s_out$compression %>% dplyr::group_split()
        c_g_names <-  c_s_out$compression %>% dplyr::group_keys()
        c_position <- (1:ncol(c_s_out$compression))[names(c_s_out$compression) %in% data_column_names]
        out_groups <- mode_clustering(g_list = c_g_list,g_names = c_g_names,
                                   position = c_position,
                                   sigma = c_s_out$sigma, maxT = .maxT, 
                                   usefrac = .use_frac,diff_eps = .diff_eps,
                                   verbose = verbose)
        
        simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
        simulation_info_df <- simulation_info_out[[1]]
        ordering_list <- simulation_info_out[[2]]
      }
    } else {
      out_groups <- sim_names %>% dplyr::mutate(groupings = 1)
      simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
    }
    
    tm_radius <- inner_convert_single_radius_to_structure(mm_delta, ordering_list)
    
    conformal_df <- inner_containment_conformal_score_mode_radius(
      df_row_group = truth_df_inner,
      simulations_group_df = simulations_group_df_inner, 
      data_column_names = data_column_names,
      simulation_info_df = simulation_info_df, 
      list_radius_info = tm_radius,
      list_grouping_id = ordering_list,
      verbose = verbose)
    
  } else {
    # varying radius
    mm_delta <- NA
    # get groupings...
    if (.mode_cluster){
      if(!.do_mode_filament_compression){
        out_groups <- mode_clustering(g_list = sim_list, g_names = sim_names,
                                      position = which(names(sim_list[[1]]) %in% data_column_names),
                                      sigma = sigma_val, maxT = .maxT,
                                      usefrac = .use_frac, diff_eps =.diff_eps,
                                      verbose = verbose)
        
        simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
        simulation_info_df <- simulation_info_out[[1]]
        ordering_list <- simulation_info_out[[2]]
      } else {
        #browser()
        c_s_out <- compression_and_sigma_estimate(sim_grouped_df = simulations_group_df_inner,
                                                  data_columns = data_column_names,
                                                  usefrac = .use_frac,
                                                  number_points = .small_size_mode_cluster,
                                                  .sigma_string = .sigma_string,
                                                  verbose = verbose)
        c_g_list <- c_s_out$compression %>% dplyr::group_split()
        c_g_names <-  c_s_out$compression %>% dplyr::group_keys()
        c_position <- (1:ncol(c_s_out$compression))[names(c_s_out$compression) %in% data_column_names]
        out_groups <- mode_clustering(g_list = c_g_list,g_names = c_g_names,
                                      position = c_position,
                                      sigma = c_s_out$sigma, maxT = .maxT, 
                                      usefrac = .use_frac,diff_eps = .diff_eps,
                                      verbose = verbose)
        
        simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
        simulation_info_df <- simulation_info_out[[1]]
        ordering_list <- simulation_info_out[[2]]
      }
    } else {
      out_groups <- sim_names %>% dplyr::mutate(groupings = 1)
      simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
    }
    
    
    tm_radius <- coverage_down_mlist(data_ll = sim_list,
                                     e_cols = data_column_names,
                                     g_order_ll =  ordering_list,
                                     names_df = sim_names,
                                     verbose = verbose)
    
    conformal_df <- inner_containment_conformal_score_mode_radius(
      df_row_group = truth_df_inner,
      simulations_group_df = simulations_group_df_inner, 
      data_column_names = data_column_names,
      simulation_info_df = simulation_info_df, 
      list_radius_info = tm_radius,
      list_grouping_id = ordering_list,
      verbose = verbose)
    
    
  }
  
  return(list(conformal_score = conformal_df, containment_df = simulation_info_df,
              mm_delta = mm_delta, tm_radius = tm_radius,
              truth_df_inner = truth_df_inner,
              simulations_group_df_inner = simulations_group_df_inner,
              parameters = c("mm_delta_prop" = .2,
                             "sigma_percentage" = percentage_inner,
                             "filament_num_points" = number_points)))
}

#' global wrapper for simulation-based conformal score calculation that allows 
#' for mode clustering and changes of radius. 
#' 
#' Does all the 4 types of approaches (mode vs non mode clustering) and (
#' change in radius vs fixed radius)
#'
#' @param truth_grouped_df data frame with filaments information (can be a 
#' grouped_df) to be assessed.
#' @param simulations_grouped_df grouped data frame with multiple simulated 
#' filaments information
#' @param data_column_names columns of \code{df_row_group} and 
#' \code{simulations_group_df} that correspond to the output space.
#' @param number_points  the number of points the filament should be compressed 
#' to (if \code{Inf} then no compression is applied.
#' @param .to_simplex boolean, if we should project points on the unit simplex.
#' @param .use_frac boolean, if distances should be defined relative to a
#' scaling for the number of points in the path.
#' @param .small_size_mode_cluster int, the number of points the filament should 
#' be compressed to when estimating the mode clustering. We recommend a smaller
#' number to speed up the mode clustering algorithm.  (if Inf, then no change in 
#' size relative to \code{number_points})
#' @param .maxT int, max number of steps for mode clustering mean-shift 
#' algorithm
#' @param .sigma_string string, the quantile of the distance matrix to define
#' the sigma (e.g. \code{"30\%"})
#' @param .diff_eps float, the error allows between mode clustered final steps
#' to still be counted as the same mode.
#' @param verbose boolean, be verbose about progression
#' @param return_min boolean, if list of information returned is mimimum (for slurm)
#'
#' @return list of information...
#' @export
simulation_based_conformal4 <- function(truth_grouped_df, simulations_grouped_df,
                                        data_column_names = c("S", "I", "R"),
                                        number_points = Inf,
                                        .to_simplex = FALSE,
                                        .use_frac = TRUE,
                                        .small_size_mode_cluster = 50,
                                        .maxT = 50,
                                        .sigma_string = "35%",
                                        .diff_eps = 1e-06,
                                        verbose = FALSE,
                                        return_min = FALSE){
  # "record keeping" (keeping track of keys for sims and new obs)
  assertthat::assert_that(inherits(simulations_grouped_df, "grouped_df"))
  sim_group_names <- names(dplyr::group_keys(simulations_grouped_df))
  simulations_group_df_inner <- simulations_grouped_df[c(sim_group_names,
                                                         data_column_names)]
  
  assertthat::assert_that(inherits(truth_grouped_df, "grouped_df"))
  group_names_containment <- names(dplyr::group_keys(truth_grouped_df))
  truth_df_inner <- truth_grouped_df[c(group_names_containment,
                                       data_column_names)]
  
  sim_list <- simulations_group_df_inner %>% dplyr::group_split()
  sim_names <- simulations_group_df_inner %>% dplyr::group_keys()
  
  all_info_list <- truth_grouped_df %>% dplyr::group_split()
  all_info_df <- do.call(rbind,all_info_list)
  all_info_nrows <- all_info_list %>% sapply(nrow)
  all_info_which <- do.call(c,lapply(1:length(all_info_nrows),
                                     function(x) rep(x, all_info_nrows[x])))
  all_info_names <- truth_grouped_df %>% dplyr::group_keys()
  g_names <- names(all_info_names)
  all_info_names <- all_info_names %>%
    dplyr::mutate(single_numeric_idx = 1:nrow(all_info_names))
  

  if (number_points == Inf){
    .do_filament_compression <- FALSE
  } else {
    .do_filament_compression <- TRUE
  }

  
  if (.to_simplex){
    truth_df_inner <- truth_df_inner %>% as.data.frame() %>%
      get_xy_coord(xyz_col = data_column_names) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_names_containment)))
    
    simulations_group_df_inner <- simulations_group_df_inner %>%
      as.data.frame() %>% get_xy_coord(xyz_col = data_column_names) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(g_names)))
    data_column_names <- c("x","y")
  }
  
  #truth_df_inner <- truth_df_inner %>%
  #  mutate_at(vars(one_of(group_names_containment)), as.character)
  
  truth_df_non_group_idx <- which(
    names(truth_df_inner) %in% data_column_names)
  
  #simulations_group_df_inner <- simulations_group_df_inner %>%
  #  mutate_at(vars(one_of(sim_group_names)), as.character)  # why is this done?
  
  simulations_group_df_non_group_idx <- which(
    names(simulations_group_df_inner) %in% data_column_names)
  
  if (.do_filament_compression){
    truth_df_inner <- truth_df_inner %>%
      filament_compression(data_columns = truth_df_non_group_idx,
                           number_points = number_points)
    
    simulations_group_df_inner <- simulations_group_df_inner %>%
      filament_compression(data_columns = simulations_group_df_non_group_idx,
                           number_points = number_points)
  }
  
  tdm_sims <- dist_matrix_innersq_direction(simulations_group_df_inner,
                                            position = simulations_group_df_non_group_idx,
                                            usefrac = .use_frac,
                                            tdm_out = T,
                                            verbose = verbose)
  
  # sigma selection
  sigma_lower <- check_character_percent(.sigma_string, ".sigma_string")
  sigma_sizes <- sapply(sigma_lower + .05*(0:5), function(v) min(v, 1))
  
  percentage_inner <- sigma_sizes[stats::quantile(as.matrix(tdm_sims), sigma_sizes) > 0][1]
  
  sigma_val <- stats::quantile(as.matrix(tdm_sims), percentage_inner)
  
  # rank_df
  pseudo_density_df <- distance_psuedo_density_function(
    tdm_sims,
    sigma = sigma_val, df_out = T) %>%
    dplyr::mutate(ranking = rank(-.data$psuedo_density,ties.method = "random")) #spelling error... :(, no ties! need ordering
  
  assertthat::assert_that(all(!is.na(pseudo_density_df$psuedo_density)),
                          msg = paste("internal error in",
                                      "distance_psuedo_density_function",
                                      "function's sigma selection."))
  
  
  if (.small_size_mode_cluster == Inf){
    .do_mode_filament_compression <- FALSE
  } else {
    .do_mode_filament_compression <- TRUE
  }
  
  top_points <- top_curves_to_points(simulations_group_df_inner,tidy_dm = tdm_sims,
                                     alpha = .2,
                                     quantile_func = distance_psuedo_density_function,
                                     sigma = sigma_val) # 80% curve remain
  
  assertthat::assert_that(nrow(top_points) > 0,
                          msg = paste("the number of simulations are so few",
                                      "thatthere are too few to estimate the",
                                      "shared radius."))
  
  mm_delta <- get_delta_nn(top_points[,data_column_names])
    
  # mode clustering
  if(!.do_mode_filament_compression){
      out_groups <- mode_clustering(g_list = sim_list, g_names = sim_names,
                                    position = which(names(sim_list[[1]]) %in% data_column_names),
                                    sigma = sigma_val, maxT = .maxT,
                                    usefrac = .use_frac, diff_eps = .diff_eps,
                                    verbose = verbose)
      
      simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
    } else {
      #browser()
      c_s_out <- compression_and_sigma_estimate(sim_grouped_df = simulations_group_df_inner,
                                                data_columns = data_column_names,
                                                usefrac = .use_frac,
                                                number_points = .small_size_mode_cluster,
                                                .sigma_string = .sigma_string,
                                                verbose = verbose)
      c_g_list <- c_s_out$compression %>% dplyr::group_split()
      c_g_names <-  c_s_out$compression %>% dplyr::group_keys()
      c_position <- (1:ncol(c_s_out$compression))[names(c_s_out$compression) %in% data_column_names]
      out_groups <- mode_clustering(g_list = c_g_list,g_names = c_g_names,
                                    position = c_position,
                                    sigma = c_s_out$sigma, maxT = .maxT, 
                                    usefrac = .use_frac,diff_eps = .diff_eps,
                                    verbose = verbose)
      
      simulation_info_out <- inner_expanding_info(pseudo_density_df, out_groups)
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
    }
  
  
  tm_radius_fixed <- inner_convert_single_radius_to_structure(mm_delta, 
                                                              ordering_list)
  
  # no mode clustering
  out_groups_nm <- sim_names %>% dplyr::mutate(groupings = 1)
  simulation_info_out_nm <- inner_expanding_info(pseudo_density_df, out_groups_nm)
  simulation_info_df_nm <- simulation_info_out_nm[[1]]
  ordering_list_nm <- simulation_info_out_nm[[2]]

    
  tm_radius_fixed_nm <- inner_convert_single_radius_to_structure(mm_delta, 
                                                                 ordering_list_nm)
  
    
  tm_radius_vary <- coverage_down_mlist(data_ll = sim_list,
                                       e_cols = data_column_names,
                                       g_order_ll =  ordering_list,
                                       names_df = sim_names,
                                       verbose = verbose)
  
  tm_radius_vary_nm <- coverage_down_mlist(data_ll = sim_list,
                                           e_cols = data_column_names,
                                           g_order_ll =  ordering_list_nm,
                                           names_df = sim_names,
                                           verbose = verbose)
    
  conformal_df_fixed <- inner_containment_conformal_score_mode_radius(
      df_row_group = truth_df_inner,
      simulations_group_df = simulations_group_df_inner, 
      data_column_names = data_column_names,
      simulation_info_df = simulation_info_df, 
      list_radius_info = tm_radius_fixed, # diff
      list_grouping_id = ordering_list, # diff
      verbose = verbose)
    
  conformal_df_fixed_nm <- inner_containment_conformal_score_mode_radius(
    df_row_group = truth_df_inner,
    simulations_group_df = simulations_group_df_inner, 
    data_column_names = data_column_names,
    simulation_info_df = simulation_info_df, 
    list_radius_info = tm_radius_fixed_nm, # diff
    list_grouping_id = ordering_list_nm, # diff
    verbose = verbose)
  
  conformal_df_vary <- inner_containment_conformal_score_mode_radius(
    df_row_group = truth_df_inner,
    simulations_group_df = simulations_group_df_inner, 
    data_column_names = data_column_names,
    simulation_info_df = simulation_info_df, 
    list_radius_info = tm_radius_vary, # diff
    list_grouping_id = ordering_list, # diff
    verbose = verbose)
  
  conformal_df_vary_nm <- inner_containment_conformal_score_mode_radius(
    df_row_group = truth_df_inner,
    simulations_group_df = simulations_group_df_inner, 
    data_column_names = data_column_names,
    simulation_info_df = simulation_info_df, 
    list_radius_info = tm_radius_vary_nm, # diff
    list_grouping_id = ordering_list_nm, # diff
    verbose = verbose)
    
  if (return_min){
    return(list(conformal_score = list(fixed = conformal_df_fixed,
                                       fixed_nm = conformal_df_fixed_nm,
                                       vary = conformal_df_vary,
                                       vary_nm = conformal_df_vary_nm), 
                parameters = c("mm_delta_prop" = .2,
                               "mm_delta" = mm_delta,
                               "sigma_percentage" = percentage_inner,
                               "filament_num_points" = number_points)))
  } else{
    return(list(conformal_score = list(fixed = conformal_df_fixed,
                                       fixed_nm = conformal_df_fixed_nm,
                                       vary = conformal_df_vary,
                                       vary_nm = conformal_df_vary_nm), 
                containment_df = simulation_info_df,
                mm_delta = mm_delta, 
                tm_radius = list(fixed = tm_radius_fixed,
                                 fixed_nm = tm_radius_fixed_nm,
                                 vary = tm_radius_vary,
                                 vary_nm = tm_radius_vary_nm),
                truth_df_inner = truth_df_inner,
                simulations_group_df_inner = simulations_group_df_inner,
                parameters = c("mm_delta_prop" = .2,
                               "sigma_percentage" = percentage_inner,
                               "filament_num_points" = number_points)))
  }

  
}


#' calculate maxmin distance between points
#'
#' Fast calculation of maxmin distance using kd trees and nearest neighbors from
#' the \code{RANN} package.
#' 
#' Taken from \code{simulationBands} package (written by Ben)
#'
#' @param df data.frame (with only columns that are needed)
#'
#' @return minimum radius for all points to be covered
#' @export
get_delta_nn <- function(df){
  check <- RANN::nn2(df, df, k = 2, treetype = "kd", eps = 0)
  mm_delta <- check$nn.dists[,2] %>% max()
  return(mm_delta)
}



#' inner l2 distance between filaments
#'
#' Basically the same as \code{l2filamentdist_df}
#'
#' @param df1 data.frame (n x p)
#' @param df2 data.frame (n x p)
#' @param usefrac if we should calculate the distance relative to a scaling of 
#' 1/nrow(\code{df1}) (before taking the square-root).
#'
#' @return distance between filaments 
#' @export
l2_filament_distance <- function(df1, df2, usefrac = F){
  difffunction(as.matrix(df1), as.matrix(df2), usefrac = usefrac)
}


#' compress filaments and estimate psuedo-density sigma in 1 go
#'
#' @param sim_grouped_df grouped_df, set of simulated filaments
#' @param data_columns string of data columns that are associated with location
#' information
#' @param usefrac boolean, if the distance between filaments should scale 
#' relative to the number of points in compression (before sqrt)
#' @param number_points number of points to compressed the data with
#' @param .sigma_string string with percentage use for sigma (relative to 
#' distance matrix), e.g. \code{"30\%"}
#' @param verbose boolean, if we should be verbose when processing...
#'
#' @return list of compression of data as well as the sigma estimate...
#' @export
compression_and_sigma_estimate <- function(sim_grouped_df,
                                           data_columns,
                                           usefrac = T,
                                           number_points = 50,
                                           .sigma_string = "30%",
                                           verbose = T){
  # compression of points ----
  compressed_points <- sim_grouped_df %>%
    filament_compression(data_columns = data_columns,
                         number_points = number_points)
  
  # estimate of sigma ----
  position <- c(1:ncol(compressed_points))[
    names(compressed_points) %in% data_columns]
  
  tdm_sims <- dist_matrix_innersq_direction(compressed_points,
                                            position = position,
                                            usefrac = usefrac,
                                            verbose = verbose,
                                            tdm_out = T)
  # sigma selection
  sigma_lower <- check_character_percent(.sigma_string, ".sigma_string")
  sigma_sizes <- sapply(sigma_lower + .05*(0:5), function(v) min(v, 1))
  
  percentage_inner <- sigma_sizes[stats::quantile(as.matrix(tdm_sims), sigma_sizes) > 0][1]
  
  sigma_val <- stats::quantile(as.matrix(tdm_sims), percentage_inner)
  
  return(list(compression = compressed_points,
              sigma = sigma_val))
  
}

#------------------------------------------------
# functions to compare to and use to test C++ versions 
#------------------------------------------------





#' An R (no cpp) version of `functional_psuedo_density_mode_cluster`
#' 
#' This is not exported and is only included to a R sleuth to see what the 
#' function (that is written in C++ looks like in R). This won't be updated
#' even if the C++ code is (as of April 7, 2021).
#' 
#' This function is also used to test out the C++ version of the function
#'
#' @param X_list data list
#' @param G_list points to move 
#' @param sigma sigma value
#' @param eps error term
#' @param maxT max number of steps
#' @param verbose boolean if we should be verbose
#'
#' @return list of updated filaments for points in G_list
functional_psuedo_density_mode_cluster_r <- function(X_list,
                                                     G_list=X_list,
                                                     sigma, # parameters of distance_kernel_function
                                                     eps = 1e-05,
                                                     maxT = 10,
                                                     verbose = TRUE){
  # usefrac is always false (in distance_kernel_function_r and diff_function_r)
  
  t <- 0
  n <- length(X_list)
  m <- length(G_list)
  error <- rep(1e+08, times = m) # initial error = massive error
  
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "taking time steps [:bar] :percent eta: :eta",
      total = maxT, clear = FALSE, width = 60)
  }
  
  while ((max(error) > eps) & (t < maxT)) {
    t = t + 1
    for (j in 1:m){
      
      if (error[j] > eps){
        c_val <- distance_kernel_function_r(G_list[[j]],X_list,sigma = sigma)
        denominator <- sum(c_val)
        
        unscaled_new <- c_val[1]*X_list[[1]]
        
        for (i in 2:n){
          unscaled_new <- unscaled_new + c_val[i]*X_list[[i]]
        }
        tmp <- G_list[[j]]
        G_list[[j]] <- 1/denominator * unscaled_new 
        error[j] <- sqrt(diff_function_r(tmp,G_list[[j]])^2)
      }
      
    }
    if (verbose){
      pb$tick()
    }
  }
  return(G_list)
}




#' l2 distance between filaments
#'
#' Also not exported, similar to `difffunction `in c++
#'
#' @param df1 data.frame / matrix
#' @param df2 data.frame / matrix (same size as m2)
#' @param usefrac bool if we should calculate the distance relative to a scaling
#' of 1/nrow(\code{df1}) (before taking the square root).
#'
#' @return numerical distance between \code{df1} and \code{df2}
diff_function_r <- function(df1, df2, usefrac = F){
  if (usefrac){
    frac <- 1/nrow(df1)
  } else {
    frac <- 1
  }
  
  sqrt(frac*sum((df1-df2)^2))
}


#' guassian kernel calculation relative to distances and scaling 
#'
#' This function is an r combination of `inner_dist_density` and `distvec` 
#'
#' @param df filament to calculate psuedo-density estimate for
#' @param df_list list of filaments that define the psuedo-density
#' @param sigma scalar, the scaling relative to distnaces
#'
#' @return numerical value of psuedo-density estimate for \code{df}.
distance_kernel_function_r <- function(df, df_list, sigma){
  assertthat::assert_that(sigma > 0, msg = "sigma should be greater than 0")
  
  dist_vec <- sapply(df_list, function(df2) diff_function_r(df,df2)) # distvec
  psuedo_density <- stats::dnorm(dist_vec/sigma) # inner_dist_density
  return(psuedo_density)
}
