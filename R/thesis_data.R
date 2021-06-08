
#' generate a set of simulations from the basic SIR hierachical example based on
#' a specific x value.
#' 
#' 
#'
#' @param x_inner scalar x value (between 0 -1)
#' @param n_sims_containment number of simulations to create
#' @param number_points number of points in the filamental compression
#' @param verbose boolean, if should report progress
#'
#' @return data frame with set of simulations after filamental compression
#' @export
data_generation <- function(x_inner, n_sims_containment = 1000, 
                            number_points = 150,
                            verbose = T){
  
  inner_truth_df <- data.frame(x = rep(x_inner, n_sims_containment))
  inner_truth_df$R0 <- sapply(inner_truth_df$x, function(x){
    1/7*simulationBands::lei_wasserman_data_conditional_simulate(2*(x-.5),
                                                                 n =1)[[1]]$sim+2})
  inner_truth_df <- inner_truth_df %>%
    dplyr::mutate(beta = .1,
           gamma = .1/R0,
           idx = paste0("inner_truth", 1:n_sims_containment))
  
  truth_paths <- data.frame(x_original = as.numeric(rep(NA, number_points*n_sims_containment)),
                            R0 = as.numeric(rep(NA, number_points*n_sims_containment)),
                            beta = as.numeric(rep(NA, number_points*n_sims_containment)),
                            gamma = as.numeric(rep(NA, number_points*n_sims_containment)),
                            idx = as.character(rep(NA, number_points*n_sims_containment)),
                            set_idx = as.character(rep(NA, number_points*n_sims_containment)),
                            x = as.numeric(rep(NA, number_points*n_sims_containment)),
                            y = as.numeric(rep(NA, number_points*n_sims_containment))) %>%
    dplyr::mutate(x_original = as.numeric(x_original),
                  R0 = as.numeric(R0),
                  beta = as.numeric(beta),
                  gamma = as.numeric(gamma),
                  idx = as.character(r_idx),
                  set_idx = as.character(set_idx),
                  x = as.numeric(x),
                  y = as.numeric(y))
  
  
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "creating data [:bar] :percent eta: :eta",
      total = n_sims_containment, clear = FALSE, width = 38)
  }
  for (r_idx in 1:n_sims_containment){
    inner_info <- inner_truth_df[r_idx, ]
    global_r_idx <- ((r_idx-1)* number_points + 1):(r_idx* number_points)
    
    inner_paths <- simulate_SIR_agents(n_sims = 1,
                                       n_time_steps = 500,
                                       inner_info$beta, inner_info$gamma,
                                       init_SIR = c(950,50,0)) %>%
      agents_to_aggregate(states = c("tI", "tR")) %>%
      as.data.frame() %>%
      get_xy_coord(xyz_col = c("X0","X1","X2")) %>%
      dplyr::mutate(idx = 1) %>%
      dplyr::group_by(idx) %>%
      filament_compression(data_columns = c("x","y"),
                           number_points = number_points) %>%
      dplyr::ungroup() %>% dplyr::select(-idx)
    
    
    inner_df <- cbind(x_original = inner_info$x,
                      R0 = inner_info$R0,
                      beta = inner_info$beta,
                      gamma = inner_info$gamma,
                      idx = r_idx,
                      set_idx = inner_info$idx,
                      inner_paths) %>% 
      dplyr::mutate(x_original = as.numeric(x_original),
                    R0 = as.numeric(R0),
                    beta = as.numeric(beta),
                    gamma = as.numeric(gamma),
                    idx = as.character(r_idx),
                    set_idx = as.character(set_idx))
    
    truth_paths[global_r_idx, ] <- inner_df
    if (verbose){
      pb$tick()
    }
  }
  
  truth_paths <- truth_paths %>% 
    dplyr::group_by(x_original, R0, beta, gamma, idx, set_idx)
  
  
  return(truth_paths)
}

#' Interior function for mode clustering across a range of eps given a single
#' maxT value
#'
#' @param x_list_inner list of curves to define psuedo-density 
#' @param g_list_inner list of curves to walk up the psuedo-density
#' @param sigma scalar for the sigma of the psuedo-density
#' @param range_eps vector and eps value (descending in size)
#' @param maxT int, max number of iterations
#'
#' @return linked list for each eps explored before maxT = 0.
interior_clustering <- function(x_list_inner, g_list_inner,
                                sigma,
                                range_eps,
                                maxT){
  
  eps <- range_eps[1]
  
  out <- functional_psuedo_density_mode_cluster2(X_list = x_list_inner,
                                                 G_list = g_list_inner,
                                                 sigma = sigma,
                                                 eps = eps,
                                                 verbose = F,
                                                 maxT = maxT,
                                                 usefrac = T
  )
  if (!out[[1]] & (length(range_eps) > 1)) {
    # then step to next eps...
    num_new_steps <- maxT - out[[2]]
    g_list_inner <- out[[3]]
    
    # check analysis for current eps (this is the final stop)
    # move to next eps
    out2 <- interior_clustering(x_list_inner, g_list_inner,
                                sigma,
                                range_eps[2:length(range_eps)],
                                num_new_steps)
    out_list <- list(out, next_element = out2) # linked list
    names(out_list)[1] <- as.character(eps)
    return(out_list)
  } else {
    out_list <- list(out, next_element = F)  # linked list
    names(out_list)[1] <- as.character(eps)
    
  }
  return(out_list)
}

#' Interior function to check number of groups from a distance matrix and cutoff
#'
#' @param dist_mat distance matrix
#' @param diff_eps threshold to suggest points are in the same group.
#'
#' @return int number of groups 
#' @export
check_number_of_groups <- function(dist_mat, diff_eps){
  adjmatrix <- dist_mat <= diff_eps
  ig <- igraph::graph_from_adjacency_matrix(adjmatrix, mode = "undirected")
  
  return(igraph::components(ig, mode = "strong")$no) # number of clusters
}


#' Title
#'
#' @param dist_mat 
#' @param diff_eps 
#'
#' @return
#' @export
get_groups <- function(dist_mat, diff_eps){
  adjmatrix <- dist_mat <= diff_eps
  ig <- igraph::graph_from_adjacency_matrix(adjmatrix, mode = "undirected")
  
  return(igraph::components(ig, mode = "strong")$membership) # number of clusters
}

#' smart model clustering that examines mode clustering across a range of 
#' paramters
#'
#' @param g_list list of initial objects
#' @param g_names data frame of names of the initial objects
#' @param position positions in g_list that define the trajectory
#' @param sigma scalar for the sigma of the psuedo-density
#' @param range_eps vector of eps value (descending in size) - to be used in the
#' mode accession
#' @param range_maxT vector of maxT values (increasing in size) 
#' @param range_diff_eps  vector of diff eps value (descending in size) - to 
#' be using in calculating the number of groups
#' @param usefrac boolean, if distance function should multiple by number of 
#' points used
#' @param verbose verbose boolean, if should report progress
#'
#' @return data frame with parameter values and number of clusters found
#' @export
mode_cluster2 <- function(g_list, g_names, position, sigma, 
                          range_eps = 10^seq(-5,-8, by = .5), 
                          range_maxT = c(seq(30, 130, by = 20), 200),
                          range_diff_eps = 10^seq(-5,-8, by = .5), 
                          usefrac = T,
                          verbose = T,
                          collect_group_info = F){
  x_list_inner <- g_list %>% 
    lapply(function(df) df[,position])
  
  g_list_inner <- x_list_inner
  
  # need to work on...
  inner_range_eps <- range_eps
  
  
  full_info <- data.frame()
  
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "creating data [:bar] :percent eta: :eta",
      total = length(range_maxT), clear = FALSE, width = 38)
  }
  
  steps_counter <- 0
  if (collect_group_info){
    groups_cluster_list <- list()
  }
  
  for (i_maxT in 1:length(range_maxT)) {
    if (i_maxT == 1){
      lower_maxT <- 0
    } else {
      lower_maxT <- range_maxT[i_maxT-1]
    }
    maxT <- range_maxT[i_maxT] - lower_maxT
    
    out_list <- interior_clustering(x_list_inner = x_list_inner,
                                    g_list_inner = g_list_inner,
                                    sigma = sigma,
                                    range_eps = inner_range_eps,
                                    maxT = maxT)
    
    
    next_element <- T
    counter <- 1
    
    
    while (next_element) {
      string_eps <- names(out_list)[1]
      # for each element calculate 
      # (1) distance matrix
      dist_mat <- dist_matrix_innersq_direction(out_list[[1]][[3]] %>%
                                                  lapply(function(x) as.data.frame(x)),
                                                position = 1:ncol(out_list[[1]][[3]][[1]]),
                                                usefrac = T,
                                                verbose = F) 
      
      # (2) examine across cutoffs the number of groups
      number_clusters <- sapply(range_diff_eps, function(diff_eps) check_number_of_groups(dist_mat, diff_eps))
      
      # (2.5) group clustering
      for (diff_eps in range_diff_eps){
        inner_grouping_vec <- get_groups(dist_mat, diff_eps)
        inner_grouping_df <- g_names %>% mutate(groupings = inner_grouping_vec)
        groups_cluster_list[[paste0(string_eps, "-", diff_eps)]] <- inner_grouping_df
      }
      
      
      # (3) storage 
      steps_counter <- steps_counter + out_list[[1]][[2]]
      inner_info <- data.frame(maxT = range_maxT[i_maxT],
                               total_steps = steps_counter,
                               eps = as.numeric(string_eps),
                               diff_eps = range_diff_eps,
                               number_clusters = number_clusters)
      full_info <- rbind(full_info, inner_info)
      
      # process for next step in while -or- for loop
      old_g_list <- out_list[[1]][[3]]
      out_list <- out_list$next_element
      if (inherits(out_list, "logical") && !out_list){
        # next step in while-loop
        next_element <- F
        if (length(inner_range_eps) == 1){
          inner_range_eps <- c()
        } else {
          inner_range_eps <- inner_range_eps[2:length(inner_range_eps)]
        }
        # next step in for-loop
        g_list_inner <- old_g_list
      }
      
      counter <- counter + 1
      if (counter == 1000){
        next_element <- F
      }
      
      
      
      if (length(inner_range_eps) == 0){
        break
      }
    }
    
    if (verbose) {
      pb$tick()
    }
    
    if (length(inner_range_eps) == 0){
      break
    }
  }
  if (collect_group_info) {
    return(list(full_info, groups_cluster_list))
  } else {
    return(list(full_info))
  }
}



#' update a tm_radius matrix by smoothing about minim covering radi values
#'
#' @param tm_radius original list of lists, each containing a 
#' \code{min_cover_vec} and a \code{dist_mat} to be updated.
#' @param func smoothing function that creates a  function that takes in 
#' \code{x} and \code{y} and returns a smoothed \code{y_hat} for the x values.
#'
#' @return updated \code{tm_radius}
#' @export
update_tm_smooth <- function(tm_radius, func){
  n_l <- length(tm_radius)
  tm_radius_out <- list()
  for (l_idx in 1:n_l){
    min_vec <- tm_radius[[l_idx]]$min_cover_vec
    n_inner <- length(min_vec)
    xx <- 1:n_inner
    smooth_min_vec <- func(xx, min_vec)

    dist_mat <- tm_radius[[l_idx]]$dist_mat
    dist_mat[1,1] <- smooth_min_vec[1]
    if (n_inner > 1){
      for (r_idx in 2:n_inner){
        dist_mat[r_idx, r_idx] <- smooth_min_vec[r_idx]
        update_bool <- dist_mat[1:(r_idx-1),(r_idx-1)] < smooth_min_vec[r_idx]
        dist_mat[which(update_bool), r_idx] <- smooth_min_vec[r_idx]
        dist_mat[which(!update_bool), r_idx] <- dist_mat[which(!update_bool), r_idx-1]
      }
  }
    inner_list <- list(min_cover_vec = smooth_min_vec,
                       dist_mat = dist_mat)
    tm_radius_out[[l_idx]] <- inner_list
  }
  return(tm_radius_out)
}



#' creates a cdf for any discrete uniform distribution (for all integer values)
#'
#' @param min lowest integer with non-zero prob
#' @param max highest integer with non-zero prob
#'
#' @return cumulative distribution function 
#' @export
create_discrete_uniform_cdf <- function(min = 0, max = 1000){
  inner_function <- function(x) {
    lower <- floor(x)
    lower[lower > max] <- max
    cdf_value <- (lower - min+1)/(max-min+1)
    cdf_value[lower < min] <- 0
    return(cdf_value)
  }
  return(inner_function)
}


#' test how close conformal scores are to discrete uniform distribution
#' across functions
#'
#' @param info_x 
#' @param sim_names 
#' @param sim_list 
#' @param info_x_test 
#' @param pseudo_density_df 
#' @param grouping_df_list 
#' @param named_smooth_function_list 
#' @param tdm_sim 
#' @param mm_delta 
#'
#' @return
#' @export
uniform_check_wrapper <- function(info_x,
                                  sim_list,
                                  sim_names, 
                                  info_x_test,
                                  pseudo_density_df,
                                  grouping_df_list,
                                  named_smooth_function_list,
                                  tdm_sim,
                                  mm_delta,
                                  dist_between_array,
                                  verbose = T){
  # across groupings (radius calculations):
  # for each grouping structure:
  # 1. calculate analysis with different grouping structures (with mode = True)
  # 
  # for all options (1 fixed_nm, x fixed, 1 vary_nm, x vary:
  # if not vary:
  #   1. create tm_matrix for mode cluster
  # if vary:
  #   1. create new tm_matrices across all smoothing parameters
  # 
  # 2. calculate  conformal scores 
  # 3. get pvalue of test against a discrete uniform 0:1000
  
  # track: tm_matrix type, mode selection eps (eps + diff), smoothing (include "not")
  # less than: 4 * 10 * 6 * (8)
  #
  n_list <- length(grouping_df_list)
  n_func <- length(smooth_function_list)
  
  p_storage <- data.frame()
  
  out_groups_nm <- sim_names %>% dplyr::mutate(groupings = 1)
  simulation_info_out_nm <- inner_expanding_info(pseudo_density_df, 
                                                 out_groups_nm)
  simulation_info_df_nm <- simulation_info_out_nm[[1]]
  ordering_list_nm <- simulation_info_out_nm[[2]]
  
  
  pdiscreteunif <- create_discrete_uniform_cdf(min = 0, max = nrow(sim_names))
  
  for (g_idx in 0:n_list){
    if (g_idx == 0){ # nm
      # build tm ------
      ## fixed_nm -----
      tm_radius_fixed_nm <- inner_convert_single_radius_to_structure(mm_delta, 
                                                                     ordering_list_nm)
      ## vary_nm ------
      tm_radius_vary_nm <- coverage_down_mlist_save(inner_dist_mat = dist_between_array,
                                                    g_order_ll =  ordering_list_nm,
                                                    verbose = verbose)
      
      # calculate conformal scores ----
      conformal_df_fixed_nm <- inner_containment_conformal_score_mode_radius(
        df_row_group = info_x_test, 
        simulations_group_df = info_x, 
        data_column_names = c("x","y"),
        simulation_info_df = simulation_info_df_nm, # TODO look at 
        list_radius_info = tm_radius_fixed_nm, # diff
        list_grouping_id = ordering_list_nm, # diff
        verbose = verbose)
      
      conformal_df_vary_nm <- inner_containment_conformal_score_mode_radius(
        df_row_group = info_x_test, # change
        simulations_group_df = info_x, 
        data_column_names =  c("x","y"),
        simulation_info_df = simulation_info_df_nm, 
        list_radius_info = tm_radius_vary_nm, # diff
        list_grouping_id = ordering_list_nm, # diff
        verbose = verbose)
      
      # get pvalue ------
      suppressWarnings(ptest_fixed_nm <- chisq.test(
        table(factor(conformal_df_fixed_nm$containment_val,
                     levels = 0:1000)))$p.value)
      # do chisq.test instead
      p_storage <- rbind(p_storage ,
                         data.frame(type = "fixed_nm",
                                    grouping_var = NA,
                                    func_name = NA,
                                    pvalue = ptest_fixed_nm))
      
      suppressWarnings(ptest_vary_nm <- chisq.test(
        table(factor(conformal_df_vary_nm$containment_val,
                     levels = 0:1000)))$p.value)
      p_storage <- rbind(p_storage ,
                         data.frame(type = "vary_nm",
                                    grouping_var = NA, 
                                    func_name = NA,
                                    pvalue = ptest_vary_nm))
      
      # then smoothings ---
      for (f_idx in 1:n_func){ # for vary only
        # update tm_matrix from above ----
        tm_radius_vary_nm_update <- update_tm_smooth(tm_radius_vary_nm, 
                                                     func = named_smooth_function_list[[f_idx]])
        
        # calculate conformal scores ---
        conformal_df_vary_nm_update <- inner_containment_conformal_score_mode_radius(
          df_row_group = info_x_test,
          simulations_group_df = info_x, 
          data_column_names = c("x","y"),
          simulation_info_df = simulation_info_df_nm, 
          list_radius_info = tm_radius_vary_nm_update, # diff
          list_grouping_id = ordering_list_nm, # diff
          verbose = verbose)
        
        # get pvalue ---
        suppressWarnings(ptest_vary_nm_update <- chisq.test(
          table(factor(conformal_df_vary_nm_update$containment_val,
                       levels = 0:1000)))$p.value)
        p_storage <- rbind(p_storage,
                           data.frame(type = "vary_nm",
                                      grouping_var = NA,
                                      func_name = names(named_smooth_function_list)[f_idx],
                                      pvalue = ptest_vary_nm_update))
      }
      
      
      
    } else { # mode variants
      
      # for first without modification (name tm_matrix different than below)
      
      # build tm ------
      simulation_info_out <- inner_expanding_info(pseudo_density_df, 
                                                  grouping_df_list[[g_idx]])
      simulation_info_df <- simulation_info_out[[1]]
      ordering_list <- simulation_info_out[[2]]
      
      tm_radius_fixed <- inner_convert_single_radius_to_structure(mm_delta, 
                                                                  ordering_list)
      
      
      tm_radius_vary <- coverage_down_mlist_save(inner_dist_mat = dist_between_array,
                                                 g_order_ll = ordering_list)
      
      # calculate conformal scores ----
      
      conformal_df_fixed <- inner_containment_conformal_score_mode_radius(
        df_row_group = info_x_test, # change
        simulations_group_df = info_x, 
        data_column_names = c("x","y"),
        simulation_info_df = simulation_info_df, 
        list_radius_info = tm_radius_fixed, # diff
        list_grouping_id = ordering_list, # diff
        verbose = verbose)
      
      conformal_df_vary <- inner_containment_conformal_score_mode_radius(
        df_row_group = info_x_test,
        simulations_group_df = info_x, 
        data_column_names = c("x","y"),
        simulation_info_df = simulation_info_df, 
        list_radius_info = tm_radius_vary, # diff
        list_grouping_id = ordering_list, # diff
        verbose = verbose)
      
      # get pvalue ----
      suppressWarnings(ptest_fixed <- chisq.test(
        table(factor(conformal_df_fixed$containment_val,
                     levels = 0:1000)))$p.value)
      
      p_storage <- rbind(p_storage,
                         data.frame(type = "fixed",
                                    grouping_var = names(grouping_df_list)[g_idx],
                                    func_name = NA,
                                    pvalue = ptest_fixed))
      
      suppressWarnings(ptest_vary <- chisq.test(
        table(factor(conformal_df_vary$containment_val,
                     levels = 0:1000)))$p.value)
      
      p_storage <- rbind(p_storage,
                         data.frame(type = "vary",
                                    grouping_var = names(grouping_df_list)[g_idx],
                                    func_name = NA,
                                    pvalue = ptest_vary))
      
      # then smoothings ---
      for (f_idx in 1:n_func){ # for vary only
        # update tm_matrix from above ----
        tm_radius_vary_update <- update_tm_smooth(tm_radius_vary, 
                                                  func = named_smooth_function_list[[f_idx]])
        
        # calculate conformal scores ---
        conformal_df_vary_update <- inner_containment_conformal_score_mode_radius(
          df_row_group = info_x_test,
          simulations_group_df = info_x, 
          data_column_names = c("x","y"),
          simulation_info_df = simulation_info_df, 
          list_radius_info = tm_radius_vary_update, # diff
          list_grouping_id = ordering_list, # diff
          verbose = verbose)
        
        # get pvalue ---
        suppressWarnings(ptest_vary_update <- chisq.test(
          table(factor(conformal_df_vary_update$containment_val,
                       levels = 0:1000)))$p.value)
        p_storage <- rbind(p_storage,
                           data.frame(type = "vary",
                                      grouping_var = names(grouping_df_list)[g_idx],
                                      func_name = names(named_smooth_function_list)[f_idx],
                                      pvalue = ptest_vary_update))
      }
      
      
    }
  }
  
  return(p_storage)
}


#' coverage_down with a saved large array matrix
#'
#' @param inner_dist_mat array distance matrix, [i,j,k] is the minimum distance 
#' from all of jth curves points to ith curves' kth point 
#' @param g_order vector order of curves 
#' @param verbose bool
#'
#' @return list of \code{min_cover_vec} and \code{dist_mat}
#' @export
coverage_down_slist_save <- function(inner_dist_mat,
                                     g_order, verbose = T){
  n_l <- length(g_order)
  
  if (n_l == 1){
    inner_list <- list(min_cover_vec = c(0),
                       dist_mat = matrix(0))
    
    return(inner_list)
  }
  
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "create coverage_down info [:bar] :percent eta: :eta",
      total = n_l-1, clear = FALSE, width = 60)
  }
  
  # create min_cover_vec -----
  min_cover_vec <- rep(NA, n_l)
  min_cover_vec[1] <- 0
  
  mat_inner <- matrix(Inf, nrow = n_l, ncol = dim(inner_dist_mat)[3])
  
  for (l_idx in 2:n_l){
    #inner_g_order <- g_order[1:l_idx]
    if (l_idx == 2){
      mat_inner[l_idx,] <- inner_dist_mat[g_order[l_idx], g_order[1:(l_idx-1)],] 
    } else {
      mat_inner[l_idx,] <- inner_dist_mat[g_order[l_idx], g_order[1:(l_idx-1)],] %>%
        apply(2, min) 
    }
    
    
    if (l_idx == 2){
      mat_inner[1,] <- inner_dist_mat[g_order[1], g_order[2],] 
    } else {
      mat_inner[1:(l_idx-1),] <- array(
        c(mat_inner[1:(l_idx-1),], inner_dist_mat[g_order[1:(l_idx-1)],g_order[l_idx],]),
        dim = c(l_idx-1,
                dim(inner_dist_mat)[3],
                2)) %>% 
        apply(1:2, min)
    }
    
    
    min_cover_vec[l_idx] <- mat_inner[1:l_idx,] %>%
      max # apply(1, max) %>% max is the same
    
    if(verbose){
      pb$tick()
    }
  }
  
  # create dist_mat ----
  
  
  
  dist_mat <- matrix(NA, nrow = n_l, ncol = n_l)
  dist_mat[1,1] <- min_cover_vec[1]
  for (r_idx in 2:n_l){
    dist_mat[r_idx, r_idx] <- min_cover_vec[r_idx]
    update_bool <- dist_mat[1:(r_idx-1),(r_idx-1)] < min_cover_vec[r_idx]
    dist_mat[which(update_bool), r_idx] <- min_cover_vec[r_idx]
    dist_mat[which(!update_bool), r_idx] <- dist_mat[which(!update_bool), r_idx-1]
  }
  inner_list <- list(min_cover_vec = min_cover_vec,
                     dist_mat = dist_mat)
  
  return(inner_list)
}


#' create an array to captures the coverage_down ideas to be used for any ordering
#'
#' @param data_list list of data frames  
#' @param e_cols columns associated with location of path points
#' @param .e_cols_string boolean if e_cols is a string (not a promise)
#' @param verbose boolean if we should be verbose (this can take a long time)
#'
#' @return array[i,j,k] the minimum distance from all of jth curves points to 
#' ith curves' kth point 
#' @export
#'
#' @examples
coverage_down_list_save <- function(data_list, 
                                    e_cols, 
                                    .e_cols_string,
                                    verbose = T){
  n_obs <- length(data_list)
  n_of_obs <- unique(sapply(data_list, nrow))
  
  if (length(n_of_obs) != 1){
    stop("data.frames in data_list are not the same length.")
  }
  
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "creating distance all [:bar] :percent eta: :eta",
      total = n_obs^2, clear = FALSE, width = 60)
  }
  
  inner_dist_mat <- array(Inf, dim = c(n_obs,n_obs,n_of_obs))
  for (r_idx in 1:n_obs){ # ignore self?
    for (c_idx in 1:n_obs){
      both_out <- inner_min_dist_per_point(data_list[[r_idx]], 
                                           data_list[c_idx], e_cols,
                                           .e_cols_string = .e_cols_string)
      inner_dist_mat[c_idx,r_idx,1:n_of_obs] <- both_out$nn_info_l %>% as.vector()
      if (verbose){
        pb$tick()
      }
    }
  }
  return(inner_dist_mat)
}


#' coverage_down with a saved large array matrix, multiple order vectors
#'
#' @param inner_dist_mat array distance matrix, [i,j,k] is the minimum distance 
#' from all of jth curves points to ith curves' kth point 
#' @param g_order_ll list of order of curves 
#'
#' @return list of \code{min_cover_vec} and \code{dist_mat}
#' @export
coverage_down_mlist_save <- function(inner_dist_mat,
                                     g_order_ll,
                                     verbose = F){
  out_list <- lapply(g_order_ll, function(g_order) coverage_down_slist_save(inner_dist_mat,
                                                                            g_order,
                                                                            verbose = verbose))
  
  return(out_list)
}
