
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
    mutate(beta = .1,
           gamma = .1/R0,
           idx = paste0("inner_truth", 1:n_sims_containment))
  
  truth_paths <- data.frame(x_original = as.numeric(rep(NA, number_points*n_sims_containment)),
                            R0 = as.numeric(rep(NA, number_points*n_sims_containment)),
                            beta = as.numeric(rep(NA, number_points*n_sims_containment)),
                            gamma = as.numeric(rep(NA, number_points*n_sims_containment)),
                            idx = as.character(rep(NA, number_points*n_sims_containment)),
                            set_idx = as.character(rep(NA, number_points*n_sims_containment)),
                            x = as.numeric(rep(NA, number_points*n_sims_containment)),
                            y = as.numeric(rep(NA, number_points*n_sims_containment)))
  
  
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
                      inner_paths)
    
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
                          verbose = T){
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
  return(full_info)
}
