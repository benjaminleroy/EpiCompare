# # run: Rscript 011-tune_hyperparameters.R [1-3] 1000 default
# 3rd option either : "default" or "5:85::5" style (start at 5, end at 85, with steps of length 5)
#
# expected run time:  a lot (~7.5 hours) - for 1000

suppressMessages(library(tidyverse))

suppressMessages(devtools::load_all())
suppressMessages(devtools::load_all("../../../simulationBands/"))

start_time <- Sys.time()
# input parameters -----------------------------
input_args <-  commandArgs(trailingOnly=TRUE) #c(1,1000, "default") 
x_value <- c(0, .64, .85)[as.numeric(input_args[1])]

seed_val <- c(707, 412, 1)[as.numeric(input_args[1])]
set.seed(seed_val)

#calibration_idx <- 1

# initial global parameters --------------------
n_simulations <- as.numeric(input_args[2]) #1000
number_points <- 150
compression_size <- 50

input_sigma_info <- input_args[3]
if (input_sigma_info == "default"){
  range_sigma <- paste0(c(1,3,5, seq(10, 85, by = 5)), "%")
  input_sigma_info_str <- input_sigma_info
} else {
  sigma_range_info <- input_sigma_info %>% stringr::str_split(":{1,2}", simplify = T) %>% 
    as.numeric()
  range_sigma <- paste0(seq(sigma_range_info[1], 
                            sigma_range_info[2], 
                            by = sigma_range_info[3]),
                        "%")
  input_sigma_info_str <- paste0(sigma_range_info, collapse = "-")
}


create_smooth_function <- function(df){
  
  eval(parse(text = paste0("smooth_function <- function(x,y){
    tryCatch({
      inner_ss <- smooth.spline(x,y, df = ",
      df,
      ")
      return(predict(inner_ss,x)$y)},
      error = function(cond){
        message(sprintf(\"returning y: error in size of x (%i < 4)\", length(x)))
        return(y)
      },
      warning = function(cond){
        message(sprintf(paste(\"returning y: error in size of x\",
                              \"relative to size of df (%i < %i)\"),
                        length(x), ", df, "))
        return(y)
      }
    )
  }")))

  return(smooth_function)
}

inner_moving <- function(x, y, window = 5, fill_left = T,
                         fun = min){
  n <- length(y)
  if (n < window){
    message(sprintf(paste("returning y: error in size of x",
                          "relative to size of window (%i < %i)"),
                    n, window))
    return(y)
  }
  
  y_smooth <- rep(NA, n)
  for (idx in window:n){
    y_smooth[idx] <- fun(y[(idx-window+1):idx])
  }
  
  if (fill_left){
    y_smooth[1:(window-1)] <- y_smooth[window]
  }
  
  return(y_smooth)
}

create_min_smooth_function <- function(window){
  eval(parse(text = paste0("out_function <- function(x,y){
    inner_moving(x,y, window = ", window,", fill_left = T, fun = min)
  }")))
  return(out_function)
}


smooth_function_list <- list()

for (df in c(5,10,25,50)){
  smooth_function_list[[paste0("df", df)]] <- create_smooth_function(df)
}
for (window in c(5,10,25,50)){
  smooth_function_list[[paste0("window",window)]] <- create_min_smooth_function(window)
}
range_eps2 <- 10^seq(-5,-10, by = -.5)
range_diff_eps <- 10^seq(-5,-8, by = -.5)
range_maxT <- c(200)


uniform_p_test <- function(cs_vec, range = c(0,1000), size_bins = c(1,10,25,50)){
  p_values <- rep(NA, length(size_bins))
  idx <- 1
  for (n_bins in size_bins){
    cs_binned <- cs_vec %>% cut(breaks = seq(range[1], range[2], by = n_bins)) %>%
      as.numeric() %>%
      factor(levels = 1:(ceiling((range[2]-range[1])/n_bins))) %>%
      table() 
    p <- rep(n_bins/(range[2]-range[1]), ceiling((range[2]-range[1])/n_bins))
    p_values[idx] <- chisq.test(cs_binned, p = p)$p.value
    idx <- idx + 1
  }
  out <- data.frame(size_bins = size_bins,
                    p_values = p_values)
  return(out)
}

uniform_kl_div <- function(cs_vec, range = c(0,1000), size_bins = c(1,10,25,50)){
  kl_values <- rep(NA, length(size_bins))
  idx <- 1
  for (n_bins in size_bins){
    cs_binned <- cs_vec %>% cut(breaks = seq(range[1], range[2], by = n_bins)) %>%
      as.numeric() %>%
      factor(levels = 1:(ceiling((range[2]-range[1])/n_bins))) %>%
      table() 
    cs_binned_prop <- cs_binned / length(cs_vec)
    if (sum(cs_binned_prop != 0) > 0){
      min_non_zero <- min(cs_binned_prop[cs_binned_prop != 0])
      added_prop <- min_non_zero/100 * length(cs_binned_prop)
      cs_binned_prop_corrected <- (cs_binned_prop + min_non_zero/100 )/ (1+added_prop)
    } else {
      cs_binned_prop_corrected <- cs_binned_prop
    }

    

    
    p <- rep(n_bins/(range[2]-range[1]), ceiling((range[2]-range[1])/n_bins))
    kl_values[idx] <- sum(cs_binned_prop_corrected * log(cs_binned_prop_corrected/p))
    idx <- idx + 1
  }
  out <- data.frame(size_bins = size_bins,
                    kl_values = kl_values)
}

uniform_multiple_check <- function(cs_vec, range = c(0,1000), size_bins = c(1,10,25,50)){
  p_values <- rep(NA, length(size_bins))
  kl_values <-  rep(NA, length(size_bins))
  l2_values <- rep(NA, length(size_bins))
  idx <- 1
  for (n_bins in size_bins){
    cs_binned <- cs_vec %>% cut(breaks = seq(range[1], range[2], by = n_bins)) %>%
      as.numeric() %>%
      factor(levels = 1:(ceiling((range[2]-range[1])/n_bins))) %>%
      table() 
    # for kl div
    cs_binned_prop <- cs_binned / length(cs_vec)
    if (sum(cs_binned_prop != 0) > 0){ # correction for log(0) occurring - a sight non-zero value (min_non_zero/100)
      min_non_zero <- min(cs_binned_prop[cs_binned_prop != 0])
      added_prop <- min_non_zero/100 * length(cs_binned_prop)
      cs_binned_prop_corrected <- (cs_binned_prop + min_non_zero/100 )/ (1+added_prop)
    } else {
      cs_binned_prop_corrected <- cs_binned_prop
    }
    
    
    p <- rep(n_bins/(range[2]-range[1]), ceiling((range[2]-range[1])/n_bins))
    
    p_values[idx] <- chisq.test(cs_binned, p = p)$p.value
    kl_values[idx] <- sum(cs_binned_prop_corrected * log(cs_binned_prop_corrected/p))
    l2_values[idx] <- mean((cs_binned_prop-p)^2)
    
    
    idx <- idx + 1
  }
  out <- data.frame(size_bins = size_bins,
                    p_values = p_values,
                    kl_values = kl_values,
                    l2_values = l2_values)
  return(out)
}

# testthat::test_that("test uniform_multiple_check", {
#   set.seed(1)
#   cs_vec <- ceiling(runif(.99, 1000, n = 1000))
#   suppressWarnings(out <- uniform_multiple_check(cs_vec))
#   
#   testthat::expect_true(all(out$p_values > 0))
#   testthat::expect_true(all(out$kl_values < 2))
#   
#   
#   cs_vec_bad <- ceiling(rbinom(size = 1000, prob = .5, n = 1000))
#   suppressWarnings(out_bad <- uniform_multiple_check(cs_vec_bad))
#   
#   testthat::expect_true(all(out_bad$kl_values > 2))
#   testthat::expect_true(all(out_bad$p_values < .001))
#   
# })


# \/ uniform function *** ***  ---------

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
uniform_check_wrapper2 <- function(info_x,
                                   sim_list,
                                   sim_names, 
                                   info_x_test_list,
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
  
  #p_storage <- data.frame()
  
  out_groups_nm <- sim_names %>% dplyr::mutate(groupings = 1)
  simulation_info_out_nm <- inner_expanding_info(pseudo_density_df, 
                                                 out_groups_nm)
  simulation_info_df_nm <- simulation_info_out_nm[[1]]
  ordering_list_nm <- simulation_info_out_nm[[2]]
  
  
  #pdiscreteunif <- create_discrete_uniform_cdf(min = 0, max = nrow(sim_names))
  

  n_sims <- n_simulations
  number_points <- 150
  n_sims_containment <- 300
  
  
  verbose = F
  # first pass (nm) --------
  p_storage_0 <- data.frame()
  
  # build tm ------
  ## fixed_nm -----
  tm_radius_fixed_nm <- inner_convert_single_radius_to_structure(mm_delta, 
                                                                 ordering_list_nm)
  ## vary_nm ------
  tm_radius_vary_nm <- coverage_down_mlist_save(inner_dist_mat = dist_between_array,
                                                g_order_ll =  ordering_list_nm,
                                                verbose = verbose)
  
  for (n_sim in 1:length(info_x_test_list)){
    info_x_test <- info_x_test_list[[n_sim]]
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
    suppressWarnings(
      ptest_fixed_nm <- uniform_multiple_check(conformal_df_fixed_nm$containment_val))
    # do chisq.test instead
    p_storage_0 <- rbind(p_storage_0 ,
                         ptest_fixed_nm %>%
                           dplyr::mutate(type = "fixed_nm",
                                    grouping_var = NA,
                                    func_name = NA,
                                    n_sim = n_sim))
    
  
    suppressWarnings(
      ptest_vary_nm <- uniform_multiple_check(conformal_df_vary_nm$containment_val))
    p_storage_0 <- rbind(p_storage_0 ,
                         ptest_vary_nm %>%
                           dplyr::mutate(type = "vary_nm",
                                         grouping_var = NA,
                                         func_name = NA,
                                         n_sim = n_sim))
    
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
      
      suppressWarnings(
        ptest_vary_nm_update <- uniform_multiple_check(conformal_df_vary_nm_update$containment_val))
      # do chisq.test instead
      p_storage_0 <- rbind(p_storage_0 ,
                           ptest_vary_nm_update %>%
                             dplyr::mutate(type = "vary_nm",
                                           grouping_var = NA,
                                           func_name = names(named_smooth_function_list)[f_idx],
                                           n_sim = n_sim))
    
  }
  }
  
  
  
  
  
  
  
  # objects to bring in ---------
  .export_string_into_function <- c("info_x", "sim_list", "sim_names", 
                                    "info_x_test_list",
                                    "pseudo_density_df", "grouping_df_list", 
                                    "named_smooth_function_list","tdm_sim", 
                                    "mm_delta", "dist_between_array")
  .export_string_interal <- c("n_sims", "number_points", "n_sims_containment",
                              "n_func", "simulation_info_df_nm", 
                              "ordering_list_nm")
  
  # not sure why this is needed
  .export_epicompare_functions <- c("inner_convert_single_radius_to_structure",
                         "coverage_down_mlist_save", "inner_containment_conformal_score_mode_radius",
                         "inner_expanding_info", "coverage_down_slist_save")
  
  .export_local_function <- c("create_smooth_function", "create_min_smooth_function")#, "inner_moving")

  
  .export_all = c(.export_string_into_function, .export_string_interal, 
                  .export_epicompare_functions, .export_local_function)
  
  
  suppressWarnings(suppressMessages(library(parallel)))
  suppressWarnings(suppressMessages(library(foreach)))
  suppressWarnings(suppressMessages(library(doSNOW)))
  max_cores <- parallel::detectCores()
  cl <- makeCluster(4)
  # https://stackoverflow.com/questions/17879766/variable-scope-in-boot-in-a-multiclustered-parallel-approach?noredirect=1&lq=1
  clusterExport(cl, c('inner_moving', "uniform_multiple_check")) # used to avoid problems with inner_moving function not being exported...
  registerDoSNOW(cl)
  iterations <- n_list
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  p_storage_global <- foreach(g_idx = 1:n_list,
                              .combine = "rbind",
                              .options.snow = opts,
                              .export = .export_all,
                              .packages = c("dplyr", "tidyr",
                                            "simulationBands",
                                            "EpiCompare")) %dopar% {
                                              # head junk ----
                                            
                                              verbose <- F
                                              p_storage <- data.frame()
                                              
                                              # mode variants ---------
                                              
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
                                              for (n_sim in 1:length(info_x_test_list)){
                                                info_x_test <- info_x_test_list[[n_sim]]
                                                
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
                                                suppressWarnings(
                                                  ptest_fixed <- uniform_multiple_check(conformal_df_fixed$containment_val))
                                                # do chisq.test instead
                                                p_storage <- rbind(p_storage ,
                                                                   ptest_fixed %>%
                                                                       dplyr::mutate(type = "fixed",
                                                                                     grouping_var = names(grouping_df_list)[g_idx],
                                                                                     func_name = NA,
                                                                                     n_sim = n_sim))
                                                
                                                suppressWarnings(
                                                  ptest_vary <- uniform_multiple_check(conformal_df_vary$containment_val))
                                                # do chisq.test instead
                                                p_storage <- rbind(p_storage ,
                                                                   ptest_vary %>%
                                                                     dplyr::mutate(type = "fixed",
                                                                                   grouping_var = names(grouping_df_list)[g_idx],
                                                                                   func_name = NA,
                                                                                   n_sim = n_sim))
                                                
  
                                                
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
                                                  suppressWarnings(
                                                    ptest_vary_update <- uniform_multiple_check(conformal_df_vary_update$containment_val))
                                                  # do chisq.test instead
                                                  p_storage <- rbind(p_storage ,
                                                                     ptest_vary_update %>%
                                                                       dplyr::mutate(type = "vary",
                                                                                     grouping_var = names(grouping_df_list)[g_idx],
                                                                                     func_name = names(named_smooth_function_list)[f_idx],
                                                                                     n_sim = n_sim))
                                                  
                                                  
                                                }
                                              }
                                              
                                              
                                              
                                              return(p_storage)
                                            }
  
  p_out <- rbind(p_storage_0,
                 p_storage_global)
  parallel::stopCluster(cl)
  return(p_out)
}

# ^ uniform function *** *** ----------------


# GENERATE DATA ----------


info_x <- data_generation(x_inner = x_value,
                          n_sims_containment = n_simulations,
                          number_points = number_points,
                          verbose = T)

# creating compressed data
compressed_info_x <- info_x %>%
  filament_compression(data_columns = 
                           c(1:ncol(info_x))[names(info_x) %in% c("x","y")],
                         number_points = compression_size)

# calculating distance matrix (for sigma values)
tdm_sims_compress_x <- dist_matrix_innersq_direction(compressed_info_x,
                                                      position = c(1:ncol(compressed_info_x))[names(compressed_info_x) %in% c("x","y")],
                                                      usefrac = T,
                                                      tdm_out = T,
                                                      verbose = T)

# calculating distance matrix (for sigma values)
# tdm_sims_x <- dist_matrix_innersq_direction(info_x,
#                                              position = c(1:ncol(info_x))[names(info_x) %in% c("x","y")],
#                                              usefrac = T,
#                                              tdm_out = T,
#                                              verbose = T)


# test info -----
compressed_info_x_test_list <- list()
for (n_sims in 1:20){
  print(paste("sim test: ", n_sims))
  info_x_test <- data_generation(x_inner = x_value,
                                 n_sims_containment = n_simulations,
                                 number_points = number_points,
                                 verbose = T)
  
  # creating compressed data
  compressed_info_x_test <- info_x_test %>%
    filament_compression(data_columns =
                           c(1:ncol(info_x_test))[names(info_x_test) %in% c("x","y")],
                         number_points = compression_size)
  
  compressed_info_x_test_list[[n_sims]] <- compressed_info_x_test
}




# \/ visualization of curves ---------------

data2d <- data.frame(x = c(1,0,0,1),
                     y = c(0,1,0,0),
                     z = c(0,0,1,0)) %>% 
  get_xy_coord(xyz_col = c("x","y", "z"))

ggvis_curves <- compressed_info_x %>% ggplot() +
  geom_path(aes(x = x, y = y, group = idx), alpha = 20/n_simulations) +
  geom_path(data = data2d, aes(x=x,y=y), color = "black") +
  theme_minimal() +
  labs(x="", y="", title = paste0("X val: ", x_value, "\nnum sims: ", n_simulations))
  
# ^ visualization of curves -----------------



g_list <- compressed_info_x %>% group_split()
g_names <-compressed_info_x %>% group_keys()

sim_list <- g_list#info_x %>% group_split()
sim_names <- g_names#info_x %>% group_keys()

#assertthat::assert_that(assertthat::are_equal(sim_names, g_names))

position <- which(names(g_list[[1]]) %in% c('x','y'))


sigma_info_list_x <- list()
ggvis_density_list <- list() 
uniform_info_list <- list()

dist_between_array <- coverage_down_list_save(sim_list, e_cols = c("x","y"),
                        .e_cols_string = F, verbose = T)


# START OF ANALYSIS -----------

for (.sigma_string in range_sigma){
  
  print(.sigma_string)
  # converting .sigma_string to numerical value.
  
  sigma_lower <- check_character_percent(.sigma_string, ".sigma_string")
  sigma_sizes <- sapply(sigma_lower + .05*(0:5), function(v) min(v, 1))
  
  percentage_inner <- sigma_sizes[stats::quantile(as.matrix(tdm_sims_compress_x), sigma_sizes) > 0][1]
  
  sigma_val <- stats::quantile(as.matrix(tdm_sims_compress_x), percentage_inner)
  sigma <- sigma_val
  
  # (1) model clustering ----------
  sigma_info_list <- mode_cluster2(g_list = g_list, g_names = g_names,
                                 sigma = sigma,
                                 position = position,
                                 range_eps = range_eps2,
                                 range_maxT = range_maxT,
                                 range_diff_eps = range_diff_eps,
                                 verbose = T, usefrac = T, 
                                 collect_group_info = T)
  
  sigma_info_df <- sigma_info_list[[1]]
  group_info_list <- sigma_info_list[[2]]
  
  # ^ need to update function so that it can return grouping information...
  
  sigma_info_df <- sigma_info_df %>% mutate(.sigma_string = .sigma_string,
                                            sigma = sigma_val,
                                            sigma_name = names(sigma_val))
  
  sigma_info_list_x[[.sigma_string]] <- sigma_info_df
  
  # (2) pseudo_density_df --------
  pd_df <- EpiCompare::distance_psuedo_density_function(tdm_sims_compress_x,
                                                       sigma = sigma, 
                                                       df_out = T)
  
  
  ggvis_density_list[[.sigma_string]] <- compressed_info_x %>%
    left_join(pd_df, 
              by = c("x_original", "R0", "beta",
                     "gamma", "idx", "set_idx")) %>% 
    ggplot() +
    geom_path(aes(x = x, y = y, group = idx, color = psuedo_density), 
              alpha = 20/n_simulations) +
    geom_path(data = data2d, aes(x=x,y=y), color = "black") +
    theme_minimal() +
    labs(x="", y="", title = paste0("X val: ", x_value, "\nnum sims: ", n_simulations,
                                    "\nsigma :", .sigma_string)) +
    scale_color_gradient2(midpoint = median(pd_df$psuedo_density))
  
  # (3) test-uniform relative to smoothness ----------
  
  # (3.1) first get mm_delta ---------
  

  
  top_points <- top_curves_to_points(compressed_info_x,
                                     tidy_dm = tdm_sims_compress_x,
                                     alpha = .2,
                                     quantile_func = distance_psuedo_density_function,
                                     sigma = sigma) # 80% curve remain
  
  assertthat::assert_that(nrow(top_points) > 0,
                          msg = paste("the number of simulations are so few",
                                      "thatthere are too few to estimate the",
                                      "shared radius."))
  
  mm_delta <- get_delta_nn(top_points[,c("x","y")])
  
  # (3.2) now get p-values across different smoothings ------
  
  # grouping_df_list <- group_info_list
  
  uniform_info_df <- uniform_check_wrapper2(compressed_info_x,
                                    sim_list = sim_list,
                                    sim_names = sim_names, 
                                    compressed_info_x_test_list,
                                    pseudo_density_df = pd_df,
                                    grouping_df_list = group_info_list,
                                    named_smooth_function_list = smooth_function_list,
                                    tdm_sim = tdm_sims_compress_x,
                                    mm_delta,
                                    dist_between_array = dist_between_array,
                                    verbose = F)
  
  
  uniform_info_df <- uniform_info_df %>% 
    mutate(.sigma_string = .sigma_string,
           sigma = sigma_val,
           sigma_name = names(sigma_val))
  
  
  uniform_info_list[[.sigma_string]] <- uniform_info_df
}


total_time <- Sys.time() - start_time
save(sigma_info_list_x, total_time, ggvis_curves, ggvis_density_list,
     uniform_info_list,
     file = paste0("data/tuning_info_x",x_value*100,
                   "_", n_simulations, "_", input_sigma_info_str,
                   "_3.Rdata"))
