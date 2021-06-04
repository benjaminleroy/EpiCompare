# # run: Rscript 01-tune_hyperparameters.R [1-3] 1000 default
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
  
  smooth_function <- function(x,y){
    tryCatch({
      inner_ss <- smooth.spline(x,y, df = df)
      return(predict(inner_ss,x)$y)},
      error = function(cond){
        message(sprintf("returning y: error in size of x (%i < 4)", length(x)))
        return(y)
      },
      warning = function(cond){
        message(sprintf(paste("returning y: error in size of x",
                              "relative to size of df (%i < %i)"),
                        length(x), df))
        return(y)
      }
    )
    

  }
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
  out_function <- function(x,y){
    inner_moving(x,y, window = window, fill_left = T, fun = min)
  }
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

# generate data


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
info_x_test <- data_generation(x_inner = x_value,
                               n_sims_containment = n_simulations,
                               number_points = number_points,
                               verbose = T)

# creating compressed data
compressed_info_x_test <- info_x_test %>%
  filament_compression(data_columns =
                         c(1:ncol(info_x_test))[names(info_x_test) %in% c("x","y")],
                       number_points = compression_size)



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
  
  uniform_info_df <- uniform_check_wrapper(compressed_info_x,
                                    sim_list = sim_list,
                                    sim_names = sim_names, 
                                    compressed_info_x_test,
                                    pseudo_density_df = pd_df,
                                    grouping_df_list = group_info_list,
                                    named_smooth_function_list = smooth_function_list,
                                    tdm_sim,
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
     file = paste0("data/tuning_info_x",x_value*100,
                   "_", n_simulations, "_", input_sigma_info_str,
                   "_3.Rdata"))
