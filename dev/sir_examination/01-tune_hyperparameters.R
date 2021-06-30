# # run: Rscript 01-tune_hyperparameters.R [1-3] 1000 
# 3rd option either : "default" or "5:85::5" style (start at 5, end at 85, with steps of length 5)
#
# expected run time:  a lot (~7.5 hours) - for 1000

suppressMessages(library(tidyverse))

suppressMessages(devtools::load_all())
suppressMessages(devtools::load_all("../../../simulationBands/"))

start_time <- Sys.time()
# input parameters -----------------------------
input_args <- c(1,1000, "default") #commandArgs(trailingOnly=TRUE)
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
  range_sigma <- paste0(c(1:5, seq(10, 85, by = 5)), "%")
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

position <- which(names(g_list[[1]]) %in% c('x','y'))


sigma_info_list_x <- list()
ggvis_density_list <- list() 

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
    left_join(pd_df, by = c() )%>% ggplot() +
    geom_path(aes(x = x, y = y, group = idx, color = psuedo_density), alpha = 20/n_simulations) +
    geom_path(data = data2d, aes(x=x,y=y), color = "black") +
    theme_minimal() +
    labs(x="", y="", title = paste0("X val: ", x_value, "\nnum sims: ", n_simulations,
                                    "\nsigma :", .sigma_string))
  
  # (3) test-uniform relative to smoothness ----------
  
  # across groupings (radius calculations):
  # 1. calculate all grouping structures
  # 2.
  
  # uniform_check_wrapper <- function(compressed_info_x,
  #                       compressed_info_x_test,
  #                       grouping_list,
  #                       pd_df = pdf_df,
  #                       group_df = ){
  #   # needs to 
  #   # 1. 
  # }
  
  
  
}


total_time <- Sys.time() - start_time
save(sigma_info_list_x, total_time, ggvis_curves, ggvis_density_list,
     file = paste0("data/sigma_info_x",x_value*100,
                   "_", n_simulations, "_", input_sigma_info_str,
                   "_3.Rdata"))
