# run: Rscript 02-visualize_hyperparameters.R 1000 default

suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))

input_args <- commandArgs(trailingOnly=TRUE)
num_sims <- as.numeric(input_args[1]) # 1000

input_sigma_info <- input_args[2]
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




data_list <- list()
vis_list <- list()
path_vis_list <- list()
total_time_all <- rep(0,3)
expected_number_clusters <- data.frame(x = c(0,.64,.85), num = c(1,2,2))
x_storage_values <- c(0,.64,.85) # needed to do .01 instead of 1 before 01-... was written
names(x_storage_values) <- as.character(c(0,.64,.85))

for (x_value in c(0,.64,.85)){ 
  # load in data
  
  x_storage_value <- x_storage_values[as.character(x_value)]
  saved_out <- load(file = paste0("data/tuning_info_x",x_storage_value*100,"_",
                                  num_sims, "_", input_sigma_info_str,
                                  "_3.Rdata"))
  eval(parse(text = paste("sigma_info_list_x <-", saved_out[1])))
  total_time_all[expected_number_clusters$x == x_value] <- total_time
  
  expected_number_clusters_inner <- expected_number_clusters$num[
    expected_number_clusters$x == x_value]
  
  data_list[[as.character(x_value)]] <- sigma_info_list_x %>%
    do.call(rbind, .) %>%
    mutate(.sigma_string = factor(.sigma_string, 
                                  levels = range_sigma),
           correct_num_clusters = number_clusters == expected_number_clusters_inner)
  

  
  vis_list[[as.character(x_value)]] <- data_list[[as.character(x_value)]] %>%
    ggplot() +
    geom_tile(aes(x = -log10(eps), y = -log10(diff_eps),
                  fill = correct_num_clusters)) +
    labs(x = "-log10(eps)", y = "-log10(diff_eps)",
         title = paste0("X = ", x_value),
         fill = paste0("correct num\nclusters: ", expected_number_clusters_inner))  +
    facet_wrap(~.sigma_string)
  
  
  data2d <- data.frame(x = c(1,0,0,1),
                       y = c(0,1,0,0),
                       z = c(0,0,1,0)) %>% 
    EpiCompare::get_xy_coord(xyz_col = c("x","y", "z"))
  
  path_vis_list[[as.character(x_value)]] <- ggvis_curves$data %>% 
    ggplot() +
    geom_path(aes(x = x, y = y, group = idx), alpha = 20/num_sims) +
    geom_path(data = data2d, aes(x=x,y=y), color = "black") +
    theme_minimal() +
    labs(x="", y="", title = paste0("X val: ", x_value, "\nnum sims: ", num_sims))
}




vis_all <- grid.arrange(grobs = vis_list, nrow = 1,
                        top = "Grid search relative to number of modes (maxT = 200)",
                        bottom = "Optimized over: sigma %, convergence eps, and clustering diff_eps")

ggsave(vis_all, filename = paste0("data/vis_all_",num_sims,"_", 
                                  input_sigma_info_str, "_3.pdf"), 
       width = 30, height = 10)


suppressMessages(devtools::load_all())
vis_all2 <- grid.arrange(grobs = path_vis_list, nrow = 1)

ggsave(vis_all2, filename = paste0("data/vis_all_paths_",num_sims,"_", 
                                  input_sigma_info_str, "_3.pdf"), 
       width = 30, height = 10)






for (x_value in c(0,.64,.85)){
  if (x_value == 0){
    data_combo <- data_list[[as.character(x_value)]] %>%
      select(eps, diff_eps, .sigma_string, correct_num_clusters)
  } else {
    data_inner <- data_list[[as.character(x_value)]] %>%
      select(eps, diff_eps, .sigma_string, correct_num_clusters)
    data_combo <- data_combo %>% full_join(data_inner, by = c("eps", "diff_eps", ".sigma_string"),
                                           suffix = c("", as.character(x_value)))
  }
}

nafalse <- function(x){
  x[is.na(x)] <- FALSE
  return(x)
}

data_combo <- data_combo %>%
  mutate(correct_num_across = nafalse(correct_num_clusters) +
           nafalse(correct_num_clusters0.64) +
           nafalse(correct_num_clusters0.85))

vis_combo <- data_combo %>%
  ggplot() +
  geom_tile(aes(x = -log10(eps), y = -log10(diff_eps),
                fill = factor(correct_num_across))) +
  labs(x = "-log10(eps)", y = "-log10(diff_eps)",
       title = "Parameters with correct amount of modes",
       subtitle = "(sum across X values, NA = False)",
       fill = "# of correct modes")  +
  facet_wrap(~.sigma_string) +
  scale_fill_manual(values = c("0"="grey70",
                               "1" = '#fee8c8',
                               "2" = '#fdbb84',
                               "3" = '#e34a33')) +
  theme(legend.position = "bottom")

ggsave(vis_combo, filename = paste0("data/vis_combo_",num_sims,"_", 
                                   input_sigma_info_str, "_3.pdf"), 
       width = 10, height = 10)

save(vis_all, vis_all2, vis_combo, data_combo,
     file = paste0("data/ggvis_",num_sims,"_", 
                       input_sigma_info_str, "_3.Rdata"))


