# run: Rscript 3-test_set_conformal_scores.R [1-101] 1000 T
# expected run time 1.6 hours (shoot for 2:30)

input_args <- commandArgs(trailingOnly=TRUE)
test_idx <- as.numeric(input_args[1])

# global parameters ----------

n_simulations <- as.numeric(input_args[2])
n_sims_containment <- 300
number_points <- 150

return_min_bool <- as.logical(input_args[3])



# loading data and libraries (minimal data) ----



suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(devtools))

# load 
user_name <- Sys.info()["user"]

if (user_name == "benjaminleroy"){ # personal computer
  if (getwd() != "/Users/benjaminleroy/Documents/CMU/research/EpiCompare/dev/sir_examination") {
    setwd("/Users/benjaminleroy/Documents/CMU/research/EpiCompare/dev/sir_examination")
  }
  
  suppressMessages(load_all("../../")) # just EpiCompare@ben_dev_and_thesis
  suppressMessages(load_all("../../../simulationBands"))
} else { # vera.psc.edu
  suppressMessages(library(EpiCompare))
  suppressMessages(library(simulationBands))
}

# hyperparameter load -----------
input_sigma_info_str <- "default"
load(paste0("data/selected_parameters_", 
            n_simulations, "_", input_sigma_info_str,
            "_3.Rdata"))

# test set conformal scores ------------------

set.seed(123)
n_test_grid <- 101
test_set_r0 <- tibble(x = seq(0,1, length.out = n_test_grid))



start_time <- Sys.time()

x_inner <- test_set_r0$x[test_idx]

                                     
# sim_df filaments ----

sim_df <- data.frame(x = rep(x_inner, n_simulations))
sim_df$R0 <- 1/7*simulationBands::lei_wasserman_data_conditional_simulate(2*(x_inner-.5),
                                                                         n = n_simulations)[[1]]$sim+2

sim_df <- sim_df %>%
 dplyr::mutate(beta = .1,
               gamma = .1/R0,
               idx = paste0("inner", 1:n_simulations))

                                     
sim_paths <- data.frame(x_original = as.numeric(rep(NA, number_points*n_simulations)),
                       R0 = as.numeric(rep(NA, number_points*n_simulations)),
                       beta = as.numeric(rep(NA, number_points*n_simulations)),
                       gamma = as.numeric(rep(NA, number_points*n_simulations)),
                       idx = as.character(rep(NA, number_points*n_simulations)),
                       set_idx = as.character(rep(NA, number_points*n_simulations)),
                       x = as.numeric(rep(NA, number_points*n_simulations)),
                       y = as.numeric(rep(NA, number_points*n_simulations)))
                                     
                                     
                                     
                                     
for (r_idx in 1:n_simulations){
 inner_info <- sim_df[r_idx, ]
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
 
 sim_paths[global_r_idx, ] <- inner_df
}

sim_paths <- sim_paths %>% 
 dplyr::group_by(x_original, R0, beta, gamma, idx, set_idx)

# inner_truth filaments ----

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
}

truth_paths <- truth_paths %>% 
 dplyr::group_by(x_original, R0, beta, gamma, idx, set_idx)

           
if (F){
   sim_paths %>% 
      mutate(obs = "sim") %>% 
      rbind(truth_paths %>%
               mutate(obs = "truth")) %>%
      ggplot() +
      geom_path(aes(x=x,y=y, group = set_idx, 
                    color = obs), alpha = .1) +
      facet_wrap(obs~.)
}                          

conformal_score <- simulation_based_conformal4(truth_grouped_df = truth_paths,
                                               simulations_grouped_df = sim_paths,
                                               data_column_names = c("x", "y"),
                                               number_points = Inf,
                                               .to_simplex = F,
                                               .maxT = selected_parameters$.maxT,
                                               .sigma_string = selected_parameters$.sigma_string,
                                               .diff_eps = selected_parameters$.diff_eps,
                                               .eps = selected_parameters$.eps,
                                               verbose = F,
                                               return_min = return_min_bool)                                   
                                     
                                  
elapsed_time <- Sys.time() - start_time
                                     
                                     
data_out <- list(test_idx = test_idx,
                 n_simulations = n_simulations,
                 elapsed_time = elapsed_time, 
                 conformal_scores = conformal_score,
                 selected_parameters = selected_parameters)

save(data_out, file = paste0("data/test-info-nsim_",
                             n_simulations, "-", test_idx,
                             ".Rdata"))

