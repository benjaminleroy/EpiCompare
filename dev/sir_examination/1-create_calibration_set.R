# global parameters -----------------------
number_points <- 150

# library loading -------------------------

library(tidyverse)
library(devtools)

# load 
if (getwd() != "/Users/benjaminleroy/Documents/CMU/research/EpiCompare/dev/sir_examination") {
  setwd("/Users/benjaminleroy/Documents/CMU/research/EpiCompare/dev/sir_examination")
}

load_all("../../") # EpiCompare@ben_dev_and_thesis
load_all("../../../simulationBands")

# calibration data generation -----------------

set.seed(5)
calibration_set_r0 <- tibble(x = runif(300))
calibration_set_r0$R0 <- sapply(calibration_set_r0$x, function(x){
  1/7*simulationBands::lei_wasserman_data_conditional_simulate(2*(x-.5),
                                                               n = 1)[[1]]$sim+2})
calibration_set_r0 <- calibration_set_r0 %>%
  mutate(beta = .1,
         gamma = .1/R0,
         idx = 1:300)

calibration_paths <- calibration_set_r0 %>% 
  dplyr::rename(x_inner = "x") %>% 
  dplyr::group_by(idx, x_inner, R0, beta, gamma) %>%
  tidyr::nest() %>%
  dplyr::mutate(sim = purrr::pmap(list(beta, gamma),
                           function(b,g){simulate_SIR_agents(n_sims = 1,
                                                             n_time_steps = 500,
                                                             b, g,
                                                             init_SIR = c(950,50,0)) %>%
                               agents_to_aggregate(states = c("tI", "tR")) %>%
                               as.data.frame() %>%
                               get_xy_coord(xyz_col = c("X0","X1","X2")) %>%
                               dplyr::mutate(idx = 1) %>%
                               dplyr::group_by(idx) %>%
                               filament_compression(data_columns = c("x","y"),
                                                    number_points = number_points) %>%
                               dplyr::ungroup() %>% dplyr::select(-idx)})) %>%
  dplyr::select(-data) %>%
  tidyr::unnest(sim) 

save(calibration_set_r0, calibration_paths, file = "data/calibration_paths.Rdata")


