testthat::test_that("test psuedo_density_mode_cluster", {
  # basic
  X_list <- list(matrix(1:12, nrow = 3), matrix(rnorm(12), nrow = 3))
  out <- psuedo_density_mode_cluster(X_list = X_list,
                                     G_list = X_list, 
                                     sigma = 1, maxT = 1, verbose = F)
  testthat::expect_equal(dim(X_list[[1]]), dim(out[[1]]))
  testthat::expect_equal(dim(X_list[[1]]), dim(out[[2]]))
  
  testthat::expect_equal(X_list[[1]], matrix(1:12, nrow = 3)) # shouldn't override current matrix values

  # check same as R version
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  all_curves <- list(curve1,curve2, curve3, curve4, curve5) %>%
    lapply(function(df) df %>% select(-id))
  
  
  
  dist_mat <- dist_matrix_innersq_direction(all_curves, position = 1:2)
  sigma <- stats::quantile(dist_mat, .4)
  
  G_out <- functional_psuedo_density_mode_cluster_r(X_list = all_curves,
                                                     G_list = all_curves,
                                                     sigma = sigma,maxT = 30,
                                                    verbose = F)
  
  G_out2 <- psuedo_density_mode_cluster(X_list = lapply(all_curves, as.matrix),
                                        G_list = lapply(all_curves, as.matrix),
                                        sigma = sigma, maxT = 30,
                                        verbose = F)
  
  if (F){
    a <- lapply(1:length(G_out2), function(i) data.frame(G_out2[[i]]) %>% mutate(id = i)) %>% 
      do.call(rbind,.) %>% 
      ggplot() +
      geom_path(aes(x=X1,y =X2, color = factor(id))) +
      geom_path(data = lapply(1:length(all_curves), function(i) data.frame(all_curves[[i]]) %>% mutate(id = i)) %>% 
                  do.call(rbind,.),
                aes(x=x,y =y, color = factor(id)), alpha = .1)
    b <- lapply(1:length(G_out), function(i) data.frame(G_out[[i]]) %>% mutate(id = i)) %>% 
      do.call(rbind,.) %>% 
      ggplot() +
      geom_path(aes(x=x,y =y, color = factor(id))) +
      geom_path(data = lapply(1:length(all_curves), function(i) data.frame(all_curves[[i]]) %>% mutate(id = i)) %>% 
                  do.call(rbind,.),
                aes(x=x,y =y, color = factor(id)), alpha = .1)
    
    grid.arrange(a,b,nrow = 1)
  }
  
  testthat::expect_true(max(dist_matrix_innersq_2d(G_out2,position = 1:2)) 
                        < 1e-10) # pretty much converged to 1 point
  
  all_curves_plus <- all_curves
  all_curves_plus[[6]] <- as.data.frame(G_out2[[1]])
  
  pd_vs_original <- sapply(all_curves_plus, distance_kernel_function_r, 
                           df_list = all_curves, sigma = sigma) %>% 
    apply(2, mean)
  testthat::expect_true(all(pd_vs_original[-6] < pd_vs_original[6]))
  
  # speed check 
  mb_out <- microbenchmark::microbenchmark(
    rcpp = psuedo_density_mode_cluster(X_list = lapply(all_curves, as.matrix),
                                       G_list = lapply(all_curves, as.matrix), sigma = sigma,
                                       maxT = 30, verbose = F),
    rwrapper_cpp = functional_psuedo_density_mode_cluster(X_list = lapply(all_curves, as.matrix),
                                                          G_list = lapply(all_curves, as.matrix), sigma = sigma,
                                                          maxT = 30, verbose = F),
    rsmart = functional_psuedo_density_mode_cluster_r(X_list = all_curves,
                                                      G_list = all_curves,
                                                      sigma = sigma, verbose = F,
                                                      maxT = 30),times = 10)
  median_info <- summary(mb_out)$median
  testthat::expect_lt(median_info[1], median_info[2])
  testthat::expect_lt(median_info[2], median_info[3])
})

testthat::test_that("test scalartimesmat", {
  rcpp_value <- scalartimesmat(matrix(1:12, ncol = 3), 2.5)
  testthat::expect_equal(2.5*matrix(1:12, ncol = 3),rcpp_value )
})

testthat::test_that("test addmats", {
  rcpp_value <- addmats(matrix(1:12, nrow = 3), matrix(2, nrow = 3, ncol = 4))
  testthat::expect_equal(matrix(1:12, nrow = 3) + 
                           matrix(2, nrow = 3, ncol = 4),
                         rcpp_value )
})

testthat::test_that("test difffunction", {
  # basics
  df1 <- matrix(1:12, nrow = 3)
  df2 <- matrix(rnorm(12), nrow = 3)
  rcpp_out <- difffunction(df1, df2)
  
  testthat::expect_equal(sqrt(sum((df1-df2)^2)),
                         rcpp_out)
  
  # speed tests
  df1 <- EpiCompare::simulate_SIR_agents(n_sims = 1,n_time_steps = 100,
                                         beta = .1, gamma = .1/2,init_SIR = c(900,100,0)) %>%
    EpiCompare::agents_to_aggregate(states = c(tI, tR)) %>%
    dplyr::select(-t)
  
  df2 <- EpiCompare::simulate_SIR_agents(n_sims = 1,n_time_steps = 100,
                                         beta = .1, gamma = .1/2,init_SIR = c(900,100,0)) %>%
    EpiCompare::agents_to_aggregate(states = c(tI, tR)) %>%
    dplyr::select(-t)
  
  
  mb_out <- microbenchmark::microbenchmark(
    rcpp_complete = distvec(as.matrix(df1), list(as.matrix(df2), as.matrix(df1))),
    rcpp_and_r = sapply(list(as.matrix(df2), as.matrix(df1)), function(df) difffunction(as.matrix(df1), df)),
    r_complete = sapply(list(df2, df1), function(df) diff_function_r(df1, df))
  )
  median_info <- summary(mb_out)$median
  testthat::expect_lt(median_info[1], median_info[2])
  testthat::expect_lt(median_info[2], median_info[3])
})


testthat::test_that("test distvec", {
  df1 <- matrix(1:12, nrow = 3)
  df2 <- matrix(rnorm(12), nrow = 3)
  rcpp_out <- distvec(df1, list(df1, df2))
  
  testthat::expect_equal(c(0,sqrt(sum((df1-df2)^2))),
                         rcpp_out)
})

testthat::test_that("test inner_dist_density", {
  dvec <- c(-1,0,1)
  sigma <- 1
  out <- inner_dist_density(dvec, sigma)
  testthat::expect_equal(out, dnorm(dvec))
})

testthat::test_that("test combo of distvec and inner_dist_density", {
  X_list <- lapply(1:4, function(i) matrix(rnorm(100), ncol = 4))
  sigma <- 1
  rversion <- distance_kernel_function_r(X_list[[1]], X_list, sigma)
  rcppversion <- inner_dist_density(distvec(X_list[[1]], X_list), sigma)
  
  testthat::expect_equal(rversion,rcppversion)
})

testthat::test_that("test mode_clustering, basic (no usefrac)", {
  # single grouping
  
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  all_curves <- list(curve1,curve2, curve3, curve4, curve5) %>%
    lapply(function(df) df %>% select(-id))
  
  
  
  dist_mat <- dist_matrix_innersq_direction(all_curves, position = 1:2)
  sigma <- stats::quantile(dist_mat, .45)
  
  m_df <- mode_clustering(g_list = all_curves, 
                          g_names = as.character(1:5),
                          position = 1:2, sigma = sigma,
                          usefrac = F,diff_eps = 1e-05,
                          maxT = 30, verbose = F)
  testthat::expect_equal(m_df$groupings, rep(1, 5))
  
  m_df2 <- mode_clustering(g_list = all_curves, 
                           g_names = data.frame(id = as.character(1:5),
                                                id2 = c(1,1,2,2,3)),
                           diff_eps = 1e-05, usefrac = F,
                           position = 1:2, sigma = sigma, maxT = 30, verbose = F)
  testthat::expect_equal(m_df2$groupings, rep(1, 5))
  
  # 2 groups
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  all_curves <- rbind(curve1, curve2, curve3, curve4, curve5)
  
  all_curves2 <- all_curves %>% 
    mutate(x = x + 7, y = y - 7,
           id = as.character(as.numeric(id) + 5))
  
  all_curves <- rbind(all_curves, all_curves2)
  
  
  if (FALSE){ # visualize
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  all_curves_list <- all_curves %>%
    mutate(id = as.numeric(id)) %>% group_split(id)
  dist_mat2 <- dist_matrix_innersq_direction(all_curves_list, position = 1:2)
  sigma2 <- stats::quantile(dist_mat2, .3)
  
  m_df2 <- mode_clustering(g_list = all_curves_list, 
                           g_names = as.character(1:10), 
                           position = 1:2, sigma = sigma2, 
                           diff_eps = 1e-05, usefrac = F,
                           maxT = 50, verbose = F)
  
  testthat::expect_equal(m_df2$groupings, rep(1:2, each = 5))
})


testthat::test_that("test coverage_down_list, basics", {
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  
  all_curves <- list(curve1,curve2, curve3, curve4, curve5) %>%
    lapply(function(df) df %>% select(-id))
  
  out <- coverage_down_list(all_curves, c(x,y))
  
  testthat::expect_equal(out$min_cover_vec, c(0,1,1,1,1))
  testthat::expect_equal(out$dist_mat[upper.tri(out$dist_mat, T)], 
                         c(0,rep(1, 14)))
  testthat::expect_true(all(is.na(out$dist_mat[lower.tri(out$dist_mat,F)])))
  
  all_curves2 <- list(curve5, curve1,curve2, curve3, curve4) %>%
    lapply(function(df) df %>% select(-id))
  
  out2 <- coverage_down_list(all_curves2, c(x,y))
  
  testthat::expect_equal(out2$min_cover_vec, c(0,.52,1,1,1))
  testthat::expect_equal(out2$dist_mat[upper.tri(out2$dist_mat, T)], 
                         c(0,rep(.52,2),rep(1, 12)))
  testthat::expect_true(all(is.na(out2$dist_mat[lower.tri(out2$dist_mat,F)])))
  
  
  all_curves3 <- list(curve1,curve3,curve5, curve2, curve4) %>%
    lapply(function(df) df %>% select(-id))
  
  out3 <- coverage_down_list(all_curves3, c(x,y))
  
  testthat::expect_equal(out3$min_cover_vec, c(0,1,.52,1,1))
  testthat::expect_equal(out3$dist_mat[upper.tri(out3$dist_mat, T)], 
                         c(0,rep(1,4),.52, rep(1, 9)))
  testthat::expect_true(all(is.na(out2$dist_mat[lower.tri(out3$dist_mat,F)])))
})


testthat::test_that("test inner_min_dist_per_point, basic", {
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  all_curves <- list(curve1,curve2, curve3, curve4, curve5) %>%
    lapply(function(df) df %>% select(-id))
  
  
  out1 <- inner_min_dist_per_point(all_curves[[1]], all_curves[2:5], c(x,y))
  
  testthat::expect_equal(nrow(out1$nn_info_l),4)
  testthat::expect_equal(ncol(out1$nn_info_l),50)
  testthat::expect_equal(length(out1$nn_info_current),50)
  
  out2 <- inner_min_dist_per_point(all_curves[[1]], all_curves[1:5], c(x,y))
  testthat::expect_equal(nrow(out2$nn_info_l),5)
  testthat::expect_equal(ncol(out2$nn_info_l),50)
  testthat::expect_equal(length(out2$nn_info_current),50)
  
  testthat::expect_equal(out2$nn_info_current, rep(0,50))
  testthat::expect_equal(out2$nn_info_l[1,], rep(0,50))
  
})

testthat::test_that("test inner_containment_conformal_score_mode_radius", {
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  all_curves <- rbind(curve1, curve2, curve3)
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  all_curves <- rbind(curve1, curve2, curve3, curve4, curve5)
  
  all_curves2 <- all_curves %>% 
    mutate(x = x + 7, y = y - 7,
           id = as.character(as.numeric(id) + 5))
  
  all_curves <- rbind(all_curves, all_curves2)
  
  
  if (FALSE){ # visualize
    #all_curves <- rbind(curve1, curve2, curve3, curve4, curve5)
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  sim_curves <- all_curves %>% filter(!(id %in% c("8","2"))) %>% 
    group_by(id)
  new_curves <- all_curves %>% filter(id %in% c("8","2")) %>% 
    group_by(id)
  
  tdmat <- dist_matrix_innersq_direction.grouped_df(sim_curves,position = 1:2,
                                                    verbose = T,tdm_out = T) 
  sigma <- as.matrix(tdmat) %>% quantile(.35)
  
  ps_d <- distance_psuedo_density_function(tdmat,sigma = sigma,df_out = T) %>%
    mutate(id = as.numeric(id))
  sim_curves_list <- sim_curves %>%
    mutate(id = as.numeric(id)) %>% group_split() 
  sim_curves_names <- sim_curves %>% 
    mutate(id = as.numeric(id)) %>% 
    group_keys()
  
  m_d <- mode_clustering(g_list = sim_curves_list, 
                         g_names = sim_curves_names, 
                         position = 1:2, diff_eps = 1e-05,usefrac = F,
                         sigma =sigma, maxT = 30,
                         verbose = F)
  
  order_group_info <- inner_expanding_info(psuedo_density_df = ps_d, 
                                           mode_grouping = m_d)
  
  simulation_info_df <- order_group_info[[1]] %>%
    mutate(id = as.character(id))
  order_listing <- order_group_info[[2]]
  
  list_radius_info <- coverage_down_mlist(data_ll = sim_curves_list, 
                                          e_cols = c(x,y),
                                          g_order_ll = order_listing,
                                          names_df = ps_d %>% select(-psuedo_density),
                                          verbose = F)
  
  out <- inner_containment_conformal_score_mode_radius(df_row_group = new_curves,
                                                       simulations_group_df = sim_curves, # need the list structure?
                                                       data_column_names = c("x", "y"),
                                                       simulation_info_df = simulation_info_df, # needs grouping info? group ranking and overall ranking? and psuedo_density
                                                       list_radius_info = list_radius_info,
                                                       list_grouping_id = order_listing,
                                                       verbose = F)
  
  testthat::expect_equal(out %>% ungroup(), 
                         tibble(id = as.character(c(2,8)), 
                                containment_val = c(7,4)))
  
  # Fixed radius
  tdm_sims <- dist_matrix_innersq_direction(sim_curves,
                                            position = 1:2,
                                            tdm_out = T)
  
  sigma_size <- c("25%" = .25, "30%" = .3,
                  "35%" = .35, "40%" = .4, "45%" = .45)
  
  percentage_str <- names(sigma_size)[stats::quantile(as.matrix(tdm_sims), sigma_size) > 0][1]
  
  percentage_inner <- check_character_percent(percentage_str, "sigma")
  sigma_val <- stats::quantile(as.matrix(tdm_sims), percentage_inner)
  
  
  ps_d <- distance_psuedo_density_function(tdmat,sigma = sigma_val,df_out = T) %>%
    mutate(id = as.numeric(id))
  sim_curves_list <- sim_curves %>%
    mutate(id = as.numeric(id)) %>% group_split() 
  sim_curves_names <- sim_curves %>% 
    mutate(id = as.numeric(id)) %>% 
    group_keys()
  
  m_d <- mode_clustering(g_list = sim_curves_list, 
                         g_names = sim_curves_names, 
                         position = 1:2, diff_eps = 1e-05, usefrac = F,
                         sigma =sigma_val, maxT = 50,
                         verbose = F)
  

  
  top_points <- top_curves_to_points(sim_curves,
                                     tidy_dm = tdm_sims,
                                     alpha = .2,
                                     quantile_func = distance_psuedo_density_function,
                                     sigma = sigma_val) # 80% curve remain
  
  mm_delta <- get_delta_nn(top_points[,1:2])

  
  for (. in 1:5){
    order_group_info <- inner_expanding_info(psuedo_density_df = ps_d, 
                                             mode_grouping = m_d)
    
    simulation_info_df <- order_group_info[[1]] %>%
      mutate(id = as.character(id))
    order_listing <- order_group_info[[2]]
    list_radius_info_fake <- inner_convert_single_radius_to_structure(mm_delta, 
                                                                      order_listing)
    
    out_fixed_r <- inner_containment_conformal_score_mode_radius(df_row_group = new_curves,
                                                         simulations_group_df = sim_curves, # need the list structure?
                                                         data_column_names = c("x", "y"),
                                                         simulation_info_df = simulation_info_df, # needs grouping info? group ranking and overall ranking? and psuedo_density
                                                         list_radius_info = list_radius_info_fake,
                                                         list_grouping_id = order_listing,
                                                         verbose = F)
    testthat::expect_equal(out_fixed_r %>% ungroup(),
                           tibble(id = as.character(c(2,8)),
                                  containment_val = c(0,3)))
  }
  
})


testthat::test_that("test simulation_based_conformal3",{
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  all_curves <- rbind(curve1, curve2, curve3)
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  if (FALSE){ # visualize
    all_curves <- rbind(curve1, curve2, curve3, curve4)
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  sim_curves <- rbind(curve1, curve2, curve3) %>%
    group_by(id)
  truth_curve <- curve4 %>% 
    group_by(id)
  
  
  # same radius, no mode clustering 
  cs4 <- simulation_based_conformal3(truth_grouped_df = truth_curve,
                                     simulations_grouped_df = sim_curves,
                                     data_column_names = c("x", "y"),
                                     .change_radius = F,
                                     .mode_cluster = F,
                                     number_points = Inf,
                                     .to_simplex = F, 
                                     verbose = F)
  testthat::expect_equal(cs4[[1]]$containment_val, 0)
  
  for (. in 1:5){
    cs1 <- simulation_based_conformal3(truth_grouped_df = curve1 %>% group_by(id),
                                       simulations_grouped_df = sim_curves,
                                       data_column_names = c("x", "y"),
                                       .change_radius = F,
                                       .mode_cluster = F,
                                       number_points = Inf,
                                       .to_simplex = F, 
                                       verbose = F)
    testthat::expect_true(cs1[[1]]$containment_val %in% 3) # shouldn't be effected by random ordering...
  }
  
  curve3.5 <- curve1 %>%
    mutate(x = x - 1.5 * sqrt(2)/2,
           y = y + 1.5 * sqrt(2)/2,
           id = "3.5")
  
  if (F){
    sim_curves %>% 
      rbind(curve3.5) %>% 
      ggplot(aes(x = x,y=y,color = id)) +
      geom_line()
  }  
  
  for (i in 1:5){
    cs1.5 <- simulation_based_conformal3(truth_grouped_df = curve3.5 %>% group_by(id),
                                         simulations_grouped_df = sim_curves,
                                         data_column_names = c("x", "y"),
                                         .change_radius = F,
                                         .mode_cluster = F,
                                         number_points = Inf,
                                         .to_simplex = F, 
                                         verbose = F)
    testthat::expect_true(cs1.5[[1]]$containment_val %in% c(2,1))
  }
  
  # different radius, no mode clustering 
  
  cs4 <- simulation_based_conformal3(truth_grouped_df = truth_curve,
                                     simulations_grouped_df = sim_curves,
                                     data_column_names = c("x", "y"),
                                     .change_radius = F,
                                     .mode_cluster = F,
                                     number_points = Inf,
                                     .to_simplex = F, 
                                     verbose = F)
  testthat::expect_equal(cs4[[1]]$containment_val, 0)
  
  for (. in 1:5){
    cs1 <- simulation_based_conformal3(truth_grouped_df = curve1 %>% group_by(id),
                                       simulations_grouped_df = sim_curves,
                                       data_column_names = c("x", "y"),
                                       .change_radius = F,
                                       .mode_cluster = F,
                                       number_points = Inf,
                                       .to_simplex = F, 
                                       verbose = F)
    testthat::expect_true(cs1[[1]]$containment_val %in% 3) 
  }
  
  curve3.5 <- curve1 %>%
    mutate(x = x - 1.5 * sqrt(2)/2,
           y = y + 1.5 * sqrt(2)/2,
           id = "3.5")
  
  if (F){
    sim_curves %>% 
      rbind(curve3.5) %>% 
      ggplot(aes(x = x,y=y,color = id)) +
      geom_line()
  }  
  
  for (i in 1:5){
    cs1.5 <- simulation_based_conformal3(truth_grouped_df = curve3.5 %>% group_by(id),
                                         simulations_grouped_df = sim_curves,
                                         data_column_names = c("x", "y"),
                                         .change_radius = F,
                                         .mode_cluster = F,
                                         number_points = Inf,
                                         .to_simplex = F, 
                                         verbose = F)
    testthat::expect_true(cs1.5[[1]]$containment_val %in% c(2,1))
  }
  
  
  # mode clustering (currently only testing fixed radius)
  # 2 groups
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  all_curves <- rbind(curve1, curve2, curve3, curve4, curve5)
  
  all_curves2 <- all_curves %>% 
    mutate(x = x + 7, y = y - 7,
           id = as.character(as.numeric(id) + 5))
  
  all_curves <- rbind(all_curves, all_curves2)
  
  
  if (FALSE){ # visualize
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  all_curves_list <- all_curves %>%
    mutate(id = as.numeric(id)) %>% group_split(id)
  dist_mat2 <- dist_matrix_innersq_direction(all_curves_list, position = 1:2)
  sigma2 <- stats::quantile(dist_mat2, .1)
  
  m_df2 <- mode_clustering(g_list = all_curves_list, 
                           g_names = as.character(1:10),
                           position = 1:2, sigma = sigma2, 
                           maxT = 20, verbose = F)
  
  truth_grouped_df <- all_curves_list[(1:10 %in% c(2,8))] %>% 
    do.call(rbind,.) %>% group_by(id)
  simulations_grouped_df <- all_curves_list[!(1:10 %in% c(2,8))] %>% 
    do.call(rbind,.) %>% group_by(id)
  
  for (. in 1:5){
    out_mult_model <- simulation_based_conformal3(
      truth_grouped_df = truth_grouped_df,
      simulations_grouped_df = simulations_grouped_df,
      data_column_names = c("x", "y"),
      .change_radius = F,
      .mode_cluster = T, # mode cluster is true...
      number_points = Inf,
      .to_simplex = F, 
      verbose = F)
    
    if(F){
      all_curves_list %>% 
        do.call(rbind,.) %>%
        mutate(use = id %in% c(2,8)) %>%
        ggplot(aes(x=x,y=y,color = factor(id), alpha = use)) + 
        geom_line()
      # visual shows that 2 is probably never contained.
      # visual shows that 8 has cs score above 2 (as it should be included)
      # in regions not defined by the lines with large kinks in them.
      # 
      # because of the inclusion of observations in the simulation set it's 
      # possible that it's actually a value of 2 (so >= 1)...
    }
    
    cs_info <- out_mult_model$conformal_score
    
    testthat::expect_equal(cs_info$containment_val[cs_info$id == 2], 0)
    testthat::expect_gt(cs_info$containment_val[cs_info$id == 8], 2)
  }
})

testthat::test_that("test simulation_based_conformal3.5 with tests from 3",{
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  all_curves <- rbind(curve1, curve2, curve3)
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  if (FALSE){ # visualize
    all_curves <- rbind(curve1, curve2, curve3, curve4)
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  sim_curves <- rbind(curve1, curve2, curve3) %>%
    group_by(id)
  truth_curve <- curve4 %>% 
    group_by(id)
  
  
  # same radius, no mode clustering 
  cs4 <- simulation_based_conformal3.5(truth_grouped_df = truth_curve,
                                     simulations_grouped_df = sim_curves,
                                     data_column_names = c("x", "y"),
                                     .change_radius = F,
                                     .mode_cluster = F,
                                     number_points = Inf,
                                     .to_simplex = F, 
                                     verbose = F)
  testthat::expect_equal(cs4[[1]]$containment_val, 0)
  
  for (. in 1:5){
    cs1 <- simulation_based_conformal3.5(truth_grouped_df = curve1 %>% group_by(id),
                                       simulations_grouped_df = sim_curves,
                                       data_column_names = c("x", "y"),
                                       .change_radius = F,
                                       .mode_cluster = F,
                                       number_points = Inf,
                                       .to_simplex = F, 
                                       verbose = F)
    testthat::expect_true(cs1[[1]]$containment_val %in% 3) # shouldn't be effected by random ordering...
  }
  
  curve3.5 <- curve1 %>%
    mutate(x = x - 1.5 * sqrt(2)/2,
           y = y + 1.5 * sqrt(2)/2,
           id = "3.5")
  
  if (F){
    sim_curves %>% 
      rbind(curve3.5) %>% 
      ggplot(aes(x = x,y=y,color = id)) +
      geom_line()
  }  
  
  for (i in 1:5){
    cs1.5 <- simulation_based_conformal3(truth_grouped_df = curve3.5 %>% group_by(id),
                                         simulations_grouped_df = sim_curves,
                                         data_column_names = c("x", "y"),
                                         .change_radius = F,
                                         .mode_cluster = F,
                                         number_points = Inf,
                                         .to_simplex = F, 
                                         verbose = F)
    testthat::expect_true(cs1.5[[1]]$containment_val %in% c(2,1))
  }
  
  # different radius, no mode clustering 
  
  cs4 <- simulation_based_conformal3.5(truth_grouped_df = truth_curve,
                                     simulations_grouped_df = sim_curves,
                                     data_column_names = c("x", "y"),
                                     .change_radius = F,
                                     .mode_cluster = F,
                                     number_points = Inf,
                                     .to_simplex = F, 
                                     verbose = F)
  testthat::expect_equal(cs4[[1]]$containment_val, 0)
  
  for (. in 1:5){
    cs1 <- simulation_based_conformal3.5(truth_grouped_df = curve1 %>% group_by(id),
                                       simulations_grouped_df = sim_curves,
                                       data_column_names = c("x", "y"),
                                       .change_radius = F,
                                       .mode_cluster = F,
                                       number_points = Inf,
                                       .to_simplex = F, 
                                       verbose = F)
    testthat::expect_true(cs1[[1]]$containment_val %in% 3) 
  }
  
  curve3.5 <- curve1 %>%
    mutate(x = x - 1.5 * sqrt(2)/2,
           y = y + 1.5 * sqrt(2)/2,
           id = "3.5")
  
  if (F){
    sim_curves %>% 
      rbind(curve3.5) %>% 
      ggplot(aes(x = x,y=y,color = id)) +
      geom_line()
  }  
  
  for (i in 1:5){
    cs1.5 <- simulation_based_conformal3.5(truth_grouped_df = curve3.5 %>% group_by(id),
                                         simulations_grouped_df = sim_curves,
                                         data_column_names = c("x", "y"),
                                         .change_radius = F,
                                         .mode_cluster = F,
                                         number_points = Inf,
                                         .to_simplex = F, 
                                         verbose = F)
    testthat::expect_true(cs1.5[[1]]$containment_val %in% c(2,1))
  }
  
  
  # mode clustering (currently only testing fixed radius)
  # 2 groups
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  all_curves <- rbind(curve1, curve2, curve3, curve4, curve5)
  
  all_curves2 <- all_curves %>% 
    mutate(x = x + 7, y = y - 7,
           id = as.character(as.numeric(id) + 5))
  
  all_curves <- rbind(all_curves, all_curves2)
  
  
  if (FALSE){ # visualize
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  all_curves_list <- all_curves %>%
    mutate(id = as.numeric(id)) %>% group_split(id)
  dist_mat2 <- dist_matrix_innersq_direction(all_curves_list, position = 1:2)
  sigma2 <- stats::quantile(dist_mat2, .1)
  
  m_df2 <- mode_clustering(g_list = all_curves_list, 
                           g_names = as.character(1:10),
                           position = 1:2, sigma = sigma2, 
                           maxT = 20, verbose = F)
  
  truth_grouped_df <- all_curves_list[(1:10 %in% c(2,8))] %>% 
    do.call(rbind,.) %>% group_by(id)
  simulations_grouped_df <- all_curves_list[!(1:10 %in% c(2,8))] %>% 
    do.call(rbind,.) %>% group_by(id)
  
  for (. in 1:5){
    out_mult_model <- simulation_based_conformal3.5(
      truth_grouped_df = truth_grouped_df,
      simulations_grouped_df = simulations_grouped_df,
      data_column_names = c("x", "y"),
      .change_radius = F,
      .mode_cluster = T, # mode cluster is true...
      number_points = Inf,
      .to_simplex = F, 
      verbose = F)
    
    if(F){
      all_curves_list %>% 
        do.call(rbind,.) %>%
        mutate(use = id %in% c(2,8)) %>%
        ggplot(aes(x=x,y=y,color = factor(id), alpha = use)) + 
        geom_line()
      # visual shows that 2 is probably never contained.
      # visual shows that 8 has cs score above 2 (as it should be included)
      # in regions not defined by the lines with large kinks in them.
      # 
      # because of the inclusion of observations in the simulation set it's 
      # possible that it's actually a value of 2 (so >= 1)...
    }
    
    cs_info <- out_mult_model$conformal_score
    
    testthat::expect_equal(cs_info$containment_val[cs_info$id == 2], 0)
    testthat::expect_gt(cs_info$containment_val[cs_info$id == 8], 2)
  }
})


testthat::test_that("test simulation_based_conformal3.5 with tests with compression and .use_frac",{
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  all_curves <- rbind(curve1, curve2, curve3)
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  if (FALSE){ # visualize
    all_curves <- rbind(curve1, curve2, curve3, curve4)
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  sim_curves <- rbind(curve1, curve2, curve3) %>%
    group_by(id)
  truth_curve <- curve4 %>% 
    group_by(id)
  
  

 
  
  # mode clustering (currently only testing fixed radius)
  # 2 groups
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  curve5 <- curve2 %>%
    mutate(x = x - 1.52 * sqrt(2)/2,
           y = y + 1.52 * sqrt(2)/2,
           id = "5")
  
  all_curves <- rbind(curve1, curve2, curve3, curve4, curve5)
  
  all_curves2 <- all_curves %>% 
    mutate(x = x + 7, y = y - 7,
           id = as.character(as.numeric(id) + 5))
  
  all_curves <- rbind(all_curves, all_curves2)
  
  
  if (FALSE){ # visualize
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  all_curves_list <- all_curves %>%
    mutate(id = as.numeric(id)) %>% group_split(id)
  dist_mat2 <- dist_matrix_innersq_direction(all_curves_list, position = 1:2)
  sigma2 <- stats::quantile(dist_mat2, .1)
  
  m_df2 <- mode_clustering(g_list = all_curves_list, 
                           g_names = as.character(1:10),
                           position = 1:2, sigma = sigma2, 
                           maxT = 20, verbose = F)
  
  truth_grouped_df <- all_curves_list[(1:10 %in% c(2,8))] %>% 
    do.call(rbind,.) %>% group_by(id)
  simulations_grouped_df <- all_curves_list[!(1:10 %in% c(2,8))] %>% 
    do.call(rbind,.) %>% group_by(id)
  
  for (. in 1:5){
    out_mult_model <- simulation_based_conformal3.5(
      truth_grouped_df = truth_grouped_df,
      simulations_grouped_df = simulations_grouped_df,
      data_column_names = c("x", "y"),
      number_points = Inf,
      .change_radius = F,
      .mode_cluster = T, # mode cluster is true...
      .small_size_mode_cluster = 25, # smaller size
      .use_frac = F,
      .to_simplex = F, 
      verbose = F)
    
    if(F){
      all_curves_list %>% 
        do.call(rbind,.) %>%
        mutate(use = id %in% c(2,8)) %>%
        ggplot(aes(x=x,y=y,color = factor(id), alpha = use)) + 
        geom_line()
      # visual shows that 2 is probably never contained.
      # visual shows that 8 has cs score above 2 (as it should be included)
      # in regions not defined by the lines with large kinks in them.
      # 
      # because of the inclusion of observations in the simulation set it's 
      # possible that it's actually a value of 2 (so >= 1)...
    }
    
    cs_info <- out_mult_model$conformal_score
    
    testthat::expect_equal(cs_info$containment_val[cs_info$id == 2], 0)
    testthat::expect_gt(cs_info$containment_val[cs_info$id == 8], 2)
  }
})

testthat::test_that("test compression_and_sigma_estimate", {
  curve1 <- data.frame(x = (1:50)/2,
                       y = (1:50)/2,
                       id = "1")
  curve2 <- curve1 %>%
    mutate(x = x + sqrt(2)/2,
           y = y - sqrt(2)/2,
           id = "2")
  curve3 <- curve1 %>%
    mutate(x = x - sqrt(2)/2,
           y = y + sqrt(2)/2,
           id = "3")
  all_curves <- rbind(curve1, curve2, curve3)
  
  curve4 <- curve1 %>% mutate(id = "4")
  curve4$index <- curve4$x > 12.5
  curve4$x <- curve4$x + sqrt(2) * c(-1,1)[curve4$index+1]
  curve4$y <- curve4$y + sqrt(2) * c(1,-1)[curve4$index+1]
  curve4 <- curve4 %>% select(-index)
  
  if (FALSE){ # visualize
    all_curves <- rbind(curve1, curve2, curve3, curve4)
    all_curves %>% ggplot() +
      geom_line(aes(x = x , y = y, color = id))
  }
  
  sim_curves <- rbind(curve1, curve2, curve3) %>%
    group_by(id)
  
  # .use_frac is False
  info_out <- compression_and_sigma_estimate(sim_grouped_df = sim_curves,
                                 data_columns = c("x","y"),
                                 usefrac = F,
                                 number_points = 25,
                                 verbose = F)
  testthat::expect_equal(dim(info_out$compression), c(75,3)) # compression done
  testthat::expect_equal(names(info_out), c("compression", "sigma"))
  
  # .use_frac is True
  info_out2 <- compression_and_sigma_estimate(sim_grouped_df = sim_curves,
                                             data_columns = c("x","y"),
                                             usefrac = T,
                                             number_points = 30,
                                             verbose = F)
  testthat::expect_equal(dim(info_out2$compression), c(90,3)) # compression done
  testthat::expect_equal(names(info_out2), c("compression", "sigma"))
  testthat::expect_lt(info_out2$sigma, info_out$sigma) # smaller radius (even with more points...)
})


testthat::test_that("test get_delta_nn", {
  df_big <- data.frame(x = rnorm(5000),
                       y = rnorm(5000))
  mm_delta_nn <- get_delta_nn(df_big)
  get_delta_simple <- function(dist_mat){
    assertthat::assert_that(assertthat::are_equal(t(dist_mat), dist_mat),
                            msg = "t(dist_mat) needs to equal dist_mat")
    
    diag(dist_mat) <- max(dist_mat) # replacing diag so that it's selected as min
    mm_delta <- apply(dist_mat, MARGIN = 1, min) %>% max
    return(mm_delta)
  }
  
  mm_delta_simple <- get_delta_simple(as.matrix(dist(df_big)))
  
  testthat::expect_equal(mm_delta_nn, mm_delta_simple)
  
  d <- data.frame(x = 1:5)
  testthat::expect_equal(get_delta_nn(d),1)
  
  d2 <- data.frame(x = c(1,3:5))
  testthat::expect_equal(get_delta_nn(d2),2)
})


