testthat::test_that("test update_tm_smooth",{
  smooth_function <- function(x,y){
    inner_ss <- smooth.spline(x,y, df = 5)
    return(predict(inner_ss,x)$y)
  }
  
  tm_radius <- list(list(min_cover_vec = c(1,1,.2,1,1,1.1),
                         dist_mat = matrix(1:36, nrow = 6)))
  
  tm_radius_update <- update_tm_smooth(tm_radius, smooth_function)
  
  testthat::expect_equal(tm_radius_update[[1]]$min_cover_vec,
                         diag(tm_radius_update[[1]]$dist_mat))
  testthat::expect_equal(tm_radius[[1]]$dist_mat[lower.tri(tm_radius[[1]]$dist_mat)],
                         tm_radius_update[[1]]$dist_mat[lower.tri(tm_radius_update[[1]]$dist_mat)])
  
  clean_storage <- tm_radius_update[[1]]$dist_mat
  clean_storage[lower.tri(clean_storage)] <- -Inf
  for (ncol in 2:ncol(clean_storage)){
    testthat::expect_true(all(clean_storage[,ncol-1] <= clean_storage[, ncol]))
  }
})




testthat::test_that("test create_discrete_uniform_cdf", {
  sillyf <- create_discrete_uniform_cdf(0,1)
  testthat::expect_equal(sillyf(c(-.1,0,.5,1,1.5, 20)),
                         c(0,.5,.5,1,1,1))
})



testthat::test_that("test coverage_down_list_save", {
  data_list <- list(data.frame(x = rep(1,10), y = 1:10),
                    data.frame(x = c(rep(1, each = 10)),
                               y = c(rep(c(.9, 10), each = 5))))
  
  dm_out <- coverage_down_list_save(data_list, e_cols = c("x","y"),
                                    .e_cols_string = F,
                                    verbose = F)
  testthat::expect_equal(dm_out[1,1,] ,rep(0, 10))
  testthat::expect_equal(dm_out[2,2,] ,rep(0, 10))
  
  testthat::expect_equal(dm_out[1,2,] ,c(.1,1.1,2.1,3.1,4.1,4,3,2,1,0))
  testthat::expect_equal(dm_out[2,1,] ,rep(c(.1,0), each = 5))
})


testthat::test_that("test coverage_down_slist_save", {
  data_list <- list(data.frame(x = rep(1,10), y = 1:10),
                    data.frame(x = c(rep(1, each = 10)),
                               y = c(rep(c(.9, 10), each = 5))),
                    data.frame(x = rep(1,10), y = 1:10+2))
  
  dm_out <- coverage_down_list_save(data_list, e_cols = c("x","y"),
                                    .e_cols_string = F,
                                    verbose = F)
  
  tm_radius <- coverage_down_slist_save(dm_out, g_order = c(1:3))
  
  testthat::expect_equal(tm_radius$min_cover_vec, c(0,4.1,2))
  testthat::expect_equal(diag(tm_radius$dist_mat), c(0,4.1,2))
  testthat::expect_equal(tm_radius$dist_mat[lower.tri(tm_radius$dist_mat)], 
                         as.numeric(rep(NA, 3)))
  testthat::expect_equal(tm_radius$dist_mat[upper.tri(tm_radius$dist_mat)], 
                         as.numeric(rep(4.1, 3)))
  
  tm_radius2 <- coverage_down_slist_save(dm_out, g_order = c(1,3,2))
  testthat::expect_equal(tm_radius2$min_cover_vec, c(0,2,2))
  testthat::expect_equal(diag(tm_radius2$dist_mat), c(0,2,2))
  testthat::expect_equal(tm_radius2$dist_mat[lower.tri(tm_radius2$dist_mat)], 
                         as.numeric(rep(NA, 3)))
  testthat::expect_equal(tm_radius2$dist_mat[upper.tri(tm_radius2$dist_mat)], 
                         as.numeric(rep(2, 3)))
  
  
  
  data_list2 <- list(data.frame(x = rep(1,10), y = 1:10),
                    data.frame(x = c(rep(1, each = 10)),
                               y = c(rep(c(.9, 10), each = 5))),
                    data.frame(x = rep(1,10), y = 1:10+2),
                    data.frame(x = rep(1.1,10), y = 1:10-1))
  
  dm_out2 <- coverage_down_list_save(data_list2, e_cols = c("x","y"),
                                    .e_cols_string = F,
                                    verbose = F)
  
  tm_radius3 <- coverage_down_slist_save(dm_out2, g_order = 1:4)
  
  testthat::expect_equal(tm_radius3$min_cover_vec, diag(tm_radius3$dist_mat))
  # TODO: probably should do more here...
  
  
  tm_radius4 <- coverage_down_slist_save(dm_out2, g_order = 1)
  
  testthat::expect_equal(tm_radius4$min_cover_vec, 0)
  testthat::expect_equal(tm_radius4$dist_mat, matrix(0))
})

