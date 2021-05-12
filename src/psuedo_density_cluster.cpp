#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

//' multiples a matrix by a scalar
//' 
//' @param m matrix
//' @param s scalar (double)
//' @return matrix who's entries are m's time s.
// [[Rcpp::export]]
NumericMatrix scalartimesmat(NumericMatrix m, double s){
  int nrow = m.nrow();
  int ncol = m.ncol();
  
  NumericMatrix out(nrow, ncol);
  
  for (int i = 0; i < nrow; ++i){
    for (int j = 0; j < ncol; ++j){
      out(i,j) = s*m(i,j);
    }
  }
  return(out);
}

//' adds two matrices together
//' 
//' @param m1 matrix
//' @param m2 matrix
//' @return matrix with elements that are pointwise addition of elements in
//' m1 and m2
// [[Rcpp::export]]
NumericMatrix addmats(NumericMatrix m1, NumericMatrix m2){
  int nrow = m1.nrow();
  int ncol = m1.ncol();
  
  if (ncol != m2.ncol()) {
    throw std::runtime_error("Incompatible number of dimensions (columns)");
  }
  
  if (nrow != m2.nrow()) {
    throw std::runtime_error("Incompatible number of dimensions (rows)");
  }
  
  NumericMatrix out(nrow, ncol);
  
  for (int i = 0; i < nrow; ++i){
    for (int j = 0; j < ncol; ++j){
      out(i,j) = m1(i,j) + m2(i,j);
    }
  }
  return(out);
}


//' calculates the l2 distance between a matrix and each matrix in a list
//' 
//' @param m1 primary matrix
//' @param l2 list of matrices of the same shape of m1
//' @param usefrac bool if we should calculate the distance relative to a scaling of 
//' 1/nrow(m1) (before taking the sqrt).
//' @return numerical vector of distance between m1 and each matrix in list l2
// [[Rcpp::export]]
NumericVector distvec(NumericMatrix m1, List l2, bool usefrac = false) {
  int nrow = m1.nrow();
  int ncol = m1.ncol();
  
  int n_list = l2.size();
  
  double frac;
  if (usefrac) {
    frac = 1.0/nrow;
  } else {
    frac = 1.0;
  }
  
  NumericVector out(n_list);
  
  for (int l = 0; l < n_list; ++l){
    double total = 0.0;
    NumericMatrix inner2 = l2[l];
    
    if (ncol != inner2.ncol()) {
      throw std::runtime_error("Incompatible number of dimensions (columns)");
    }
    
    if (nrow != inner2.nrow()) {
      throw std::runtime_error("Incompatible number of dimensions (rows)");
    }
    
    for (int i = 0; i < nrow; ++i){
      for (int j = 0; j < ncol; ++j){
        total += pow(m1(i,j) - inner2(i,j), 2.0);
      }
    }
    out[l] = sqrt(frac * total);
  }
  return out;
}

//' calculates the l2 distance two matrices
//' 
//' almost  same as EpiCompare::l2filamentdist...
//' 
//' @param m1 matrix
//' @param m2 matrix (same size as m2)
//' @param usefrac bool if we should calculate the distance relative to a scaling of 
//' 1/nrow(m1) (before taking the sqrt).
//' @return numerical distance between m1 and m2
// [[Rcpp::export]]
double difffunction(NumericMatrix m1, NumericMatrix m2, bool usefrac = false) {
  int nrow = m1.nrow();
  int ncol = m1.ncol();
  
  if (ncol != m2.ncol()) {
    throw std::runtime_error("Incompatible number of dimensions (columns)");
  }
  
  if (nrow != m2.nrow()) {
    throw std::runtime_error("Incompatible number of dimensions (rows)");
  }
  
  double frac;
  if (usefrac) {
    frac = 1.0/nrow;
  } else {
    frac = 1.0;
  }
  
  double total = 0.0;
  
  for (int i = 0; i < nrow; ++i){
    for (int j = 0; j < ncol; ++j){
      total += pow(m1(i,j) - m2(i,j), 2.0);
    }
  }
  double out = sqrt(frac * total);
  
  return out;
}


// [[Rcpp::export]]
NumericVector inner_dist_density(NumericVector distance, double sigma){
  int vlength = distance.size();
  NumericVector psd(vlength);
  
  if (sigma <= 0) {
    throw std::runtime_error("sigma needs to be strictly positive.");
  }
  
  for (int v_idx = 0; v_idx < vlength; ++v_idx) {
    psd[v_idx] = 1.0/sqrt(2.0 * atan(1)*4.0) * exp(-1.0/2.0 * pow(distance[v_idx]/sigma,2));
  }
  return(psd);
}

//' psuedo density walking towards the mode
//' 
//' @param X_list list of matrices (same size), these define the psuedo density
//' @param G_list list of matrices to walk up (can be the same as X_list, but 
//' need not be).
//' @param sigma double, sigma for distance psuedo density calculation
//' @param eps double, distance between iterative steps when to stop progressing 
//' (aka stop at a mode)
//'@param maxT int, max number of interations
//'@param verbose boolean, if we should use a progress bar
//'@param usefrac boolean, if we should calculate the distance relative to this 
//'scaling (should also impact the earlier calculation of sigma)
//'
//'@return Updated step of G_list
// [[Rcpp::export]]
List psuedo_density_mode_cluster(List X_list, List G_list, double sigma,
                                    double eps = 1E-06, int maxT = 10,
                                    bool verbose = true, bool usefrac = false) {
  
  int n_listX = X_list.size();
  int n_listG = G_list.size();
  int t = 0;
  
  List G_list_out(n_listG);
  
  Progress p(maxT, verbose);
  
  
  NumericVector error (n_listG, 1E08);  // initial error = massive error

  double max_error = 1E08;
  
  while (max_error > eps && t < maxT){
    if (Progress::check_abort() )
      return -1.0;
    
    for (int g_index = 0; g_index < n_listG; ++g_index){
      if (error[g_index] > eps){
        
        NumericVector dist_vals(n_listX);
        if (t == 0) {
          dist_vals = inner_dist_density(distvec(G_list[g_index],X_list, usefrac), sigma);
        } else {
          dist_vals = inner_dist_density(distvec(G_list_out[g_index],X_list, usefrac), sigma);
        }
      
        double denominator = std::accumulate(dist_vals.begin(), dist_vals.end(), 0.0);
  
        if (denominator == 0) {
          throw std::runtime_error("Error: no distance between a observation in G_list and all observations in X_list > 0.");
        }
  
        NumericMatrix unscaled_new = scalartimesmat(X_list[0],
                                                    dist_vals[0]);
  
        for (int x_index = 1; x_index < n_listX; ++x_index){
          unscaled_new = addmats(unscaled_new,
                                 scalartimesmat(X_list[x_index],dist_vals[x_index]));
        }
        
        NumericMatrix tmp = scalartimesmat(unscaled_new, 1.0/denominator);
        if (t == 0){
          error[g_index] = difffunction(tmp,G_list[g_index]);
        } else {
          error[g_index] = difffunction(tmp,G_list_out[g_index]);
        }
        G_list_out[g_index] = tmp;
      }
    }
    p.increment(); 
    max_error = *std::max_element(error.begin(), error.end());
    t += 1;
  }
  return(G_list_out);
}



//' psuedo density walking towards the mode
//' 
//' @param X_list list of matrices (same size), these define the psuedo density
//' @param G_list list of matrices to walk up (can be the same as X_list, but 
//' need not be).
//' @param sigma double, sigma for distance psuedo density calculation
//' @param eps double, distance between iterative steps when to stop progressing 
//' (aka stop at a mode)
//'@param maxT int, max number of interations
//'@param verbose boolean, if we should use a progress bar
//'@param usefrac boolean, if we should calculate the distance relative to this 
//'scaling (should also impact the earlier calculation of sigma)
//'
//'@return Updated step of G_list
// [[Rcpp::export]]
List psuedo_density_mode_cluster2(List X_list, List G_list, double sigma,
                                 double eps = 1E-06, int maxT = 10,
                                 bool verbose = true, bool usefrac = false) {
  
  int n_listX = X_list.size();
  int n_listG = G_list.size();
  int t = 0;
  
  List G_list_out(n_listG);
  
  Progress p(maxT, verbose);
  
  
  NumericVector error (n_listG, 1E08);  // initial error = massive error
  
  double max_error = 1E08;
  
  while (max_error > eps && t < maxT){
    if (Progress::check_abort() )
      return -1.0;
    
    for (int g_index = 0; g_index < n_listG; ++g_index){
      if (error[g_index] > eps){
        
        NumericVector dist_vals(n_listX);
        if (t == 0) {
          dist_vals = inner_dist_density(distvec(G_list[g_index],X_list, usefrac), sigma);
        } else {
          dist_vals = inner_dist_density(distvec(G_list_out[g_index],X_list, usefrac), sigma);
        }
        
        double denominator = std::accumulate(dist_vals.begin(), dist_vals.end(), 0.0);
        
        if (denominator == 0) {
          throw std::runtime_error("Error: no distance between a observation in G_list and all observations in X_list > 0.");
        }
        
        NumericMatrix unscaled_new = scalartimesmat(X_list[0],
                                                    dist_vals[0]);
        
        for (int x_index = 1; x_index < n_listX; ++x_index){
          unscaled_new = addmats(unscaled_new,
                                 scalartimesmat(X_list[x_index],dist_vals[x_index]));
        }
        
        NumericMatrix tmp = scalartimesmat(unscaled_new, 1.0/denominator);
        if (t == 0){
          error[g_index] = difffunction(tmp,G_list[g_index]);
        } else {
          error[g_index] = difffunction(tmp,G_list_out[g_index]);
        }
        G_list_out[g_index] = tmp;
      }
    }
    p.increment(); 
    max_error = *std::max_element(error.begin(), error.end());
    t += 1;
  }
  
  
  List out(3);
  if (t < maxT){
    out[0] = false;
    out[1] = t;
    out[2] = G_list_out;
    return(out);
  } else {
    out[0] = true;
    out[1] = t;
    out[2] = G_list_out;
    return(out);
  }
}