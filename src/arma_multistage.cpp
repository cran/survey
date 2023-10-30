//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(rng=FALSE)]]
arma::mat arma_onestage(const arma::mat& Y,
                        const arma::colvec& samp_unit_ids,
                        const arma::colvec& strata_ids,
                        const arma::colvec& strata_samp_sizes,
                        const arma::colvec& strata_pop_sizes,
                        const Rcpp::CharacterVector& singleton_method,
                        const Rcpp::LogicalVector& use_singleton_method_for_domains,
                        const int& stage) { 
  
  size_t n_col_y = Y.n_cols;
  arma::uword number_of_data_rows = Y.n_rows;
  
  arma::mat result(n_col_y, n_col_y, arma::fill::zeros);
  
  if (number_of_data_rows == 0) {
    return result;
  }
  
  // Get distinct strata ids and their length, H
  arma::colvec distinct_strata_ids = unique(strata_ids);
  arma::uword H = distinct_strata_ids.n_elem;
  
  // Initialize vectors giving start and end rows for each stratum
  arma::uvec strata_ends;
  arma::uvec strata_starts;
  strata_ends.set_size(H);
  strata_starts.set_size(H);
  
  // Initialize checks for end of stratum
  bool at_end_of_stratum = true;
  
  // Prepare to iterate over list of strata
  arma::uword h = 0;
  
  // Get end index of each stratum
  arma::uvec row_indices = arma::linspace<arma::uvec>(0L, number_of_data_rows - 1L, number_of_data_rows);
  arma::uword next_row_index;
  for (arma::uvec::iterator row_index = row_indices.begin(); row_index != row_indices.end(); ++row_index) {
    
    next_row_index = (*(row_index+1));
    
    if ((*row_index) == (number_of_data_rows-1)) {
      at_end_of_stratum = true;
    } else {
      at_end_of_stratum = strata_ids(*row_index) != strata_ids(next_row_index);
    }
    
    if (at_end_of_stratum) {
      strata_ends(h) = *(row_index);
      h += 1L;
    }
  }
  
  // Get start index of each stratum
  if (H > 1) {
    h = 1;
    while (h < H) {
      strata_starts(h) = strata_ends((h-1)) + 1;
      h += 1;
    }
  }
  strata_starts(0) = 0;
  
  // Initialize strata summaries
  arma::mat cov_h(n_col_y, n_col_y, arma::fill::zeros);
  double n_h;
  double N_h;
  int n_obs_h;
  double f_h;
  double scale;
  
  // Initialize checks for end of sampling unit
  bool at_end_of_samp_unit = true;
  
  // Initialize sampling unit total
  arma::rowvec Yhi;
  Yhi = Yhi.zeros(n_col_y);
  
  // Get count of sampling units from all strata,
  // and check for singleton strata
  double n_samp_units_all_strata = 0;
  arma::rowvec singleton_indicators;
  singleton_indicators.zeros(H);
  for (arma::uword h = 0; h < H; ++h) {
    arma::uword stratum_start = strata_starts(h);
    n_samp_units_all_strata += strata_samp_sizes(stratum_start);
    // Check if there is only one sampling unit sampled in the stratum
    if (strata_samp_sizes(stratum_start) < 2) {
      // If only one sampling unit, check that it wasn't sampled with certainty
      if ((strata_pop_sizes(stratum_start) > 1)) {
        singleton_indicators(h) = 1;
        if ((singleton_method[0] == "fail")) {
          Rcpp::String error_msg("At least one stratum contains only one PSU at stage ");
          error_msg += stage;
          Rcpp::stop(error_msg);
        }
      }
      break;
    }
    
    if (use_singleton_method_for_domains[0] && ((singleton_method[0] == "adjust") || (singleton_method[0] == "average"))) {
      arma::uword h_row_index = strata_starts(h);
      // Set singleton indicator equal to 1 unless there are multiple observed sampling units in the stratum
      singleton_indicators(h) = 1;
      while ((singleton_indicators(h) == 1) && (h_row_index <= strata_ends(h))) {
        if (samp_unit_ids(h_row_index) != samp_unit_ids(strata_starts(h))) {
          singleton_indicators(h) = 0;
        }
        h_row_index += 1;
      }
    }
    
    if ((static_cast<int>(strata_ends[h] - strata_starts[h] + 1)) < 2) {
      if (use_singleton_method_for_domains[0] && ((singleton_method[0] == "adjust") || (singleton_method[0] == "average"))) {
        singleton_indicators(h) = 1;
      }
    }
  }
  int n_singleton_strata = sum(singleton_indicators);
  bool any_singleton_strata = n_singleton_strata > 0;
  
  arma::rowvec Y_means;
  Y_means = Y_means.zeros(n_col_y);
  if (any_singleton_strata & (singleton_method[0] == "adjust")) {
    Y_means = sum(Y, 0) / n_samp_units_all_strata;
  }
  
  // Iterate over each stratum
  for (arma::uword h = 0; h < H; ++h) {
    
    if ((singleton_indicators(h) == 1) && (singleton_method[0] != "adjust")) {
      break;
    }
    
    arma::uvec h_row_indices = arma::linspace<arma::uvec>(strata_starts[h], strata_ends[h], strata_ends[h] - strata_starts[h] + 1);
    
    // Get stratum-wide summary
    arma::uword stratum_start = strata_starts(h);
    
    n_h = static_cast<long double>(strata_samp_sizes(stratum_start));
    N_h = static_cast<double>(strata_pop_sizes(stratum_start));
    n_obs_h = 0;
    
    if (arma::is_finite(N_h)) {
      f_h = static_cast<double>(n_h) / static_cast<double>(N_h);
    } else {
      f_h = 0.0;
    }
    
    arma::rowvec Ybar_h;
    if ((singleton_indicators(h) == 1) && (singleton_method[0] == "adjust")) {
      Ybar_h = Y_means;
    } else {
      Ybar_h = arma::sum(Y.rows(h_row_indices), 0) / n_h;
    }
    
    // Iterate over each row in the stratum and its following row
    arma::uword h_next_row_index;
    for (arma::uvec::iterator h_row_index = h_row_indices.begin(); h_row_index != h_row_indices.end(); ++h_row_index) {
      
      h_next_row_index = (*(h_row_index+1));
      
      // Determine whether current row is last for the current sampling unit
      if ((*h_row_index) < strata_ends(h)) {
        at_end_of_samp_unit = samp_unit_ids(*h_row_index) != samp_unit_ids(h_next_row_index);
      } else {
        at_end_of_samp_unit = true;
      }
      
      // Get contribution to sampling unit\'s total
      Yhi += Y.row(*h_row_index);
      
      if (at_end_of_samp_unit) {
        n_obs_h += 1;
        Yhi.each_row() -= Ybar_h;
        cov_h += (arma::trans(Yhi)*Yhi);
        Yhi = Yhi.zeros();
      }
    }
    
    // If the data were subsetted, some sampling units
    // may not have rows of data that appear in inputs.
    // Make sure these units contribute to the variance.
    double n_h_missing = n_h - n_obs_h;
    if (n_h_missing > 0) {
      cov_h += n_h_missing*(arma::trans(Ybar_h)*Ybar_h);
    }
    
    // Determine scaling factor to use for normalizing sum of squares
    // and handling finite population correction
    scale = ((1.0 - f_h) * n_h);
    if (n_h > 1) {
      scale /= (n_h - 1);
    }
    
    result += (scale * cov_h);
    cov_h = cov_h.zeros();
  }
  
  
  // For the 'average' method for handling singleton strata,
  // the total variance contribution of each singleton stratum
  // is assumed to be equal to the variance contribution of the average nonsingleton stratum
  if ((singleton_method[0] == "average") && any_singleton_strata) {
    int n_nonsingleton_strata = H - n_singleton_strata;
    double scaling_factor;
    if (n_nonsingleton_strata > 0) {
      scaling_factor = static_cast<double>(H)/static_cast<double>(n_nonsingleton_strata);
    } else {
      scaling_factor = R_NaN;
    }
    result *= scaling_factor;
  }
  
  return result;
}

// [[Rcpp::export(rng=FALSE)]]
arma::mat arma_multistage(arma::mat Y,
                          arma::mat samp_unit_ids,
                          arma::mat strata_ids,
                          arma::mat strata_samp_sizes,
                          arma::mat strata_pop_sizes,
                          Rcpp::CharacterVector singleton_method,
                          Rcpp::LogicalVector use_singleton_method_for_domains,
                          Rcpp::LogicalVector use_only_first_stage,
                          int stage) {
  
  size_t n_stages = samp_unit_ids.n_cols;
  
  if (Y.n_rows == 0) {
    arma::mat V(Y.n_cols, Y.n_cols, arma::fill::zeros);
    return V;
  }
  
  // First reorder inputs by first-stage sample unit IDs
  arma::uvec samp_unit_id_order = arma::stable_sort_index(samp_unit_ids.col(0), "ascend");
  Y = Y.rows(samp_unit_id_order);
  samp_unit_ids = samp_unit_ids.rows(samp_unit_id_order);
  strata_ids = strata_ids.rows(samp_unit_id_order);
  strata_samp_sizes = strata_samp_sizes.rows(samp_unit_id_order);
  strata_pop_sizes = strata_pop_sizes.rows(samp_unit_id_order);
  
  // Next reorder inputs by first-stage strata IDs
  arma::uvec strata_id_order = arma::stable_sort_index(strata_ids.col(0), "ascend");
  Y = Y.rows(strata_id_order);
  samp_unit_ids = samp_unit_ids.rows(strata_id_order);
  strata_ids = strata_ids.rows(strata_id_order);
  strata_samp_sizes = strata_samp_sizes.rows(strata_id_order);
  strata_pop_sizes = strata_pop_sizes.rows(strata_id_order);
  
  // Obtain first stage information
  arma::colvec first_stage_ids = samp_unit_ids.col(0);
  arma::colvec first_stage_strata = strata_ids.col(0);
  arma::colvec first_stage_strata_samp_sizes = strata_samp_sizes.col(0);
  arma::colvec first_stage_strata_pop_sizes = strata_pop_sizes.col(0);
  
  // If there are later stages of sampling,
  // obtain the necessary columns from inputs,
  // which will be used recursively
  
  arma::mat later_stage_ids;
  arma::mat later_stage_strata;
  arma::mat later_stage_strata_samp_sizes;
  arma::mat later_stage_strata_pop_sizes;
  
  if ((n_stages > 1) && !use_only_first_stage[0]) {
    later_stage_ids = samp_unit_ids.tail_cols(n_stages - 1);
    later_stage_strata = strata_ids.tail_cols(n_stages - 1);
    later_stage_strata_samp_sizes = strata_samp_sizes.tail_cols(n_stages-1);
    later_stage_strata_pop_sizes = strata_pop_sizes.tail_cols(n_stages-1);
  }
  
  // Calculate first-stage variance
  arma::mat V = arma_onestage(Y,
                              first_stage_ids,
                              first_stage_strata,
                              first_stage_strata_samp_sizes,
                              first_stage_strata_pop_sizes,
                              singleton_method,
                              use_singleton_method_for_domains,
                              stage);
  
  // For each first-stage unit, get variance contribution from next stage
  if ((n_stages > 1) && !use_only_first_stage[0]) {
    
    // Get distinct first-stage strata ids and their length, H
    arma::colvec distinct_strata_ids = unique(first_stage_strata);
    arma::uword H = distinct_strata_ids.n_elem;
    
    for (arma::uword h = 0; h < H; ++h) {
      
      // Determine which rows of data correspond to the current first-stage stratum
      arma::uvec h_indices = arma::find(first_stage_strata==distinct_strata_ids(h));
      
      // Get submatrices of inputs corresponding to the current first-stage stratum
      arma::mat Y_h = Y.rows(h_indices);
      
      arma::mat h_samp_unit_ids = later_stage_ids.rows(h_indices);
      
      arma::mat h_strata = later_stage_strata.rows(h_indices);
      
      arma::mat h_strata_samp_sizes = later_stage_strata_samp_sizes.rows(h_indices);
      arma::mat h_strata_pop_sizes = later_stage_strata_pop_sizes.rows(h_indices);
      
      // Get count of first-stage sampling units in first-stage stratum
      // and finite population correction, based on all first-stage units in sample design
      arma::uword n_h = min(first_stage_strata_samp_sizes.elem(h_indices));
      double N_h = static_cast<double>(min(strata_pop_sizes.elem(h_indices)));
      
      double f_h;
      if (arma::is_finite(N_h)) {
        f_h = static_cast<double>(n_h) /  N_h;
      } else {
        f_h = 0.0;
        continue; // Skip to next stratum, since variance contribution is 0
      }
      
      // Get list of first-stage units in the current subset of data, and count them
      arma::colvec h_first_stage_units = first_stage_ids.elem(h_indices);
      arma::colvec h_unique_first_stage_units = unique(h_first_stage_units);
      arma::uword n_h_subset = h_unique_first_stage_units.n_elem;
      
      for (arma::uword i=0; i < n_h_subset; ++i ) {
        // Create subsets of inputs specific to current first-stage sampling unit
        arma::uvec unit_indices = arma::find(first_stage_ids.elem(h_indices) == h_unique_first_stage_units(i));
        arma::mat Y_hi = Y_h.rows(unit_indices);
        arma::mat hi_samp_unit_ids = h_samp_unit_ids.rows(unit_indices);
        arma::mat hi_strata = h_strata.rows(unit_indices);
        arma::mat hi_strata_samp_sizes = h_strata_samp_sizes.rows(unit_indices);
        arma::mat hi_strata_pop_sizes = h_strata_pop_sizes.rows(unit_indices);
        
        // Estimate later-stage variance contribution
        arma::mat V_hi = f_h * arma_multistage(Y_hi,
                                               hi_samp_unit_ids,
                                               hi_strata,
                                               hi_strata_samp_sizes,
                                               hi_strata_pop_sizes,
                                               singleton_method,
                                               use_singleton_method_for_domains,
                                               use_only_first_stage,
                                               stage = stage + 1);
        V += V_hi;
        
      }
    }
  }
  // Deal with potential floating-point error
  V.diag() = arma::abs(V.diag());
  for (arma::uword i=0; i < V.n_cols; ++i ) {
    if (std::abs(V(i,i)) < 1e-16) {
      V.row(i).zeros();
      V.col(i).zeros();
    }
  }
  
  // Return result
  return V;
}
