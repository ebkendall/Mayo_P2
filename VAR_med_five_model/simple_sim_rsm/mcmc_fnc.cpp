#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

NumericVector csample_num( NumericVector x, int size, bool replace,
                           NumericVector prob = NumericVector::create()) {
    NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
    return ret;
}

// Recursive functions to find all possible state sequences --------------------
void findStateSequences_backward(int currentState,
                                 int endState,
                                 const arma::imat& adjacencyMatrix,
                                 int sequenceLength,
                                 std::vector<std::vector<int>>& allSequences,
                                 std::vector<int>& currentSequence) {
    
    // Add current state to the current sequence
    currentSequence.push_back(currentState);
    
    // If current sequence length matches desired length and current state is endState, add to results
    if (currentSequence.size() == sequenceLength && currentState == endState) {
        allSequences.push_back(currentSequence);
    }
    
    // If current sequence length is less than desired length, continue exploring
    if (currentSequence.size() < sequenceLength) {
        // Find next states
        arma::uvec nextStates = find(adjacencyMatrix.row(currentState) == 1);
        
        // Recursively find sequences of desired length ending in endState starting from each next state
        for (size_t i = 0; i < nextStates.n_elem; ++i) {
            findStateSequences_backward(nextStates(i), endState, adjacencyMatrix, sequenceLength, allSequences, currentSequence);
        }
    }
    
    // Remove current state from the current sequence (backtracking)
    currentSequence.pop_back();
}

List getAllStateSequences_backward(int end_state, const arma::imat& adjacency_matrix, int sequence_length) {
    int nStates = adjacency_matrix.n_rows;
    
    // Check validity of end_state
    if (end_state < 0 || end_state >= nStates) {
        stop("Invalid ending state index.");
    }
    
    std::vector<std::vector<int>> allSequences;
    std::vector<int> currentSequence;
    
    // Find all sequences of specified length ending in end_state
    for (int start_state = 0; start_state < nStates; ++start_state) {
        findStateSequences_backward(start_state, end_state, adjacency_matrix, sequence_length, allSequences, currentSequence);
    }
    
    // // Convert C++ vectors to R list of integer vectors
    List result(allSequences.size());
    for (size_t i = 0; i < allSequences.size(); ++i) {
        result[i] = wrap(allSequences[i]);
    }
    
    return result;
}

void findStateSequences_forward(int currentState,
                                const arma::imat& adjacencyMatrix,
                                int sequenceLength,
                                std::vector<std::vector<int>>& allSequences,
                                std::vector<int>& currentSequence) {
    
    // Add current state to the current sequence
    currentSequence.push_back(currentState);
    
    // If current sequence length matches desired length, add to results
    if (currentSequence.size() == sequenceLength) {
        allSequences.push_back(currentSequence);
    }
    
    // If current sequence length is less than desired length, continue exploring
    if (currentSequence.size() < sequenceLength) {
        // Find next states
        arma::uvec nextStates = find(adjacencyMatrix.row(currentState) == 1);
        
        // Recursively find sequences of desired length starting from each next state
        for (size_t i = 0; i < nextStates.n_elem; ++i) {
            findStateSequences_forward(nextStates(i), adjacencyMatrix, sequenceLength, allSequences, currentSequence);
        }
    }
    
    // Remove current state from the current sequence (backtracking)
    currentSequence.pop_back();
}

List getAllStateSequences_forward(int start_state, const arma::imat& adjacency_matrix, int sequence_length) {
    int nStates = adjacency_matrix.n_rows;
    
    // Check validity of start_state
    if (start_state < 0 || start_state >= nStates) {
        stop("Invalid starting state index.");
    }
    
    std::vector<std::vector<int>> allSequences;
    std::vector<int> currentSequence;
    
    // Find all sequences of specified length starting from start_state
    findStateSequences_forward(start_state, adjacency_matrix, sequence_length, allSequences, currentSequence);
    
    // Convert C++ vectors to R list of integer vectors
    List result(allSequences.size());
    for (size_t i = 0; i < allSequences.size(); ++i) {
        result[i] = wrap(allSequences[i]);
    }
    
    return result;
}

void findStateSequences_both(int currentState,
                             int endState,
                             const arma::imat& adjacencyMatrix,
                             int sequenceLength,
                             std::vector<std::vector<int>>& allSequences,
                             std::vector<int>& currentSequence) {
    
    // Add current state to the current sequence
    currentSequence.push_back(currentState);
    
    // If current sequence length matches desired length and current state is endState, add to results
    if (currentSequence.size() == sequenceLength && currentState == endState) {
        allSequences.push_back(currentSequence);
    }
    
    // If current sequence length is less than desired length, continue exploring
    if (currentSequence.size() < sequenceLength) {
        // Find next states
        arma::uvec nextStates = find(adjacencyMatrix.row(currentState) == 1);
        
        // Recursively find sequences of desired length starting from each next state
        for (size_t i = 0; i < nextStates.n_elem; ++i) {
            findStateSequences_both(nextStates(i), endState, adjacencyMatrix, sequenceLength, allSequences, currentSequence);
        }
    }
    
    // Remove current state from the current sequence (backtracking)
    currentSequence.pop_back();
}

List getAllStateSequences_both(int start_state, int end_state, const arma::imat& adjacency_matrix, int sequence_length) {
    int nStates = adjacency_matrix.n_rows;
    
    // Check validity of start_state and end_state
    if (start_state < 0 || start_state >= nStates || end_state < 0 || end_state >= nStates) {
        stop("Invalid starting or ending state index.");
    }
    
    std::vector<std::vector<int>> allSequences;
    std::vector<int> currentSequence;
    
    // Find all sequences of specified length starting from start_state and ending in end_state
    findStateSequences_both(start_state, end_state, adjacency_matrix, sequence_length, allSequences, currentSequence);
    
    // Convert C++ vectors to R list of integer vectors
    List result(allSequences.size());
    for (size_t i = 0; i < allSequences.size(); ++i) {
        result[i] = wrap(allSequences[i]);
    }
    
    return result;
}

arma::field<arma::field<arma::mat>> get_Omega_list(const arma::imat &a_mat, int s) {
    
    int N = a_mat.n_cols; // dimension of adj matrix
    
    arma::field<arma::mat> c(N);
    arma::field<arma::mat> b(N, N);
    arma::field<arma::mat> a(N);
    
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            b(i, j) = arma::mat(1, s, arma::fill::zeros);
        }
    }
    
    for(int i = 0; i < N; i++) {
        // (a) -----------------------------------------------------------------
        Rcpp::List a_temp_list = getAllStateSequences_forward(i, a_mat, s+1);
        arma::mat a_i;
        if(a_temp_list.length() > 0) {
            arma::mat a_temp_mat(a_temp_list.length(), s+1);
            for(int k = 0; k < a_temp_list.length(); k++) {
                arma::rowvec a_i_temp = a_temp_list[k];
                a_temp_mat.row(k) = a_i_temp;
            }
            a_i = a_temp_mat.cols(1, s); // ignore the first column
        } else{
            arma::mat a_temp_mat(1, s, arma::fill::ones);
            a_temp_mat = -2 * a_temp_mat;
            a_i = a_temp_mat;
        }
        a(i) = a_i + 1;
        
        // (c) -----------------------------------------------------------------
        Rcpp::List c_temp_list = getAllStateSequences_backward(i, a_mat, s+1);
        arma::mat c_i;
        if(c_temp_list.length() > 0) {
            arma::mat c_temp_mat(c_temp_list.length(), s+1);
            for(int k = 0; k < c_temp_list.length(); k++) {
                arma::rowvec c_i_temp = c_temp_list[k];
                c_temp_mat.row(k) = c_i_temp;
            }
            c_i = c_temp_mat.cols(0, s-1); // ignore the last column
        } else{ 
            arma::mat c_temp_mat(1, s, arma::fill::ones);
            c_temp_mat = -2 * c_temp_mat;
            c_i = c_temp_mat;
        } 
        c(i) = c_i + 1; 
        
        // (b) -----------------------------------------------------------------
        for(int j = 0; j < N; j++) {
            Rcpp::List b_temp_list = getAllStateSequences_both(i, j, a_mat, s+2);
            arma::mat b_i;
            if(b_temp_list.length() > 0) {
                arma::mat b_temp_mat(b_temp_list.length(), s+2);
                for(int k = 0; k < b_temp_list.length(); k++) {
                    arma::rowvec b_i_temp = b_temp_list[k];
                    b_temp_mat.row(k) = b_i_temp;
                }
                b_i = b_temp_mat.cols(1, s); // ignore the first & last columns
            } else{ 
                arma::mat b_temp_mat(1, s, arma::fill::ones);
                b_temp_mat = -2 * b_temp_mat;
                b_i = b_temp_mat;
            } 
            b(i,j) = b_i + 1;
        }
        
    }
    
    arma::field<arma::field<arma::mat>> Omega_List(3);
    Omega_List(0) = c; Omega_List(1) = b; Omega_List(2) = a;
    
    return Omega_List;
}
// -----------------------------------------------------------------------------

arma::imat adj_mat_GLOBAL;
int states_per_step_GLOBAL;
arma::field<arma::field<arma::mat>> Omega_List_GLOBAL_multi;
arma::vec P_init;

// [[Rcpp::export]]
void initialize_cpp(arma::imat a_mat, int s_per_s) {
    adj_mat_GLOBAL = a_mat;
    states_per_step_GLOBAL = s_per_s;
    Omega_List_GLOBAL_multi = get_Omega_list(adj_mat_GLOBAL, states_per_step_GLOBAL);

    arma::vec init_global_temp(a_mat.n_rows, arma::fill::ones);
    P_init = init_global_temp / arma::accu(init_global_temp);
}

arma::mat get_omega_list(const int k, const int n_i, const arma::vec &b_i, int states_per_step) {
    
    arma::mat Omega_set;
    
    // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
    if (k == 0) {
        // () -> () -> 1-5
        Omega_set = Omega_List_GLOBAL_multi(0)(b_i(k + states_per_step) - 1);
    } else if (k == n_i - states_per_step) {
        // 1-5 -> () -> ()
        Omega_set = Omega_List_GLOBAL_multi(2)(b_i(n_i - states_per_step - 1) - 1);
    } else {
        // 1-5 -> () -> () -> 1-5
        Omega_set = Omega_List_GLOBAL_multi(1)(b_i(k - 1) - 1, 
                                            b_i(k + states_per_step) - 1);
    } 
    
    return Omega_set;
    
}

arma::mat Omega_fun_cpp_new_multi(const int k, const int n_i, const arma::vec &b_i,
                                  int t_pt_length) {
    
    arma::mat Omega_set;
    
    // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
    if (k == 1) {
        // () -> () -> 1-5
        Omega_set = Omega_List_GLOBAL_multi(0)(b_i(k + t_pt_length - 1) - 1);
    } else if (k <= n_i - t_pt_length) {
        // 1-5 -> () -> () -> 1-5
        Omega_set = Omega_List_GLOBAL_multi(1)(b_i(k - 2) - 1, 
                                            b_i(k + t_pt_length - 1) - 1);
    } else if (k == n_i - t_pt_length + 1) {
        // 1-5 -> () -> ()
        Omega_set = Omega_List_GLOBAL_multi(2)(b_i(n_i - t_pt_length - 1) - 1);
    }
    
    return Omega_set;
    
} 

// -----------------------------------------------------------------------------

double log_f_i_cpp_total(const arma::vec &EIDs, const arma::vec &par, 
                         const arma::field<arma::uvec> &par_index, 
                         const arma::field <arma::vec> &B, const arma::vec &y, 
                         const arma::vec &eids, int n_cores) {
    
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        double like_comp_transition = 0;
        double like_comp = 0;
        
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::vec b_i = B(ii);
        arma::vec y_i = y.rows(sub_ind);
        
        arma::vec ones(b_i.n_elem, arma::fill::ones);
        arma::vec twos(b_i.n_elem, arma::fill::zeros);
        twos.elem(arma::find(b_i == 2)) += 1;
        
        arma::mat D_b = arma::join_rows(ones, arma::cumsum(twos));
        arma::vec mean_b = D_b * mu;
        
        for(int jj = 0; jj < n_i; jj++) {
            
            // (1) Transition probabilities ------------------------------------
            if(jj == 0) {
                like_comp_transition = like_comp_transition + log(P_init(b_i(jj) - 1));
            } else{ 
                arma::vec qz = exp(zeta);
                arma::mat Q = { {    1,   qz(0)},
                                {qz(1),       1}};
                
                arma::vec q_row_sums = arma::sum(Q, 1);
                arma::mat P_i = Q.each_col() / q_row_sums;
                
                int curr_b_k_1 = b_i(jj-1);
                int curr_b_k   = b_i(jj);
                
                like_comp_transition = like_comp_transition + 
                    log(P_i(curr_b_k_1 - 1, curr_b_k - 1));
            } 
            
            // (2) Likelihood component ------------------------------------
            like_comp = like_comp + R::dnorm(y_i(jj), mean_b(jj), 1, true);
        } 
        
        in_vals(ii) = like_comp + like_comp_transition;
    } 
    
    double in_value = arma::accu(in_vals);
    
    return in_value;
}

// [[Rcpp::export]]
double log_post_cpp(const arma::vec &EIDs, const arma::vec &par, 
                    const arma::field<arma::uvec> &par_index,
                    const arma::field<arma::vec> &B, const arma::vec &y, 
                    const arma::vec &eids, int n_cores) {

  // Compute the likelihood ----------------------------------------------------
  double value = log_f_i_cpp_total(EIDs, par, par_index, B, y, eids, n_cores);
  
  // Prior densities
  arma::vec prior_mean(par.n_elem, arma::fill::zeros);
  
  arma::vec prior_var_diag(par.n_elem, arma::fill::ones);
  prior_var_diag = 20 * prior_var_diag;
  arma::mat prior_sd = arma::diagmat(prior_var_diag);
  
  arma::vec prior_dens = dmvnorm(par.t(), prior_mean, prior_sd, true);
  double prior_dens_val = arma::as_scalar(prior_dens);

  value = value + prior_dens_val;
  
  return value;
}

// [[Rcpp::export]]
arma::field<arma::vec> gibbs_up(const arma::vec EIDs, const arma::vec &par,
                                const arma::field<arma::uvec> &par_index,
                                arma::field <arma::vec> &B,
                                const arma::vec &y, const arma::vec &eids,
                                int n_cores, int states_per_step) {
    // par_index KEY: (0) mu, (1) zeta;
    // "i" is the numeric EID number; "ii" is the index of the EID
    arma::field<arma::vec> B_return(EIDs.n_elem);

    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);

    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};

    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        arma::vec b_i = B(ii);
        arma::vec y_i = y.rows(sub_ind);

        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - states_per_step + 1; k++) {

            // All possible state transitions given current time point -----
            arma::mat Omega_set = get_omega_list(k, n_i, b_i, states_per_step);
            
            if(Omega_set.n_rows > 1) {
                
                arma::vec ss_prob(Omega_set.n_rows, arma::fill::zeros);
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
                
                for(int jj = 0; jj < Omega_set.n_rows; jj++) {
                    
                    double log_like_val = 0;
                    arma::vec ss_jj = b_i;
                    ss_jj.rows(k, k + states_per_step - 1) = Omega_set.row(jj).t();
                    
                    arma::vec ones(b_i.n_elem, arma::fill::ones);
                    arma::vec twos_s(b_i.n_elem, arma::fill::zeros);
                    twos_s.elem(arma::find(ss_jj == 2)) += 1;
                    arma::mat D_s = arma::join_rows(ones, arma::cumsum(twos_s));
                    arma::vec mean_s = D_s * mu;
                    
                    for(int kk = k; kk < n_i; kk++) {
                        if(kk == 0) {
                            log_like_val = log_like_val + 
                                           log(P_init(ss_jj(kk) - 1)) + 
                                           R::dnorm(y_i(kk), mean_s(kk), 1, true);
                        } else if(kk <= n_i - states_per_step) {
                            log_like_val = log_like_val + 
                                           log(P_i(ss_jj(kk-1) - 1, ss_jj(kk) - 1)) + 
                                           R::dnorm(y_i(kk), mean_s(kk), 1, true);
                        }
                    }
                    
                    ss_prob(jj) = log_like_val;
                }
                
                double prob_log_max = ss_prob.max();
                ss_prob = ss_prob - prob_log_max;
                ss_prob = exp(ss_prob);
                ss_prob = (1/arma::accu(ss_prob)) * ss_prob;
                
                arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                b_i.rows(k, k+states_per_step-1) = Omega_set.row(row_ind(0)).t();
            }
        }
        
        B_return(ii) = b_i;
    }

    return B_return;
}

// [[Rcpp::export]]
arma::field<arma::vec> almost_gibbs_up(const arma::vec EIDs, const arma::vec &par,
                                       const arma::field<arma::uvec> &par_index,
                                       arma::field <arma::vec> &B,
                                       const arma::vec &y, const arma::vec &eids,
                                       int n_cores, int states_per_step) {
    // par_index KEY: (0) mu, (1) zeta;
    // "i" is the numeric EID number; "ii" is the index of the EID
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    
    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};
    
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::vec b_i = B(ii);
        arma::vec y_i = y.rows(sub_ind);

        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - states_per_step + 1; k++) {
        
            // All possible state transitions given current time point -----
            arma::mat Omega_set = get_omega_list(k, n_i, b_i, states_per_step);
            
            if(Omega_set.n_rows > 1) {
                
                // Learn state-sampling distribution via "Almost-Gibbs" --------
                arma::vec ss_prob(Omega_set.n_rows, arma::fill::zeros);
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
                
                for(int jj = 0; jj < Omega_set.n_rows; jj++) {
                    
                    double log_like_val = 0;
                    arma::vec ss_jj = b_i;
                    ss_jj.rows(k, k + states_per_step - 1) = Omega_set.row(jj).t();
                    
                    arma::vec ones(b_i.n_elem, arma::fill::ones);
                    arma::vec twos_s(b_i.n_elem, arma::fill::zeros);
                    twos_s.elem(arma::find(ss_jj == 2)) += 1;
                    arma::mat D_s = arma::join_rows(ones, arma::cumsum(twos_s));
                    arma::vec mean_s = D_s * mu;
                    
                    for(int kk = k; kk < k + states_per_step + 1; kk++) {
                        if(kk == 0) {
                            log_like_val = log_like_val + 
                                log(P_init(ss_jj(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_s(kk), 1, true);
                        } else if(kk <= n_i - states_per_step) {
                            log_like_val = log_like_val + 
                                log(P_i(ss_jj(kk-1) - 1, ss_jj(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_s(kk), 1, true);
                        }
                    }
                    
                    ss_prob(jj) = log_like_val;
                } 
                
                double prob_log_max = ss_prob.max();
                ss_prob = ss_prob - prob_log_max;
                ss_prob = exp(ss_prob);
                ss_prob = (1/arma::accu(ss_prob)) * ss_prob;
                
                
                // Compute the MH Ratio to accept or reject --------------------
                double log_like_b = 0;
                double log_like_s = 0;
                arma::vec s_i = b_i;
                
                arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                s_i.rows(k, k+states_per_step-1) = Omega_set.row(row_ind(0)).t();
                
                if(arma::approx_equal(s_i, b_i, "absdiff", 0.001)) {
                    b_i = s_i;
                } else {
                    arma::vec ones(b_i.n_elem, arma::fill::ones);
                    arma::vec twos_b(b_i.n_elem, arma::fill::zeros);
                    arma::vec twos_s(b_i.n_elem, arma::fill::zeros);
                    
                    twos_b.elem(arma::find(b_i == 2)) += 1;
                    twos_s.elem(arma::find(s_i == 2)) += 1;
                    
                    arma::mat D_b = arma::join_rows(ones, arma::cumsum(twos_b));
                    arma::mat D_s = arma::join_rows(ones, arma::cumsum(twos_s));
                    
                    arma::vec mean_b = D_b * mu;
                    arma::vec mean_s = D_s * mu;
                    
                    for(int kk = k; kk < n_i; kk++) {
                        if(kk == 0) {
                            log_like_b = log_like_b + 
                                log(P_init(b_i(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_b(kk), 1, true);
                            
                            log_like_s = log_like_s + 
                                log(P_init(s_i(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_s(kk), 1, true);
                            
                        } else if(kk <= n_i - states_per_step) {
                            log_like_b = log_like_b + 
                                log(P_i(b_i(kk-1) - 1, b_i(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_b(kk), 1, true);
                            
                            log_like_s = log_like_s + 
                                log(P_i(s_i(kk-1) - 1, s_i(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_s(kk), 1, true);
                        }
                    }
                    double diff_check = log_like_s - log_like_b;
                    double min_log = log(arma::randu(arma::distr_param(0,1)));
                    if(diff_check > min_log){b_i = s_i;}    
                }
            }
        }
        
        B_return(ii) = b_i;
    }
    
    return B_return;
}

// [[Rcpp::export]]
arma::field<arma::vec> mh_up(const arma::vec EIDs, const arma::vec &par,
                                const arma::field<arma::uvec> &par_index,
                                arma::field <arma::vec> &B,
                                const arma::vec &y, const arma::vec &eids,
                                int n_cores, int states_per_step) {
    
    // par_index KEY: (0) mu, (1) zeta;
    // "i" is the numeric EID number; "ii" is the index of the EID
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);

    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};
    
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
    
        arma::vec b_i = B(ii);
        arma::vec y_i = y.rows(sub_ind);
    
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - states_per_step + 1; k++) {
        
            // All possible state transitions given current time point -----
            arma::mat Omega_set = get_omega_list(k, n_i, b_i, states_per_step);
            
            if(Omega_set.n_rows > 1) {
                
                arma::vec ss_prob(Omega_set.n_rows, arma::fill::ones);
                ss_prob = (1/arma::accu(ss_prob)) * ss_prob;
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
                arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                    
                double log_like_b = 0;
                double log_like_s = 0;
                arma::vec s_i = b_i;
                s_i.rows(k, k+states_per_step-1) = Omega_set.row(row_ind(0)).t();
                
                if(arma::approx_equal(s_i, b_i, "absdiff", 0.001)) {
                    b_i = s_i;
                } else {
                    arma::vec ones(b_i.n_elem, arma::fill::ones);
                    arma::vec twos_b(b_i.n_elem, arma::fill::zeros);
                    arma::vec twos_s(b_i.n_elem, arma::fill::zeros);
                    
                    twos_b.elem(arma::find(b_i == 2)) += 1;
                    twos_s.elem(arma::find(s_i == 2)) += 1;
                    
                    arma::mat D_b = arma::join_rows(ones, arma::cumsum(twos_b));
                    arma::mat D_s = arma::join_rows(ones, arma::cumsum(twos_s));
                    
                    arma::vec mean_b = D_b * mu;
                    arma::vec mean_s = D_s * mu;
                        
                    for(int kk = k; kk < n_i; kk++) {
                        if(kk == 0) {
                            log_like_b = log_like_b + 
                                log(P_init(b_i(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_b(kk), 1, true);
                            
                            log_like_s = log_like_s + 
                                log(P_init(s_i(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_s(kk), 1, true);
                            
                        } else if(kk <= n_i - states_per_step) {
                            log_like_b = log_like_b + 
                                log(P_i(b_i(kk-1) - 1, b_i(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_b(kk), 1, true);
                            
                            log_like_s = log_like_s + 
                                log(P_i(s_i(kk-1) - 1, s_i(kk) - 1)) + 
                                R::dnorm(y_i(kk), mean_s(kk), 1, true);
                        }
                    }
                    double diff_check = log_like_s - log_like_b;
                    double min_log = log(arma::randu(arma::distr_param(0,1)));
                    if(diff_check > min_log){b_i = s_i;} 
                }
            }
        }
        
        B_return(ii) = b_i;
    }

    return B_return;
    
}

// [[Rcpp::export]]
arma::field<arma::vec> mh_up_all(const arma::vec EIDs, const arma::vec &par,
                                 const arma::field<arma::uvec> &par_index,
                                 arma::field <arma::vec> &B,
                                 const arma::vec &y, const arma::vec &eids,
                                 int n_cores) {
    
    // par_index KEY: (0) mu, (1) zeta;
    // "i" is the numeric EID number; "ii" is the index of the EID
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    
    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};
    
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for(int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::vec b_i = B(ii);
        arma::vec s_i(n_i, arma::fill::zeros);
        arma::vec y_i = y.rows(sub_ind);
        
        arma::vec all_like_vals_b(n_i, arma::fill::zeros);
        arma::vec all_like_vals_s(n_i, arma::fill::zeros);
    
        // Propose a full state sequence ---------------------------------------
        for(int k = 0; k < n_i; k++) {
            
            arma::vec like_vals_s(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
            arma::vec like_vals_b(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
            
            for(int m = 0; m < adj_mat_GLOBAL.n_cols; m++) {
                
                // Necessary computations for mean of candidate state seq (s_i)
                arma::vec s_temp = s_i.subvec(0, k);
                s_temp(k) = m+1;
                arma::vec twos_s(s_temp.n_elem, arma::fill::zeros);
                twos_s.elem(arma::find(s_temp == 2)) += 1;
                arma::vec D_s = {1, arma::accu(twos_s)};
                double mean_s = arma::as_scalar(D_s.t() * mu);
                
                // Necessary computations for mean of current state seq (b_i)
                arma::vec b_temp = b_i.subvec(0, k);
                arma::vec twos_b(b_temp.n_elem, arma::fill::zeros);
                twos_b.elem(arma::find(b_temp == 2)) += 1;
                arma::vec D_b = {1, arma::accu(twos_b)};
                double mean_b = arma::as_scalar(D_b.t() * mu);
                
                if(k == 0) {
                    like_vals_s(m) = P_init(m) * R::dnorm(y_i(k), mean_s, 1, false);
                    like_vals_b(m) = like_vals_s(m);
                } else {
                    like_vals_s(m) = P_i(s_i(k-1) - 1, m) * R::dnorm(y_i(k), mean_s, 1, false);
                    like_vals_b(m) = P_i(b_i(k-1) - 1, m) * R::dnorm(y_i(k), mean_b, 1, false);
                }
            }
            
            all_like_vals_s(k) = arma::accu(like_vals_s);
            all_like_vals_b(k) = arma::accu(like_vals_b);
            
            // Determine sampling distribution for s_i -------------------------
            arma::vec ss_ind = arma::linspace(0, adj_mat_GLOBAL.n_cols-1, 
                                              adj_mat_GLOBAL.n_cols);
            
            double prob_max = like_vals_s.max();
            like_vals_s = like_vals_s / prob_max;
            arma::vec ss_prob = (1/arma::accu(like_vals_s)) * like_vals_s;
            
            arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
            s_i(k) = row_ind(0) + 1;
        }
        
        double diff_check = arma::accu(log(all_like_vals_s)) - arma::accu(log(all_like_vals_b));
        double min_log = log(arma::randu(arma::distr_param(0,1)));
        if(diff_check > min_log){b_i = s_i;} 
        
        B_return(ii) = b_i;
    }
    
    return B_return;
}

// [[Rcpp::export]]
arma::field<arma::vec> mle_state_seq(const arma::vec &EIDs, const arma::vec &par, 
                                     const arma::field<arma::uvec> &par_index,
                                     const arma::vec &y, const arma::vec &eids, 
                                     int n_cores) {
    
    arma::field<arma::vec> B_mle(EIDs.n_elem);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    
    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};
    
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;
    
    // Find the MLE state sequence given the current parameters ----------------
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        arma::vec y_i = y.rows(sub_ind);
        
        arma::vec b_i(n_i, arma::fill::zeros);
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i; k++) {
            if(k == 0) {
                // Initial state
                arma::vec init_vals(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
                for(int jj = 0; jj < init_vals.n_elem; jj++) {
                    b_i(k) = jj + 1;
                    
                    arma::vec b_sub = b_i.subvec(0,k);
                    
                    arma::vec ones(b_sub.n_elem, arma::fill::ones);
                    arma::vec twos(b_sub.n_elem, arma::fill::zeros);
                    twos.elem(arma::find(b_sub == 2)) += 1;
                    
                    arma::vec D = {1, arma::accu(twos)};
                    
                    double mean_b = arma::as_scalar(D.t() * mu);
                    
                    init_vals(jj) = log(P_init(jj)) + R::dnorm(y_i(k), mean_b, 1, true);
                } 
                
                b_i(k) = arma::index_max(init_vals) + 1;
            } else {
                // All others
                int prev_state = b_i(k-1);
                arma::vec poss_next_state(arma::accu(adj_mat_GLOBAL.row(prev_state-1)), arma::fill::zeros);
                arma::vec poss_state_like(arma::accu(adj_mat_GLOBAL.row(prev_state-1)), arma::fill::zeros);
                
                int w_ind = 0;
                for(int jj = 0; jj < adj_mat_GLOBAL.n_cols; jj++) {
                    if(adj_mat_GLOBAL(prev_state-1, jj) != 0) {
                        poss_next_state(w_ind) = jj + 1;
                        
                        b_i(k) = jj + 1;
                        
                        arma::vec b_sub = b_i.subvec(0,k);
                        
                        arma::vec twos(b_sub.n_elem, arma::fill::zeros);
                        twos.elem(arma::find(b_sub == 2)) += 1;
                        
                        arma::vec D = {1, arma::accu(twos)};
                        
                        double mean_b = arma::as_scalar(D.t() * mu);
                        
                        poss_state_like(w_ind) = log(P_i(prev_state-1, jj)) + R::dnorm(y_i(k), mean_b, 1, true);
                        
                        w_ind += 1;
                    }
                }
            
                b_i(k) = poss_next_state(poss_state_like.index_max());
            }
        }
        
        B_mle(ii) = b_i;
    }
    
    return B_mle;
}

// [[Rcpp::export]]
arma::field<arma::vec> viterbi_alg(const arma::vec &EIDs, const arma::vec &par, 
                                   const arma::field<arma::uvec> &par_index,
                                   const arma::vec &y, const arma::vec &eids, 
                                   int n_cores) {
    
    arma::field<arma::vec> B_mle(EIDs.n_elem);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    
    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};
    
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;
    
    // Find the MLE state sequence given the current parameters ----------------
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        arma::vec y_i = y.rows(sub_ind);
        
        arma::vec b_i(n_i, arma::fill::zeros);
        
        arma::mat viterbi_P(n_i, adj_mat_GLOBAL.n_cols, arma::fill::zeros);
        arma::mat viterbi_Q(n_i, adj_mat_GLOBAL.n_cols, arma::fill::zeros);
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i; k++) {
            if(k == 0) {
                // Initial state
                arma::vec p_temp(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
                arma::vec q_temp(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
                
                for(int jj = 0; jj < adj_mat_GLOBAL.n_cols; jj++) {
                    p_temp(jj) = log(P_init(jj)) + R::dnorm(y_i(k), mu(jj), 1, true);
                } 
                
                viterbi_P.row(k) = p_temp.t();
            } else {
                // All others
                arma::vec P_t_1_dot = viterbi_P.row(k-1).t();
                
                for(int s = 0; s < adj_mat_GLOBAL.n_cols; s++) {
                    arma::vec a_dot_s = P_i.col(s);
                    a_dot_s = log(a_dot_s);
                    
                    double b_s_o = R::dnorm(y_i(k), mu(s), 1, true);
                    arma::vec b_s_dot(adj_mat_GLOBAL.n_cols); b_s_dot.fill(b_s_o);
                    
                    arma::vec prod_max_over = P_t_1_dot + a_dot_s + b_s_dot;
                    viterbi_P(k, s) = prod_max_over.max();
                    viterbi_Q(k, s) = prod_max_over.index_max();
                }
            }
        }
        
        for(int k = n_i-1; k > 0; k--) {
            if(k == n_i - 1) {
                b_i(k) = viterbi_P.row(n_i-1).index_max();
            }
            b_i(k-1) = viterbi_Q(k, b_i(k));
        } 
        
        B_mle(ii) = b_i + 1;
    }
    
    return B_mle;
    
}

// [[Rcpp::export]]
void test_fnc() {
    // Rcpp::Rcout << Omega_List_GLOBAL_multi << std::endl;
    // Rcpp::Rcout << adj_mat_GLOBAL << std::endl;
    // 
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(1) << std::endl;
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(1)(0, 1) << std::endl;
    
    arma::vec test = {0,1,2,3,0};
    arma::vec test2 = log(test);
    Rcpp::Rcout << test2.t() << std::endl;
    Rcpp::Rcout << test2.max() << std::endl;
    Rcpp::Rcout << test2.index_max() << std::endl;
    Rcpp::Rcout << test2.min() << std::endl;
    Rcpp::Rcout << test2.index_min() << std::endl;
    
}



