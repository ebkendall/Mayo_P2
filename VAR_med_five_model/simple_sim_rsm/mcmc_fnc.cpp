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
arma::field<arma::field<arma::mat>> Omega_List_GLOBAL_multi;

// [[Rcpp::export]]
void initialize_cpp(arma::imat a_mat) {
    adj_mat_GLOBAL = a_mat;
    Omega_List_GLOBAL_multi = get_Omega_list(adj_mat_GLOBAL, 2);
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

double log_f_i_cpp_total(const arma::vec &EIDs, const arma::vec &par,
                         const arma::field<arma::uvec> &par_index,
                         const arma::field <arma::vec> &B, const arma::mat &y,
                         const arma::vec &eids, int n_cores, bool dgm, 
                         const arma::mat &gamma_i) {
    
    // (0) alpha, (1) zeta, (2) R, (3) init, (4) G -----------------------------

    // Parameter initialization ------------------------------------------------
    arma::mat alpha_miss_y;
    if(dgm) {
        alpha_miss_y = arma::reshape(par.elem(par_index(0) - 1), 3, 4);
    } else {
        alpha_miss_y = arma::reshape(par.elem(par_index(0) - 1), 2, 4);   
    }
    
    arma::vec zeta = par.elem(par_index(1) - 1);
    
    arma::mat R = arma::diagmat(exp(par.elem(par_index(2) - 1)));
    
    arma::mat G(R.n_rows, R.n_cols, arma::fill::zeros);
    if(!dgm) {
        G = arma::diagmat(exp(par.elem(par_index(4) - 1)));
    }

    arma::vec lp_temp = exp(par.elem(par_index(3) - 1));
    arma::vec logit_prob = {1, lp_temp(0), lp_temp(1)};
    arma::vec init_prob = logit_prob / arma::accu(logit_prob);

    arma::vec qz = exp(zeta);
    arma::mat Q = { {1, qz(0), 0}, {0, 1, qz(1)}, {qz(2), qz(3), 1}};

    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P = Q.each_col() / q_row_sums;
    // -------------------------------------------------------------------------
    
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        double like_comp = 0;
        double like_comp_transition = 0;

        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        arma::vec b_i = B(ii);
        arma::mat y_i = y.rows(sub_ind);

        arma::mat alpha_i(3, 4);
        if(dgm) {
            alpha_i = alpha_miss_y;
        } else {
            alpha_i.row(0) = gamma_i.row(ii); 
            alpha_i.row(1) = alpha_miss_y.row(0);
            alpha_i.row(2) = alpha_miss_y.row(1);    
        }

        for(int jj = 0; jj < n_i; jj++) {
            
            if(jj == 0) {
                
                arma::vec log_y_pdf = dmvnorm(y_i.row(jj), alpha_i.row(0).t(), R, true);
                
                double log_g_val = 0.0;
                if(!dgm) {
                    arma::vec log_g_pdf = dmvnorm(alpha_i.row(0), y_i.row(jj).t(), G, true);
                    log_g_val = arma::as_scalar(log_g_pdf);
                }
                
                like_comp = like_comp + log_g_val + arma::as_scalar(log_y_pdf);
                like_comp_transition = like_comp_transition + log(init_prob(b_i(jj) - 1));
                
            } else {
                
                int curr_b = b_i(jj);
                int curr_b_1 = b_i(jj-1);
                
                arma::vec b_sub = b_i.subvec(1, jj); // start = 1
                arma::vec twos(b_sub.n_elem, arma::fill::zeros);
                arma::vec threes(b_sub.n_elem, arma::fill::zeros);
                twos.elem(arma::find(b_sub == 2)) += 1;
                threes.elem(arma::find(b_sub == 3)) += 1;
                
                arma::vec mean_b = alpha_i.row(0).t() + 
                    arma::accu(twos) * alpha_i.row(1).t() + 
                    arma::accu(threes) * alpha_i.row(2).t();
                
                arma::vec log_y_pdf = dmvnorm(y_i.row(jj), mean_b, R, true);
                
                like_comp = like_comp + arma::as_scalar(log_y_pdf);
                like_comp_transition = like_comp_transition + log(P(curr_b_1 - 1, curr_b - 1));   
            }
        }

        in_vals(ii) = like_comp + like_comp_transition;
    }

    double in_value = arma::accu(in_vals);

    return in_value;
}

// [[Rcpp::export]]
double log_post_cpp(const arma::vec &EIDs, const arma::vec &par,
                    const arma::field<arma::uvec> &par_index,
                    const arma::field<arma::vec> &B, const arma::mat &y,
                    const arma::vec &eids, int n_cores, bool dgm, 
                    const arma::mat &gamma_i) {
    
    // (0) alpha, (1) zeta, (2) R, (3) init, (4) G -----------------------------

    // Compute the likelihood --------------------------------------------------
    double value = log_f_i_cpp_total(EIDs, par, par_index, B, y, eids, n_cores, dgm, gamma_i);

    // Prior densities
    arma::vec prior_mean;
    if(dgm) {
        prior_mean = {50,  -5,   5, 100,  10, -10, 100, -10,  10, 50,   5,  -5,
                      -2, -2, -1.5, -1.5,
                      0, 0, 0, 0,
                      0, 0};
    } else {
        prior_mean = {-5, 5, 10, -10, -10, 10, 5, -5,
                      -2, -2, -1.5, -1.5,
                      0, 0, 0, 0,
                      0, 0,
                      0, 0, 0, 0};
    }

    arma::vec prior_var_diag(prior_mean.n_elem, arma::fill::ones);
    prior_var_diag = 100 * prior_var_diag;
    arma::mat prior_var = arma::diagmat(prior_var_diag);

    arma::vec prior_dens = dmvnorm(par.t(), prior_mean, prior_var, true);
    double prior_dens_val = arma::as_scalar(prior_dens);
    
    value = value + prior_dens_val;
    
    return value;
}

arma::vec p_2_sampler(int n_i, arma::mat &y_i, arma::imat adj_mat_i,
                      arma::vec &b_i, arma::mat alpha_i, arma::mat R,
                      arma::mat P, arma::vec init_prob) {
    
    // Because we ignore the first time point, the length of the state sequence
    // we care about is 1 less than the original.
    for(int k = 0; k < n_i - 1; k++) {

        arma::vec s_i = b_i;
        
        // arma::mat omega_set = get_omega_list(k-1, n_i - 1, b_i.subvec(1, n_i-1), 2);
        arma::mat omega_set = get_omega_list(k, n_i, b_i, 2);
        
        arma::vec prob_omega(omega_set.n_rows, arma::fill::ones);
        prob_omega = (1/arma::accu(prob_omega)) * prob_omega;
        arma::vec ind_omega = arma::linspace(0, omega_set.n_rows-1, omega_set.n_rows);
        arma::vec row_omega = RcppArmadillo::sample(ind_omega, 1, false, prob_omega);

        s_i.subvec(k, k + 1) = omega_set.row(row_omega(0)).t();

        // Step 3: compute MH-ratio to accept/reject -------------------
        if(arma::accu(arma::abs(s_i - b_i)) != 0) {

            double log_like_s = 0;
            double log_like_b = 0;

            arma::vec twos_s(s_i.n_elem, arma::fill::zeros);
            arma::vec threes_s(s_i.n_elem, arma::fill::zeros);
            twos_s.elem(arma::find(s_i == 2)) += 1;
            threes_s.elem(arma::find(s_i == 3)) += 1;

            arma::vec twos_b(b_i.n_elem, arma::fill::zeros);
            arma::vec threes_b(b_i.n_elem, arma::fill::zeros);
            twos_b.elem(arma::find(b_i == 2)) += 1;
            threes_b.elem(arma::find(b_i == 3)) += 1;

            for(int t = k; t < n_i; t++) {
                
                if(t == 0) {
                    log_like_s = log_like_s + log(init_prob(s_i(t) - 1));
                    log_like_b = log_like_b + log(init_prob(b_i(t) - 1));
                } else {
                    // INDEXING STARTS AT 1 (NOT 0)
                    // Computations for mean of candidate (s_i) ------------
                    arma::vec s_2 = twos_s.subvec(1, t);
                    arma::vec s_3 = threes_s.subvec(1, t);
                    
                    arma::vec mean_s = alpha_i.row(0).t() + 
                        arma::accu(s_2) * alpha_i.row(1).t() +
                        arma::accu(s_3) * alpha_i.row(2).t();
                    
                    // Computations for mean of current (b_i) --------------
                    arma::vec b_2 = twos_b.subvec(1, t);
                    arma::vec b_3 = threes_b.subvec(1, t);
                    
                    arma::vec mean_b = alpha_i.row(0).t() +
                        arma::accu(b_2) * alpha_i.row(1).t() +
                        arma::accu(b_3) * alpha_i.row(2).t();
                    
                    arma::vec log_y_pdf_s = dmvnorm(y_i.row(t), mean_s, R, true);
                    arma::vec log_y_pdf_b = dmvnorm(y_i.row(t), mean_b, R, true); 
                    
                    if(t <= k + 2) {
                        log_like_s = log_like_s + log(P(s_i(t-1)-1, s_i(t)-1));
                        log_like_b = log_like_b + log(P(b_i(t-1)-1, b_i(t)-1));
                    }
                    
                    log_like_s = log_like_s + arma::as_scalar(log_y_pdf_s);
                    log_like_b = log_like_b + arma::as_scalar(log_y_pdf_b);    
                }
            }

            double diff_check = log_like_s - log_like_b;

            double min_log = log(arma::randu(arma::distr_param(0,1)));
            if(diff_check > min_log){b_i = s_i;}
        }
    }
    
    return b_i;
}

arma::vec p_full_sampler(int n_i, arma::mat &y_i, arma::imat adj_mat_i,
                         arma::vec &b_i, arma::mat alpha_i, arma::mat R,
                         arma::mat P, arma::vec init_prob) {
    
    // Because we ignore the first time point, the length of the state sequence
    // we care about is 1 less than the original.
    
    arma::vec all_like_vals_b(n_i, arma::fill::zeros);
    arma::vec all_like_vals_s(n_i, arma::fill::zeros);
    
    arma::vec s_i(n_i, arma::fill::zeros);

    // Propose a full state sequence -------------------------------------------
    for (int k = 0; k < n_i; k++) {

        arma::vec like_vals_s;
        arma::vec like_vals_b;
        arma::vec ss_ind;

        if(k == 0) {
            like_vals_s = log(init_prob);
            like_vals_b = log(init_prob);

            ss_ind = arma::linspace(0, adj_mat_i.n_cols-1, adj_mat_i.n_cols);

        } else {
            
            int prev_state_s = s_i(k-1);
            int prev_state_b = b_i(k-1);

            like_vals_s.zeros(arma::accu(adj_mat_i.row(prev_state_s-1)));
            like_vals_b.zeros(arma::accu(adj_mat_i.row(prev_state_b-1)));

            arma::vec state_vals_s(arma::accu(adj_mat_i.row(prev_state_s-1)), arma::fill::zeros);
            arma::vec state_vals_b(arma::accu(adj_mat_i.row(prev_state_b-1)), arma::fill::zeros);

            int s_ind = 0;
            int b_ind = 0;
            
            // INDEXING STARTS AT 1 (NOT 0)
            for(int m = 0; m < adj_mat_i.n_cols; m++) {

                // Computations for mean of candidate (s_i) --------------------
                if(adj_mat_i(prev_state_s-1, m) != 0) {
                    
                    arma::vec s_temp = s_i.subvec(1, k);
                    s_temp(k-1) = m+1;
                    arma::vec twos_s(s_temp.n_elem, arma::fill::zeros);
                    arma::vec threes_s(s_temp.n_elem, arma::fill::zeros);
                    twos_s.elem(arma::find(s_temp == 2)) += 1;
                    threes_s.elem(arma::find(s_temp == 3)) += 1;
                    
                    arma::vec mean_s = alpha_i.row(0).t() + 
                        arma::accu(twos_s) * alpha_i.row(1).t() + 
                        arma::accu(threes_s) * alpha_i.row(2).t();
                    
                    arma::vec log_y_pdf = dmvnorm(y_i.row(k), mean_s, R, true);
                    
                    like_vals_s(s_ind) = log(P(prev_state_s-1, m)) + arma::as_scalar(log_y_pdf);
                    state_vals_s(s_ind) = m;
                    s_ind = s_ind + 1;
                }

                // Computations for mean of current (b_i) ----------------------
                if(adj_mat_i(prev_state_b-1, m) != 0) {
                    
                    arma::vec b_temp = b_i.subvec(1, k);
                    b_temp(k-1) = m+1;
                    arma::vec twos_b(b_temp.n_elem, arma::fill::zeros);
                    arma::vec threes_b(b_temp.n_elem, arma::fill::zeros);
                    twos_b.elem(arma::find(b_temp == 2)) += 1;
                    threes_b.elem(arma::find(b_temp == 3)) += 1;
                    
                    arma::vec mean_b = alpha_i.row(0).t() + 
                        arma::accu(twos_b) * alpha_i.row(1).t() + 
                        arma::accu(threes_b) * alpha_i.row(2).t();
                    
                    arma::vec log_y_pdf = dmvnorm(y_i.row(k), mean_b, R, true);
                    
                    like_vals_b(b_ind) = log(P(prev_state_b-1, m)) + arma::as_scalar(log_y_pdf);
                    state_vals_b(b_ind) = m;
                    b_ind = b_ind + 1;
                }
            }

            ss_ind = state_vals_s;
        }

        all_like_vals_s(k) = arma::accu(exp(like_vals_s));
        all_like_vals_b(k) = arma::accu(exp(like_vals_b));

        // Determine sampling distribution for s_i -------------------------
        double prob_log_max = like_vals_s.max();
        like_vals_s = like_vals_s - prob_log_max;
        like_vals_s = exp(like_vals_s);
        arma::vec ss_prob = (1/arma::accu(like_vals_s)) * like_vals_s;

        arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
        s_i(k) = row_ind(0) + 1;
    }

    double diff_check = arma::accu(log(all_like_vals_s)) - arma::accu(log(all_like_vals_b));
    double min_log = log(arma::randu(arma::distr_param(0,1)));
    if(diff_check > min_log){ b_i = s_i; }
    
    return b_i;
}

arma::vec p_flex_sampler(int n_i, arma::mat &y_i, arma::imat adj_mat_i,
                         arma::vec &b_i, arma::mat alpha_i, arma::mat R, 
                         arma::mat P, arma::vec init_prob, int sps) {
    
    int k = 0;
    
    while(k < n_i - 2) {
        
        arma::vec s_i = b_i;
        
        arma::vec all_like_vals_b(sps - 2, arma::fill::zeros);
        arma::vec all_like_vals_s(sps - 2, arma::fill::zeros);
        
        // Step 1: sample 1-at-a-time ----------------------------------
        arma::ivec t_lim = {k + sps - 2, n_i - 2};
        int t_max = arma::min(t_lim);
        for(int t = k; t < t_max; t++) {
            
            arma::vec like_vals_s;
            arma::vec like_vals_b;
            arma::vec ss_ind;
            
            if(t == 0) {
                
                like_vals_s = log(init_prob);
                like_vals_b = log(init_prob);
                
                ss_ind = arma::linspace(0, adj_mat_i.n_cols-1, adj_mat_i.n_cols);
                
            } else {
                int prev_state_s = s_i(t-1);
                int prev_state_b = b_i(t-1);
                
                like_vals_s.zeros(arma::accu(adj_mat_i.row(prev_state_s-1)));
                like_vals_b.zeros(arma::accu(adj_mat_i.row(prev_state_b-1)));
                
                arma::vec state_vals_s(arma::accu(adj_mat_i.row(prev_state_s-1)), arma::fill::zeros);
                arma::vec state_vals_b(arma::accu(adj_mat_i.row(prev_state_b-1)), arma::fill::zeros);
                
                int s_ind = 0;
                int b_ind = 0;
                
                // INDEXING STARTS AT 1 (NOT 0)
                for(int m = 0; m < adj_mat_i.n_cols; m++) {
                    
                    // Computations for mean of candidate (s_i) --------
                    if(adj_mat_i(prev_state_s-1, m) != 0) {
                        
                        arma::vec s_temp = s_i.subvec(1, t); 
                        s_temp(t-1) = m+1;
                        arma::vec twos_s(s_temp.n_elem, arma::fill::zeros);
                        arma::vec threes_s(s_temp.n_elem, arma::fill::zeros);
                        twos_s.elem(arma::find(s_temp == 2)) += 1;
                        threes_s.elem(arma::find(s_temp == 3)) += 1;
                        
                        arma::vec mean_s = alpha_i.row(0).t() + 
                            arma::accu(twos_s) * alpha_i.row(1).t() + 
                            arma::accu(threes_s) * alpha_i.row(2).t();
                        
                        arma::vec log_y_pdf = dmvnorm(y_i.row(t), mean_s, R, true);
                        
                        like_vals_s(s_ind) = log(P(prev_state_s-1, m)) + arma::as_scalar(log_y_pdf);
                        state_vals_s(s_ind) = m;
                        s_ind = s_ind + 1;
                    }
                    
                    // Computations for mean of current (b_i) ----------
                    if(adj_mat_i(prev_state_b-1, m) != 0) {
                        
                        arma::vec b_temp = b_i.subvec(1, t);
                        b_temp(t-1) = m+1;
                        arma::vec twos_b(b_temp.n_elem, arma::fill::zeros);
                        arma::vec threes_b(b_temp.n_elem, arma::fill::zeros);
                        twos_b.elem(arma::find(b_temp == 2)) += 1;
                        threes_b.elem(arma::find(b_temp == 3)) += 1;
                        
                        arma::vec mean_b = alpha_i.row(0).t() + 
                            arma::accu(twos_b) * alpha_i.row(1).t() + 
                            arma::accu(threes_b) * alpha_i.row(2).t();
                        
                        arma::vec log_y_pdf = dmvnorm(y_i.row(t), mean_b, R, true);
                        
                        like_vals_b(b_ind) = log(P(prev_state_b-1, m)) + arma::as_scalar(log_y_pdf);
                        state_vals_b(b_ind) = m;
                        b_ind = b_ind + 1;
                    }
                }
                
                ss_ind = state_vals_s;
            }
            
            all_like_vals_s(t-k) = arma::accu(exp(like_vals_s));
            all_like_vals_b(t-k) = arma::accu(exp(like_vals_b));
            
            // Determine sampling distribution for s_i -----------------
            double prob_log_max = like_vals_s.max();
            like_vals_s = like_vals_s - prob_log_max;
            like_vals_s = exp(like_vals_s);
            arma::vec ss_prob = (1/arma::accu(like_vals_s)) * like_vals_s;
            
            arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
            s_i(t) = row_ind(0) + 1;
        }
        
        // Step 2: sample the last 2 times together --------------------
        arma::mat Omega_set_s = get_omega_list(t_max, n_i, s_i, 2);
        arma::mat Omega_set_b = get_omega_list(t_max, n_i, b_i, 2);
        
        arma::vec prob_omega_s(Omega_set_s.n_rows, arma::fill::ones);
        prob_omega_s = (1/arma::accu(prob_omega_s)) * prob_omega_s;
        arma::vec ind_omega_s = arma::linspace(0, Omega_set_s.n_rows-1, Omega_set_s.n_rows);
        arma::vec row_omega_s = RcppArmadillo::sample(ind_omega_s, 1, false, prob_omega_s);
        
        s_i.subvec(t_max, t_max + 1) = Omega_set_s.row(row_omega_s(0)).t();
        
        // Step 3: compute MH-ratio to accept/reject -------------------
        if(arma::accu(arma::abs(s_i - b_i)) != 0) {
            
            double log_prob_diff = 0;
            double log_like_diff = 0;
            
            arma::vec twos_s(s_i.n_elem, arma::fill::zeros);
            arma::vec threes_s(s_i.n_elem, arma::fill::zeros);
            twos_s.elem(arma::find(s_i == 2)) += 1;
            threes_s.elem(arma::find(s_i == 3)) += 1;
            
            arma::vec twos_b(b_i.n_elem, arma::fill::zeros);
            arma::vec threes_b(b_i.n_elem, arma::fill::zeros);
            twos_b.elem(arma::find(b_i == 2)) += 1;
            threes_b.elem(arma::find(b_i == 3)) += 1;
            
            // INDEXING STARTS AT 1 (NOT 0)
            for(int t = t_max; t < n_i; t++) {
                
                // Computations for mean of candidate (s_i) ------------
                arma::vec s_2 = twos_s.subvec(1, t);
                arma::vec s_3 = threes_s.subvec(1, t);
                
                arma::vec mean_s = alpha_i.row(0).t() +
                    arma::accu(s_2) * alpha_i.row(1).t() +
                    arma::accu(s_3) * alpha_i.row(2).t();
                
                // Computations for mean of current (b_i) --------------
                arma::vec b_2 = twos_b.subvec(1, t);
                arma::vec b_3 = threes_b.subvec(1, t);
                
                arma::vec mean_b = alpha_i.row(0).t() +
                    arma::accu(b_2) * alpha_i.row(1).t() +
                    arma::accu(b_3) * alpha_i.row(2).t();
                
                arma::vec log_y_pdf_s = dmvnorm(y_i.row(t), mean_s, R, true);
                arma::vec log_y_pdf_b = dmvnorm(y_i.row(t), mean_b, R, true);
                
                if(t < k + sps + 1) {
                    log_prob_diff = log_prob_diff +
                        log(P(s_i(t-1)-1, s_i(t)-1)) - log(P(b_i(t-1)-1, b_i(t)-1));
                }
                
                log_like_diff = log_like_diff +
                    arma::as_scalar(log_y_pdf_s) - arma::as_scalar(log_y_pdf_b);
            }
            
            double log_diff_1 = arma::accu(log(all_like_vals_s)) - arma::accu(log(all_like_vals_b));
            double log_diff_2 = log(Omega_set_s.n_rows) - log(Omega_set_b.n_rows);
            
            double diff_check = log_prob_diff + log_like_diff + log_diff_1 + log_diff_2;
            
            double min_log = log(arma::randu(arma::distr_param(0,1)));
            if(diff_check > min_log){b_i = s_i;}
        }
        
        k = k + sps - 2;
    }
    
    return b_i;
}

// [[Rcpp::export]]
arma::field<arma::vec> state_sampler(const arma::vec EIDs, const arma::vec &par,
                                     const arma::field<arma::uvec> &par_index,
                                     arma::field <arma::vec> &B,
                                     const arma::mat &y, const arma::vec &eids,
                                     int n_cores, int states_per_step, bool dgm,
                                     const arma::mat &gamma_i) {
    
    // (0) alpha, (1) zeta, (2) R, (3) init, (4) G -----------------------------
    
    // "i" is the numeric EID number; "ii" is the index of the EID
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    // Parameter initialization ------------------------------------------------
    arma::mat alpha_miss_y;
    if(dgm) {
        alpha_miss_y = arma::reshape(par.elem(par_index(0) - 1), 3, 4);
    } else {
        alpha_miss_y = arma::reshape(par.elem(par_index(0) - 1), 2, 4);   
    }
    
    arma::vec zeta = par.elem(par_index(1) - 1);
    
    arma::mat R = arma::diagmat(exp(par.elem(par_index(2) - 1)));
    
    arma::vec lp_temp = exp(par.elem(par_index(3) - 1));
    arma::vec logit_prob = {1, lp_temp(0), lp_temp(1)};
    arma::vec init_prob = logit_prob / arma::accu(logit_prob);

    arma::vec qz = exp(zeta);
    arma::mat Q = { {1, qz(0), 0}, {0, 1, qz(1)}, {qz(2), qz(3), 1}};

    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P = Q.each_col() / q_row_sums;
    // -------------------------------------------------------------------------
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::vec b_i = B(ii);
        arma::mat y_i = y.rows(sub_ind);
        
        arma::imat adj_mat_i = adj_mat_GLOBAL;

        arma::mat alpha_i(3, 4);
        if(dgm) {
            alpha_i = alpha_miss_y;
        } else {
            alpha_i.row(0) = gamma_i.row(ii);
            alpha_i.row(1) = alpha_miss_y.row(0);
            alpha_i.row(2) = alpha_miss_y.row(1);
        }
        
        arma::vec new_b_i;
        if(states_per_step == 2) {
            // p = 2
            new_b_i = p_2_sampler(n_i, y_i, adj_mat_i, b_i, alpha_i,
                                  R, P, init_prob);
            
        } else if(states_per_step >= n_i - 1) {
            // p >= n_i - 1 (first time point doesn't count)
            new_b_i = p_full_sampler(n_i, y_i, adj_mat_i, b_i, alpha_i,
                                     R, P, init_prob);
            
        } else {
            // 2 < p < n_i - 1
            new_b_i = p_flex_sampler(n_i, y_i, adj_mat_i, b_i, alpha_i,
                                     R, P, init_prob, states_per_step);
        }
        
        B_return(ii) = new_b_i;
    }
    
    return B_return;
}

// [[Rcpp::export]]
arma::field<arma::vec> mle_state_seq(const arma::vec &EIDs, const arma::vec &par,
                                     const arma::field<arma::uvec> &par_index,
                                     const arma::mat &y, const arma::vec &eids,
                                     const bool dgm) {

    // (0) alpha, (1) zeta, (2) R, (3) init, (4) G -----------------------------
    arma::field<arma::vec> B_mle(EIDs.n_elem);

    // Parameter initialization ------------------------------------------------
    arma::mat alpha_miss_y;
    if(dgm) {
        alpha_miss_y = arma::reshape(par.elem(par_index(0) - 1), 3, 4);
    } else {
        alpha_miss_y = arma::reshape(par.elem(par_index(0) - 1), 2, 4);   
    }
    
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::mat R = arma::diagmat(exp(par.elem(par_index(2) - 1)));
    
    arma::mat G(R.n_rows, R.n_cols, arma::fill::zeros);
    if(!dgm) {
        G = arma::diagmat(exp(par.elem(par_index(4) - 1)));
    }

    arma::vec lp_temp = exp(par.elem(par_index(3) - 1));
    arma::vec logit_prob = {1, lp_temp(0), lp_temp(1)};
    arma::vec init_prob = logit_prob / arma::accu(logit_prob);

    arma::vec qz = exp(zeta);
    arma::mat Q = { {1, qz(0), 0}, {0, 1, qz(1)}, {qz(2), qz(3), 1}};

    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P = Q.each_col() / q_row_sums;

    // Find the MLE state sequence given the current parameters ----------------
    for(int ii = 0; ii < EIDs.n_elem; ii++) {

        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        arma::mat y_i = y.rows(sub_ind);

        arma::vec b_i(n_i, arma::fill::zeros);
        
        arma::mat alpha_i(3, 4);
        if(dgm) {
            alpha_i = alpha_miss_y;
        } else {
            alpha_i.row(0) = y_i.row(0); 
            alpha_i.row(1) = alpha_miss_y.row(0);
            alpha_i.row(2) = alpha_miss_y.row(1);    
        }

        // Looping through subject state space ---------------------------------
        for(int k = 0; k < n_i; k++) {
            
            if(k == 0) {
                
                b_i(k) = arma::index_max(init_prob) + 1;
                
            } else {
                
                int prev_state = b_i(k-1);
                arma::vec poss_next_state(arma::accu(adj_mat_GLOBAL.row(prev_state-1)), arma::fill::zeros);
                arma::vec poss_state_like(arma::accu(adj_mat_GLOBAL.row(prev_state-1)), arma::fill::zeros);

                int w_ind = 0;
                for(int jj = 0; jj < adj_mat_GLOBAL.n_cols; jj++) {
                    if(adj_mat_GLOBAL(prev_state-1, jj) != 0) {
                        
                        b_i(k) = jj + 1;
                        
                        arma::vec b_sub = b_i.subvec(1, k); // start = 1
                        arma::vec twos(b_sub.n_elem, arma::fill::zeros);
                        arma::vec threes(b_sub.n_elem, arma::fill::zeros);
                        twos.elem(arma::find(b_sub == 2)) += 1;
                        threes.elem(arma::find(b_sub == 3)) += 1;
                        
                        arma::vec mean_b = alpha_i.row(0).t() + arma::accu(twos) * alpha_i.row(1).t() + 
                            arma::accu(threes) * alpha_i.row(2).t();

                        arma::vec log_y_pdf = dmvnorm(y_i.row(k), mean_b, R, true);
                        
                        poss_next_state(w_ind) = jj + 1;
                        poss_state_like(w_ind) = log(P(prev_state-1, jj)) + arma::as_scalar(log_y_pdf);

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
arma::mat gamma_i_sample(const arma::vec &EIDs, const arma::vec &par,
                         const arma::field<arma::uvec> &par_index,
                         const arma::field <arma::vec> &B, const arma::mat &y,
                         const arma::vec &eids, int n_cores) {

    // (0) alpha, (1) zeta, (2) R, (3) init, (4) G -----------------------------

    // Parameter initialization ------------------------------------------------
    arma::mat alpha_miss_y = arma::reshape(par.elem(par_index(0) - 1), 2, 4);
    
    arma::vec vec_R = exp(par.elem(par_index(2) - 1));
    arma::mat R = arma::diagmat(vec_R);
    arma::mat inv_R = arma::diagmat(1 / vec_R);

    arma::vec vec_G = exp(par.elem(par_index(4) - 1));
    arma::mat G = arma::diagmat(vec_G); 
    arma::mat inv_G = arma::diagmat(1 / vec_G);
    // -------------------------------------------------------------------------

    arma::mat gamma_i(EIDs.n_elem, 4, arma::fill::zeros);

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        arma::vec b_i = B(ii);
        arma::mat y_i = y.rows(sub_ind);

        arma::vec hold1(y_i.n_cols, arma::fill::zeros);

        for(int jj = 1; jj < n_i; jj++) {

            arma::vec b_sub = b_i.subvec(1, jj); // start = 1
            arma::vec twos(b_sub.n_elem, arma::fill::zeros);
            arma::vec threes(b_sub.n_elem, arma::fill::zeros);
            twos.elem(arma::find(b_sub == 2)) += 1;
            threes.elem(arma::find(b_sub == 3)) += 1;

            arma::vec mean_b = arma::accu(twos) * alpha_miss_y.row(0).t() +
                arma::accu(threes) * alpha_miss_y.row(1).t();

            hold1 = hold1 + (y_i.row(jj).t() - mean_b);
        }

        arma::mat inv_W_i = inv_G + n_i * inv_R;
        arma::mat W_i = arma::inv_sympd(inv_W_i);

        arma::vec V_i = (inv_R + inv_G) * y_i.row(0).t() + inv_R * hold1;
        
        arma::vec gamma_i_mu = W_i * V_i;

        gamma_i.row(ii) = rmvnorm(1, gamma_i_mu, W_i);
    }

    return gamma_i;
}

// [[Rcpp::export]]
void test_fnc(const arma::vec &par,
              const arma::field<arma::uvec> &par_index) {
    
    arma::mat R = arma::diagmat(exp(par.elem(par_index(2) - 1)));
    Rcpp::Rcout << R << std::endl;
    
    arma::vec test = {1,2,3,4,5,6};
    Rcpp::Rcout << test << std::endl;
    test.subvec(1,2) = {0,0};
    Rcpp::Rcout << test << std::endl;
    
    // Rcpp::Rcout << "() -> () -> 1" << std::endl;
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(0)(0) << std::endl;
    // Rcpp::Rcout << "() -> () -> 2" << std::endl;
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(0)(1) << std::endl;
    // Rcpp::Rcout << "() -> () -> 3" << std::endl;
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(0)(2) << std::endl;
    // 
    // for(int i = 0; i < 3; i++){
    //     for(int j = 0; j < 3; j++){
    //         Rcpp::Rcout << i + 1 << " -> () -> () -> " << j + 1 << std::endl;
    //         Rcpp::Rcout << Omega_List_GLOBAL_multi(1)(i, j) << std::endl;
    //     } 
    // }
    // 
    // Rcpp::Rcout << "1 -> () -> ()" << std::endl;
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(2)(0) << std::endl;
    // Rcpp::Rcout << "2 -> () -> ()" << std::endl;
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(2)(1) << std::endl;
    // Rcpp::Rcout << "3 -> () -> ()" << std::endl;
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(2)(2) << std::endl;
    // 
    // Rcpp::Rcout << adj_mat_GLOBAL << std::endl;

    // Rcpp::Rcout << Omega_List_GLOBAL_multi(1) << std::endl;
    // Rcpp::Rcout << Omega_List_GLOBAL_multi(1)(0, 1) << std::endl;
    
    // arma::vec test = {0,1,2,3,0};
    // arma::vec test2 = log(test);
    // Rcpp::Rcout << test2.t() << std::endl;
    // Rcpp::Rcout << test2.max() << std::endl;
    // Rcpp::Rcout << test2.index_max() << std::endl;
    // Rcpp::Rcout << test2.min() << std::endl;
    // Rcpp::Rcout << test2.index_min() << std::endl;
    
    // int n_i = 10;
    // int states_per_step = 4;
    // for (int k = 0; k < n_i - states_per_step + 1; k++) {
    //     
    //     double log_prob_diff = 0;
    //     double log_like_diff = 0;
    //     
    //     for(int kk = k + states_per_step; kk < n_i; kk++) {
    //         Rcpp::Rcout << "(k, kk) = (" << k << ", " << kk << ")" << std::endl;
    //     }
    //     double diff_check = log_prob_diff + log_like_diff;
    //     double min_log = log(arma::randu(arma::distr_param(0,1)));
    //     if(diff_check > min_log){Rcpp::Rcout << "accept" << std::endl;}
    // }
}
