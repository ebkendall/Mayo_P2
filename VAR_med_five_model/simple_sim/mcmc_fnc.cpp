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

// [[Rcpp::export]]
void initialize_cpp(arma::imat a_mat, int s_per_s) {
    adj_mat_GLOBAL = a_mat;
    states_per_step_GLOBAL = s_per_s;
    Omega_List_GLOBAL_multi = get_Omega_list(adj_mat_GLOBAL, states_per_step_GLOBAL);
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

// [[Rcpp::export]]
arma::vec full_like_MH(const int k, const int n_i, int t_pt_length, const arma::vec &par, 
                       const arma::field<arma::uvec> &par_index, const arma::mat &y_i, 
                       arma::vec &pr_b_i, arma::vec &curr_b_i) {
    // Parameter initialization ------------------------------------------------
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::vec P_init = {0.5,0.5}; 
    
    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;
    
    arma::vec t_pts;
    if (k == 1) {
        // () -> () -> 1-5
        t_pts = arma::linspace(0, n_i - 1, n_i);
        
    } else if (k <= n_i - t_pt_length) {
        // 1-5 -> () -> () -> 1-5
        t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
        
    } else if (k == n_i - t_pt_length + 1) { 
        // 1-5 -> () -> ()
        t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
    }
    // ------------------------------------------------------------------------
    
    // Likelihood computations ------------------------------------------------
    double like_comp_prop_transition = 0;
    double like_comp_curr_transition = 0;
    double like_comp_prop = 0;
    double like_comp_curr = 0;
    int jj_start = 0;
    
    for(int jj = jj_start; jj < t_pts.n_elem; jj++) {
        
        int t_j = t_pts(jj);
        
        // (1) Transition probabilities ------------------------------------
        if(jj < t_pt_length + 1) {
            if(t_j == 0) {
                // Current -----------
                like_comp_curr_transition = like_comp_curr_transition + log(P_init(curr_b_i(t_j) - 1));
                
                // Proposal  ---------
                like_comp_prop_transition = like_comp_prop_transition + log(P_init(pr_b_i(t_j) - 1));
                
            } else{ 
                // Current -----------
                int curr_b_k_1 = curr_b_i(t_j-1);
                int curr_b_k   = curr_b_i(t_j);
                
                like_comp_curr_transition = like_comp_curr_transition + 
                    log(P_i(curr_b_k_1 - 1, curr_b_k - 1));
                
                // Proposal  ---------
                int prop_b_k_1 = pr_b_i(t_j-1);
                int prop_b_k   = pr_b_i(t_j);
                
                like_comp_prop_transition = like_comp_prop_transition + 
                    log(P_i(prop_b_k_1 - 1, prop_b_k - 1));
            } 
        }
        
        // (2) Likelihood component ------------------------------------
        // Current ---------
        like_comp_curr = like_comp_curr + R::dnorm(y_i(t_j), mu(curr_b_i(t_j) - 1), 1, true);
        
        // Proposal  ---------
        like_comp_prop = like_comp_prop + R::dnorm(y_i(t_j), mu(pr_b_i(t_j) - 1), 1, true);
    } 
    
    arma::vec mh_like = {like_comp_prop + like_comp_prop_transition, 
                         like_comp_curr + like_comp_curr_transition};
    return mh_like;
    
}

// [[Rcpp::export]]
arma::vec sub_like_MH(const int k, const int n_i, int t_pt_length, const arma::vec &par, 
                      const arma::field<arma::uvec> &par_index, const arma::vec &y_i, 
                      arma::vec &pr_b_i, arma::vec &curr_b_i, int option_num) {
    // Parameter initialization ------------------------------------------------
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::vec P_init = {0.5,0.5}; 
    
    arma::vec t_pts;
    if (k == 1) {
        // () -> () -> 1-5
        t_pts = arma::linspace(0, n_i - 1, n_i);
        
    } else if (k <= n_i - t_pt_length) {
        // 1-5 -> () -> () -> 1-5
        t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
        
    } else if (k == n_i - t_pt_length + 1) { 
        // 1-5 -> () -> ()
        t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
    }
    
    // The last two time points are Gibbs updates -----------------------------
    if(k >= n_i - t_pt_length) {
        arma::vec gibbs_done = {0,0};
        return gibbs_done;
    }
    // ------------------------------------------------------------------------
    
    // Likelihood computations ------------------------------------------------
    double like_comp_prop = 0;
    double like_comp_curr = 0;
    int jj_start = t_pt_length+1;
    
    if(option_num == 2) { 
        jj_start = 0; // proposal distribution is only transition probabilities
    }
    
    for(int jj = jj_start; jj < t_pts.n_elem; jj++) {
        
        int t_j = t_pts(jj);
        
        // Current ---------
        like_comp_curr = like_comp_curr + R::dnorm(y_i(t_j), mu(curr_b_i(t_j) - 1), 1, true);
        
        // Proposal  ---------
        like_comp_prop = like_comp_prop + R::dnorm(y_i(t_j), mu(pr_b_i(t_j) - 1), 1, true);
    } 
    
    arma::vec mh_like = {like_comp_prop, like_comp_curr};
    return mh_like;
    
}

arma::mat state_prob_MH(const int k, const int n_i, int t_pt_length, 
                        const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                        const arma::vec &y_i, const arma::vec &b_i,
                        arma::mat omega_set) {
    
    // Parameter initialization ------------------------------------------------
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::vec P_init = {0.5,0.5}; 
    
    arma::vec curr_b_i_k;
    if(k == 1) {
        curr_b_i_k = b_i.rows(0, omega_set.n_cols - 1);
    } else if (k <= n_i - t_pt_length) {
        curr_b_i_k = b_i.rows(k - 1, k + t_pt_length - 2);
    } else if (k == n_i - t_pt_length + 1) {
        curr_b_i_k = b_i.rows(k - 1, k + t_pt_length - 2);
    }
    // ------------------------------------------------------------------------
    
    // Columns: proposal distribution, remaining likelihood, indicator
    arma::mat prob_dist(omega_set.n_rows, 3, arma::fill::zeros);
    
    for(int j = 0; j < omega_set.n_rows; j++) {
        
        // Initialize the new state sequence and the points to evaluate likelihood
        arma::vec t_pts;
        arma::vec ss_j = b_i;
        if (k == 1) {
            // () -> () -> 1-5
            ss_j.rows(0, omega_set.n_cols - 1) = omega_set.row(j).t();
            t_pts = arma::linspace(0, n_i - 1, n_i);
            
        } else if (k <= n_i - t_pt_length) {
            // 1-5 -> () -> () -> 1-5
            ss_j.rows(k - 1, k + t_pt_length - 2) = omega_set.row(j).t();
            t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
            
        } else if (k == n_i - t_pt_length + 1) { 
            // 1-5 -> () -> ()
            ss_j.rows(k - 1, k + t_pt_length - 2) = omega_set.row(j).t();
            t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
        } 
        
        // Check if the new state sequence is the same as b_i -----------------
        arma::vec state_diff = curr_b_i_k - omega_set.row(j).t();
        if(arma::any(state_diff != 0)) {
            prob_dist(j,2) = 0;
        } else {
            prob_dist(j,2) = 1;
        }
        
        // Likelihood computations ---------------------------------------------
        double like_comp_prob = 0;
        double like_comp_resp = 0;
        
        int loop_max = t_pt_length + 1;
        if(k == n_i - t_pt_length + 1) {
            loop_max = t_pt_length;
        }
        
        for(int jj = 0; jj < loop_max; jj++) {
            
            int t_j = t_pts(jj);
            
            // (1) Transition probabilities ------------------------------------
            if(t_j == 0) {
                like_comp_prob = like_comp_prob + log(P_init(ss_j(t_j) - 1));
            } else{ 
                
                arma::vec qz = exp(zeta);
                arma::mat Q = { {    1,   qz(0)},
                                {qz(1),       1}};
                
                arma::vec q_row_sums = arma::sum(Q, 1);
                arma::mat P_i = Q.each_col() / q_row_sums;
                int b_k_1 = ss_j(t_j-1);
                int b_k = ss_j(t_j);
                
                like_comp_prob = like_comp_prob + log(P_i(b_k_1 - 1, b_k - 1));
            } 
            
            // (2) Response outcome --------------------------------------------
            like_comp_resp = like_comp_resp + R::dnorm(y_i(t_j), mu(ss_j(t_j) - 1), 1, true);
        } 
        
        prob_dist(j,0) = like_comp_prob + like_comp_resp;
    }
    
    // For numerical stability, we will scale on the log-scale
    double prob_log_max = prob_dist.col(0).max();
    
    prob_dist.col(0) = prob_dist.col(0) - prob_log_max;
    prob_dist.col(0) = exp(prob_dist.col(0));
    prob_dist.col(0)= (1/arma::accu(prob_dist.col(0))) * prob_dist.col(0);
    
    arma::uvec num_non_zero = arma::find(prob_dist.col(2) == 1);
    if(num_non_zero.n_elem > 1) {
        Rcpp::Rcout << "bad indicator" << std::endl;
    }
    
    return prob_dist;
}
arma::mat state_prob_only(const int k, const int n_i, int t_pt_length, 
                          const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                          const arma::vec &y_i, const arma::vec &b_i,
                          arma::mat omega_set) {
    
    // Parameter initialization ------------------------------------------------
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::vec P_init = {0.5,0.5}; 
    
    arma::vec curr_b_i_k;
    if(k == 1) {
        curr_b_i_k = b_i.rows(0, omega_set.n_cols - 1);
    } else if (k <= n_i - t_pt_length) {
        curr_b_i_k = b_i.rows(k - 1, k + t_pt_length - 2);
    } else if (k == n_i - t_pt_length + 1) {
        curr_b_i_k = b_i.rows(k - 1, k + t_pt_length - 2);
    }
    // ------------------------------------------------------------------------
    
    // Columns: proposal distribution, remaining likelihood, indicator
    arma::mat prob_dist(omega_set.n_rows, 3, arma::fill::zeros);
    
    for(int j = 0; j < omega_set.n_rows; j++) {
        
        // Initialize the new state sequence and the points to evaluate likelihood
        arma::vec t_pts;
        arma::vec ss_j = b_i;
        if (k == 1) {
            // () -> () -> 1-5
            ss_j.rows(0, omega_set.n_cols - 1) = omega_set.row(j).t();
            t_pts = arma::linspace(0, n_i - 1, n_i);
            
        } else if (k <= n_i - t_pt_length) {
            // 1-5 -> () -> () -> 1-5
            ss_j.rows(k - 1, k + t_pt_length - 2) = omega_set.row(j).t();
            t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
            
        } else if (k == n_i - t_pt_length + 1) {  
            // 1-5 -> () -> ()
            ss_j.rows(k - 1, k + t_pt_length - 2) = omega_set.row(j).t();
            t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
        }  
        
        // Check if the new state sequence is the same as b_i -----------------
        arma::vec state_diff = curr_b_i_k - omega_set.row(j).t();
        if(arma::any(state_diff != 0)) {
            prob_dist(j,2) = 0;
        } else { 
            prob_dist(j,2) = 1;
        } 
        
        // State transition probability ---------------------------------------
        double like_comp_prob = 0;
        
        int loop_max = t_pt_length + 1;
        if(k == n_i - t_pt_length + 1) {
            loop_max = t_pt_length;
        } 
        
        for(int jj = 0; jj < loop_max; jj++) {
            
            int t_j = t_pts(jj);
            
            // (1) Transition probabilities ------------------------------------
            if(t_j == 0) {
                like_comp_prob = like_comp_prob + log(P_init(ss_j(t_j) - 1));
            } else{ 
                
                arma::vec qz = exp(zeta);
                arma::mat Q = { {    1,   qz(0)},
                                {qz(1),       1}};
                
                arma::vec q_row_sums = arma::sum(Q, 1);
                arma::mat P_i = Q.each_col() / q_row_sums;
                int b_k_1 = ss_j(t_j-1);
                int b_k = ss_j(t_j);
                
                like_comp_prob = like_comp_prob + log(P_i(b_k_1 - 1, b_k - 1));
            } 
        } 
        
        prob_dist(j,0) = like_comp_prob;
    } 
    
    // For numerical stability, we will scale on the log-scale
    double prob_log_max = prob_dist.col(0).max();
    
    prob_dist.col(0) = prob_dist.col(0) - prob_log_max;
    prob_dist.col(0) = exp(prob_dist.col(0));
    prob_dist.col(0)= (1/arma::accu(prob_dist.col(0))) * prob_dist.col(0);
    arma::uvec num_non_zero = arma::find(prob_dist.col(2) == 1);
    if(num_non_zero.n_elem > 1) {
        Rcpp::Rcout << "bad indicator" << std::endl;
    }
    
    return prob_dist;
}

// [[Rcpp::export]]
arma::field<arma::vec> update_b_i_MH(const arma::vec EIDs, const arma::vec &par, 
                                     const arma::field<arma::uvec> &par_index, 
                                     arma::field <arma::vec> &B, 
                                     const arma::vec &y, const arma::vec &eids, 
                                     int n_cores, int t_pt_length, int option_num) {
    
    // par_index KEY: (0) mu, (1) zeta; 
    // "i" is the numeric EID number; "ii" is the index of the EID
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::vec B_temp = B(ii);
        arma::vec y_i = y.rows(sub_ind);
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - (t_pt_length - 1); k++) {
            
            arma::vec pr_B = B_temp;
            
            // All possible state transitions given current time point ---------
            arma::mat Omega_set = Omega_fun_cpp_new_multi(k + 1, n_i, B_temp, t_pt_length);
            
            // Learn the proposal distribution ---------------------------------
            int row_ind;
            if(option_num == 1) {
                row_ind = arma::randi(arma::distr_param(1, Omega_set.n_rows));
            } else { 
                arma::mat ss_prob;
                if(option_num == 2) {
                    ss_prob = state_prob_only(k+1, n_i, t_pt_length, par,
                                            par_index, y_i, B_temp, Omega_set);
                } else if(option_num == 3) {
                    ss_prob = state_prob_MH(k+1, n_i, t_pt_length, par,
                                            par_index, y_i, B_temp, Omega_set);
                }
                
                arma::vec x_sample = arma::linspace(1, Omega_set.n_rows, Omega_set.n_rows);
                arma::vec row_ind_vec = RcppArmadillo::sample(x_sample, 1, false, ss_prob.col(0));
                row_ind = row_ind_vec(0);
            }
            
            pr_B.rows(k, k+t_pt_length-1) = Omega_set.row(row_ind-1).t();
            
            // log MH-ratio ----------------------------------------------------
            arma::vec mh_like_val;
            if(option_num == 1) {
                mh_like_val = full_like_MH(k+1, n_i, t_pt_length, par, 
                                           par_index, y_i, pr_B, B_temp);
            } else {
                mh_like_val = sub_like_MH(k+1, n_i, t_pt_length, par, 
                                          par_index, y_i, pr_B, B_temp, option_num);   
            }
            
            double diff_check = mh_like_val(0) - mh_like_val(1);
            double min_log = log(arma::randu(arma::distr_param(0,1)));
            if(diff_check > min_log){B_temp = pr_B;} 
        }
        B_return(ii) = B_temp;
    } 
    
    return B_return;
} 

// [[Rcpp::export]]
arma::vec state_prob_gibbs(const int k, const int n_i, int t_pt_length, 
                           const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                           const arma::mat &y_i, const arma::vec &b_i,
                           arma::mat omega_set) {
    
    // Parameter initialization ------------------------------------------------
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::vec P_init = {0.5,0.5}; 
    // ------------------------------------------------------------------------
    
    arma::vec prob_dist(omega_set.n_rows, arma::fill::ones);
    
    for(int j = 0; j < omega_set.n_rows; j++) {
        
        // Initialize the new state sequence and the points to evaluate likelihood
        arma::vec t_pts;
        arma::vec ss_j = b_i;
        if (k == 1) {
            // () -> () -> 1-5
            ss_j.rows(0, omega_set.n_cols - 1) = omega_set.row(j).t();
            t_pts = arma::linspace(0, n_i - 1, n_i);
            
        } else if (k <= n_i - t_pt_length) {
            // 1-5 -> () -> () -> 1-5
            ss_j.rows(k - 1, k + t_pt_length - 2) = omega_set.row(j).t();
            t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
            
        } else if (k == n_i - t_pt_length + 1) { 
            // 1-5 -> () -> ()
            ss_j.rows(k - 1, k + t_pt_length - 2) = omega_set.row(j).t();
            t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
        } 
        
        // Likelihood computations ---------------------------------------------
        double like_comp_prob = 0;
        double like_comp_resp = 0;
        for(int jj = 0; jj < t_pts.n_elem; jj++) {
            
            int t_j = t_pts(jj);
            
            // (1) Transition probabilities ------------------------------------
            if(jj < t_pt_length + 1) {
                if(t_j == 0) {
                    like_comp_prob = like_comp_prob + log(P_init(ss_j(t_j) - 1));
                } else{ 
                    
                    arma::vec qz = exp(zeta);
                    arma::mat Q = { {    1,   qz(0)},
                                    {qz(1),       1}};
                    
                    arma::vec q_row_sums = arma::sum(Q, 1);
                    arma::mat P_i = Q.each_col() / q_row_sums;
                    int b_k_1 = ss_j(t_j-1);
                    int b_k = ss_j(t_j);
                    
                    like_comp_prob = like_comp_prob + log(P_i(b_k_1 - 1, b_k - 1));
                }    
            }
            
            // (2) Response outcome --------------------------------------------
            like_comp_resp = like_comp_resp + R::dnorm(y_i(t_j), mu(ss_j(t_j) - 1), 1, true);
        } 
        
        prob_dist(j) = like_comp_prob + like_comp_resp;
    }
    
    // For numerical stability, we will scale on the log-scale
    double prob_log_max = prob_dist.max();
    prob_dist = prob_dist - prob_log_max;
    prob_dist = exp(prob_dist);
    prob_dist= (1/arma::accu(prob_dist)) * prob_dist;
    
    return prob_dist;
}

// [[Rcpp::export]]
arma::field<arma::vec> update_b_i_gibbs( const arma::vec EIDs, const arma::vec &par, 
                                         const arma::field<arma::uvec> &par_index, 
                                         arma::field <arma::vec> &B, 
                                         const arma::vec &y, const arma::vec &eids,
                                         int n_cores, int t_pt_length) {
    
    // par_index KEY: (0) mu, (1) zeta; 
    // "i" is the numeric EID number; "ii" is the index of the EID
    arma::field<arma::vec> B_return(EIDs.n_elem);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::vec B_temp = B(ii);
        arma::vec y_i = y.rows(sub_ind);
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - (t_pt_length - 1); k++) {
            
            // All possible state transitions given current time point ---------
            arma::mat Omega_set = Omega_fun_cpp_new_multi(k + 1, n_i, B_temp, t_pt_length);
            
            // Learn the proposal distribution ---------------------------------
            arma::vec ss_prob = state_prob_gibbs(k+1, n_i, t_pt_length, par,
                                                 par_index, y_i, B_temp, Omega_set);
            
            arma::vec x_sample = arma::linspace(1, Omega_set.n_rows, Omega_set.n_rows);
            arma::vec row_ind = RcppArmadillo::sample(x_sample, 1, false, ss_prob);
            
            // Gibbs update ----------------------------------------------------
            B_temp.rows(k, k+t_pt_length-1) = Omega_set.row(row_ind(0)-1).t();
            
        }
        B_return(ii) = B_temp;
    } 
    
    return B_return;
} 

double log_f_i_cpp_total(const arma::vec &EIDs, const arma::vec &par, 
                         const arma::field<arma::uvec> &par_index, 
                         const arma::field <arma::vec> &B, const arma::vec &y, 
                         const arma::vec &eids, int n_cores) {
    
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::vec P_init = {0.5,0.5}; 

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
        
        for(int jj = 0; jj < n_i; jj++) {
            
            int t_j = jj;
            
            // (1) Transition probabilities ------------------------------------
            if(t_j == 0) {
                like_comp_transition = like_comp_transition + log(P_init(b_i(t_j) - 1));
            } else{ 
                arma::vec qz = exp(zeta);
                arma::mat Q = { {    1,   qz(0)},
                                {qz(1),       1}};
                
                arma::vec q_row_sums = arma::sum(Q, 1);
                arma::mat P_i = Q.each_col() / q_row_sums;
                
                int curr_b_k_1 = b_i(t_j-1);
                int curr_b_k   = b_i(t_j);
                
                like_comp_transition = like_comp_transition + 
                    log(P_i(curr_b_k_1 - 1, curr_b_k - 1));
            } 
            
            // (2) Likelihood component ------------------------------------
            like_comp = like_comp + R::dnorm(y_i(t_j), mu(b_i(t_j) - 1), 1, true);
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
double log_post_cpp_no_b(const arma::vec &EIDs, const arma::vec &par, 
                         const arma::field<arma::uvec> &par_index,
                         const arma::vec &y, const arma::vec &eids, int n_cores) {
    // Compute the likelihood while integrating out the states
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::rowvec P_init = {0.5,0.5}; 
    
    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};
                    
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        arma::rowvec f_i(P_init.n_elem, arma::fill::ones);
        arma::rowvec val(P_init.n_elem, arma::fill::ones);
        double log_norm = 0;
        
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::vec y_i = y.rows(sub_ind);
        
        for(int jj = 0; jj < n_i; jj++) {
            
            if(jj == 0) {
                arma::vec resp_like(P_init.n_elem, arma::fill::zeros);
                for(int ll = 0; ll < P_init.n_elem; ll++) {
                    resp_like(ll) = R::dnorm(y_i(jj), mu(ll), 1, false);
                }
                arma::mat D = arma::diagmat(resp_like);
                f_i = P_init * D;
            } else{ 
                
                arma::vec resp_like(P_init.n_elem, arma::fill::zeros);
                for(int ll = 0; ll < P_init.n_elem; ll++) {
                    resp_like(ll) = R::dnorm(y_i(jj), mu(ll), 1, false);
                }
                arma::mat D = arma::diagmat(resp_like);
                
                val = f_i * P_i * D;
            } 
            
            double norm_val = arma::norm(val, 2);
            f_i = val / norm_val;
            log_norm = log_norm + log(norm_val);
        } 
        
        in_vals(ii) = log(arma::accu(f_i)) + log_norm;
    } 
    
    double in_value = arma::accu(in_vals);
    
    
    // Prior densities
    arma::vec prior_mean(par.n_elem, arma::fill::zeros);
    
    arma::vec prior_var_diag(par.n_elem, arma::fill::ones);
    prior_var_diag = 20 * prior_var_diag;
    arma::mat prior_sd = arma::diagmat(prior_var_diag);
    
    arma::vec prior_dens = dmvnorm(par.t(), prior_mean, prior_sd, true);
    double prior_dens_val = arma::as_scalar(prior_dens);
    
    in_value = in_value + prior_dens_val;
    
    return in_value;
}

// [[Rcpp::export]]
double pseudo_like(const arma::vec &EIDs, const arma::vec &par, 
                      const arma::field<arma::uvec> &par_index,
                      const arma::vec &y, const arma::vec &eids, int n_cores) {
    
    arma::field<arma::vec> B_mle(EIDs.n_elem);
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::vec P_init = {0.5,0.5}; 
    
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
        double log_like_i = 0;
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i; k++) {
            arma::vec like_vals(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
            if(k == 0) {
                // Initial state
                for(int jj = 0; jj < like_vals.n_elem; jj++) {
                    like_vals(jj) = log(P_init(jj)) + R::dnorm(y_i(k), mu(jj), 1, true);
                }
                
                b_i(k) = arma::index_max(like_vals) + 1;
                log_like_i = log_like_i + arma::max(like_vals);
            } else {
                // All others
                for(int jj = 0; jj < like_vals.n_elem; jj++) {
                    like_vals(jj) = log(P_i(b_i(k-1) - 1, jj)) + R::dnorm(y_i(k), mu(jj), 1, true);
                }
                
                b_i(k) = arma::index_max(like_vals) + 1;
                log_like_i = log_like_i + arma::max(like_vals);
            }
        }
        
        B_mle(ii) = b_i;
        in_vals(ii) = log_like_i;
    } 
    
    double in_value = arma::accu(in_vals);
    
    // Prior densities
    arma::vec prior_mean(par.n_elem, arma::fill::zeros);
    
    arma::vec prior_var_diag(par.n_elem, arma::fill::ones);
    prior_var_diag = 20 * prior_var_diag;
    arma::mat prior_sd = arma::diagmat(prior_var_diag);
    
    arma::vec prior_dens = dmvnorm(par.t(), prior_mean, prior_sd, true);
    double prior_dens_val = arma::as_scalar(prior_dens);
    
    in_value = in_value + prior_dens_val;
    
    return in_value;
}

// [[Rcpp::export]]
double pseudo_like2(const arma::vec &EIDs, const arma::vec &par, 
                    const arma::field<arma::uvec> &par_index,
                    const arma::vec &y, const arma::vec &eids, int n_cores) {
    
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    arma::vec mu = par.elem(par_index(0) - 1);
    arma::vec zeta = par.elem(par_index(1) - 1);
    arma::vec P_init = {0.5,0.5}; 
    
    arma::vec qz = exp(zeta);
    arma::mat Q = { {    1,   qz(0)},
                    {qz(1),       1}};
    
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;
    
    arma::vec state_samp = arma::linspace(1, adj_mat_GLOBAL.n_cols, adj_mat_GLOBAL.n_cols);
    
    // Find the MLE state sequence given the current parameters ----------------
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        arma::vec y_i = y.rows(sub_ind);
        
        int M = 20;
        arma::vec emp_mean_vec(M, arma::fill::zeros);
        
        // Empirical mean estimate ---------------------------------------------
        for(int m = 0; m < M; m++) {
            
            arma::vec b_i(n_i, arma::fill::zeros);
            double log_like_i = 0;
            double log_q_b_i = 0;
            
            // Looping through subject state space -----------------------------
            for (int k = 0; k < n_i; k++) {
                arma::vec like_vals(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
                if(k == 0) {
                    // Initial state
                    for(int jj = 0; jj < like_vals.n_elem; jj++) {
                        like_vals(jj) = P_init(jj) * R::dnorm(y_i(k), mu(jj), 1, false);
                    }

                    arma::vec sample_prob(adj_mat_GLOBAL.n_cols, arma::fill::ones);
                    sample_prob = like_vals / arma::accu(like_vals);
                    // arma::vec sample_prob(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
                    // sample_prob.elem(arma::find(P_init > 0)) += 1;
                    // sample_prob = sample_prob / arma::accu(sample_prob);
                    
                    arma::vec row_ind_vec = RcppArmadillo::sample(state_samp, 1, false, sample_prob);
                    b_i(k) = row_ind_vec(0);
                    
                    log_like_i = log_like_i + log(P_init(b_i(k) - 1)) + 
                        R::dnorm(y_i(k), mu(b_i(k) - 1), 1, true);
                    log_q_b_i = log_q_b_i + log(sample_prob(b_i(k) - 1));
                } else {
                    // All others
                    for(int jj = 0; jj < like_vals.n_elem; jj++) {
                        like_vals(jj) = P_i(b_i(k-1) - 1, jj) * R::dnorm(y_i(k), mu(jj), 1, false);
                    }

                    arma::vec sample_prob(adj_mat_GLOBAL.n_cols, arma::fill::ones);
                    sample_prob = like_vals / arma::accu(like_vals);
                    // arma::ivec possible_states = adj_mat_GLOBAL.row(b_i(k-1) - 1).t();
                    // arma::vec sample_prob(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
                    // sample_prob.elem(arma::find(possible_states > 0)) += 1;
                    // sample_prob = sample_prob / arma::accu(sample_prob);
                    
                    arma::vec row_ind_vec = RcppArmadillo::sample(state_samp, 1, false, sample_prob);
                    b_i(k) = row_ind_vec(0);
                    
                    log_like_i = log_like_i + log(P_i(b_i(k-1) - 1, b_i(k) - 1)) +
                        R::dnorm(y_i(k), mu(b_i(k) - 1), 1, true);
                    log_q_b_i = log_q_b_i + log(sample_prob(b_i(k) - 1));
                }
            }
            emp_mean_vec(m) = log_like_i - log_q_b_i;
        }
        
        double max_log_emp = emp_mean_vec.max();
        emp_mean_vec = emp_mean_vec - max_log_emp;
        arma::vec exp_emp_mean_vec = exp(emp_mean_vec);
        
        in_vals(ii) = log(arma::mean(exp_emp_mean_vec)) + max_log_emp;
    }
    
    double in_value = arma::accu(in_vals);
    
    // Prior densities
    arma::vec prior_mean(par.n_elem, arma::fill::zeros);
    
    arma::vec prior_var_diag(par.n_elem, arma::fill::ones);
    prior_var_diag = 20 * prior_var_diag;
    arma::mat prior_sd = arma::diagmat(prior_var_diag);
    
    arma::vec prior_dens = dmvnorm(par.t(), prior_mean, prior_sd, true);
    double prior_dens_val = arma::as_scalar(prior_dens);
    
    in_value = in_value + prior_dens_val;
    
    return in_value;
} 

// [[Rcpp::export]]
void test_fnc() {
    Rcpp::Rcout << Omega_List_GLOBAL_multi << std::endl;
    Rcpp::Rcout << adj_mat_GLOBAL << std::endl;
    
    Rcpp::Rcout << Omega_List_GLOBAL_multi(1) << std::endl;
    Rcpp::Rcout << Omega_List_GLOBAL_multi(1)(0, 1) << std::endl;
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
    arma::vec P_init = {0.5,0.5};

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
                    
                    for(int kk = k; kk < k + states_per_step + 1; kk++) {
                        if(kk == 0) {
                            log_like_val = log_like_val + 
                                           log(P_init(ss_jj(kk) - 1)) + 
                                           R::dnorm(y_i(kk), mu(ss_jj(kk) - 1), 1, true);
                        } else if(kk <= n_i - states_per_step) {
                            log_like_val = log_like_val + 
                                           log(P_i(ss_jj(kk-1) - 1, ss_jj(kk) - 1)) + 
                                           R::dnorm(y_i(kk), mu(ss_jj(kk) - 1), 1, true);
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
    arma::vec P_init = {0.5,0.5};

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
                    
                for(int kk = k; kk < k + states_per_step + 1; kk++) {
                    if(kk == 0) {
                        log_like_b = log_like_b + 
                            log(P_init(b_i(kk) - 1)) + 
                            R::dnorm(y_i(kk), mu(b_i(kk) - 1), 1, true);
                        
                        log_like_s = log_like_s + 
                            log(P_init(s_i(kk) - 1)) + 
                            R::dnorm(y_i(kk), mu(s_i(kk) - 1), 1, true);
                        
                    } else if(kk <= n_i - states_per_step) {
                        log_like_b = log_like_b + 
                            log(P_i(b_i(kk-1) - 1, b_i(kk) - 1)) + 
                            R::dnorm(y_i(kk), mu(b_i(kk) - 1), 1, true);
                        
                        log_like_s = log_like_s + 
                            log(P_i(s_i(kk-1) - 1, s_i(kk) - 1)) + 
                            R::dnorm(y_i(kk), mu(s_i(kk) - 1), 1, true);
                    }
                }
                double diff_check = log_like_s - log_like_b;
                double min_log = log(arma::randu(arma::distr_param(0,1)));
                if(diff_check > min_log){b_i = s_i;} 
            }
        }
        
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
    arma::vec P_init = {0.5,0.5}; 
    
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
            arma::vec like_vals(adj_mat_GLOBAL.n_cols, arma::fill::zeros);
            if(k == 0) {
                // Initial state
                for(int jj = 0; jj < like_vals.n_elem; jj++) {
                    like_vals(jj) = log(P_init(jj)) + R::dnorm(y_i(k), mu(jj), 1, true);
                } 
                
                b_i(k) = arma::index_max(like_vals) + 1;
            } else {
                // All others
                for(int jj = 0; jj < like_vals.n_elem; jj++) {
                    like_vals(jj) = log(P_i(b_i(k-1) - 1, jj)) + R::dnorm(y_i(k), mu(jj), 1, true);
                }
            
                b_i(k) = arma::index_max(like_vals) + 1;
            }
        }
        
        B_mle(ii) = b_i;
    } 
    
    
    return B_mle;
}