#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

NumericVector csample_num( NumericVector x, int size, bool replace,
                           NumericVector prob = NumericVector::create()) {
    NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
    return ret;
}

// Combinatoric functions for finding possible state sequences -----------------
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

arma::field<arma::field<arma::mat>> get_Omega_list(const arma::imat &adj_mat, int s) {
    
    int N = adj_mat.n_cols; // dimension of adj matrix
    
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
        Rcpp::List a_temp_list = getAllStateSequences_forward(i, adj_mat, s+1);
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
        Rcpp::List c_temp_list = getAllStateSequences_backward(i, adj_mat, s+1);
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
            Rcpp::List b_temp_list = getAllStateSequences_both(i, j, adj_mat, s+2);
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
arma::imat adj_mat_sub_GLOBAL;

arma::field<arma::field<arma::mat>> Omega_List_GLOBAL_multi;
arma::field<arma::field<arma::mat>> Omega_List_GLOBAL_sub_multi;

int states_per_step_GLOBAL;

// [[Rcpp::export]]
void initialize_cpp(arma::imat a_mat, arma::imat a_mat_sub, int s_per_s) {
    adj_mat_GLOBAL = a_mat;
    adj_mat_sub_GLOBAL = a_mat_sub;
    states_per_step_GLOBAL = s_per_s;
    
    Omega_List_GLOBAL_multi = get_Omega_list(adj_mat_GLOBAL, states_per_step_GLOBAL);
    Omega_List_GLOBAL_sub_multi = get_Omega_list(adj_mat_sub_GLOBAL, states_per_step_GLOBAL);
}

arma::mat Omega_fun_cpp_new(const int k, const int n_i, const arma::vec &b_i,
                            const bool sub) {
    
    arma::mat Omega_set;
    // We are in a reduced dimension case
    if(sub) { 
        // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
        if (k == 1){
            // () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_sub_multi(0)(b_i(2) - 1);
        } else if (k <= n_i - 2) {
            // 1-5 -> () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_sub_multi(1)(b_i(k - 2) - 1, b_i(k + 1) - 1);
        } else if (k == n_i - 1) {
            // 1-5 -> () -> ()
            Omega_set = Omega_List_GLOBAL_sub_multi(2)(b_i(n_i - 3) - 1);
        }
    } else {
        // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
        if (k == 1) {
            // () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_multi(0)(b_i(2) - 1);
        } else if (k <= n_i - 2) {
            // 1-5 -> () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_multi(1)(b_i(k - 2) - 1, b_i(k + 1) - 1);
        } else if (k == n_i - 1) {
            // 1-5 -> () -> ()
            Omega_set = Omega_List_GLOBAL_multi(2)(b_i(n_i - 3) - 1);
        }
    }
    
    return Omega_set;
    
}

arma::mat Omega_fun_cpp_new_multi(const int k, const int n_i, const arma::vec &b_i,
                                  const bool sub, int t_pt_length) {
    
    arma::mat Omega_set;
    // We are in a reduced dimension case
    if(sub) { 
        // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
        if (k == 1){
            // () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_sub_multi(0)(b_i(t_pt_length) - 1);
        } else if (k <= n_i - t_pt_length) {
            // 1-5 -> () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_sub_multi(1)(b_i(k - 2) - 1, 
                                                    b_i(k + t_pt_length - 1) - 1);
        } else if (k == n_i - t_pt_length + 1) {
            // 1-5 -> () -> ()
            Omega_set = Omega_List_GLOBAL_sub_multi(2)(b_i(n_i - t_pt_length - 1) - 1);
        }
    } else {
        // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
        if (k == 1) {
            // () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_multi(0)(b_i(t_pt_length) - 1);
        } else if (k <= n_i - t_pt_length) {
            // 1-5 -> () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_multi(1)(b_i(k - 2) - 1, 
                                                b_i(k + t_pt_length - 1) - 1);
        } else if (k == n_i - t_pt_length + 1) {
            // 1-5 -> () -> ()
            Omega_set = Omega_List_GLOBAL_multi(2)(b_i(n_i - t_pt_length - 1) - 1);
        }
    }
    
    return Omega_set;
    
}

arma::mat get_omega_list(const int k, const int n_i, const arma::vec &b_i, 
                         int states_per_step, const bool sub) {
    
    arma::mat Omega_set;
    
    if(sub) {
        // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
        if (k == 0) {
            // () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_sub_multi(0)(b_i(k + states_per_step) - 1);
        } else if (k == n_i - states_per_step) {
            // 1-5 -> () -> ()
            Omega_set = Omega_List_GLOBAL_sub_multi(2)(b_i(n_i - states_per_step - 1) - 1);
        } else {
            // 1-5 -> () -> () -> 1-5
            Omega_set = Omega_List_GLOBAL_sub_multi(1)(b_i(k - 1) - 1, 
                                                        b_i(k + states_per_step) - 1);
        }
    } else {
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
    }
    
    return Omega_set;
    
}
// -----------------------------------------------------------------------------

// *** Using ***
// [[Rcpp::export]]
double log_f_i_cpp(const int i, const int ii, arma::vec t_pts, const arma::vec &par, 
                   const arma::field<arma::uvec> &par_index, const arma::vec &A, 
                   const arma::vec &B, const arma::mat &Y, const arma::mat &z, 
                   const arma::field<arma::mat> &Dn, const arma::field<arma::mat> &Xn, 
                   const arma::field<arma::mat> &Dn_omega, const arma::vec &W) {
  
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    double in_value = 0;
    
    // Initializing parameters ---------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                          exp(vec_init_content(2)), exp(vec_init_content(3))}; // THREE STATE
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    // ---------------------------------------------------------------------------
    
    arma::vec eids = Y.col(0);
    arma::uvec sub_ind = arma::find(eids == i);
    int n_i = sub_ind.n_elem;
    
    arma::mat Y_temp = Y.rows(sub_ind);
    arma::mat Y_i = Y_temp.cols(1, 4);
    Y_i = Y_i.t();
    arma::mat z_i = z.rows(sub_ind);
    
    arma::vec b_i = B;
    
    arma::field<arma::mat> Dn_alpha_full = Dn;
    arma::field<arma::mat> Dn_omega_full = Dn_omega;
    arma::field<arma::mat> Xn_full = Xn;
    arma::vec vec_alpha_ii = A;
    arma::vec vec_omega_ii = W;
    
    if (any(t_pts == -1)) { t_pts = arma::linspace(1, n_i, n_i);}
    
    for(int w = 0; w < t_pts.n_elem; w++){
        
        int k = t_pts(w);
        
        if(k==1){
            
            int b_0 = b_i(0);
            
            arma::vec vec_A = A_all_state.col(b_0 - 1);
            
            arma::mat Gamma = { {R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                 R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                 R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                 R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                 {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                  R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                  R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                  R(1,3) / (1 - vec_A(1) * vec_A(3))},
                                  {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                   R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                   R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                   R(2,3) / (1 - vec_A(2) * vec_A(3))},
                                   {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                    R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                    R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                    R(3,3) / (1 - vec_A(3) * vec_A(3))}};
          
            arma::vec y_1 = Y_i.col(0);
            arma::vec nu_1 = Dn_alpha_full(0) * vec_alpha_ii + 
                                Dn_omega_full(0) * vec_omega_ii +
                                  Xn_full(0) * vec_beta;
            arma::vec log_y_pdf = dmvnorm(y_1.t(), nu_1, Gamma, true);
            
            in_value = in_value + log(P_init(b_0 - 1)) + arma::as_scalar(log_y_pdf);
            
        } else{
            
            double q1_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(0));
            double q1 = exp(q1_sub);
            double q2_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(1));
            double q2 = exp(q2_sub);
            double q3_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(2));
            double q3 = exp(q3_sub);
            double q4_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(3));
            double q4 = exp(q4_sub);
            
            double q5_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(4));
            double q5 = exp(q5_sub);
            double q6_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(5));
            double q6 = exp(q6_sub);
            double q7_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(6));
            double q7 = exp(q7_sub);
            double q8_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(7));
            double q8 = exp(q8_sub);
            
            double q9_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(8));
            double q9 = exp(q9_sub);
            double q10_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(9));
            double q10 = exp(q10_sub);
            double q11_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(10));
            double q11 = exp(q11_sub);
            double q12_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(11));
            double q12 = exp(q12_sub);
            
            arma::mat Q = {{   1,   q1,  0,  q2,  0},
                           {   0,    1, q3,  q4,  0},
                           {  q5,   q6,  1,  q7,  0},
                           {   0,   q8,  0,   1, q9},
                           { q10,  q11,  0, q12,  1}}; // THREE STATE
            
            arma::vec q_row_sums = arma::sum(Q, 1);
            arma::mat P_i = Q.each_col() / q_row_sums;
            
            int b_k_1 = b_i(k-2);
            int b_k = b_i(k-1);
            
            arma::vec vec_A = A_all_state.col(b_k - 1);
            arma::mat A_1 = arma::diagmat(vec_A);
            
            arma::vec y_k_1 = Y_i.col(k-2);
            arma::vec y_k = Y_i.col(k-1);
            arma::vec nu_k_1 = Dn_alpha_full(k-2) * vec_alpha_ii + 
                                    Dn_omega_full(k-2) * vec_omega_ii +
                                        Xn_full(k-2) * vec_beta;
            arma::vec nu_k = Dn_alpha_full(k-1) * vec_alpha_ii + 
                                Dn_omega_full(k-1) * vec_omega_ii +
                                    Xn_full(k-1) * vec_beta;
            
            arma::vec mean_k = nu_k + A_1 * (y_k_1 - nu_k_1);
            
            arma::vec log_y_k_pdf = dmvnorm(y_k.t(), mean_k, R, true);
            
            in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1)) + arma::as_scalar(log_y_k_pdf);
        }
    }
    
    return in_value;
}

// *** Using ***
double log_f_i_cpp_total(const arma::vec &EIDs, arma::vec t_pts, const arma::vec &par, 
                         const arma::field<arma::uvec> &par_index, 
                         const arma::field <arma::vec> &A, const arma::field <arma::vec> &B, 
                         const arma::mat &Y, const arma::mat &z, 
                         const arma::field<arma::field<arma::mat>> &Dn, 
                         const arma::field <arma::field<arma::mat>> &Xn, 
                         const arma::field<arma::field<arma::mat>> &Dn_omega, 
                         const arma::field <arma::vec> &W, int n_cores) {

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
  //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID
  
  arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);

    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        in_vals(ii) = log_f_i_cpp(i, ii, t_pts, par, par_index, A(ii), B(ii), 
                                  Y, z, Dn(ii), Xn(ii), Dn_omega(ii), W(ii));
    }
    
    double in_value = arma::accu(in_vals);
    
    return in_value;
}

// *** Using ***
// [[Rcpp::export]]
double log_post_cpp(const arma::vec &EIDs, const arma::vec &par, 
                    const arma::field<arma::uvec> &par_index,
                    const arma::field<arma::vec> &A, const arma::field<arma::vec> &B,
                    const arma::mat &Y, const arma::mat &z, 
                    const arma::field<arma::field<arma::mat>> &Dn,
                    const arma::field<arma::field<arma::mat>> &Xn, 
                    const arma::field<arma::field<arma::mat>> &Dn_omega, 
                    const arma::field<arma::vec> &W, int n_cores) {

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Likelihood --------------------------------------------------------------
    double value;
    arma::vec t_pts = {-1};
    value = log_f_i_cpp_total(EIDs, t_pts, par, par_index, A, B, Y, z, Dn, 
                              Xn, Dn_omega, W, n_cores);
    
    // Autocorrelation prior ---------------------------------------------------
    arma::vec vec_A_content = par.elem(par_index(3) - 1);
    
    arma::vec vec_A_mean = {1.5,  1.5,  1.5,  1.5, -1.0, -1.0, -1.0, -1.0,
                            0.1,  0.1,  0.1,  0.1,    0,    0,    0,    0,
                            0.1,  0.1,  0.1,  0.1};
    arma::vec scalar_A(vec_A_content.n_elem, arma::fill::ones);
    scalar_A = 20 * scalar_A;
    arma::mat A_var = arma::diagmat(scalar_A);
    
    arma::vec prior_A = dmvnorm(vec_A_content.t(), vec_A_mean, A_var, true);
    double prior_A_val = arma::as_scalar(prior_A);
    
    // Error-variance prior ----------------------------------------------------
    arma::vec vec_R_content = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R_content, 4, 4);
    
    int nu_R = 1000;
    arma::vec scalar_vec_R = {9, 9, 9, 9};
    scalar_vec_R = (nu_R - 4 - 1) * scalar_vec_R;
    arma::mat Psi_R = arma::diagmat(scalar_vec_R);
    
    double prior_R_val = diwish(R, nu_R, Psi_R, true);
    
    // Zeta prior --------------------------------------------------------------
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    
    arma::vec vec_zeta_mean = {-4.7405, 4.5, -5.2152,   1, -3.6473,-0.5, -3.1475, -0.2, 
                               -6.4459,  -1, -3.9404,   2, -4.2151,   1, -4.1778, 2.5, 
                               -3.0523,   0, -6.4459,-0.2, -4.2404, 3.5, -4.2151,   1};
    arma::vec scalar_zeta(vec_zeta_mean.n_elem, arma::fill::ones);
    scalar_zeta = 20 * scalar_zeta;
    arma::mat zeta_var = arma::diagmat(scalar_zeta);
    
    arma::vec prior_zeta = dmvnorm(vec_zeta_content.t(), vec_zeta_mean, zeta_var, true);
    double prior_zeta_val = arma::as_scalar(prior_zeta);
    
    // Initial Probability prior -----------------------------------------------
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    
    arma::vec vec_init_mean = {0, 0, 0, 0}; 
    arma::vec scalar_init(vec_init_content.n_elem, arma::fill::ones); 
    scalar_init = 20 * scalar_init;
    arma::mat init_var = arma::diagmat(scalar_init);
    
    arma::vec prior_init = dmvnorm(vec_init_content.t(), vec_init_mean, init_var, true);
    double prior_init_val = arma::as_scalar(prior_init);
    
    // Upsilon omega priors ----------------------------------------------------
    arma::vec vec_up_omega_content = par.elem(par_index(8) - 1);
    arma::vec omega_mean(vec_up_omega_content.n_elem, arma::fill::zeros);
    
    arma::vec scalar_omega(vec_up_omega_content.n_elem, arma::fill::ones);
    scalar_omega = 20 * scalar_omega;
    arma::mat omega_var = arma::diagmat(scalar_omega);
    
    arma::vec prior_omega = dmvnorm(vec_up_omega_content.t(), omega_mean, omega_var, true);
    double prior_omega_val = arma::as_scalar(prior_omega);
    
    // Full log-posterior ------------------------------------------------------
    value = value + prior_A_val + prior_R_val + prior_zeta_val + prior_init_val + prior_omega_val;
    
    return value;
}

// [[Rcpp::export]]
arma::vec full_like_MH(const int k, const int n_i, int t_pt_length, const arma::vec &par, 
                       const arma::field<arma::uvec> &par_index, const arma::mat &y_i, 
                       const arma::mat &z_i, const arma::vec &alpha_i, 
                       const arma::field<arma::mat> &x_i, 
                       const arma::field<arma::mat> &Dn_omega, const arma::vec &w_i,
                       const arma::field<arma::mat> &pr_Dn, arma::vec &pr_b_i,
                       const arma::field<arma::mat> &curr_Dn, arma::vec &curr_b_i) {
    // Parameter initialization ------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    
    arma::vec vec_A_scale = { (exp(vec_A_total(0)) - 1) / (1+exp(vec_A_total(0))),
                              (exp(vec_A_total(1)) - 1) / (1+exp(vec_A_total(1))),
                              (exp(vec_A_total(2)) - 1) / (1+exp(vec_A_total(2))),
                              (exp(vec_A_total(3)) - 1) / (1+exp(vec_A_total(3))),
                              (exp(vec_A_total(4)) - 1) / (1+exp(vec_A_total(4))),
                              (exp(vec_A_total(5)) - 1) / (1+exp(vec_A_total(5))),
                              (exp(vec_A_total(6)) - 1) / (1+exp(vec_A_total(6))),
                              (exp(vec_A_total(7)) - 1) / (1+exp(vec_A_total(7))),
                              (exp(vec_A_total(8)) - 1) / (1+exp(vec_A_total(8))),
                              (exp(vec_A_total(9)) - 1) / (1+exp(vec_A_total(9))),
                              (exp(vec_A_total(10)) - 1) / (1+exp(vec_A_total(10))),
                              (exp(vec_A_total(11)) - 1) / (1+exp(vec_A_total(11))),
                              (exp(vec_A_total(12)) - 1) / (1+exp(vec_A_total(12))),
                              (exp(vec_A_total(13)) - 1) / (1+exp(vec_A_total(13))),
                              (exp(vec_A_total(14)) - 1) / (1+exp(vec_A_total(14))),
                              (exp(vec_A_total(15)) - 1) / (1+exp(vec_A_total(15))),
                              (exp(vec_A_total(16)) - 1) / (1+exp(vec_A_total(16))),
                              (exp(vec_A_total(17)) - 1) / (1+exp(vec_A_total(17))),
                              (exp(vec_A_total(18)) - 1) / (1+exp(vec_A_total(18))),
                              (exp(vec_A_total(19)) - 1) / (1+exp(vec_A_total(19)))};
    
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    
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
                // State space component
                double q1_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(0));
                double q1 = exp(q1_sub);
                double q2_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(1));
                double q2 = exp(q2_sub);
                double q3_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(2));
                double q3 = exp(q3_sub);
                double q4_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(3));
                double q4 = exp(q4_sub);
                
                double q5_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(4));
                double q5 = exp(q5_sub);
                double q6_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(5));
                double q6 = exp(q6_sub);
                double q7_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(6));
                double q7 = exp(q7_sub);
                double q8_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(7));
                double q8 = exp(q8_sub);
                
                double q9_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(8));
                double q9 = exp(q9_sub);
                double q10_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(9));
                double q10 = exp(q10_sub);
                double q11_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(10));
                double q11 = exp(q11_sub);
                double q12_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(11));
                double q12 = exp(q12_sub);
                
                arma::mat Q = { {   1,   q1,  0,  q2,  0},
                                {   0,    1, q3,  q4,  0},
                                {  q5,   q6,  1,  q7,  0},
                                {   0,   q8,  0,   1, q9},
                                { q10,  q11,  0, q12,  1}}; 
                
                arma::vec q_row_sums = arma::sum(Q, 1);
                arma::mat P_i = Q.each_col() / q_row_sums;
                
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
        if(t_j == 0) {
            
            arma::vec y_1 = y_i.col(t_j);
            
            // Current ---------
            arma::vec vec_A_curr = A_all_state.col(curr_b_i(t_j) - 1);
            
            arma::mat Gamma_curr = {{R(0,0) / (1 - vec_A_curr(0) * vec_A_curr(0)), 
                                     R(0,1) / (1 - vec_A_curr(0) * vec_A_curr(1)), 
                                     R(0,2) / (1 - vec_A_curr(0) * vec_A_curr(2)), 
                                     R(0,3) / (1 - vec_A_curr(0) * vec_A_curr(3))},
                                     {R(1,0) / (1 - vec_A_curr(1) * vec_A_curr(0)), 
                                      R(1,1) / (1 - vec_A_curr(1) * vec_A_curr(1)), 
                                      R(1,2) / (1 - vec_A_curr(1) * vec_A_curr(2)), 
                                      R(1,3) / (1 - vec_A_curr(1) * vec_A_curr(3))},
                                      {R(2,0) / (1 - vec_A_curr(2) * vec_A_curr(0)), 
                                       R(2,1) / (1 - vec_A_curr(2) * vec_A_curr(1)), 
                                       R(2,2) / (1 - vec_A_curr(2) * vec_A_curr(2)), 
                                       R(2,3) / (1 - vec_A_curr(2) * vec_A_curr(3))},
                                       {R(3,0) / (1 - vec_A_curr(3) * vec_A_curr(0)), 
                                        R(3,1) / (1 - vec_A_curr(3) * vec_A_curr(1)), 
                                        R(3,2) / (1 - vec_A_curr(3) * vec_A_curr(2)), 
                                        R(3,3) / (1 - vec_A_curr(3) * vec_A_curr(3))}};
            
            arma::vec nu_1_curr = curr_Dn(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
            arma::vec like_y_curr = dmvnorm(y_1.t(), nu_1_curr, Gamma_curr, true);
            
            like_comp_curr = like_comp_curr + arma::as_scalar(like_y_curr);
            
            // Proposal  ---------
            arma::vec vec_A_prop = A_all_state.col(pr_b_i(t_j) - 1);
            
            arma::mat Gamma_prop = {{R(0,0) / (1 - vec_A_prop(0) * vec_A_prop(0)), 
                                     R(0,1) / (1 - vec_A_prop(0) * vec_A_prop(1)), 
                                     R(0,2) / (1 - vec_A_prop(0) * vec_A_prop(2)), 
                                     R(0,3) / (1 - vec_A_prop(0) * vec_A_prop(3))},
                                     {R(1,0) / (1 - vec_A_prop(1) * vec_A_prop(0)), 
                                      R(1,1) / (1 - vec_A_prop(1) * vec_A_prop(1)), 
                                      R(1,2) / (1 - vec_A_prop(1) * vec_A_prop(2)), 
                                      R(1,3) / (1 - vec_A_prop(1) * vec_A_prop(3))},
                                      {R(2,0) / (1 - vec_A_prop(2) * vec_A_prop(0)), 
                                       R(2,1) / (1 - vec_A_prop(2) * vec_A_prop(1)), 
                                       R(2,2) / (1 - vec_A_prop(2) * vec_A_prop(2)), 
                                       R(2,3) / (1 - vec_A_prop(2) * vec_A_prop(3))},
                                       {R(3,0) / (1 - vec_A_prop(3) * vec_A_prop(0)), 
                                        R(3,1) / (1 - vec_A_prop(3) * vec_A_prop(1)), 
                                        R(3,2) / (1 - vec_A_prop(3) * vec_A_prop(2)), 
                                        R(3,3) / (1 - vec_A_prop(3) * vec_A_prop(3))}};
            
            arma::vec nu_1_prop = pr_Dn(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
            arma::vec like_y_prop = dmvnorm(y_1.t(), nu_1_prop, Gamma_prop, true);
            
            like_comp_prop = like_comp_prop + arma::as_scalar(like_y_prop);
            
        } else {
            
            arma::vec y_k_1 = y_i.col(t_j - 1);
            arma::vec y_k = y_i.col(t_j);
            
            // Current ---------
            arma::vec vec_A_curr = A_all_state.col(curr_b_i(t_j) - 1);
            arma::mat A_1_curr = arma::diagmat(vec_A_curr);
            
            arma::vec nu_k_1_curr = curr_Dn(t_j - 1) * alpha_i 
                                    + Dn_omega(t_j - 1) * w_i 
                                    + x_i(t_j - 1) * vec_beta;
            arma::vec nu_k_curr = curr_Dn(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
            
            arma::vec mean_k_curr = nu_k_curr + A_1_curr * (y_k_1 - nu_k_1_curr);
            arma::vec like_y_k_curr = dmvnorm(y_k.t(), mean_k_curr, R, true);
            
            like_comp_curr = like_comp_curr + arma::as_scalar(like_y_k_curr);
            
            // Proposal ---------
            arma::vec vec_A_prop = A_all_state.col(pr_b_i(t_j) - 1);
            arma::mat A_1_prop = arma::diagmat(vec_A_prop);
            
            arma::vec nu_k_1_prop = pr_Dn(t_j - 1) * alpha_i 
                                    + Dn_omega(t_j - 1) * w_i 
                                    + x_i(t_j - 1) * vec_beta;
            arma::vec nu_k_prop = pr_Dn(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
            
            arma::vec mean_k_prop = nu_k_prop + A_1_prop * (y_k_1 - nu_k_1_prop);
            arma::vec like_y_k_prop = dmvnorm(y_k.t(), mean_k_prop, R, true);
            
            like_comp_prop = like_comp_prop + arma::as_scalar(like_y_k_prop);
        }
    } 
    
    arma::vec mh_like = {like_comp_prop + like_comp_prop_transition, 
                         like_comp_curr + like_comp_curr_transition};
    return mh_like;
    
}

// [[Rcpp::export]]
arma::vec sub_like_MH(const int k, const int n_i, int t_pt_length, const arma::vec &par, 
                      const arma::field<arma::uvec> &par_index, const arma::mat &y_i, 
                      const arma::mat &z_i, const arma::vec &alpha_i, 
                      const arma::field<arma::mat> &x_i, 
                      const arma::field<arma::mat> &Dn_omega, const arma::vec &w_i,
                      const arma::field<arma::mat> &pr_Dn, arma::vec &pr_b_i,
                      const arma::field<arma::mat> &curr_Dn, arma::vec &curr_b_i,
                      int option_num) {
    // Parameter initialization ------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    // Transform to support (-1,1)
    arma::vec vec_A_scale = { (exp(vec_A_total(0)) - 1) / (1+exp(vec_A_total(0))),
                              (exp(vec_A_total(1)) - 1) / (1+exp(vec_A_total(1))),
                              (exp(vec_A_total(2)) - 1) / (1+exp(vec_A_total(2))),
                              (exp(vec_A_total(3)) - 1) / (1+exp(vec_A_total(3))),
                              (exp(vec_A_total(4)) - 1) / (1+exp(vec_A_total(4))),
                              (exp(vec_A_total(5)) - 1) / (1+exp(vec_A_total(5))),
                              (exp(vec_A_total(6)) - 1) / (1+exp(vec_A_total(6))),
                              (exp(vec_A_total(7)) - 1) / (1+exp(vec_A_total(7))),
                              (exp(vec_A_total(8)) - 1) / (1+exp(vec_A_total(8))),
                              (exp(vec_A_total(9)) - 1) / (1+exp(vec_A_total(9))),
                              (exp(vec_A_total(10)) - 1) / (1+exp(vec_A_total(10))),
                              (exp(vec_A_total(11)) - 1) / (1+exp(vec_A_total(11))),
                              (exp(vec_A_total(12)) - 1) / (1+exp(vec_A_total(12))),
                              (exp(vec_A_total(13)) - 1) / (1+exp(vec_A_total(13))),
                              (exp(vec_A_total(14)) - 1) / (1+exp(vec_A_total(14))),
                              (exp(vec_A_total(15)) - 1) / (1+exp(vec_A_total(15))),
                              (exp(vec_A_total(16)) - 1) / (1+exp(vec_A_total(16))),
                              (exp(vec_A_total(17)) - 1) / (1+exp(vec_A_total(17))),
                              (exp(vec_A_total(18)) - 1) / (1+exp(vec_A_total(18))),
                              (exp(vec_A_total(19)) - 1) / (1+exp(vec_A_total(19)))};
    
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    
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
        
        if(t_j == 0) {
            
            arma::vec y_1 = y_i.col(t_j);
            
            // Current ---------
            arma::vec vec_A_curr = A_all_state.col(curr_b_i(t_j) - 1);
            
            arma::mat Gamma_curr = {{R(0,0) / (1 - vec_A_curr(0) * vec_A_curr(0)), 
                                     R(0,1) / (1 - vec_A_curr(0) * vec_A_curr(1)), 
                                     R(0,2) / (1 - vec_A_curr(0) * vec_A_curr(2)), 
                                     R(0,3) / (1 - vec_A_curr(0) * vec_A_curr(3))},
                                    {R(1,0) / (1 - vec_A_curr(1) * vec_A_curr(0)), 
                                     R(1,1) / (1 - vec_A_curr(1) * vec_A_curr(1)), 
                                     R(1,2) / (1 - vec_A_curr(1) * vec_A_curr(2)), 
                                     R(1,3) / (1 - vec_A_curr(1) * vec_A_curr(3))},
                                    {R(2,0) / (1 - vec_A_curr(2) * vec_A_curr(0)), 
                                     R(2,1) / (1 - vec_A_curr(2) * vec_A_curr(1)), 
                                     R(2,2) / (1 - vec_A_curr(2) * vec_A_curr(2)), 
                                     R(2,3) / (1 - vec_A_curr(2) * vec_A_curr(3))},
                                    {R(3,0) / (1 - vec_A_curr(3) * vec_A_curr(0)), 
                                     R(3,1) / (1 - vec_A_curr(3) * vec_A_curr(1)), 
                                     R(3,2) / (1 - vec_A_curr(3) * vec_A_curr(2)), 
                                     R(3,3) / (1 - vec_A_curr(3) * vec_A_curr(3))}};
            
            arma::vec nu_1_curr = curr_Dn(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
            arma::vec like_y_curr = dmvnorm(y_1.t(), nu_1_curr, Gamma_curr, true);
            
            like_comp_curr = like_comp_curr + arma::as_scalar(like_y_curr);
            
            // Proposal  ---------
            arma::vec vec_A_prop = A_all_state.col(pr_b_i(t_j) - 1);
            
            arma::mat Gamma_prop = {{R(0,0) / (1 - vec_A_prop(0) * vec_A_prop(0)), 
                                     R(0,1) / (1 - vec_A_prop(0) * vec_A_prop(1)), 
                                     R(0,2) / (1 - vec_A_prop(0) * vec_A_prop(2)), 
                                     R(0,3) / (1 - vec_A_prop(0) * vec_A_prop(3))},
                                    {R(1,0) / (1 - vec_A_prop(1) * vec_A_prop(0)), 
                                     R(1,1) / (1 - vec_A_prop(1) * vec_A_prop(1)), 
                                     R(1,2) / (1 - vec_A_prop(1) * vec_A_prop(2)), 
                                     R(1,3) / (1 - vec_A_prop(1) * vec_A_prop(3))},
                                    {R(2,0) / (1 - vec_A_prop(2) * vec_A_prop(0)), 
                                     R(2,1) / (1 - vec_A_prop(2) * vec_A_prop(1)), 
                                     R(2,2) / (1 - vec_A_prop(2) * vec_A_prop(2)), 
                                     R(2,3) / (1 - vec_A_prop(2) * vec_A_prop(3))},
                                    {R(3,0) / (1 - vec_A_prop(3) * vec_A_prop(0)), 
                                     R(3,1) / (1 - vec_A_prop(3) * vec_A_prop(1)), 
                                     R(3,2) / (1 - vec_A_prop(3) * vec_A_prop(2)), 
                                     R(3,3) / (1 - vec_A_prop(3) * vec_A_prop(3))}};
            
            arma::vec nu_1_prop = pr_Dn(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
            arma::vec like_y_prop = dmvnorm(y_1.t(), nu_1_prop, Gamma_prop, true);
            
            like_comp_prop = like_comp_prop + arma::as_scalar(like_y_prop);
            
        } else {
            
            arma::vec y_k_1 = y_i.col(t_j - 1);
            arma::vec y_k = y_i.col(t_j);
            
            // Current ---------
            arma::vec vec_A_curr = A_all_state.col(curr_b_i(t_j) - 1);
            arma::mat A_1_curr = arma::diagmat(vec_A_curr);
            
            arma::vec nu_k_1_curr = curr_Dn(t_j - 1) * alpha_i 
                                    + Dn_omega(t_j - 1) * w_i 
                                    + x_i(t_j - 1) * vec_beta;
            arma::vec nu_k_curr = curr_Dn(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
            
            arma::vec mean_k_curr = nu_k_curr + A_1_curr * (y_k_1 - nu_k_1_curr);
            arma::vec like_y_k_curr = dmvnorm(y_k.t(), mean_k_curr, R, true);
            
            like_comp_curr = like_comp_curr + arma::as_scalar(like_y_k_curr);
            
            // Proposal ---------
            arma::vec vec_A_prop = A_all_state.col(pr_b_i(t_j) - 1);
            arma::mat A_1_prop = arma::diagmat(vec_A_prop);
            
            arma::vec nu_k_1_prop = pr_Dn(t_j - 1) * alpha_i 
                                    + Dn_omega(t_j - 1) * w_i 
                                    + x_i(t_j - 1) * vec_beta;
            arma::vec nu_k_prop = pr_Dn(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
            
            arma::vec mean_k_prop = nu_k_prop + A_1_prop * (y_k_1 - nu_k_1_prop);
            arma::vec like_y_k_prop = dmvnorm(y_k.t(), mean_k_prop, R, true);
            
            like_comp_prop = like_comp_prop + arma::as_scalar(like_y_k_prop);
        }
    } 
    
    arma::vec mh_like = {like_comp_prop, like_comp_curr};
    return mh_like;
    
}


// [[Rcpp::export]]
arma::mat state_prob_MH(const int k, const int n_i, int t_pt_length, 
                        const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                        const arma::mat &y_i, const arma::mat &z_i, const arma::vec &b_i,
                        const arma::vec &alpha_i, const arma::field<arma::mat> &x_i,
                        const arma::field<arma::mat> &Dn_omega, const arma::vec &w_i,
                        arma::mat omega_set) {
    
    // Parameter initialization ------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    // Transform to support (-1,1)
    arma::vec vec_A_scale = { (exp(vec_A_total(0)) - 1) / (1+exp(vec_A_total(0))),
                              (exp(vec_A_total(1)) - 1) / (1+exp(vec_A_total(1))),
                              (exp(vec_A_total(2)) - 1) / (1+exp(vec_A_total(2))),
                              (exp(vec_A_total(3)) - 1) / (1+exp(vec_A_total(3))),
                              (exp(vec_A_total(4)) - 1) / (1+exp(vec_A_total(4))),
                              (exp(vec_A_total(5)) - 1) / (1+exp(vec_A_total(5))),
                              (exp(vec_A_total(6)) - 1) / (1+exp(vec_A_total(6))),
                              (exp(vec_A_total(7)) - 1) / (1+exp(vec_A_total(7))),
                              (exp(vec_A_total(8)) - 1) / (1+exp(vec_A_total(8))),
                              (exp(vec_A_total(9)) - 1) / (1+exp(vec_A_total(9))),
                              (exp(vec_A_total(10)) - 1) / (1+exp(vec_A_total(10))),
                              (exp(vec_A_total(11)) - 1) / (1+exp(vec_A_total(11))),
                              (exp(vec_A_total(12)) - 1) / (1+exp(vec_A_total(12))),
                              (exp(vec_A_total(13)) - 1) / (1+exp(vec_A_total(13))),
                              (exp(vec_A_total(14)) - 1) / (1+exp(vec_A_total(14))),
                              (exp(vec_A_total(15)) - 1) / (1+exp(vec_A_total(15))),
                              (exp(vec_A_total(16)) - 1) / (1+exp(vec_A_total(16))),
                              (exp(vec_A_total(17)) - 1) / (1+exp(vec_A_total(17))),
                              (exp(vec_A_total(18)) - 1) / (1+exp(vec_A_total(18))),
                              (exp(vec_A_total(19)) - 1) / (1+exp(vec_A_total(19)))};
    
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    
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
        
        // Create new design matrix -------------------------------------------
        arma::field<arma::mat> Dn_ss_j(ss_j.n_elem);
        
        arma::vec twos(ss_j.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(ss_j == 2)) += 1;
        threes.elem(arma::find(ss_j == 3)) += 1; 
        fours.elem(arma::find(ss_j == 4)) += 1;
        fives.elem(arma::find(ss_j == 5)) += 1;
        
        arma::vec ones(ss_j.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes));
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            Dn_ss_j(jj) = arma::kron(I, bigB.row(jj));
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
                // State space component
                double q1_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(0));
                double q1 = exp(q1_sub);
                double q2_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(1));
                double q2 = exp(q2_sub);
                double q3_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(2));
                double q3 = exp(q3_sub);
                double q4_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(3));
                double q4 = exp(q4_sub);
                
                double q5_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(4));
                double q5 = exp(q5_sub);
                double q6_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(5));
                double q6 = exp(q6_sub);
                double q7_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(6));
                double q7 = exp(q7_sub);
                double q8_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(7));
                double q8 = exp(q8_sub);
                
                double q9_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(8));
                double q9 = exp(q9_sub);
                double q10_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(9));
                double q10 = exp(q10_sub);
                double q11_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(10));
                double q11 = exp(q11_sub);
                double q12_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(11));
                double q12 = exp(q12_sub);
                
                arma::mat Q = { {   1,   q1,  0,  q2,  0},
                                {   0,    1, q3,  q4,  0},
                                {  q5,   q6,  1,  q7,  0},
                                {   0,   q8,  0,   1, q9},
                                { q10,  q11,  0, q12,  1}}; 
                
                arma::vec q_row_sums = arma::sum(Q, 1);
                arma::mat P_i = Q.each_col() / q_row_sums;
                int b_k_1 = ss_j(t_j-1);
                int b_k = ss_j(t_j);
                
                like_comp_prob = like_comp_prob + log(P_i(b_k_1 - 1, b_k - 1));
            } 
            
            // (2) Response outcome --------------------------------------------
            if(t_j == 0) {
                arma::vec vec_A = A_all_state.col(ss_j(t_j) - 1);
                
                arma::mat Gamma = { {R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                     R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                     R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                     R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                    {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                     R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                     R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                     R(1,3) / (1 - vec_A(1) * vec_A(3))},
                                    {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                     R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                     R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                     R(2,3) / (1 - vec_A(2) * vec_A(3))},
                                    {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                     R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                     R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                     R(3,3) / (1 - vec_A(3) * vec_A(3))}};
                
                arma::vec y_1 = y_i.col(t_j);
                arma::vec nu_1 = Dn_ss_j(t_j) * alpha_i 
                + Dn_omega(t_j) * w_i 
                + x_i(t_j) * vec_beta;
                arma::vec like_y = dmvnorm(y_1.t(), nu_1, Gamma, true);
                
                like_comp_resp = like_comp_resp + arma::as_scalar(like_y);
                
            } else {
                arma::vec vec_A = A_all_state.col(ss_j(t_j) - 1);
                arma::mat A_1 = arma::diagmat(vec_A);
                
                arma::vec y_k_1 = y_i.col(t_j - 1);
                arma::vec y_k = y_i.col(t_j);
                arma::vec nu_k_1 = Dn_ss_j(t_j - 1) * alpha_i 
                + Dn_omega(t_j - 1) * w_i 
                + x_i(t_j - 1) * vec_beta;
                arma::vec nu_k = Dn_ss_j(t_j) * alpha_i 
                + Dn_omega(t_j) * w_i 
                + x_i(t_j) * vec_beta;
                
                arma::vec mean_k = nu_k + A_1 * (y_k_1 - nu_k_1);
                arma::vec like_y_k = dmvnorm(y_k.t(), mean_k, R, true);
                
                like_comp_resp = like_comp_resp + arma::as_scalar(like_y_k);
            }
        } 
        
        // Rcpp::Rcout << "Total: " << like_comp_prob + like_comp_resp << std::endl;
        // Rcpp::Rcout << omega_set.row(j) << std::endl;
        
        prob_dist(j,0) = like_comp_prob + like_comp_resp;
        // prob_dist(j,1) = like_comp_MH;
    }
    
    // For numerical stability, we will scale on the log-scale
    double prob_log_max = prob_dist.col(0).max();
    // Rcpp::Rcout << "Max: " << prob_log_max << std::endl;
    // double prob_log_avg = arma::mean(prob_dist.col(0));
    // Rcpp::Rcout << "Avg: " << prob_log_avg << std::endl;
    
    prob_dist.col(0) = prob_dist.col(0) - prob_log_max;
    // Rcpp::Rcout << "Centered: " << std::endl;
    // Rcpp::Rcout << prob_dist.col(0) << std::endl;
    
    prob_dist.col(0) = exp(prob_dist.col(0));
    
    // Rcpp::Rcout << "Expit" << std::endl;
    // Rcpp::Rcout << prob_dist.col(0) << std::endl;
    
    prob_dist.col(0)= (1/arma::accu(prob_dist.col(0))) * prob_dist.col(0);
    // Rcpp::Rcout << "Prob dist" << std::endl;
    // Rcpp::Rcout << prob_dist << std::endl;
    
    arma::uvec num_non_zero = arma::find(prob_dist.col(2) == 1);
    if(num_non_zero.n_elem > 1) {
        Rcpp::Rcout << "bad indicator" << std::endl;
    }
    
    return prob_dist;
}

arma::mat state_prob_only(const int k, const int n_i, int t_pt_length, 
                          const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                          const arma::mat &y_i, const arma::mat &z_i, const arma::vec &b_i,
                          arma::mat omega_set) {
    
    // Parameter initialization -----------------------------------------------
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    
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
                // State space component
                double q1_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(0));
                double q1 = exp(q1_sub);
                double q2_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(1));
                double q2 = exp(q2_sub);
                double q3_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(2));
                double q3 = exp(q3_sub);
                double q4_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(3));
                double q4 = exp(q4_sub);
                
                double q5_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(4));
                double q5 = exp(q5_sub);
                double q6_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(5));
                double q6 = exp(q6_sub);
                double q7_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(6));
                double q7 = exp(q7_sub);
                double q8_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(7));
                double q8 = exp(q8_sub);
                
                double q9_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(8));
                double q9 = exp(q9_sub);
                double q10_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(9));
                double q10 = exp(q10_sub);
                double q11_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(10));
                double q11 = exp(q11_sub);
                double q12_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(11));
                double q12 = exp(q12_sub);
                
                arma::mat Q = { {   1,   q1,  0,  q2,  0},
                                {   0,    1, q3,  q4,  0},
                                {  q5,   q6,  1,  q7,  0},
                                {   0,   q8,  0,   1, q9},
                                { q10,  q11,  0, q12,  1}}; 
                
                arma::vec q_row_sums = arma::sum(Q, 1);
                arma::mat P_i = Q.each_col() / q_row_sums;
                int b_k_1 = ss_j(t_j-1);
                int b_k = ss_j(t_j);
                
                like_comp_prob = like_comp_prob + log(P_i(b_k_1 - 1, b_k - 1));
            }
        } 
        
        // Rcpp::Rcout << "Total: " << like_comp_prob + like_comp_resp << std::endl;
        // Rcpp::Rcout << omega_set.row(j) << std::endl;
        
        prob_dist(j,0) = like_comp_prob;
    }
    
    // For numerical stability, we will scale on the log-scale
    double prob_log_max = prob_dist.col(0).max();
    // Rcpp::Rcout << "Max: " << prob_log_max << std::endl;
    // double prob_log_avg = arma::mean(prob_dist.col(0));
    // Rcpp::Rcout << "Avg: " << prob_log_avg << std::endl;
    
    prob_dist.col(0) = prob_dist.col(0) - prob_log_max;
    // Rcpp::Rcout << "Centered: " << std::endl;
    // Rcpp::Rcout << prob_dist.col(0) << std::endl;
    
    prob_dist.col(0) = exp(prob_dist.col(0));
    
    // Rcpp::Rcout << "Expit" << std::endl;
    // Rcpp::Rcout << prob_dist.col(0) << std::endl;
    
    prob_dist.col(0)= (1/arma::accu(prob_dist.col(0))) * prob_dist.col(0);
    // Rcpp::Rcout << "Prob dist" << std::endl;
    // Rcpp::Rcout << prob_dist << std::endl;
    
    arma::uvec num_non_zero = arma::find(prob_dist.col(2) == 1);
    if(num_non_zero.n_elem > 1) {
        Rcpp::Rcout << "bad indicator" << std::endl;
    }
    
    return prob_dist;
}

// [[Rcpp::export]]
Rcpp::List update_b_i_MH(const arma::vec EIDs, const arma::vec &par, 
                         const arma::field<arma::uvec> &par_index, 
                         const arma::field <arma::vec> &A, 
                         arma::field <arma::vec> &B, 
                         const arma::mat &Y, const arma::mat &z, 
                         arma::field<arma::field<arma::mat>> &Dn, 
                         const arma::field <arma::field<arma::mat>> &Xn, 
                         const arma::field<arma::field<arma::mat>> &Dn_omega, 
                         const arma::field <arma::vec> &W,
                         const arma::vec &bleed_indicator, int n_cores,
                         int t_pt_length, int option_num) {
    
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::vec eids = Y.col(0); 
    // arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6); 
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        // int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        // Sub-setting fields --------------------------------------------------
        arma::vec B_temp = B(ii);
        arma::vec A_temp = A(ii);
        arma::vec W_temp = W(ii);
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
        
        arma::field<arma::mat> Dn_temp = Dn(ii);
        arma::field<arma::mat> Dn_omega_temp = Dn_omega(ii);
        arma::field<arma::mat> Xn_temp = Xn(ii);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat z_temp = z.rows(sub_ind);
        
        arma::mat y_i = Y_temp.cols(1, 4);
        y_i = y_i.t();
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - (t_pt_length - 1); k++) {
            
            arma::vec pr_B = B_temp;
            arma::field<arma::mat> pr_Dn(Dn_temp.n_elem);
            
            // All possible state transitions given current time point ---------
            arma::mat Omega_set;
            if (clinic_rule >= 0) {
                Omega_set = Omega_fun_cpp_new_multi(k + 1, n_i, B_temp, false, 
                                                    t_pt_length);
            } else { 
                Omega_set = Omega_fun_cpp_new_multi(k + 1, n_i, B_temp, true, 
                                                    t_pt_length);
            } 
            
            // Learn the proposal distribution ---------------------------------
            int row_ind;
            if(option_num == 1) {
                row_ind = arma::randi(arma::distr_param(1, Omega_set.n_rows));
            } else {
                arma::mat ss_prob;
                if(option_num == 2) {
                    ss_prob = state_prob_only(k+1, n_i, t_pt_length, par,
                                              par_index, y_i, z_temp, B_temp,
                                              Omega_set);
                } else if(option_num == 3) {
                    ss_prob = state_prob_MH(k+1, n_i, t_pt_length, par,
                                            par_index, y_i, z_temp, B_temp,
                                            A_temp, Xn_temp, Dn_omega_temp,
                                            W_temp, Omega_set);
                }
                
                
                arma::vec x_sample = arma::linspace(1, Omega_set.n_rows, Omega_set.n_rows);
                arma::vec row_ind_vec = RcppArmadillo::sample(x_sample, 1, false, ss_prob.col(0));
                row_ind = row_ind_vec(0);
            }
            
            // Format ----------------------------------------------------------
            pr_B.rows(k, k+t_pt_length-1) = Omega_set.row(row_ind-1).t();
            
            arma::vec twos(pr_B.n_elem, arma::fill::zeros);
            arma::vec threes = twos; // THREE STATE
            arma::vec fours = twos;
            arma::vec fives = twos;
            
            twos.elem(arma::find(pr_B == 2)) += 1;
            threes.elem(arma::find(pr_B == 3)) += 1; // THREE STATE
            fours.elem(arma::find(pr_B == 4)) += 1;
            fives.elem(arma::find(pr_B == 5)) += 1;
            
            arma::vec ones(pr_B.n_elem, arma::fill::ones);
            
            arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
            bigB = arma::join_rows(bigB, arma::cumsum(threes)); // THREE STATE
            bigB = arma::join_rows(bigB, arma::cumsum(fours));
            bigB = arma::join_rows(bigB, arma::cumsum(fives));
            
            arma::mat I = arma::eye(4,4);
            for(int jj = 0; jj < n_i; jj++) {
                pr_Dn(jj) = arma::kron(I, bigB.row(jj));
            } 
            
            // log MH-ratio ----------------------------------------------------
            arma::vec mh_like_val;
            if(option_num == 1) {
                mh_like_val = full_like_MH(k+1, n_i, t_pt_length, par, 
                                           par_index, y_i, z_temp, A_temp, 
                                           Xn_temp, Dn_omega_temp, W_temp,
                                           pr_Dn, pr_B, Dn_temp, B_temp);
            } else {
                mh_like_val = sub_like_MH(k+1, n_i, t_pt_length, par, 
                                          par_index, y_i, z_temp, A_temp, 
                                          Xn_temp, Dn_omega_temp, W_temp,
                                          pr_Dn, pr_B, Dn_temp, B_temp, option_num);   
            }
            
            double diff_check = mh_like_val(0) - mh_like_val(1);
            double min_log = log(arma::randu(arma::distr_param(0,1)));
            if(diff_check > min_log){
                B_temp = pr_B;
                Dn_temp = pr_Dn;
            } 
            
        }
        B_return(ii) = B_temp;
        Dn_return(ii) = Dn_temp;
    } 
    List B_Dn = List::create(B_return, Dn_return);
    
    return B_Dn;
} 


// [[Rcpp::export]]
arma::vec state_prob_gibbs(const int k, const int n_i, int t_pt_length, 
                          const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                          const arma::mat &y_i, const arma::mat &z_i, const arma::vec &b_i,
                          const arma::vec &alpha_i, const arma::field<arma::mat> &x_i,
                          const arma::field<arma::mat> &Dn_omega, const arma::vec &w_i,
                          arma::mat omega_set) {
    
    // Parameter initialization ------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    // Transform to support (-1,1)
    arma::vec vec_A_scale = { (exp(vec_A_total(0)) - 1) / (1+exp(vec_A_total(0))),
                              (exp(vec_A_total(1)) - 1) / (1+exp(vec_A_total(1))),
                              (exp(vec_A_total(2)) - 1) / (1+exp(vec_A_total(2))),
                              (exp(vec_A_total(3)) - 1) / (1+exp(vec_A_total(3))),
                              (exp(vec_A_total(4)) - 1) / (1+exp(vec_A_total(4))),
                              (exp(vec_A_total(5)) - 1) / (1+exp(vec_A_total(5))),
                              (exp(vec_A_total(6)) - 1) / (1+exp(vec_A_total(6))),
                              (exp(vec_A_total(7)) - 1) / (1+exp(vec_A_total(7))),
                              (exp(vec_A_total(8)) - 1) / (1+exp(vec_A_total(8))),
                              (exp(vec_A_total(9)) - 1) / (1+exp(vec_A_total(9))),
                              (exp(vec_A_total(10)) - 1) / (1+exp(vec_A_total(10))),
                              (exp(vec_A_total(11)) - 1) / (1+exp(vec_A_total(11))),
                              (exp(vec_A_total(12)) - 1) / (1+exp(vec_A_total(12))),
                              (exp(vec_A_total(13)) - 1) / (1+exp(vec_A_total(13))),
                              (exp(vec_A_total(14)) - 1) / (1+exp(vec_A_total(14))),
                              (exp(vec_A_total(15)) - 1) / (1+exp(vec_A_total(15))),
                              (exp(vec_A_total(16)) - 1) / (1+exp(vec_A_total(16))),
                              (exp(vec_A_total(17)) - 1) / (1+exp(vec_A_total(17))),
                              (exp(vec_A_total(18)) - 1) / (1+exp(vec_A_total(18))),
                              (exp(vec_A_total(19)) - 1) / (1+exp(vec_A_total(19)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
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
        
        // Create new design matrix -------------------------------------------
        arma::field<arma::mat> Dn_ss_j(ss_j.n_elem);
        
        arma::vec twos(ss_j.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(ss_j == 2)) += 1;
        threes.elem(arma::find(ss_j == 3)) += 1; 
        fours.elem(arma::find(ss_j == 4)) += 1;
        fives.elem(arma::find(ss_j == 5)) += 1;
        
        arma::vec ones(ss_j.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes));
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            Dn_ss_j(jj) = arma::kron(I, bigB.row(jj));
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
                    // State space component
                    double q1_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(0));
                    double q1 = exp(q1_sub);
                    double q2_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(1));
                    double q2 = exp(q2_sub);
                    double q3_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(2));
                    double q3 = exp(q3_sub);
                    double q4_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(3));
                    double q4 = exp(q4_sub);
                    
                    double q5_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(4));
                    double q5 = exp(q5_sub);
                    double q6_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(5));
                    double q6 = exp(q6_sub);
                    double q7_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(6));
                    double q7 = exp(q7_sub);
                    double q8_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(7));
                    double q8 = exp(q8_sub);
                    
                    double q9_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(8));
                    double q9 = exp(q9_sub);
                    double q10_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(9));
                    double q10 = exp(q10_sub);
                    double q11_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(10));
                    double q11 = exp(q11_sub);
                    double q12_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(11));
                    double q12 = exp(q12_sub);
                    
                    arma::mat Q = { {   1,   q1,  0,  q2,  0},
                                    {   0,    1, q3,  q4,  0},
                                    {  q5,   q6,  1,  q7,  0},
                                    {   0,   q8,  0,   1, q9},
                                    { q10,  q11,  0, q12,  1}}; 
                    
                    arma::vec q_row_sums = arma::sum(Q, 1);
                    arma::mat P_i = Q.each_col() / q_row_sums;
                    int b_k_1 = ss_j(t_j-1);
                    int b_k = ss_j(t_j);
                    
                    like_comp_prob = like_comp_prob + log(P_i(b_k_1 - 1, b_k - 1));
                } 
            }
            
            // (2) Response outcome --------------------------------------------
            if(t_j == 0) {
                arma::vec vec_A = A_all_state.col(ss_j(t_j) - 1);
                
                arma::mat Gamma = { {R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                     R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                     R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                     R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                     {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                      R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                      R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                      R(1,3) / (1 - vec_A(1) * vec_A(3))},
                                      {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                       R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                       R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                       R(2,3) / (1 - vec_A(2) * vec_A(3))},
                                       {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                        R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                        R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                        R(3,3) / (1 - vec_A(3) * vec_A(3))}};
                
                arma::vec y_1 = y_i.col(t_j);
                arma::vec nu_1 = Dn_ss_j(t_j) * alpha_i 
                                + Dn_omega(t_j) * w_i 
                                + x_i(t_j) * vec_beta;
                arma::vec like_y = dmvnorm(y_1.t(), nu_1, Gamma, true);
                
                like_comp_resp = like_comp_resp + arma::as_scalar(like_y);
                
            } else { 
                arma::vec vec_A = A_all_state.col(ss_j(t_j) - 1);
                arma::mat A_1 = arma::diagmat(vec_A);
                
                arma::vec y_k_1 = y_i.col(t_j - 1);
                arma::vec y_k = y_i.col(t_j);
                arma::vec nu_k_1 = Dn_ss_j(t_j - 1) * alpha_i 
                                    + Dn_omega(t_j - 1) * w_i 
                                    + x_i(t_j - 1) * vec_beta;
                arma::vec nu_k = Dn_ss_j(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;
                
                arma::vec mean_k = nu_k + A_1 * (y_k_1 - nu_k_1);
                arma::vec like_y_k = dmvnorm(y_k.t(), mean_k, R, true);
                
                like_comp_resp = like_comp_resp + arma::as_scalar(like_y_k);
            }
            // Rcpp::Rcout << "prob: " << like_comp_prob << std::endl;
            // Rcpp::Rcout << "resp: " << like_comp_resp << std::endl;
        } 
        
        // Rcpp::Rcout << "Total: " << like_comp_prob + like_comp_resp << std::endl;
        // Rcpp::Rcout << omega_set.row(j) << std::endl;
        
        prob_dist(j) = like_comp_prob + like_comp_resp;
        // Rcpp::Rcout << ss_j.t() << std::endl;
        // Rcpp::Rcout << "prob = " << prob_dist(j) << std::endl;
    }
    
    // For numerical stability, we will scale on the log-scale
    double prob_log_max = prob_dist.max();
    // Rcpp::Rcout << "Max: " << prob_log_max << std::endl;
    // double prob_log_avg = arma::mean(prob_dist);
    // Rcpp::Rcout << "Avg: " << prob_log_avg << std::endl;
    
    prob_dist = prob_dist - prob_log_max;
    // Rcpp::Rcout << "Centered: " << std::endl;
    // Rcpp::Rcout << prob_dist << std::endl;
    
    prob_dist = exp(prob_dist);
    
    // Rcpp::Rcout << "Expit" << std::endl;
    // Rcpp::Rcout << prob_dist << std::endl;
    
    prob_dist= (1/arma::accu(prob_dist)) * prob_dist;
    // Rcpp::Rcout << "Prob dist" << std::endl;
    // Rcpp::Rcout << prob_dist << std::endl;
    
    return prob_dist;
}

// [[Rcpp::export]]
Rcpp::List update_b_i_gibbs( const arma::vec EIDs, const arma::vec &par, 
                             const arma::field<arma::uvec> &par_index, 
                             const arma::field <arma::vec> &A, 
                             arma::field <arma::vec> &B, 
                             const arma::mat &Y, const arma::mat &z, 
                             arma::field<arma::field<arma::mat>> &Dn, 
                             const arma::field <arma::field<arma::mat>> &Xn, 
                             const arma::field<arma::field<arma::mat>> &Dn_omega, 
                             const arma::field <arma::vec> &W,
                             const arma::vec &bleed_indicator, int n_cores,
                             int t_pt_length) {
    
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::vec eids = Y.col(0); 
    // arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6); 
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        // int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        // Sub-setting fields --------------------------------------------------
        arma::vec B_temp = B(ii);
        arma::vec A_temp = A(ii);
        arma::vec W_temp = W(ii);
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
        
        arma::field<arma::mat> Dn_temp = Dn(ii);
        arma::field<arma::mat> Dn_omega_temp = Dn_omega(ii);
        arma::field<arma::mat> Xn_temp = Xn(ii);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat z_temp = z.rows(sub_ind);
        
        arma::mat y_i = Y_temp.cols(1, 4);
        y_i = y_i.t();
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - (t_pt_length - 1); k++) {
        // for (int k = n_i - t_pt_length; k >= 0; k--) {
            
            // All possible state transitions given current time point ---------
            arma::mat Omega_set;
            if (clinic_rule >= 0) {
                Omega_set = Omega_fun_cpp_new_multi(k + 1, n_i, B_temp, false, 
                                                    t_pt_length);
            } else {
                Omega_set = Omega_fun_cpp_new_multi(k + 1, n_i, B_temp, true, 
                                                    t_pt_length);
            } 
            
            // Learn the proposal distribution ---------------------------------
            arma::vec ss_prob = state_prob_gibbs(k+1, n_i, t_pt_length, par,
                                                 par_index, y_i, z_temp, B_temp,
                                                 A_temp, Xn_temp, Dn_omega_temp,
                                                 W_temp, Omega_set);
            
            arma::vec x_sample = arma::linspace(1, Omega_set.n_rows, Omega_set.n_rows);
            arma::vec row_ind = RcppArmadillo::sample(x_sample, 1, false, ss_prob);
            
            // Gibbs update ----------------------------------------------------
            B_temp.rows(k, k+t_pt_length-1) = Omega_set.row(row_ind(0)-1).t();
            
            arma::vec twos(B_temp.n_elem, arma::fill::zeros);
            arma::vec threes = twos; // THREE STATE
            arma::vec fours = twos;
            arma::vec fives = twos;
            
            twos.elem(arma::find(B_temp == 2)) += 1;
            threes.elem(arma::find(B_temp == 3)) += 1; // THREE STATE
            fours.elem(arma::find(B_temp == 4)) += 1;
            fives.elem(arma::find(B_temp == 5)) += 1;
            
            arma::vec ones(B_temp.n_elem, arma::fill::ones);
            
            arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
            bigB = arma::join_rows(bigB, arma::cumsum(threes)); // THREE STATE
            bigB = arma::join_rows(bigB, arma::cumsum(fours));
            bigB = arma::join_rows(bigB, arma::cumsum(fives));
            
            arma::mat I = arma::eye(4,4);
            for(int jj = 0; jj < n_i; jj++) {
                Dn_temp(jj) = arma::kron(I, bigB.row(jj));
            } 
        }
        B_return(ii) = B_temp;
        Dn_return(ii) = Dn_temp;
    } 
    List B_Dn = List::create(B_return, Dn_return);
    
    return B_Dn;
} 

// *** Using ***
// [[Rcpp::export]]
arma::field <arma::vec> update_alpha_i_cpp( const arma::vec &EIDs, const arma::vec &par, 
                                            const arma::field<arma::uvec> &par_index,
                                            const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn, 
                                            const arma::field <arma::field<arma::mat>> &Xn,
                                            const arma::field<arma::field<arma::mat>> &Dn_omega, 
                                            const arma::field <arma::vec> &W, 
                                            arma::field <arma::vec> &B, int n_cores){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_alpha_tilde = par.elem(par_index(1) - 1);
    
    arma::vec sigma_upsilon_vec = par.elem(par_index(2) - 1);
    arma::mat Upsilon = arma::reshape(sigma_upsilon_vec, vec_alpha_tilde.n_elem,
                                      vec_alpha_tilde.n_elem);
    arma::mat inv_Upsilon = arma::inv_sympd(Upsilon);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    // -------------------------------------------------------------------------
    
    arma::vec eids = Y.col(0);
    
    arma::field<arma::vec> A(EIDs.n_elem);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::field<arma::mat> Dn_alpha_full = Dn(ii);
        arma::field<arma::mat> Dn_omega_full = Dn_omega(ii);
        arma::field<arma::mat> Xn_full = Xn(ii);
        arma::vec vec_omega_ii = W(ii);
        arma::vec b_i = B(ii);
        
        arma::mat interm_W(vec_alpha_tilde.n_elem, vec_alpha_tilde.n_elem, 
                           arma::fill::zeros); 
        arma::vec interm_V(vec_alpha_tilde.n_elem, arma::fill::zeros); 
        
        // Skip k = 0 because we deal with it later
        for(int k = 1; k < n_i; k++) {
            
            arma::mat D_k = Dn_alpha_full(k);
            arma::mat D_k_o = Dn_omega_full(k);
            arma::mat X_k = Xn_full(k);
            
            arma::mat D_k_1 = Dn_alpha_full(k-1);
            arma::mat D_k_o_1 = Dn_omega_full(k-1);
            arma::mat X_k_1 = Xn_full(k-1);
            
            arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            arma::mat A_1_k = arma::diagmat(vec_A_k);
            
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            arma::mat diff_hold_y = Y_i.col(k) - A_1_k * Y_i.col(k-1);
            
            arma::vec m_k = diff_hold_y - diff_hold_d_o*vec_omega_ii - diff_hold_x * vec_beta;
            
            interm_W = interm_W + diff_hold_d.t() * invR * diff_hold_d;
            
            interm_V = interm_V + diff_hold_d.t() * invR * m_k;
        }
        
        arma::vec vec_A = A_all_state.col(b_i(0) - 1);
        arma::mat Gamma = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                            R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                            R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                            R(0,3) / (1 - vec_A(0) * vec_A(3))},
                           {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                            R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                            R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                            R(1,3) / (1 - vec_A(1) * vec_A(3))},
                           {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                            R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                            R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                            R(2,3) / (1 - vec_A(2) * vec_A(3))},
                           {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                            R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                            R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                            R(3,3) / (1 - vec_A(3) * vec_A(3))}};
        arma::mat inv_Gamma = arma::inv_sympd(Gamma);
        
        arma::mat W_i_inv = inv_Upsilon + Dn_alpha_full(0).t() * inv_Gamma * Dn_alpha_full(0) + interm_W;
        
        arma::mat W_i = arma::inv_sympd(W_i_inv);
        
        arma::vec y_diff_0 = Y_i.col(0) - Dn_omega_full(0) * vec_omega_ii - Xn_full(0) * vec_beta;
        arma::mat V_i = inv_Upsilon * vec_alpha_tilde + Dn_alpha_full(0).t() * inv_Gamma * y_diff_0 + interm_V;
        
        arma::vec mu = W_i * V_i;
        
        arma::mat alpha_i = rmvnorm(1, mu, W_i);
        arma::vec vec_alpha_i = alpha_i.t();
        
        A(ii) = vec_alpha_i;
    }
    
    return A;
  
}

// *** Using ***
// [[Rcpp::export]]
arma::field <arma::vec> update_omega_i_cpp( const arma::vec &EIDs, const arma::vec &par, 
                                            const arma::field<arma::uvec> &par_index,
                                            const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn, 
                                            const arma::field <arma::field<arma::mat>> &Xn,
                                            const arma::field<arma::field<arma::mat>> &Dn_omega, 
                                            const arma::field <arma::vec> &A, 
                                            arma::field <arma::vec> &B, int n_cores){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_omega_tilde = par.elem(par_index(7) - 1);
    
    arma::vec vec_upsilon_omega = par.elem(par_index(8) - 1);
    arma::vec vec_upsilon_omega_inv = 1 / exp(vec_upsilon_omega);
    arma::mat inv_Upsilon = arma::diagmat(vec_upsilon_omega_inv);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    // -------------------------------------------------------------------------
    
    arma::vec eids = Y.col(0);
    
    arma::field<arma::vec> W(EIDs.n_elem);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::field<arma::mat> Dn_alpha_full = Dn(ii);
        arma::field<arma::mat> Dn_omega_full = Dn_omega(ii);
        arma::field<arma::mat> Xn_full = Xn(ii);
        arma::vec vec_alpha_ii = A(ii);
        arma::vec b_i = B(ii);

        arma::mat interm_W(vec_omega_tilde.n_elem, vec_omega_tilde.n_elem, arma::fill::zeros);
        arma::vec interm_V(vec_omega_tilde.n_elem, arma::fill::zeros);
        for(int k = 1; k < n_i; k++) {
            
            arma::mat D_k = Dn_alpha_full(k);
            arma::mat D_k_o = Dn_omega_full(k);
            arma::mat X_k = Xn_full(k);
            
            arma::mat D_k_1 = Dn_alpha_full(k-1);
            arma::mat D_k_o_1 = Dn_omega_full(k-1);
            arma::mat X_k_1 = Xn_full(k-1);
            
            arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            arma::mat A_1_k = arma::diagmat(vec_A_k);
            
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            arma::mat diff_hold_y = Y_i.col(k) - A_1_k * Y_i.col(k-1);
            
            arma::vec m_k = diff_hold_y - diff_hold_d*vec_alpha_ii - diff_hold_x*vec_beta;
            
            interm_W = interm_W + diff_hold_d_o.t() * invR * diff_hold_d_o;
            
            interm_V = interm_V + diff_hold_d_o.t() * invR * m_k;
        }
        
        arma::vec vec_A = A_all_state.col(b_i(0) - 1);
        arma::mat Gamma = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                            R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                            R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                            R(0,3) / (1 - vec_A(0) * vec_A(3))},
                           {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                            R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                            R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                            R(1,3) / (1 - vec_A(1) * vec_A(3))},
                           {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                            R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                            R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                            R(2,3) / (1 - vec_A(2) * vec_A(3))},
                           {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                            R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                            R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                            R(3,3) / (1 - vec_A(3) * vec_A(3))}};
        arma::mat inv_Gamma = arma::inv_sympd(Gamma);
        
        arma::mat W_i_inv = inv_Upsilon + Dn_omega_full(0).t() * inv_Gamma * Dn_omega_full(0) + interm_W;
        arma::mat W_i = arma::inv_sympd(W_i_inv);
        
        arma::vec y_diff_0 = Y_i.col(0) - Dn_alpha_full(0) * vec_alpha_ii - Xn_full(0) * vec_beta;
        arma::mat V_i = inv_Upsilon * vec_omega_tilde + Dn_omega_full(0).t() * inv_Gamma * y_diff_0 + interm_V;
        
        arma::vec mu = W_i * V_i;
        
        arma::mat omega_i = rmvnorm(1, mu, W_i);
        
        arma::vec vec_omega_i = omega_i.t();
        
        W(ii) = vec_omega_i;
    }
    
    return W;
}

// *** Using ***
// [[Rcpp::export]]
arma::vec update_alpha_tilde_cpp( const arma::vec EIDs, arma::vec par, 
                                  const arma::field<arma::uvec> par_index,
                                  const arma::field <arma::vec> A, const arma::mat Y){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    // Prior MEAN (vec_alpha_tilde)
    arma::vec vec_alpha_tilde_0 = { 9, -1,  1, 0, 0,
                                    85,  5, -5, 0, 0,
                                    75, -5,  5, 0, 0,
                                    5,  1, -1, 0, 0};
    
    // Prior PRECISION (vec_alpha_tilde)
    arma::vec inv_Sigma_alpha_diag(vec_alpha_tilde_0.n_elem, arma::fill::ones);
    inv_Sigma_alpha_diag = 0.05 * inv_Sigma_alpha_diag;
    arma::mat inv_Sigma_alpha = arma::diagmat(inv_Sigma_alpha_diag);
    
    arma::vec sigma_upsilon_vec = par.elem(par_index(2) - 1);
    arma::mat Upsilon = arma::reshape(sigma_upsilon_vec, 
                                      vec_alpha_tilde_0.n_elem, 
                                      vec_alpha_tilde_0.n_elem);
    arma::mat inv_Upsilon = arma::inv_sympd(Upsilon);
    // -------------------------------------------------------------------------
    
    int N_id = EIDs.n_elem;
    
    arma::mat inv_W = inv_Sigma_alpha + N_id * inv_Upsilon;
    arma::mat W = arma::inv_sympd(inv_W);
    
    arma::vec total_alpha = A(0);
    for (int i = 1; i < A.n_elem; i++){
      total_alpha = total_alpha + A(i);
    }
    
    arma::vec V = inv_Sigma_alpha * vec_alpha_tilde_0 + inv_Upsilon * total_alpha;
    
    arma::vec mu = W * V;
    
    arma::vec alpha_tilde_temp = arma::mvnrnd(mu, W, 1);
    par.elem(par_index(1) - 1) = alpha_tilde_temp;
    
    // // Ensuring that the parent means for HR and MAP obey the directional assumption
    // int count_while_loop = 0;
    // int count_while_loop_big = 0;
    // while(alpha_tilde_temp(1) > alpha_tilde_temp(2)) {
    //     alpha_tilde_temp = arma::mvnrnd(mu, U, 1);
    //     count_while_loop += 1;
    //     if(count_while_loop > 1000) {
    //         count_while_loop_big += 1;
    //         Rcpp::Rcout << "stuck in alpha tilde " << count_while_loop_big << std::endl;
    //         // Rcpp::Rcout << "hr bleed: " << alpha_tilde_temp(4) << ", hr recov: " << alpha_tilde_temp(5) <<
    //         //     ", MAP bleed: " << alpha_tilde_temp(7) << ", MAP recov: " << alpha_tilde_temp(8) << std::endl;
    //         count_while_loop = 0;
    //     }
    // }
    
    return par;
}

// *** Using ***
// [[Rcpp::export]]
arma::vec update_omega_tilde_cpp( const arma::vec EIDs, arma::vec par, 
                                  const arma::field<arma::uvec> par_index,
                                  const arma::field <arma::vec> W, const arma::mat Y){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    // Prior MEAN (vec_omega_tilde)
    arma::vec vec_omega_tilde_0 = {-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1,
                                    1,-1,-1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1,-1,-1,
                                    1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,-1,-1,-1,
                                    1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,
                                   -1,-1,-1,-1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,
                                   -1,-1,-1, 1,-1,-1,-1,-1,-1};
    
    // Prior PRECISION (vec_alpha_tilde)
    arma::vec inv_sig_omega_vec(vec_omega_tilde_0.n_elem, arma::fill::ones);
    inv_sig_omega_vec = 0.05 * inv_sig_omega_vec;
    arma::mat inv_Sigma_omega = arma::diagmat(inv_sig_omega_vec);
    
    arma::vec vec_upsilon_omega = par.elem(par_index(8) - 1);
    arma::vec vec_upsilon_omega_inv = 1 / exp(vec_upsilon_omega);
    arma::mat inv_Upsilon = arma::diagmat(vec_upsilon_omega_inv);
    // -------------------------------------------------------------------------
    
    int N_id = EIDs.n_elem;
    
    arma::mat inv_W = inv_Sigma_omega + N_id * inv_Upsilon;
    arma::mat W_var = arma::inv_sympd(inv_W);
    
    arma::mat total_omega = W(0);
    for (int i = 1; i < W.n_elem; i++) {
        total_omega = total_omega + W(i);
    }
    
    arma::vec V = inv_Sigma_omega * vec_omega_tilde_0 + inv_Upsilon * total_omega;
    
    arma::vec mu = W_var * V;
    
    par.elem(par_index(7) - 1) = rmvnorm(1, mu, W_var);
    
    return par;
}

// *** Using ***
// [[Rcpp::export]]
arma::vec update_beta_upsilon_cpp(const arma::vec &EIDs, arma::vec par, 
                                  const arma::field<arma::uvec> &par_index,
                                  const arma::field <arma::vec> &A, const arma::mat &Y,
                                  const arma::field<arma::field<arma::mat>> &Dn, 
                                  const arma::field <arma::field<arma::mat>> &Xn, 
                                  const arma::field<arma::field<arma::mat>> &Dn_omega, 
                                  const arma::field <arma::vec> &W, arma::field <arma::vec> &B,
                                  int n_cores) {
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    // Prior (beta)
    arma::uvec vec_beta_ind = par_index(0);
    arma::vec vec_beta_0 = {0.25, -2, 2, -0.25};
    arma::vec scalar_mult_beta(vec_beta_ind.n_elem, arma::fill::ones);
    scalar_mult_beta = 0.05 * scalar_mult_beta;
    scalar_mult_beta(0) = 4;
    arma::mat inv_Sigma_beta = arma::diagmat(scalar_mult_beta);
    
    // Prior (upsilon_alpha)
    arma::uvec vec_sigma_upsilon_ind = par_index(2);
    int nu_Upsilon = 22;
    arma::vec scalar_mult2 = {  4, 0.01, 0.01, 0.25, 0.25, 
                              100,    1,    1,   25,   25, 
                              100,    1,    1,   25,   25, 
                                1, 0.01, 0.01, 0.25, 0.25};
    scalar_mult2 = (nu_Upsilon - 20 - 1) * scalar_mult2;
    arma::mat Psi_Upsilon = arma::diagmat(scalar_mult2);
    
    // Other parameters
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    
    arma::vec vec_alpha_tilde = par.elem(par_index(1) - 1);
    // -------------------------------------------------------------------------
    
    arma::vec eids = Y.col(0);

    // CHANGE THIS!
    arma::vec sum_in_V(vec_beta_ind.n_elem, arma::fill::zeros);
    arma::mat sum_in_inv_W(vec_beta_ind.n_elem, vec_beta_ind.n_elem, 
                           arma::fill::zeros);
    arma::mat sum_in_Upsilon_cov(vec_alpha_tilde.n_elem, vec_alpha_tilde.n_elem,
                                 arma::fill::zeros);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
    
        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_omega_i = W(ii);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::field<arma::mat> Dn_alpha_full = Dn(ii);
        arma::field<arma::mat> Dn_omega_full = Dn_omega(ii);
        arma::field<arma::mat> Xn_full = Xn(ii);
        
        arma::mat interm_W(vec_beta_ind.n_elem, vec_beta_ind.n_elem, 
                           arma::fill::zeros);
        arma::vec interm_V(vec_beta_ind.n_elem, arma::fill::zeros);
        
        for(int k = 1; k < n_i; k++) {
            arma::mat D_k = Dn_alpha_full(k);
            arma::mat D_k_o = Dn_omega_full(k);
            arma::mat X_k = Xn_full(k);
            
            arma::mat D_k_1 = Dn_alpha_full(k-1);
            arma::mat D_k_o_1 = Dn_omega_full(k-1);
            arma::mat X_k_1 = Xn_full(k-1);
            
            arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            arma::mat A_1_k = arma::diagmat(vec_A_k);
            
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            
            arma::vec m_k = Y_i.col(k) - A_1_k * Y_i.col(k-1) 
                                - diff_hold_d * vec_alpha_i - diff_hold_d_o * vec_omega_i;
            
            interm_W = interm_W + diff_hold_x.t() * invR * diff_hold_x;
            interm_V = interm_V + diff_hold_x.t() * invR * m_k;
        }
        
        arma::vec vec_A = A_all_state.col(b_i(0) - 1);
        arma::mat Gamma = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                            R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                            R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                            R(0,3) / (1 - vec_A(0) * vec_A(3))},
                           {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                            R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                            R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                            R(1,3) / (1 - vec_A(1) * vec_A(3))},
                           {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                            R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                            R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                            R(2,3) / (1 - vec_A(2) * vec_A(3))},
                           {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                            R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                            R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                            R(3,3) / (1 - vec_A(3) * vec_A(3))}};
        arma::mat inv_Gamma = arma::inv_sympd(Gamma);
        
        arma::mat W_i_inv = Xn_full(0).t() * inv_Gamma * Xn_full(0) + interm_W;
        
        arma::mat V_i = Xn_full(0).t() * inv_Gamma * 
                (Y_i.col(0) - Dn_alpha_full(0)*vec_alpha_i - Dn_omega_full(0)*vec_omega_i) + interm_V;
        
        arma::mat hold2 = vec_alpha_i - vec_alpha_tilde;
        arma::mat in_Upsilon_cov = hold2 * hold2.t();
        
        sum_in_V = sum_in_V + V_i;
        sum_in_inv_W = sum_in_inv_W + W_i_inv;
        sum_in_Upsilon_cov = sum_in_Upsilon_cov + in_Upsilon_cov;
    }
    
    arma::vec V = inv_Sigma_beta * vec_beta_0 + sum_in_V;
    
    arma::mat inv_W = inv_Sigma_beta + sum_in_inv_W;
    arma::mat W_b = arma::inv_sympd(inv_W);
    
    arma::mat Upsilon_cov = Psi_Upsilon + sum_in_Upsilon_cov;
    
    int N_id = EIDs.n_elem;
    
    par.elem(vec_beta_ind - 1) = arma::mvnrnd(W_b * V, W_b);
    par.elem(vec_sigma_upsilon_ind - 1) = arma::vectorise(riwish(nu_Upsilon + N_id, Upsilon_cov));
    
    return par;
}

// [[Rcpp::export]]
arma::mat update_Y_i_cpp( const arma::vec &EIDs, const arma::vec &par, 
                          const arma::field<arma::uvec> &par_index, 
                          const arma::field <arma::vec> &A, arma::mat Y,
                          arma::field <arma::field<arma::mat>> &Dn, 
                          const arma::field <arma::field<arma::mat>> &Xn, const arma::mat &otype,
                          const arma::field<arma::field<arma::mat>> &Dn_omega,
                          const arma::field <arma::vec> &W, arma::field <arma::vec> &B,
                          int n_cores) {

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID

    arma::mat newY(Y.n_rows, 4); 
    arma::vec eids = Y.col(0);
    
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);

    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    // Transform to support (-1,1)
    arma::vec vec_A_scale = { (exp(vec_A_total(0)) - 1) / (1+exp(vec_A_total(0))),
                              (exp(vec_A_total(1)) - 1) / (1+exp(vec_A_total(1))),
                              (exp(vec_A_total(2)) - 1) / (1+exp(vec_A_total(2))),
                              (exp(vec_A_total(3)) - 1) / (1+exp(vec_A_total(3))),
                              (exp(vec_A_total(4)) - 1) / (1+exp(vec_A_total(4))),
                              (exp(vec_A_total(5)) - 1) / (1+exp(vec_A_total(5))),
                              (exp(vec_A_total(6)) - 1) / (1+exp(vec_A_total(6))),
                              (exp(vec_A_total(7)) - 1) / (1+exp(vec_A_total(7))),
                              (exp(vec_A_total(8)) - 1) / (1+exp(vec_A_total(8))),
                              (exp(vec_A_total(9)) - 1) / (1+exp(vec_A_total(9))),
                              (exp(vec_A_total(10)) - 1) / (1+exp(vec_A_total(10))),
                              (exp(vec_A_total(11)) - 1) / (1+exp(vec_A_total(11))),
                              (exp(vec_A_total(12)) - 1) / (1+exp(vec_A_total(12))),
                              (exp(vec_A_total(13)) - 1) / (1+exp(vec_A_total(13))),
                              (exp(vec_A_total(14)) - 1) / (1+exp(vec_A_total(14))),
                              (exp(vec_A_total(15)) - 1) / (1+exp(vec_A_total(15))),
                              (exp(vec_A_total(16)) - 1) / (1+exp(vec_A_total(16))),
                              (exp(vec_A_total(17)) - 1) / (1+exp(vec_A_total(17))),
                              (exp(vec_A_total(18)) - 1) / (1+exp(vec_A_total(18))),
                              (exp(vec_A_total(19)) - 1) / (1+exp(vec_A_total(19)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE
    
      // omp_set_num_threads(n_cores);
      // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {	
        
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);

        arma::vec b_i = B(ii);

        arma::field<arma::mat> Dn_ii = Dn(ii);
        arma::field<arma::mat> Dn_omega_ii = Dn_omega(ii);
        arma::field<arma::mat> Xn_ii = Xn(ii);
        arma::vec vec_alpha_ii = A(ii);
        arma::vec vec_omega_ii = W(ii);

        // Index of observed versus missing data
        // 1 = observed, 0 = missing
        arma::mat otype_i = otype.rows(sub_ind);
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        otype_i = otype_i.t();
        arma::mat Y_i_new = Y_i;

        for(int k = 0; k < Y_i.n_cols; k++) {
            if(all(otype_i.col(k) == 1)) {
                Y_i_new.col(k) = Y_i.col(k);
            } else {
                if(k == 0) {
                    arma::vec vec_A = A_all_state.col(b_i(k) - 1);
                    arma::mat Gamma = { {R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                         R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                         R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                         R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                         {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                          R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                          R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                          R(1,3) / (1 - vec_A(1) * vec_A(3))},
                                          {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                           R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                           R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                           R(2,3) / (1 - vec_A(2) * vec_A(3))},
                                           {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                            R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                            R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                            R(3,3) / (1 - vec_A(3) * vec_A(3))}};
                    arma::mat inv_Gamma = arma::inv_sympd(Gamma);
                    
                    arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                    arma::mat A_p1 = arma::diagmat(vec_A_p1);
                    
                    arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                    arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + Dn_omega_ii(k+1) * vec_omega_ii + Xn_ii(k+1) * vec_beta;
                    
                    arma::vec y_val_kp1 = Y_i_new.col(k+1);

                    arma::mat inv_W_i = inv_Gamma + A_p1.t() * invR * A_p1;
                    arma::mat W_i = inv(inv_W_i);
                    arma::vec V_i = inv_Gamma*nu_k + A_p1.t()*invR*(y_val_kp1 - nu_k_p1 + A_p1*nu_k);

                    arma::vec y_i_mean = W_i * V_i;

                    arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    // int count_while_loop = 0;
                    // int count_while_loop_big = 0;
                    // while(arma::any(update_value <= 0)) {
                    //         new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                    //         update_value = Y_i_new.col(k);
                    //         update_value.elem(ind_replace) = new_value.elem(ind_replace);

                    //         count_while_loop += 1;
                    //         if(count_while_loop > 10000) {
                    //             count_while_loop_big += 1;
                    //             Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                    //             count_while_loop = 0;
                    //         }
                    //         if(count_while_loop_big > 1000) {
                    //             break;
                    //         }
                    // }

                    Y_i_new.col(k) = update_value;
                } else if(k == Y_i.n_cols - 1) {
                    
                    arma::vec vec_A = A_all_state.col(b_i(k) - 1);
                    arma::mat A_k = arma::diagmat(vec_A);
                    
                    arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                    arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + Dn_omega_ii(k-1) * vec_omega_ii +  Xn_ii(k-1) * vec_beta;
                    
                    arma::vec y_val_km1 = Y_i_new.col(k-1);

                    arma::vec y_i_mean = nu_k + A_k * (y_val_km1 - nu_k_m1);

                    arma::vec new_value = arma::mvnrnd(y_i_mean, R, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    // int count_while_loop = 0;
                    // int count_while_loop_big = 0;
                    // while(arma::any(update_value <= 0)) {
                    //     new_value = arma::mvnrnd(y_i_mean, R, 1);
                    //     update_value = Y_i_new.col(k);
                    //     update_value.elem(ind_replace) = new_value.elem(ind_replace);

                    //     count_while_loop += 1;
                    //     if(count_while_loop > 10000) {
                    //         count_while_loop_big += 1;
                    //         Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                    //         count_while_loop = 0;
                    //     }
                    //     if(count_while_loop_big > 1000) {
                    //         break;
                    //     }
                    // }

                    Y_i_new.col(k) = update_value;
                } else {

                    arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
                    arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                    arma::mat A_k = arma::diagmat(vec_A_k);                    
                    arma::mat A_p1 = arma::diagmat(vec_A_p1);

                    arma::vec nu_k    = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                    arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + Dn_omega_ii(k-1) * vec_omega_ii + Xn_ii(k-1) * vec_beta;
                    arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + Dn_omega_ii(k+1) * vec_omega_ii + Xn_ii(k+1) * vec_beta;
                    
                    arma::vec y_val_km1 = Y_i_new.col(k-1);
                    arma::vec y_val_kp1 = Y_i_new.col(k+1);

                    arma::mat inv_W_i = invR + A_p1.t() * invR * A_p1;
                    arma::mat W_i = inv(inv_W_i);

                    arma::vec V_i = invR * (nu_k + A_k * (y_val_km1 - nu_k_m1)) + 
                                        A_p1.t() * invR * (y_val_kp1 - nu_k_p1 + A_p1 * nu_k);

                    arma::vec y_i_mean = W_i * V_i;

                    arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    // int count_while_loop = 0;
                    // int count_while_loop_big = 0;
                    // while(arma::any(update_value <= 0)) {
                    //     new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                    //     update_value = Y_i_new.col(k);
                    //     update_value.elem(ind_replace) = new_value.elem(ind_replace);

                    //     count_while_loop += 1;
                    //     if(count_while_loop > 10000) {
                    //         count_while_loop_big += 1;
                    //         Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                    //         count_while_loop = 0;
                    //     }
                    //     if(count_while_loop_big > 1000) {
                    //         break;
                    //     }
                    // }

                    Y_i_new.col(k) = update_value;
                }
            }
        }

        Y_i_new = Y_i_new.t();
        newY.rows(sub_ind) = Y_i_new;
    }
    
    Y.cols(1,4) = newY;
    return Y;
}

// *** Using ***
// [[Rcpp::export]]
Rcpp::List proposal_R_cpp_new(const int nu_R, const arma::mat psi_R,const arma::mat &Y, 
                              arma::field<arma::field<arma::mat>> &Dn, 
                              const arma::field<arma::field<arma::mat>> &Xn, 
                              const arma::field <arma::vec> A, const arma::vec par, 
                              const arma::field<arma::uvec> par_index, 
                              const arma::vec EIDs, arma::field <arma::vec> B,
                              const arma::field<arma::field<arma::mat>> Dn_omega,
                              const arma::field <arma::vec> W) {
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    // -------------------------------------------------------------------------
    
    arma::vec eids = Y.col(0);
    arma::mat psi_q(4, 4, arma::fill::zeros);
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        
        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_omega_i = W(ii);

        arma::field<arma::mat> Xn_i = Xn(ii);
        arma::field<arma::mat> Dn_i = Dn(ii);
        arma::field<arma::mat> Dn_omega_i = Dn_omega(ii);
        
        for(int k = 0; k < Y_i.n_cols; k++) {
            
            arma::vec nu_k = Dn_i(k) * vec_alpha_i + Dn_omega_i(k) * vec_omega_i + Xn_i(k) * vec_beta;
            arma::vec y_diff_k = Y_i.col(k) - nu_k;
            arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            arma::mat A_k = arma::diagmat(vec_A_k);
            
            if(k == 0) {
                psi_q = psi_q + y_diff_k * y_diff_k.t();
            } else {
                arma::vec nu_k_1 = Dn_i(k-1) * vec_alpha_i + 
                                   Dn_omega_i(k-1) * vec_omega_i + 
                                   Xn_i(k-1) * vec_beta;
                arma::vec y_diff_k_1 = Y_i.col(k-1) - nu_k_1;
                
                arma::vec hold_k = y_diff_k - A_k * y_diff_k_1;
                
                psi_q = psi_q + hold_k * hold_k.t();
            }
        }
    }

    psi_q = psi_q + psi_R;
    int nu_q = Y.n_rows + nu_R;

    List nu_psi_R = List::create(psi_q, nu_q);
    return nu_psi_R;
}

// [[Rcpp::export]]
arma::mat small_impute_Y_i_cpp( const int i, const int ii, const arma::vec &par, 
                                const arma::field<arma::uvec> &par_index, 
                                const arma::vec &A, arma::mat Y,
                                const arma::field<arma::mat> &Dn, 
                                const arma::field<arma::mat> &Xn, 
                                const arma::vec &B, const arma::vec &W, 
                                const arma::field<arma::mat> &Dn_omega,
                                const arma::mat &otype, arma::vec t_pts) {
    
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameter values -------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    // Transform to support (-1,1)
    arma::vec vec_A_scale = { (exp(vec_A_total(0)) - 1) / (1+exp(vec_A_total(0))),
                              (exp(vec_A_total(1)) - 1) / (1+exp(vec_A_total(1))),
                              (exp(vec_A_total(2)) - 1) / (1+exp(vec_A_total(2))),
                              (exp(vec_A_total(3)) - 1) / (1+exp(vec_A_total(3))),
                              (exp(vec_A_total(4)) - 1) / (1+exp(vec_A_total(4))),
                              (exp(vec_A_total(5)) - 1) / (1+exp(vec_A_total(5))),
                              (exp(vec_A_total(6)) - 1) / (1+exp(vec_A_total(6))),
                              (exp(vec_A_total(7)) - 1) / (1+exp(vec_A_total(7))),
                              (exp(vec_A_total(8)) - 1) / (1+exp(vec_A_total(8))),
                              (exp(vec_A_total(9)) - 1) / (1+exp(vec_A_total(9))),
                              (exp(vec_A_total(10)) - 1) / (1+exp(vec_A_total(10))),
                              (exp(vec_A_total(11)) - 1) / (1+exp(vec_A_total(11))),
                              (exp(vec_A_total(12)) - 1) / (1+exp(vec_A_total(12))),
                              (exp(vec_A_total(13)) - 1) / (1+exp(vec_A_total(13))),
                              (exp(vec_A_total(14)) - 1) / (1+exp(vec_A_total(14))),
                              (exp(vec_A_total(15)) - 1) / (1+exp(vec_A_total(15))),
                              (exp(vec_A_total(16)) - 1) / (1+exp(vec_A_total(16))),
                              (exp(vec_A_total(17)) - 1) / (1+exp(vec_A_total(17))),
                              (exp(vec_A_total(18)) - 1) / (1+exp(vec_A_total(18))),
                              (exp(vec_A_total(19)) - 1) / (1+exp(vec_A_total(19)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE
        
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))}; // THREE STATE
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    
    // -------------------------------------------------------------------------
    
    // Subsetting the data -----------------------------------------------------
    
    arma::vec b_i = B;
    arma::field<arma::mat> Dn_ii = Dn;
    arma::field<arma::mat> Dn_omega_ii = Dn_omega;
    arma::field<arma::mat> Xn_ii = Xn;
    arma::vec vec_alpha_ii = A;
    arma::vec vec_omega_ii = W;
    
    // Index of observed versus missing data
    // 1 = observed, 0 = missing
    arma::mat otype_i = otype;
    arma::mat Y_temp = Y;
    arma::mat Y_i = Y_temp.cols(1,4);
    Y_i = Y_i.t();
    otype_i = otype_i.t();
    arma::mat Y_i_new = Y_i;
    
    // Need to impute all time points because for a given subject --------------
    // for(int k = 0; k < Y_i.n_cols; k++) {
    for(int w=0; w < t_pts.n_elem; w++) {
        int k = t_pts(w) - 1;
        
        if(all(otype_i.col(k) == 1)) {
            Y_i_new.col(k) = Y_i.col(k);
        } else {
            if(k == 0) {
                arma::vec vec_A = A_all_state.col(b_i(k) - 1);
                arma::mat Gamma = { {R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                     R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                     R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                     R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                     {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                      R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                      R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                      R(1,3) / (1 - vec_A(1) * vec_A(3))},
                                      {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                       R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                       R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                       R(2,3) / (1 - vec_A(2) * vec_A(3))},
                                       {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                        R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                        R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                        R(3,3) / (1 - vec_A(3) * vec_A(3))}};
                arma::mat inv_Gamma = arma::inv_sympd(Gamma);
                
                arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                arma::mat A_p1 = arma::diagmat(vec_A_p1);
                
                arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + Dn_omega_ii(k+1) * vec_omega_ii + Xn_ii(k+1) * vec_beta;
                
                arma::vec y_val_kp1 = Y_i_new.col(k+1);
                
                arma::mat inv_W_i = inv_Gamma + A_p1.t() * invR * A_p1;
                arma::mat W_i = inv(inv_W_i);
                arma::vec V_i = inv_Gamma*nu_k + A_p1.t()*invR*(y_val_kp1 - nu_k_p1 + A_p1*nu_k);
                
                arma::vec y_i_mean = W_i * V_i;
                
                arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                arma::vec update_value = Y_i_new.col(k);
                arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                update_value.elem(ind_replace) = new_value.elem(ind_replace);
                
                // // Prevent negatives
                // int count_while_loop = 0;
                // int count_while_loop_big = 0;
                // while(arma::any(update_value <= 0)) {
                //     new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                //     update_value = Y_i_new.col(k);
                //     update_value.elem(ind_replace) = new_value.elem(ind_replace);
                //     
                //     count_while_loop += 1;
                //     if(count_while_loop > 10000) {
                //         count_while_loop_big += 1;
                //         Rcpp::Rcout << "stuck in small impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                //         count_while_loop = 0;
                //     }
                //     if(count_while_loop_big > 1000) {
                //         break;
                //     }
                // }
                
                Y_i_new.col(k) = update_value;
            } else if(k == Y_i.n_cols - 1) {
                
                arma::vec vec_A = A_all_state.col(b_i(k) - 1);
                arma::mat A_k = arma::diagmat(vec_A);
                
                arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + Dn_omega_ii(k-1) * vec_omega_ii +  Xn_ii(k-1) * vec_beta;
                
                arma::vec y_val_km1 = Y_i_new.col(k-1);
                
                arma::vec y_i_mean = nu_k + A_k * (y_val_km1 - nu_k_m1);
                
                arma::vec new_value = arma::mvnrnd(y_i_mean, R, 1);
                arma::vec update_value = Y_i_new.col(k);
                arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                update_value.elem(ind_replace) = new_value.elem(ind_replace);
                
                // // Prevent negatives
                // int count_while_loop = 0;
                // int count_while_loop_big = 0;
                // while(arma::any(update_value <= 0)) {
                //     new_value = arma::mvnrnd(y_i_mean, R, 1);
                //     update_value = Y_i_new.col(k);
                //     update_value.elem(ind_replace) = new_value.elem(ind_replace);
                //     
                //     count_while_loop += 1;
                //     if(count_while_loop > 10000) {
                //         count_while_loop_big += 1;
                //         Rcpp::Rcout << "stuck in small impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                //         count_while_loop = 0;
                //     }
                //     if(count_while_loop_big > 1000) {
                //         break;
                //     }
                // }
                
                Y_i_new.col(k) = update_value;
            } else {
                
                arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
                arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                arma::mat A_k = arma::diagmat(vec_A_k);                    
                arma::mat A_p1 = arma::diagmat(vec_A_p1);
                
                arma::vec nu_k    = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + Dn_omega_ii(k-1) * vec_omega_ii + Xn_ii(k-1) * vec_beta;
                arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + Dn_omega_ii(k+1) * vec_omega_ii + Xn_ii(k+1) * vec_beta;
                
                arma::vec y_val_km1 = Y_i_new.col(k-1);
                arma::vec y_val_kp1 = Y_i_new.col(k+1);
                
                arma::mat inv_W_i = invR + A_p1.t() * invR * A_p1;
                arma::mat W_i = inv(inv_W_i);
                
                arma::vec V_i = invR * (nu_k + A_k * (y_val_km1 - nu_k_m1)) + 
                    A_p1.t() * invR * (y_val_kp1 - nu_k_p1 + A_p1 * nu_k);
                
                arma::vec y_i_mean = W_i * V_i;
                
                arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                arma::vec update_value = Y_i_new.col(k);
                arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                update_value.elem(ind_replace) = new_value.elem(ind_replace);
                
                // // Prevent negatives
                // int count_while_loop = 0;
                // int count_while_loop_big = 0;
                // while(arma::any(update_value <= 0)) {
                //     new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                //     update_value = Y_i_new.col(k);
                //     update_value.elem(ind_replace) = new_value.elem(ind_replace);
                //     
                //     count_while_loop += 1;
                //     if(count_while_loop > 10000) {
                //         count_while_loop_big += 1;
                //         Rcpp::Rcout << "stuck in small impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                //         count_while_loop = 0;
                //     }
                //     if(count_while_loop_big > 1000) {
                //         break;
                //     }
                // }
                
                Y_i_new.col(k) = update_value;
            }
        }
    }
    
    Y_i_new = Y_i_new.t();
    
    Y.cols(1,4) = Y_i_new;
    return Y;
}

// [[Rcpp::export]] 
Rcpp::List update_b_i_impute_cpp( const arma::vec EIDs, const arma::vec &par, 
                                  const arma::field<arma::uvec> &par_index, 
                                  const arma::field <arma::vec> &A, 
                                  arma::field <arma::vec> &B, 
                                  arma::mat &Y, const arma::mat &z, 
                                  arma::field<arma::field<arma::mat>> &Dn, 
                                  const arma::field <arma::field<arma::mat>> &Xn, 
                                  const arma::field<arma::field<arma::mat>> &Dn_omega, 
                                  const arma::field <arma::vec> &W,
                                  const arma::vec &bleed_indicator, int n_cores,
                                  const arma::mat &otype) {
    
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::vec eids = Y.col(0); 
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6); 
    arma::mat Y_return = Y;
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        
        int n_i = sub_ind.n_elem;
        
        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        // Subsetting fields
        arma::vec B_temp = B(ii);
        arma::vec A_temp = A(ii);
        arma::vec W_temp = W(ii);
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
        
        arma::field<arma::mat> Dn_temp = Dn(ii);
        arma::field<arma::mat> Dn_omega_temp = Dn_omega(ii);
        arma::field<arma::mat> Xn_temp = Xn(ii);
        
        // Subsetting the remaining data
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat otype_i = otype.rows(sub_ind);
        arma::mat z_temp = z.rows(sub_ind);
        
        for (int k = 0; k < n_i - 1; k++) {
            
            arma::vec t_pts;
            if (k == n_i - 2) {
                t_pts = arma::linspace(k+1, k+2, 2);
            } else {
                t_pts = arma::linspace(k+1, k+3, 3);
            }
            
            arma::vec pr_B = B_temp;
            arma::field<arma::mat> pr_Dn = Dn_temp;
            
            // Sample and update the two neighboring states
            arma::mat Omega_set;
            if (clinic_rule >= 0) {
                Omega_set = Omega_fun_cpp_new(k + 1, n_i, B_temp, false);
            } else {
                Omega_set = Omega_fun_cpp_new(k + 1, n_i, B_temp, true);
            }
            
            int sampled_index = arma::randi(arma::distr_param(1, Omega_set.n_rows));
            
            pr_B.rows(k, k+1) = Omega_set.row(sampled_index-1).t();
            
            // State sampling using RBC and Clinic information -----------------------
            bool valid_prop = false;
            bool b_i_rule = arma::any(arma::vectorise(pr_B)==2);
            
            if(clinic_rule == 1) {
                if(rbc_rule == 1) {
                    // clinic=1, rbc=1 -> NEED S2, *yes* restriction on time of S2
                    int pos_bleed = arma::as_scalar(arma::find(bleed_ind_i == 1));
                    arma::uvec b_i_time = arma::find(pr_B == 2);
                    if(arma::any(b_i_time <= pos_bleed)) {
                        valid_prop = true;
                    }
                } else {
                    // clinic=1, rbc=0 -> NEED S2, *no* restriction on time of S2
                    if(b_i_rule) {
                        valid_prop = true;
                    }
                }
            } else if(clinic_rule == 0) {
                if(rbc_rule == 1) {
                    // clinic=0, rbc=1 -> NEED S2, *yes* restriction on time of S2
                    int pos_bleed = arma::as_scalar(arma::find(bleed_ind_i == 1));
                    arma::uvec b_i_time = arma::find(pr_B == 2);
                    if(arma::any(b_i_time <= pos_bleed)) {
                        valid_prop = true;
                    }
                } else {
                    // clinic=0, rbc=0 -> No restrictions, consider all state seq.
                    valid_prop = true;
                }
            } else {
                // clinic=-1, rbc=1 -> evaluate likelihood anyways because S1,S4,S5
                // clinic=-1, rbc=0 -> evaluate likelihood anyways because S1,S4,S5
                valid_prop = true; 
            }
            
            // If the proposed state sequence is the same, then we do not need to 
            // evaluate the likelihood. Thus valid_prop = false can be set.
            if(arma::accu(pr_B == B_temp) == pr_B.n_elem) {
                valid_prop = false;
            }
            // -----------------------------------------------------------------------
            
            if(valid_prop) {
                
                // Assuming a valid proposal, let's compare likelihoods when imputation occurs
                arma::vec log_prev_vec = log_f_i_cpp(i, ii, t_pts, par, par_index,A_temp,
                                                     B_temp,Y_temp,z_temp,Dn_temp,Xn_temp,
                                                     Dn_omega_temp, W_temp);
                
                double log_target_prev = log_prev_vec(0);
                
                arma::vec twos(pr_B.n_elem, arma::fill::zeros);
                arma::vec threes = twos; // THREE STATE
                arma::vec fours = twos;
                arma::vec fives = twos;
                
                twos.elem(arma::find(pr_B == 2)) += 1;
                threes.elem(arma::find(pr_B == 3)) += 1; // THREE STATE
                fours.elem(arma::find(pr_B == 4)) += 1;
                fives.elem(arma::find(pr_B == 5)) += 1;
                
                arma::vec ones(pr_B.n_elem, arma::fill::ones);
                
                arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
                bigB = arma::join_rows(bigB, arma::cumsum(threes)); // THREE STATE
                bigB = arma::join_rows(bigB, arma::cumsum(fours));
                bigB = arma::join_rows(bigB, arma::cumsum(fives));
                
                arma::mat I = arma::eye(4,4);
                for(int jj = 0; jj < n_i; jj++) {
                    pr_Dn(jj) = arma::kron(I, bigB.row(jj));
                }
                
                // Impute new Y values based on the new state sequence 
                arma::mat pr_Y = small_impute_Y_i_cpp(i, ii, par, par_index, 
                                                      A_temp, Y_temp, pr_Dn, 
                                                      Xn_temp, pr_B, W_temp, 
                                                      Dn_omega_temp, otype_i,
                                                      t_pts);
                
                arma::vec log_target_vec = log_f_i_cpp(i,ii,t_pts,par,par_index,A_temp,
                                                       pr_B, pr_Y,z_temp,pr_Dn,Xn_temp,
                                                       Dn_omega_temp, W_temp);
                double log_target = log_target_vec(0);
                
                // Note that the proposal probs cancel in the MH ratio
                double diff_check = log_target - log_target_prev;
                double min_log = log(arma::randu(arma::distr_param(0,1)));
                if(diff_check > min_log){
                    B_temp = pr_B;
                    Dn_temp = pr_Dn;
                    Y_temp = pr_Y;
                }
            }
        }
        B_return(ii) = B_temp;
        Dn_return(ii) = Dn_temp;
        Y_return.rows(sub_ind) = Y_temp;
    }
    List B_Dn = List::create(B_return, Dn_return, Y_return);
    
    return B_Dn;
}

// NEW ADDITIONS ---------------------------------------------------------------
// *** Using ***
bool rule_check(int clinic_rule, int rbc_rule, arma::vec bleed_ind_i, 
                arma::vec s_i) {
    
    bool eval_like;
    
    if(clinic_rule == -1) {
        // clinic = -1, rbc = 0,1
        eval_like = true;
    } else if(clinic_rule == 1) {
        if(rbc_rule > 0) {
            // clinic = 1, rbc = 1
            if(arma::any(s_i == 2)) {
                arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
                int first_bleed_ind = arma::as_scalar(bleed_ind_ind);
                arma::uvec sim_bleed_ind = arma::find(s_i == 2);
                arma::vec check_vec = {0, s_i(first_bleed_ind)};
                if(first_bleed_ind > 0) {
                    check_vec(0) = s_i(first_bleed_ind-1);
                }
                
                if(arma::any(sim_bleed_ind <= first_bleed_ind) && 
                   arma::any(check_vec == 2)) {
                    eval_like = true;
                } else {
                    eval_like = false;
                }
            } else {
                eval_like = false;
            }
        } else {
            // clinic = 1, rbc = 0
            if(arma::any(s_i == 2)) {
                eval_like = true;
            } else {
                eval_like = false;
            }
        }
    } else {
        if(rbc_rule > 0) {
            // clinic = 0, rbc = 1
            if(arma::any(s_i == 2)) {
                arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
                int first_bleed_ind = arma::as_scalar(bleed_ind_ind);
                arma::uvec sim_bleed_ind = arma::find(s_i == 2);
                arma::vec check_vec = {0, s_i(first_bleed_ind)};
                if(first_bleed_ind > 0) {
                    check_vec(0) = s_i(first_bleed_ind-1);
                }
                
                if(arma::any(sim_bleed_ind <= first_bleed_ind) && 
                   arma::any(check_vec == 2)) {
                    eval_like = true;
                } else {
                    eval_like = false;
                }
            } else {
                eval_like = false;
            }
        } else {
            // clinic = 0, rbc = 0
            eval_like = true;
        }    
    }
    
    return eval_like;
}

// *** Using ***
// [[Rcpp::export]]
arma::field<arma::field<arma::mat>> initialize_X(const arma::vec EIDs,
                                                 const arma::mat &Y,
                                                 const arma::mat &x) {
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::vec eids = Y.col(0);
    arma::field<arma::field<arma::mat>> Xn(EIDs.n_elem);
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::field<arma::mat> Xn_i(n_i);
        
        arma::mat x_i = x.rows(sub_ind);
        
        arma::mat I = arma::eye(4,4);
        
        for(int jj = 0; jj < n_i; jj++) {
            Xn_i(jj) = arma::kron(I, x_i.row(jj));
        }
        Xn(ii) = Xn_i;
    }
    
    return Xn;
}

// *** Using ***
double log_like_i(int k, arma::mat y_i, arma::mat z_i, arma::vec b_i, 
                  arma::mat alpha_i, arma::mat omega_i, 
                  arma::mat X_beta, arma::vec D_alpha, arma::mat D_omega, 
                  arma::mat X_beta_1, arma::vec D_alpha_1, arma::mat D_omega_1,
                  arma::vec vec_beta, arma::mat A_all_state, arma::mat R, 
                  arma::mat zeta, arma::vec P_init) {
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    int b_i_k = b_i(k);
    arma::vec y_i_k = y_i.col(k);
    arma::vec nu_k = (D_alpha.t() * alpha_i).t() + 
                        D_omega * omega_i + 
                        X_beta * vec_beta;
    arma::vec vec_A = A_all_state.col(b_i_k - 1);
    
    double log_like_i_k;
    
    // alpha_i is a 5x4 matrix, omega_i is a 84x1 matrix
    if(k == 0) {
        arma::mat Gamma = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                            R(0,2) / (1 - vec_A(0) * vec_A(2)), R(0,3) / (1 - vec_A(0) * vec_A(3))},
                           {R(1,0) / (1 - vec_A(1) * vec_A(0)), R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                            R(1,2) / (1 - vec_A(1) * vec_A(2)), R(1,3) / (1 - vec_A(1) * vec_A(3))},
                           {R(2,0) / (1 - vec_A(2) * vec_A(0)), R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                            R(2,2) / (1 - vec_A(2) * vec_A(2)), R(2,3) / (1 - vec_A(2) * vec_A(3))},
                           {R(3,0) / (1 - vec_A(3) * vec_A(0)), R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                            R(3,2) / (1 - vec_A(3) * vec_A(2)), R(3,3) / (1 - vec_A(3) * vec_A(3))}};
        
        arma::vec log_y_pdf = dmvnorm(y_i_k.t(), nu_k, Gamma, true);
        log_like_i_k = log(P_init(b_i_k - 1)) + arma::as_scalar(log_y_pdf);
        
    } else {
        // State space component
        double q1_sub = arma::as_scalar(z_i.row(k) * zeta.col(0));
        double q1 = exp(q1_sub);
        double q2_sub = arma::as_scalar(z_i.row(k) * zeta.col(1));
        double q2 = exp(q2_sub);
        double q3_sub = arma::as_scalar(z_i.row(k) * zeta.col(2));
        double q3 = exp(q3_sub);
        double q4_sub = arma::as_scalar(z_i.row(k) * zeta.col(3));
        double q4 = exp(q4_sub);
        
        double q5_sub = arma::as_scalar(z_i.row(k) * zeta.col(4));
        double q5 = exp(q5_sub);
        double q6_sub = arma::as_scalar(z_i.row(k) * zeta.col(5));
        double q6 = exp(q6_sub);
        double q7_sub = arma::as_scalar(z_i.row(k) * zeta.col(6));
        double q7 = exp(q7_sub);
        double q8_sub = arma::as_scalar(z_i.row(k) * zeta.col(7));
        double q8 = exp(q8_sub);
        
        double q9_sub = arma::as_scalar(z_i.row(k) * zeta.col(8));
        double q9 = exp(q9_sub);
        double q10_sub = arma::as_scalar(z_i.row(k) * zeta.col(9));
        double q10 = exp(q10_sub);
        double q11_sub = arma::as_scalar(z_i.row(k) * zeta.col(10));
        double q11 = exp(q11_sub);
        double q12_sub = arma::as_scalar(z_i.row(k) * zeta.col(11));
        double q12 = exp(q12_sub);
        
        arma::mat Q = { {   1,   q1,  0,  q2,  0},
                        {   0,    1, q3,  q4,  0},
                        {  q5,   q6,  1,  q7,  0},
                        {   0,   q8,  0,   1, q9},
                        { q10,  q11,  0, q12,  1}};
        
        arma::vec q_row_sums = arma::sum(Q, 1);
        arma::mat P_i = Q.each_col() / q_row_sums;
        int b_i_k_1 = b_i(k-1);
        
        arma::mat A_1 = arma::diagmat(vec_A);
        
        arma::vec y_i_k_1 = y_i.col(k-1);
        arma::vec nu_k_1 = (D_alpha_1.t() * alpha_i).t() + 
                            D_omega_1 * omega_i + 
                            X_beta_1 * vec_beta;
        
        arma::vec mean_k = nu_k + A_1 * (y_i_k_1 - nu_k_1);
        
        arma::vec log_y_pdf = dmvnorm(y_i_k.t(), mean_k, R, true);
        log_like_i_k = log(P_i(b_i_k_1 - 1, b_i_k - 1)) + arma::as_scalar(log_y_pdf);
    }
    
    return log_like_i_k;
}

// *** Using ***
// [[Rcpp::export]]
Rcpp::List mle_state_seq(const arma::vec EIDs, const arma::vec &par,
                         const arma::field<arma::uvec> &par_index,
                         const arma::field <arma::vec> &A,
                         const arma::mat &Y, const arma::mat &z,
                         const arma::field <arma::field<arma::mat>> &Xn,
                         const arma::field<arma::field<arma::mat>> &Dn_omega,
                         const arma::field <arma::vec> &W, int n_cores) {
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID

    // Initializing parameters -------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);

    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);

    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);

    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12);

    arma::vec vec_init = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init(0)), exp(vec_init(1)),
                            exp(vec_init(2)), exp(vec_init(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    // -------------------------------------------------------------------------

    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);

    arma::vec eids = Y.col(0);
    arma::vec clinic_rule_vec = Y.col(6);

    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        int clinic_rule = clinic_rule_vec(sub_ind.min());

        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat y_i = Y_temp.cols(1, 4);
        y_i = y_i.t();

        arma::mat z_i = z.rows(sub_ind);
        arma::field<arma::mat> X_i = Xn(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::vec alpha_i_vec = A(ii);
        arma::mat alpha_i = arma::reshape(alpha_i_vec, 5, 4);
        arma::vec omega_i = W(ii);

        arma::vec b_i(n_i, arma::fill::zeros);
        arma::imat adj_mat_i;
        if(clinic_rule >= 0) {
            adj_mat_i = adj_mat_GLOBAL;
        } else {
            adj_mat_i = adj_mat_sub_GLOBAL;
        }

        for(int k = 0; k < n_i; k++) {

            // Shared information across k ------------
            arma::mat X_beta = X_i(k);
            arma::mat X_beta_1(X_beta.n_rows, X_beta.n_cols, arma::fill::zeros);
            arma::mat D_omega = D_omega_i(k);
            arma::mat D_omega_1(D_omega.n_rows, D_omega.n_cols, arma::fill::zeros);

            if(k==0) {
                // Initial state --------------------
                arma::vec init_vals(adj_mat_i.n_cols, arma::fill::zeros);
                for(int jj = 0; jj < init_vals.n_elem; jj++) {
                    b_i(k) = jj + 1;

                    arma::vec b_sub = b_i.subvec(0,k);

                    arma::vec twos(b_sub.n_elem, arma::fill::zeros);
                    arma::vec threes(b_sub.n_elem, arma::fill::zeros);
                    arma::vec fours(b_sub.n_elem, arma::fill::zeros);
                    arma::vec fives(b_sub.n_elem, arma::fill::zeros);

                    twos.elem(arma::find(b_sub == 2)) += 1;
                    threes.elem(arma::find(b_sub == 3)) += 1;
                    fours.elem(arma::find(b_sub == 4)) += 1;
                    fives.elem(arma::find(b_sub == 5)) += 1;

                    arma::vec D_alpha = {1, arma::accu(twos), arma::accu(threes),
                                         arma::accu(fours), arma::accu(fives)};
                    arma::vec D_alpha_1(D_alpha.n_elem, arma::fill::zeros);

                    init_vals(jj) = log_like_i(k, y_i, z_i, b_i, alpha_i,
                                               omega_i, X_beta, D_alpha,
                                               D_omega, X_beta_1, D_alpha_1,
                                               D_omega_1, vec_beta, A_all_state,
                                               R, zeta, P_init);
                }

                int selected_state = arma::index_max(init_vals) + 1;

                // Double check that bleeding isn't selected when clinic_rule=-1
                if(clinic_rule < 0) {
                    if(selected_state == 2) {
                        arma::uvec state_order = arma::sort_index(init_vals, "descend");
                        // Take second largest
                        int new_state = state_order(1)+1;
                        selected_state = new_state;
                    }
                }

                b_i(k) = selected_state;

            } else {
                // All other states -----------------
                int prev_state = b_i(k-1);
                arma::vec poss_next_state(arma::accu(adj_mat_i.row(prev_state-1)), arma::fill::zeros);
                arma::vec poss_state_like(arma::accu(adj_mat_i.row(prev_state-1)), arma::fill::zeros);
                
                int w_ind = 0;
                for(int jj = 0; jj < adj_mat_i.n_cols; jj++) {
                    if(adj_mat_i(prev_state-1, jj) != 0) {
                        poss_next_state(w_ind) = jj + 1;
                        
                        b_i(k) = jj + 1;
                        
                        arma::vec b_sub = b_i.subvec(0,k);
                        
                        arma::vec twos(b_sub.n_elem, arma::fill::zeros);
                        arma::vec threes(b_sub.n_elem, arma::fill::zeros);
                        arma::vec fours(b_sub.n_elem, arma::fill::zeros);
                        arma::vec fives(b_sub.n_elem, arma::fill::zeros);
                        
                        twos.elem(arma::find(b_sub == 2)) += 1;
                        threes.elem(arma::find(b_sub == 3)) += 1;
                        fours.elem(arma::find(b_sub == 4)) += 1;
                        fives.elem(arma::find(b_sub == 5)) += 1;
                        
                        arma::vec D_alpha = {1, arma::accu(twos), arma::accu(threes),
                                             arma::accu(fours), arma::accu(fives)};
                        arma::vec D_alpha_1 = {1, arma::accu(twos.subvec(0,k-1)), 
                                                  arma::accu(threes.subvec(0,k-1)),
                                                  arma::accu(fours.subvec(0,k-1)),
                                                  arma::accu(fives.subvec(0,k-1))};
                        X_beta_1 = X_i(k-1);
                        D_omega_1 = D_omega_i(k-1);
                        
                        
                        poss_state_like(w_ind) = log_like_i(k, y_i, z_i, b_i, alpha_i,
                                                            omega_i, X_beta, D_alpha,
                                                            D_omega, X_beta_1, D_alpha_1,
                                                            D_omega_1, vec_beta, A_all_state,
                                                            R, zeta, P_init);
                        
                        w_ind += 1;
                    }
                }
                
                b_i(k) = poss_next_state(poss_state_like.index_max());
            }
        }
        
        // Format the Dn_alpha list --------------------------------------------
        arma::field<arma::mat> Dn_temp(n_i);
        arma::vec twos(b_i.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(b_i == 2)) += 1;
        threes.elem(arma::find(b_i == 3)) += 1;
        fours.elem(arma::find(b_i == 4)) += 1;
        fives.elem(arma::find(b_i == 5)) += 1;
        
        arma::vec ones(b_i.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes)); 
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            Dn_temp(jj) = arma::kron(I, bigB.row(jj));
        }
        
        B_return(ii) = b_i;
        Dn_return(ii) = Dn_temp;
    }
    
    List B_Dn = List::create(B_return, Dn_return);
    return B_Dn;
}

// *** Using ***
// [[Rcpp::export]]
Rcpp::List mh_up(const arma::vec EIDs, const arma::vec &par, 
                 const arma::field<arma::uvec> &par_index, 
                 const arma::field <arma::vec> &A, arma::field <arma::vec> &B, 
                 const arma::mat &Y, const arma::mat &z, 
                 const arma::field <arma::field<arma::mat>> &Xn, 
                 const arma::field<arma::field<arma::mat>> &Dn_omega, 
                 const arma::field <arma::vec> &W, const arma::vec &bleed_indicator, 
                 int n_cores, int states_per_step) {
    
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init(0)), exp(vec_init(1)),
                            exp(vec_init(2)), exp(vec_init(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    // -------------------------------------------------------------------------
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    arma::vec eids = Y.col(0); 
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat y_i = Y_temp.cols(1, 4);
        y_i = y_i.t();
        
        arma::vec b_i = B(ii);
        arma::mat z_i = z.rows(sub_ind);
        arma::field<arma::mat> X_i = Xn(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::vec alpha_i_vec = A(ii);
        arma::mat alpha_i = arma::reshape(alpha_i_vec, 5, 4);
        arma::vec omega_i = W(ii);
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - states_per_step + 1; k++) {
            
            // All possible state transitions given current time point ---------
            arma::mat Omega_set;
            if(clinic_rule >= 0) {
                Omega_set = get_omega_list(k, n_i, b_i, states_per_step, false);
            } else {
                Omega_set = get_omega_list(k, n_i, b_i, states_per_step, true);
            }
            
            if(Omega_set.n_rows > 1) {
                
                // State-space proposal prob. distribution ---------------------
                arma::vec ss_prob(Omega_set.n_rows, arma::fill::ones);
                ss_prob = (1/arma::accu(ss_prob)) * ss_prob;
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
                arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                
                arma::vec s_i = b_i;
                s_i.rows(k, k+states_per_step-1) = Omega_set.row(row_ind(0)).t();
                
                double log_like_b = 0;
                double log_like_s = 0;
                
                if(arma::accu(arma::abs(s_i - b_i)) == 0) {
                    b_i = s_i;
                } else {
                    
                    // Incorporate RBC_rule & clinic_rule ----------------------
                    bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i);
                    
                    if(eval_like) {
                        // Computations for mean of current (b_i) --------------
                        arma::vec ones_b(b_i.n_elem, arma::fill::ones);
                        arma::vec twos_b(b_i.n_elem, arma::fill::zeros);
                        arma::vec threes_b(b_i.n_elem, arma::fill::zeros);
                        arma::vec fours_b(b_i.n_elem, arma::fill::zeros);
                        arma::vec fives_b(b_i.n_elem, arma::fill::zeros);
                        
                        twos_b.elem(arma::find(b_i == 2)) += 1;
                        threes_b.elem(arma::find(b_i == 3)) += 1;
                        fours_b.elem(arma::find(b_i == 4)) += 1;
                        fives_b.elem(arma::find(b_i == 5)) += 1;
                        
                        arma::mat D_alpha_b = arma::join_rows(ones_b, arma::cumsum(twos_b));
                        D_alpha_b = arma::join_rows(D_alpha_b, arma::cumsum(threes_b)); 
                        D_alpha_b = arma::join_rows(D_alpha_b, arma::cumsum(fours_b));
                        D_alpha_b = arma::join_rows(D_alpha_b, arma::cumsum(fives_b));
                        
                        // Computations for mean of candidate (s_i) ------------
                        arma::vec ones_s(s_i.n_elem, arma::fill::ones);
                        arma::vec twos_s(s_i.n_elem, arma::fill::zeros);
                        arma::vec threes_s(s_i.n_elem, arma::fill::zeros);
                        arma::vec fours_s(s_i.n_elem, arma::fill::zeros);
                        arma::vec fives_s(s_i.n_elem, arma::fill::zeros);
                        
                        twos_s.elem(arma::find(s_i == 2)) += 1;
                        threes_s.elem(arma::find(s_i == 3)) += 1;
                        fours_s.elem(arma::find(s_i == 4)) += 1;
                        fives_s.elem(arma::find(s_i == 5)) += 1;
                        
                        arma::mat D_alpha_s = arma::join_rows(ones_s, arma::cumsum(twos_s));
                        D_alpha_s = arma::join_rows(D_alpha_s, arma::cumsum(threes_s)); 
                        D_alpha_s = arma::join_rows(D_alpha_s, arma::cumsum(fours_s));
                        D_alpha_s = arma::join_rows(D_alpha_s, arma::cumsum(fives_s));
                        
                        for(int kk = k; kk < n_i; kk++) {
                            
                            // Shared components for the mean ------------------
                            arma::vec D_alpha_b_k = D_alpha_b.row(kk).t();
                            arma::vec D_alpha_b_k_1(D_alpha_b_k.n_elem, arma::fill::zeros);
                            arma::vec D_alpha_s_k = D_alpha_s.row(kk).t();
                            arma::vec D_alpha_s_k_1(D_alpha_s_k.n_elem, arma::fill::zeros);
                            
                            arma::mat X_beta = X_i(kk);
                            arma::mat X_beta_1(X_beta.n_rows, X_beta.n_cols, arma::fill::zeros);
                            arma::mat D_omega = D_omega_i(kk);
                            arma::mat D_omega_1(D_omega.n_rows, D_omega.n_cols, arma::fill::zeros);
                            
                            if(kk > 0) {
                                D_alpha_b_k_1 = D_alpha_b.row(kk-1).t();
                                D_alpha_s_k_1 = D_alpha_s.row(kk-1).t();
                                X_beta_1 = X_i(kk-1);
                                D_omega_1 = D_omega_i(kk-1);
                            }
                            
                            log_like_b = log_like_b + log_like_i(kk, y_i, z_i, b_i, alpha_i, 
                                                                 omega_i, X_beta, D_alpha_b_k, 
                                                                 D_omega, X_beta_1, D_alpha_b_k_1,
                                                                 D_omega_1, vec_beta, A_all_state,
                                                                 R, zeta, P_init);
                            
                            log_like_s = log_like_s + log_like_i(kk, y_i, z_i, s_i, alpha_i, 
                                                                 omega_i, X_beta, D_alpha_s_k, 
                                                                 D_omega, X_beta_1, D_alpha_s_k_1,
                                                                 D_omega_1, vec_beta, A_all_state,
                                                                 R, zeta, P_init);
                        }
                        double diff_check = log_like_s - log_like_b;
                        double min_log = log(arma::randu(arma::distr_param(0,1)));
                        if(diff_check > min_log){b_i = s_i;} 
                    }
                }
            }
        }
        
        // Format the Dn_alpha list --------------------------------------------
        arma::field<arma::mat> Dn_temp(n_i);
        arma::vec twos(b_i.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(b_i == 2)) += 1;
        threes.elem(arma::find(b_i == 3)) += 1;
        fours.elem(arma::find(b_i == 4)) += 1;
        fives.elem(arma::find(b_i == 5)) += 1;
        
        arma::vec ones(b_i.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes)); 
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            Dn_temp(jj) = arma::kron(I, bigB.row(jj));
        }
        
        B_return(ii) = b_i;
        Dn_return(ii) = Dn_temp;
    }
    
    List B_Dn = List::create(B_return, Dn_return);
    return B_Dn;
}

// *** Using ***
// [[Rcpp::export]]
Rcpp::List almost_gibbs_up(const arma::vec EIDs, const arma::vec &par, 
                           const arma::field<arma::uvec> &par_index, 
                           const arma::field <arma::vec> &A, arma::field <arma::vec> &B, 
                           const arma::mat &Y, const arma::mat &z, 
                           const arma::field <arma::field<arma::mat>> &Xn, 
                           const arma::field<arma::field<arma::mat>> &Dn_omega, 
                           const arma::field <arma::vec> &W, const arma::vec &bleed_indicator, 
                           int n_cores, int states_per_step) {
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init(0)), exp(vec_init(1)),
                            exp(vec_init(2)), exp(vec_init(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    // -------------------------------------------------------------------------
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    arma::vec eids = Y.col(0); 
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat y_i = Y_temp.cols(1, 4);
        y_i = y_i.t();
        
        arma::vec b_i = B(ii);
        arma::mat z_i = z.rows(sub_ind);
        arma::field<arma::mat> X_i = Xn(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::vec alpha_i_vec = A(ii);
        arma::mat alpha_i = arma::reshape(alpha_i_vec, 5, 4);
        arma::vec omega_i = W(ii);
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - states_per_step + 1; k++) {
            
            // All possible state transitions given current time point ---------
            arma::mat Omega_set;
            if(clinic_rule >= 0) {
                Omega_set = get_omega_list(k, n_i, b_i, states_per_step, false);
            } else {
                Omega_set = get_omega_list(k, n_i, b_i, states_per_step, true);
            }
            
            if(Omega_set.n_rows > 1) {
                
                // Compute proposal distribution ("Almost Gibbs") --------------
                arma::vec ss_prob(Omega_set.n_rows, arma::fill::ones);
                for(int jj = 0; jj < Omega_set.n_rows; jj++) {
                    double log_like_val = 0;
                    arma::vec ss_jj = b_i;
                    ss_jj.rows(k, k + states_per_step - 1) = Omega_set.row(jj).t();
                    
                    arma::vec ones_b(ss_jj.n_elem, arma::fill::ones);
                    arma::vec twos_b(ss_jj.n_elem, arma::fill::zeros);
                    arma::vec threes_b(ss_jj.n_elem, arma::fill::zeros);
                    arma::vec fours_b(ss_jj.n_elem, arma::fill::zeros);
                    arma::vec fives_b(ss_jj.n_elem, arma::fill::zeros);
                    
                    twos_b.elem(arma::find(ss_jj == 2)) += 1;
                    threes_b.elem(arma::find(ss_jj == 3)) += 1;
                    fours_b.elem(arma::find(ss_jj == 4)) += 1;
                    fives_b.elem(arma::find(ss_jj == 5)) += 1;
                    
                    arma::mat D_alpha = arma::join_rows(ones_b, arma::cumsum(twos_b));
                    D_alpha = arma::join_rows(D_alpha, arma::cumsum(threes_b)); 
                    D_alpha = arma::join_rows(D_alpha, arma::cumsum(fours_b));
                    D_alpha = arma::join_rows(D_alpha, arma::cumsum(fives_b));
                    
                    for(int kk = k; kk <= k + states_per_step; kk++) {
                        if(kk < n_i) {
                            arma::vec D_alpha_k = D_alpha.row(kk).t();
                            arma::vec D_alpha_k_1(D_alpha_k.n_elem, arma::fill::zeros);
                            arma::mat X_beta = X_i(kk);
                            arma::mat X_beta_1(X_beta.n_rows, X_beta.n_cols, arma::fill::zeros);
                            arma::mat D_omega = D_omega_i(kk);
                            arma::mat D_omega_1(D_omega.n_rows, D_omega.n_cols, arma::fill::zeros);
                            
                            if(kk > 0) {
                                D_alpha_k_1 = D_alpha.row(kk-1).t();
                                X_beta_1 = X_i(kk-1);
                                D_omega_1 = D_omega_i(kk-1);
                            }
                            
                            log_like_val = log_like_val + log_like_i(kk, y_i, z_i, b_i, alpha_i, 
                                                                     omega_i, X_beta, D_alpha_k, 
                                                                     D_omega, X_beta_1, D_alpha_k_1,
                                                                     D_omega_1, vec_beta, A_all_state,
                                                                     R, zeta, P_init);    
                        }
                    }
                    
                    ss_prob(jj) = log_like_val;
                }
                
                double prob_log_max = ss_prob.max();
                ss_prob = ss_prob - prob_log_max;
                ss_prob = exp(ss_prob);
                ss_prob = (1/arma::accu(ss_prob)) * ss_prob;
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
                arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                // -------------------------------------------------------------
                
                arma::vec s_i = b_i;
                s_i.rows(k, k+states_per_step-1) = Omega_set.row(row_ind(0)).t();
                
                double log_like_b = 0;
                double log_like_s = 0;
                
                if(arma::accu(arma::abs(s_i - b_i)) == 0) {
                    b_i = s_i;
                } else {
                    
                    // Incorporate RBC_rule & clinic_rule ----------------------
                    bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i);
                    
                    if(eval_like) {
                        if(k >= n_i - states_per_step - 1) {
                            // Gibbs update
                            b_i = s_i;
                        } else {
                            // Computations for mean of current (b_i) --------------
                            arma::vec ones_b(b_i.n_elem, arma::fill::ones);
                            arma::vec twos_b(b_i.n_elem, arma::fill::zeros);
                            arma::vec threes_b(b_i.n_elem, arma::fill::zeros);
                            arma::vec fours_b(b_i.n_elem, arma::fill::zeros);
                            arma::vec fives_b(b_i.n_elem, arma::fill::zeros);
                            
                            twos_b.elem(arma::find(b_i == 2)) += 1;
                            threes_b.elem(arma::find(b_i == 3)) += 1;
                            fours_b.elem(arma::find(b_i == 4)) += 1;
                            fives_b.elem(arma::find(b_i == 5)) += 1;
                            
                            arma::mat D_alpha_b = arma::join_rows(ones_b, arma::cumsum(twos_b));
                            D_alpha_b = arma::join_rows(D_alpha_b, arma::cumsum(threes_b)); 
                            D_alpha_b = arma::join_rows(D_alpha_b, arma::cumsum(fours_b));
                            D_alpha_b = arma::join_rows(D_alpha_b, arma::cumsum(fives_b));
                            
                            // Computations for mean of candidate (s_i) ------------
                            arma::vec ones_s(s_i.n_elem, arma::fill::ones);
                            arma::vec twos_s(s_i.n_elem, arma::fill::zeros);
                            arma::vec threes_s(s_i.n_elem, arma::fill::zeros);
                            arma::vec fours_s(s_i.n_elem, arma::fill::zeros);
                            arma::vec fives_s(s_i.n_elem, arma::fill::zeros);
                            
                            twos_s.elem(arma::find(s_i == 2)) += 1;
                            threes_s.elem(arma::find(s_i == 3)) += 1;
                            fours_s.elem(arma::find(s_i == 4)) += 1;
                            fives_s.elem(arma::find(s_i == 5)) += 1;
                            
                            arma::mat D_alpha_s = arma::join_rows(ones_s, arma::cumsum(twos_s));
                            D_alpha_s = arma::join_rows(D_alpha_s, arma::cumsum(threes_s)); 
                            D_alpha_s = arma::join_rows(D_alpha_s, arma::cumsum(fours_s));
                            D_alpha_s = arma::join_rows(D_alpha_s, arma::cumsum(fives_s));
                            
                            for(int kk = k + states_per_step + 1; kk < n_i; kk++) {
                                
                                // Shared components for the mean ------------------
                                arma::vec D_alpha_b_k = D_alpha_b.row(kk).t();
                                arma::vec D_alpha_b_k_1(D_alpha_b_k.n_elem, arma::fill::zeros);
                                arma::vec D_alpha_s_k = D_alpha_s.row(kk).t();
                                arma::vec D_alpha_s_k_1(D_alpha_s_k.n_elem, arma::fill::zeros);
                                
                                arma::mat X_beta = X_i(kk);
                                arma::mat X_beta_1(X_beta.n_rows, X_beta.n_cols, arma::fill::zeros);
                                arma::mat D_omega = D_omega_i(kk);
                                arma::mat D_omega_1(D_omega.n_rows, D_omega.n_cols, arma::fill::zeros);
                                
                                if(kk > 0) {
                                    D_alpha_b_k_1 = D_alpha_b.row(kk-1).t();
                                    D_alpha_s_k_1 = D_alpha_s.row(kk-1).t();
                                    X_beta_1 = X_i(kk-1);
                                    D_omega_1 = D_omega_i(kk-1);
                                }
                                
                                log_like_b = log_like_b + log_like_i(kk, y_i, z_i, b_i, alpha_i, 
                                                                     omega_i, X_beta, D_alpha_b_k, 
                                                                     D_omega, X_beta_1, D_alpha_b_k_1,
                                                                     D_omega_1, vec_beta, A_all_state,
                                                                     R, zeta, P_init);
                                
                                log_like_s = log_like_s + log_like_i(kk, y_i, z_i, s_i, alpha_i, 
                                                                     omega_i, X_beta, D_alpha_s_k, 
                                                                     D_omega, X_beta_1, D_alpha_s_k_1,
                                                                     D_omega_1, vec_beta, A_all_state,
                                                                     R, zeta, P_init);
                            }
                            double diff_check = log_like_s - log_like_b;
                            double min_log = log(arma::randu(arma::distr_param(0,1)));
                            if(diff_check > min_log){b_i = s_i;} 
                        }
                    }
                }
            }
        }
        
        // Format the Dn_alpha list --------------------------------------------
        arma::field<arma::mat> Dn_temp(n_i);
        arma::vec twos(b_i.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(b_i == 2)) += 1;
        threes.elem(arma::find(b_i == 3)) += 1;
        fours.elem(arma::find(b_i == 4)) += 1;
        fives.elem(arma::find(b_i == 5)) += 1;
        
        arma::vec ones(b_i.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes)); 
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            Dn_temp(jj) = arma::kron(I, bigB.row(jj));
        }
        
        B_return(ii) = b_i;
        Dn_return(ii) = Dn_temp;
    }
    
    List B_Dn = List::create(B_return, Dn_return);
    return B_Dn;
}

// *** Using ***
// [[Rcpp::export]]
Rcpp::List gibbs_up(const arma::vec EIDs, const arma::vec &par, 
                    const arma::field<arma::uvec> &par_index, 
                    const arma::field <arma::vec> &A, arma::field <arma::vec> &B, 
                    const arma::mat &Y, const arma::mat &z, 
                    const arma::field <arma::field<arma::mat>> &Xn, 
                    const arma::field<arma::field<arma::mat>> &Dn_omega, 
                    const arma::field <arma::vec> &W, const arma::vec &bleed_indicator, 
                    int n_cores, int states_per_step) {
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init(0)), exp(vec_init(1)),
                            exp(vec_init(2)), exp(vec_init(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    // -------------------------------------------------------------------------
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    arma::vec eids = Y.col(0); 
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat y_i = Y_temp.cols(1, 4);
        y_i = y_i.t();
        
        arma::vec b_i = B(ii);
        arma::mat z_i = z.rows(sub_ind);
        arma::field<arma::mat> X_i = Xn(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::vec alpha_i_vec = A(ii);
        arma::mat alpha_i = arma::reshape(alpha_i_vec, 5, 4);
        arma::vec omega_i = W(ii);
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - states_per_step + 1; k++) {
            
            // All possible state transitions given current time point ---------
            arma::mat Omega_set;
            if(clinic_rule >= 0) {
                Omega_set = get_omega_list(k, n_i, b_i, states_per_step, false);
            } else {
                Omega_set = get_omega_list(k, n_i, b_i, states_per_step, true);
            } 
            
            if(Omega_set.n_rows > 1) {
                
                // Compute proposal distribution ("Gibbs") ---------------------
                arma::vec ss_prob(Omega_set.n_rows, arma::fill::ones);
                for(int jj = 0; jj < Omega_set.n_rows; jj++) {
                    double log_like_val = 0;
                    arma::vec ss_jj = b_i;
                    ss_jj.rows(k, k + states_per_step - 1) = Omega_set.row(jj).t();
                    
                    arma::vec ones_b(ss_jj.n_elem, arma::fill::ones);
                    arma::vec twos_b(ss_jj.n_elem, arma::fill::zeros);
                    arma::vec threes_b(ss_jj.n_elem, arma::fill::zeros);
                    arma::vec fours_b(ss_jj.n_elem, arma::fill::zeros);
                    arma::vec fives_b(ss_jj.n_elem, arma::fill::zeros);
                    
                    twos_b.elem(arma::find(ss_jj == 2)) += 1;
                    threes_b.elem(arma::find(ss_jj == 3)) += 1;
                    fours_b.elem(arma::find(ss_jj == 4)) += 1;
                    fives_b.elem(arma::find(ss_jj == 5)) += 1;
                    
                    arma::mat D_alpha = arma::join_rows(ones_b, arma::cumsum(twos_b));
                    D_alpha = arma::join_rows(D_alpha, arma::cumsum(threes_b)); 
                    D_alpha = arma::join_rows(D_alpha, arma::cumsum(fours_b));
                    D_alpha = arma::join_rows(D_alpha, arma::cumsum(fives_b));
                    
                    for(int kk = k; kk < n_i; kk++) {
                        arma::vec D_alpha_k = D_alpha.row(kk).t();
                        arma::vec D_alpha_k_1(D_alpha_k.n_elem, arma::fill::zeros);
                        arma::mat X_beta = X_i(kk);
                        arma::mat X_beta_1(X_beta.n_rows, X_beta.n_cols, arma::fill::zeros);
                        arma::mat D_omega = D_omega_i(kk);
                        arma::mat D_omega_1(D_omega.n_rows, D_omega.n_cols, arma::fill::zeros);
                        
                        if(kk > 0) {
                            D_alpha_k_1 = D_alpha.row(kk-1).t();
                            X_beta_1 = X_i(kk-1);
                            D_omega_1 = D_omega_i(kk-1);
                        }
                        
                        log_like_val = log_like_val + log_like_i(kk, y_i, z_i, b_i, alpha_i, 
                                                                 omega_i, X_beta, D_alpha_k, 
                                                                 D_omega, X_beta_1, D_alpha_k_1,
                                                                 D_omega_1, vec_beta, A_all_state,
                                                                 R, zeta, P_init);    
                    }
                    
                    ss_prob(jj) = log_like_val;
                } 
                
                double prob_log_max = ss_prob.max();
                ss_prob = ss_prob - prob_log_max;
                ss_prob = exp(ss_prob);
                ss_prob = (1/arma::accu(ss_prob)) * ss_prob;
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
                arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                // -------------------------------------------------------------
                
                arma::vec s_i = b_i;
                s_i.rows(k, k+states_per_step-1) = Omega_set.row(row_ind(0)).t();
                
                if(arma::accu(arma::abs(s_i - b_i)) == 0) {
                    b_i = s_i;
                } else {
                    // Incorporate RBC_rule & clinic_rule ----------------------
                    bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i);
                    if(eval_like) {b_i = s_i;}
                }
            }
        }
        
        // Format the Dn_alpha list --------------------------------------------
        arma::field<arma::mat> Dn_temp(n_i);
        arma::vec twos(b_i.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(b_i == 2)) += 1;
        threes.elem(arma::find(b_i == 3)) += 1;
        fours.elem(arma::find(b_i == 4)) += 1;
        fives.elem(arma::find(b_i == 5)) += 1;
        
        arma::vec ones(b_i.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes)); 
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            Dn_temp(jj) = arma::kron(I, bigB.row(jj));
        }
        
        B_return(ii) = b_i;
        Dn_return(ii) = Dn_temp;
    }
    
    List B_Dn = List::create(B_return, Dn_return);
    return B_Dn;
}

// *** Using ***
// [[Rcpp::export]]
Rcpp::List mh_up_all(const arma::vec EIDs, const arma::vec &par, 
                     const arma::field<arma::uvec> &par_index, 
                     const arma::field <arma::vec> &A, arma::field <arma::vec> &B, 
                     const arma::mat &Y, const arma::mat &z, 
                     const arma::field <arma::field<arma::mat>> &Xn, 
                     const arma::field<arma::field<arma::mat>> &Dn_omega, 
                     const arma::field <arma::vec> &W, const arma::vec &bleed_indicator, 
                     int n_cores) {
    
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init(0)), exp(vec_init(1)),
                            exp(vec_init(2)), exp(vec_init(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    // -------------------------------------------------------------------------
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    arma::vec eids = Y.col(0); 
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6);
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat y_i = Y_temp.cols(1, 4);
        y_i = y_i.t();
        
        arma::vec b_i = B(ii);
        arma::vec s_i(n_i, arma::fill::zeros);
        arma::mat z_i = z.rows(sub_ind);
        arma::field<arma::mat> X_i = Xn(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::vec alpha_i_vec = A(ii);
        arma::mat alpha_i = arma::reshape(alpha_i_vec, 5, 4);
        arma::vec omega_i = W(ii);
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
        
        arma::imat adj_mat_i;
        if(clinic_rule >= 0) {
            adj_mat_i = adj_mat_GLOBAL;
        } else {
            adj_mat_i = adj_mat_sub_GLOBAL;
        } 
        
        arma::vec all_like_vals_b(n_i, arma::fill::zeros);
        arma::vec all_like_vals_s(n_i, arma::fill::zeros);
        
        // Propose a full state sequence ---------------------------------------
        for (int k = 0; k < n_i; k++) {
            
            arma::vec like_vals_s(adj_mat_i.n_cols, arma::fill::zeros);
            arma::vec like_vals_b(adj_mat_i.n_cols, arma::fill::zeros);
            
            for(int m = 0; m < adj_mat_i.n_cols; m++) {
                
                // Computations for mean of candidate (s_i) ------------
                arma::vec s_temp = s_i.subvec(0, k);
                s_temp(k) = m+1;
                arma::vec ones_s(s_temp.n_elem, arma::fill::ones);
                arma::vec twos_s(s_temp.n_elem, arma::fill::zeros);
                arma::vec threes_s(s_temp.n_elem, arma::fill::zeros);
                arma::vec fours_s(s_temp.n_elem, arma::fill::zeros);
                arma::vec fives_s(s_temp.n_elem, arma::fill::zeros);
                
                twos_s.elem(arma::find(s_temp == 2)) += 1;
                threes_s.elem(arma::find(s_temp == 3)) += 1;
                fours_s.elem(arma::find(s_temp == 4)) += 1;
                fives_s.elem(arma::find(s_temp == 5)) += 1;
                
                // Computations for mean of current (b_i) --------------
                arma::vec b_temp = b_i.subvec(0, k);
                arma::vec ones_b(b_temp.n_elem, arma::fill::ones);
                arma::vec twos_b(b_temp.n_elem, arma::fill::zeros);
                arma::vec threes_b(b_temp.n_elem, arma::fill::zeros);
                arma::vec fours_b(b_temp.n_elem, arma::fill::zeros);
                arma::vec fives_b(b_temp.n_elem, arma::fill::zeros);
                
                twos_b.elem(arma::find(b_temp == 2)) += 1;
                threes_b.elem(arma::find(b_temp == 3)) += 1;
                fours_b.elem(arma::find(b_temp == 4)) += 1;
                fives_b.elem(arma::find(b_temp == 5)) += 1;
                
                arma::vec D_alpha_s_k = {1, arma::accu(twos_s), arma::accu(threes_s),
                                         arma::accu(fours_s), arma::accu(fives_s)};
                arma::vec D_alpha_s_k_1(D_alpha_s_k.n_elem, arma::fill::zeros);
                
                arma::vec D_alpha_b_k = {1, arma::accu(twos_b), arma::accu(threes_b),
                                         arma::accu(fours_b), arma::accu(fives_b)};
                arma::vec D_alpha_b_k_1(D_alpha_b_k.n_elem, arma::fill::zeros);
                
                arma::mat X_beta = X_i(k);
                arma::mat X_beta_1(X_beta.n_rows, X_beta.n_cols, arma::fill::zeros);
                arma::mat D_omega = D_omega_i(k);
                arma::mat D_omega_1(D_omega.n_rows, D_omega.n_cols, arma::fill::zeros);
                
                if(k > 0) {
                    D_alpha_s_k_1 = {1,  arma::accu(twos_s.subvec(0,k-1)), 
                                         arma::accu(threes_s.subvec(0,k-1)),
                                         arma::accu(fours_s.subvec(0,k-1)), 
                                         arma::accu(fives_s.subvec(0,k-1))};
                    D_alpha_b_k_1 = {1,  arma::accu(twos_b.subvec(0,k-1)), 
                                         arma::accu(threes_b.subvec(0,k-1)),
                                         arma::accu(fours_b.subvec(0,k-1)), 
                                         arma::accu(fives_b.subvec(0,k-1))};
                    X_beta_1 = X_i(k-1);
                    D_omega_1 = D_omega_i(k-1);
                }
                
                arma::vec s_i_curr = s_i;
                s_i_curr(k) = m+1;
                
                if(k == 0) {
                    double log_like_s = log_like_i(k, y_i, z_i, s_i_curr, alpha_i, 
                                                   omega_i, X_beta, D_alpha_s_k, 
                                                   D_omega, X_beta_1, D_alpha_s_k_1,
                                                   D_omega_1, vec_beta, A_all_state,
                                                   R, zeta, P_init);
                    like_vals_s(m) = log_like_s;
                    like_vals_b(m) = like_vals_s(m);
                } else {
                    double log_like_s = log_like_i(k, y_i, z_i, s_i_curr, alpha_i, 
                                                   omega_i, X_beta, D_alpha_s_k, 
                                                   D_omega, X_beta_1, D_alpha_s_k_1,
                                                   D_omega_1, vec_beta, A_all_state,
                                                   R, zeta, P_init);
                    double log_like_b = log_like_i(k, y_i, z_i, b_i, alpha_i, 
                                                   omega_i, X_beta, D_alpha_b_k, 
                                                   D_omega, X_beta_1, D_alpha_b_k_1,
                                                   D_omega_1, vec_beta, A_all_state,
                                                   R, zeta, P_init);
                    like_vals_s(m) = log_like_s;
                    like_vals_b(m) = log_like_b;
                }
            }
            
            all_like_vals_s(k) = arma::accu(exp(like_vals_s));
            all_like_vals_b(k) = arma::accu(exp(like_vals_b));
            
            // Determine sampling distribution for s_i -------------------------
            arma::vec ss_ind = arma::linspace(0, adj_mat_i.n_cols-1, adj_mat_i.n_cols);
            
            double prob_log_max = like_vals_s.max();
            like_vals_s = like_vals_s - prob_log_max;
            like_vals_s = exp(like_vals_s);
            arma::vec ss_prob = (1/arma::accu(like_vals_s)) * like_vals_s;
            
            arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
            s_i(k) = row_ind(0) + 1;
        }
        
        double diff_check = arma::accu(log(all_like_vals_s)) - arma::accu(log(all_like_vals_b));
        double min_log = log(arma::randu(arma::distr_param(0,1)));
        if(diff_check > min_log){
            bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i);
            if(eval_like) {
                b_i = s_i;   
            }
        }
        
        // Format the Dn_alpha list --------------------------------------------
        arma::field<arma::mat> Dn_temp(n_i);
        arma::vec twos(b_i.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(b_i == 2)) += 1;
        threes.elem(arma::find(b_i == 3)) += 1;
        fours.elem(arma::find(b_i == 4)) += 1;
        fives.elem(arma::find(b_i == 5)) += 1;
        
        arma::vec ones(b_i.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes)); 
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            Dn_temp(jj) = arma::kron(I, bigB.row(jj));
        }
        
        B_return(ii) = b_i;
        Dn_return(ii) = Dn_temp;
    }
    
    List B_Dn = List::create(B_return, Dn_return);
    return B_Dn;
}

// [[Rcpp::export]]
void test_fnc() {
    
    arma::field<arma::vec> test_A(3);
    test_A(0) = {1,2,3,4,5};
    test_A(1) = {6,7,8,9,10};
    test_A(2) = {0,3,6,9,12};
    
    Rcpp::Rcout << test_A << std::endl;
    

    // arma::vec vec_A_mean = {1.5,  1.5,  1.5,  1.5, -1.0, -1.0, -1.0, -1.0,
    //                         0.1,  0.1,  0.1,  0.1,    0,    0,    0,    0,
    //                         0.1,  0.1,  0.1,  0.1};
    // 
    // arma::vec A_scale = exp(vec_A_mean) / (1 + exp(vec_A_mean));
    // Rcpp::Rcout << A_scale << std::endl;
    // 
    // arma::vec a = {1,2,-1,3,0,8};
    // arma::vec b = arma::sort(a, "descend");
    // arma::uvec c = arma::sort_index(a, "descend");
    // 
    // Rcpp::Rcout << a << std::endl;
    // Rcpp::Rcout << b << std::endl;
    // Rcpp::Rcout << c+1 << std::endl;
    // 
    // int clinic_rule = 0;
    // int rbc_rule = 1;
    // 
    // arma::vec s_i = {4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,2,2,2,2,2,2,2,2,2,
    //                  2,2,2,2,4,4,4,4,5,5,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
    //                  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,
    //                  3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
    //                  2,2,2,2,4,4,4,4,4,4,4,4,5,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
    //                  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
    //                  3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4};
    // arma::vec bleed_ind_i = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
    //                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    // 
    // if((clinic_rule < 0) | (rbc_rule >= 1)) {
    //     Rcpp::Rcout << "hello worlds " << std::endl;
    // }
    // 
    // 
    // bool eval_like;
    // 
    // // Incorporate RBC_rule & clinic_rule ----------------------
    // if(clinic_rule == -1) {
    //     // clinic = -1, rbc = 0,1
    //     eval_like = true;
    // } else if(clinic_rule == 1) {
    //     if(rbc_rule > 0) {
    //         // clinic = 1, rbc = 1
    //         if(arma::any(s_i == 2)) {
    //             arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
    //             int first_bleed_ind = arma::as_scalar(bleed_ind_ind);
    //             arma::uvec sim_bleed_ind = arma::find(s_i == 2);
    //             
    //             arma::vec check_vec = {0, s_i(first_bleed_ind)};
    //             if(first_bleed_ind > 0) {
    //                 check_vec(0) = s_i(first_bleed_ind-1);
    //             }
    //             
    //             if(arma::any(sim_bleed_ind <= first_bleed_ind) && arma::any(check_vec == 2)) {
    //                 eval_like = true;
    //             } else {
    //                 eval_like = false;
    //             }
    //         } else {
    //             eval_like = false;
    //         }
    //     } else {
    //         // clinic = 1, rbc = 0
    //         if(arma::any(s_i == 2)) {
    //             eval_like = true;
    //         } else {
    //             eval_like = false;
    //         }
    //     }
    // } else {
    //     if(rbc_rule > 0) {
    //         Rcpp::Rcout << "here1" << std::endl;
    //         // clinic = 0, rbc = 1
    //         if(arma::any(s_i == 2)) {
    //             Rcpp::Rcout << "here2" << std::endl;
    //             arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
    //             int first_bleed_ind = arma::as_scalar(bleed_ind_ind);
    //             arma::uvec sim_bleed_ind = arma::find(s_i == 2);
    //             
    //             arma::vec check_vec = {0, s_i(first_bleed_ind)};
    //             if(first_bleed_ind > 0) {
    //                 check_vec(0) = s_i(first_bleed_ind-1);
    //             }
    //             
    //             Rcpp::Rcout << sim_bleed_ind.t() << std::endl;
    //             Rcpp::Rcout << first_bleed_ind << std::endl;
    //             Rcpp::Rcout << check_vec << std::endl;
    //             Rcpp::Rcout << arma::any(sim_bleed_ind <= first_bleed_ind) << std::endl;
    //             Rcpp::Rcout << arma::any(check_vec == 2) << std::endl;
    //             
    //             if(arma::any(sim_bleed_ind <= first_bleed_ind) && arma::any(check_vec == 2)) {
    //                 Rcpp::Rcout << "here3" << std::endl;
    //                 eval_like = true;
    //             } else {
    //                 eval_like = false;
    //             }
    //         } else {
    //             eval_like = false;
    //         }
    //     } else {
    //         // clinic = 0, rbc = 0
    //         eval_like = true;
    //     }    
    // }
    // 
    // Rcpp::Rcout << eval_like << std::endl;
    
} 
