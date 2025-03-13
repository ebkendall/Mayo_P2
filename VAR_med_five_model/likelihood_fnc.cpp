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
void initialize_cpp(arma::imat a_mat, arma::imat a_mat_sub, int s_per_s, 
                    int sampling_num) {
    adj_mat_GLOBAL = a_mat;
    adj_mat_sub_GLOBAL = a_mat_sub;
    states_per_step_GLOBAL = s_per_s;
    
    if(sampling_num == 4 || sampling_num == 5) {
        Omega_List_GLOBAL_multi = get_Omega_list(adj_mat_GLOBAL, 2);    
        Omega_List_GLOBAL_sub_multi = get_Omega_list(adj_mat_sub_GLOBAL, 2);
    } else {
        Omega_List_GLOBAL_multi = get_Omega_list(adj_mat_GLOBAL, states_per_step_GLOBAL);    
        Omega_List_GLOBAL_sub_multi = get_Omega_list(adj_mat_sub_GLOBAL, states_per_step_GLOBAL);
    }
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
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
            
            // arma::vec vec_A = A_all_state.col(b_0 - 1);
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
            
            // arma::vec vec_A = A_all_state.col(b_k - 1);
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

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
    arma::vec vec_A_mean(vec_A_content.n_elem, arma::fill::zeros);
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    // -------------------------------------------------------------------------
    
    arma::vec eids = Y.col(0);
    
    arma::field<arma::vec> A(EIDs.n_elem);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
            
            // arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            // arma::mat A_1_k = arma::diagmat(vec_A_k);
            arma::mat A_1_k = arma::diagmat(vec_A);
            
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            arma::mat diff_hold_y = Y_i.col(k) - A_1_k * Y_i.col(k-1);
            
            arma::vec m_k = diff_hold_y - diff_hold_d_o*vec_omega_ii - diff_hold_x * vec_beta;
            
            interm_W = interm_W + diff_hold_d.t() * invR * diff_hold_d;
            
            interm_V = interm_V + diff_hold_d.t() * invR * m_k;
        }
        
        // arma::vec vec_A = A_all_state.col(b_i(0) - 1);
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    // -------------------------------------------------------------------------
    
    arma::vec eids = Y.col(0);
    
    arma::field<arma::vec> W(EIDs.n_elem);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
            
            // arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            // arma::mat A_1_k = arma::diagmat(vec_A_k);
            arma::mat A_1_k = arma::diagmat(vec_A);
            
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            arma::mat diff_hold_y = Y_i.col(k) - A_1_k * Y_i.col(k-1);
            
            arma::vec m_k = diff_hold_y - diff_hold_d*vec_alpha_ii - diff_hold_x*vec_beta;
            
            interm_W = interm_W + diff_hold_d_o.t() * invR * diff_hold_d_o;
            
            interm_V = interm_V + diff_hold_d_o.t() * invR * m_k;
        }
        
        // arma::vec vec_A = A_all_state.col(b_i(0) - 1);
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    
    arma::vec vec_alpha_tilde = par.elem(par_index(1) - 1);
    // -------------------------------------------------------------------------
    
    arma::vec eids = Y.col(0);

    arma::field<arma::vec> V_list(EIDs.n_elem);
    arma::field<arma::mat> inv_W_list(EIDs.n_elem);
    arma::field<arma::mat> in_Up_list(EIDs.n_elem);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
            
            // arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            // arma::mat A_1_k = arma::diagmat(vec_A_k);
            arma::mat A_1_k = arma::diagmat(vec_A);
            
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            
            arma::vec m_k = Y_i.col(k) - A_1_k * Y_i.col(k-1) 
                                - diff_hold_d * vec_alpha_i - diff_hold_d_o * vec_omega_i;
            
            interm_W = interm_W + diff_hold_x.t() * invR * diff_hold_x;
            interm_V = interm_V + diff_hold_x.t() * invR * m_k;
        }
        
        // arma::vec vec_A = A_all_state.col(b_i(0) - 1);
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
        
        V_list(ii) = V_i;
        inv_W_list(ii) = W_i_inv;
        in_Up_list(ii) = in_Upsilon_cov;
    }
    
    arma::vec sum_in_V = V_list(0);
    arma::mat sum_in_inv_W = inv_W_list(0);
    arma::mat sum_in_Upsilon_cov = in_Up_list(0);
    
    for(int ii = 1; ii < EIDs.n_elem; ii++){
        sum_in_V = sum_in_V + V_list(ii);
        sum_in_inv_W = sum_in_inv_W + inv_W_list(ii);
        sum_in_Upsilon_cov = sum_in_Upsilon_cov + in_Up_list(ii);
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

// *** Using ***
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

    // Initializing parameters -------------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    // -------------------------------------------------------------------------
    
    arma::mat newY(Y.n_rows, 4); 
    arma::vec eids = Y.col(0);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {	
        
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);

        arma::field<arma::mat> Dn_ii = Dn(ii);
        arma::field<arma::mat> Dn_omega_ii = Dn_omega(ii);
        arma::field<arma::mat> Xn_ii = Xn(ii);
        arma::vec vec_alpha_ii = A(ii);
        arma::vec vec_omega_ii = W(ii);
        arma::vec b_i = B(ii);

        // Index of observed versus missing data
        // 1 = observed, 0 = missing
        arma::mat otype_i = otype.rows(sub_ind);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        otype_i = otype_i.t();
        
        arma::mat Y_i_new = Y_i;

        for(int k = 0; k < Y_i.n_cols; k++) {
            if(arma::all(otype_i.col(k) == 1)) {
                Y_i_new.col(k) = Y_i.col(k);
            } else {
                if(k == 0) {
                    // arma::vec vec_A = A_all_state.col(b_i(k) - 1);
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
                    
                    // arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                    // arma::mat A_p1 = arma::diagmat(vec_A_p1);
                    arma::mat A_p1 = arma::diagmat(vec_A);
                    
                    arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + 
                        Dn_omega_ii(k) * vec_omega_ii + 
                        Xn_ii(k) * vec_beta;
                    arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + 
                        Dn_omega_ii(k+1) * vec_omega_ii + 
                        Xn_ii(k+1) * vec_beta;
                    
                    arma::vec y_val_kp1 = Y_i_new.col(k+1);

                    arma::mat inv_W_i = inv_Gamma + A_p1.t() * invR * A_p1;
                    arma::mat W_i = arma::inv_sympd(inv_W_i);
                    
                    arma::vec V_i = inv_Gamma*nu_k + A_p1.t()*invR*(y_val_kp1 - nu_k_p1 + A_p1*nu_k);

                    arma::vec y_i_mean = W_i * V_i;

                    arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    int count_while_loop = 0;
                    int count_while_loop_big = 0;
                    while(arma::any(update_value <= 0)) {
                            new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                            update_value = Y_i_new.col(k);
                            update_value.elem(ind_replace) = new_value.elem(ind_replace);
                            count_while_loop += 1;
                            if(count_while_loop > 10000) {
                                count_while_loop_big += 1;
                                Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                                count_while_loop = 0;
                            }
                            if(count_while_loop_big > 1000) {
                                break;
                            }
                    }
                    
                    Y_i_new.col(k) = update_value;
                    
                } else if(k == Y_i.n_cols - 1) {
                    
                    // arma::vec vec_A = A_all_state.col(b_i(k) - 1);
                    arma::mat A_k = arma::diagmat(vec_A);
                    
                    arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + 
                        Dn_omega_ii(k) * vec_omega_ii + 
                        Xn_ii(k) * vec_beta;
                    arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + 
                        Dn_omega_ii(k-1) * vec_omega_ii +  
                        Xn_ii(k-1) * vec_beta;
                    
                    arma::vec y_val_km1 = Y_i_new.col(k-1);

                    arma::vec y_i_mean = nu_k + A_k * (y_val_km1 - nu_k_m1);

                    arma::vec new_value = arma::mvnrnd(y_i_mean, R, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    int count_while_loop = 0;
                    int count_while_loop_big = 0;
                    while(arma::any(update_value <= 0)) {
                        new_value = arma::mvnrnd(y_i_mean, R, 1);
                        update_value = Y_i_new.col(k);
                        update_value.elem(ind_replace) = new_value.elem(ind_replace);
                        count_while_loop += 1;
                        if(count_while_loop > 10000) {
                            count_while_loop_big += 1;
                            Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                            count_while_loop = 0;
                        }
                        if(count_while_loop_big > 1000) {
                            break;
                        }
                    }
                    
                    Y_i_new.col(k) = update_value;
                    
                } else {
                    
                    // arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
                    // arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                    // arma::mat A_k = arma::diagmat(vec_A_k);                    
                    // arma::mat A_p1 = arma::diagmat(vec_A_p1);
                    arma::mat A_k = arma::diagmat(vec_A);
                    arma::mat A_p1 = arma::diagmat(vec_A);

                    arma::vec nu_k    = Dn_ii(k) * vec_alpha_ii + 
                        Dn_omega_ii(k) * vec_omega_ii + 
                        Xn_ii(k) * vec_beta;
                    arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + 
                        Dn_omega_ii(k-1) * vec_omega_ii + 
                        Xn_ii(k-1) * vec_beta;
                    arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + 
                        Dn_omega_ii(k+1) * vec_omega_ii + 
                        Xn_ii(k+1) * vec_beta;
                    
                    arma::vec y_val_km1 = Y_i_new.col(k-1);
                    arma::vec y_val_kp1 = Y_i_new.col(k+1);

                    arma::mat inv_W_i = invR + A_p1.t() * invR * A_p1;
                    arma::mat W_i = arma::inv_sympd(inv_W_i);

                    arma::vec V_i = invR * (nu_k + A_k * (y_val_km1 - nu_k_m1)) + 
                                        A_p1.t() * invR * (y_val_kp1 - nu_k_p1 + A_p1 * nu_k);

                    arma::vec y_i_mean = W_i * V_i;

                    arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    int count_while_loop = 0;
                    int count_while_loop_big = 0;
                    while(arma::any(update_value <= 0)) {
                        new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                        update_value = Y_i_new.col(k);
                        update_value.elem(ind_replace) = new_value.elem(ind_replace);
                        count_while_loop += 1;
                        if(count_while_loop > 10000) {
                            count_while_loop_big += 1;
                            Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                            count_while_loop = 0;
                        }
                        if(count_while_loop_big > 1000) {
                            break;
                        }
                    }
                    
                    Y_i_new.col(k) = update_value;
                }
            }
        }

        newY.rows(sub_ind) = Y_i_new.t();
    }
    
    Y.cols(1,4) = newY;
    return Y;
}

// *** Using ***
// [[Rcpp::export]]
Rcpp::List proposal_R_cpp_new(const int nu_R, const arma::mat psi_R, arma::mat curr_R,
                              const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn, 
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::mat curr_R_sqrt = arma::sqrtmat_sympd(curr_R);
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
            // arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            // arma::mat A_k = arma::diagmat(vec_A_k);
            arma::mat A_k = arma::diagmat(vec_A);
            
            if(k == 0) {
                arma::mat Gamma = {{curr_R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                    curr_R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                    curr_R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                    curr_R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                   {curr_R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                    curr_R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                    curr_R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                    curr_R(1,3) / (1 - vec_A(1) * vec_A(3))},
                                   {curr_R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                    curr_R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                    curr_R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                    curr_R(2,3) / (1 - vec_A(2) * vec_A(3))},
                                   {curr_R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                    curr_R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                    curr_R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                    curr_R(3,3) / (1 - vec_A(3) * vec_A(3))}};
                
                arma::mat inv_Gamma = arma::inv_sympd(Gamma);
                arma::mat inv_gamma_sqrt = arma::sqrtmat_sympd(inv_Gamma);
                
                arma::vec y_diff_temp = curr_R_sqrt * inv_gamma_sqrt * y_diff_k;
                
                psi_q = psi_q + y_diff_temp * y_diff_temp.t();

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

// NEW ADDITIONS ---------------------------------------------------------------
// *** Using ***
bool rule_check(int clinic_rule, int rbc_rule, arma::vec bleed_ind_i, 
                arma::vec s_i, int states_per_step, int k) {
    
    bool eval_like;
    
    if(states_per_step == -1) {
        // Full state sequence update ------------------------------------------
        if(clinic_rule == -1) {
            // clinic = -1, rbc = 0,1
            eval_like = true;
        } else if(clinic_rule == 1) {
            if(rbc_rule > 0) {
                // clinic = 1, rbc = 1
                arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
                int first_bleed_ind = arma::as_scalar(bleed_ind_ind);
                arma::vec check_vec = {0, s_i(first_bleed_ind)};
                if(first_bleed_ind > 0) {
                    check_vec(0) = s_i(first_bleed_ind-1);
                }
                
                if(arma::any(check_vec == 2)) {
                    eval_like = true;
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
                arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
                int first_bleed_ind = arma::as_scalar(bleed_ind_ind);
                arma::vec check_vec = {0, s_i(first_bleed_ind)};
                if(first_bleed_ind > 0) {
                    check_vec(0) = s_i(first_bleed_ind-1);
                }
                
                if(arma::any(check_vec == 2)) {
                    eval_like = true;
                } else { 
                    eval_like = false;
                }
            } else {
                // clinic = 0, rbc = 0
                eval_like = true;
            }    
        } 
    } else {
        // states_per_step update ----------------------------------------------
        if(clinic_rule == -1) {
            // clinic = -1, rbc = 0,1
            eval_like = true;
        } else if(clinic_rule == 1) {
            if(rbc_rule > 0) {
                // clinic = 1, rbc = 1
                arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
                double first_bleed_ind = arma::as_scalar(bleed_ind_ind);
                arma::vec check_vec = {0, s_i(first_bleed_ind)};
                arma::vec check_vec_ind = {0, first_bleed_ind};
                if(first_bleed_ind > 0) {
                    check_vec(0) = s_i(first_bleed_ind-1);
                    check_vec_ind(0) = first_bleed_ind - 1;
                }
                
                arma::vec k_inds = arma::linspace(k, k + states_per_step - 1, states_per_step);
                arma::vec bleed_k_int = arma::intersect(k_inds, check_vec_ind);
                
                if(bleed_k_int.n_elem > 0) {
                    if(arma::any(check_vec == 2)) {
                        eval_like = true;
                    } else {
                        eval_like = false;
                    }
                } else {
                    eval_like = true;
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
                arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
                double first_bleed_ind = arma::as_scalar(bleed_ind_ind);
                arma::vec check_vec = {0, s_i(first_bleed_ind)};
                arma::vec check_vec_ind = {0, first_bleed_ind};
                if(first_bleed_ind > 0) {
                    check_vec(0) = s_i(first_bleed_ind-1);
                    check_vec_ind(0) = first_bleed_ind - 1;
                }
                
                arma::vec k_inds = arma::linspace(k, k + states_per_step - 1, states_per_step);
                arma::vec bleed_k_int = arma::intersect(k_inds, check_vec_ind);
                
                if(bleed_k_int.n_elem > 0) {
                    if(arma::any(check_vec == 2)) {
                        eval_like = true;
                    } else {
                        eval_like = false;
                    }
                } else {
                    eval_like = true;
                }
            } else {
                // clinic = 0, rbc = 0
                eval_like = true;
            }    
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
arma::vec log_like_i(int k, arma::mat y_i, arma::mat z_i, arma::vec b_i, 
                     arma::mat alpha_i, arma::mat omega_i, 
                     arma::mat X_beta, arma::vec D_alpha, arma::mat D_omega, 
                     arma::mat X_beta_1, arma::vec D_alpha_1, arma::mat D_omega_1,
                     arma::vec vec_beta, arma::mat A_all_state, arma::mat R, 
                     arma::mat zeta, arma::vec P_init, bool log_l) {
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    int b_i_k = b_i(k);
    arma::vec y_i_k = y_i.col(k);
    arma::vec nu_k = (D_alpha.t() * alpha_i).t() + 
                        D_omega * omega_i + 
                        X_beta * vec_beta;
    // arma::vec vec_A = A_all_state.col(b_i_k - 1);
    arma::vec vec_A = A_all_state;
    
    arma::vec log_like_i_k(2, arma::fill::zeros);
    
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
        
        if(log_l) {
            arma::vec log_y_pdf = dmvnorm(y_i_k.t(), nu_k, Gamma, true);
            log_like_i_k(0) = log(P_init(b_i_k - 1));
            log_like_i_k(1) = arma::as_scalar(log_y_pdf);    
        } else {
            arma::vec log_y_pdf = dmvnorm(y_i_k.t(), nu_k, Gamma, false);
            log_like_i_k(0) = P_init(b_i_k - 1);
            log_like_i_k(1) = arma::as_scalar(log_y_pdf); 
        }
        
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
        
        if(log_l) {
            arma::vec log_y_pdf = dmvnorm(y_i_k.t(), mean_k, R, true);
            log_like_i_k(0) = log(P_i(b_i_k_1 - 1, b_i_k - 1));
            log_like_i_k(1) = arma::as_scalar(log_y_pdf);    
        } else {
            arma::vec log_y_pdf = dmvnorm(y_i_k.t(), mean_k, R, false);
            log_like_i_k(0) = P_i(b_i_k_1 - 1, b_i_k - 1);
            log_like_i_k(1) = arma::as_scalar(log_y_pdf);
        }
    }
    
    return log_like_i_k;
}

// *** Using ***
double log_like_i_sub(int k, arma::mat y_i, arma::mat z_i, arma::vec b_i,
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
    y_i_k = y_i_k.subvec(1,2);
    
    arma::vec nu_k = (D_alpha.t() * alpha_i).t() +
        D_omega * omega_i +
        X_beta * vec_beta;
    nu_k = nu_k.subvec(1,2);
    
    // arma::vec vec_A = A_all_state.col(b_i_k - 1);
    arma::vec vec_A = A_all_state;
    
    double log_like_i_k;
    

    // alpha_i is a 5x4 matrix, omega_i is a 84x1 matrix
    if(k == 0) {
        arma::mat Gamma = {{R(1,1) / (1 - vec_A(1) * vec_A(1)),R(1,2) / (1 - vec_A(1) * vec_A(2))},
                           {R(2,1) / (1 - vec_A(2) * vec_A(1)),R(2,2) / (1 - vec_A(2) * vec_A(2))}};

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
        
        vec_A = vec_A.subvec(1,2);
        arma::mat A_1 = arma::diagmat(vec_A);
        
        R = R.submat(1,1,2,2);

        arma::vec y_i_k_1 = y_i.col(k-1);
        y_i_k_1 = y_i_k_1.subvec(1,2);
        
        arma::vec nu_k_1 = (D_alpha_1.t() * alpha_i).t() +
            D_omega_1 * omega_i +
            X_beta_1 * vec_beta;
        nu_k_1 = nu_k_1.subvec(1,2);

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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);

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

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
                    
                    arma::vec log_like_temp = log_like_i(k, y_i, z_i, b_i, alpha_i,
                                                         omega_i, X_beta, D_alpha,
                                                         D_omega, X_beta_1, D_alpha_1,
                                                         D_omega_1, vec_beta, vec_A,
                                                         R, zeta, P_init, true);

                    init_vals(jj) = arma::accu(log_like_temp);
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
                        
                        arma::vec log_like_temp2 = log_like_i(k, y_i, z_i, b_i, alpha_i,
                                                              omega_i, X_beta, D_alpha,
                                                              D_omega, X_beta_1, D_alpha_1,
                                                              D_omega_1, vec_beta, vec_A,
                                                              R, zeta, P_init, true);
                        poss_state_like(w_ind) = arma::accu(log_like_temp2);
                        
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
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
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
                
                if(arma::accu(arma::abs(s_i - b_i)) != 0) {
                    // Local RBC_rule & clinic_rule ----------------------------
                    bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i,
                                                states_per_step, k);

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
                        
                        double log_prob_diff = 0;
                        double log_like_diff = 0;

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
                            
                            arma::vec log_like_b_temp = log_like_i(kk, y_i, z_i, b_i, alpha_i,
                                                                   omega_i, X_beta, D_alpha_b_k,
                                                                   D_omega, X_beta_1, D_alpha_b_k_1,
                                                                   D_omega_1, vec_beta, vec_A,
                                                                   R, zeta, P_init, true);
                            arma::vec log_like_s_temp = log_like_i(kk, y_i, z_i, s_i, alpha_i,
                                                                   omega_i, X_beta, D_alpha_s_k,
                                                                   D_omega, X_beta_1, D_alpha_s_k_1,
                                                                   D_omega_1, vec_beta, vec_A,
                                                                   R, zeta, P_init, true);

                            if(kk <= k + states_per_step) {
                                log_prob_diff = log_prob_diff + 
                                    log_like_s_temp(0) - log_like_b_temp(0);
                            }
                            
                            log_like_diff = log_like_diff + 
                                log_like_s_temp(1) - log_like_b_temp(1);
                        }
                        
                        double diff_check = log_prob_diff + log_like_diff;
                        
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
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
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
                
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
                    
                    for(int kk = k; kk < k + states_per_step; kk++) {
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
                        
                        arma::vec log_like_temp = log_like_i(kk, y_i, z_i, ss_jj, alpha_i, 
                                                             omega_i, X_beta, D_alpha_k, 
                                                             D_omega, X_beta_1, D_alpha_k_1,
                                                             D_omega_1, vec_beta, vec_A,
                                                             R, zeta, P_init, true);
                        
                        log_like_val = log_like_val + arma::accu(log_like_temp);
                    }
                    
                    ss_prob(jj) = log_like_val;
                }
                
                double prob_log_max = ss_prob.max();
                ss_prob = ss_prob - prob_log_max;
                ss_prob = exp(ss_prob);
                ss_prob = (1/arma::accu(ss_prob)) * ss_prob;
                
                // Compute the MH Ratio to accept or reject --------------------
                arma::vec s_i = b_i;
                
                arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                s_i.rows(k, k+states_per_step-1) = Omega_set.row(row_ind(0)).t();
                
                if(arma::accu(arma::abs(s_i - b_i)) != 0) {
                    // Local RBC_rule & clinic_rule ----------------------------
                    bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i,
                                                states_per_step, k);

                    if(eval_like) {
                        // Computations for mean of current (b_i) ----------
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
                        
                        double log_prob_diff = 0;
                        double log_like_diff = 0;
                        
                        for(int kk = k + states_per_step; kk < n_i; kk++) {
                            
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
                            
                            arma::vec log_like_b_temp = log_like_i(kk, y_i, z_i, b_i, alpha_i,
                                                                   omega_i, X_beta, D_alpha_b_k,
                                                                   D_omega, X_beta_1, D_alpha_b_k_1,
                                                                   D_omega_1, vec_beta, vec_A,
                                                                   R, zeta, P_init, true);
                            arma::vec log_like_s_temp = log_like_i(kk, y_i, z_i, s_i, alpha_i,
                                                                   omega_i, X_beta, D_alpha_s_k,
                                                                   D_omega, X_beta_1, D_alpha_s_k_1,
                                                                   D_omega_1, vec_beta, vec_A,
                                                                   R, zeta, P_init, true);
                            
                            if(kk == k + states_per_step) {
                                log_prob_diff = log_prob_diff +
                                    log_like_s_temp(0) - log_like_b_temp(0);
                            }
                            
                            log_like_diff = log_like_diff +
                                log_like_s_temp(1) - log_like_b_temp(1);
                            
                        }
                        
                        double diff_check = log_prob_diff + log_like_diff;
                        
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
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
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
                
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
                        
                        arma::vec log_like_temp = log_like_i(kk, y_i, z_i, ss_jj, alpha_i, 
                                                             omega_i, X_beta, D_alpha_k, 
                                                             D_omega, X_beta_1, D_alpha_k_1,
                                                             D_omega_1, vec_beta, vec_A,
                                                             R, zeta, P_init, true);
                        
                        log_like_val = log_like_val + arma::accu(log_like_temp);
                    }
                    
                    ss_prob(jj) = log_like_val;
                } 
                
                double prob_log_max = ss_prob.max();
                ss_prob = ss_prob - prob_log_max;
                ss_prob = exp(ss_prob);
                ss_prob = (1/arma::accu(ss_prob)) * ss_prob;
                
                arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                // -------------------------------------------------------------
                
                arma::vec s_i = b_i;
                s_i.rows(k, k+states_per_step-1) = Omega_set.row(row_ind(0)).t();
                
                if(arma::accu(arma::abs(s_i - b_i)) != 0) {
                    // Local RBC_rule & clinic_rule ----------------------------
                    bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i,
                                                states_per_step, k);
                    
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

arma::vec full_seq_update(int n_i, arma::mat y_i, arma::imat adj_mat_i, 
                          arma::vec b_i, arma::vec s_i, arma::mat z_i,
                          arma::field<arma::mat> X_i, arma::field<arma::mat> D_omega_i,
                          arma::mat alpha_i, arma::vec omega_i, arma::vec bleed_ind_i,
                          arma::mat vec_beta, arma::mat A_all_state, arma::mat R,
                          arma::mat zeta, arma::vec P_init, int clinic_rule, int rbc_rule) {
    
    arma::vec all_like_vals_b(n_i, arma::fill::zeros);
    arma::vec all_like_vals_s(n_i, arma::fill::zeros);
    
    // Propose a full state sequence ---------------------------------------
    for (int k = 0; k < n_i; k++) {
        
        arma::vec ss_ind;
        arma::vec like_vals_s;
        arma::vec like_vals_b;
        
        if(k == 0) {
            like_vals_s.zeros(adj_mat_i.n_cols);
            like_vals_b.zeros(adj_mat_i.n_cols);
            
            for(int m = 0; m < adj_mat_i.n_cols; m++) {
                
                // Computations for mean of candidate (s_i) ------------
                arma::vec s_temp = s_i.subvec(0, k);
                s_temp(k) = m+1;
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
                b_temp(k) = m+1;
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
                
                arma::vec s_i_curr = s_i;
                arma::vec b_i_curr = b_i;
                s_i_curr(k) = m+1;
                b_i_curr(k) = m+1;
                
                arma::vec log_like_s_temp = log_like_i(k, y_i, z_i, s_i_curr, alpha_i,
                                                       omega_i, X_beta, D_alpha_s_k,
                                                       D_omega, X_beta_1, D_alpha_s_k_1,
                                                       D_omega_1, vec_beta, A_all_state,
                                                       R, zeta, P_init, true);
                arma::vec log_like_b_temp = log_like_i(k, y_i, z_i, b_i_curr, alpha_i,
                                                       omega_i, X_beta, D_alpha_b_k,
                                                       D_omega, X_beta_1, D_alpha_b_k_1,
                                                       D_omega_1, vec_beta, A_all_state,
                                                       R, zeta, P_init, true);
                
                like_vals_s(m) = arma::accu(log_like_s_temp);
                like_vals_b(m) = arma::accu(log_like_b_temp);
            }
            
            ss_ind = arma::linspace(0, adj_mat_i.n_cols-1, adj_mat_i.n_cols);
            
            if(clinic_rule < 0) {
                // Remove state 2 as a possibility
                like_vals_s = {like_vals_s(0), like_vals_s(2), 
                               like_vals_s(3), like_vals_s(4)};
                like_vals_b = {like_vals_b(0), like_vals_b(2), 
                               like_vals_b(3), like_vals_b(4)};
                ss_ind = {ss_ind(0), ss_ind(2), ss_ind(3), ss_ind(4)};
            }
            
        } else {
            int prev_state_s = s_i(k-1);
            int prev_state_b = b_i(k-1);
            
            like_vals_s.zeros(arma::accu(adj_mat_i.row(prev_state_s-1)));
            like_vals_b.zeros(arma::accu(adj_mat_i.row(prev_state_b-1)));
            
            arma::vec state_vals_s(arma::accu(adj_mat_i.row(prev_state_s-1)), arma::fill::zeros);
            arma::vec state_vals_b(arma::accu(adj_mat_i.row(prev_state_b-1)), arma::fill::zeros);
            
            int s_ind = 0;
            int b_ind = 0;
            
            for(int m = 0; m < adj_mat_i.n_cols; m++) {
                
                // Computations for mean of candidate (s_i) ------------
                arma::vec s_temp = s_i.subvec(0, k);
                s_temp(k) = m+1;
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
                b_temp(k) = m+1;
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
                arma::vec D_alpha_s_k_1= {1,  arma::accu(twos_s.subvec(0,k-1)),
                                          arma::accu(threes_s.subvec(0,k-1)),
                                          arma::accu(fours_s.subvec(0,k-1)),
                                          arma::accu(fives_s.subvec(0,k-1))};
                
                arma::vec D_alpha_b_k = {1, arma::accu(twos_b), arma::accu(threes_b),
                                         arma::accu(fours_b), arma::accu(fives_b)};
                arma::vec D_alpha_b_k_1 = {1,  arma::accu(twos_b.subvec(0,k-1)),
                                           arma::accu(threes_b.subvec(0,k-1)),
                                           arma::accu(fours_b.subvec(0,k-1)),
                                           arma::accu(fives_b.subvec(0,k-1))};
                
                arma::mat X_beta = X_i(k);
                arma::mat X_beta_1 = X_i(k-1);
                arma::mat D_omega = D_omega_i(k);
                arma::mat D_omega_1 = D_omega_i(k-1);
                
                arma::vec s_i_curr = s_i;
                arma::vec b_i_curr = b_i;
                s_i_curr(k) = m+1;
                b_i_curr(k) = m+1;
                
                if(adj_mat_i(prev_state_s-1, m) != 0) {
                    arma::vec log_like_s_temp = log_like_i(k, y_i, z_i, s_i_curr, alpha_i,
                                                           omega_i, X_beta, D_alpha_s_k,
                                                           D_omega, X_beta_1, D_alpha_s_k_1,
                                                           D_omega_1, vec_beta, A_all_state,
                                                           R, zeta, P_init, true);    
                    like_vals_s(s_ind) = arma::accu(log_like_s_temp);
                    state_vals_s(s_ind) = m;
                    s_ind = s_ind + 1;
                }
                
                if(adj_mat_i(prev_state_b-1, m) != 0) {
                    arma::vec log_like_b_temp = log_like_i(k, y_i, z_i, b_i_curr, alpha_i,
                                                           omega_i, X_beta, D_alpha_b_k,
                                                           D_omega, X_beta_1, D_alpha_b_k_1,
                                                           D_omega_1, vec_beta, A_all_state,
                                                           R, zeta, P_init, true);
                    like_vals_b(b_ind) = arma::accu(log_like_b_temp);    
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
    if(diff_check > min_log){
        bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i, -1, -1);
        if(eval_like) {
            b_i = s_i;   
        }
    }
    
    return b_i;
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); 
    
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
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
        if(clinic_rule < 0) {
            adj_mat_i = adj_mat_sub_GLOBAL;
        } else {
            adj_mat_i = adj_mat_GLOBAL;
        }
        
        arma::vec new_b_i = full_seq_update(n_i, y_i, adj_mat_i, b_i, s_i, z_i,
                                            X_i, D_omega_i, alpha_i, omega_i, 
                                            bleed_ind_i, vec_beta, vec_A, 
                                            R, zeta, P_init, clinic_rule, rbc_rule);
        b_i = new_b_i;
        
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
Rcpp::List almost_gibbs_fast_b(const arma::vec EIDs, const arma::vec &par, 
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);

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

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
        
        arma::imat adj_mat_i;
        bool clinic_bool;
        if(clinic_rule < 0) {
            adj_mat_i = adj_mat_sub_GLOBAL;
            clinic_bool = true;
        } else {
            adj_mat_i = adj_mat_GLOBAL;
            clinic_bool = false;
        }

        if(states_per_step >= n_i) {
            arma::vec s_i(n_i, arma::fill::zeros);
            arma::vec new_b_i = full_seq_update(n_i, y_i, adj_mat_i, b_i, s_i, z_i,
                                                X_i, D_omega_i, alpha_i, omega_i, 
                                                bleed_ind_i, vec_beta, vec_A, 
                                                R, zeta, P_init, clinic_rule, rbc_rule);
            b_i = new_b_i;
        } else {
            int k = 0;
            while(k < n_i - 2) {
                arma::vec s_i = b_i;
                
                arma::vec all_like_vals_b(states_per_step - 2, arma::fill::zeros);
                arma::vec all_like_vals_s(states_per_step - 2, arma::fill::zeros);
                
                // Step 1: sample 1-at-a-time --------------------------------------
                arma::ivec t_lim = {k + states_per_step - 2, n_i - 2};
                int t_max = arma::min(t_lim);
                for(int t = k; t < t_max; t++) {
                    arma::vec ss_ind;
                    arma::vec like_vals_s;
                    arma::vec like_vals_b;
                    
                    if(t == 0) {
                        like_vals_s.zeros(adj_mat_i.n_cols);
                        like_vals_b.zeros(adj_mat_i.n_cols);
                        
                        for(int m = 0; m < adj_mat_i.n_cols; m++) {
                            
                            // Computations for mean of candidate (s_i) ------------
                            arma::vec s_temp = s_i.subvec(0, t);
                            s_temp(t) = m+1;
                            arma::vec twos_s(s_temp.n_elem, arma::fill::zeros);
                            arma::vec threes_s(s_temp.n_elem, arma::fill::zeros);
                            arma::vec fours_s(s_temp.n_elem, arma::fill::zeros);
                            arma::vec fives_s(s_temp.n_elem, arma::fill::zeros);
                            
                            twos_s.elem(arma::find(s_temp == 2)) += 1;
                            threes_s.elem(arma::find(s_temp == 3)) += 1;
                            fours_s.elem(arma::find(s_temp == 4)) += 1;
                            fives_s.elem(arma::find(s_temp == 5)) += 1;
                            
                            // Computations for mean of current (b_i) --------------
                            arma::vec b_temp = b_i.subvec(0, t);
                            b_temp(t) = m+1;
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
                            
                            arma::mat X_beta = X_i(t);
                            arma::mat X_beta_1(X_beta.n_rows, X_beta.n_cols, arma::fill::zeros);
                            arma::mat D_omega = D_omega_i(t);
                            arma::mat D_omega_1(D_omega.n_rows, D_omega.n_cols, arma::fill::zeros);
                            
                            arma::vec s_i_curr = s_i;
                            arma::vec b_i_curr = b_i;
                            s_i_curr(t) = m+1;
                            b_i_curr(t) = m+1;
                            
                            arma::vec log_like_s_temp = log_like_i(t, y_i, z_i, s_i_curr, alpha_i,
                                                                   omega_i, X_beta, D_alpha_s_k,
                                                                   D_omega, X_beta_1, D_alpha_s_k_1,
                                                                   D_omega_1, vec_beta, vec_A,
                                                                   R, zeta, P_init, true);
                            arma::vec log_like_b_temp = log_like_i(t, y_i, z_i, b_i_curr, alpha_i,
                                                                   omega_i, X_beta, D_alpha_b_k,
                                                                   D_omega, X_beta_1, D_alpha_b_k_1,
                                                                   D_omega_1, vec_beta, vec_A,
                                                                   R, zeta, P_init, true);
                            
                            like_vals_s(m) = arma::accu(log_like_s_temp);
                            like_vals_b(m) = arma::accu(log_like_b_temp);
                        }
                        
                        ss_ind = arma::linspace(0, adj_mat_i.n_cols-1, adj_mat_i.n_cols);
                        
                        if(clinic_rule < 0) {
                            // Remove state 2 as a possibility
                            like_vals_s = {like_vals_s(0), like_vals_s(2), 
                                           like_vals_s(3), like_vals_s(4)};
                            like_vals_b = {like_vals_b(0), like_vals_b(2), 
                                           like_vals_b(3), like_vals_b(4)};
                            ss_ind = {ss_ind(0), ss_ind(2), ss_ind(3), ss_ind(4)};
                        }
                    } else {
                        int prev_state_s = s_i(t-1);
                        int prev_state_b = b_i(t-1);
                        
                        like_vals_s.zeros(arma::accu(adj_mat_i.row(prev_state_s-1)));
                        like_vals_b.zeros(arma::accu(adj_mat_i.row(prev_state_b-1)));
                        
                        arma::vec state_vals_s(arma::accu(adj_mat_i.row(prev_state_s-1)), arma::fill::zeros);
                        arma::vec state_vals_b(arma::accu(adj_mat_i.row(prev_state_b-1)), arma::fill::zeros);
                        
                        int s_ind = 0;
                        int b_ind = 0;
                        
                        for(int m = 0; m < adj_mat_i.n_cols; m++) {
                            
                            // Computations for mean of candidate (s_i) ------------
                            arma::vec s_temp = s_i.subvec(0, t);
                            s_temp(t) = m+1;
                            arma::vec twos_s(s_temp.n_elem, arma::fill::zeros);
                            arma::vec threes_s(s_temp.n_elem, arma::fill::zeros);
                            arma::vec fours_s(s_temp.n_elem, arma::fill::zeros);
                            arma::vec fives_s(s_temp.n_elem, arma::fill::zeros);
                            
                            twos_s.elem(arma::find(s_temp == 2)) += 1;
                            threes_s.elem(arma::find(s_temp == 3)) += 1;
                            fours_s.elem(arma::find(s_temp == 4)) += 1;
                            fives_s.elem(arma::find(s_temp == 5)) += 1;
                            
                            // Computations for mean of current (b_i) --------------
                            arma::vec b_temp = b_i.subvec(0, t);
                            b_temp(t) = m+1;
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
                            arma::vec D_alpha_s_k_1= {1,  arma::accu(twos_s.subvec(0,t-1)),
                                                      arma::accu(threes_s.subvec(0,t-1)),
                                                      arma::accu(fours_s.subvec(0,t-1)),
                                                      arma::accu(fives_s.subvec(0,t-1))};
                            
                            arma::vec D_alpha_b_k = {1, arma::accu(twos_b), arma::accu(threes_b),
                                                     arma::accu(fours_b), arma::accu(fives_b)};
                            arma::vec D_alpha_b_k_1 = {1,  arma::accu(twos_b.subvec(0,t-1)),
                                                       arma::accu(threes_b.subvec(0,t-1)),
                                                       arma::accu(fours_b.subvec(0,t-1)),
                                                       arma::accu(fives_b.subvec(0,t-1))};
                            
                            arma::mat X_beta = X_i(t);
                            arma::mat X_beta_1 = X_i(t-1);
                            arma::mat D_omega = D_omega_i(t);
                            arma::mat D_omega_1 = D_omega_i(t-1);
                            
                            arma::vec s_i_curr = s_i;
                            arma::vec b_i_curr = b_i;
                            s_i_curr(t) = m+1;
                            b_i_curr(t) = m+1;
                            
                            if(adj_mat_i(prev_state_s-1, m) != 0) {
                                arma::vec log_like_s_temp = log_like_i(t, y_i, z_i, s_i_curr, alpha_i,
                                                                       omega_i, X_beta, D_alpha_s_k,
                                                                       D_omega, X_beta_1, D_alpha_s_k_1,
                                                                       D_omega_1, vec_beta, vec_A,
                                                                       R, zeta, P_init, true);    
                                like_vals_s(s_ind) = arma::accu(log_like_s_temp);
                                state_vals_s(s_ind) = m;
                                s_ind = s_ind + 1;
                            }
                            
                            if(adj_mat_i(prev_state_b-1, m) != 0) {
                                arma::vec log_like_b_temp = log_like_i(t, y_i, z_i, b_i_curr, alpha_i,
                                                                       omega_i, X_beta, D_alpha_b_k,
                                                                       D_omega, X_beta_1, D_alpha_b_k_1,
                                                                       D_omega_1, vec_beta, vec_A,
                                                                       R, zeta, P_init, true);
                                like_vals_b(b_ind) = arma::accu(log_like_b_temp);    
                                state_vals_b(b_ind) = m;
                                b_ind = b_ind + 1;
                            }
                        }
                        
                        ss_ind = state_vals_s;
                    }
                    
                    all_like_vals_s(t-k) = arma::accu(exp(like_vals_s));
                    all_like_vals_b(t-k) = arma::accu(exp(like_vals_b));
                    
                    // Determine sampling distribution for s_i ---------------------
                    double prob_log_max = like_vals_s.max();
                    like_vals_s = like_vals_s - prob_log_max;
                    like_vals_s = exp(like_vals_s);
                    arma::vec ss_prob = (1/arma::accu(like_vals_s)) * like_vals_s;
                    
                    arma::vec row_ind = RcppArmadillo::sample(ss_ind, 1, false, ss_prob);
                    s_i(t) = row_ind(0) + 1;   
                }
                
                // Step 2: sample the last 2 times together ------------------------
                arma::mat Omega_set_s = get_omega_list(t_max, n_i, s_i, 2, clinic_bool);
                arma::mat Omega_set_b = get_omega_list(t_max, n_i, b_i, 2, clinic_bool);
                
                arma::vec prob_omega_s(Omega_set_s.n_rows, arma::fill::ones);
                prob_omega_s = (1/arma::accu(prob_omega_s)) * prob_omega_s;
                arma::vec ind_omega_s = arma::linspace(0, Omega_set_s.n_rows-1, Omega_set_s.n_rows);
                arma::vec row_omega_s = RcppArmadillo::sample(ind_omega_s, 1, false, prob_omega_s);
                
                s_i.rows(t_max, t_max + 1) = Omega_set_s.row(row_omega_s(0)).t();
                
                // Step 3: compute MH-ratio to accept/reject -----------------------
                if(arma::accu(arma::abs(s_i - b_i)) != 0) {
                    // Local RBC_rule & clinic_rule ----------------------------
                    bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i,
                                                states_per_step, k);
                    
                    if(eval_like) {
                        // Computations for mean of current (b_i) ----------
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
                        
                        double log_prob_diff = 0;
                        double log_like_diff = 0;
                        
                        for(int t = k + states_per_step - 2; t < n_i; t++) {
                            
                            // Shared components for the mean ------------------
                            arma::vec D_alpha_b_k = D_alpha_b.row(t).t();
                            arma::vec D_alpha_b_k_1(D_alpha_b_k.n_elem, arma::fill::zeros);
                            arma::vec D_alpha_s_k = D_alpha_s.row(t).t();
                            arma::vec D_alpha_s_k_1(D_alpha_s_k.n_elem, arma::fill::zeros);
                            
                            arma::mat X_beta = X_i(t);
                            arma::mat X_beta_1(X_beta.n_rows, X_beta.n_cols, arma::fill::zeros);
                            arma::mat D_omega = D_omega_i(t);
                            arma::mat D_omega_1(D_omega.n_rows, D_omega.n_cols, arma::fill::zeros);
                            
                            if(t > 0) {
                                D_alpha_b_k_1 = D_alpha_b.row(t-1).t();
                                D_alpha_s_k_1 = D_alpha_s.row(t-1).t();
                                X_beta_1 = X_i(t-1);
                                D_omega_1 = D_omega_i(t-1);
                            }
                            
                            arma::vec log_like_b_temp = log_like_i(t, y_i, z_i, b_i, alpha_i,
                                                                   omega_i, X_beta, D_alpha_b_k,
                                                                   D_omega, X_beta_1, D_alpha_b_k_1,
                                                                   D_omega_1, vec_beta, vec_A,
                                                                   R, zeta, P_init, true);
                            arma::vec log_like_s_temp = log_like_i(t, y_i, z_i, s_i, alpha_i,
                                                                   omega_i, X_beta, D_alpha_s_k,
                                                                   D_omega, X_beta_1, D_alpha_s_k_1,
                                                                   D_omega_1, vec_beta, vec_A,
                                                                   R, zeta, P_init, true);
                            if(t < k + states_per_step + 1) {
                                log_prob_diff = log_prob_diff +
                                    log_like_s_temp(0) - log_like_b_temp(0);
                            }
                            
                            log_like_diff = log_like_diff +
                                log_like_s_temp(1) - log_like_b_temp(1);
                            
                        }
                        
                        double log_diff_1 = arma::accu(log(all_like_vals_s)) - arma::accu(log(all_like_vals_b));
                        double log_diff_2 = log(Omega_set_s.n_rows) - log(Omega_set_b.n_rows);
                        
                        double diff_check = log_prob_diff + log_like_diff + log_diff_1 + log_diff_2;
                        
                        double min_log = log(arma::randu(arma::distr_param(0,1)));
                        if(diff_check > min_log){b_i = s_i;}
                    }
                }
                
                k = k + states_per_step - 2;
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
Rcpp::List initialize_Y(const arma::vec &EIDs, const arma::vec &par, 
                        const arma::field<arma::uvec> &par_index, 
                        const arma::field <arma::vec> &A, arma::mat Y,
                        const arma::mat &z, const arma::field <arma::field<arma::mat>> &Xn, 
                        const arma::field<arma::field<arma::mat>> &Dn_omega,
                        const arma::field <arma::vec> &W, const arma::mat &otype,
                        int n_cores, arma::vec vital_means) {
    
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);

    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);

    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::vec std_R = sqrt(arma::diagvec(R));

    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12);

    arma::vec vec_init = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init(0)), exp(vec_init(1)),
                            exp(vec_init(2)), exp(vec_init(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    // -------------------------------------------------------------------------

    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    arma::mat newY(Y.n_rows, 4); 
    arma::vec eids = Y.col(0);
    arma::vec clinic_rule_vec = Y.col(6);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
    
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        arma::mat z_i = z.rows(sub_ind);
        arma::field<arma::mat> X_i = Xn(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::vec alpha_i_vec = A(ii);
        arma::mat alpha_i = arma::reshape(alpha_i_vec, 5, 4);
        arma::vec omega_i = W(ii);
        arma::vec b_i(n_i, arma::fill::zeros);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        
        arma::mat Y_i_new = Y_i;
        
        arma::imat adj_mat_i;
        if(clinic_rule >= 0) {
            adj_mat_i = adj_mat_GLOBAL;
        } else { 
            adj_mat_i = adj_mat_sub_GLOBAL;
        } 
        
        // Index of observed versus missing data
        // 1 = observed, 0 = missing
        arma::mat otype_i = otype.rows(sub_ind);
        otype_i = otype_i.t();
        
        // STEP 1: -------------------------------------------------------------
        // Linearly interpolate HEM, HR, MAP, LACT w/ noise
        for(int h = 0; h <= 3; h++) {
            arma::uvec obs_indices = arma::find(otype_i.row(h) == 1);
            arma::uvec miss_indices = arma::find(otype_i.row(h) == 0);
            
            if(miss_indices.n_elem == n_i) {
                // No observations
                arma::vec h_vec(n_i, arma::fill::value(vital_means(h)));
                // h_vec = h_vec + std_R(h) * arma::randn<arma::vec>(n_i);
                Y_i_new.row(h) = h_vec.t();
            } else if(obs_indices.n_elem == 1){
                // One observation
                double obs_val = Y_i_new(h, obs_indices(0));
                
                arma::vec h_vec(n_i, arma::fill::value(obs_val));
                // h_vec = h_vec + std_R(h) * arma::randn<arma::vec>(n_i);
                // h_vec(obs_indices(0)) = obs_val;
                
                Y_i_new.row(h) = h_vec.t();
            } else {
                for(int j = 0; j < miss_indices.n_elem; j++) {
                    int ind_j = miss_indices(j); 
                    if(ind_j < obs_indices(0)) {
                        
                        Y_i_new(h, ind_j) = Y_i_new(h, obs_indices(0));// + R::rnorm(0, std_R(h));
                        
                    } else if(ind_j > obs_indices(obs_indices.n_elem - 1)) {
                        
                        Y_i_new(h, ind_j) = Y_i_new(h, obs_indices(obs_indices.n_elem - 1));// + R::rnorm(0, std_R(h));
                        
                    } else {
                        arma::uvec end_pts_ind = {arma::find(obs_indices < ind_j).max(),
                                                  arma::find(obs_indices > ind_j).min()};
                        arma::uvec end_pts = obs_indices(end_pts_ind);
                        
                        double slope = (Y_i_new(h,end_pts(1)) - Y_i_new(h,end_pts(0))) / (end_pts(1) - end_pts(0));
                        Y_i_new(h,ind_j) = slope * (ind_j - end_pts(0)) + Y_i_new(h,end_pts(0));// + R::rnorm(0, std_R(h));
                    }
                }
            }
        }
        
        // STEP 2: Use HR and MAP to find "most likely state seq" --------------
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

                    init_vals(jj) = log_like_i_sub(k, Y_i_new, z_i, b_i, alpha_i,
                                                   omega_i, X_beta, D_alpha,
                                                   D_omega, X_beta_1, D_alpha_1,
                                                   D_omega_1, vec_beta, vec_A,
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


                        poss_state_like(w_ind) = log_like_i_sub(k, Y_i_new, z_i, b_i, alpha_i,
                                                                omega_i, X_beta, D_alpha,
                                                                D_omega, X_beta_1, D_alpha_1,
                                                                D_omega_1, vec_beta, vec_A,
                                                                R, zeta, P_init);

                        w_ind += 1;
                    }
                }

                b_i(k) = poss_next_state(poss_state_like.index_max());
            }
        }

        // Format the Dn_alpha list 
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
        newY.rows(sub_ind) = Y_i_new.t();
    }
    
    List Y_B_Dn = List::create(newY, B_return, Dn_return);
    return Y_B_Dn;
}

// MISC ------------------------------------------------------------------------
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
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::vec vec_A_scale = exp(vec_A_total) / (1 + exp(vec_A_total));
    // arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);

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
                // arma::vec vec_A = A_all_state.col(b_i(k) - 1);
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

                // arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                // arma::mat A_p1 = arma::diagmat(vec_A_p1);
                arma::mat A_p1 = arma::diagmat(vec_A); 

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

                // arma::vec vec_A = A_all_state.col(b_i(k) - 1);
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

                // arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
                // arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                // arma::mat A_k = arma::diagmat(vec_A_k);
                // arma::mat A_p1 = arma::diagmat(vec_A_p1);
                arma::mat A_k = arma::diagmat(vec_A);
                arma::mat A_p1 = arma::diagmat(vec_A);

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

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
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
                double log_target_prev = log_f_i_cpp(i, ii, t_pts, par, par_index,A_temp,
                                                     B_temp,Y_temp,z_temp,Dn_temp,Xn_temp,
                                                     Dn_omega_temp, W_temp);

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

                double log_target = log_f_i_cpp(i,ii,t_pts,par,par_index,A_temp,
                                                pr_B, pr_Y,z_temp,pr_Dn,Xn_temp,
                                                Dn_omega_temp, W_temp);

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

// [[Rcpp::export]]
void test_fnc(int states_per_step) {

    arma::vec s_i = {4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,2,2,2,2,2,2,2,2,2,
                     2,2,2,2,4,4,4,4,5,5,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                     2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,
                     3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                     2,2,2,2,4,4,4,4,4,4,4,4,5,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                     3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                     3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4};
    arma::vec bleed_ind_i = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
    double first_bleed_ind = arma::as_scalar(bleed_ind_ind);
    arma::vec bleed_ind_checks = {0, first_bleed_ind};
    if(first_bleed_ind > 0) {
        bleed_ind_checks(0) = first_bleed_ind - 1;
    }
    
    int N = 5;
    
    Rcpp::Rcout << "Case (c) Full" << std::endl;
    for(int w=0; w < N; w++) {
        Rcpp::Rcout << "() -> () -> " << w+1 << std::endl;
        Rcpp::Rcout << Omega_List_GLOBAL_multi(0)(w).n_rows << " combos" << std::endl;
        Rcpp::Rcout << Omega_List_GLOBAL_multi(0)(w) << std::endl;
    }
    
    Rcpp::Rcout << "Case (b) Full" << std::endl;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            Rcpp::Rcout << i+1 << "-->" << j+1 << std::endl;
            Rcpp::Rcout << Omega_List_GLOBAL_multi(1)(i, j).n_rows << " combos" << std::endl;
            Rcpp::Rcout << Omega_List_GLOBAL_multi(1)(i, j) << std::endl;
        }
    }
    
    Rcpp::Rcout << "Case (a) Full" << std::endl;
    for(int w=0; w < N; w++) {
        Rcpp::Rcout << w + 1 << " -> () -> ()" << std::endl;
        Rcpp::Rcout << Omega_List_GLOBAL_multi(2)(w).n_rows << " combos" << std::endl;
        Rcpp::Rcout << Omega_List_GLOBAL_multi(2)(w) << std::endl;
    }
    
} 
