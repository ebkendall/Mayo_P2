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

// [[Rcpp::export]]
void initialize_cpp(arma::imat a_mat, arma::imat a_mat_sub) {
    adj_mat_GLOBAL = a_mat;
    adj_mat_sub_GLOBAL = a_mat_sub;
    
    Omega_List_GLOBAL_multi = get_Omega_list(adj_mat_GLOBAL, 2);    
    Omega_List_GLOBAL_sub_multi = get_Omega_list(adj_mat_sub_GLOBAL, 2);
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

double log_like(const arma::vec &EIDs, const arma::vec &par, 
                const arma::field<arma::uvec> &par_index, 
                const arma::field <arma::vec> &A, const arma::field <arma::vec> &B, 
                const arma::field<arma::field<arma::mat>> &Dn, const arma::mat &Y, 
                int n_cores, const arma::mat &gamma_1) {

    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::vec alpha_tilde = par.elem(par_index(0) - 1);
    
    arma::vec vec_upsilon = par.elem(par_index(1) - 1);
    arma::mat upsilon_alpha = arma::reshape(vec_upsilon, alpha_tilde.n_elem, 
                                            alpha_tilde.n_elem);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(3) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    // arma::vec vec_R = par.elem(par_index(3) - 1);
    // arma::mat R = arma::reshape(vec_R, 4, 4);
    
    arma::vec vec_A_total = par.elem(par_index(2) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    
    arma::vec vec_zeta_content = par.elem(par_index(4) - 1);
    arma::vec qz = exp(vec_zeta_content);
    
    arma::vec vec_init_content = par.elem(par_index(5) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    
    arma::mat gamma_var = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
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
    
    arma::mat G_sqrt_t = arma::reshape(par.elem(par_index(6) - 1), 4, 4);
    arma::mat G_sqrt = G_sqrt_t.t() * G_sqrt_t;
    arma::mat G = G_sqrt * G_sqrt;
    
    arma::vec eids = Y.col(0);
    // -------------------------------------------------------------------------
    
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::field<arma::mat> D_alpha = Dn(ii);
        
        double like_comp_i = 0;
        
        for(int jj = 0; jj < n_i; jj++) {
            
            if(jj == 0) {
                
                arma::vec gamma_1_mu = Y_i.col(jj);
                arma::mat gamma_1_v = gamma_var + G;
                // arma::mat gamma_1_v = G;
                
                arma::vec log_g_pdf = dmvnorm(gamma_1.row(ii), gamma_1_mu, gamma_1_v, true);
                
                arma::vec log_a_pdf = dmvnorm(vec_alpha_i.t(), alpha_tilde, upsilon_alpha, true);
                
                like_comp_i = like_comp_i + arma::as_scalar(log_g_pdf) + arma::as_scalar(log_a_pdf);
                
            } else {
                
                arma::mat Q = { {    1,  qz(0),     0,  qz(1),     0},
                                {    0,      1, qz(2),  qz(3),     0},
                                {qz(4),  qz(5),     1,  qz(6),     0},
                                {    0,  qz(7),     0,      1, qz(8)},
                                {qz(9), qz(10),     0, qz(11),     1}};
                
                arma::vec q_row_sums = arma::sum(Q, 1);
                arma::mat P_i = Q.each_col() / q_row_sums;
                
                arma::vec y_j = Y_i.col(jj);
                arma::vec y_j_1 = Y_i.col(jj-1); 
                arma::vec nu_j = gamma_1.row(ii).t() + D_alpha(jj) * vec_alpha_i;
                arma::vec nu_j_1 = gamma_1.row(ii).t() + D_alpha(jj-1) * vec_alpha_i;
                
                arma::vec y_mean = nu_j + A_1 * (y_j_1 - nu_j_1);
                
                arma::vec log_y_pdf = dmvnorm(y_j.t(), y_mean, R, true);
                
                like_comp_i = like_comp_i + arma::as_scalar(log_y_pdf);
                
                if(jj == 1) {
                    like_comp_i = like_comp_i + log(P_init(b_i(jj) - 1));
                } else {
                    like_comp_i = like_comp_i + log(P_i(b_i(jj-1) - 1, b_i(jj) - 1));
                }
            }
        }
        
        in_vals(ii) = like_comp_i;
    }
    
    double in_value = arma::accu(in_vals);
    
    return in_value;
}

// [[Rcpp::export]]
double log_post(const arma::vec &EIDs, const arma::vec &par, 
                const arma::field<arma::uvec> &par_index, 
                const arma::field <arma::vec> &A, 
                const arma::field <arma::vec> &B, 
                const arma::field<arma::field<arma::mat>> &Dn, 
                const arma::mat &Y, int n_cores,
                const arma::mat &gamma_1) {

    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Likelihood --------------------------------------------------------------
    double value = log_like(EIDs, par, par_index, A, B, Dn, Y, n_cores, gamma_1);
    
    // Autocorrelation prior ---------------------------------------------------
    arma::vec vec_A_content = par.elem(par_index(2) - 1);
    arma::vec vec_A_mean(vec_A_content.n_elem, arma::fill::zeros);
    arma::vec scalar_A(vec_A_content.n_elem, arma::fill::ones);
    scalar_A = 100 * scalar_A;
    arma::mat A_var = arma::diagmat(scalar_A);
    
    arma::vec prior_A = dmvnorm(vec_A_content.t(), vec_A_mean, A_var, true);
    double prior_A_val = arma::as_scalar(prior_A);
    
    // Error-variance prior ----------------------------------------------------
    arma::vec vec_R_content = par.elem(par_index(3) - 1);
    arma::vec vec_R_mean(vec_R_content.n_elem, arma::fill::zeros);
    arma::vec scalar_R(vec_R_content.n_elem, arma::fill::ones);
    scalar_R = 100 * scalar_R;
    arma::mat R_var = arma::diagmat(scalar_R);
    
    arma::vec prior_R = dmvnorm(vec_R_content.t(), vec_R_mean, R_var, true);
    double prior_R_val = arma::as_scalar(prior_R);
    // arma::vec vec_R_content = par.elem(par_index(3) - 1);
    // arma::mat R = arma::reshape(vec_R_content, 4, 4);
    // 
    // int nu_R = 8;
    // arma::vec scalar_vec_R = {9, 9, 9, 9};
    // scalar_vec_R = (nu_R - 4 - 1) * scalar_vec_R;
    // arma::mat Psi_R = arma::diagmat(scalar_vec_R);
    // double prior_R_val = diwish(R, nu_R, Psi_R, true);
    
    // Zeta prior --------------------------------------------------------------
    arma::vec vec_zeta_content = par.elem(par_index(4) - 1);
    arma::vec vec_zeta_mean = {-3.7405, -4.2152, -2.6473, -2.1475, 
                               -3.4459, -2.9404, -3.2151, -3.1778, 
                               -2.0523, -3.4459, -3.2404, -3.2151};
    arma::vec scalar_zeta(vec_zeta_mean.n_elem, arma::fill::ones);
    scalar_zeta = 100 * scalar_zeta;
    arma::mat zeta_var = arma::diagmat(scalar_zeta);
    
    arma::vec prior_zeta = dmvnorm(vec_zeta_content.t(), vec_zeta_mean, zeta_var, true);
    double prior_zeta_val = arma::as_scalar(prior_zeta);
    
    // Initial Probability prior -----------------------------------------------
    arma::vec vec_init_content = par.elem(par_index(5) - 1);
    arma::vec vec_init_mean = {0, 0, 0, 0};
    arma::vec scalar_init(vec_init_content.n_elem, arma::fill::ones);
    scalar_init = 100 * scalar_init;
    arma::mat init_var = arma::diagmat(scalar_init);

    arma::vec prior_init = dmvnorm(vec_init_content.t(), vec_init_mean, init_var, true);
    double prior_init_val = arma::as_scalar(prior_init);
    
    // Gamma RE variance prior -------------------------------------------------
    arma::vec vec_G_content = par.elem(par_index(6) - 1);
    arma::vec vec_G_mean(vec_G_content.n_elem, arma::fill::zeros);
    arma::vec scalar_G(vec_G_content.n_elem, arma::fill::ones);
    scalar_G = 100 * scalar_G;
    arma::mat G_var = arma::diagmat(scalar_G);
    
    arma::vec prior_G = dmvnorm(vec_G_content.t(), vec_G_mean, G_var, true);
    double prior_G_val = arma::as_scalar(prior_G);
    
    // Full log-posterior ------------------------------------------------------
    value = value + prior_A_val + prior_R_val + prior_zeta_val + prior_init_val + prior_G_val;
    
    return value;
}

// [[Rcpp::export]]
arma::mat update_gamma_i(const arma::vec &EIDs, const arma::vec &par,
                         const arma::field<arma::uvec> &par_index,
                         const arma::field <arma::vec> &A, 
                         const arma::field <arma::vec> &B,
                         const arma::field<arma::field<arma::mat>> &Dn,
                         const arma::mat &Y, int n_cores){

    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID

    // Initializing parameters -------------------------------------------------
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(3) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    // arma::vec vec_R = par.elem(par_index(3) - 1);
    // arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat inv_R = arma::inv_sympd(R);

    arma::vec vec_A_total = par.elem(par_index(2) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    arma::mat I(A_1.n_rows, A_1.n_cols, arma::fill::eye);
    
    arma::mat gamma_var = {{R(0,0) / (1 - vec_A(0) * vec_A(0)),
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
    
    arma::mat G_sqrt_t = arma::reshape(par.elem(par_index(6) - 1), 4, 4);
    arma::mat G_sqrt = G_sqrt_t.t() * G_sqrt_t;
    arma::mat G = G_sqrt * G_sqrt;
    arma::mat g_var_combo = gamma_var + G;
    arma::mat g_var_combo_inv = arma::inv_sympd(g_var_combo);

    arma::vec eids = Y.col(0);
    // -------------------------------------------------------------------------

    arma::mat gamma_1(EIDs.n_elem, 4, arma::fill::zeros);

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();

        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::field<arma::mat> D_alpha = Dn(ii);

        arma::mat diff_g = I - A_1;
        arma::mat inv_W_i = g_var_combo_inv + (n_i - 1) * (diff_g.t() * inv_R * diff_g);
        arma::mat W_i = arma::inv_sympd(inv_W_i);
        
        arma::vec V_i = g_var_combo_inv * Y_i.col(0);

        for(int k = 1; k < n_i; k++) {

            arma::mat D_k = D_alpha(k);
            arma::mat D_k_1 = D_alpha(k-1);

            arma::mat diff_d = D_k - A_1 * D_k_1;
            arma::mat diff_y = Y_i.col(k) - A_1 * Y_i.col(k-1);
            
            arma::vec m_k = diff_y - diff_d * vec_alpha_i;

            V_i = V_i + diff_g.t() * inv_R * m_k;
        }

        arma::vec mu_gamma_i = W_i * V_i;

        arma::vec gamma_i = arma::mvnrnd(mu_gamma_i, W_i, 1);

        gamma_1.row(ii) = gamma_i.t();
    }

    return gamma_1;

}

// [[Rcpp::export]]
arma::field<arma::vec> update_alpha_i(const arma::vec &EIDs, const arma::vec &par, 
                                      const arma::field<arma::uvec> &par_index, 
                                      const arma::field <arma::vec> &B, 
                                      const arma::field<arma::field<arma::mat>> &Dn, 
                                      const arma::mat &Y, int n_cores,
                                      const arma::mat &gamma_1){

    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::vec alpha_tilde = par.elem(par_index(0) - 1);
    
    arma::vec vec_upsilon = par.elem(par_index(1) - 1);
    arma::mat upsilon_alpha = arma::reshape(vec_upsilon, alpha_tilde.n_elem, 
                                            alpha_tilde.n_elem);
    arma::mat inv_ups_alpha = arma::inv_sympd(upsilon_alpha);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(3) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    // arma::vec vec_R = par.elem(par_index(3) - 1);
    // arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat inv_R = arma::inv_sympd(R);
    
    arma::vec vec_A_total = par.elem(par_index(2) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    arma::mat I(A_1.n_rows, A_1.n_cols, arma::fill::eye);
    
    arma::vec eids = Y.col(0);
    // -------------------------------------------------------------------------
    
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
        
        arma::vec b_i = B(ii);
        arma::field<arma::mat> D_alpha = Dn(ii);
        
        arma::mat inv_W_i = inv_ups_alpha;
        arma::vec V_i     = inv_ups_alpha * alpha_tilde;
        
        for(int k = 1; k < n_i; k++) {
            
            arma::mat D_k = D_alpha(k);
            arma::mat D_k_1 = D_alpha(k-1);
            
            arma::mat diff_d = D_k - A_1 * D_k_1;
            arma::mat diff_y = Y_i.col(k) - A_1 * Y_i.col(k-1);
            arma::mat diff_g = I - A_1;
            
            arma::vec m_k = diff_y - diff_g * gamma_1.row(ii).t();
            
            inv_W_i = inv_W_i + diff_d.t() * inv_R * diff_d;
            
            V_i = V_i + diff_d.t() * inv_R * m_k;
        }
        
        arma::mat W_i = arma::inv_sympd(inv_W_i);
        
        arma::vec mu_alpha_i = W_i * V_i;
        
        arma::vec alpha_i = arma::mvnrnd(mu_alpha_i, W_i, 1);
        
        A(ii) = alpha_i;
    }
    
    return A;
  
}

// [[Rcpp::export]]
arma::vec update_alpha_tilde(const arma::vec EIDs, arma::vec par, 
                             const arma::field<arma::uvec> par_index,
                             const arma::field <arma::vec> A, const arma::mat Y){

    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::vec alpha_0(par_index(0).n_elem, arma::fill::zeros);
    arma::vec inv_sigma_alpha_diag(alpha_0.n_elem, arma::fill::ones);
    inv_sigma_alpha_diag = 0.01 * inv_sigma_alpha_diag;
    arma::mat inv_sigma_alpha = arma::diagmat(inv_sigma_alpha_diag);
    
    arma::vec vec_upsilon = par.elem(par_index(1) - 1);
    arma::mat upsilon_alpha = arma::reshape(vec_upsilon, par_index(0).n_elem, 
                                            par_index(0).n_elem);
    arma::mat inv_ups_alpha = arma::inv_sympd(upsilon_alpha);
    // -------------------------------------------------------------------------
    
    int N_id = EIDs.n_elem;
    
    arma::mat inv_W = inv_sigma_alpha + N_id * inv_ups_alpha;
    arma::mat W = arma::inv_sympd(inv_W);
    
    arma::vec total_alpha = A(0);
    for (int ii = 1; ii < N_id; ii++){
      total_alpha = total_alpha + A(ii);
    }
    
    arma::vec V = inv_sigma_alpha * alpha_0 + inv_ups_alpha * total_alpha;
    
    arma::vec mu = W * V;
    
    arma::vec alpha_tilde_temp = arma::mvnrnd(mu, W, 1);
    par.elem(par_index(0) - 1) = alpha_tilde_temp;
    
    return par;
}

// [[Rcpp::export]]
arma::vec update_beta_upsilon(const arma::vec &EIDs, arma::vec par, 
                              const arma::field<arma::uvec> &par_index,
                              const arma::field <arma::vec> &A,int n_cores) {
    
    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_upsilon = par.elem(par_index(1) - 1);
    int nu_ups = 20;
    arma::vec scalar_mult = {4, 4, 1, 1, 
                             4, 4, 1, 1, 
                             4, 4, 1, 1, 
                             4, 4, 1, 1};
    scalar_mult = (nu_ups - 16 - 1) * scalar_mult;
    arma::mat psi_ups = arma::diagmat(scalar_mult);
    
    arma::vec vec_alpha_tilde = par.elem(par_index(0) - 1);
    
    int N_id = EIDs.n_elem;
    // -------------------------------------------------------------------------
    
    arma::field<arma::mat> ups_list(N_id);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < N_id; ii++) {
        
        arma::vec vec_alpha_i = A(ii);
        
        arma::mat hold2 = vec_alpha_i - vec_alpha_tilde;
        arma::mat in_Upsilon_cov = hold2 * hold2.t();
        
        ups_list(ii) = in_Upsilon_cov;
    }
    
    arma::mat sum_ups_cov = psi_ups;
    
    for(int ii = 0; ii < N_id; ii++){
        sum_ups_cov += ups_list(ii);
    }
    
    par.elem(par_index(1) - 1) = arma::vectorise(riwish(nu_ups + N_id, sum_ups_cov));
    
    return par;
}

// [[Rcpp::export]]
Rcpp::List proposal_R(const int nu_R, const arma::mat psi_R, arma::mat curr_R,
                      const arma::vec &EIDs, const arma::vec &par,
                      const arma::field<arma::uvec> &par_index,
                      const arma::field <arma::vec> &A, const arma::field <arma::vec> &B,
                      const arma::field<arma::field<arma::mat>> &Dn,
                      const arma::mat &Y, int n_cores,
                      const arma::mat &gamma_1) {
    
    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_A_total = par.elem(par_index(2) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    
    arma::mat curr_R_sqrt = arma::sqrtmat_sympd(curr_R);
    
    arma::vec qz = exp(par.elem(par_index(4) - 1));
    
    arma::vec vec_init_content = par.elem(par_index(5) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    
    // arma::mat gamma_var = {{curr_R(0,0) / (1 - vec_A(0) * vec_A(0)), 
    //                         curr_R(0,1) / (1 - vec_A(0) * vec_A(1)), 
    //                         curr_R(0,2) / (1 - vec_A(0) * vec_A(2)), 
    //                         curr_R(0,3) / (1 - vec_A(0) * vec_A(3))},
    //                        {curr_R(1,0) / (1 - vec_A(1) * vec_A(0)),
    //                         curr_R(1,1) / (1 - vec_A(1) * vec_A(1)),
    //                         curr_R(1,2) / (1 - vec_A(1) * vec_A(2)),
    //                         curr_R(1,3) / (1 - vec_A(1) * vec_A(3))},
    //                        {curr_R(2,0) / (1 - vec_A(2) * vec_A(0)),
    //                         curr_R(2,1) / (1 - vec_A(2) * vec_A(1)),
    //                         curr_R(2,2) / (1 - vec_A(2) * vec_A(2)),
    //                         curr_R(2,3) / (1 - vec_A(2) * vec_A(3))},
    //                        {curr_R(3,0) / (1 - vec_A(3) * vec_A(0)),
    //                         curr_R(3,1) / (1 - vec_A(3) * vec_A(1)),
    //                         curr_R(3,2) / (1 - vec_A(3) * vec_A(2)),
    //                         curr_R(3,3) / (1 - vec_A(3) * vec_A(3))}};
    // 
    // arma::mat g_diag = arma::diagmat(exp(par.elem(par_index(6) - 1)));
    // 
    // arma::mat gamma_var_combo = gamma_var + g_diag;
    
    arma::mat G_sqrt_t = arma::reshape(par.elem(par_index(6) - 1), 4, 4);
    arma::mat G_sqrt = G_sqrt_t.t() * G_sqrt_t;
    arma::mat G = G_sqrt * G_sqrt;
    arma::mat gamma_var_combo = G;
    
    arma::mat inv_gamma_var = arma::inv_sympd(gamma_var_combo);
    arma::mat inv_gamma_sqrt = arma::sqrtmat_sympd(inv_gamma_var);
    
    arma::vec eids = Y.col(0);
    // -------------------------------------------------------------------------
    
    arma::field<arma::mat> psi_q_list(EIDs.n_elem);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        
        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::field<arma::mat> Dn_i = Dn(ii);
        
        arma::mat psi_q_i(4, 4, arma::fill::zeros);
        
        for(int k = 0; k < Y_i.n_cols; k++) {
            if(k == 0) {
                
                arma::vec diff_k = curr_R_sqrt * inv_gamma_sqrt * 
                                    (gamma_1.row(ii).t() - Y_i.col(k));
                
                psi_q_i = psi_q_i + diff_k * diff_k.t();

            } else {
                
                arma::vec hold_k_1 = Y_i.col(k-1) - gamma_1.row(ii).t() - Dn_i(k-1) * vec_alpha_i;
                arma::vec mu_k = gamma_1.row(ii).t() + Dn_i(k) * vec_alpha_i + A_1 * hold_k_1;
                
                arma::vec diff_k = Y_i.col(k) - mu_k;
                
                psi_q_i = psi_q_i + diff_k * diff_k.t();
            }
        }
        
        psi_q_list(ii) = psi_q_i;
    }
    
    arma::mat psi_q = psi_R;
    for(int ii = 0; ii < EIDs.n_elem; ii++) {
        psi_q += psi_q_list(ii);
    }

    int nu_q = Y.n_rows + nu_R;

    List nu_psi_R = List::create(psi_q, nu_q);
    return nu_psi_R;
}

// [[Rcpp::export]]
arma::field<arma::field<arma::mat>> initialize_Dn(const arma::vec EIDs,
                                                  arma::field <arma::vec> &B) {
    
    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);

    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        // Subject-specific information ----------------------------------------
        arma::vec b_i = B(ii);
        
        int n_i = b_i.n_elem;
        
        arma::field<arma::mat> Dn_temp(n_i);
        
        arma::vec twos(n_i, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(b_i == 2)) += 1;
        threes.elem(arma::find(b_i == 3)) += 1;
        fours.elem(arma::find(b_i == 4)) += 1;
        fives.elem(arma::find(b_i == 5)) += 1;
        twos(0) = 0; threes(0) = 0; fours(0) = 0; fives(0) = 0;
        
        arma::mat bigB = arma::join_rows(arma::cumsum(twos), arma::cumsum(threes));
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            // Guarantee Dn(0) = 0
            if(jj == 0) {
                arma::mat zero_jj(1, bigB.n_cols, arma::fill::zeros);
                Dn_temp(jj) = arma::kron(I, zero_jj);
            } else {
                Dn_temp(jj) = arma::kron(I, bigB.row(jj));    
            }
        }
        
        Dn_return(ii) = Dn_temp;
    }
    
    return Dn_return;
}

arma::vec p_2_sampler(int n_i, arma::mat &Y_i, arma::imat adj_mat_i, 
                      arma::vec &b_i, arma::mat alpha_i, arma::mat A_1, 
                      arma::mat R, arma::mat P_i, arma::vec P_init) {

    // Because we ignore the first time point, the length of the state sequence
    // we care about is 1 less than the original.
    for(int k = 1; k < n_i - 1; k++) {

        arma::vec s_i = b_i;

        arma::mat omega_set = get_omega_list(k-1, n_i - 1, b_i.subvec(1, n_i-1), 2, false);
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
            arma::vec threes_s = twos_s;
            arma::vec fours_s = twos_s;
            arma::vec fives_s = twos_s;
            twos_s.elem(arma::find(s_i == 2)) += 1;
            threes_s.elem(arma::find(s_i == 3)) += 1;
            fours_s.elem(arma::find(s_i == 4)) += 1;
            fives_s.elem(arma::find(s_i == 5)) += 1;

            arma::vec twos_b(b_i.n_elem, arma::fill::zeros);
            arma::vec threes_b = twos_b;
            arma::vec fours_b = twos_b;
            arma::vec fives_b = twos_b;
            twos_b.elem(arma::find(b_i == 2)) += 1;
            threes_b.elem(arma::find(b_i == 3)) += 1;
            fours_b.elem(arma::find(b_i == 4)) += 1;
            fives_b.elem(arma::find(b_i == 5)) += 1;

            for(int t = k; t < n_i; t++) {

                // INDEXING STARTS AT 1 (NOT 0)
                // Computations for mean of candidate (s_i) ------------
                arma::vec s_2 = twos_s.subvec(1, t);
                arma::vec s_3 = threes_s.subvec(1, t);
                arma::vec s_4 = fours_s.subvec(1, t);
                arma::vec s_5 = fives_s.subvec(1, t);

                arma::vec nu_s = alpha_i.row(0).t() + 
                    arma::accu(s_2) * alpha_i.row(1).t() + 
                    arma::accu(s_3) * alpha_i.row(2).t() + 
                    arma::accu(s_4) * alpha_i.row(3).t() + 
                    arma::accu(s_5) * alpha_i.row(4).t();
                
                arma::vec nu_s_1;
                if(t == 1) {
                    nu_s_1 = alpha_i.row(0).t();
                } else {
                    arma::vec s_2_1 = twos_s.subvec(1, t-1);
                    arma::vec s_3_1 = threes_s.subvec(1, t-1);
                    arma::vec s_4_1 = fours_s.subvec(1, t-1);
                    arma::vec s_5_1 = fives_s.subvec(1, t-1);
                    
                    nu_s_1 = alpha_i.row(0).t() + 
                        arma::accu(s_2_1) * alpha_i.row(1).t() + 
                        arma::accu(s_3_1) * alpha_i.row(2).t() + 
                        arma::accu(s_4_1) * alpha_i.row(3).t() + 
                        arma::accu(s_5_1) * alpha_i.row(4).t();
                }
                
                arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);

                // Computations for mean of current (b_i) --------------
                arma::vec b_2 = twos_b.subvec(1, t);
                arma::vec b_3 = threes_b.subvec(1, t);
                arma::vec b_4 = fours_b.subvec(1, t);
                arma::vec b_5 = fives_b.subvec(1, t);
                
                arma::vec nu_b = alpha_i.row(0).t() + 
                    arma::accu(b_2) * alpha_i.row(1).t() + 
                    arma::accu(b_3) * alpha_i.row(2).t() + 
                    arma::accu(b_4) * alpha_i.row(3).t() + 
                    arma::accu(b_5) * alpha_i.row(4).t();
                
                arma::vec nu_b_1;
                if(t == 1) {
                    nu_b_1 = alpha_i.row(0).t();
                } else {
                    arma::vec b_2_1 = twos_b.subvec(1, t-1);
                    arma::vec b_3_1 = threes_b.subvec(1, t-1);
                    arma::vec b_4_1 = fours_b.subvec(1, t-1);
                    arma::vec b_5_1 = fives_b.subvec(1, t-1);
            
                    nu_b_1 = alpha_i.row(0).t() + 
                        arma::accu(b_2_1) * alpha_i.row(1).t() + 
                        arma::accu(b_3_1) * alpha_i.row(2).t() + 
                        arma::accu(b_4_1) * alpha_i.row(3).t() + 
                        arma::accu(b_5_1) * alpha_i.row(4).t();
                }
                
                arma::vec mean_b = nu_b + A_1 * (Y_i.col(t-1) - nu_b_1);

                arma::vec log_y_pdf_s = dmvnorm(Y_i.col(t).t(), mean_s, R, true);
                arma::vec log_y_pdf_b = dmvnorm(Y_i.col(t).t(), mean_b, R, true);

                if(t <= k + 2) {
                    if(t == 1) {
                        log_like_s = log_like_s + log(P_init(s_i(t) - 1));
                        log_like_b = log_like_b + log(P_init(b_i(t) - 1));
                    } else {
                        log_like_s = log_like_s + log(P_i(s_i(t-1)-1, s_i(t)-1));
                        log_like_b = log_like_b + log(P_i(b_i(t-1)-1, b_i(t)-1));
                    }
                }

                log_like_s = log_like_s + arma::as_scalar(log_y_pdf_s);
                log_like_b = log_like_b + arma::as_scalar(log_y_pdf_b);
            }

            double diff_check = log_like_s - log_like_b;

            double min_log = log(arma::randu(arma::distr_param(0,1)));
            if(diff_check > min_log){b_i = s_i;}
        }
    }

    return b_i;
}

arma::vec p_full_sampler(int n_i, arma::mat &Y_i, arma::imat adj_mat_i,
                         arma::vec &b_i, arma::mat alpha_i, arma::mat A_1, 
                         arma::mat R, arma::mat P_i, arma::vec P_init) {

    // Because we ignore the first time point, the length of the state sequence
    // we care about is 1 less than the original.

    arma::vec all_like_vals_b(n_i - 1, arma::fill::zeros);
    arma::vec all_like_vals_s(n_i - 1, arma::fill::zeros);

    arma::vec s_i(n_i, arma::fill::zeros);

    // Propose a full state sequence -------------------------------------------
    for (int k = 1; k < n_i; k++) {

        arma::vec like_vals_s;
        arma::vec like_vals_b;
        arma::vec ss_ind;

        // INDEXING STARTS AT 1 (NOT 0)
        if(k == 1) {
            like_vals_s.zeros(adj_mat_i.n_cols);
            like_vals_b.zeros(adj_mat_i.n_cols);

            for(int m = 0; m < adj_mat_i.n_cols; m++) {

                // Computations for mean of candidate (s_i) --------------------
                arma::vec s_sub = s_i.subvec(1, k); // start = 1
                s_sub(k-1) = m+1;
                arma::vec twos_s(s_sub.n_elem, arma::fill::zeros);
                arma::vec threes_s(s_sub.n_elem, arma::fill::zeros);
                arma::vec fours_s(s_sub.n_elem, arma::fill::zeros);
                arma::vec fives_s(s_sub.n_elem, arma::fill::zeros);
                twos_s.elem(arma::find(s_sub == 2)) += 1;
                threes_s.elem(arma::find(s_sub == 3)) += 1;
                fours_s.elem(arma::find(s_sub == 4)) += 1;
                fives_s.elem(arma::find(s_sub == 5)) += 1;
                
                arma::vec nu_s = alpha_i.row(0).t() + 
                    arma::accu(twos_s)   * alpha_i.row(1).t() + 
                    arma::accu(threes_s) * alpha_i.row(2).t() + 
                    arma::accu(fours_s)  * alpha_i.row(3).t() + 
                    arma::accu(fives_s)  * alpha_i.row(4).t();
                arma::vec nu_s_1 = alpha_i.row(0).t();
                
                arma::vec mean_s = nu_s + A_1 * (Y_i.col(k-1) - nu_s_1);

                // Computations for mean of current (b_i) ----------------------
                arma::vec b_sub = b_i.subvec(1, k); // start = 1
                b_sub(k-1) = m+1;
                arma::vec twos_b(b_sub.n_elem, arma::fill::zeros);
                arma::vec threes_b(b_sub.n_elem, arma::fill::zeros);
                arma::vec fours_b(b_sub.n_elem, arma::fill::zeros);
                arma::vec fives_b(b_sub.n_elem, arma::fill::zeros);
                twos_b.elem(arma::find(b_sub == 2)) += 1;
                threes_b.elem(arma::find(b_sub == 3)) += 1;
                fours_b.elem(arma::find(b_sub == 4)) += 1;
                fives_b.elem(arma::find(b_sub == 5)) += 1;
                
                arma::vec nu_b = alpha_i.row(0).t() + 
                    arma::accu(twos_b)   * alpha_i.row(1).t() + 
                    arma::accu(threes_b) * alpha_i.row(2).t() + 
                    arma::accu(fours_b)  * alpha_i.row(3).t() + 
                    arma::accu(fives_b)  * alpha_i.row(4).t();
                arma::vec nu_b_1 = alpha_i.row(0).t();
                
                arma::vec mean_b = nu_b + A_1 * (Y_i.col(k-1) - nu_b_1);

                arma::vec log_y_pdf_s = dmvnorm(Y_i.col(k).t(), mean_s, R, true);
                arma::vec log_y_pdf_b = dmvnorm(Y_i.col(k).t(), mean_b, R, true);

                like_vals_s(m) = log(P_init(m)) + arma::as_scalar(log_y_pdf_s);
                like_vals_b(m) = log(P_init(m)) + arma::as_scalar(log_y_pdf_b);
            }

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

            for(int m = 0; m < adj_mat_i.n_cols; m++) {

                // Computations for mean of candidate (s_i) --------------------
                if(adj_mat_i(prev_state_s-1, m) != 0) {
                    
                    arma::vec s_sub = s_i.subvec(1, k); // start = 1
                    s_sub(k-1) = m+1;
                    arma::vec twos_s(s_sub.n_elem, arma::fill::zeros);
                    arma::vec threes_s(s_sub.n_elem, arma::fill::zeros);
                    arma::vec fours_s(s_sub.n_elem, arma::fill::zeros);
                    arma::vec fives_s(s_sub.n_elem, arma::fill::zeros);
                    twos_s.elem(arma::find(s_sub == 2)) += 1;
                    threes_s.elem(arma::find(s_sub == 3)) += 1;
                    fours_s.elem(arma::find(s_sub == 4)) += 1;
                    fives_s.elem(arma::find(s_sub == 5)) += 1;
                    
                    arma::vec nu_s = alpha_i.row(0).t() + 
                        arma::accu(twos_s)   * alpha_i.row(1).t() + 
                        arma::accu(threes_s) * alpha_i.row(2).t() + 
                        arma::accu(fours_s)  * alpha_i.row(3).t() + 
                        arma::accu(fives_s)  * alpha_i.row(4).t();
                    arma::vec nu_s_1 = alpha_i.row(0).t() + 
                        arma::accu(twos_s.subvec(0,k-2))   * alpha_i.row(1).t() + 
                        arma::accu(threes_s.subvec(0,k-2)) * alpha_i.row(2).t() + 
                        arma::accu(fours_s.subvec(0,k-2))  * alpha_i.row(3).t() + 
                        arma::accu(fives_s.subvec(0,k-2))  * alpha_i.row(4).t();
                    
                    arma::vec mean_s = nu_s + A_1 * (Y_i.col(k-1) - nu_s_1);

                    arma::vec log_y_pdf = dmvnorm(Y_i.col(k).t(), mean_s, R, true);

                    like_vals_s(s_ind) = log(P_i(prev_state_s-1, m)) + arma::as_scalar(log_y_pdf);
                    state_vals_s(s_ind) = m;
                    s_ind = s_ind + 1;
                }

                // Computations for mean of current (b_i) ----------------------
                if(adj_mat_i(prev_state_b-1, m) != 0) {

                    arma::vec b_sub = b_i.subvec(1, k); // start = 1
                    b_sub(k-1) = m+1;
                    arma::vec twos_b(b_sub.n_elem, arma::fill::zeros);
                    arma::vec threes_b(b_sub.n_elem, arma::fill::zeros);
                    arma::vec fours_b(b_sub.n_elem, arma::fill::zeros);
                    arma::vec fives_b(b_sub.n_elem, arma::fill::zeros);
                    twos_b.elem(arma::find(b_sub == 2)) += 1;
                    threes_b.elem(arma::find(b_sub == 3)) += 1;
                    fours_b.elem(arma::find(b_sub == 4)) += 1;
                    fives_b.elem(arma::find(b_sub == 5)) += 1;
                    
                    arma::vec nu_b = alpha_i.row(0).t() + 
                        arma::accu(twos_b)   * alpha_i.row(1).t() + 
                        arma::accu(threes_b) * alpha_i.row(2).t() + 
                        arma::accu(fours_b)  * alpha_i.row(3).t() + 
                        arma::accu(fives_b)  * alpha_i.row(4).t();
                    arma::vec nu_b_1 = alpha_i.row(0).t() + 
                        arma::accu(twos_b.subvec(0,k-2))   * alpha_i.row(1).t() + 
                        arma::accu(threes_b.subvec(0,k-2)) * alpha_i.row(2).t() + 
                        arma::accu(fours_b.subvec(0,k-2))  * alpha_i.row(3).t() + 
                        arma::accu(fives_b.subvec(0,k-2))  * alpha_i.row(4).t();
                    
                    arma::vec mean_b = nu_b + A_1 * (Y_i.col(k-1) - nu_b_1);

                    arma::vec log_y_pdf = dmvnorm(Y_i.col(k).t(), mean_b, R, true);

                    like_vals_b(b_ind) = log(P_i(prev_state_b-1, m)) + arma::as_scalar(log_y_pdf);
                    state_vals_b(b_ind) = m;
                    b_ind = b_ind + 1;
                }
            }

            ss_ind = state_vals_s;
        }

        all_like_vals_s(k-1) = arma::accu(exp(like_vals_s));
        all_like_vals_b(k-1) = arma::accu(exp(like_vals_b));

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

arma::vec p_flex_sampler(int n_i, arma::mat &Y_i, arma::imat adj_mat_i,
                         arma::vec &b_i, arma::mat alpha_i, arma::mat A_1,
                         arma::mat R, arma::mat P_i, arma::vec P_init, int sps) {

    // Because we ignore the first time point, the length of the state sequence
    // we care about is 1 less than the original.

    int k = 1;

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

            // INDEXING STARTS AT 1 (NOT 0)
            if(t == 1) {
                like_vals_s.zeros(adj_mat_i.n_cols);
                like_vals_b.zeros(adj_mat_i.n_cols);

                for(int m = 0; m < adj_mat_i.n_cols; m++) {

                    // Computations for mean of candidate (s_i) --------
                    arma::vec s_sub = s_i.subvec(1, t); // start = 1
                    s_sub(t-1) = m+1;
                    arma::vec twos_s(s_sub.n_elem, arma::fill::zeros);
                    arma::vec threes_s(s_sub.n_elem, arma::fill::zeros);
                    arma::vec fours_s(s_sub.n_elem, arma::fill::zeros);
                    arma::vec fives_s(s_sub.n_elem, arma::fill::zeros);
                    twos_s.elem(arma::find(s_sub == 2)) += 1;
                    threes_s.elem(arma::find(s_sub == 3)) += 1;
                    fours_s.elem(arma::find(s_sub == 4)) += 1;
                    fives_s.elem(arma::find(s_sub == 5)) += 1;
                    
                    arma::vec nu_s = alpha_i.row(0).t() + 
                        arma::accu(twos_s)   * alpha_i.row(1).t() + 
                        arma::accu(threes_s) * alpha_i.row(2).t() + 
                        arma::accu(fours_s)  * alpha_i.row(3).t() + 
                        arma::accu(fives_s)  * alpha_i.row(4).t();
                    arma::vec nu_s_1 = alpha_i.row(0).t();
                    
                    arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);

                    // Computations for mean of current (b_i) ----------
                    arma::vec b_sub = b_i.subvec(1, t); // start = 1
                    b_sub(t-1) = m+1;
                    arma::vec twos_b(b_sub.n_elem, arma::fill::zeros);
                    arma::vec threes_b(b_sub.n_elem, arma::fill::zeros);
                    arma::vec fours_b(b_sub.n_elem, arma::fill::zeros);
                    arma::vec fives_b(b_sub.n_elem, arma::fill::zeros);
                    twos_b.elem(arma::find(b_sub == 2)) += 1;
                    threes_b.elem(arma::find(b_sub == 3)) += 1;
                    fours_b.elem(arma::find(b_sub == 4)) += 1;
                    fives_b.elem(arma::find(b_sub == 5)) += 1;
                    
                    arma::vec nu_b = alpha_i.row(0).t() + 
                        arma::accu(twos_b)   * alpha_i.row(1).t() + 
                        arma::accu(threes_b) * alpha_i.row(2).t() + 
                        arma::accu(fours_b)  * alpha_i.row(3).t() + 
                        arma::accu(fives_b)  * alpha_i.row(4).t();
                    arma::vec nu_b_1 = alpha_i.row(0).t();
                    
                    arma::vec mean_b = nu_b + A_1 * (Y_i.col(t-1) - nu_b_1);
                    
                    arma::vec log_y_pdf_s = dmvnorm(Y_i.col(t).t(), mean_s, R, true);
                    arma::vec log_y_pdf_b = dmvnorm(Y_i.col(t).t(), mean_b, R, true);

                    like_vals_s(m) = log(P_init(m)) + arma::as_scalar(log_y_pdf_s);
                    like_vals_b(m) = log(P_init(m)) + arma::as_scalar(log_y_pdf_b);
                }

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

                for(int m = 0; m < adj_mat_i.n_cols; m++) {

                    // Computations for mean of candidate (s_i) --------
                    if(adj_mat_i(prev_state_s-1, m) != 0) {
                        
                        arma::vec s_sub = s_i.subvec(1, t); // start = 1
                        s_sub(t-1) = m+1;
                        arma::vec twos_s(s_sub.n_elem, arma::fill::zeros);
                        arma::vec threes_s(s_sub.n_elem, arma::fill::zeros);
                        arma::vec fours_s(s_sub.n_elem, arma::fill::zeros);
                        arma::vec fives_s(s_sub.n_elem, arma::fill::zeros);
                        twos_s.elem(arma::find(s_sub == 2)) += 1;
                        threes_s.elem(arma::find(s_sub == 3)) += 1;
                        fours_s.elem(arma::find(s_sub == 4)) += 1;
                        fives_s.elem(arma::find(s_sub == 5)) += 1;
                        
                        arma::vec nu_s = alpha_i.row(0).t() + 
                            arma::accu(twos_s)   * alpha_i.row(1).t() + 
                            arma::accu(threes_s) * alpha_i.row(2).t() + 
                            arma::accu(fours_s)  * alpha_i.row(3).t() + 
                            arma::accu(fives_s)  * alpha_i.row(4).t();
                        arma::vec nu_s_1 = alpha_i.row(0).t() + 
                            arma::accu(twos_s.subvec(0,t-2))   * alpha_i.row(1).t() + 
                            arma::accu(threes_s.subvec(0,t-2)) * alpha_i.row(2).t() + 
                            arma::accu(fours_s.subvec(0,t-2))  * alpha_i.row(3).t() + 
                            arma::accu(fives_s.subvec(0,t-2))  * alpha_i.row(4).t();
                        
                        arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);

                        arma::vec log_y_pdf = dmvnorm(Y_i.col(t).t(), mean_s, R, true);

                        like_vals_s(s_ind) = log(P_i(prev_state_s-1, m)) + arma::as_scalar(log_y_pdf);
                        state_vals_s(s_ind) = m;
                        s_ind = s_ind + 1;
                    }

                    // Computations for mean of current (b_i) ----------
                    if(adj_mat_i(prev_state_b-1, m) != 0) {
                        
                        arma::vec b_sub = b_i.subvec(1, t); // start = 1
                        b_sub(t-1) = m+1;
                        arma::vec twos_b(b_sub.n_elem, arma::fill::zeros);
                        arma::vec threes_b(b_sub.n_elem, arma::fill::zeros);
                        arma::vec fours_b(b_sub.n_elem, arma::fill::zeros);
                        arma::vec fives_b(b_sub.n_elem, arma::fill::zeros);
                        twos_b.elem(arma::find(b_sub == 2)) += 1;
                        threes_b.elem(arma::find(b_sub == 3)) += 1;
                        fours_b.elem(arma::find(b_sub == 4)) += 1;
                        fives_b.elem(arma::find(b_sub == 5)) += 1;
                        
                        arma::vec nu_b = alpha_i.row(0).t() + 
                            arma::accu(twos_b)   * alpha_i.row(1).t() + 
                            arma::accu(threes_b) * alpha_i.row(2).t() + 
                            arma::accu(fours_b)  * alpha_i.row(3).t() + 
                            arma::accu(fives_b)  * alpha_i.row(4).t();
                        arma::vec nu_b_1 = alpha_i.row(0).t() + 
                            arma::accu(twos_b.subvec(0,t-2))   * alpha_i.row(1).t() + 
                            arma::accu(threes_b.subvec(0,t-2)) * alpha_i.row(2).t() + 
                            arma::accu(fours_b.subvec(0,t-2))  * alpha_i.row(3).t() + 
                            arma::accu(fives_b.subvec(0,t-2))  * alpha_i.row(4).t();
                        
                        arma::vec mean_b = nu_b + A_1 * (Y_i.col(t-1) - nu_b_1);
                        
                        arma::vec log_y_pdf = dmvnorm(Y_i.col(t).t(), mean_b, R, true);

                        like_vals_b(b_ind) = log(P_i(prev_state_b-1, m)) + arma::as_scalar(log_y_pdf);
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
        arma::mat Omega_set_s = get_omega_list(t_max, n_i, s_i, 2, false);
        arma::mat Omega_set_b = get_omega_list(t_max, n_i, b_i, 2, false);

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
            arma::vec threes_s = twos_s;
            arma::vec fours_s = twos_s;
            arma::vec fives_s = twos_s;
            twos_s.elem(arma::find(s_i == 2)) += 1;
            threes_s.elem(arma::find(s_i == 3)) += 1;
            fours_s.elem(arma::find(s_i == 4)) += 1;
            fives_s.elem(arma::find(s_i == 5)) += 1;
            
            arma::vec twos_b(b_i.n_elem, arma::fill::zeros);
            arma::vec threes_b = twos_b;
            arma::vec fours_b = twos_b;
            arma::vec fives_b = twos_b;
            twos_b.elem(arma::find(b_i == 2)) += 1;
            threes_b.elem(arma::find(b_i == 3)) += 1;
            fours_b.elem(arma::find(b_i == 4)) += 1;
            fives_b.elem(arma::find(b_i == 5)) += 1;

            for(int t = t_max; t < n_i; t++) {
                
                // INDEXING STARTS AT 1 (NOT 0)
                // Computations for mean of candidate (s_i) ------------
                arma::vec s_2 = twos_s.subvec(1, t);
                arma::vec s_3 = threes_s.subvec(1, t);
                arma::vec s_4 = fours_s.subvec(1, t);
                arma::vec s_5 = fives_s.subvec(1, t);
                
                arma::vec s_2_1 = twos_s.subvec(1, t-1);
                arma::vec s_3_1 = threes_s.subvec(1, t-1);
                arma::vec s_4_1 = fours_s.subvec(1, t-1);
                arma::vec s_5_1 = fives_s.subvec(1, t-1);

                arma::vec nu_s = alpha_i.row(0).t() +
                    arma::accu(s_2) * alpha_i.row(1).t() +
                    arma::accu(s_3) * alpha_i.row(2).t() +
                    arma::accu(s_4) * alpha_i.row(3).t() +
                    arma::accu(s_5) * alpha_i.row(4).t();
                arma::vec nu_s_1 = alpha_i.row(0).t() +
                    arma::accu(s_2_1) * alpha_i.row(1).t() +
                    arma::accu(s_3_1) * alpha_i.row(2).t() +
                    arma::accu(s_4_1) * alpha_i.row(3).t() +
                    arma::accu(s_5_1) * alpha_i.row(4).t();
                
                arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);

                // Computations for mean of current (b_i) --------------
                arma::vec b_2 = twos_b.subvec(1, t);
                arma::vec b_3 = threes_b.subvec(1, t);
                arma::vec b_4 = fours_b.subvec(1, t);
                arma::vec b_5 = fives_b.subvec(1, t);
                
                arma::vec b_2_1 = twos_b.subvec(1, t-1);
                arma::vec b_3_1 = threes_b.subvec(1, t-1);
                arma::vec b_4_1 = fours_b.subvec(1, t-1);
                arma::vec b_5_1 = fives_b.subvec(1, t-1);

                arma::vec nu_b = alpha_i.row(0).t() +
                    arma::accu(b_2) * alpha_i.row(1).t() +
                    arma::accu(b_3) * alpha_i.row(2).t() +
                    arma::accu(b_4) * alpha_i.row(3).t() +
                    arma::accu(b_5) * alpha_i.row(4).t();
                arma::vec nu_b_1 = alpha_i.row(0).t() +
                    arma::accu(b_2_1) * alpha_i.row(1).t() +
                    arma::accu(b_3_1) * alpha_i.row(2).t() +
                    arma::accu(b_4_1) * alpha_i.row(3).t() +
                    arma::accu(b_5_1) * alpha_i.row(4).t();
                
                arma::vec mean_b = nu_b + A_1 * (Y_i.col(t-1) - nu_b_1);

                arma::vec log_y_pdf_s = dmvnorm(Y_i.col(t).t(), mean_s, R, true);
                arma::vec log_y_pdf_b = dmvnorm(Y_i.col(t).t(), mean_b, R, true);

                if(t < k + sps + 1) {
                    log_prob_diff = log_prob_diff +
                        log(P_i(s_i(t-1)-1, s_i(t)-1)) - log(P_i(b_i(t-1)-1, b_i(t)-1));
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
Rcpp::List state_sampler(const arma::vec EIDs, const arma::vec &par,
                         const arma::field<arma::uvec> &par_index,
                         const arma::field <arma::vec> A,
                         arma::field <arma::vec> &B,
                         const arma::mat &Y, int states_per_step,
                         const arma::mat &gamma_1, int n_cores) {
    
    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);

    // Initializing parameters -------------------------------------------------
    arma::vec vec_A_total = par.elem(par_index(2) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(3) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    // arma::vec vec_R = par.elem(par_index(3) - 1);
    // arma::mat R = arma::reshape(vec_R, 4, 4);
    
    arma::vec qz = exp(par.elem(par_index(4) - 1));
    
    arma::vec vec_init_content = par.elem(par_index(5) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    
    arma::vec eids = Y.col(0);
    // -------------------------------------------------------------------------

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::mat Q = { {    1,  qz(0),     0,  qz(1),     0},
                        {    0,      1, qz(2),  qz(3),     0},
                        {qz(4),  qz(5),     1,  qz(6),     0},
                        {    0,  qz(7),     0,      1, qz(8)},
                        {qz(9), qz(10),     0, qz(11),     1}};
        
        arma::vec q_row_sums = arma::sum(Q, 1);
        arma::mat P_i = Q.each_col() / q_row_sums;
        
        arma::vec b_i = B(ii);

        arma::vec vec_alpha_i_no_base = A(ii);
        arma::mat alpha_i_no_base = arma::reshape(vec_alpha_i_no_base, 4, 4);
        
        arma::mat alpha_i(5, 4, arma::fill::zeros);
        alpha_i.row(0) = gamma_1.row(ii);
        alpha_i.row(1) = alpha_i_no_base.row(0);
        alpha_i.row(2) = alpha_i_no_base.row(1);
        alpha_i.row(3) = alpha_i_no_base.row(2);
        alpha_i.row(4) = alpha_i_no_base.row(3);
        
        arma::imat adj_mat_i = adj_mat_GLOBAL;

        arma::vec new_b_i;
        
        if(states_per_step == 2) {
            // p = 2
            new_b_i = p_2_sampler(n_i, Y_i, adj_mat_i, b_i, alpha_i, A_1, R, 
                                  P_i, P_init);

        } else if(states_per_step >= n_i - 1) {
            // p >= n_i - 1 (first time point doesn't count)
            new_b_i = p_full_sampler(n_i, Y_i, adj_mat_i, b_i, alpha_i, A_1, R,
                                     P_i, P_init);

        } else {
            // 2 < p < n_i - 1
            new_b_i = p_flex_sampler(n_i, Y_i, adj_mat_i, b_i, alpha_i, A_1, R,
                                     P_i, P_init, states_per_step);
        }
        
        // Format the Dn_alpha list --------------------------------------------
        arma::field<arma::mat> Dn_temp(n_i);
        arma::vec twos(new_b_i.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(new_b_i == 2)) += 1;
        threes.elem(arma::find(new_b_i == 3)) += 1;
        fours.elem(arma::find(new_b_i == 4)) += 1;
        fives.elem(arma::find(new_b_i == 5)) += 1;
        twos(0) = 0; threes(0) = 0; fours(0) = 0; fives(0) = 0;
        
        arma::mat bigB = arma::join_rows(arma::cumsum(twos), arma::cumsum(threes));
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            // Guarantee Dn(0) = 0
            if(jj == 0) {
                arma::mat zero_jj(1, bigB.n_cols, arma::fill::zeros);
                Dn_temp(jj) = arma::kron(I, zero_jj);
            } else {
                Dn_temp(jj) = arma::kron(I, bigB.row(jj));    
            }
        }

        B_return(ii) = new_b_i;
        Dn_return(ii) = Dn_temp;
    }

    List B_Dn = List::create(B_return, Dn_return);
    return B_Dn;
}

// [[Rcpp::export]]
Rcpp::List mle_state_seq(const arma::vec EIDs, const arma::vec &par,
                         const arma::field<arma::uvec> &par_index,
                         const arma::field <arma::vec> &A, const arma::mat &Y, 
                         int n_cores) {
    
    // par_index: (0) alpha_tilde, (1) upsilon, (2) A, (3) R, (4) zeta, (5) init, (6) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);

    // Initializing parameters -------------------------------------------------
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(3) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    // arma::vec vec_R = par.elem(par_index(3) - 1);
    // arma::mat R = arma::reshape(vec_R, 4, 4);
    
    arma::vec vec_A_total = par.elem(par_index(2) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    
    arma::vec qz = exp(par.elem(par_index(4) - 1));
    
    arma::vec vec_init_content = par.elem(par_index(5) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    
    arma::mat gamma_var = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
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
    
    arma::vec eids = Y.col(0);
    // -------------------------------------------------------------------------

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::vec b_i(n_i, arma::fill::zeros);

        arma::vec vec_alpha_i_no_base = A(ii);
        arma::mat alpha_i_no_base = arma::reshape(vec_alpha_i_no_base, 4, 4);
        
        arma::mat alpha_i(5, 4, arma::fill::zeros);
        alpha_i.row(0) = Y_i.col(0).t();
        alpha_i.row(1) = alpha_i_no_base.row(0);
        alpha_i.row(2) = alpha_i_no_base.row(1);
        alpha_i.row(3) = alpha_i_no_base.row(2);
        alpha_i.row(4) = alpha_i_no_base.row(3);

        arma::imat adj_mat_i = adj_mat_GLOBAL;
        
        arma::mat Q = { {    1,  qz(0),     0,  qz(1),     0},
                        {    0,      1, qz(2),  qz(3),     0},
                        {qz(4),  qz(5),     1,  qz(6),     0},
                        {    0,  qz(7),     0,      1, qz(8)},
                        {qz(9), qz(10),     0, qz(11),     1}};
        
        arma::vec q_row_sums = arma::sum(Q, 1);
        arma::mat P_i = Q.each_col() / q_row_sums;

        // Looping through subject state space ---------------------------------
        // NOTE: start with k = 1 (not k = 0)!
        for(int k = 1; k < n_i; k++) {
            
            if(k == 1) {
                arma::vec init_vals(adj_mat_i.n_cols, arma::fill::zeros);
                
                // Consider all possible initial states
                for(int jj = 0; jj < init_vals.n_elem; jj++) {
                    
                    b_i(k) = jj + 1;

                    arma::vec b_sub = b_i.subvec(1, k); // start = 1

                    arma::vec twos(b_sub.n_elem, arma::fill::zeros);
                    arma::vec threes(b_sub.n_elem, arma::fill::zeros);
                    arma::vec fours(b_sub.n_elem, arma::fill::zeros);
                    arma::vec fives(b_sub.n_elem, arma::fill::zeros);
                    twos.elem(arma::find(b_sub == 2)) += 1;
                    threes.elem(arma::find(b_sub == 3)) += 1;
                    fours.elem(arma::find(b_sub == 4)) += 1;
                    fives.elem(arma::find(b_sub == 5)) += 1;
                    
                    arma::vec nu_k = alpha_i.row(0).t() + 
                        arma::accu(twos)   * alpha_i.row(1).t() + 
                        arma::accu(threes) * alpha_i.row(2).t() + 
                        arma::accu(fours)  * alpha_i.row(3).t() + 
                        arma::accu(fives)  * alpha_i.row(4).t();
                    
                    arma::vec nu_k_1 = alpha_i.row(0).t();
                    
                    arma::vec mean_b = nu_k + A_1 * (Y_i.col(k-1) - nu_k_1);
                    
                    arma::vec log_y_pdf = dmvnorm(Y_i.col(k).t(), mean_b, R, true);
                    init_vals(jj) = log(P_init(jj)) + arma::as_scalar(log_y_pdf);
                }

                b_i(k) = arma::index_max(init_vals) + 1;

            } else {
                
                // All other states -----------------
                int prev_state = b_i(k-1);
                arma::vec poss_next_state(arma::accu(adj_mat_i.row(prev_state-1)), arma::fill::zeros);
                arma::vec poss_state_like(arma::accu(adj_mat_i.row(prev_state-1)), arma::fill::zeros);
                
                int w_ind = 0;
                for(int jj = 0; jj < adj_mat_i.n_cols; jj++) {
                    if(adj_mat_i(prev_state-1, jj) != 0) {
                        
                        b_i(k) = jj + 1;
                        
                        arma::vec b_sub = b_i.subvec(1, k); // start = 1
                        
                        arma::vec twos(b_sub.n_elem, arma::fill::zeros);
                        arma::vec threes(b_sub.n_elem, arma::fill::zeros);
                        arma::vec fours(b_sub.n_elem, arma::fill::zeros);
                        arma::vec fives(b_sub.n_elem, arma::fill::zeros);
                        twos.elem(arma::find(b_sub == 2)) += 1;
                        threes.elem(arma::find(b_sub == 3)) += 1;
                        fours.elem(arma::find(b_sub == 4)) += 1;
                        fives.elem(arma::find(b_sub == 5)) += 1;
                        
                        arma::vec nu_k = alpha_i.row(0).t() + 
                            arma::accu(twos)   * alpha_i.row(1).t() + 
                            arma::accu(threes) * alpha_i.row(2).t() + 
                            arma::accu(fours)  * alpha_i.row(3).t() + 
                            arma::accu(fives)  * alpha_i.row(4).t();
                        
                        arma::vec nu_k_1 = alpha_i.row(0).t() + 
                            arma::accu(twos.subvec(0,k-2))   * alpha_i.row(1).t() + 
                            arma::accu(threes.subvec(0,k-2)) * alpha_i.row(2).t() + 
                            arma::accu(fours.subvec(0,k-2))  * alpha_i.row(3).t() + 
                            arma::accu(fives.subvec(0,k-2))  * alpha_i.row(4).t();
                        
                        arma::vec mean_b = nu_k + A_1 * (Y_i.col(k-1) - nu_k_1);
                        
                        arma::vec log_y_pdf = dmvnorm(Y_i.col(k).t(), mean_b, R, true);
                        
                        poss_next_state(w_ind) = jj + 1;
                        poss_state_like(w_ind) = log(P_i(prev_state-1, jj)) + arma::as_scalar(log_y_pdf);
                        
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
        twos(0) = 0; threes(0) = 0; fours(0) = 0; fives(0) = 0;
        
        arma::mat bigB = arma::join_rows(arma::cumsum(twos), arma::cumsum(threes));
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            // Guarantee Dn(0) = 0
            if(jj == 0) {
                arma::mat zero_jj(1, bigB.n_cols, arma::fill::zeros);
                Dn_temp(jj) = arma::kron(I, zero_jj);
            } else {
                Dn_temp(jj) = arma::kron(I, bigB.row(jj));    
            }
        }
        
        B_return(ii) = b_i;
        Dn_return(ii) = Dn_temp;
    }
    
    List B_Dn = List::create(B_return, Dn_return);
    return B_Dn;
}

// [[Rcpp::export]]
void test_fnc(int states_per_step) {
    
    arma::vec temp = {1,2,3,4,5};

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
    // arma::uvec bleed_ind_ind = arma::find(bleed_ind_i == 1);
    // double first_bleed_ind = arma::as_scalar(bleed_ind_ind);
    // arma::vec bleed_ind_checks = {0, first_bleed_ind};
    // if(first_bleed_ind > 0) {
    //     bleed_ind_checks(0) = first_bleed_ind - 1;
    // }
    // 
    // int N = 5;
    // 
    // Rcpp::Rcout << "Case (c) Full" << std::endl;
    // for(int w=0; w < N; w++) {
    //     Rcpp::Rcout << "() -> () -> " << w+1 << std::endl;
    //     Rcpp::Rcout << Omega_List_GLOBAL_multi(0)(w).n_rows << " combos" << std::endl;
    //     Rcpp::Rcout << Omega_List_GLOBAL_multi(0)(w) << std::endl;
    // }
    // 
    // Rcpp::Rcout << "Case (b) Full" << std::endl;
    // for(int i = 0; i < N; i++) {
    //     for(int j = 0; j < N; j++) {
    //         Rcpp::Rcout << i+1 << "-->" << j+1 << std::endl;
    //         Rcpp::Rcout << Omega_List_GLOBAL_multi(1)(i, j).n_rows << " combos" << std::endl;
    //         Rcpp::Rcout << Omega_List_GLOBAL_multi(1)(i, j) << std::endl;
    //     }
    // }
    // 
    // Rcpp::Rcout << "Case (a) Full" << std::endl;
    // for(int w=0; w < N; w++) {
    //     Rcpp::Rcout << w + 1 << " -> () -> ()" << std::endl;
    //     Rcpp::Rcout << Omega_List_GLOBAL_multi(2)(w).n_rows << " combos" << std::endl;
    //     Rcpp::Rcout << Omega_List_GLOBAL_multi(2)(w) << std::endl;
    // }
    
} 
