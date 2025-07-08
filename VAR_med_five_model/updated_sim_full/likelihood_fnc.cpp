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

arma::mat get_P_i(int jj, arma::mat &z_i, arma::mat &zeta) {
    
    double q1_sub = arma::as_scalar(z_i.row(jj) * zeta.col(0));
    double q1 = exp(q1_sub);
    double q2_sub = arma::as_scalar(z_i.row(jj) * zeta.col(1));
    double q2 = exp(q2_sub);
    double q3_sub = arma::as_scalar(z_i.row(jj) * zeta.col(2));
    double q3 = exp(q3_sub);
    double q4_sub = arma::as_scalar(z_i.row(jj) * zeta.col(3));
    double q4 = exp(q4_sub);
    
    double q5_sub = arma::as_scalar(z_i.row(jj) * zeta.col(4));
    double q5 = exp(q5_sub);
    double q6_sub = arma::as_scalar(z_i.row(jj) * zeta.col(5));
    double q6 = exp(q6_sub);
    double q7_sub = arma::as_scalar(z_i.row(jj) * zeta.col(6));
    double q7 = exp(q7_sub);
    double q8_sub = arma::as_scalar(z_i.row(jj) * zeta.col(7));
    double q8 = exp(q8_sub);
    
    double q9_sub = arma::as_scalar(z_i.row(jj) * zeta.col(8));
    double q9 = exp(q9_sub);
    double q10_sub = arma::as_scalar(z_i.row(jj) * zeta.col(9));
    double q10 = exp(q10_sub);
    double q11_sub = arma::as_scalar(z_i.row(jj) * zeta.col(10));
    double q11 = exp(q11_sub);
    double q12_sub = arma::as_scalar(z_i.row(jj) * zeta.col(11));
    double q12 = exp(q12_sub);
    
    arma::mat Q = { {   1,   q1,  0,  q2,  0},
                    {   0,    1, q3,  q4,  0},
                    {  q5,   q6,  1,  q7,  0},
                    {   0,   q8,  0,   1, q9},
                    { q10,  q11,  0, q12,  1}}; 
    
    arma::vec q_row_sums = arma::sum(Q, 1);
    arma::mat P_i = Q.each_col() / q_row_sums;
    
    return P_i;
}

double log_like(const arma::vec &EIDs, const arma::vec &par, 
                const arma::field<arma::uvec> &par_index, 
                const arma::field <arma::vec> &B, const arma::field <arma::vec> &A,
                const arma::field <arma::vec> &W, const arma::mat &Y, const arma::mat &z,
                const arma::field<arma::field<arma::mat>> &Dn_alpha, 
                const arma::field<arma::field<arma::mat>> &Dn_omega,
                const arma::field<arma::field<arma::mat>> &Xn, 
                const arma::mat &gamma_1, int n_cores){
  
    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec alpha_tilde = par.elem(par_index(1) - 1);
    
    arma::vec vec_upsilon = par.elem(par_index(2) - 1);
    arma::mat upsilon_alpha = arma::reshape(vec_upsilon, alpha_tilde.n_elem, 
                                            alpha_tilde.n_elem);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(4) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    
    arma::vec omega_tilde = par.elem(par_index(7) - 1);
    
    arma::vec eta_omega = exp(par.elem(par_index(8) - 1));
    arma::mat upsilon_omega = arma::diagmat(eta_omega);
    
    arma::mat G_sqrt_t = arma::reshape(par.elem(par_index(9) - 1), 4, 4);
    arma::mat G_sqrt = G_sqrt_t.t() * G_sqrt_t;
    arma::mat G = G_sqrt * G_sqrt;
    
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
    
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // "i" is the numeric EID number
        // "ii" is the index of the EID
        
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        arma::mat z_i = z.rows(sub_ind);
        
        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_omega_i = W(ii);
        arma::field<arma::mat> D_alpha_i = Dn_alpha(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::field<arma::mat> X_i = Xn(ii);
        
        double like_comp_i = 0;
        
        for(int jj = 0; jj < n_i; jj++) {
            
            if(jj == 0) {
                
                arma::vec gamma_1_mu = Y_i.col(jj) - D_omega_i(jj) * vec_omega_i  - X_i(jj) * vec_beta;
                arma::mat gamma_1_v = gamma_var + G;
                
                arma::vec log_g_pdf = dmvnorm(gamma_1.row(ii), gamma_1_mu, gamma_1_v, true);
                arma::vec log_a_pdf = dmvnorm(vec_alpha_i.t(), alpha_tilde, upsilon_alpha, true);
                arma::vec log_o_pdf = dmvnorm(vec_omega_i.t(), omega_tilde, upsilon_omega, true);
                
                double log_g_val = arma::as_scalar(log_g_pdf);
                double log_a_val = arma::as_scalar(log_a_pdf);
                double log_o_val = arma::as_scalar(log_o_pdf);
                
                like_comp_i = like_comp_i + log_g_val + log_a_val + log_o_val;
                
            } else {
                
                arma::vec y_j = Y_i.col(jj);
                arma::vec y_j_1 = Y_i.col(jj-1); 
                arma::vec nu_j = gamma_1.row(ii).t() + D_alpha_i(jj) * vec_alpha_i + 
                                    D_omega_i(jj) * vec_omega_i + X_i(jj) * vec_beta;
                arma::vec nu_j_1 = gamma_1.row(ii).t() + D_alpha_i(jj-1) * vec_alpha_i + 
                                    D_omega_i(jj-1) * vec_omega_i + X_i(jj-1) * vec_beta;
                
                arma::vec y_mean = nu_j + A_1 * (y_j_1 - nu_j_1);
                
                arma::vec log_y_pdf = dmvnorm(y_j.t(), y_mean, R, true);
                
                like_comp_i = like_comp_i + arma::as_scalar(log_y_pdf);
                
                if(jj == 1) {
                    like_comp_i = like_comp_i + log(P_init(b_i(jj) - 1));
                } else {
                    arma::mat P_i = get_P_i(jj, z_i, zeta);
                    
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
                const arma::field <arma::vec> &B, const arma::field <arma::vec> &A,
                const arma::field <arma::vec> &W, const arma::mat &Y, const arma::mat &z,
                const arma::field<arma::field<arma::mat>> &Dn_alpha, 
                const arma::field<arma::field<arma::mat>> &Dn_omega,
                const arma::field<arma::field<arma::mat>> &Xn, 
                const arma::mat &gamma_1, int n_cores){

    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
    // Likelihood --------------------------------------------------------------
    double value = log_like(EIDs, par, par_index, B, A, W, Y, z, Dn_alpha, Dn_omega,
                            Xn, gamma_1, n_cores);
    
    // Autocorrelation prior ---------------------------------------------------
    arma::vec vec_A_content = par.elem(par_index(3) - 1);
    arma::vec vec_A_mean(vec_A_content.n_elem, arma::fill::zeros);
    arma::vec scalar_A(vec_A_content.n_elem, arma::fill::ones);
    scalar_A = 100 * scalar_A;
    arma::mat A_var = arma::diagmat(scalar_A);
    
    arma::vec prior_A = dmvnorm(vec_A_content.t(), vec_A_mean, A_var, true);
    double prior_A_val = arma::as_scalar(prior_A);
    
    // Error-variance prior ----------------------------------------------------
    arma::vec vec_R_content = par.elem(par_index(4) - 1);
    arma::vec vec_R_mean(vec_R_content.n_elem, arma::fill::zeros);
    arma::vec scalar_R(vec_R_content.n_elem, arma::fill::ones);
    scalar_R = 100 * scalar_R;
    arma::mat R_var = arma::diagmat(scalar_R);
    
    arma::vec prior_R = dmvnorm(vec_R_content.t(), vec_R_mean, R_var, true);
    double prior_R_val = arma::as_scalar(prior_R);
    
    // Zeta prior --------------------------------------------------------------
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::vec vec_zeta_mean = {-3.7405, 2.5, -4.2152,   1, -2.6473,-0.5, -2.1475, -0.2, 
                               -3.4459,  -1, -2.9404,   1, -3.2151,   1, -3.1778,  1.5, 
                               -2.0523,   0, -3.4459,-0.2, -3.2404, 2.5, -3.2151,    1};
    arma::vec scalar_zeta(vec_zeta_mean.n_elem, arma::fill::ones);
    scalar_zeta = 100 * scalar_zeta;
    arma::mat zeta_var = arma::diagmat(scalar_zeta);
    
    arma::vec prior_zeta = dmvnorm(vec_zeta_content.t(), vec_zeta_mean, zeta_var, true);
    double prior_zeta_val = arma::as_scalar(prior_zeta);
    
    // Initial Probability prior -----------------------------------------------
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec vec_init_mean(vec_init_content.n_elem, arma::fill::zeros);
    arma::vec scalar_init(vec_init_content.n_elem, arma::fill::ones);
    scalar_init = 100 * scalar_init;
    arma::mat init_var = arma::diagmat(scalar_init);

    arma::vec prior_init = dmvnorm(vec_init_content.t(), vec_init_mean, init_var, true);
    double prior_init_val = arma::as_scalar(prior_init);
    
    // Omega RE variance prior -------------------------------------------------
    arma::vec vec_up_omega_content = par.elem(par_index(8) - 1);
    arma::vec omega_mean(vec_up_omega_content.n_elem, arma::fill::zeros);
    arma::vec scalar_omega(vec_up_omega_content.n_elem, arma::fill::ones);
    scalar_omega = 100 * scalar_omega;
    arma::mat omega_var = arma::diagmat(scalar_omega);
    
    arma::vec prior_omega = dmvnorm(vec_up_omega_content.t(), omega_mean, omega_var, true);
    double prior_omega_val = arma::as_scalar(prior_omega);
    
    // Gamma RE variance prior -------------------------------------------------
    arma::vec vec_G_content = par.elem(par_index(9) - 1);
    arma::vec vec_G_mean(vec_G_content.n_elem, arma::fill::zeros);
    arma::vec scalar_G(vec_G_content.n_elem, arma::fill::ones);
    scalar_G = 100 * scalar_G;
    arma::mat G_var = arma::diagmat(scalar_G);
    
    arma::vec prior_G = dmvnorm(vec_G_content.t(), vec_G_mean, G_var, true);
    double prior_G_val = arma::as_scalar(prior_G);
    
    // Full log-posterior ------------------------------------------------------
    value = value + prior_A_val + prior_R_val + prior_zeta_val + prior_init_val + prior_omega_val + prior_G_val;
    
    return value;
}

// [[Rcpp::export]]
arma::mat update_gamma_i(const arma::vec &EIDs, const arma::vec &par,
                         const arma::field<arma::uvec> &par_index,
                         const arma::field <arma::vec> &A,
                         const arma::field <arma::vec> &W, const arma::mat &Y,
                         const arma::field<arma::field<arma::mat>> &Dn_alpha, 
                         const arma::field<arma::field<arma::mat>> &Dn_omega,
                         const arma::field<arma::field<arma::mat>> &Xn, int n_cores){

    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    arma::mat I(A_1.n_rows, A_1.n_cols, arma::fill::eye);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(4) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    arma::mat inv_R = arma::inv_sympd(R);
    
    arma::mat G_sqrt_t = arma::reshape(par.elem(par_index(9) - 1), 4, 4);
    arma::mat G_sqrt = G_sqrt_t.t() * G_sqrt_t;
    arma::mat G = G_sqrt * G_sqrt;
    
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
        
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_omega_i = W(ii);
        arma::field<arma::mat> D_alpha_i = Dn_alpha(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::field<arma::mat> X_i = Xn(ii);

        arma::mat diff_g = I - A_1;
        arma::mat inv_W_i = g_var_combo_inv + (n_i - 1) * (diff_g.t() * inv_R * diff_g);
        arma::mat W_i = arma::inv_sympd(inv_W_i);

        arma::vec V_i = g_var_combo_inv * (Y_i.col(0) - D_omega_i(0) * vec_omega_i  - X_i(0) * vec_beta);
        
        arma::vec int_sum(Y_i.n_rows, arma::fill::zeros);
        for(int k = 1; k < n_i; k++) {
            arma::mat diff_a = D_alpha_i(k) - A_1 * D_alpha_i(k-1);
            arma::mat diff_o = D_omega_i(k) - A_1 * D_omega_i(k-1);
            arma::mat diff_x = X_i(k) - A_1 * X_i(k-1);
            arma::mat diff_y = Y_i.col(k) - A_1 * Y_i.col(k-1);

            arma::vec m_k = diff_y - diff_a*vec_alpha_i - diff_o*vec_omega_i - diff_x*vec_beta;
            int_sum = int_sum + m_k;
        }

        V_i = V_i + diff_g.t() * inv_R * int_sum;
        
        arma::vec mu_gamma_i = W_i * V_i;

        arma::vec gamma_i = arma::mvnrnd(mu_gamma_i, W_i, 1);

        gamma_1.row(ii) = gamma_i.t();
    }

    return gamma_1;
}


// [[Rcpp::export]]
arma::field <arma::vec> update_alpha_i(const arma::vec &EIDs, const arma::vec &par, 
                                       const arma::field<arma::uvec> &par_index,
                                       const arma::field <arma::vec> &W, const arma::mat &Y,
                                       const arma::field<arma::field<arma::mat>> &Dn_alpha, 
                                       const arma::field<arma::field<arma::mat>> &Dn_omega,
                                       const arma::field<arma::field<arma::mat>> &Xn, 
                                       const arma::mat &gamma_1, int n_cores){

    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec alpha_tilde = par.elem(par_index(1) - 1);
    
    arma::vec vec_upsilon = par.elem(par_index(2) - 1);
    arma::mat upsilon_alpha = arma::reshape(vec_upsilon, alpha_tilde.n_elem, 
                                            alpha_tilde.n_elem);
    arma::mat inv_ups_alpha = arma::inv_sympd(upsilon_alpha);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    arma::mat I(A_1.n_rows, A_1.n_cols, arma::fill::eye);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(4) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    arma::mat inv_R = arma::inv_sympd(R);
    
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
        
        arma::vec vec_omega_i = W(ii);
        arma::field<arma::mat> D_alpha_i = Dn_alpha(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::field<arma::mat> X_i = Xn(ii);
        
        arma::mat inv_W_i = inv_ups_alpha;
        arma::vec V_i     = inv_ups_alpha * alpha_tilde;
        
        arma::mat diff_g = I - A_1;
        
        for(int k = 1; k < n_i; k++) {
            
            arma::mat diff_a = D_alpha_i(k) - A_1 * D_alpha_i(k-1);
            arma::mat diff_o = D_omega_i(k) - A_1 * D_omega_i(k-1);
            arma::mat diff_x = X_i(k) - A_1 * X_i(k-1);
            arma::mat diff_y = Y_i.col(k) - A_1 * Y_i.col(k-1);
            
            arma::vec m_k = diff_y - diff_o*vec_omega_i - diff_x*vec_beta - diff_g*gamma_1.row(ii).t();
            
            inv_W_i = inv_W_i + diff_a.t() * inv_R * diff_a;
            
            V_i = V_i + diff_a.t() * inv_R * m_k;
        }
        
        arma::mat W_i = arma::inv_sympd(inv_W_i);
        
        arma::vec mu_alpha_i = W_i * V_i;
        
        arma::vec vec_alpha_i = arma::mvnrnd(mu_alpha_i, W_i, 1);
        
        A(ii) = vec_alpha_i;
    }
    
    return A;
  
}

// *** Using ***
// [[Rcpp::export]]
arma::field <arma::vec> update_omega_i(const arma::vec &EIDs, const arma::vec &par, 
                                       const arma::field<arma::uvec> &par_index,
                                       const arma::field <arma::vec> &A, const arma::mat &Y,
                                       const arma::field<arma::field<arma::mat>> &Dn_alpha, 
                                       const arma::field<arma::field<arma::mat>> &Dn_omega,
                                       const arma::field<arma::field<arma::mat>> &Xn, 
                                       const arma::mat &gamma_1, int n_cores){

    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    arma::mat I(A_1.n_rows, A_1.n_cols, arma::fill::eye);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(4) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    arma::mat inv_R = arma::inv_sympd(R);
    
    arma::vec omega_tilde = par.elem(par_index(7) - 1);
    
    arma::vec eta_omega = exp(par.elem(par_index(8) - 1));
    arma::vec inv_eta_omega = 1 / eta_omega;
    arma::mat inv_ups_omega = arma::diagmat(inv_eta_omega);
    
    arma::mat G_sqrt_t = arma::reshape(par.elem(par_index(9) - 1), 4, 4);
    arma::mat G_sqrt = G_sqrt_t.t() * G_sqrt_t;
    arma::mat G = G_sqrt * G_sqrt;
    
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
    
    arma::mat g_var_combo = gamma_var + G;
    arma::mat g_var_combo_inv = arma::inv_sympd(g_var_combo);
    
    arma::vec eids = Y.col(0);
    // -------------------------------------------------------------------------
    
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
        
        arma::vec vec_alpha_i = A(ii);
        arma::field<arma::mat> D_alpha_i = Dn_alpha(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::field<arma::mat> X_i = Xn(ii);
        
        arma::mat inv_W_i = inv_ups_omega + D_omega_i(0).t() * g_var_combo_inv * D_omega_i(0);
        arma::vec V_i     = inv_ups_omega * omega_tilde + 
            D_omega_i(0).t() * g_var_combo_inv * (Y_i.col(0) - gamma_1.row(ii).t() - X_i(0) * vec_beta);
        
        arma::mat diff_g = I - A_1;

        for(int k = 1; k < n_i; k++) {
            
            arma::mat diff_a = D_alpha_i(k) - A_1 * D_alpha_i(k-1);
            arma::mat diff_o = D_omega_i(k) - A_1 * D_omega_i(k-1);
            arma::mat diff_x = X_i(k) - A_1 * X_i(k-1);
            arma::mat diff_y = Y_i.col(k) - A_1 * Y_i.col(k-1);
            
            arma::vec m_k = diff_y - diff_a*vec_alpha_i - diff_x*vec_beta - diff_g*gamma_1.row(ii).t();
            
            inv_W_i = inv_W_i + diff_o.t() * inv_R * diff_o;
            
            V_i = V_i + diff_o.t() * inv_R * m_k;
        }
        
        arma::mat W_i = arma::inv_sympd(inv_W_i);
        
        arma::vec mu_omega_i = W_i * V_i;
        
        arma::vec vec_omega_i = arma::mvnrnd(mu_omega_i, W_i, 1);
        
        W(ii) = vec_omega_i;
    }
    
    return W;
}

// [[Rcpp::export]]
arma::vec update_beta_upsilon(const arma::vec &EIDs, arma::vec &par, 
                              const arma::field<arma::uvec> &par_index,
                              const arma::field <arma::vec> &A, 
                              const arma::field <arma::vec> &W, const arma::mat &Y,
                              const arma::field<arma::field<arma::mat>> &Dn_alpha, 
                              const arma::field<arma::field<arma::mat>> &Dn_omega,
                              const arma::field<arma::field<arma::mat>> &Xn, 
                              const arma::mat &gamma_1, int n_cores){
    
    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic

    // Initializing parameters -------------------------------------------------
    arma::vec alpha_tilde = par.elem(par_index(1) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    arma::mat I(A_1.n_rows, A_1.n_cols, arma::fill::eye);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(4) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    arma::mat inv_R = arma::inv_sympd(R);
    
    arma::mat G_sqrt_t = arma::reshape(par.elem(par_index(9) - 1), 4, 4);
    arma::mat G_sqrt = G_sqrt_t.t() * G_sqrt_t;
    arma::mat G = G_sqrt * G_sqrt;
    
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
    
    arma::mat g_var_combo = gamma_var + G;
    arma::mat g_var_combo_inv = arma::inv_sympd(g_var_combo);
    
    arma::vec eids = Y.col(0);
    int N_id = EIDs.n_elem;
    
    // Priors ------------------------------------------------------------------
    arma::vec vec_beta_0 = {0.25, -2, 2, -0.25};
    arma::vec scalar_beta = {1, 0.01, 0.01, 0.01};
    arma::mat inv_sigma_beta = arma::diagmat(scalar_beta);
    
    int nu_ups = 20;
    arma::vec scalar_ups = {4, 4, 1, 1, 
                            4, 4, 1, 1, 
                            4, 4, 1, 1, 
                            4, 4, 1, 1};
    scalar_ups = (nu_ups - 16 - 1) * scalar_ups;
    arma::mat psi_ups = arma::diagmat(scalar_ups);
    
    // -------------------------------------------------------------------------

    arma::field<arma::mat> beta_W_list(N_id);
    arma::field<arma::vec> beta_V_list(N_id);
    arma::field<arma::mat> ups_list(N_id);

    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < N_id; ii++) {
        
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_omega_i = W(ii);
        arma::field<arma::mat> D_alpha_i = Dn_alpha(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::field<arma::mat> X_i = Xn(ii);

        arma::mat inv_W_i = X_i(0).t() * g_var_combo_inv * X_i(0);
        arma::vec V_i     = X_i(0).t() * g_var_combo_inv * (Y_i.col(0) - gamma_1.row(ii).t() - D_omega_i(0) * vec_omega_i);
        
        arma::mat diff_g = I - A_1;

        for(int k = 1; k < n_i; k++) {
            arma::mat diff_a = D_alpha_i(k) - A_1 * D_alpha_i(k-1);
            arma::mat diff_o = D_omega_i(k) - A_1 * D_omega_i(k-1);
            arma::mat diff_x = X_i(k) - A_1 * X_i(k-1);
            arma::mat diff_y = Y_i.col(k) - A_1 * Y_i.col(k-1);
            
            arma::vec m_k = diff_y - diff_a*vec_alpha_i - diff_o*vec_omega_i - diff_g*gamma_1.row(ii).t();
            
            inv_W_i = inv_W_i + diff_x.t() * inv_R * diff_x;
            
            V_i = V_i + diff_x.t() * inv_R * m_k;
        }

        arma::mat hold_ups = vec_alpha_i - alpha_tilde;
        arma::mat in_Upsilon_cov = hold_ups * hold_ups.t();

        beta_W_list(ii) = inv_W_i;
        beta_V_list(ii) = V_i;
        ups_list(ii) = in_Upsilon_cov;
    }

    arma::mat sum_beta_W = inv_sigma_beta; 
    arma::vec sum_beta_V = inv_sigma_beta * vec_beta_0;
    arma::mat sum_ups_cov = psi_ups;

    for(int ii = 0; ii < N_id; ii++){
        sum_beta_W += beta_W_list(ii);
        sum_beta_V += beta_V_list(ii);
        sum_ups_cov += ups_list(ii);
    }

    arma::mat W_b = arma::inv_sympd(sum_beta_W);
    arma::vec V_b = sum_beta_V;

    par.elem(par_index(0) - 1) = arma::mvnrnd(W_b*V_b, W_b, 1); 
    par.elem(par_index(2) - 1) = arma::vectorise(riwish(nu_ups + N_id, sum_ups_cov));

    return par;
}


// [[Rcpp::export]]
arma::vec update_alpha_tilde(const arma::vec &EIDs, arma::vec &par, 
                             const arma::field<arma::uvec> &par_index,
                             const arma::field <arma::vec> &A, const arma::mat &Y){

    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_upsilon = par.elem(par_index(2) - 1);
    arma::mat upsilon_alpha = arma::reshape(vec_upsilon, par_index(1).n_elem, 
                                            par_index(1).n_elem);
    arma::mat inv_ups_alpha = arma::inv_sympd(upsilon_alpha);
    
    // Priors ------------------------------------------------------------------
    arma::vec alpha_0(par_index(1).n_elem, arma::fill::zeros);
    arma::vec inv_sigma_alpha_diag(alpha_0.n_elem, arma::fill::ones);
    inv_sigma_alpha_diag = 0.01 * inv_sigma_alpha_diag;
    arma::mat inv_sigma_alpha = arma::diagmat(inv_sigma_alpha_diag);
    
    // -------------------------------------------------------------------------
    
    int N_id = EIDs.n_elem;
    
    arma::mat inv_W = inv_sigma_alpha + N_id * inv_ups_alpha;
    arma::mat W_a = arma::inv_sympd(inv_W);
    
    arma::vec total_alpha = A(0);
    for (int ii = 1; ii < N_id; ii++){
        total_alpha = total_alpha + A(ii);
    }
    
    arma::vec V = inv_sigma_alpha * alpha_0 + inv_ups_alpha * total_alpha;
    
    arma::vec mu = W_a * V;
    
    par.elem(par_index(1) - 1) = arma::mvnrnd(mu, W_a, 1);
    
    return par;
}

// [[Rcpp::export]]
arma::vec update_omega_tilde(const arma::vec &EIDs, arma::vec &par, 
                             const arma::field<arma::uvec> &par_index,
                             const arma::field <arma::vec> &W, const arma::mat &Y){

    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
    // Initializing parameters -------------------------------------------------
    arma::vec eta_omega = exp(par.elem(par_index(8) - 1));
    arma::vec inv_eta_omega = 1 / eta_omega;
    arma::mat inv_ups_omega = arma::diagmat(inv_eta_omega);
    
    // Priors ------------------------------------------------------------------
    arma::vec omega_0 = {-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                         -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                         -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                         -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1};
    omega_0 = 2 * omega_0;
    
    arma::vec inv_sigma_omega_diag(omega_0.n_elem, arma::fill::ones);
    inv_sigma_omega_diag = 0.01 * inv_sigma_omega_diag;
    arma::mat inv_sigma_omega = arma::diagmat(inv_sigma_omega_diag);
    
    // -------------------------------------------------------------------------
    
    int N_id = EIDs.n_elem;
    
    arma::mat inv_W = inv_sigma_omega + N_id * inv_ups_omega;
    arma::mat W_o = arma::inv_sympd(inv_W);
    
    arma::vec total_omega = W(0);
    for (int ii = 1; ii < N_id; ii++){
        total_omega = total_omega + W(ii);
    }
    
    arma::vec V = inv_sigma_omega * omega_0 + inv_ups_omega * total_omega;
    
    arma::vec mu = W_o * V;
    
    par.elem(par_index(7) - 1) = arma::mvnrnd(mu, W_o, 1);
    
    return par;
}

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

// [[Rcpp::export]]
arma::field<arma::field<arma::mat>> initialize_X(const arma::vec EIDs,
                                                 const arma::mat &Y,
                                                 const arma::mat &x) {
    
    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
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

// [[Rcpp::export]]
arma::field<arma::field<arma::mat>> initialize_Dn(const arma::vec EIDs,
                                                  arma::field <arma::vec> &B) {
    
    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    
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

    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic

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
Rcpp::List mle_state_seq(const arma::vec &EIDs, const arma::vec &par, 
                         const arma::field<arma::uvec> &par_index, 
                         const arma::field <arma::vec> &A,
                         const arma::field <arma::vec> &W, 
                         const arma::mat &Y, const arma::mat &z,
                         const arma::field<arma::field<arma::mat>> &Dn_omega,
                         const arma::field<arma::field<arma::mat>> &Xn, 
                         const arma::mat &gamma_1, int n_cores) {
    
    // par_index: (0) beta, (1) alpha_tilde, (2) upsilon, (3) A, (4) R, (5) zeta, 
    //            (6) init, (7) omega_tilde, (8) eta_omega, (9) G
    // Y: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic

    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    // Initializing parameters -------------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));
    arma::mat A_1 = arma::diagmat(vec_A);
    
    arma::mat R_sqrt_t = arma::reshape(par.elem(par_index(4) - 1), 4, 4);
    arma::mat R_sqrt = R_sqrt_t.t() * R_sqrt_t;
    arma::mat R = R_sqrt * R_sqrt;
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    
    arma::vec eids = Y.col(0);
    arma::vec clinic_rule_vec = Y.col(6);
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
        arma::mat z_i = z.rows(sub_ind);
        
        arma::vec vec_omega_i = W(ii);
        arma::field<arma::mat> D_omega_i = Dn_omega(ii);
        arma::field<arma::mat> X_i = Xn(ii);
        
        arma::vec vec_alpha_i_no_base = A(ii);
        arma::mat alpha_i_no_base = arma::reshape(vec_alpha_i_no_base, 4, 4);
        
        arma::mat alpha_i(5, 4, arma::fill::zeros);
        alpha_i.row(0) = Y_i.col(0).t();
        alpha_i.row(1) = alpha_i_no_base.row(0);
        alpha_i.row(2) = alpha_i_no_base.row(1);
        alpha_i.row(3) = alpha_i_no_base.row(2);
        alpha_i.row(4) = alpha_i_no_base.row(3);
        
        arma::vec b_i(n_i, arma::fill::zeros);

        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        arma::imat adj_mat_i;
        if(clinic_rule < 0) {
            adj_mat_i = adj_mat_sub_GLOBAL;
        } else {
            adj_mat_i = adj_mat_GLOBAL;
        }
        
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
                    
                    arma::vec g_plus_a = alpha_i.row(0).t() + 
                        arma::accu(twos)   * alpha_i.row(1).t() + 
                        arma::accu(threes) * alpha_i.row(2).t() + 
                        arma::accu(fours)  * alpha_i.row(3).t() + 
                        arma::accu(fives)  * alpha_i.row(4).t();
                    
                    arma::vec g_plus_a_1 = alpha_i.row(0).t();
                    
                    arma::vec nu_k = g_plus_a + D_omega_i(k)*vec_omega_i + X_i(k)*vec_beta;
                    arma::vec nu_k_1 = g_plus_a_1 + D_omega_i(k-1)*vec_omega_i + X_i(k-1)*vec_beta;
                    
                    arma::vec mean_b = nu_k + A_1 * (Y_i.col(k-1) - nu_k_1);
                    
                    arma::vec log_y_pdf = dmvnorm(Y_i.col(k).t(), mean_b, R, true);
                    
                    init_vals(jj) = log(P_init(jj)) + arma::as_scalar(log_y_pdf);
                }
                
                int selected_state;
                
                if(clinic_rule < 0) {
                    arma::vec poss_states = {1,4,5};
                    init_vals = {init_vals(0), init_vals(3), init_vals(4)};
                    
                    selected_state = poss_states(arma::index_max(init_vals));
                } else {
                    selected_state = arma::index_max(init_vals) + 1;    
                }

                b_i(k) = selected_state;

            } else {
                
                arma::mat P_i = get_P_i(k, z_i, zeta);
                
                // All other states --------------------------------------------
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
                        
                        arma::vec g_plus_a = alpha_i.row(0).t() + 
                            arma::accu(twos)   * alpha_i.row(1).t() + 
                            arma::accu(threes) * alpha_i.row(2).t() + 
                            arma::accu(fours)  * alpha_i.row(3).t() + 
                            arma::accu(fives)  * alpha_i.row(4).t();
                        
                        arma::vec g_plus_a_1 = alpha_i.row(0).t() + 
                            arma::accu(twos.subvec(0,k-2))   * alpha_i.row(1).t() + 
                            arma::accu(threes.subvec(0,k-2)) * alpha_i.row(2).t() + 
                            arma::accu(fours.subvec(0,k-2))  * alpha_i.row(3).t() + 
                            arma::accu(fives.subvec(0,k-2))  * alpha_i.row(4).t();
                        
                        arma::vec nu_k = g_plus_a + D_omega_i(k)*vec_omega_i + X_i(k)*vec_beta;
                        arma::vec nu_k_1 = g_plus_a_1 + D_omega_i(k-1)*vec_omega_i + X_i(k-1)*vec_beta;
                        
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

// Miscellaneous ---------------------------------------------------------------

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

    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::vec std_R = sqrt(arma::diagvec(R));

    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12);

    arma::vec vec_init = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init(0)), exp(vec_init(1)),
                            exp(vec_init(2)), exp(vec_init(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);
    // arma::vec P_init = {1, 0, 0, 0, 0};
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

                int selected_state;
                
                if(clinic_rule < 0) {
                    arma::vec poss_states = {1,4,5};
                    init_vals = {init_vals(0), init_vals(3), init_vals(4)};
                    
                    selected_state = poss_states(arma::index_max(init_vals));
                } else {
                    selected_state = arma::index_max(init_vals) + 1;    
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

                    // // Prevent negatives
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

                    // // Prevent negatives
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

                    // // Prevent negatives
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

        newY.rows(sub_ind) = Y_i_new.t();
    }

    Y.cols(1,4) = newY;
    return Y;
}

// [[Rcpp::export]]
Rcpp::List proposal_R_cpp_new(const int nu_R, const arma::mat psi_R, arma::mat curr_R,
                              const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn,
                              const arma::field<arma::field<arma::mat>> &Xn,
                              const arma::field <arma::vec> A, const arma::vec par,
                              const arma::field<arma::uvec> par_index,
                              const arma::vec EIDs, arma::field <arma::vec> B,
                              const arma::field<arma::field<arma::mat>> Dn_omega,
                              const arma::field <arma::vec> W, int n_cores) {
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega

    // Initializing parameters -------------------------------------------------
    arma::vec vec_beta = par.elem(par_index(0) - 1);

    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A = exp(vec_A_total) / (1 + exp(vec_A_total));

    arma::mat curr_R_sqrt = arma::sqrtmat_sympd(curr_R);
    // -------------------------------------------------------------------------

    arma::vec eids = Y.col(0);
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
        arma::vec vec_omega_i = W(ii);

        arma::field<arma::mat> Xn_i = Xn(ii);
        arma::field<arma::mat> Dn_i = Dn(ii);
        arma::field<arma::mat> Dn_omega_i = Dn_omega(ii);

        arma::mat psi_q_i(4, 4, arma::fill::zeros);

        for(int k = 0; k < Y_i.n_cols; k++) {

            arma::vec nu_k = Dn_i(k) * vec_alpha_i + Dn_omega_i(k) * vec_omega_i + Xn_i(k) * vec_beta;
            arma::vec y_diff_k = Y_i.col(k) - nu_k;
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

                // psi_q = psi_q + y_diff_temp * y_diff_temp.t();
                psi_q_i = psi_q_i + y_diff_temp * y_diff_temp.t();

            } else {
                arma::vec nu_k_1 = Dn_i(k-1) * vec_alpha_i +
                    Dn_omega_i(k-1) * vec_omega_i +
                    Xn_i(k-1) * vec_beta;
                arma::vec y_diff_k_1 = Y_i.col(k-1) - nu_k_1;

                arma::vec hold_k = y_diff_k - A_k * y_diff_k_1;

                // psi_q = psi_q + hold_k * hold_k.t();
                psi_q_i = psi_q_i + hold_k * hold_k.t();
            }
        }
        psi_q_list(ii) = psi_q_i;
    }

    arma::mat psi_q = psi_R;
    for(int ii = 0; ii < EIDs.n_elem; ii++) {
        psi_q += psi_q_list(ii);
    }
    // psi_q = psi_q + psi_R;

    int nu_q = Y.n_rows + nu_R;

    List nu_psi_R = List::create(psi_q, nu_q);
    return nu_psi_R;
}

// [[Rcpp::export]]
void test_fnc(int states_per_step) {

    arma::vec temp = {log(2),log(2),log(2),log(2)};
    arma::vec eta_omega = exp(temp);
    arma::vec inv_eta_omega = 1 / eta_omega;
    arma::mat inv_ups_omega = arma::diagmat(inv_eta_omega);
    Rcpp::Rcout << inv_eta_omega << std::endl;
    Rcpp::Rcout << inv_ups_omega << std::endl;
    
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
