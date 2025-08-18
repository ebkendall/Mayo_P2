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

// [[Rcpp::export]]
void initialize_cpp(arma::imat a_mat, arma::imat a_mat_sub, int sps) {
    adj_mat_GLOBAL = a_mat;
    adj_mat_sub_GLOBAL = a_mat_sub;

    Omega_List_GLOBAL_multi = get_Omega_list(adj_mat_GLOBAL, sps);
    Omega_List_GLOBAL_sub_multi = get_Omega_list(adj_mat_sub_GLOBAL, sps);
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

// [[Rcpp::export]]
arma::field<arma::field<arma::mat>> initialize_Xn(const arma::vec EIDs,
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
            // Guarantee Dn_alpha(0) = 0
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
// -----------------------------------------------------------------------------

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

                arma::vec check_vec = s_i.subvec(0, first_bleed_ind);

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

                arma::vec check_vec = s_i.subvec(0, first_bleed_ind);

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

                arma::vec check_vec = s_i.subvec(0, first_bleed_ind);
                arma::vec check_vec_ind = arma::linspace(0, first_bleed_ind, first_bleed_ind + 1);
                
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

                arma::vec check_vec = s_i.subvec(0, first_bleed_ind);
                arma::vec check_vec_ind = arma::linspace(0, first_bleed_ind, first_bleed_ind + 1);
                
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

arma::vec p_2_sampler(int n_i, arma::mat &Y_i, arma::mat z_i, arma::imat adj_mat_i,
                      arma::vec &b_i, arma::mat alpha_i, arma::vec vec_omega_i,
                      arma::vec vec_beta, arma::mat A_1, arma::mat R,
                      arma::field<arma::mat> D_omega_i, arma::field<arma::mat> X_i,
                      arma::mat zeta, arma::vec P_init, int clinic_rule,
                      int rbc_rule, arma::vec bleed_ind_i) {

    bool clinic_sub;
    if(clinic_rule < 0) {
        clinic_sub = true;
    } else {
        clinic_sub = false;
    }

    for(int k = 0; k < n_i - 1; k++) {

        arma::vec s_i = b_i;

        // arma::mat omega_set = get_omega_list(k-1, n_i - 1, b_i.subvec(1, n_i-1), 2, clinic_sub);
        arma::mat omega_set = get_omega_list(k, n_i, b_i, 2, clinic_sub);

        arma::vec prob_omega(omega_set.n_rows, arma::fill::ones);
        prob_omega = (1/arma::accu(prob_omega)) * prob_omega;
        arma::vec ind_omega = arma::linspace(0, omega_set.n_rows-1, omega_set.n_rows);
        arma::vec row_omega = RcppArmadillo::sample(ind_omega, 1, false, prob_omega);

        s_i.subvec(k, k + 1) = omega_set.row(row_omega(0)).t();

        // Step 3: compute MH-ratio to accept/reject -------------------
        if(arma::accu(arma::abs(s_i - b_i)) != 0) {

            bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i, 2, k);

            if(eval_like) {
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

                    if(t == 0) {
                        log_like_s = log_like_s + log(P_init(s_i(t) - 1));
                        log_like_b = log_like_b + log(P_init(b_i(t) - 1));
                    } else {
                        // INDEXING STARTS AT 1 (NOT 0)
                        // Computations for mean of candidate (s_i) ------------
                        arma::vec s_2 = twos_s.subvec(1, t);
                        arma::vec s_3 = threes_s.subvec(1, t);
                        arma::vec s_4 = fours_s.subvec(1, t);
                        arma::vec s_5 = fives_s.subvec(1, t);

                        arma::vec g_p_a_s = alpha_i.row(0).t() +
                            arma::accu(s_2) * alpha_i.row(1).t() +
                            arma::accu(s_3) * alpha_i.row(2).t() +
                            arma::accu(s_4) * alpha_i.row(3).t() +
                            arma::accu(s_5) * alpha_i.row(4).t();

                        arma::vec g_p_a_s_1;
                        if(t == 1) {
                            g_p_a_s_1 = alpha_i.row(0).t();
                        } else {
                            arma::vec s_2_1 = twos_s.subvec(1, t-1);
                            arma::vec s_3_1 = threes_s.subvec(1, t-1);
                            arma::vec s_4_1 = fours_s.subvec(1, t-1);
                            arma::vec s_5_1 = fives_s.subvec(1, t-1);

                            g_p_a_s_1 = alpha_i.row(0).t() +
                                arma::accu(s_2_1) * alpha_i.row(1).t() +
                                arma::accu(s_3_1) * alpha_i.row(2).t() +
                                arma::accu(s_4_1) * alpha_i.row(3).t() +
                                arma::accu(s_5_1) * alpha_i.row(4).t();
                        }

                        arma::vec nu_s = g_p_a_s + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                        arma::vec nu_s_1 = g_p_a_s_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;

                        arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);

                        // Computations for mean of current (b_i) --------------
                        arma::vec b_2 = twos_b.subvec(1, t);
                        arma::vec b_3 = threes_b.subvec(1, t);
                        arma::vec b_4 = fours_b.subvec(1, t);
                        arma::vec b_5 = fives_b.subvec(1, t);

                        arma::vec g_p_a_b = alpha_i.row(0).t() +
                            arma::accu(b_2) * alpha_i.row(1).t() +
                            arma::accu(b_3) * alpha_i.row(2).t() +
                            arma::accu(b_4) * alpha_i.row(3).t() +
                            arma::accu(b_5) * alpha_i.row(4).t();

                        arma::vec g_p_a_b_1;
                        if(t == 1) {
                            g_p_a_b_1 = alpha_i.row(0).t();
                        } else {
                            arma::vec b_2_1 = twos_b.subvec(1, t-1);
                            arma::vec b_3_1 = threes_b.subvec(1, t-1);
                            arma::vec b_4_1 = fours_b.subvec(1, t-1);
                            arma::vec b_5_1 = fives_b.subvec(1, t-1);

                            g_p_a_b_1 = alpha_i.row(0).t() +
                                arma::accu(b_2_1) * alpha_i.row(1).t() +
                                arma::accu(b_3_1) * alpha_i.row(2).t() +
                                arma::accu(b_4_1) * alpha_i.row(3).t() +
                                arma::accu(b_5_1) * alpha_i.row(4).t();
                        }

                        arma::vec nu_b = g_p_a_b + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                        arma::vec nu_b_1 = g_p_a_b_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;

                        arma::vec mean_b = nu_b + A_1 * (Y_i.col(t-1) - nu_b_1);

                        arma::vec log_y_pdf_s = dmvnorm(Y_i.col(t).t(), mean_s, R, true);
                        arma::vec log_y_pdf_b = dmvnorm(Y_i.col(t).t(), mean_b, R, true);

                        log_like_s = log_like_s + arma::as_scalar(log_y_pdf_s);
                        log_like_b = log_like_b + arma::as_scalar(log_y_pdf_b);

                        if(t <= k + 2) {
                            arma::mat P_i = get_P_i(t, z_i, zeta);

                            log_like_s = log_like_s + log(P_i(s_i(t-1)-1, s_i(t)-1));
                            log_like_b = log_like_b + log(P_i(b_i(t-1)-1, b_i(t)-1));
                        }
                    }
                }

                double diff_check = log_like_s - log_like_b;

                double min_log = log(arma::randu(arma::distr_param(0,1)));
                if(diff_check > min_log){b_i = s_i;}
            }
        }
    }

    return b_i;
}

arma::vec p_full_sampler(int n_i, arma::mat &Y_i, arma::mat z_i, arma::imat adj_mat_i,
                         arma::vec &b_i, arma::mat alpha_i, arma::vec vec_omega_i,
                         arma::vec vec_beta, arma::mat A_1, arma::mat R,
                         arma::field<arma::mat> D_omega_i, arma::field<arma::mat> X_i,
                         arma::mat zeta, arma::vec P_init, int clinic_rule,
                         int rbc_rule, arma::vec bleed_ind_i) {

    arma::vec all_like_vals_b(n_i, arma::fill::zeros);
    arma::vec all_like_vals_s(n_i, arma::fill::zeros);

    arma::vec s_i(n_i, arma::fill::zeros);

    // Propose a full state sequence -------------------------------------------
    for (int k = 0; k < n_i; k++) {

        arma::vec like_vals_s;
        arma::vec like_vals_b;
        arma::vec ss_ind;

        // INDEXING STARTS AT 1 (NOT 0)
        if(k == 0) {
            like_vals_s = log(P_init);
            like_vals_b = log(P_init);

            ss_ind = arma::linspace(0, adj_mat_i.n_cols-1, adj_mat_i.n_cols);

            if(clinic_rule < 0) {
                // Remove state 2 & 3 as a possibility
                like_vals_s = {like_vals_s(0), like_vals_s(3), like_vals_s(4)};
                like_vals_b = {like_vals_b(0), like_vals_b(3), like_vals_b(4)};
                ss_ind = {ss_ind(0), ss_ind(3), ss_ind(4)};
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

            arma::mat P_i = get_P_i(k, z_i, zeta);

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

                    arma::vec g_p_a_s = alpha_i.row(0).t() +
                        arma::accu(twos_s)   * alpha_i.row(1).t() +
                        arma::accu(threes_s) * alpha_i.row(2).t() +
                        arma::accu(fours_s)  * alpha_i.row(3).t() +
                        arma::accu(fives_s)  * alpha_i.row(4).t();

                    arma::vec g_p_a_s_1;
                    if(k == 1) {
                        g_p_a_s_1 = alpha_i.row(0).t();
                    } else {
                        g_p_a_s_1 = alpha_i.row(0).t() +
                            arma::accu(twos_s.subvec(0,k-2))   * alpha_i.row(1).t() +
                            arma::accu(threes_s.subvec(0,k-2)) * alpha_i.row(2).t() +
                            arma::accu(fours_s.subvec(0,k-2))  * alpha_i.row(3).t() +
                            arma::accu(fives_s.subvec(0,k-2))  * alpha_i.row(4).t();
                    }

                    arma::vec nu_s = g_p_a_s + D_omega_i(k)*vec_omega_i + X_i(k)*vec_beta;
                    arma::vec nu_s_1 = g_p_a_s_1 + D_omega_i(k-1)*vec_omega_i + X_i(k-1)*vec_beta;

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

                    arma::vec g_p_a_b = alpha_i.row(0).t() +
                        arma::accu(twos_b)   * alpha_i.row(1).t() +
                        arma::accu(threes_b) * alpha_i.row(2).t() +
                        arma::accu(fours_b)  * alpha_i.row(3).t() +
                        arma::accu(fives_b)  * alpha_i.row(4).t();

                    arma::vec g_p_a_b_1;
                    if(k == 1) {
                        g_p_a_b_1 = alpha_i.row(0).t();
                    } else {
                        g_p_a_b_1 = alpha_i.row(0).t() +
                            arma::accu(twos_b.subvec(0,k-2))   * alpha_i.row(1).t() +
                            arma::accu(threes_b.subvec(0,k-2)) * alpha_i.row(2).t() +
                            arma::accu(fours_b.subvec(0,k-2))  * alpha_i.row(3).t() +
                            arma::accu(fives_b.subvec(0,k-2))  * alpha_i.row(4).t();
                    }

                    arma::vec nu_b = g_p_a_b + D_omega_i(k)*vec_omega_i + X_i(k)*vec_beta;
                    arma::vec nu_b_1 = g_p_a_b_1 + D_omega_i(k-1)*vec_omega_i + X_i(k-1)*vec_beta;

                    arma::vec mean_b = nu_b + A_1 * (Y_i.col(k-1) - nu_b_1);

                    arma::vec log_y_pdf = dmvnorm(Y_i.col(k).t(), mean_b, R, true);

                    like_vals_b(b_ind) = log(P_i(prev_state_b-1, m)) + arma::as_scalar(log_y_pdf);
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
    if(diff_check > min_log) {
        bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i, -1, -1);
        if(eval_like) {
            b_i = s_i;
        }
    }

    return b_i;
}

arma::vec p_flex_sampler(int n_i, arma::mat &Y_i, arma::mat z_i, arma::imat adj_mat_i,
                         arma::vec &b_i, arma::mat alpha_i, arma::vec vec_omega_i,
                         arma::vec vec_beta, arma::mat A_1, arma::mat R,
                         arma::field<arma::mat> D_omega_i, arma::field<arma::mat> X_i,
                         arma::mat zeta, arma::vec P_init, int clinic_rule,
                         int rbc_rule, arma::vec bleed_ind_i, int sps) {

    bool clinic_sub;
    if(clinic_rule < 0) {
        clinic_sub = true;
    } else {
        clinic_sub = false;
    }

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

                like_vals_s = log(P_init);
                like_vals_b = log(P_init);

                ss_ind = arma::linspace(0, adj_mat_i.n_cols-1, adj_mat_i.n_cols);

                if(clinic_rule < 0) {
                    // Remove state 2 & 3 as a possibility
                    like_vals_s = {like_vals_s(0), like_vals_s(3), like_vals_s(4)};
                    like_vals_b = {like_vals_b(0), like_vals_b(3), like_vals_b(4)};
                    ss_ind = {ss_ind(0), ss_ind(3), ss_ind(4)};
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

                arma::mat P_i = get_P_i(t, z_i, zeta);

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

                        arma::vec g_p_a_s = alpha_i.row(0).t() +
                            arma::accu(twos_s)   * alpha_i.row(1).t() +
                            arma::accu(threes_s) * alpha_i.row(2).t() +
                            arma::accu(fours_s)  * alpha_i.row(3).t() +
                            arma::accu(fives_s)  * alpha_i.row(4).t();

                        arma::vec g_p_a_s_1;
                        if(t == 1) {
                            g_p_a_s_1 = alpha_i.row(0).t();
                        } else {
                            g_p_a_s_1 = alpha_i.row(0).t() +
                                arma::accu(twos_s.subvec(0,t-2))   * alpha_i.row(1).t() +
                                arma::accu(threes_s.subvec(0,t-2)) * alpha_i.row(2).t() +
                                arma::accu(fours_s.subvec(0,t-2))  * alpha_i.row(3).t() +
                                arma::accu(fives_s.subvec(0,t-2))  * alpha_i.row(4).t();
                        }

                        arma::vec nu_s = g_p_a_s + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                        arma::vec nu_s_1 = g_p_a_s_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;

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

                        arma::vec g_p_a_b = alpha_i.row(0).t() +
                            arma::accu(twos_b)   * alpha_i.row(1).t() +
                            arma::accu(threes_b) * alpha_i.row(2).t() +
                            arma::accu(fours_b)  * alpha_i.row(3).t() +
                            arma::accu(fives_b)  * alpha_i.row(4).t();

                        arma::vec g_p_a_b_1;
                        if(t == 1) {
                            g_p_a_b_1 = alpha_i.row(0).t();
                        } else {
                            g_p_a_b_1 = alpha_i.row(0).t() +
                                arma::accu(twos_b.subvec(0,t-2))   * alpha_i.row(1).t() +
                                arma::accu(threes_b.subvec(0,t-2)) * alpha_i.row(2).t() +
                                arma::accu(fours_b.subvec(0,t-2))  * alpha_i.row(3).t() +
                                arma::accu(fives_b.subvec(0,t-2))  * alpha_i.row(4).t();
                        }

                        arma::vec nu_b = g_p_a_b + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                        arma::vec nu_b_1 = g_p_a_b_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;

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
        arma::mat Omega_set_s = get_omega_list(t_max, n_i, s_i, 2, clinic_sub);
        arma::mat Omega_set_b = get_omega_list(t_max, n_i, b_i, 2, clinic_sub);

        arma::vec prob_omega_s(Omega_set_s.n_rows, arma::fill::ones);
        prob_omega_s = (1/arma::accu(prob_omega_s)) * prob_omega_s;
        arma::vec ind_omega_s = arma::linspace(0, Omega_set_s.n_rows-1, Omega_set_s.n_rows);
        arma::vec row_omega_s = RcppArmadillo::sample(ind_omega_s, 1, false, prob_omega_s);

        s_i.subvec(t_max, t_max + 1) = Omega_set_s.row(row_omega_s(0)).t();

        // Step 3: compute MH-ratio to accept/reject -------------------
        if(arma::accu(arma::abs(s_i - b_i)) != 0) {

            bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i, 2, t_max);

            if(eval_like) {

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

                    arma::mat P_i = get_P_i(t, z_i, zeta);

                    // INDEXING STARTS AT 1 (NOT 0)
                    // Computations for mean of candidate (s_i) ------------
                    arma::vec s_2 = twos_s.subvec(1, t);
                    arma::vec s_3 = threes_s.subvec(1, t);
                    arma::vec s_4 = fours_s.subvec(1, t);
                    arma::vec s_5 = fives_s.subvec(1, t);

                    arma::vec g_p_a_s = alpha_i.row(0).t() +
                        arma::accu(s_2) * alpha_i.row(1).t() +
                        arma::accu(s_3) * alpha_i.row(2).t() +
                        arma::accu(s_4) * alpha_i.row(3).t() +
                        arma::accu(s_5) * alpha_i.row(4).t();

                    arma::vec g_p_a_s_1;
                    if(t == 1) {
                        g_p_a_s_1 = alpha_i.row(0).t();
                    } else {
                        arma::vec s_2_1 = twos_s.subvec(1, t-1);
                        arma::vec s_3_1 = threes_s.subvec(1, t-1);
                        arma::vec s_4_1 = fours_s.subvec(1, t-1);
                        arma::vec s_5_1 = fives_s.subvec(1, t-1);

                        g_p_a_s_1 = alpha_i.row(0).t() +
                            arma::accu(s_2_1) * alpha_i.row(1).t() +
                            arma::accu(s_3_1) * alpha_i.row(2).t() +
                            arma::accu(s_4_1) * alpha_i.row(3).t() +
                            arma::accu(s_5_1) * alpha_i.row(4).t();
                    }

                    arma::vec nu_s = g_p_a_s + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                    arma::vec nu_s_1 = g_p_a_s_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;

                    arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);

                    // Computations for mean of current (b_i) --------------
                    arma::vec b_2 = twos_b.subvec(1, t);
                    arma::vec b_3 = threes_b.subvec(1, t);
                    arma::vec b_4 = fours_b.subvec(1, t);
                    arma::vec b_5 = fives_b.subvec(1, t);

                    arma::vec g_p_a_b = alpha_i.row(0).t() +
                        arma::accu(b_2) * alpha_i.row(1).t() +
                        arma::accu(b_3) * alpha_i.row(2).t() +
                        arma::accu(b_4) * alpha_i.row(3).t() +
                        arma::accu(b_5) * alpha_i.row(4).t();

                    arma::vec g_p_a_b_1;
                    if(t == 1) {
                        g_p_a_b_1 = alpha_i.row(0).t();
                    } else {
                        arma::vec b_2_1 = twos_b.subvec(1, t-1);
                        arma::vec b_3_1 = threes_b.subvec(1, t-1);
                        arma::vec b_4_1 = fours_b.subvec(1, t-1);
                        arma::vec b_5_1 = fives_b.subvec(1, t-1);

                        g_p_a_b_1 = alpha_i.row(0).t() +
                            arma::accu(b_2_1) * alpha_i.row(1).t() +
                            arma::accu(b_3_1) * alpha_i.row(2).t() +
                            arma::accu(b_4_1) * alpha_i.row(3).t() +
                            arma::accu(b_5_1) * alpha_i.row(4).t();
                    }

                    arma::vec nu_b = g_p_a_b + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                    arma::vec nu_b_1 = g_p_a_b_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;

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
        }

        k = k + sps - 2;
    }

    return b_i;
}

// [[Rcpp::export]]
Rcpp::List state_sampler(const arma::vec EIDs, const arma::vec &par,
                         const arma::field<arma::uvec> &par_index,
                         const arma::field <arma::vec> &B, const arma::field <arma::vec> &A,
                         const arma::field <arma::vec> &W, const arma::mat &Y, const arma::mat &z,
                         const arma::field<arma::field<arma::mat>> &Dn_omega,
                         const arma::field<arma::field<arma::mat>> &Xn,
                         const arma::vec &bleed_indicator, 
                         int states_per_step, int n_cores) {

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

    arma::mat R = arma::reshape(par.elem(par_index(4) - 1), 4, 4);

    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12);

    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);

    arma::vec eids = Y.col(0);
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6);
    // -------------------------------------------------------------------------

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
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);

        arma::vec vec_alpha_i_no_base = A(ii);
        arma::mat alpha_i_no_base = arma::reshape(vec_alpha_i_no_base, 4, 4);
        
        arma::vec alpha_base = Y_i.col(0) - D_omega_i(0)*vec_omega_i - X_i(0)*vec_beta;

        arma::mat alpha_i(5, 4, arma::fill::zeros);
        alpha_i.row(0) = alpha_base.t(); // NOT gamma_i
        alpha_i.row(1) = alpha_i_no_base.row(0);
        alpha_i.row(2) = alpha_i_no_base.row(1);
        alpha_i.row(3) = alpha_i_no_base.row(2);
        alpha_i.row(4) = alpha_i_no_base.row(3);

        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());

        arma::imat adj_mat_i;
        if(clinic_rule < 0) {
            adj_mat_i = adj_mat_sub_GLOBAL;
        } else {
            adj_mat_i = adj_mat_GLOBAL;
        }

        arma::vec b_i = B(ii);
        arma::vec new_b_i;

        if(states_per_step == 2) {
            // p = 2
            new_b_i = p_2_sampler(n_i, Y_i, z_i, adj_mat_i, b_i, alpha_i,
                                  vec_omega_i, vec_beta, A_1, R, D_omega_i, X_i,
                                  zeta, P_init, clinic_rule, rbc_rule, bleed_ind_i);

        } else if(states_per_step >= n_i) {
            new_b_i = p_full_sampler(n_i, Y_i, z_i, adj_mat_i, b_i, alpha_i,
                                     vec_omega_i, vec_beta, A_1, R, D_omega_i,
                                     X_i, zeta, P_init, clinic_rule, rbc_rule, bleed_ind_i);

        } else {
            // 2 < p < n_i
            new_b_i = p_flex_sampler(n_i, Y_i, z_i, adj_mat_i, b_i, alpha_i,
                                     vec_omega_i, vec_beta, A_1, R, D_omega_i,
                                     X_i, zeta, P_init, clinic_rule, rbc_rule, bleed_ind_i,
                                     states_per_step);
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
            // Guarantee Dn_alpha(0) = 0
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
Rcpp::List state_coin_flip(const arma::vec EIDs, const arma::vec &par,
                           const arma::field<arma::uvec> &par_index,
                           const arma::field <arma::vec> &B, const arma::field <arma::vec> &A,
                           const arma::field <arma::vec> &W, const arma::mat &Y, const arma::mat &z,
                           const arma::field<arma::field<arma::mat>> &Dn_omega,
                           const arma::field<arma::field<arma::mat>> &Xn,
                           const arma::vec &bleed_indicator,
                           int states_per_step, int n_cores) {

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

    arma::mat R = arma::reshape(par.elem(par_index(4) - 1), 4, 4);

    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12);

    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);

    arma::vec eids = Y.col(0);
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6);
    // -------------------------------------------------------------------------

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
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);

        arma::vec vec_alpha_i_no_base = A(ii);
        arma::mat alpha_i_no_base = arma::reshape(vec_alpha_i_no_base, 4, 4);

        arma::vec alpha_base = Y_i.col(0) - D_omega_i(0)*vec_omega_i - X_i(0)*vec_beta;

        arma::mat alpha_i(5, 4, arma::fill::zeros);
        alpha_i.row(0) = alpha_base.t(); // NOT gamma_i
        alpha_i.row(1) = alpha_i_no_base.row(0);
        alpha_i.row(2) = alpha_i_no_base.row(1);
        alpha_i.row(3) = alpha_i_no_base.row(2);
        alpha_i.row(4) = alpha_i_no_base.row(3);

        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());

        bool clinic_sub;
        if(clinic_rule < 0) {
            clinic_sub = true;
        } else {
            clinic_sub = false;
        }
        
        arma::vec b_i = B(ii);

        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - states_per_step + 1; k++) {
            
            arma::vec s_i = b_i;

            arma::mat omega_set = get_omega_list(k, n_i, b_i, states_per_step, clinic_sub);
            
            arma::vec prob_omega(omega_set.n_rows, arma::fill::ones);
            prob_omega = (1/arma::accu(prob_omega)) * prob_omega;
            arma::vec ind_omega = arma::linspace(0, omega_set.n_rows-1, omega_set.n_rows);
            arma::vec row_omega = RcppArmadillo::sample(ind_omega, 1, false, prob_omega);

            s_i.rows(k, k+states_per_step-1) = omega_set.row(row_omega(0)).t();

            if(arma::accu(arma::abs(s_i - b_i)) != 0) {
                // Local RBC_rule & clinic_rule ----------------------------
                bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i,
                                            states_per_step, k);

                if(eval_like) {
                    
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
                        
                        if(t == 0) {
                            log_like_s = log_like_s + log(P_init(s_i(t) - 1));
                            log_like_b = log_like_b + log(P_init(b_i(t) - 1));
                        } else {
                            // INDEXING STARTS AT 1 (NOT 0)
                            // Computations for mean of candidate (s_i) ------------
                            arma::vec s_2 = twos_s.subvec(1, t);
                            arma::vec s_3 = threes_s.subvec(1, t);
                            arma::vec s_4 = fours_s.subvec(1, t);
                            arma::vec s_5 = fives_s.subvec(1, t);
                            
                            arma::vec g_p_a_s = alpha_i.row(0).t() +
                                arma::accu(s_2) * alpha_i.row(1).t() +
                                arma::accu(s_3) * alpha_i.row(2).t() +
                                arma::accu(s_4) * alpha_i.row(3).t() +
                                arma::accu(s_5) * alpha_i.row(4).t();
                            
                            arma::vec g_p_a_s_1;
                            if(t == 1) {
                                g_p_a_s_1 = alpha_i.row(0).t();
                            } else {
                                arma::vec s_2_1 = twos_s.subvec(1, t-1);
                                arma::vec s_3_1 = threes_s.subvec(1, t-1);
                                arma::vec s_4_1 = fours_s.subvec(1, t-1);
                                arma::vec s_5_1 = fives_s.subvec(1, t-1);
                                
                                g_p_a_s_1 = alpha_i.row(0).t() +
                                    arma::accu(s_2_1) * alpha_i.row(1).t() +
                                    arma::accu(s_3_1) * alpha_i.row(2).t() +
                                    arma::accu(s_4_1) * alpha_i.row(3).t() +
                                    arma::accu(s_5_1) * alpha_i.row(4).t();
                            }
                            
                            arma::vec nu_s = g_p_a_s + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                            arma::vec nu_s_1 = g_p_a_s_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;
                            
                            arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);
                            
                            // Computations for mean of current (b_i) --------------
                            arma::vec b_2 = twos_b.subvec(1, t);
                            arma::vec b_3 = threes_b.subvec(1, t);
                            arma::vec b_4 = fours_b.subvec(1, t);
                            arma::vec b_5 = fives_b.subvec(1, t);
                            
                            arma::vec g_p_a_b = alpha_i.row(0).t() +
                                arma::accu(b_2) * alpha_i.row(1).t() +
                                arma::accu(b_3) * alpha_i.row(2).t() +
                                arma::accu(b_4) * alpha_i.row(3).t() +
                                arma::accu(b_5) * alpha_i.row(4).t();
                            
                            arma::vec g_p_a_b_1;
                            if(t == 1) {
                                g_p_a_b_1 = alpha_i.row(0).t();
                            } else {
                                arma::vec b_2_1 = twos_b.subvec(1, t-1);
                                arma::vec b_3_1 = threes_b.subvec(1, t-1);
                                arma::vec b_4_1 = fours_b.subvec(1, t-1);
                                arma::vec b_5_1 = fives_b.subvec(1, t-1);
                                
                                g_p_a_b_1 = alpha_i.row(0).t() +
                                    arma::accu(b_2_1) * alpha_i.row(1).t() +
                                    arma::accu(b_3_1) * alpha_i.row(2).t() +
                                    arma::accu(b_4_1) * alpha_i.row(3).t() +
                                    arma::accu(b_5_1) * alpha_i.row(4).t();
                            }
                            
                            arma::vec nu_b = g_p_a_b + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                            arma::vec nu_b_1 = g_p_a_b_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;
                            
                            arma::vec mean_b = nu_b + A_1 * (Y_i.col(t-1) - nu_b_1);
                            
                            arma::vec log_y_pdf_s = dmvnorm(Y_i.col(t).t(), mean_s, R, true);
                            arma::vec log_y_pdf_b = dmvnorm(Y_i.col(t).t(), mean_b, R, true);
                            
                            log_like_s = log_like_s + arma::as_scalar(log_y_pdf_s);
                            log_like_b = log_like_b + arma::as_scalar(log_y_pdf_b);
                            
                            if(t <= k + states_per_step) {
                                arma::mat P_i = get_P_i(t, z_i, zeta);
                                
                                log_like_s = log_like_s + log(P_i(s_i(t-1)-1, s_i(t)-1));
                                log_like_b = log_like_b + log(P_i(b_i(t-1)-1, b_i(t)-1));
                            }    
                        }
                    }
                    
                    double diff_check = log_like_s - log_like_b;
                    
                    double min_log = log(arma::randu(arma::distr_param(0,1)));
                    if(diff_check > min_log){b_i = s_i;} 
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
        twos(0) = 0; threes(0) = 0; fours(0) = 0; fives(0) = 0;

        arma::mat bigB = arma::join_rows(arma::cumsum(twos), arma::cumsum(threes));
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            // Guarantee Dn_alpha(0) = 0
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
Rcpp::List state_almost_gibbs(const arma::vec EIDs, const arma::vec &par,
                              const arma::field<arma::uvec> &par_index,
                              const arma::field <arma::vec> &B, const arma::field <arma::vec> &A,
                              const arma::field <arma::vec> &W, const arma::mat &Y, const arma::mat &z,
                              const arma::field<arma::field<arma::mat>> &Dn_omega,
                              const arma::field<arma::field<arma::mat>> &Xn,
                              const arma::vec &bleed_indicator,
                              int states_per_step, int n_cores) {
    
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

    arma::mat R = arma::reshape(par.elem(par_index(4) - 1), 4, 4);

    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12);

    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))};
    arma::vec P_init = init_logit / arma::accu(init_logit);

    arma::vec eids = Y.col(0);
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6);
    // -------------------------------------------------------------------------

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
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);

        arma::vec vec_alpha_i_no_base = A(ii);
        arma::mat alpha_i_no_base = arma::reshape(vec_alpha_i_no_base, 4, 4);

        arma::vec alpha_base = Y_i.col(0) - D_omega_i(0)*vec_omega_i - X_i(0)*vec_beta;

        arma::mat alpha_i(5, 4, arma::fill::zeros);
        alpha_i.row(0) = alpha_base.t(); // NOT gamma_i
        alpha_i.row(1) = alpha_i_no_base.row(0);
        alpha_i.row(2) = alpha_i_no_base.row(1);
        alpha_i.row(3) = alpha_i_no_base.row(2);
        alpha_i.row(4) = alpha_i_no_base.row(3);

        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());

        bool clinic_sub;
        if(clinic_rule < 0) {
            clinic_sub = true;
        } else {
            clinic_sub = false;
        }

        arma::vec b_i = B(ii);

        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - states_per_step + 1; k++) {

            arma::mat Omega_set = get_omega_list(k, n_i, b_i, states_per_step, clinic_sub);

            if(Omega_set.n_rows > 1) {

                // Compute proposal distribution ("Almost Gibbs") --------------
                arma::vec ss_prob(Omega_set.n_rows, arma::fill::ones);
                arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);

                for(int jj = 0; jj < Omega_set.n_rows; jj++) {

                    double log_like_val = 0;
                    
                    arma::vec ss_jj = b_i;
                    ss_jj.subvec(k, k + states_per_step - 1) = Omega_set.row(jj).t();
                    
                    arma::vec twos_s(ss_jj.n_elem, arma::fill::zeros);
                    arma::vec threes_s = twos_s;
                    arma::vec fours_s = twos_s;
                    arma::vec fives_s = twos_s;
                    twos_s.elem(arma::find(ss_jj == 2)) += 1;
                    threes_s.elem(arma::find(ss_jj == 3)) += 1;
                    fours_s.elem(arma::find(ss_jj == 4)) += 1;
                    fives_s.elem(arma::find(ss_jj == 5)) += 1;

                    for(int t = k; t < k + states_per_step; t++) {
                        
                        if(t == 0) {
                            log_like_val = log_like_val + log(P_init(ss_jj(t) - 1));
                        } else {
                            // INDEXING STARTS AT 1 (NOT 0)
                            // Computations for mean of candidate (s_i) ------------
                            arma::vec s_2 = twos_s.subvec(1, t);
                            arma::vec s_3 = threes_s.subvec(1, t);
                            arma::vec s_4 = fours_s.subvec(1, t);
                            arma::vec s_5 = fives_s.subvec(1, t);
                            
                            arma::vec g_p_a_s = alpha_i.row(0).t() +
                                arma::accu(s_2) * alpha_i.row(1).t() +
                                arma::accu(s_3) * alpha_i.row(2).t() +
                                arma::accu(s_4) * alpha_i.row(3).t() +
                                arma::accu(s_5) * alpha_i.row(4).t();
                            
                            arma::vec g_p_a_s_1;
                            if(t == 1) {
                                g_p_a_s_1 = alpha_i.row(0).t();
                            } else {
                                arma::vec s_2_1 = twos_s.subvec(1, t-1);
                                arma::vec s_3_1 = threes_s.subvec(1, t-1);
                                arma::vec s_4_1 = fours_s.subvec(1, t-1);
                                arma::vec s_5_1 = fives_s.subvec(1, t-1);
                                
                                g_p_a_s_1 = alpha_i.row(0).t() +
                                    arma::accu(s_2_1) * alpha_i.row(1).t() +
                                    arma::accu(s_3_1) * alpha_i.row(2).t() +
                                    arma::accu(s_4_1) * alpha_i.row(3).t() +
                                    arma::accu(s_5_1) * alpha_i.row(4).t();
                            }
                            
                            arma::vec nu_s = g_p_a_s + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                            arma::vec nu_s_1 = g_p_a_s_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;
                            
                            arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);
                            
                            arma::vec log_y_pdf_s = dmvnorm(Y_i.col(t).t(), mean_s, R, true);
                            
                            arma::mat P_i = get_P_i(t, z_i, zeta);
                            
                            log_like_val = log_like_val + arma::as_scalar(log_y_pdf_s) + 
                                                log(P_i(ss_jj(t-1)-1, ss_jj(t)-1));
                        }
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
                s_i.subvec(k, k + states_per_step - 1) = Omega_set.row(row_ind(0)).t();

                if(arma::accu(arma::abs(s_i - b_i)) != 0) {
                    // Local RBC_rule & clinic_rule ----------------------------
                    bool eval_like = rule_check(clinic_rule, rbc_rule, bleed_ind_i, s_i,
                                                states_per_step, k);

                    if(eval_like) {
                        
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
                        
                        for(int t = k + states_per_step; t < n_i; t++) {
                            
                            if(t == 0) {
                                log_like_s = log_like_s + log(P_init(s_i(t) - 1));
                                log_like_b = log_like_b + log(P_init(b_i(t) - 1));
                            } else {
                                // INDEXING STARTS AT 1 (NOT 0)
                                // Computations for mean of candidate (s_i) ------------
                                arma::vec s_2 = twos_s.subvec(1, t);
                                arma::vec s_3 = threes_s.subvec(1, t);
                                arma::vec s_4 = fours_s.subvec(1, t);
                                arma::vec s_5 = fives_s.subvec(1, t);
                                
                                arma::vec g_p_a_s = alpha_i.row(0).t() +
                                    arma::accu(s_2) * alpha_i.row(1).t() +
                                    arma::accu(s_3) * alpha_i.row(2).t() +
                                    arma::accu(s_4) * alpha_i.row(3).t() +
                                    arma::accu(s_5) * alpha_i.row(4).t();
                                
                                arma::vec g_p_a_s_1;
                                if(t == 1) {
                                    g_p_a_s_1 = alpha_i.row(0).t();
                                } else {
                                    arma::vec s_2_1 = twos_s.subvec(1, t-1);
                                    arma::vec s_3_1 = threes_s.subvec(1, t-1);
                                    arma::vec s_4_1 = fours_s.subvec(1, t-1);
                                    arma::vec s_5_1 = fives_s.subvec(1, t-1);
                                    
                                    g_p_a_s_1 = alpha_i.row(0).t() +
                                        arma::accu(s_2_1) * alpha_i.row(1).t() +
                                        arma::accu(s_3_1) * alpha_i.row(2).t() +
                                        arma::accu(s_4_1) * alpha_i.row(3).t() +
                                        arma::accu(s_5_1) * alpha_i.row(4).t();
                                }
                                
                                arma::vec nu_s = g_p_a_s + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                                arma::vec nu_s_1 = g_p_a_s_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;
                                
                                arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);
                                
                                // Computations for mean of current (b_i) --------------
                                arma::vec b_2 = twos_b.subvec(1, t);
                                arma::vec b_3 = threes_b.subvec(1, t);
                                arma::vec b_4 = fours_b.subvec(1, t);
                                arma::vec b_5 = fives_b.subvec(1, t);
                                
                                arma::vec g_p_a_b = alpha_i.row(0).t() +
                                    arma::accu(b_2) * alpha_i.row(1).t() +
                                    arma::accu(b_3) * alpha_i.row(2).t() +
                                    arma::accu(b_4) * alpha_i.row(3).t() +
                                    arma::accu(b_5) * alpha_i.row(4).t();
                                
                                arma::vec g_p_a_b_1;
                                if(t == 1) {
                                    g_p_a_b_1 = alpha_i.row(0).t();
                                } else {
                                    arma::vec b_2_1 = twos_b.subvec(1, t-1);
                                    arma::vec b_3_1 = threes_b.subvec(1, t-1);
                                    arma::vec b_4_1 = fours_b.subvec(1, t-1);
                                    arma::vec b_5_1 = fives_b.subvec(1, t-1);
                                    
                                    g_p_a_b_1 = alpha_i.row(0).t() +
                                        arma::accu(b_2_1) * alpha_i.row(1).t() +
                                        arma::accu(b_3_1) * alpha_i.row(2).t() +
                                        arma::accu(b_4_1) * alpha_i.row(3).t() +
                                        arma::accu(b_5_1) * alpha_i.row(4).t();
                                }
                                
                                arma::vec nu_b = g_p_a_b + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                                arma::vec nu_b_1 = g_p_a_b_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;
                                
                                arma::vec mean_b = nu_b + A_1 * (Y_i.col(t-1) - nu_b_1);
                                
                                arma::vec log_y_pdf_s = dmvnorm(Y_i.col(t).t(), mean_s, R, true);
                                arma::vec log_y_pdf_b = dmvnorm(Y_i.col(t).t(), mean_b, R, true);
                                
                                log_like_s = log_like_s + arma::as_scalar(log_y_pdf_s);
                                log_like_b = log_like_b + arma::as_scalar(log_y_pdf_b);
                                
                                if(t == k + states_per_step) {
                                    arma::mat P_i = get_P_i(t, z_i, zeta);
                                    
                                    log_like_s = log_like_s + log(P_i(s_i(t-1)-1, s_i(t)-1));
                                    log_like_b = log_like_b + log(P_i(b_i(t-1)-1, b_i(t)-1));
                                }    
                            }
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
        twos(0) = 0; threes(0) = 0; fours(0) = 0; fives(0) = 0;

        arma::mat bigB = arma::join_rows(arma::cumsum(twos), arma::cumsum(threes));
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));

        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            // Guarantee Dn_alpha(0) = 0
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
Rcpp::List state_gibbs(const arma::vec EIDs, const arma::vec &par,
                       const arma::field<arma::uvec> &par_index,
                       const arma::field <arma::vec> &B, const arma::field <arma::vec> &A,
                       const arma::field <arma::vec> &W, const arma::mat &Y, const arma::mat &z,
                       const arma::field<arma::field<arma::mat>> &Dn_omega,
                       const arma::field<arma::field<arma::mat>> &Xn,
                       const arma::vec &bleed_indicator,
                       int states_per_step, int n_cores) {
        
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
    
        arma::mat R = arma::reshape(par.elem(par_index(4) - 1), 4, 4);
    
        arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
        arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12);
    
        arma::vec vec_init_content = par.elem(par_index(6) - 1);
        arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                                exp(vec_init_content(2)), exp(vec_init_content(3))};
        arma::vec P_init = init_logit / arma::accu(init_logit);
    
        arma::vec eids = Y.col(0);
        arma::vec rbc_rule_vec = Y.col(5);
        arma::vec clinic_rule_vec = Y.col(6);
        // -------------------------------------------------------------------------
    
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
            arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
    
            arma::vec vec_alpha_i_no_base = A(ii);
            arma::mat alpha_i_no_base = arma::reshape(vec_alpha_i_no_base, 4, 4);
    
            arma::vec alpha_base = Y_i.col(0) - D_omega_i(0)*vec_omega_i - X_i(0)*vec_beta;
    
            arma::mat alpha_i(5, 4, arma::fill::zeros);
            alpha_i.row(0) = alpha_base.t(); // NOT gamma_i
            alpha_i.row(1) = alpha_i_no_base.row(0);
            alpha_i.row(2) = alpha_i_no_base.row(1);
            alpha_i.row(3) = alpha_i_no_base.row(2);
            alpha_i.row(4) = alpha_i_no_base.row(3);
    
            int rbc_rule = rbc_rule_vec(sub_ind.min());
            int clinic_rule = clinic_rule_vec(sub_ind.min());
    
            bool clinic_sub;
            if(clinic_rule < 0) {
                clinic_sub = true;
            } else {
                clinic_sub = false;
            }
    
            arma::vec b_i = B(ii);
    
            // Looping through subject state space ---------------------------------
            for (int k = 0; k < n_i - states_per_step + 1; k++) {
    
                arma::mat Omega_set = get_omega_list(k, n_i, b_i, states_per_step, clinic_sub);
    
                if(Omega_set.n_rows > 1) {
    
                    // Compute proposal distribution ("Gibbs") ---------------------
                    arma::vec ss_prob(Omega_set.n_rows, arma::fill::ones);
                    arma::vec ss_ind = arma::linspace(0, Omega_set.n_rows-1, Omega_set.n_rows);
    
                    for(int jj = 0; jj < Omega_set.n_rows; jj++) {
    
                        double log_like_val = 0;
                        
                        arma::vec ss_jj = b_i;
                        ss_jj.subvec(k, k + states_per_step - 1) = Omega_set.row(jj).t();
                        
                        arma::vec twos_s(ss_jj.n_elem, arma::fill::zeros);
                        arma::vec threes_s = twos_s;
                        arma::vec fours_s = twos_s;
                        arma::vec fives_s = twos_s;
                        twos_s.elem(arma::find(ss_jj == 2)) += 1;
                        threes_s.elem(arma::find(ss_jj == 3)) += 1;
                        fours_s.elem(arma::find(ss_jj == 4)) += 1;
                        fives_s.elem(arma::find(ss_jj == 5)) += 1;
                        
                        for(int t = k; t < n_i; t++) {
    
                            if(t == 0) {
                                log_like_val = log_like_val + log(P_init(ss_jj(t) - 1));
                            } else {
                                // INDEXING STARTS AT 1 (NOT 0)
                                // Computations for mean of candidate (s_i) ------------
                                arma::vec s_2 = twos_s.subvec(1, t);
                                arma::vec s_3 = threes_s.subvec(1, t);
                                arma::vec s_4 = fours_s.subvec(1, t);
                                arma::vec s_5 = fives_s.subvec(1, t);
    
                                arma::vec g_p_a_s = alpha_i.row(0).t() +
                                    arma::accu(s_2) * alpha_i.row(1).t() +
                                    arma::accu(s_3) * alpha_i.row(2).t() +
                                    arma::accu(s_4) * alpha_i.row(3).t() +
                                    arma::accu(s_5) * alpha_i.row(4).t();
    
                                arma::vec g_p_a_s_1;
                                if(t == 1) {
                                    g_p_a_s_1 = alpha_i.row(0).t();
                                } else {
                                    arma::vec s_2_1 = twos_s.subvec(1, t-1);
                                    arma::vec s_3_1 = threes_s.subvec(1, t-1);
                                    arma::vec s_4_1 = fours_s.subvec(1, t-1);
                                    arma::vec s_5_1 = fives_s.subvec(1, t-1);
    
                                    g_p_a_s_1 = alpha_i.row(0).t() +
                                        arma::accu(s_2_1) * alpha_i.row(1).t() +
                                        arma::accu(s_3_1) * alpha_i.row(2).t() +
                                        arma::accu(s_4_1) * alpha_i.row(3).t() +
                                        arma::accu(s_5_1) * alpha_i.row(4).t();
                                }
    
                                arma::vec nu_s = g_p_a_s + D_omega_i(t)*vec_omega_i + X_i(t)*vec_beta;
                                arma::vec nu_s_1 = g_p_a_s_1 + D_omega_i(t-1)*vec_omega_i + X_i(t-1)*vec_beta;
    
                                arma::vec mean_s = nu_s + A_1 * (Y_i.col(t-1) - nu_s_1);
    
                                arma::vec log_y_pdf_s = dmvnorm(Y_i.col(t).t(), mean_s, R, true);
    
                                arma::mat P_i = get_P_i(t, z_i, zeta);
    
                                log_like_val = log_like_val + arma::as_scalar(log_y_pdf_s) +
                                    log(P_i(ss_jj(t-1)-1, ss_jj(t)-1));
                            }
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
                    s_i.subvec(k, k + states_per_step - 1) = Omega_set.row(row_ind(0)).t();
    
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
            twos(0) = 0; threes(0) = 0; fours(0) = 0; fives(0) = 0;
    
            arma::mat bigB = arma::join_rows(arma::cumsum(twos), arma::cumsum(threes));
            bigB = arma::join_rows(bigB, arma::cumsum(fours));
            bigB = arma::join_rows(bigB, arma::cumsum(fives));
    
            arma::mat I = arma::eye(4,4);
            for(int jj = 0; jj < n_i; jj++) {
                // Guarantee Dn_alpha(0) = 0
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
void get_dimension(){
    
    int N = 5;
    
    Rcpp::Rcout << "Case (c) Full" << std::endl;
    for(int w=0; w < N; w++) {
        Rcpp::Rcout << "() -> () -> " << w+1 << " : " << Omega_List_GLOBAL_multi(0)(w).n_rows << " combos" << std::endl;
    }

    Rcpp::Rcout << "Case (b) Full" << std::endl;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            Rcpp::Rcout << i+1 << " --> " << j+1 << " : " << Omega_List_GLOBAL_multi(1)(i, j).n_rows << " combos" << std::endl;
        }
    }

    Rcpp::Rcout << "Case (a) Full" << std::endl;
    for(int w=0; w < N; w++) {
        Rcpp::Rcout << w + 1 << " -> () -> () : " << Omega_List_GLOBAL_multi(2)(w).n_rows << " combos" << std::endl;
    }
}
