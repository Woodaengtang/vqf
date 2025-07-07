#include "rls_mag.h"

// 내부 초기화 함수 선언
void perform_batch_lsm(RLS_State* state);

int main() {

    RLS_State* rls_state = rls_init();
    if (!rls_state) return 1;

    while (/**mag raw값을 calibration으로 사용하는 동안~ */) {
        rls_update(rls_state, mag);
    }
    double final_hard[3], final_soft[9], final_v[9];
    extract_iron_parameters(rls_state->est_param, final_hard, final_soft, final_v);
    rls_destroy(rls_state);
    return 0;
}

// =============================================================================
// ========================= RLS Filter Management =============================
// =============================================================================

RLS_State* rls_init() {
    RLS_State* state = (RLS_State*)malloc(sizeof(RLS_State));
    if (!state) return NULL;
    state->P_cov = (double*)malloc(N * N * sizeof(double));
    state->est_param = (double*)calloc(N, sizeof(double));
    state->P_init = (double*)malloc(INITIAL_DATA_LEN * N * sizeof(double));
    state->Y_init = (double*)malloc(INITIAL_DATA_LEN * sizeof(double));
    if (!state->P_cov || !state->est_param || !state->P_init || !state->Y_init) {
        free(state->P_cov); free(state->est_param); free(state->P_init); free(state->Y_init); free(state);
        return NULL;
    }
    state->data_count = 0;
    state->is_initialized = 0; // false
    return state;
}

int rls_update(RLS_State* state, const vqf_real_t mag[3]) {
    if (!state) return -1;

    if (!state->is_initialized) {
        double input[N];
        create_input_vector(mag[0], mag[1], mag[2], input);
        memcpy(&state->P_init[IDX(state->data_count, 0, N)], input, N * sizeof(double));
        state->Y_init[state->data_count] = -1.0;
        state->data_count++;
        if (state->data_count >= INITIAL_DATA_LEN) {
            perform_batch_lsm(state);
            state->is_initialized = 1;
            return 1; // 초기화 완료
        }
        return 0; // 버퍼링 중
    }

    double input[N];
    create_input_vector(mag[0], mag[1], mag[2], input);
    double e_k = -1.0;
    for (int j = 0; j < N; ++j) e_k -= input[j] * state->est_param[j];
    double P_input_T[N];
    matrix_multiply(state->P_cov, input, P_input_T, N, N, 1);
    double input_P_input_T = 0;
    for (int j = 0; j < N; ++j) input_P_input_T += input[j] * P_input_T[j];
    double denominator = 1.0 + input_P_input_T;
    double numerator_term[N * N];
    matrix_multiply(P_input_T, input, numerator_term, N, 1, N);
    for (int j = 0; j < N * N; ++j) {
        state->P_cov[j] -= numerator_term[j] / denominator;
    }
    for (int j = 0; j < N; ++j) {
        state->est_param[j] += P_input_T[j] * e_k;
    }
    return 2; // 업데이트됨
}

void rls_destroy(RLS_State* state) {
    if (state) {
        free(state->P_cov);
        free(state->est_param);
        free(state->P_init); // is_initialized가 false일 경우를 대비해 호출
        free(state->Y_init); // NULL 포인터에 free를 호출하는 것은 안전함
        free(state);
    }
}

// =============================================================================
// ========================= INTERNAL HELPER FUNCTIONS =========================
// =============================================================================

void perform_batch_lsm(RLS_State* state) {
    double* P_T = (double*)malloc(N * INITIAL_DATA_LEN * sizeof(double));
    double* P_T_P = (double*)malloc(N * N * sizeof(double));
    double* P_T_P_inv = (double*)malloc(N * N * sizeof(double));
    double* P_T_Y = (double*)malloc(N * sizeof(double));

    matrix_transpose(state->P_init, P_T, INITIAL_DATA_LEN, N);
    matrix_multiply(P_T, state->P_init, P_T_P, N, INITIAL_DATA_LEN, N);
    if (matrix_inverse(P_T_P, P_T_P_inv, N)) {
        memcpy(state->P_cov, P_T_P_inv, N * N * sizeof(double));
        matrix_multiply(P_T, state->Y_init, P_T_Y, N, INITIAL_DATA_LEN, 1);
        matrix_multiply(state->P_cov, P_T_Y, state->est_param, N, N, 1);
    } else {
        fprintf(stderr, "Error: Initial matrix is singular. Calibration may fail.\n");
    }

    free(P_T); free(P_T_P); free(P_T_P_inv); free(P_T_Y);
    
    // 핵심: 사용이 끝난 초기화 버퍼 메모리 해제
    free(state->P_init);
    free(state->Y_init);
    state->P_init = NULL;
    state->Y_init = NULL;
}

void create_input_vector(vqf_real_t mx, vqf_real_t my, vqf_real_t mz, double* input_vec) {
    input_vec[0] = mx * mx;
    input_vec[1] = my * my;
    input_vec[2] = mz * mz;
    input_vec[3] = 2.0 * mx * my;
    input_vec[4] = 2.0 * mx * mz;
    input_vec[5] = 2.0 * my * mz;
    input_vec[6] = mx;
    input_vec[7] = my;
    input_vec[8] = mz;
}


void extract_iron_parameters(const double* est_param, double* hard_iron, double* soft_iron, double* v_matrix) {
    double A[9];
    A[IDX(0,0,3)] = est_param[0]; A[IDX(0,1,3)] = est_param[3]; A[IDX(0,2,3)] = est_param[4];
    A[IDX(1,0,3)] = est_param[3]; A[IDX(1,1,3)] = est_param[1]; A[IDX(1,2,3)] = est_param[5];
    A[IDX(2,0,3)] = est_param[4]; A[IDX(2,1,3)] = est_param[5]; A[IDX(2,2,3)] = est_param[2];
    const double* b = &est_param[6];
    double D_diag[3];
    eigen_decomposition_3x3_symmetric(A, D_diag, v_matrix);
    double b_new[3];
    matrix_multiply(b, v_matrix, b_new, 1, 3, 3);
    double temp_sum = 0.0;
    for (int i = 0; i < 3; ++i) {
        if (fabs(D_diag[i]) > 1e-9) { temp_sum += (b_new[i] * b_new[i]) / (4 * D_diag[i]); }
    }
    double scale_q = sqrt(fmax(0, temp_sum - 1.0));
    double scale_matrix[9] = {0};
    if (fabs(scale_q) > 1e-9) {
        for (int i=0; i<3; ++i) { scale_matrix[IDX(i,i,3)] = sqrt(fmax(0, D_diag[i])) / scale_q; }
    }
    for (int i=0; i<3; ++i) {
        hard_iron[i] = (fabs(D_diag[i]) > 1e-9) ? -b_new[i] / (2 * D_diag[i]) : 0;
    }
    double hard_temp[3];
    memcpy(hard_temp, hard_iron, 3*sizeof(double));
    matrix_multiply(hard_temp, v_matrix, hard_iron, 1, 3, 3);
    matrix_multiply(v_matrix, scale_matrix, soft_iron, 3, 3, 3);
}

void mag_raw_to_cal(double* hard_iron, double* soft_iron, double* v_matrix, vqf_real_t* mag){
    if (!hard_iron || !soft_iron || !v_matrix || !mag) return;

    double v_T[9];
    matrix_transpose(v_matrix, v_T, 3, 3);
    vqf_real_t temp_mag[3];
    temp_mag[0] = mag[0] * v_T[IDX(0, 0, 3)] + mag[1] * v_T[IDX(0, 1, 3)] + mag[2] * v_T[IDX(0, 2, 3)];
    temp_mag[1] = mag[0] * v_T[IDX(1, 0, 3)] + mag[1] * v_T[IDX(1, 1, 3)] + mag[2] * v_T[IDX(1, 2, 3)];
    temp_mag[2] = mag[0] * v_T[IDX(2, 0, 3)] + mag[1] * v_T[IDX(2, 1, 3)] + mag[2] * v_T[IDX(2, 2, 3)];

    temp_mag[0] += hard_iron[0];
    temp_mag[1] += hard_iron[1];
    temp_mag[2] += hard_iron[2];

    mag[0] = temp_mag[0] * soft_iron[IDX(0, 0, 3)] + temp_mag[1] * soft_iron[IDX(0, 1, 3)] + temp_mag[2] * soft_iron[IDX(0, 2, 3)];
    mag[1] = temp_mag[0] * soft_iron[IDX(1, 0, 3)] + temp_mag[1] * soft_iron[IDX(1, 1, 3)] + temp_mag[2] * soft_iron[IDX(1, 2, 3)];
    mag[2] = temp_mag[0] * soft_iron[IDX(2, 0, 3)] + temp_mag[1] * soft_iron[IDX(2, 1, 3)] + temp_mag[2] * soft_iron[IDX(2, 2, 3)];
}

void matrix_print(const char* title, const double* A, int rows, int cols) {
    printf("%s\n", title);
    for (int i = 0; i < rows; ++i) {
        printf("  [ ");
        for (int j = 0; j < cols; ++j) { printf("%9.6f ", A[IDX(i, j, cols)]); }
        printf("]\n");
    }
    printf("\n");
}

void matrix_multiply(const double* A, const double* B, double* C, int r1, int c1, int c2) {
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c2; ++j) {
            C[IDX(i, j, c2)] = 0;
            for (int k = 0; k < c1; ++k) { C[IDX(i, j, c2)] += A[IDX(i, k, c1)] * B[IDX(k, j, c2)]; }
        }
    }
}

void matrix_transpose(const double* A, double* A_T, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) { A_T[IDX(j, i, rows)] = A[IDX(i, j, cols)]; }
    }
}

int matrix_inverse(const double* A, double* A_inv, int n) {
    double* temp = (double*)malloc(n * n * sizeof(double));
    if (!temp) return 0;
    memcpy(temp, A, n * n * sizeof(double));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) { 
            A_inv[IDX(i, j, n)] = (i == j) ? 1.0 : 0.0; 
        }
    }
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(temp[IDX(j, i, n)]) > fabs(temp[IDX(pivot, i, n)])) { pivot = j; }
        }
        for (int k = 0; k < n; ++k) {
            double t = temp[IDX(i, k, n)]; 
            temp[IDX(i, k, n)] = temp[IDX(pivot, k, n)]; 
            temp[IDX(pivot, k, n)] = t;
            t = A_inv[IDX(i, k, n)];
            A_inv[IDX(i, k, n)] = A_inv[IDX(pivot, k, n)]; 
            A_inv[IDX(pivot, k, n)] = t;
        }
        if (fabs(temp[IDX(i, i, n)]) < 1e-12) { free(temp); return 0; }
        double div = temp[IDX(i, i, n)];
        for (int j = 0; j < n; ++j) {
            temp[IDX(i, j, n)] /= div; 
            A_inv[IDX(i, j, n)] /= div;
        }
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double mult = temp[IDX(j, i, n)];
                for (int k = 0; k < n; ++k) {
                    temp[IDX(j, k, n)] -= mult * temp[IDX(i, k, n)];
                    A_inv[IDX(j, k, n)] -= mult * A_inv[IDX(i, k, n)];
                }
            }
        }
    }
    free(temp);
    return 1;
}

void eigen_decomposition_3x3_symmetric(const double* A, double* eigenvalues, double* eigenvectors) {
    double a = A[IDX(0,0,3)], b = A[IDX(0,1,3)], c = A[IDX(0,2,3)];
    double d = A[IDX(1,1,3)], e = A[IDX(1,2,3)];
    double f = A[IDX(2,2,3)];
    double trA = a + d + f;
    double detA = a*(d*f - e*e) - b*(b*f - c*e) + c*(b*e - c*d);
    double term2 = (a*d - b*b) + (a*f - c*c) + (d*f - e*e);
    double p1 = term2;
    double q1 = -detA;
    double p = p1 - (trA*trA)/3.0;
    double q = q1 - (trA*p1)/3.0 + (2.0*trA*trA*trA)/27.0;
    double R = sqrt(fmax(0,-p*p*p/27.0));
    double phi = acos(fmin(1.0, fmax(-1.0, -q / (2.0 * R))));
    eigenvalues[0] = 2.0 * R * cos(phi/3.0) + trA/3.0;
    eigenvalues[1] = 2.0 * R * cos((phi + 2*M_PI)/3.0) + trA/3.0;
    eigenvalues[2] = 2.0 * R * cos((phi + 4*M_PI)/3.0) + trA/3.0;
    for (int i = 0; i < 3; ++i) {
        double lambda = eigenvalues[i];
        double B[9];
        memcpy(B, A, 9*sizeof(double));
        B[IDX(0,0,3)] -= lambda; B[IDX(1,1,3)] -= lambda; B[IDX(2,2,3)] -= lambda;
        double v[3];
        v[0] = B[IDX(0,1,3)] * B[IDX(1,2,3)] - B[IDX(0,2,3)] * B[IDX(1,1,3)];
        v[1] = B[IDX(0,2,3)] * B[IDX(1,0,3)] - B[IDX(0,0,3)] * B[IDX(1,2,3)];
        v[2] = B[IDX(0,0,3)] * B[IDX(1,1,3)] - B[IDX(0,1,3)] * B[IDX(1,0,3)];
        double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        if (norm < 1e-9) {
             v[0] = B[IDX(1,1,3)] * B[IDX(2,2,3)] - B[IDX(1,2,3)] * B[IDX(2,1,3)];
             v[1] = B[IDX(1,2,3)] * B[IDX(2,0,3)] - B[IDX(1,0,3)] * B[IDX(2,2,3)];
             v[2] = B[IDX(1,0,3)] * B[IDX(2,1,3)] - B[IDX(1,1,3)] * B[IDX(2,0,3)];
             norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        }
        if (norm > 1e-9) {
            eigenvectors[IDX(0,i,3)] = v[0] / norm;
            eigenvectors[IDX(1,i,3)] = v[1] / norm;
            eigenvectors[IDX(2,i,3)] = v[2] / norm;
        }
    }
}