#ifndef RLS_MAG_H
#define RLS_MAG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ==========================================================
// ==                  CONSTANTS & MACROS                  ==
// ==========================================================

#define N 9 // 상태 벡터의 파라미터 수
#define L 5 // 초기 데이터 길이 계산을 위한 승수
#define INITIAL_DATA_LEN (N * L) // 초기화에 필요한 데이터 총량

// 1차원 배열을 2D 행렬처럼 사용하기 위한 매크로
#define IDX(row, col, width) ((row) * (width) + (col))

// ==========================================================
// ==                  DATA STRUCTURES                     ==
// ==========================================================

// 센서 데이터 타입을 double로 정의 (float를 사용하려면 float로 변경)
typedef double vqf_real_t;

// RLS 필터의 모든 상태를 관리하는 구조체
typedef struct {
    double* P_cov;      // 공분산 행렬 (N x N)
    double* est_param;  // 추정 파라미터 벡터 (N x 1)
    double* P_init;     // 초기 입력 데이터 버퍼 (INITIAL_DATA_LEN x N)
    double* Y_init;     // 초기 출력 데이터 버퍼 (INITIAL_DATA_LEN x 1)
    int data_count;     // 현재까지 수집된 데이터 개수
    int is_initialized; // RLS 필터 초기화 완료 여부 플래그
} RLS_State;

// ==========================================================
// ==                FUNCTION PROTOTYPES                   ==
// ==========================================================

// --- RLS Filter Management ---
RLS_State* rls_init();
int rls_update(RLS_State* state, const vqf_real_t mag[3]);
void rls_destroy(RLS_State* state);

// --- Core Algorithm Functions ---
void create_input_vector(vqf_real_t mx, vqf_real_t my, vqf_real_t mz, double* input_vec);
void extract_iron_parameters(const double* est_param, double* hard_iron, double* soft_iron, double* v_matrix);
void mag_raw_to_cal(double* hard_iron, double* soft_iron, double* v_matrix, vqf_real_t* mag);

// --- Linear Algebra Utility Functions ---
void matrix_print(const char* title, const double* A, int rows, int cols);
void matrix_multiply(const double* A, const double* B, double* C, int r1, int c1, int c2);
void matrix_transpose(const double* A, double* A_T, int rows, int cols);
int matrix_inverse(const double* A, double* A_inv, int n);
void eigen_decomposition_3x3_symmetric(const double* A, double* eigenvalues, double* eigenvectors);

#endif // RLS_MAG_H