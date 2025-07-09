#include <stdio.h>
#include <string.h>
#include <math.h>

#define N 5     // Number of parameters
#define LEN 100 // Data length

// mat : matrix
// vec : vector
// out : output vector

void mat_vec_mul(vqf_real_t *mat, vqf_real_t *vec, vqf_real_t *out, int n)
{
    for (int i = 0; i < n; i++)
    {
        out[i] = 0;
        for (int j = 0; j < n; j++)
            out[i] += mat[i * n + j] * vec[j];
    }
}

void vec_outer(vqf_real_t *a, vqf_real_t *b, vqf_real_t *out, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            out[i * n + j] = a[i] * b[j];
}

void mat_add(vqf_real_t *a, vqf_real_t *b, vqf_real_t *out, int n)
{
    for (int i = 0; i < n * n; i++)
        out[i] = a[i] + b[i];
}

void mat_sub(vqf_real_t *a, vqf_real_t *b, vqf_real_t *out, int n)
{
    for (int i = 0; i < n * n; i++)
        out[i] = a[i] - b[i];
}

// Inverse of a square matrix using Gauss-Jordan elimination
void mat_inv(vqf_real_t *A, vqf_real_t *A_inv, int n)
{
    vqf_real_t aug[n * 2 * n];
    // Form [A | I]
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            aug[i * 2 * n + j] = A[i * n + j];
            aug[i * 2 * n + (j + n)] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Gauss-Jordan elimination
    for (int i = 0; i < n; i++)
    {
        vqf_real_t piv = aug[i * 2 * n + i];
        if (piv == 0)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (aug[j * 2 * n + i] != 0)
                {
                    for (int k = 0; k < 2 * n; k++)
                    {
                        vqf_real_t tmp = aug[i * 2 * n + k];
                        aug[i * 2 * n + k] = aug[j * 2 * n + k];
                        aug[j * 2 * n + k] = tmp;
                    }
                    piv = aug[i * 2 * n + i];
                    break;
                }
            }
        }
        // Normalize row
        for (int k = 0; k < 2 * n; k++)
            aug[i * 2 * n + k] /= piv;
        // Eliminate other rows
        for (int j = 0; j < n; j++)
        {
            if (j == i)
                continue;
            vqf_real_t f = aug[j * 2 * n + i];
            for (int k = 0; k < 2 * n; k++)
                aug[j * 2 * n + k] -= f * aug[i * 2 * n + k];
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A_inv[i * n + j] = aug[i * 2 * n + (j + n)];
}

void mat_copy(vqf_real_t *src, vqf_real_t *dst, int n)
{
    for (int i = 0; i < n * n; i++)
        dst[i] = src[i];
}

void eye(vqf_real_t *mat, int n)
{
    memset(mat, 0, sizeof(vqf_real_t) * n * n);
    for (int i = 0; i < n; i++)
        mat[i * n + i] = 1.0;
}

void get_cali_param(const vqf_real_t *est_param, vqf_real_t *cali_data)
{
    vqf_real_t b = est_param[0], c = est_param[1], d = est_param[2], e = est_param[3], f = est_param[4];
    vqf_real_t denom = b * b - 4.0 * c;
    vqf_real_t delta = 0.5 * atan2(b, 1 - c);
    vqf_real_t cx = (2.0 * c * d - b * e) / denom;
    vqf_real_t cy = (2.0 * b * e - b * d) / denom;
    vqf_real_t cosd = cos(delta), sind = sin(delta);
    vqf_real_t numer = cx * cx + b * cx * cy + c * cy * cy - f;
    vqf_real_t w = sqrt(numer / (cosd * cosd + b * cosd * sind + c * sind * sind));
    vqf_real_t h = sqrt(numer / (sind * sind + b * cosd * sind + c * cosd * cosd));
    cali_data[0] = delta;
    cali_data[1] = cx;
    cali_data[2] = cy;
    cali_data[3] = h / w;
}

int main()
{
    // vqf_real_t raw_data_mx[LEN];
    // vqf_real_t raw_data_my[LEN];
    // 실제로 raw data가 와야 함 (아래는 예시 데이터)
    // for (int i = 0; i < LEN; i++){
    //     raw_data_mx[i] = i * 0.1;
    //     raw_data_my[i] = i * 0.2;
    // }

    vqf_real_t P[N * N] = {0};
    vqf_real_t Y[N] = {0};
    vqf_real_t est_param[N] = {0};
    int count = 0;
    int min_len = N;
    int data_len = LEN;
    int init_flag = 1;

    while (count < data_len)
    {
        // 실제로 raw_data_mx, raw_data_my에 실제 mx, my 데이터가 와야 함 (아래는 예시)
        // vqf_real_t new_mag_x = raw_data_mx[count];
        // vqf_real_t new_mag_y = raw_data_my[count];
        vqf_real_t input[N] = {
            new_mag_x * new_mag_y,
            new_mag_y * new_mag_y,
            new_mag_x,
            new_mag_y,
            1.0};
        vqf_real_t output = -new_mag_x * new_mag_x;

        if (count >= min_len)
        {
            if (init_flag)
            {
                // est_param = inv(P' * P) * P' * Y;
                vqf_real_t PtP[N * N] = {0};
                vqf_real_t Pt[N * N] = {0};
                // Compute P'
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++)
                        Pt[i * N + j] = P[j * N + i];
                // PtP = P' * P
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++)
                    {
                        PtP[i * N + j] = 0;
                        for (int k = 0; k < N; k++)
                            PtP[i * N + j] += Pt[i * N + k] * P[k * N + j];
                    }
                // inv(PtP)
                vqf_real_t inv_PtP[N * N] = {0};
                mat_inv(PtP, inv_PtP, N);
                // PtY = P' * Y
                vqf_real_t PtY[N] = {0};
                for (int i = 0; i < N; i++)
                {
                    PtY[i] = 0;
                    for (int j = 0; j < N; j++)
                        PtY[i] += Pt[i * N + j] * Y[j];
                }
                // est_param = inv(PtP) * PtY
                for (int i = 0; i < N; i++)
                {
                    est_param[i] = 0;
                    for (int j = 0; j < N; j++)
                        est_param[i] += inv_PtP[i * N + j] * PtY[j];
                }
                init_flag = 0;
            }
            // e_k = output - input * est_param;
            vqf_real_t pred = 0;
            for (int i = 0; i < N; i++)
                pred += input[i] * est_param[i];
            vqf_real_t e_k = output - pred;

            // P = P - (P * input') / (1 + input * P * input') * (input * P);
            vqf_real_t Pin[N];
            mat_vec_mul(P, input, Pin, N); // P * input'
            vqf_real_t denom = 1.0;
            for (int i = 0; i < N; i++)
                denom += input[i] * Pin[i];

            vqf_real_t outer[N * N];
            vec_outer(Pin, Pin, outer, N);

            vqf_real_t update[N * N];
            for (int i = 0; i < N * N; i++)
                update[i] = outer[i] / denom;

            vqf_real_t P_new[N * N];
            mat_sub(P, update, P_new, N);
            mat_copy(P_new, P, N);

            // est_param = est_param + P * input' * e_k;
            mat_vec_mul(P, input, Pin, N);
            for (int i = 0; i < N; i++)
                est_param[i] += Pin[i] * e_k;

            vqf_real_t cali_data[4] = {0};
            get_cali_param(est_param, cali_data);
        }
        else
        {
            for (int i = 0; i < N; i++)
                P[count * N + i] = input[i];
            Y[count] = output;
        }
        count++;
    }

    return 0;
}
