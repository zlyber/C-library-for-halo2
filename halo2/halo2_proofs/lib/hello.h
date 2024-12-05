#include <stdint.h>

void NTT(uint64_t vector[][4], uint64_t* omega, uint32_t k);

void MSM(uint64_t** coeffs, uint64_t** bases_x, uint64_t** bases_y, uint64_t** bases_z, uint64_t* acc_x, uint64_t* acc_y, uint64_t* acc_z, uint64_t len);

void SETUP(uint32_t k, uint64_t** g_x, uint64_t** g_y, uint64_t** g_lagrange_x, uint64_t** g_lagrange_y, uint64_t g2_x[8], uint64_t g2_y[8], uint64_t s_g2_x[8], uint64_t s_g2_y[8]);