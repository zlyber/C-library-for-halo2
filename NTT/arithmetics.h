#ifndef __ARITHMETICS_H__
#define __ARITHMETICS_H__

#include <cstdint>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include <sys/timeb.h>
#include "fr-bn256.h"


// using namespace std;

void u64_to_u64(uint64_t* u641, uint64_t* u642);
void swap(uint64_t* a, uint64_t* b) ;
uint64_t reverse_bits(uint64_t operand,int bit_count);
void group_scale(uint64_t* self, uint64_t* rhs, uint64_t*res);
void group_add(uint64_t* self, uint64_t* rhs, uint64_t* res);
void group_sub(uint64_t* self, uint64_t* rhs, uint64_t* res);
void scalar_one(uint64_t* result);

// precompute w^n
void precompute_omega(uint64_t* root_of_unity, int k, uint64_t* omega, int j = 3);
// precompute twiddle factors
void precompute_twiddles(uint64_t** twiddles, uint64_t *omega, uint64_t n);
//recursive version
void recursive_butterfly(uint64_t** vector,uint64_t N, uint64_t twiddle_chunk, uint64_t** twiddles);
//k=log_n
void fft(uint64_t** vector, uint64_t* omega, int k);



#endif