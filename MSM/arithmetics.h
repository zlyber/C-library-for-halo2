#ifndef __ARITHMETICS_H__
#define __ARITHMETICS_H__

#include "curve.h"
#include "math.h"
#include "fields.h"
#include "fq-bn256.h"
#include "fr.h"

//get the target bit string's in decimalism
uint64_t get_at(int segment, int c, uint8_t* bytes);

//judge whether the bucket is null now
bool is_null(uint64_t* x);

//MSM
void multiexp_serial(uint64_t** coeffs, uint64_t** bases_x, uint64_t** bases_y, uint64_t** bases_z, uint64_t* acc_x, uint64_t* acc_y, uint64_t* acc_z, uint64_t len);


#endif