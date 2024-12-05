#ifndef __CURVE_H__
#define __CURVE_H__
#include <stdint.h>
#include <stdio.h>
#include "fq-bn256.h"
#include "fr-bn256.h"
#include "fq2-bn256.h"

//judge the affine coordinate is the curve's identity?
bool is_identity(uint64_t* x,uint64_t* y);

//judge the projective coordinate is the curve's identity?
bool is_identity_project(uint64_t* x,uint64_t* y,uint64_t* z);

//affine to projection
void to_curve(uint64_t* x,uint64_t* y,uint64_t* z);

//projection to affine
void to_affine(uint64_t* x,uint64_t* y,uint64_t* z);

void batch_normalize(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* affine_x,uint64_t* affine_y);

//point add in affine form : (x,y)+(X,Y) and convert to projection
void affine_add_affine(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* X,uint64_t* Y,uint64_t* Z);

//projective point add affine point : (x,y,z)+(X,Y) and convert to projection
void project_add_affine(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* X,uint64_t* Y,uint64_t* Z);

//point add in Jacobian form : (x,y,z)+(X,Y,Z)
void project_add_project(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* X,uint64_t* Y,uint64_t* Z);

//point double in Jacobian form: 2(x,y,z)
void point_double(uint64_t* x,uint64_t* y,uint64_t* z);

//return the identity of the curve in projective form
void identity(uint64_t* x,uint64_t* y,uint64_t* z);

//return the identity of the curve in affine form
void affine_identity(uint64_t* x,uint64_t* y);

/// Returns a fixed generator of the prime-order subgroup.
void G1_affine_generator(uint64_t* x,uint64_t* y);

////double-and-add implementation of point multiplication
void PMULT(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* scalar,uint64_t* acc_x,uint64_t* acc_y,uint64_t* acc_z);

// PMULT where P is in affine mode
void Affine_PMULT(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* scalar,uint64_t* acc_x,uint64_t* acc_y,uint64_t* acc_z);

//get the G2 generator in affine mode
void G2_affine_generator(uint64_t* x,uint64_t* y);

// G2 affine to projective 
void G2_to_curve(uint64_t* x,uint64_t* y,uint64_t* z); 

//G2 projection to affine
void G2_to_affine(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* affine_x,uint64_t* affine_y);

//G2 point double in Jacobian form: 2(x,y,z)
void G2_point_double(uint64_t* x,uint64_t* y,uint64_t* z);

//projective G2 point add affine point : (x,y,z)+(X,Y) and convert to projection
void G2_project_add_affine(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* X,uint64_t* Y,uint64_t* Z);

// PMULT where P is a G2 point and in affine mode
void G2_Affine_PMULT(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* scalar,uint64_t* acc_x,uint64_t* acc_y,uint64_t* acc_z);

#endif