#ifndef __CURVE_H__
#define __CURVE_H__
#include <stdint.h>
#include <stdio.h>
#include "fq-bn256.h"
#include "fields.h"
#include "fr.h"

//judge the affine coordinate is the curve's identity?
bool is_identity(uint64_t* x,uint64_t* y);

//judge the projective coordinate is the curve's identity?
bool is_identity_project(uint64_t* x,uint64_t* y,uint64_t* z);

//affine to projective
void to_curve(uint64_t* x,uint64_t* y,uint64_t* z);

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

#endif