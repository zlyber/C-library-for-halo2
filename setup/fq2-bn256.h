#include"fields.h"

//redefine field's operators in extend Field Fq2

void fq2_one(uint64_t* self,const uint64_t* R);

void fq2_zero(uint64_t* self);

void fq2_sub(uint64_t* self, uint64_t* rhs, uint64_t* res, const uint64_t* MODULUS);

void fq2_add(uint64_t* self, uint64_t* rhs, uint64_t* res, const uint64_t* MODULUS);

void fq2_mul(uint64_t* self, uint64_t* other, uint64_t* res, const uint64_t INV, const uint64_t* MODULUS);

void fq2_square(uint64_t* self, uint64_t* res, const uint64_t INV, const uint64_t* MODULUS);

void fq2_invert(uint64_t* self, const uint64_t* exp, uint64_t* res, const uint64_t* R, const uint64_t INV, const uint64_t* MODULUS);

// void fq2_pow(uint64_t* self, const uint64_t* by, uint64_t* res, const uint64_t* R, const uint64_t INV, const uint64_t* MODULUS);