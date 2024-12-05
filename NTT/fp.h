#ifndef __FP_H__
#define __FP_H__
#include "fields.h"

// using namespace std;

// This represents an element of $\mathbb{F}_p$ where
//
// `p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001`
//
// is the base field of the Pallas curve.
// The internal representation of this type is four 64-bit unsigned
// integers in little-endian order. `Fp` values are always in
// Montgomery form; i.e., Fp(a) = aR mod p, with R = 2^256.

/// Constant representing the modulus
/// p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
const uint64_t MODULUS[4] = { 0x992d30ed00000001, 0x224698fc094cf91b, 0x0000000000000000, 0x4000000000000000 };

/// INV = -(p^{-1} mod 2^64) mod 2^64
const uint64_t INV= 0x992d30ecffffffff;

/// R = 2^256 mod p
const uint64_t R[4]={
    0x34786d38fffffffd,
    0x992c350be41914ad,
    0xffffffffffffffff,
    0x3fffffffffffffff
};

/// R^2 = 2^512 mod p
const uint64_t R2[4]={
    0x8c78ecb30000000f,
    0xd7d30dbd8b0de0e7,
    0x7797a99bc3c95d18,
    0x096d41af7b9cb714
};

/// R^3 = 2^768 mod p
const uint64_t R3[4]={
    0xf185a5993a9e10f9,
    0xf6a68f3b6ac5b1d1,
    0xdf8d1014353fd42c,
    0x2ae309222d2d9910,
};

//raw 2^s root of unity
const uint64_t ROU[4]={0xbdad6fabd87ea32f,0xea322bf2b7bb7584,0x362120830561f81a,0x2bce74deac30ebda};

const uint64_t ROOT_OF_UNITY_INV[4]={
    0xf0b87c7db2ce91f6,
    0x84a0a1d8859f066f,
    0xb4ed8e647196dad1,
    0x2cd5282c53116b5c
};

void SUB(uint64_t* self, const uint64_t* rhs, uint64_t* result);

void ADD(uint64_t* self, const uint64_t* rhs, uint64_t* result);

void montgomery_reduce(
        uint64_t r0,
        uint64_t r1,
        uint64_t r2,
        uint64_t r3,
        uint64_t r4,
        uint64_t r5,
        uint64_t r6,
        uint64_t r7,
        uint64_t* result
    );


// Schoolbook multiplication
void MUL(uint64_t* self, const uint64_t* rhs, uint64_t* result);

/// Squares this element.
void SQUARE(uint64_t* self);

void from(uint64_t val,uint64_t* result);

/// Returns zero, the additive identity.
// uint64_t* zero();

/// Returns one, the multiplicative identity.
void one(uint64_t* res);

void from_raw(uint64_t* val, uint64_t* result);

void from_u512(uint64_t arr[8], uint64_t* res);

void Random(uint64_t* res);

//Return 2^s root of unity
void root_of_unity(uint64_t* result);

//uint64_t* pow_vartime(uint64_t* self, uint64_t* exp);

/// Computes the multiplicative inverse of this element,
/// failing if the element is zero.
//Fp invert(Fp* self);

#endif