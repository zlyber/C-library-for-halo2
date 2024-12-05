#ifndef __FR_H__
#define __FR_H__
#include "fields.h"

/// This represents an element of $\mathbb{F}_r$ where
///
/// `r = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001`
///
/// is the scalar field of the BN254 curve.
// The internal representation of this type is four 64-bit unsigned
// integers in little-endian order. `Fr` values are always in
// Montgomery form; i.e., Fr(a) = aR mod r, with R = 2^256.

const uint64_t MODULUS[4] = { 0x43e1f593f0000001,
    0x2833e84879b97091,
    0xb85045b68181585d,
    0x30644e72e131a029};

/// INV = -(p^{-1} mod 2^64) mod 2^64
const uint64_t INV= 0xc2e1f593efffffff;

/// R = 2^256 mod p
const uint64_t R[4]={
    0xac96341c4ffffffb,
    0x36fc76959f60cd29,
    0x666ea36f7879462e,
    0x0e0a77c19a07df2f
};

/// R^2 = 2^512 mod p
const uint64_t R2[4]={
    0x1bb8e645ae216da7,
    0x53fe3ab1e35c59e3,
    0x8c49833d53bb8085,
    0x0216d0b17f4e44a5
};

/// R^3 = 2^768 mod p
const uint64_t R3[4]={
    0x5e94d8e1b4bf0040,
    0x2a489cbe1cfbb6b8,
    0x893cc664a19fcfed,
    0x0cf8594b7fcc657c
};

//raw 2^s root of unity
const uint64_t ROU[4]={0xd34f1ed960c37c9c,
    0x3215cf6dd39329c8,
    0x98865ea93dd31f74,
    0x03ddb9f5166d18b7};

///raw 1 / ROOT_OF_UNITY mod r
const uint64_t ROOT_OF_UNITY_INV[4]={
    0x0ed3e50a414e6dba,
    0xb22625f59115aba7,
    0x1bbe587180f34361,
    0x048127174daabc26
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