#ifndef __FQ_H__
#define __FQ_H__
#include<stdio.h>
#include<stdint.h>
#include <endian.h>
#include <string.h>

/// This represents an element of $\mathbb{F}_q$ where
///
/// `p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47`
///
/// is the base field of the BN254 curve.
// The internal representation of this type is four 64-bit unsigned
// integers in little-endian order. `Fq` values are always in
// Montgomery form; i.e., Fq(a) = aR mod q, with R = 2^256.

///we need define parameters for different field
///Fq parameters
//Fq modulus 
const uint64_t MODULUS_q[4] = { 0x3c208c16d87cfd47,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029};

//Fq INV = -(q^{-1} mod 2^64) mod 2^64
const uint64_t INV_q= 0x87d20782e4866389;

//Fq R = 2^256 mod q
const uint64_t R_q[4]={
    0xd35d438dc58f0d9d,
    0x0a78eb28f5c70b3d,
    0x666ea36f7879462c,
    0x0e0a77c19a07df2f
};

/// R^2 = 2^512 mod q
const uint64_t R2_q[4]={
    0xf32cfc5b538afa89,
    0xb5e71911d44501fb,
    0x47ab1eff0a417ff6,
    0x06d89f71cab8351f
};

/// R^3 = 2^768 mod q
const uint64_t R3_q[4]={
    0xb1cd6dafda1530df,
    0x62f210e6a7283db6,
    0xef7f0b0c0ada0afb,
    0x20fd6e902d592544
};



#endif
