#ifndef __FR_H__
#define __FR_H__
#include<stdio.h>
#include<stdint.h>
#include <endian.h>
#include <string.h>
#include "fields.h"

///Fr parameters
//Fr modulus

const uint64_t MODULUS_r[4]={
    0x43e1f593f0000001,
    0x2833e84879b97091,
    0xb85045b68181585d,
    0x30644e72e131a029
};

//Fr INV = -(r^{-1} mod 2^64) mod 2^64
const uint64_t INV_r= 0xc2e1f593efffffff;

const uint64_t R_r[4]={
    0xac96341c4ffffffb,
    0x36fc76959f60cd29,
    0x666ea36f7879462e,
    0x0e0a77c19a07df2f
};

const uint64_t R2_r[4]={
    0x1bb8e645ae216da7,
    0x53fe3ab1e35c59e3,
    0x8c49833d53bb8085,
    0x0216d0b17f4e44a5
};

const uint64_t R3_r[4]={
    0x5e94d8e1b4bf0040,
    0x2a489cbe1cfbb6b8,
    0x893cc664a19fcfed,
    0x0cf8594b7fcc657c
};

const uint64_t exp_r[4]={
    0x43e1f593efffffff,
    0x2833e84879b97091,
    0xb85045b68181585d,
    0x30644e72e131a029,
};



#endif