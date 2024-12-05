#ifndef __FIELDS_H__
#define __FIELDS_H__
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#define ARRAY_SIZE 8

// using namespace std;

typedef __uint128_t uint128_t;


/// Compute a + b + carry, returning the result and the new carry over.
inline void adc(uint64_t a, uint64_t b, uint64_t carry, uint64_t* result, uint64_t* new_carry) {
    uint128_t ret = (uint128_t)a + (uint128_t)b + (uint128_t)carry;
    *result = (uint64_t)ret;
    *new_carry = (uint64_t)(ret >> 64);
}


/// Compute a - (b + borrow), returning the result and the new borrow.
inline void sbb(uint64_t a, uint64_t b, uint64_t borrow, uint64_t* result, uint64_t* new_borrow) {
    uint128_t ret = (uint128_t)a - (uint128_t)b - (uint128_t)(borrow >> 63);
    *result = (uint64_t)ret;
    *new_borrow = (uint64_t)(ret >> 64);
}

// Compute a + (b * c) + carry, returning the result and the new carry over.
inline void mac(uint64_t a, uint64_t b, uint64_t c, uint64_t carry, uint64_t* result, uint64_t* new_carry) {
    uint128_t ret = (uint128_t)a + ((uint128_t)b * (uint128_t)c) + (uint128_t)carry;
    *result = (uint64_t)ret;
    *new_carry = (uint64_t)(ret >> 64);
}

inline uint64_t next_u64(void) {
    uint64_t r = 0;
    for (int i = 0; i < 64; i += 8) {
        r |= ((uint64_t)rand() & 0xff) << i;
    }
    return r;
}

#endif