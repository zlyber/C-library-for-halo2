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
inline void adc(uint64_t a, uint64_t b, uint64_t carry, uint64_t* result, uint64_t* new_carry);


/// Compute a - (b + borrow), returning the result and the new borrow.
inline void sbb(uint64_t a, uint64_t b, uint64_t borrow, uint64_t* result, uint64_t* new_borrow);

// Compute a + (b * c) + carry, returning the result and the new carry over.
inline void mac(uint64_t a, uint64_t b, uint64_t c, uint64_t carry, uint64_t* result, uint64_t* new_carry);

inline uint64_t next_u64(void);

void SUB(uint64_t* self, const uint64_t* rhs,uint64_t* result,const uint64_t* MODULUS);

void ADD(uint64_t* self, const uint64_t* rhs, uint64_t* result, const uint64_t* MODULUS);

void montgomery_reduce(
        uint64_t r0,
        uint64_t r1,
        uint64_t r2,
        uint64_t r3,
        uint64_t r4,
        uint64_t r5,
        uint64_t r6,
        uint64_t r7,
        uint64_t* result,
        const uint64_t INV,
        const uint64_t* MODULUS
    ); 



// Schoolbook multiplication
void MUL(uint64_t* self, const uint64_t* rhs, uint64_t* result, const uint64_t INV, const uint64_t* MODULUS); 

/// Squares this element.
void SQUARE(uint64_t* self,uint64_t* result, const uint64_t INV, const uint64_t* MODULUS); 

//judge two 256bit field element's equality
bool is_equal(uint64_t* a, uint64_t* b);

//copy a 256bit field element 
void u64_to_u64(uint64_t* a,uint64_t* b);

#endif