#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <endian.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

typedef __uint128_t uint128_t;
#define ARRAY_SIZE 8
//we need define parameters for different field
//fp parameters
const uint64_t MODULUS_p[4] = { 
    0x43e1f593f0000001,
    0x2833e84879b97091,
    0xb85045b68181585d,
    0x30644e72e131a029};

const uint64_t INV_p= 0xc2e1f593efffffff;

const uint64_t R_p[4]={
    0xac96341c4ffffffb,
    0x36fc76959f60cd29,
    0x666ea36f7879462e,
    0x0e0a77c19a07df2f
};

const uint64_t R2_p[4]={
    0x1bb8e645ae216da7,
    0x53fe3ab1e35c59e3,
    0x8c49833d53bb8085,
    0x0216d0b17f4e44a5
};

const uint64_t R3_p[4]={
    0x5e94d8e1b4bf0040,
    0x2a489cbe1cfbb6b8,
    0x893cc664a19fcfed,
    0x0cf8594b7fcc657c
};

const uint64_t ROU_p[4]={
    0xd34f1ed960c37c9c,
    0x3215cf6dd39329c8,
    0x98865ea93dd31f74,
    0x03ddb9f5166d18b7};

const uint64_t ROOT_OF_UNITY_INV_p[4]={
    0x0ed3e50a414e6dba,
    0xb22625f59115aba7,
    0x1bbe587180f34361,
    0x048127174daabc26
};

///Fq parameters
//Fq modulus 
const uint64_t MODULUS_q[4] = { 
    0x3c208c16d87cfd47,
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
const uint64_t exp_q[4] = {
    0x3c208c16d87cfd45,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029
};
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

void sbb(uint64_t a, uint64_t b, uint64_t borrow, uint64_t* result, uint64_t* new_borrow) {
    uint128_t ret = (uint128_t)a - (uint128_t)b - (uint128_t)(borrow >> 63);
    *result = (uint64_t)ret;
    *new_borrow = (uint64_t)(ret >> 64);
}
void adc(uint64_t a, uint64_t b, uint64_t carry, uint64_t* result, uint64_t* new_carry) {
    uint128_t ret = (uint128_t)a + (uint128_t)b + (uint128_t)carry;
    *result = (uint64_t)ret;
    *new_carry = (uint64_t)(ret >> 64);
}
void mac(uint64_t a, uint64_t b, uint64_t c, uint64_t carry, uint64_t* result, uint64_t* new_carry) {
    uint128_t ret = (uint128_t)a + ((uint128_t)b * (uint128_t)c) + (uint128_t)carry;
    *result = (uint64_t)ret;
    *new_carry = (uint64_t)(ret >> 64);
}

uint64_t next_u64(void) {
    uint64_t r = 0;
    for (int i = 0; i < 64; i += 8) {
        r |= ((uint64_t)rand() & 0xff) << i;
    }
    return r;
}
void SUB(uint64_t* self, const uint64_t* rhs,uint64_t* result,const uint64_t* MODULUS) {
        uint64_t borrow=0;
        uint64_t d0,d1,d2,d3;
        sbb(self[0], rhs[0], borrow, &d0, &borrow);
        sbb(self[1], rhs[1], borrow, &d1, &borrow);
        sbb(self[2], rhs[2], borrow, &d2, &borrow);
        sbb(self[3], rhs[3], borrow, &d3, &borrow);
        uint64_t carry=0;
        uint64_t Ins;   
        adc(d0, MODULUS[0] & borrow, carry, &d0, &carry);
        adc(d1, MODULUS[1] & borrow, carry, &d1, &carry);
        adc(d2, MODULUS[2] & borrow, carry, &d2, &carry);
        adc(d3, MODULUS[3] & borrow, carry, &d3, &Ins);
        result[0]=d0;
        result[1]=d1;
        result[2]=d2;
        result[3]=d3;
}
//get the negative self
void NEG(uint64_t* self, uint64_t* result,const uint64_t* MODULUS){
        // Subtract `self` from `MODULUS` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        uint64_t borrow=0;
        uint64_t d0,d1,d2,d3;
        sbb(MODULUS[0],self[0],borrow,&d0,&borrow);
        sbb(MODULUS[1],self[1],borrow,&d1,&borrow);
        sbb(MODULUS[2],self[2],borrow,&d2,&borrow);
        sbb(MODULUS[3],self[3],borrow,&d3,&borrow);
        // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        uint64_t mask = ((self[0] | self[1] | self[2] | self[3]) == 0) ? 0 : UINT64_MAX;

        result[0]=d0 & mask;
        result[1]=d1 & mask;
        result[2]=d2 & mask;
        result[3]=d3 & mask;
}
void ADD(uint64_t* self, const uint64_t* rhs, uint64_t* result, const uint64_t* MODULUS) {
        uint64_t carry=0;
        uint64_t Ins;
        uint64_t d0,d1,d2,d3;
        adc(self[0], rhs[0], carry, &d0, &carry);
        adc(self[1], rhs[1], carry, &d1, &carry);
        adc(self[2], rhs[2], carry, &d2, &carry);
        adc(self[3], rhs[3], carry, &d3, &Ins);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        uint64_t f[4]={d0, d1, d2, d3};
        SUB(f, MODULUS, result, MODULUS);
}
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
    ) {
        // The Montgomery reduction here is based on Algorithm 14.32 in
        // Handbook of Applied Cryptography
        // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

        uint64_t k = r0 * INV;   //INV=m', b=2^64 , correspond to 2.1 : u_i=a_i*m' mod b
        uint64_t carry=0;
        uint64_t carry2;
        uint64_t Ins;
        mac(r0, k, MODULUS[0], carry, &Ins, &carry); //correspond to 2.2 : A=A+u_i*m*b^i
        mac(r1, k, MODULUS[1], carry, &r1, &carry);
        mac(r2, k, MODULUS[2], carry, &r2, &carry);
        mac(r3, k, MODULUS[3], carry, &r3, &carry);
        adc(r4, 0, carry, &r4, &carry2);

        k = r1 * INV;
        carry=0;
        mac(r1, k, MODULUS[0], carry, &Ins, &carry);
        mac(r2, k, MODULUS[1], carry, &r2, &carry);
        mac(r3, k, MODULUS[2], carry, &r3, &carry);
        mac(r4, k, MODULUS[3], carry, &r4, &carry);
        adc(r5, carry2, carry, &r5, &carry2);

        k = r2 * INV;
        carry=0;
        mac(r2, k, MODULUS[0], carry, &Ins, &carry);
        mac(r3, k, MODULUS[1], carry, &r3, &carry);
        mac(r4, k, MODULUS[2], carry, &r4, &carry);
        mac(r5, k, MODULUS[3], carry, &r5, &carry);
        adc(r6, carry2, carry, &r6, &carry2);

        k = r3 * INV;
        carry=0;
        mac(r3, k, MODULUS[0], carry, &Ins, &carry);
        mac(r4, k, MODULUS[1], carry, &r4, &carry);
        mac(r5, k, MODULUS[2], carry, &r5, &carry);
        mac(r6, k, MODULUS[3], carry, &r6, &carry);
        adc(r7, carry2, carry, &r7, &Ins);
        
        
        // Result may be within MODULUS of the correct value
        uint64_t f[4]={r4,r5,r6,r7}; // correspond to 3 : A=A/b^n ,b=2^64,n=4
        SUB(f, MODULUS, result, MODULUS);
    }
void MUL(uint64_t* self, const uint64_t* rhs, uint64_t* result, const uint64_t INV, const uint64_t* MODULUS) {
    uint64_t r0, r1, r2, r3, r4, r5, r6, r7, carry;
    carry=0;
    mac(0, self[0], rhs[0], carry, &r0, &carry);
    mac(0, self[0], rhs[1], carry, &r1, &carry);
    mac(0, self[0], rhs[2], carry, &r2, &carry);
    mac(0, self[0], rhs[3], carry, &r3, &r4);
    carry=0;
    mac(r1, self[1], rhs[0], carry, &r1, &carry);
    mac(r2, self[1], rhs[1], carry, &r2, &carry);
    mac(r3, self[1], rhs[2], carry, &r3, &carry);
    mac(r4, self[1], rhs[3], carry, &r4, &r5);
    carry=0;
    mac(r2, self[2], rhs[0], carry, &r2, &carry);
    mac(r3, self[2], rhs[1], carry, &r3, &carry);
    mac(r4, self[2], rhs[2], carry, &r4, &carry);
    mac(r5, self[2], rhs[3], carry, &r5, &r6);
    carry=0;
    mac(r3, self[3], rhs[0], carry, &r3, &carry);
    mac(r4, self[3], rhs[1], carry, &r4, &carry);
    mac(r5, self[3], rhs[2], carry, &r5, &carry);
    mac(r6, self[3], rhs[3], carry, &r6, &r7);

    montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7, result, INV, MODULUS);
}

/// Squares this element.
void SQUARE(uint64_t* self,uint64_t* result, const uint64_t INV, const uint64_t* MODULUS) {
    uint64_t r0,r1,r2,r3,r4,r5,r6,r7;
    uint64_t carry=0;
    mac(0, self[0], self[1], carry, &r1, &carry);
    mac(0, self[0], self[2], carry, &r2, &carry);
    mac(0, self[0], self[3], carry, &r3, &r4);

    carry=0;
    mac(r3, self[1], self[2], carry, &r3, &carry);
    mac(r4, self[1], self[3], carry, &r4, &r5);
    
    carry=0;
    mac(r5, self[2], self[3], carry, &r5, &r6);

    r7 = r6 >> 63;
    r6 = (r6 << 1) | (r5 >> 63);
    r5 = (r5 << 1) | (r4 >> 63);
    r4 = (r4 << 1) | (r3 >> 63);
    r3 = (r3 << 1) | (r2 >> 63);
    r2 = (r2 << 1) | (r1 >> 63);
    r1 = r1 << 1;

    carry=0;
    uint64_t Ins;
    mac(0, self[0], self[0], carry, &r0, &carry);
    adc(0, r1, carry, &r1, &carry);
    mac(r2, self[1], self[1], carry, &r2, &carry);
    adc(0, r3, carry, &r3, &carry);
    mac(r4, self[2], self[2], carry, &r4, &carry);
    adc(0, r5, carry, &r5, &carry);
    mac(r6, self[3], self[3], carry, &r6, &carry);
    adc(0, r7, carry, &r7, &Ins);

    montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7, result, INV, MODULUS);
}

//judge two 256bit field element's equality
bool is_equal(uint64_t* a, uint64_t* b){
    bool equal=1;
    for(int i=0;i<4;i++){
        equal=equal & (a[i]==b[i]);
    }
    return equal;
}

//copy a 256bit field element 
void u64_to_u64(uint64_t* a,uint64_t* b){
    for(int i=0;i<4;i++){
        a[i]=b[i];
    }
}

void swap(uint64_t* a, uint64_t* b) {
    uint64_t temp[4];
    u64_to_u64(temp,a);
    u64_to_u64(a,b);
    u64_to_u64(b,temp);
}

uint64_t reverse_bits(uint64_t operand,int bit_count){
    uint64_t acc = 0;
    for (int i = 0; i < bit_count; i++) {
      acc = acc << 1;
      acc |= (operand >> i) & 1;
    }
    return acc;
}

void fp_one(uint64_t* res){
    for(int i=0;i<4;i++){
        res[i]=R_p[i];
    }
}

void group_scale(uint64_t* self, uint64_t* rhs, uint64_t*res){
    MUL(self, rhs, res, INV_p, MODULUS_p); 
}
void group_add(uint64_t* self, uint64_t* rhs, uint64_t* res){
    ADD(self, rhs, res, MODULUS_p);
}
void group_sub(uint64_t* self, uint64_t* rhs, uint64_t* res){
    SUB(self, rhs, res, MODULUS_p);
}
void scalar_one(uint64_t* result){
    fp_one(result);
}

void precompute_twiddles(uint64_t** twiddles, uint64_t *omega, uint64_t n) {
    uint64_t w[4];
    scalar_one(w);
    for (uint64_t i = 0; i < n / 2; i++) {
        for(int j=0;j<4;j++){
            twiddles[i][j]=w[j];
        }
        group_scale(w, omega, w);
    }
}
void recursive_butterfly(uint64_t** vector,uint64_t N, uint64_t twiddle_chunk, uint64_t** twiddles){
    if (N==2){
        uint64_t t[4];
        u64_to_u64(t,vector[1]);
        u64_to_u64(vector[1],vector[0]);
        group_add(vector[0],t,vector[0]);
        group_sub(vector[1],t,vector[1]);
    }
    else{
        uint64_t** Left;
        uint64_t** Right; 
        Left = (uint64_t**)malloc(sizeof(uint64_t*)*(N/2));
        Right = (uint64_t**)malloc(sizeof(uint64_t*)*(N/2));
        for (uint64_t i=0; i<N/2; i++)
	    {
		    Left[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
            Right[i]= (uint64_t*)malloc(sizeof(uint64_t)*4); 
	    }
        
        for(uint64_t i=0;i<N/2;i++){
            u64_to_u64(Left[i],vector[i]);
            u64_to_u64(Right[i],vector[i+N/2]);
        }
        uint64_t t[4];
        u64_to_u64(t,Right[0]);      
        u64_to_u64(Right[0],Left[0]);
        group_add(Left[0],t,Left[0]);
        group_sub(Right[0],t,Right[0]);
        for(uint64_t m=0;m<N/2-1;m++){
            uint64_t* Left_addr=Left[m+1];
            uint64_t* Right_addr=Right[m+1];
            uint64_t* twiddle=twiddles[(m+1)*twiddle_chunk];
            uint64_t t1[4];
            group_scale(Right_addr,twiddle,t1); 
            u64_to_u64(Right_addr,Left_addr);
            group_add(Left_addr,t1,Left_addr); 
            group_sub(Right_addr,t1,Right_addr); 
        }
        for(uint64_t i=0;i<N/2;i++){
            u64_to_u64(vector[i],Left[i]);
            u64_to_u64(vector[i+N/2],Right[i]);
        }
        for (uint64_t i = 0; i < N/2; i++){
            free(Left[i]);
            free(Right[i]);
        }
        free(Left);
        free(Right);
    }
}
void NTT(uint64_t vector[][4], uint64_t* omega, uint32_t k){
    uint64_t N=1<<k;
    for(uint64_t i=0; i<N; i++){
        uint64_t rk = reverse_bits(i, k);
        if (i < rk){
            swap(vector[rk],vector[i]);
        }
    }
    uint64_t** twiddles;
    twiddles = (uint64_t**)malloc(sizeof(uint64_t*)*(N/2));
    for (uint64_t i=0; i<N/2; i++)
	{
		twiddles[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
	}  
    precompute_twiddles(twiddles,omega,N);
    uint64_t chunk=2;
    uint64_t twiddle_chunk = N / 2;
    for(uint32_t i=0;i<k;i++){
        for(uint64_t j=0;j<N;j+=chunk){
            uint64_t* Left=&vector[j][0];
            uint64_t* Right=&vector[j+chunk/2][0];
            uint64_t t[4];
            u64_to_u64(t,Right);    
            u64_to_u64(Right,Left); 
            group_add(Left,t,Left);
            group_sub(Right,t,Right);
            for(uint64_t m=0;m<chunk/2-1;m++){   
                uint64_t* Left_addr=vector[j+m+1];
                uint64_t* Right_addr=vector[j+chunk/2+m+1];
                uint64_t* twiddle=twiddles[(m+1)*twiddle_chunk];
                uint64_t t1[4];
                
                group_scale(Right_addr,twiddle,t1); 
                u64_to_u64(Right_addr,Left_addr);
                group_add(Left_addr,t1,Left_addr); 
                group_sub(Right_addr,t1,Right_addr);  
            }  
        }
    chunk *= 2; // 向上合并
    twiddle_chunk /= 2;
    }
    for (uint64_t i = 0; i < N/2; i++){
        free(twiddles[i]);
    }
    free(twiddles);
}

/// Returns one, the multiplicative identity.
void one(uint64_t* res,const uint64_t* R){
    for(int i=0;i<4;i++){
        res[i]=R[i];
    }
}
/// Exponentiates `self` by `exp`, where `exp` is a little-endian order
/// integer exponent.
///
/// **This operation is variable time with respect to the exponent.** If the
/// exponent is fixed, this operation is effectively constant time.
void pow_vartime(uint64_t* self, const uint64_t exp, uint64_t* res, const uint64_t* R, const uint64_t INV, const uint64_t* MODULUS)
{
        one(res,R);
        for (int i=63;i>=0;i--)
        {
            SQUARE(res,res,INV,MODULUS);

            if (((exp >> i) & 1) == 1) 
            {
                MUL(res,self,res,INV,MODULUS);
            }
        }
}

/// Exponentiates `self` by `by`, where `by` is a little-endian order
/// integer exponent.
void POW(uint64_t* self, const uint64_t* by, uint64_t* res, const uint64_t* R, const uint64_t INV, const uint64_t* MODULUS) 
{
        one(res,R);
        for (int i=3;i>=0;i--){
            for (int j=63;j>=0;j--){
                SQUARE(res,res,INV,MODULUS);
                uint64_t tmp[4];
                u64_to_u64(tmp,res);
                MUL(tmp,self,tmp,INV,MODULUS);
                if (((by[i] >> j) & 0x1)==1)
                {
                    u64_to_u64(res,tmp);
                }
            }
        }
}
/// Computes the multiplicative inverse of this element,
/// failing if the element is zero.
void invert(uint64_t* self, const uint64_t* exp, uint64_t* res, const uint64_t* R, const uint64_t INV, const uint64_t* MODULUS)
 {
    uint64_t tmp[4];
    POW(self,exp,tmp,R,INV,MODULUS);
    if((self[0]||self[1]||self[2]||self[3])!=0){
        u64_to_u64(res,tmp);
    }
    else{
        printf("inversion not exist\n");
    }
}


//rescale to Mongomery form
void from_raw(uint64_t* val, uint64_t* result, const uint64_t* R2, const uint64_t INV, const uint64_t* MODULUS){
    MUL(val, R2, result, INV, MODULUS);
}

void from(uint64_t val, uint64_t* result, const uint64_t* R2, const uint64_t INV, const uint64_t* MODULUS) {
    uint64_t f[4]={val, 0, 0, 0};   
    MUL(f, R2, result, INV, MODULUS);
}

void from_u512(uint64_t arr[8], uint64_t* res){
    uint64_t d0[4]={arr[0],arr[1],arr[2],arr[3]};
    uint64_t d1[4]={arr[4],arr[5],arr[6],arr[7]};
    MUL(d0,R2_p,d0,INV_p,MODULUS_p);
    MUL(d1,R3_p,d1,INV_p,MODULUS_p);
    ADD(d0,d1,res,MODULUS_p);
}

void Random(uint64_t* res){
    uint64_t arr[8];
    for(int i=0;i<8;i++){
        arr[i]=next_u64();
    }
    return from_u512(arr, res);
}

//Return 2^s root of unity
void root_of_unity(uint64_t* result){
    from_raw((uint64_t*)ROU_p,result,R2_p,INV_p,MODULUS_p);
}

void to_repr(uint64_t* self,uint8_t* res,const uint64_t INV, const uint64_t* MODULUS) {
    // Turn into canonical form by computing
    // (a.R) / R = a
    uint64_t tmp[4];
    montgomery_reduce(self[0],self[1],self[2],self[3],0,0,0,0,tmp,INV,MODULUS);

    res[0] = tmp[0] & 0xff;
    res[1] = (tmp[0] >> 8) & 0xff;
    res[2] = (tmp[0] >> 16) & 0xff;
    res[3] = (tmp[0] >> 24) & 0xff;
    res[4] = (tmp[0] >> 32) & 0xff;
    res[5] = (tmp[0] >> 40) & 0xff;
    res[6] = (tmp[0] >> 48) & 0xff;
    res[7] = (tmp[0] >> 56) & 0xff;
    res[8] = tmp[1] & 0xff;
    res[9] = (tmp[1] >> 8) & 0xff;
    res[10] = (tmp[1] >> 16) & 0xff;
    res[11] = (tmp[1] >> 24) & 0xff;
    res[12] = (tmp[1] >> 32) & 0xff;
    res[13] = (tmp[1] >> 40) & 0xff;
    res[14] = (tmp[1] >> 48) & 0xff;
    res[15] = (tmp[1] >> 56) & 0xff;
    res[16] = tmp[2] & 0xff;
    res[17] = (tmp[2] >> 8) & 0xff;
    res[18] = (tmp[2] >> 16) & 0xff;
    res[19] = (tmp[2] >> 24) & 0xff;
    res[20] = (tmp[2] >> 32) & 0xff;
    res[21] = (tmp[2] >> 40) & 0xff;
    res[22] = (tmp[2] >> 48) & 0xff;
    res[23] = (tmp[2] >> 56) & 0xff;
    res[24] = tmp[3] & 0xff;
    res[25] = (tmp[3] >> 8) & 0xff;
    res[26] = (tmp[3] >> 16) & 0xff;
    res[27] = (tmp[3] >> 24) & 0xff;
    res[28] = (tmp[3] >> 32) & 0xff;
    res[29] = (tmp[3] >> 40) & 0xff;
    res[30] = (tmp[3] >> 48) & 0xff;
    res[31] = (tmp[3] >> 56) & 0xff;

}

void affine_identity(uint64_t* x,uint64_t* y){
    for(int i=0;i<4;i++){
        x[i]=0;
        y[i]=0;
    }
}

//return the identity of the curve in projective form
void identity(uint64_t* x,uint64_t* y,uint64_t* z){
    for(int i=0;i<4;i++){
        x[i]=0;
        y[i]=0;
        z[i]=0;
    }
}

//judge the affine coordinate is the curve's identity?
bool is_identity(uint64_t* x,uint64_t* y) {
    bool is_zero=1;
    for(int i=0;i<4;i++){
        is_zero=is_zero & (x[i]==0) & (y[i]==0);
    }
    return is_zero;
}

//judge the projective coordinate is the curve's identity?
bool is_identity_project(uint64_t* x,uint64_t* y,uint64_t* z) {
    bool is_zero=1;
    for(int i=0;i<4;i++){
        is_zero=is_zero & (x[i]==0) & (y[i]==0) & (z[i]==0);
    }
    return is_zero;
}

//affine to projective 
void to_curve(uint64_t* x,uint64_t* y,uint64_t* z) {
    bool ct=is_identity(x,y);
    for(int i=0;i<4;i++){
        z[i]=(ct)?0:R_q[i];
    }
}

//point double in Jacobian form: 2(x,y,z)
void point_double(uint64_t* x,uint64_t* y,uint64_t* z){
    uint64_t a[4],b[4],c[4],d[4],e[4],f[4];
    SQUARE(x,a,INV_q,MODULUS_q);
    SQUARE(y,b,INV_q,MODULUS_q);
    SQUARE(b,c,INV_q,MODULUS_q);
    ADD(x,b,d,MODULUS_q);
    SQUARE(d,d,INV_q,MODULUS_q);
    SUB(d,a,d,MODULUS_q);
    SUB(d,c,d,MODULUS_q);
    ADD(d,d,d,MODULUS_q);
    ADD(a,a,e,MODULUS_q);
    ADD(e,a,e,MODULUS_q);
    SQUARE(e,f,INV_q,MODULUS_q);
    
    uint64_t x3[4],y3[4],z3[4];
    uint64_t mid1[4],mid2[4]; //存储中间变量
    MUL(z,y,z3,INV_q,MODULUS_q);
    ADD(z3,z3,z3,MODULUS_q);
    ADD(d,d,mid1,MODULUS_q);
    SUB(f,mid1,x3,MODULUS_q);
    ADD(c,c,c,MODULUS_q);
    ADD(c,c,c,MODULUS_q);
    ADD(c,c,c,MODULUS_q);
    SUB(d,x3,mid1,MODULUS_q);
    MUL(e,mid1,mid2,INV_q,MODULUS_q);
    SUB(mid2,c,y3,MODULUS_q);

    bool ct=is_identity_project(x,y,z);
    for(int i=0;i<4;i++){
        x[i]=(ct)?0:x3[i];
        y[i]=(ct)?0:y3[i];
        z[i]=(ct)?0:z3[i];
    }
}

//point add in affine form : (x,y)+(X,Y) and convert to projection
void affine_add_affine(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* X,uint64_t* Y,uint64_t* Z){
    if(is_identity(x,y)){
        to_curve(X,Y,Z);
        u64_to_u64(x,X);
        u64_to_u64(y,Y);
        u64_to_u64(z,Z);
    }
    else if(is_identity(X,Y)){
        to_curve(x,y,z);
    }
    else
    {
        if(is_equal(x,X)){
            if(is_equal(y,Y)){
                to_curve(x,y,z);
                point_double(x,y,z);
            }
            else{
                identity(x,y,z);
            }
        }
        else
        {
            uint64_t h[4],hh[4],i[4],j[4],r[4],v[4];
            SUB(X,x,h,MODULUS_q);
            SQUARE(h,hh,INV_q,MODULUS_q);
            ADD(hh,hh,i,MODULUS_q);
            ADD(i,i,i,MODULUS_q);
            MUL(h,i,j,INV_q,MODULUS_q);
            SUB(Y,y,r,MODULUS_q);
            ADD(r,r,r,MODULUS_q);
            MUL(x,i,v,INV_q,MODULUS_q);

            uint64_t mid1[4],mid2[4]; //存储中间结果
            uint64_t x3[4],y3[4],z3[4];
            SQUARE(r,mid1,INV_q,MODULUS_q);
            SUB(mid1,j,mid1,MODULUS_q);
            SUB(mid1,v,mid2,MODULUS_q);
            SUB(mid2,v,x3,MODULUS_q);

            MUL(y,j,j,INV_q,MODULUS_q);
            ADD(j,j,j,MODULUS_q);
            SUB(v,x3,mid1,MODULUS_q);
            MUL(r,mid1,mid2,INV_q,MODULUS_q);
            SUB(mid2,j,y3,MODULUS_q);

            ADD(h,h,z3,MODULUS_q);
            u64_to_u64(x,x3);
            u64_to_u64(y,y3);
            u64_to_u64(z,z3);
        }
    }
}

//projective point add affine point : (x,y,z)+(X,Y) and convert to projection
void project_add_affine(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* X,uint64_t* Y,uint64_t* Z){
    if(is_identity_project(x,y,z)){
        to_curve(X,Y,Z);
        u64_to_u64(x,X);
        u64_to_u64(y,Y);
        u64_to_u64(z,Z);
    }
    else if (is_identity(X,Y))
    {
        return;
    }
    else{
        uint64_t z1z1[4],u2[4],s2[4];
        SQUARE(z,z1z1,INV_q,MODULUS_q);
        MUL(X,z1z1,u2,INV_q,MODULUS_q);
        MUL(Y,z1z1,s2,INV_q,MODULUS_q);
        MUL(s2,z,s2,INV_q,MODULUS_q);

        if(is_equal(x,u2)){
            if(is_equal(y,s2)){
                point_double(x,y,z);
            }
            else{
                identity(x,y,z);
            }
        }
        else
        {
            uint64_t h[4],hh[4],i[4],j[4],r[4],v[4];
            SUB(u2,x,h,MODULUS_q);
            SQUARE(h,hh,INV_q,MODULUS_q);
            ADD(hh,hh,i,MODULUS_q);
            ADD(i,i,i,MODULUS_q);
            MUL(h,i,j,INV_q,MODULUS_q);
            SUB(s2,y,r,MODULUS_q);
            ADD(r,r,r,MODULUS_q);
            MUL(x,i,v,INV_q,MODULUS_q);

            uint64_t mid1[4],mid2[4];
            uint64_t x3[4],y3[4],z3[4];
            SQUARE(r,mid1,INV_q,MODULUS_q);
            SUB(mid1,j,mid1,MODULUS_q);
            SUB(mid1,v,mid2,MODULUS_q);
            SUB(mid2,v,x3,MODULUS_q);

            MUL(y,j,j,INV_q,MODULUS_q);
            ADD(j,j,j,MODULUS_q);
            SUB(v,x3,mid1,MODULUS_q);
            MUL(mid1,r,mid2,INV_q,MODULUS_q);
            SUB(mid2,j,y3,MODULUS_q);

            ADD(z,h,mid1,MODULUS_q);
            SQUARE(mid1,mid1,INV_q,MODULUS_q);
            SUB(mid1,z1z1,mid2,MODULUS_q);
            SUB(mid2,hh,z3,MODULUS_q);

            u64_to_u64(x,x3);
            u64_to_u64(y,y3);
            u64_to_u64(z,z3);
        }
    }
    
}

//two projective point add : (x,y,z)+(X,Y,Z)
void project_add_project(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* X,uint64_t* Y,uint64_t* Z){
    if(is_identity_project(x,y,z)){
        u64_to_u64(x,X);
        u64_to_u64(y,Y);
        u64_to_u64(z,Z);
    }
    else if(is_identity_project(X,Y,Z)){
        return;
    }
    else{
        uint64_t z1z1[4],z2z2[4];
        SQUARE(z,z1z1,INV_q,MODULUS_q);         //z^2
        SQUARE(Z,z2z2,INV_q,MODULUS_q);         //Z^2

        uint64_t u1[4],u2[4];
        MUL(x,z2z2,u1,INV_q,MODULUS_q);
        MUL(X,z1z1,u2,INV_q,MODULUS_q);

        uint64_t s1[4],s2[4];
        MUL(y,z2z2,s1,INV_q,MODULUS_q);
        MUL(s1,Z,s1,INV_q,MODULUS_q);
        MUL(Y,z1z1,s2,INV_q,MODULUS_q);
        MUL(s2,z,s2,INV_q,MODULUS_q);

        if(is_equal(u1,u2)){
            if(is_equal(s1,s2)){
                point_double(x,y,z); 
            }
            else{
                identity(x,y,z); //return identity
            }
        }
        else{
            uint64_t h[4],i[4],j[4],r[4],v[4];
            SUB(u2,u1,h,MODULUS_q);
            ADD(h,h,i,MODULUS_q);
            SQUARE(i,i,INV_q,MODULUS_q);
            MUL(h,i,j,INV_q,MODULUS_q);
            SUB(s2,s1,r,MODULUS_q);
            ADD(r,r,r,MODULUS_q);
            MUL(u1,i,v,INV_q,MODULUS_q);

            uint64_t x3[4],y3[4],z3[4],rr[4];
            SQUARE(r,rr,INV_q,MODULUS_q); //r^2

            uint64_t mid1[4],mid2[4]; //用来存储中间变量
            SUB(rr,j,mid1,MODULUS_q);
            SUB(mid1,v,mid2,MODULUS_q);
            SUB(mid2,v,x3,MODULUS_q);

            MUL(s1,j,s1,INV_q,MODULUS_q);
            ADD(s1,s1,s1,MODULUS_q);
            SUB(v,x3,mid1,MODULUS_q);
            MUL(mid1,r,mid2,INV_q,MODULUS_q);
            SUB(mid2,s1,y3,MODULUS_q);

            ADD(z,Z,mid1,MODULUS_q);
            SQUARE(mid1,mid1,INV_q,MODULUS_q);
            SUB(mid1,z1z1,mid2,MODULUS_q);
            SUB(mid2,z2z2,z3,MODULUS_q);
            MUL(z3,h,z3,INV_q,MODULUS_q);

            u64_to_u64(x,x3);
            u64_to_u64(y,y3);
            u64_to_u64(z,z3);
        }
    }
}



uint64_t get_at(int segment, int c, uint8_t* bytes) {
    int skip_bits = segment * c;
    int skip_bytes = skip_bits / 8;

    if (skip_bytes >= 32) {
        return 0;
    }
    
    uint8_t v[8] = {0};
    int len = 32 - skip_bytes;
    if (len > 8) {
        len = 8;
    }
    memcpy(v, bytes + skip_bytes, len);

    uint64_t tmp = *(uint64_t*)v;
    tmp >>= skip_bits - (skip_bytes * 8);
    tmp = tmp % (1 << c);

    return tmp;
}

//judge whether the bucket is null now
bool is_null(uint64_t* x){
    bool is_one=1;
    for(int i=0;i<4;i++){
        is_one=is_one&(x[i]==1);
    }
    return is_one;
}

//MSM
void MSM(uint64_t** coeffs, uint64_t** bases_x, uint64_t** bases_y, uint64_t** bases_z, uint64_t* acc_x, uint64_t* acc_y, uint64_t* acc_z, uint64_t len){
    

    uint8_t** coeffs_repr;
    coeffs_repr = (uint8_t**)malloc(sizeof(uint8_t*)*len);
    for (uint64_t i=0; i<len; i++)
	{
		coeffs_repr[i] = (uint8_t*)malloc(sizeof(uint8_t)*32); 
	} 
    for (uint64_t i = 0; i < len; i++) {
        to_repr(coeffs[i],coeffs_repr[i],INV_r,MODULUS_r);
    }

    int c;
    if(len<4){
        c=1;
    }
    else if(len<32){
        c=3;
    }
    else{
        c = (int)ceil(log(len));
    }
    
    int segments = (256 / c) + 1;

    for (int current_segment = segments - 1; current_segment >= 0; current_segment--) {
        for (int i = 0; i < c; i++) {
            point_double(acc_x,acc_y,acc_z); 
        }

        //为buckets创建空间，要创建三个buckets来代表x,y,z三维
        uint64_t** buckets_x; 
        uint64_t** buckets_y;
        uint64_t** buckets_z;
        uint64_t bucket_len=(1<<c)-1;
        buckets_x=(uint64_t**)malloc(sizeof(uint64_t*)*bucket_len);
        buckets_y=(uint64_t**)malloc(sizeof(uint64_t*)*bucket_len);
        buckets_z=(uint64_t**)malloc(sizeof(uint64_t*)*bucket_len);
        for (uint64_t i=0; i<bucket_len; i++)
        {
           //为每列分配4个uint64_t大小的空间，并将初始值都设为0
           buckets_x[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
           buckets_y[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
           buckets_z[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 

           //初始化为1，视为空
           for(int j=0;j<4;j++){
            buckets_x[i][j]=1;
            buckets_y[i][j]=1;
            buckets_z[i][j]=1;
           }
        }
        for (uint64_t i = 0; i < len; i++) {
            uint64_t coeff = get_at(current_segment, c, coeffs_repr[i]);
            //以下分支达到和add_assign函数一样的效果
            if (coeff != 0) {
                //若当前桶为空，则先把坐标值放进去
                if(is_null(buckets_x[coeff - 1])&is_null(buckets_y[coeff - 1])&is_null(buckets_z[coeff - 1]))
                {
                    u64_to_u64(buckets_x[coeff - 1],bases_x[i]);
                    u64_to_u64(buckets_y[coeff - 1],bases_y[i]);
                    
                }
                //当前是仿射坐标形式，使用仿射坐标相加后转换到投影坐标
                else if(is_null(buckets_z[coeff - 1]))
                {
                    affine_add_affine(buckets_x[coeff - 1],buckets_y[coeff - 1],buckets_z[coeff - 1], bases_x[i],bases_y[i],bases_z[i]);
                    
                }
                //当前已经是投影坐标形式，使用投影坐标加仿射坐标，之后转换到仿射坐标
                else 
                {
                    project_add_affine(buckets_x[coeff - 1],buckets_y[coeff - 1],buckets_z[coeff - 1], bases_x[i],bases_y[i],bases_z[i]);
                    
                }
            }
        }
        
        // Summation by parts
        // e.g. 3a + 2b + 1c = a +
        //                    (a) + b +
        //                    ((a) + b) + c
        uint64_t sum_x[4],sum_y[4],sum_z[4];
        identity(sum_x,sum_y,sum_z); 
        for (signed long long i = bucket_len - 1; i >= 0; i--) //注意这里的i被定义为i64类型而不是u64，因为减到0后u64会自动变成2^64-1
        {
            //以下分支达到和add函数一样的效果，和add_assign的区别是add函数里sum是主体,buckets是客体
            //若当前桶为空，则跳过不管
            if(is_null(buckets_x[i])&is_null(buckets_y[i])&is_null(buckets_z[i]))
            {
                project_add_project(acc_x, acc_y, acc_z, sum_x, sum_y, sum_z);
                continue;
            }
            //若当前桶内是仿射坐标形式，使用投影坐标加仿射坐标的方法
            else if(is_null(buckets_z[i]))
            {
                project_add_affine(sum_x, sum_y, sum_z, buckets_x[i],buckets_y[i],buckets_z[i]);
            }
            //若当前已经是投影坐标形式，直接使用两个投影坐标相加的方式
            else 
            {
                project_add_project(sum_x, sum_y, sum_z, buckets_x[i],buckets_y[i],buckets_z[i]);
            }
            //*acc = *acc + &running_sum;
            project_add_project(acc_x, acc_y, acc_z, sum_x, sum_y, sum_z);
            
        }
        for (uint64_t i = 0; i < bucket_len; i++){
            free(buckets_x[i]);
            free(buckets_y[i]);
            free(buckets_z[i]);
        }
        free(buckets_x);
        free(buckets_y);
        free(buckets_z);
        }
    for (uint64_t i = 0; i < len; i++){
        free(coeffs_repr[i]);
    }
    free(coeffs_repr);
}
void fq2_one(uint64_t* self,const uint64_t* R){
    one(self,R);
    for(int i=0;i<4;i++){
        self[i+4]=0;
    }
}

void fq2_zero(uint64_t* self){
    for(int i=0;i<8;i++){
        self[i]=0;
    }
}

void fq2_sub(uint64_t* self, uint64_t* rhs, uint64_t* res, const uint64_t* MODULUS){
    SUB(self,rhs,res,MODULUS);
    SUB(self+4,rhs+4,res+4,MODULUS);
}

void fq2_add(uint64_t* self, uint64_t* rhs, uint64_t* res, const uint64_t* MODULUS){
    ADD(self,rhs,res,MODULUS);
    ADD(self+4,rhs+4,res+4,MODULUS);
}

void fq2_mul(uint64_t* self, uint64_t* other, uint64_t* res, const uint64_t INV, const uint64_t* MODULUS) {
    
    uint64_t t1[4];
    uint64_t t0[4];
    MUL(self,other,t1,INV,MODULUS);
    ADD(self,self+4,t0,MODULUS);

    uint64_t t2[4];
    MUL(self+4,other+4,t2,INV,MODULUS);

    ADD(other,other+4,res+4,MODULUS);
    SUB(t1,t2,res,MODULUS);
    ADD(t1,t2,t1,MODULUS);
    MUL(t0,res+4,t0,INV,MODULUS);
    SUB(t0,t1,res+4,MODULUS);
}

void fq2_square(uint64_t* self, uint64_t* res, const uint64_t INV, const uint64_t* MODULUS) {
    uint64_t ab[4];
    MUL(self,self+4,ab,INV,MODULUS);

    uint64_t c0c1[4];
    ADD(self,self+4,c0c1,MODULUS);

    uint64_t c0[4];
    NEG(self+4,c0,MODULUS);
    ADD(c0,self,c0,MODULUS);
    MUL(c0,c0c1,c0,INV,MODULUS);
    SUB(c0,ab,c0,MODULUS);
    
    ADD(ab,ab,res+4,MODULUS); //self.c1 = ab.double();
    ADD(c0,ab,res,MODULUS);
    
}

/// Computes the multiplicative inverse of this element,
/// failing if the element is zero.
void fq2_invert(uint64_t* self, const uint64_t* exp, uint64_t* res, const uint64_t* R, const uint64_t INV, const uint64_t* MODULUS)
 {
    uint64_t t1[4];
    u64_to_u64(t1,self+4);
    SQUARE(t1,t1,INV,MODULUS);

    uint64_t t0[4];
    u64_to_u64(t0,self);
    SQUARE(t0,t0,INV,MODULUS);
    ADD(t0,t1,t0,MODULUS);

    uint64_t t[4];
    invert(t0,exp,t,R,INV,MODULUS);
    
    u64_to_u64(res,self);
    u64_to_u64(res+4,self+4);
    MUL(res,t,res,INV,MODULUS);
    MUL(res+4,t,res+4,INV,MODULUS);
    NEG(res+4,res+4,MODULUS);
}

void batch_normalize(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* affine_x,uint64_t* affine_y){
    uint64_t acc[4];
    one(acc,R_q);
    u64_to_u64(affine_x,acc);
    if(!(is_identity_project(x,y,z))){
        MUL(acc,z,acc,INV_q,MODULUS_q);
    }
    uint64_t pow[4]={
    0x3c208c16d87cfd45,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029};
    invert(acc,pow,acc,R_q,INV_q,MODULUS_q);
    uint64_t tmp[4];
    MUL(affine_x,acc,tmp,INV_q,MODULUS_q);
    if(!(is_identity_project(x,y,z))){
        MUL(acc,z,acc,INV_q,MODULUS_q);
    }
    uint64_t tmp2[4];
    uint64_t tmp3[4];
    SQUARE(tmp,tmp2,INV_q,MODULUS_q);
    MUL(tmp2,tmp,tmp3,INV_q,MODULUS_q);
    if(!(is_identity_project(x,y,z))){
    MUL(x,tmp2,affine_x,INV_q,MODULUS_q);
    MUL(y,tmp3,affine_y,INV_q,MODULUS_q);
    }
    else{
        affine_identity(affine_x,affine_y);
    }

}

//get the G1 generator in affine mode
void G1_affine_generator(uint64_t* x,uint64_t* y){
    one(x,R_q); //G1_GENERATOR_X
    uint64_t val[4]={2,0,0,0};
    from_raw(val,y,R2_q,INV_q,MODULUS_q); //G1_GENERATOR_Y
}

//This is a simple double-and-add implementation of point
//multiplication, moving from most significant to least
//significant bit of the scalar.
void PMULT(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* scalar,uint64_t* acc_x,uint64_t* acc_y,uint64_t* acc_z){
    identity(acc_x,acc_y,acc_z);
    uint8_t scalar_byte[32];
    to_repr(scalar,scalar_byte,INV_r,MODULUS_r);
    for(int i=31;i>=0;i--){
        uint8_t byte=scalar_byte[i];
        for(int j=7;j>=0;j--){
            bool bit=(byte>>j)&1;
            point_double(acc_x,acc_y,acc_z);
            if(bit){
                project_add_project(acc_x,acc_y,acc_z,x,y,z);
            }
        }
    }
}

//PMULT where P is in affine mode
void Affine_PMULT(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* scalar,uint64_t* acc_x,uint64_t* acc_y,uint64_t* acc_z){
    identity(acc_x,acc_y,acc_z);
    uint8_t scalar_byte[32];
    to_repr(scalar,scalar_byte,INV_r,MODULUS_r);
    for(int i=31;i>=0;i--){
        uint8_t byte=scalar_byte[i];
        for(int j=7;j>=0;j--){
            bool bit=(byte>>j)&1;
            point_double(acc_x,acc_y,acc_z);
            if(bit){
                project_add_affine(acc_x,acc_y,acc_z,x,y,z);
            }
        }
    }
}

//get the G2 generator in affine mode
void G2_affine_generator(uint64_t* x,uint64_t* y){
    uint64_t x_c0[4];
    uint64_t x_c1[4];
    uint64_t y_c0[4];
    uint64_t y_c1[4];
    
    uint64_t val_x_c0[4]={
        0x46debd5cd992f6ed,
        0x674322d4f75edadd,
        0x426a00665e5c4479,
        0x1800deef121f1e76};
    uint64_t val_x_c1[4]={
        0x97e485b7aef312c2,
        0xf1aa493335a9e712,
        0x7260bfb731fb5d25,
        0x198e9393920d483a};
    uint64_t val_y_c0[4]={
        0x4ce6cc0166fa7daa,
        0xe3d1e7690c43d37b,
        0x4aab71808dcb408f,
        0x12c85ea5db8c6deb};
    uint64_t val_y_c1[4]={
        0x55acdadcd122975b,
        0xbc4b313370b38ef3,
        0xec9e99ad690c3395,
        0x090689d0585ff075};
    
    from_raw(val_x_c0,x_c0,R2_q,INV_q,MODULUS_q);
    from_raw(val_x_c1,x_c1,R2_q,INV_q,MODULUS_q);
    from_raw(val_y_c0,y_c0,R2_q,INV_q,MODULUS_q);
    from_raw(val_y_c1,y_c1,R2_q,INV_q,MODULUS_q);
    u64_to_u64(x,x_c0);
    u64_to_u64(x+4,x_c1);
    u64_to_u64(y,y_c0);
    u64_to_u64(y+4,y_c1);
    
}

//G2 affine to projective 
void G2_to_curve(uint64_t* x,uint64_t* y,uint64_t* z){
    bool ct1=is_identity(x,y);
    bool ct2=is_identity(x+4,y+4);

    if (ct1&ct2){
        fq2_zero(z);
    }
    else{
        fq2_one(z,R_q);
    }
}

//projection to affine
void to_affine(uint64_t* x,uint64_t* y,uint64_t* z){
    uint64_t zinv[4];
    invert(z,exp_q,zinv,R_q,INV_q,MODULUS_q);
    if((zinv[0]||zinv[1]||zinv[2]||zinv[3])==0){
        affine_identity(x,y);
    }
    else
    {
    uint64_t zinv2[4];
    uint64_t zinv3[4];
    SQUARE(zinv,zinv2,INV_q,MODULUS_q);
    MUL(x,zinv2,x,INV_q,MODULUS_q);
    MUL(zinv2,zinv,zinv3,INV_q,MODULUS_q);
    MUL(y,zinv3,y,INV_q,MODULUS_q);
    }
}

//G2 projection to affine
void G2_to_affine(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* affine_x,uint64_t* affine_y){
    uint64_t zinv[8];
    fq2_invert(z,exp_q,zinv,R_q,INV_q,MODULUS_q);
    if((zinv[0]||zinv[1]||zinv[2]||zinv[3]||zinv[4]||zinv[5]||zinv[6]||zinv[7])==0){
        affine_identity(affine_x,affine_y);
        affine_identity(affine_x+4,affine_y+4);
    }
    else
    {
    uint64_t zinv2[8];
    uint64_t zinv3[8];
    fq2_square(zinv,zinv2,INV_q,MODULUS_q);
    fq2_mul(zinv2,x,affine_x,INV_q,MODULUS_q);
    fq2_mul(zinv,zinv2,zinv3,INV_q,MODULUS_q);
    fq2_mul(zinv3,y,affine_y,INV_q,MODULUS_q);
    }
}

//G2 point double in Jacobian form: 2(x,y,z)
void G2_point_double(uint64_t* x,uint64_t* y,uint64_t* z){
    uint64_t a[8],b[8],c[8],d[8],e[8],f[8];
    fq2_square(x,a,INV_q,MODULUS_q);
    fq2_square(y,b,INV_q,MODULUS_q);
    fq2_square(b,c,INV_q,MODULUS_q);
    fq2_add(x,b,d,MODULUS_q);
    fq2_square(d,d,INV_q,MODULUS_q);
    fq2_sub(d,a,d,MODULUS_q);
    fq2_sub(d,c,d,MODULUS_q);
    fq2_add(d,d,d,MODULUS_q);
    fq2_add(a,a,e,MODULUS_q);
    fq2_add(e,a,e,MODULUS_q);
    fq2_square(e,f,INV_q,MODULUS_q);
    
    uint64_t x3[8],y3[8],z3[8];
    uint64_t mid1[8],mid2[8]; //存储中间变量
    fq2_mul(y,z,z3,INV_q,MODULUS_q);
    fq2_add(z3,z3,z3,MODULUS_q);
    fq2_add(d,d,mid1,MODULUS_q);
    fq2_sub(f,mid1,x3,MODULUS_q);
    fq2_add(c,c,c,MODULUS_q);
    fq2_add(c,c,c,MODULUS_q);
    fq2_add(c,c,c,MODULUS_q);
    fq2_sub(d,x3,mid1,MODULUS_q);
    fq2_mul(mid1,e,mid2,INV_q,MODULUS_q);
    fq2_sub(mid2,c,y3,MODULUS_q);

    bool ct1=is_identity_project(x,y,z);
    bool ct2=is_identity_project(x+4,y+4,z+4);
    bool ct=ct1&ct2;
    for(int i=0;i<8;i++){
        x[i]=(ct)?0:x3[i];
        y[i]=(ct)?0:y3[i];
        z[i]=(ct)?0:z3[i];
    }

}

//projective G2 point add affine point : (x,y,z)+(X,Y) and convert to projection
void G2_project_add_affine(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* X,uint64_t* Y,uint64_t* Z){
    if(is_identity_project(x,y,z)&is_identity_project(x+4,y+4,z+4)){
        G2_to_curve(X,Y,Z);
        u64_to_u64(x,X);
        u64_to_u64(x+4,X+4);
        u64_to_u64(y,Y);
        u64_to_u64(y+4,Y+4);
        u64_to_u64(z,Z);
        u64_to_u64(z+4,Z+4);
    }
    else if (is_identity(X,Y)&is_identity(X+4,Y+4))
    {
        return;
    }
    else{
        uint64_t z1z1[8],u2[8],s2[8];
        fq2_square(z,z1z1,INV_q,MODULUS_q);
        fq2_mul(z1z1,X,u2,INV_q,MODULUS_q);
        fq2_mul(z1z1,Y,s2,INV_q,MODULUS_q);
        fq2_mul(z,s2,s2,INV_q,MODULUS_q);

        if(is_equal(x,u2)&is_equal(x+4,u2+4)){
            if(is_equal(y,s2)&is_equal(y+4,s2+4)){
                G2_point_double(x,y,z);
            }
            else{
                identity(x,y,z);
                identity(x+4,y+4,z+4);
            }
        }
        else
        {
            uint64_t h[8],hh[8],i[8],j[8],r[8],v[8];
            fq2_sub(u2,x,h,MODULUS_q);
            fq2_square(h,hh,INV_q,MODULUS_q);
            fq2_add(hh,hh,i,MODULUS_q);
            fq2_add(i,i,i,MODULUS_q);
            fq2_mul(i,h,j,INV_q,MODULUS_q);
            fq2_sub(s2,y,r,MODULUS_q);
            fq2_add(r,r,r,MODULUS_q);
            fq2_mul(i,x,v,INV_q,MODULUS_q);

            uint64_t mid1[8],mid2[8];
            uint64_t x3[8],y3[8],z3[8];
            fq2_square(r,mid1,INV_q,MODULUS_q);
            fq2_sub(mid1,j,mid1,MODULUS_q);
            fq2_sub(mid1,v,mid2,MODULUS_q);
            fq2_sub(mid2,v,x3,MODULUS_q);

            fq2_mul(j,y,j,INV_q,MODULUS_q);
            fq2_add(j,j,j,MODULUS_q);
            fq2_sub(v,x3,mid1,MODULUS_q);
            fq2_mul(r,mid1,mid2,INV_q,MODULUS_q);
            fq2_sub(mid2,j,y3,MODULUS_q);

            fq2_add(z,h,mid1,MODULUS_q);
            fq2_square(mid1,mid1,INV_q,MODULUS_q);
            fq2_sub(mid1,z1z1,mid2,MODULUS_q);
            fq2_sub(mid2,hh,z3,MODULUS_q);

            u64_to_u64(x,x3);
            u64_to_u64(x+4,x3+4);
            u64_to_u64(y,y3);
            u64_to_u64(y+4,y3+4);
            u64_to_u64(z,z3);
            u64_to_u64(z+4,z3+4);
        }
    }
    
}

//PMULT where P is a G2 point and in affine mode
void G2_Affine_PMULT(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* scalar,uint64_t* acc_x,uint64_t* acc_y,uint64_t* acc_z){
    identity(acc_x,acc_y,acc_z);
    identity(acc_x+4,acc_y+4,acc_z+4);
    uint8_t scalar_byte[32];
    to_repr(scalar,scalar_byte,INV_r,MODULUS_r);
    for(int i=31;i>=0;i--){
        uint8_t byte=scalar_byte[i];
        for(int j=7;j>=0;j--){
            bool bit=(byte>>j)&1;
            G2_point_double(acc_x,acc_y,acc_z);
            if(bit){
                G2_project_add_affine(acc_x,acc_y,acc_z,x,y,z);
            }
        }
    }
}

void SETUP(uint32_t k, uint64_t** g_x, uint64_t** g_y, uint64_t** g_lagrange_x, uint64_t** g_lagrange_y, uint64_t g2_x[8], uint64_t g2_y[8], uint64_t s_g2_x[8], uint64_t s_g2_y[8])
{
    uint64_t n=1<<k;
    uint64_t g1_x[4];
    uint64_t g1_y[4];
    uint64_t g1_z[4];
    uint64_t s[4]; //toxic value
    G1_affine_generator(g1_x,g1_y);
    to_curve(g1_x,g1_y,g1_z);
    //Random(s,R2_q,R3_q,INV_q,MODULUS_q);
    uint64_t ROOT_OF_UNITY_RAW[4]={
    0xd34f1ed960c37c9c,
    0x3215cf6dd39329c8,
    0x98865ea93dd31f74,
    0x03ddb9f5166d18b7};
    from_raw(ROOT_OF_UNITY_RAW,s,R2_r,INV_r,MODULUS_r);
    // Calculate g = [G1, [s] G1, [s^2] G1, ..., [s^(n-1)] G1] in parallel.
    // Calculate g's Jacobian Projection coordinate 
    uint64_t** g_projective_x; 
    uint64_t** g_projective_y;
    uint64_t** g_projective_z;
    g_projective_x=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    g_projective_y=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    g_projective_z=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    for (uint64_t i=0; i<n; i++)
    {
        //为每列分配4个uint64_t大小的空间
        g_projective_x[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
        g_projective_y[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
        g_projective_z[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
        
    }

    uint64_t current_g_x[4];
    uint64_t current_g_y[4];
    uint64_t current_g_z[4];
    uint64_t current_s[4];
    u64_to_u64(g_projective_x[0],g1_x);
    u64_to_u64(g_projective_y[0],g1_y);
    u64_to_u64(current_g_x,g1_x);
    u64_to_u64(current_g_y,g1_y);
    u64_to_u64(current_s,s);
    to_curve(g_projective_x[0],g_projective_y[0],g_projective_z[0]);
    to_curve(current_g_x,current_g_y,current_g_z);
    for (uint64_t i=1; i<n; i++)
    { 
        u64_to_u64(current_g_x,g1_x);
        u64_to_u64(current_g_y,g1_y);
        to_curve(current_g_x,current_g_y,current_g_z);
        PMULT(current_g_x,current_g_y,current_g_z,current_s,g_projective_x[i],g_projective_y[i],g_projective_z[i]);
        MUL(current_s,s,current_s,INV_r,MODULUS_r);
    }
    // turn Projection coordinate to Affine coordinate
    for (uint64_t i=0;i<n;i++){
        batch_normalize(g_projective_x[i],g_projective_y[i],g_projective_z[i],g_x[i],g_y[i]);
    }

    //Calculate lagrange basis function point over an order-n multiplicative subgroup 
    //X coordinate are all n^th root of unity
    uint64_t** g_projective_lagrange_x; 
    uint64_t** g_projective_lagrange_y;
    uint64_t** g_projective_lagrange_z;
    g_projective_lagrange_x=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    g_projective_lagrange_y=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    g_projective_lagrange_z=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    for (uint64_t i=0; i<n; i++)
    {
        //为每列分配4个uint64_t大小的空间
        g_projective_lagrange_x[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
        g_projective_lagrange_y[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
        g_projective_lagrange_z[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
        
    }

    uint64_t ROOT_OF_UNITY_INV_r[4];
    uint64_t root[4];
    /// 1 / ROOT_OF_UNITY mod r
    uint64_t ROOT_OF_UNITY_INV_r_RAW[4]={
    0x0ed3e50a414e6dba,
    0xb22625f59115aba7,
    0x1bbe587180f34361,
    0x048127174daabc26
    };
    from_raw(ROOT_OF_UNITY_INV_r_RAW,ROOT_OF_UNITY_INV_r,R2_r,INV_r,MODULUS_r);
    invert(ROOT_OF_UNITY_INV_r,exp_r,root,R_r,INV_r,MODULUS_r); //这个地方可以优化成直接预存root

    //Calculate one n^th root of unity w and use it as the generator of the order-n multiplicative subgroup
    //S=28
    for (int i=k; i<28;i++)
    {
        SQUARE(root,root,INV_r,MODULUS_r);
    }
    //n_inv correspond to the barycentric weight c_i when i=0
    uint64_t n_inv[4];  //这个也可以预存
    from(n,n_inv,R2_r,INV_r,MODULUS_r);
    invert(n_inv,exp_r,n_inv,R_r,INV_r,MODULUS_r);
    //multiplier correspond to c_i*(x^n-1) where x=s and basis func L_i(x)=multiplier*w^i/(x-w^i) 
    uint64_t multiplier[4]; //这个也可以预存
    pow_vartime(s,n,multiplier,R_r,INV_r,MODULUS_r);
    uint64_t ONE[4];
    one(ONE,R_r);
    SUB(multiplier,ONE,multiplier,MODULUS_r);
    MUL(multiplier,n_inv,multiplier,INV_r,MODULUS_r);

    for(uint64_t i=0;i<n;i++)
    {
        uint64_t root_pow[4];
        pow_vartime(root,i,root_pow,R_r,INV_r,MODULUS_r); //compute w^i
        uint64_t scalar[4];
        uint64_t denominator[4];
        uint64_t numerator[4];
        SUB(s,root_pow,denominator,MODULUS_r);
        invert(denominator,exp_r,denominator,R_r,INV_r,MODULUS_r);
        MUL(multiplier,root_pow,numerator,INV_r,MODULUS_r);
        MUL(numerator,denominator,scalar,INV_r,MODULUS_r); //compute L_i(x)
        Affine_PMULT(g1_x,g1_y,g1_z,scalar,g_projective_lagrange_x[i],g_projective_lagrange_y[i],g_projective_lagrange_z[i]);
    }

    // turn Projection coordinate to Affine coordinate
    for (uint64_t i=0;i<n;i++){
        batch_normalize(g_projective_lagrange_x[i],g_projective_lagrange_y[i],g_projective_lagrange_z[i],g_lagrange_x[i],g_lagrange_y[i]);
    } 


    uint64_t g2_z[8];
    uint64_t g2_projective_x[8];
    uint64_t g2_projective_y[8];
    uint64_t g2_projective_z[8];

    G2_affine_generator(g2_x,g2_y);
    G2_Affine_PMULT(g2_x,g2_y,g2_z,s,g2_projective_x,g2_projective_y,g2_projective_z);
    G2_to_affine(g2_projective_x,g2_projective_y,g2_projective_z,s_g2_x,s_g2_y);
}