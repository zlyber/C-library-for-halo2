#include "fields.h"

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

void SUB(uint64_t* self, const uint64_t* rhs,uint64_t* result,const uint64_t* MODULUS) {
        uint64_t borrow=0;
        uint64_t d0,d1,d2,d3;
        sbb(self[0], rhs[0], borrow, &d0, &borrow);
        sbb(self[1], rhs[1], borrow, &d1, &borrow);
        sbb(self[2], rhs[2], borrow, &d2, &borrow);
        sbb(self[3], rhs[3], borrow, &d3, &borrow);

        // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
        // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the modulus.
        uint64_t carry=0;
        uint64_t Ins;   //insignificant...
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
        SUB(f, MODULUS, result,MODULUS);
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



// Schoolbook multiplication
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

/// Returns one, the multiplicative identity.
void one(uint64_t* res,const uint64_t* R){
    for(int i=0;i<4;i++){
        res[i]=R[i];
    }
}

void from(uint64_t val, uint64_t* result, const uint64_t* R2, const uint64_t INV, const uint64_t* MODULUS) {
    uint64_t f[4]={val, 0, 0, 0};   
    MUL(f, R2, result, INV, MODULUS);
}

//rescale to Mongomery form
void from_raw(uint64_t* val, uint64_t* result, const uint64_t* R2, const uint64_t INV, const uint64_t* MODULUS){
    MUL(val, R2, result, INV, MODULUS);
}

/// We reduce an arbitrary 512-bit number by decomposing it into two 256-bit digits
// with the higher bits multiplied by 2^256. Thus, we perform two reductions
//
// 1. the lower bits are multiplied by R^2, as normal
// 2. the upper bits are multiplied by R^2 * 2^256 = R^3
//
// and computing their sum in the field. It remains to see that arbitrary 256-bit
// numbers can be placed into Montgomery form safely using the reduction. The
// reduction works so long as the product is less than R=2^256 multiplied by
// the modulus. This holds because for any `c` smaller than the modulus, we have
// that (2^256 - 1)*c is an acceptable product for the reduction. Therefore, the
// reduction always works so long as `c` is in the field; in this case it is either the
// constant `R2` or `R3`.
void from_u512(uint64_t arr[8], uint64_t* res, const uint64_t*R2, const uint64_t*R3, const uint64_t INV, const uint64_t* MODULUS){
    uint64_t d0[4]={arr[0],arr[1],arr[2],arr[3]};
    uint64_t d1[4]={arr[4],arr[5],arr[6],arr[7]};
    
    MUL(d0,R2,d0,INV,MODULUS);
    MUL(d1,R3,d1,INV,MODULUS);
    ADD(d0,d1,res,MODULUS);
}

/// Converts from an integer represented in little endian
/// into its (congruent) `$field` representation.
// pub const fn from_raw(val: [u64; 4]) -> Self {
//     let (r0, carry) = mac(0, val[0], $r2.0[0], 0);
//     let (r1, carry) = mac(0, val[0], $r2.0[1], carry);
//     let (r2, carry) = mac(0, val[0], $r2.0[2], carry);
//     let (r3, r4) = mac(0, val[0], $r2.0[3], carry);

//     let (r1, carry) = mac(r1, val[1], $r2.0[0], 0);
//     let (r2, carry) = mac(r2, val[1], $r2.0[1], carry);
//     let (r3, carry) = mac(r3, val[1], $r2.0[2], carry);
//     let (r4, r5) = mac(r4, val[1], $r2.0[3], carry);

//     let (r2, carry) = mac(r2, val[2], $r2.0[0], 0);
//     let (r3, carry) = mac(r3, val[2], $r2.0[1], carry);
//     let (r4, carry) = mac(r4, val[2], $r2.0[2], carry);
//     let (r5, r6) = mac(r5, val[2], $r2.0[3], carry);

//     let (r3, carry) = mac(r3, val[3], $r2.0[0], 0);
//     let (r4, carry) = mac(r4, val[3], $r2.0[1], carry);
//     let (r5, carry) = mac(r5, val[3], $r2.0[2], carry);
//     let (r6, r7) = mac(r6, val[3], $r2.0[3], carry);

//     // Montgomery reduction (first part)
//     let k = r0.wrapping_mul($inv);
//     let (_, carry) = mac(r0, k, $modulus.0[0], 0);
//     let (r1, carry) = mac(r1, k, $modulus.0[1], carry);
//     let (r2, carry) = mac(r2, k, $modulus.0[2], carry);
//     let (r3, carry) = mac(r3, k, $modulus.0[3], carry);
//     let (r4, carry2) = adc(r4, 0, carry);

//     let k = r1.wrapping_mul($inv);
//     let (_, carry) = mac(r1, k, $modulus.0[0], 0);
//     let (r2, carry) = mac(r2, k, $modulus.0[1], carry);
//     let (r3, carry) = mac(r3, k, $modulus.0[2], carry);
//     let (r4, carry) = mac(r4, k, $modulus.0[3], carry);
//     let (r5, carry2) = adc(r5, carry2, carry);

//     let k = r2.wrapping_mul($inv);
//     let (_, carry) = mac(r2, k, $modulus.0[0], 0);
//     let (r3, carry) = mac(r3, k, $modulus.0[1], carry);
//     let (r4, carry) = mac(r4, k, $modulus.0[2], carry);
//     let (r5, carry) = mac(r5, k, $modulus.0[3], carry);
//     let (r6, carry2) = adc(r6, carry2, carry);

//     let k = r3.wrapping_mul($inv);
//     let (_, carry) = mac(r3, k, $modulus.0[0], 0);
//     let (r4, carry) = mac(r4, k, $modulus.0[1], carry);
//     let (r5, carry) = mac(r5, k, $modulus.0[2], carry);
//     let (r6, carry) = mac(r6, k, $modulus.0[3], carry);
//     let (r7, _) = adc(r7, carry2, carry);

//     // Montgomery reduction (sub part)
//     let (d0, borrow) = sbb(r4, $modulus.0[0], 0);
//     let (d1, borrow) = sbb(r5, $modulus.0[1], borrow);
//     let (d2, borrow) = sbb(r6, $modulus.0[2], borrow);
//     let (d3, borrow) = sbb(r7, $modulus.0[3], borrow);

//     let (d0, carry) = adc(d0, $modulus.0[0] & borrow, 0);
//     let (d1, carry) = adc(d1, $modulus.0[1] & borrow, carry);
//     let (d2, carry) = adc(d2, $modulus.0[2] & borrow, carry);
//     let (d3, _) = adc(d3, $modulus.0[3] & borrow, carry);

//     $field([d0, d1, d2, d3])

//     &$field(val)).mul(&$r2)

// }

//get a random element in field
void Random(uint64_t* res, const uint64_t*R2, const uint64_t*R3,const uint64_t INV, const uint64_t* MODULUS){
    uint64_t arr[8];
    for(int i=0;i<8;i++){
        arr[i]=next_u64();
    }
    return from_u512(arr, res, R2, R3, INV, MODULUS);
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

/// Computes the multiplicative inverse of this element,
/// failing if the element is zero.
void invert(uint64_t* self, const uint64_t* exp, uint64_t* res, const uint64_t* R, const uint64_t INV, const uint64_t* MODULUS)
 {
    uint64_t tmp[4];
    pow(self,exp,tmp,R,INV,MODULUS);
    if((self[0]||self[1]||self[2]||self[3])!=0){
        u64_to_u64(res,tmp);
    }
    else{
        printf("inversion not exist\n");
    }
}

/// Exponentiates `self` by `by`, where `by` is a little-endian order
/// integer exponent.
void pow(uint64_t* self, const uint64_t* by, uint64_t* res, const uint64_t* R, const uint64_t INV, const uint64_t* MODULUS) 
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