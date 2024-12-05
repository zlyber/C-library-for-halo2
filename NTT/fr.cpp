#include "fr-bn256.h"


void SUB(uint64_t* self, const uint64_t* rhs,uint64_t* result) {
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

void ADD(uint64_t* self, const uint64_t* rhs, uint64_t* result) {
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
        SUB(f, MODULUS, result);
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
        uint64_t* result
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
        SUB(f, MODULUS, result);
    }



// Schoolbook multiplication
void MUL(uint64_t* self, const uint64_t* rhs, uint64_t* result) {
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

    montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7, result);
    //montgomery_reduce(result[0],result[1],result[2],result[3], r4, r5, r6, r7, result);
}

/// Squares this element.
void SQUARE(uint64_t* self) {
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

    montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7, self);
}

void from(uint64_t val,uint64_t* result) {
    uint64_t f[4]={val, 0, 0, 0};   
    MUL(f, R2, result);
}

/// Returns zero, the additive identity.
// uint64_t* zero(){
//     uint64_t zero[4]={0,0,0,0};
//     return zero;
// }

/// Returns one, the multiplicative identity.
void one(uint64_t* res){
    for(int i=0;i<4;i++){
        res[i]=R[i];
    }
}

//rescale to Mongomery form
void from_raw(uint64_t* val, uint64_t* result){
    MUL(val, R2, result);
}


void from_u512(uint64_t arr[8], uint64_t* res){
    uint64_t d0[4]={arr[0],arr[1],arr[2],arr[3]};
    uint64_t d1[4]={arr[4],arr[5],arr[6],arr[7]};
    MUL(d0,R2,d0);
    MUL(d1,R3,d1);
    ADD(d0,d1,res);
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
    from_raw((uint64_t*)ROU,result);
}

