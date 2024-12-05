#include"fq2-bn256.h"

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
