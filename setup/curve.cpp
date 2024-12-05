#include "curve.h"

//return the identity of the curve in projective form
void identity(uint64_t* x,uint64_t* y,uint64_t* z){
    for(int i=0;i<4;i++){
        x[i]=0;
        y[i]=0;
        z[i]=0;
    }
}

void affine_identity(uint64_t* x,uint64_t* y){
    for(int i=0;i<4;i++){
        x[i]=0;
        y[i]=0;
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

// This is a simple double-and-add implementation of point
// multiplication, moving from most significant to least
// significant bit of the scalar.
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

// PMULT where P is in affine mode
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

// G2 affine to projective 
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

// PMULT where P is a G2 point and in affine mode
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
            printf("acc_x_c1[0]=%llu\n",acc_x[4]);
            printf("acc_y_c1[0]=%llu\n",acc_y[4]);
            printf("acc_z_c1[0]=%llu\n",acc_z[4]);
            if(bit){
                G2_project_add_affine(acc_x,acc_y,acc_z,x,y,z);
                printf("acc_x_c1[0]=%llu\n",acc_x[4]);
                printf("acc_y_c1[0]=%llu\n",acc_y[4]);
                printf("acc_z_c1[0]=%llu\n",acc_z[4]);
            }
        }
    }
}

// void G2_batch_normalize(uint64_t* x,uint64_t* y,uint64_t* z,uint64_t* affine_x,uint64_t* affine_y){
//     uint64_t acc[8];
//     fq2_one(acc,R_q);
//     u64_to_u64(affine_x,acc);
//     u64_to_u64(affine_x+4,acc+4);
//     if(! ((is_identity_project(x,y,z)) & (is_identity_project(x+4,y+4,z+4))))
//     {
//         fq2_mul(z,acc,acc,INV_q,MODULUS_q);
//     }
//     uint64_t pow[4]={
//     0x3c208c16d87cfd45,
//     0x97816a916871ca8d,
//     0xb85045b68181585d,
//     0x30644e72e131a029};
//     invert(acc,pow,acc,R_q,INV_q,MODULUS_q);
//     uint64_t tmp[4];
//     MUL(affine_x,acc,tmp,INV_q,MODULUS_q);
//     if(!(is_identity_project(x,y,z))){
//         MUL(acc,z,acc,INV_q,MODULUS_q);
//     }
//     uint64_t tmp2[4];
//     uint64_t tmp3[4];
//     SQUARE(tmp,tmp2,INV_q,MODULUS_q);
//     MUL(tmp2,tmp,tmp3,INV_q,MODULUS_q);
//     if(!(is_identity_project(x,y,z))){
//     MUL(x,tmp2,affine_x,INV_q,MODULUS_q);
//     MUL(y,tmp3,affine_y,INV_q,MODULUS_q);
//     }
//     else{
//         affine_identity(affine_x,affine_y);
//     }

// }