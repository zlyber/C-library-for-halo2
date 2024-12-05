#include "curve.h"

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
