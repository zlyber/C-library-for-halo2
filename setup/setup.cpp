#include "setup.h"

void setup(int k)
{
    uint64_t n=1<<k;
    FILE* fp;
    fp = fopen("g and g_lagrange.txt", "w");
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
    uint64_t** g_x; 
    uint64_t** g_y;
    g_x=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    g_y=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    for (uint64_t i=0; i<n; i++)
    {
        //为每列分配4个uint64_t大小的空间
        g_x[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
        g_y[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 

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
    fprintf(fp, "g_projective_x[0]=%llu,%llu,%llu,%llu\n",g_projective_x[0][0],g_projective_x[0][1],g_projective_x[0][2],g_projective_x[0][3]);
    fprintf(fp, "g_projective_y[0]=%llu,%llu,%llu,%llu\n",g_projective_y[0][0],g_projective_y[0][1],g_projective_y[0][2],g_projective_y[0][3]);
    fprintf(fp, "g_projective_z[0]=%llu,%llu,%llu,%llu\n",g_projective_z[0][0],g_projective_z[0][1],g_projective_z[0][2],g_projective_z[0][3]);
    for (uint64_t i=1; i<n; i++)
    { 
        u64_to_u64(current_g_x,g1_x);
        u64_to_u64(current_g_y,g1_y);
        to_curve(current_g_x,current_g_y,current_g_z);
        PMULT(current_g_x,current_g_y,current_g_z,current_s,g_projective_x[i],g_projective_y[i],g_projective_z[i]);
        MUL(current_s,s,current_s,INV_r,MODULUS_r);
        fprintf(fp, "g_projective_x[%llu]=%llu,%llu,%llu,%llu\n",i,g_projective_x[i][0],g_projective_x[i][1],g_projective_x[i][2],g_projective_x[i][3]);
        fprintf(fp, "g_projective_y[%llu]=%llu,%llu,%llu,%llu\n",i,g_projective_y[i][0],g_projective_y[i][1],g_projective_y[i][2],g_projective_y[i][3]);
        fprintf(fp, "g_projective_z[%llu]=%llu,%llu,%llu,%llu\n",i,g_projective_z[i][0],g_projective_z[i][1],g_projective_z[i][2],g_projective_z[i][3]);
        fprintf(fp,"\n");
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

    uint64_t** g_lagrange_x; 
    uint64_t** g_lagrange_y;
    g_lagrange_x=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    g_lagrange_y=(uint64_t**)malloc(sizeof(uint64_t*)*n);
    for (uint64_t i=0; i<n; i++)
    {
        //为每列分配4个uint64_t大小的空间
        g_lagrange_x[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
        g_lagrange_y[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 

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
        fprintf(fp, "g_projective__lagrange_x[%llu]=%llu,%llu,%llu,%llu\n",i,g_projective_lagrange_x[i][0],g_projective_lagrange_x[i][1],g_projective_lagrange_x[i][2],g_projective_lagrange_x[i][3]);
        fprintf(fp, "g_projective__lagrange_y[%llu]=%llu,%llu,%llu,%llu\n",i,g_projective_lagrange_y[i][0],g_projective_lagrange_y[i][1],g_projective_lagrange_y[i][2],g_projective_lagrange_y[i][3]);
        fprintf(fp,"\n");
    }

    // turn Projection coordinate to Affine coordinate
    for (uint64_t i=0;i<n;i++){
        batch_normalize(g_projective_lagrange_x[i],g_projective_lagrange_y[i],g_projective_lagrange_z[i],g_lagrange_x[i],g_lagrange_y[i]);
    }
    
    for (uint64_t i=0;i<n;i++){
        fprintf(fp, "g_x[%llu]=%llu,%llu,%llu,%llu\n",i,g_x[i][0],g_x[i][1],g_x[i][2],g_x[i][3]);
        fprintf(fp, "g_y[%llu]=%llu,%llu,%llu,%llu\n",i,g_y[i][0],g_y[i][1],g_y[i][2],g_y[i][3]);
        fprintf(fp,"\n");
    }
    fprintf(fp,"\n=====================\n");
    for (uint64_t i=0;i<n;i++){
        fprintf(fp, "g_lagrange_x[%llu]=%llu,%llu,%llu,%llu\n",i,g_lagrange_x[i][0],g_lagrange_x[i][1],g_lagrange_x[i][2],g_lagrange_x[i][3]);
        fprintf(fp, "g_lagrange_y[%llu]=%llu,%llu,%llu,%llu\n",i,g_lagrange_y[i][0],g_lagrange_y[i][1],g_lagrange_y[i][2],g_lagrange_y[i][3]);
        fprintf(fp,"\n");
    }
    

    uint64_t g2_x[8];
    uint64_t g2_y[8];
    uint64_t g2_z[8];
    uint64_t g2_projective_x[8];
    uint64_t g2_projective_y[8];
    uint64_t g2_projective_z[8];
    uint64_t s_g2_x[8];
    uint64_t s_g2_y[8];
    G2_affine_generator(g2_x,g2_y);
    G2_Affine_PMULT(g2_x,g2_y,g2_z,s,g2_projective_x,g2_projective_y,g2_projective_z);
    G2_to_affine(g2_projective_x,g2_projective_y,g2_projective_z,s_g2_x,s_g2_y);
    fclose(fp);
}