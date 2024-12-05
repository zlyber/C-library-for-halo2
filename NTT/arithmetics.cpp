#include "arithmetics.h"

using namespace std;

void u64_to_u64(uint64_t* u641, uint64_t* u642){
    for(int i=0;i<4;i++){
        u641[i]=u642[i];
    }
}

//swap 2 u64 array in a 2-D array
void swap(uint64_t* a, uint64_t* b) {
    uint64_t temp[4];
    u64_to_u64(temp,a);
    u64_to_u64(a,b);
    u64_to_u64(b,temp);
}

//return the bit_reverse of operand
uint64_t reverse_bits(uint64_t operand,int bit_count){

    uint64_t acc = 0;
    for (int i = 0; i < bit_count; i++) {
      acc = acc << 1;
      acc |= (operand >> i) & 1;
    }
    return acc;
}

void group_scale(uint64_t* self, uint64_t* rhs, uint64_t*res){
    MUL(self, rhs, res); 
}

void group_add(uint64_t* self, uint64_t* rhs, uint64_t* res){
    ADD(self, rhs, res);
}

void group_sub(uint64_t* self, uint64_t* rhs, uint64_t* res){
    SUB(self, rhs, res);
}

void scalar_one(uint64_t* result){
    one(result);
}

//precompute omega 实现src/poly/domain文件中new函数的部分功能
//N=2^k ,j=cs.degree()
void precompute_omega(uint64_t* root_of_unity,int k,int j,
int* extended_k, uint64_t* omega, uint64_t* extended_omega)
{
    const int S=32;
    uint64_t n=1ull<<k;
    // quotient_poly_degree * params.n - 1 is the degree of the quotient polynomial
    int quotient_poly_degree = j - 1;
    *extended_k = k;
    while ((1ull << *extended_k) < (n * quotient_poly_degree)) {
        *extended_k += 1;
    }
    // Get extended_omega, the 2^{extended_k}'th root of unity
    // The loop computes extended_omega = omega^{2 ^ (S - extended_k)}
    // Notice that extended_omega ^ {2 ^ extended_k} = omega ^ {2^S} = 1.
    u64_to_u64(extended_omega,root_of_unity);
    for (int i=*extended_k;i<S;i++) {
        SQUARE(extended_omega);
    }

    // Get omega, the 2^{k}'th root of unity (i.e. n'th root of unity)
    // The loop computes omega = extended_omega ^ {2 ^ (extended_k - k)}
    //           = (omega^{2 ^ (S - extended_k)})  ^ {2 ^ (extended_k - k)}
    //           = omega ^ {2 ^ (S - k)}.
    // Notice that omega ^ {2^k} = omega ^ {2^S} = 1.
    u64_to_u64(omega,extended_omega);
    for (int i=k;i<*extended_k;i++) {
        SQUARE(omega);
    }
    /*先直接填数据，计算过程之后再写

    uint64_t omega_inv[4];
    u64_to_u64(omega_inv,omega); // Inversion computed later

    uint64_t extended_omega_inv[4];
    u64_to_u64(extended_omega_inv,extended_omega);// Inversion computed later

    uint64_t ifft_divisor[4];
    uint64_t extended_ifft_divisor[4];
    Fp fp1,fp2;
    from(1 << k,&fp1);
    from(1 << extended_k,&fp2);

    //compute inversion
    
    Fp fp_omega_inv,fp_extended_omega_inv;
    u64_to_Fp(omega_inv,&fp_omega_inv);
    u64_to_Fp(extended_omega_inv,&fp_extended_omega_inv);
    fp_extended_omega_inv=invert(&fp_extended_omega_inv);
    fp_omega_inv=invert(&fp_omega_inv);
    invert(&fp1);
    invert(&fp2);
    Fp_to_u64(&fp_extended_omega_inv,extended_omega_inv);
    Fp_to_u64(&fp_omega_inv,omega_inv);
    Fp_to_u64(&fp1,ifft_divisor);
    Fp_to_u64(&fp2,extended_ifft_divisor);
    */
    
}

// precompute twiddle factors
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

//递归版本
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
        Left = (uint64_t**)malloc(sizeof(uint64_t*)*(N/2));//为二维数组分配n行
        Right = (uint64_t**)malloc(sizeof(uint64_t*)*(N/2));//为二维数组分配n行
        for (uint64_t i=0; i<N/2; i++)
	    {
		    //为每列分配4个uint64_t大小的空间
		    Left[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
            Right[i]= (uint64_t*)malloc(sizeof(uint64_t)*4); 
	    }
        
        for(uint64_t i=0;i<N/2;i++){
            u64_to_u64(Left[i],vector[i]);
            u64_to_u64(Right[i],vector[i+N/2]);
        }
        
        #pragma omp parallel sections
        {
            #pragma omp section
            recursive_butterfly(Left, N / 2, twiddle_chunk * 2, twiddles);
            #pragma omp section
            recursive_butterfly(Right, N / 2, twiddle_chunk * 2, twiddles);
        }
        uint64_t t[4];
        u64_to_u64(t,Right[0]);      //取出Right[0]
        u64_to_u64(Right[0],Left[0]);//Right[0]=Left[0]
        group_add(Left[0],t,Left[0]);
        group_sub(Right[0],t,Right[0]);
        for(uint64_t m=0;m<N/2-1;m++){
            uint64_t* Left_addr=Left[m+1];
            uint64_t* Right_addr=Right[m+1];
            uint64_t* twiddle=twiddles[(m+1)*twiddle_chunk];

            uint64_t t1[4];
            group_scale(Right_addr,twiddle,t1); //b*w
            u64_to_u64(Right_addr,Left_addr);
            group_add(Left_addr,t1,Left_addr); //a+b*w
            group_sub(Right_addr,t1,Right_addr); //a-b*w     
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

//非递归版本 k=log_n
void fft(uint64_t** vector, uint64_t* omega, int k){
    int threads=64;
    int log_threads = (int)log2(threads);
    uint64_t N=1<<k;

    for(uint64_t i=0; i<N; i++){
        uint64_t rk = reverse_bits(i, k);
        if (i < rk){
            swap(vector[rk],vector[i]);
        }
    }
    
    uint64_t** twiddles;
    twiddles = (uint64_t**)malloc(sizeof(uint64_t*)*(N/2));//为二维数组分配n行
    for (uint64_t i=0; i<N/2; i++)
	{
		//为每列分配4个uint64_t大小的空间
		twiddles[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
        
	} 
    precompute_twiddles(twiddles,omega,N);
   
    if(k<=log_threads){
        uint64_t chunk=2;
        uint64_t twiddle_chunk = N / 2;
        for(int i=0;i<k;i++){
            for(uint64_t j=0;j<N;j+=chunk){
                uint64_t* Left=&vector[j][0];
                uint64_t* Right=&vector[j+chunk/2][0];
                uint64_t t[4];
                u64_to_u64(t,Right);    //取出Right[0]
                u64_to_u64(Right,Left); //Right[0]=Left[0]
                group_add(Left,t,Left);
                group_sub(Right,t,Right);
                for(uint64_t m=0;m<chunk/2-1;m++){   //最底层，chunk=2时不进入该循环
                    uint64_t* Left_addr=vector[j+m+1];
                    uint64_t* Right_addr=vector[j+chunk/2+m+1];
                    uint64_t* twiddle=twiddles[(m+1)*twiddle_chunk];
                    uint64_t t1[4];
                    group_scale(Right_addr,twiddle,t1); //b*w
                    u64_to_u64(Right_addr,Left_addr);
                    group_add(Left_addr,t1,Left_addr); //a+b*w
                    group_sub(Right_addr,t1,Right_addr); //a-b*w     
                }  
            }
        chunk *= 2; // 向上合并
        twiddle_chunk /= 2;
        }
    }
    else{
        recursive_butterfly(vector,N,1,twiddles);
    }
    for (uint64_t i = 0; i < N/2; i++){
        free(twiddles[i]);
    }
    free(twiddles);
}



                
