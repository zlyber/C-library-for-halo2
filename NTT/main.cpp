#include "arithmetics.h"
#include <stdlib.h>
using namespace std;

void writeToFile(uint64_t** a, uint64_t N,char* filename)
{
    FILE* file = fopen(filename, "a+"); // 打开文件，以写入模式打开

    if (file == NULL) {
        printf("can't open the file\n");
        return;
    }

    for (uint64_t i = 0; i < N; i++) {
        for (uint64_t j = 0; j < 4; j++) {
            fprintf(file, "%llu ", a[i][j]); 
        }
        fprintf(file, "\n"); 
    }

    fclose(file); 
}

int main(){
    int k;
    printf("please input k:");
    scanf("%d",&k);
    uint64_t r_o_u[4];
    uint64_t omega[4];
    root_of_unity(r_o_u);
    precompute_omega(r_o_u, k, omega);
    uint64_t **a;
    uint64_t N=1<<k;
    a = (uint64_t**)malloc(sizeof(uint64_t*)*N);
    for (uint64_t i=0; i<N; i++)
	{
		a[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
        Random(a[i]);
	} 

    struct timeb t1,t2;
    long t;
    int i;

    ftime(&t1);  

    fft(a,omega,k);

    ftime(&t2); 
    t=(t2.time-t1.time)*1000+(t2.millitm-t1.millitm); 
    printf("NTT time is %ld ms\n",t);
    writeToFile(a, N, "NTT result.txt");

    for (uint64_t i = 0; i < N; i++){
        free(a[i]);
    }
    free(a);
    
}


