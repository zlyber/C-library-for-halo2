#include "arithmetics.h"
#include <stdlib.h>
using namespace std;

void writeToFile(uint64_t** a, uint64_t N,char* filename)
{
    FILE* file = fopen(filename, "w"); // 打开文件，以写入模式打开

    if (file == NULL) {
        printf("无法打开文件。\n");
        return;
    }

    for (uint64_t i = 0; i < N; i++) {
        for (uint64_t j = 0; j < 4; j++) {
            fprintf(file, "%llu ", a[i][j]); // 将数组元素写入文件，使用空格分隔
        }
        fprintf(file, "\n"); // 写入换行符，表示一行结束
    }

    fclose(file); // 关闭文件
}

int main(){
    int k;
    printf("please input k:");
    scanf("%d",&k);
    uint64_t omega24[4]={17104276958738319052, 7189382688254064745,1604580299050760569,2533986721453659612}; //2^24
    uint64_t omega12[4]={11791621636447142361, 3213422488462342693,7137044954843475233,1120461048903492910}; //2^12

    uint64_t **a;
    uint64_t N=1<<k;
    a = (uint64_t**)malloc(sizeof(uint64_t*)*N);//为二维数组分配n行
    for (uint64_t i=0; i<N; i++)
	{
		//为每列分配4个uint64_t大小的空间
		a[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
        Random(a[i]);
	} 
    for (int i=0;i<12;i++){
        SQUARE(omega12);
        printf("%llu,%llu,%llu,%llu\n",omega12[0],omega12[1],omega12[2],omega12[3]);
    }

    struct timeb t1,t2;
    long t;
    int i;

    ftime(&t1);  /* 获取当前时间 */

    fft(a,omega12,k);

    ftime(&t2); /* 获取当前时间 */
    t=(t2.time-t1.time)*1000+(t2.millitm-t1.millitm); /* 计算毫秒级的时间 */
    printf("NTT time is %ld ms\n",t);


    for (uint64_t i = 0; i < N; i++){
        free(a[i]);
    }
    free(a);
    
}


