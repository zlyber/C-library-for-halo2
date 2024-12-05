#include "arithmetics.h"

#define ROWS 32
#define COLS 4


void read_file_into_array(const char* filename, uint64_t** array) {
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            fscanf(fp, "%llu", &array[i][j]);
        }
    }
    fclose(fp);
}

int main(){
    uint64_t** coeffs;
    uint64_t** base_x; 
    uint64_t** base_y;
    uint64_t** base_z;
    uint64_t base_len=ROWS;
    coeffs=(uint64_t**)malloc(sizeof(uint64_t*)*base_len);
    base_x=(uint64_t**)malloc(sizeof(uint64_t*)*base_len);
    base_y=(uint64_t**)malloc(sizeof(uint64_t*)*base_len);
    base_z=(uint64_t**)malloc(sizeof(uint64_t*)*base_len);
    for (uint64_t i=0; i<base_len; i++)
    {
        //为每列分配4个uint64_t大小的空间
        coeffs[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 

        base_x[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
        base_y[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
        base_z[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
        memset(base_z[i], 0, sizeof(uint64_t)*4);
    }
    read_file_into_array("coeffs.txt", coeffs);
    read_file_into_array("basex.txt", base_x);
    read_file_into_array("basey.txt", base_y);

    uint64_t acc_x[4],acc_y[4],acc_z[4];
    identity(acc_x,acc_y,acc_z);
    multiexp_serial(coeffs,base_x,base_y,base_z,acc_x,acc_y,acc_z,ROWS);
    //print the result
    printf("%llu\t%llu\t%llu\t%llu\n",acc_x[0],acc_x[1],acc_x[2],acc_x[3]);
    printf("%llu\t%llu\t%llu\t%llu\n",acc_y[0],acc_y[1],acc_y[2],acc_y[3]);
    printf("%llu\t%llu\t%llu\t%llu\n",acc_z[0],acc_z[1],acc_z[2],acc_z[3]);
}