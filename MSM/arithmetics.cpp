#include "arithmetics.h"



uint64_t get_at(int segment, int c, uint8_t* bytes) {
    int skip_bits = segment * c;
    int skip_bytes = skip_bits / 8;

    if (skip_bytes >= 32) {
        return 0;
    }

    uint8_t v[8] = {0};
    memcpy(v, bytes + skip_bytes, 8);

    uint64_t tmp = *(uint64_t*)v;
    tmp >>= skip_bits - (skip_bytes * 8);
    tmp = tmp % (1 << c);

    return tmp;
}

//judge whether the bucket is null now
bool is_null(uint64_t* x){
    bool is_one=1;
    for(int i=0;i<4;i++){
        is_one=is_one&(x[i]==1);
    }
    return is_one;
}

//MSM
void multiexp_serial(uint64_t** coeffs, uint64_t** bases_x, uint64_t** bases_y, uint64_t** bases_z, uint64_t* acc_x, uint64_t* acc_y, uint64_t* acc_z, uint64_t len){
    uint8_t** coeffs_repr;
    coeffs_repr = (uint8_t**)malloc(sizeof(uint8_t*)*len);//为二维数组分配行
    for (uint64_t i=0; i<len; i++)
	{
		//为每列分配32个u8大小的空间
		coeffs_repr[i] = (uint8_t*)malloc(sizeof(uint8_t)*32); 
	} 
    for (uint64_t i = 0; i < len; i++) {
        to_repr(coeffs[i],coeffs_repr[i]);
    }

    int c;
    if(len<4){
        c=1;
    }
    else if(len<32){
        c=3;
    }
    else{
        c = (int)ceil(log(len));
    }
    
    int segments = (256 / c) + 1;

    for (int current_segment = segments - 1; current_segment >= 0; current_segment--) {
        for (int i = 0; i < c; i++) {
            point_double(acc_x,acc_y,acc_z); 
        }

        //为buckets创建空间，要创建三个buckets来代表x,y,z三维
        uint64_t** buckets_x; 
        uint64_t** buckets_y;
        uint64_t** buckets_z;
        uint64_t bucket_len=(1<<c)-1;
        buckets_x=(uint64_t**)malloc(sizeof(uint64_t*)*bucket_len);
        buckets_y=(uint64_t**)malloc(sizeof(uint64_t*)*bucket_len);
        buckets_z=(uint64_t**)malloc(sizeof(uint64_t*)*bucket_len);
        for (uint64_t i=0; i<bucket_len; i++)
        {
           //为每列分配4个uint64_t大小的空间，并将初始值都设为0
           buckets_x[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
           buckets_y[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 
    
           buckets_z[i] = (uint64_t*)malloc(sizeof(uint64_t)*4); 

           //初始化为1，视为空
           for(int j=0;j<4;j++){
            buckets_x[i][j]=1;
            buckets_y[i][j]=1;
            buckets_z[i][j]=1;
           }
        }
        
        for (uint64_t i = 0; i < len; i++) {          
            uint64_t coeff = get_at(current_segment, c, coeffs_repr[i]);
    
            //以下分支达到和add_assign函数一样的效果
            if (coeff != 0) {
                //若当前桶为空，则先把坐标值放进去
                if(is_null(buckets_x[coeff - 1])&is_null(buckets_y[coeff - 1])&is_null(buckets_z[coeff - 1]))
                {
                    u64_to_u64(buckets_x[coeff - 1],bases_x[i]);
                    u64_to_u64(buckets_y[coeff - 1],bases_y[i]);
                }
                //当前是仿射坐标形式，使用仿射坐标相加后转换到投影坐标
                else if(is_null(buckets_z[coeff - 1]))
                {
                    affine_add_affine(buckets_x[coeff - 1],buckets_y[coeff - 1],buckets_z[coeff - 1], bases_x[i],bases_y[i],bases_z[i]);
        
                }
                //当前已经是投影坐标形式，使用投影坐标加仿射坐标，之后转换到仿射坐标
                else 
                {
                    project_add_affine(buckets_x[coeff - 1],buckets_y[coeff - 1],buckets_z[coeff - 1], bases_x[i],bases_y[i],bases_z[i]);
        
                }
            }
        }
        
        // Summation by parts
        // e.g. 3a + 2b + 1c = a +
        //                    (a) + b +
        //                    ((a) + b) + c
        uint64_t sum_x[4],sum_y[4],sum_z[4];
        identity(sum_x,sum_y,sum_z); 
        for (int i = bucket_len - 1; i >= 0; i--) {
            //以下分支达到和add函数一样的效果，和add_assign的区别是add函数里sum是主体,buckets是客体
            //若当前桶为空，则跳过不管
            if(is_null(buckets_x[i])&is_null(buckets_y[i])&is_null(buckets_z[i]))
            {
                project_add_project(acc_x, acc_y, acc_z, sum_x, sum_y, sum_z);
                continue;
            }
            //若当前桶内是仿射坐标形式，使用投影坐标加仿射坐标的方法
            else if(is_null(buckets_z[i]))
            {
                project_add_affine(sum_x, sum_y, sum_z, buckets_x[i],buckets_y[i],buckets_z[i]);
            }
            //若当前已经是投影坐标形式，直接使用两个投影坐标相加的方式
            else 
            {
                project_add_project(sum_x, sum_y, sum_z, buckets_x[i],buckets_y[i],buckets_z[i]);
            }
            //*acc = *acc + &running_sum;
            project_add_project(acc_x, acc_y, acc_z, sum_x, sum_y, sum_z);
            
        }
    }
}