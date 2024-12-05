#include "fr.h"

void to_repr(uint64_t* self,uint8_t* res) {
    // Turn into canonical form by computing
    // (a.R) / R = a
    uint64_t tmp[4];
    montgomery_reduce(self[0],self[1],self[2],self[3],0,0,0,0,tmp,INV_r,MODULUS_r);

    res[0] = tmp[0] & 0xff;
    res[1] = (tmp[0] >> 8) & 0xff;
    res[2] = (tmp[0] >> 16) & 0xff;
    res[3] = (tmp[0] >> 24) & 0xff;
    res[4] = (tmp[0] >> 32) & 0xff;
    res[5] = (tmp[0] >> 40) & 0xff;
    res[6] = (tmp[0] >> 48) & 0xff;
    res[7] = (tmp[0] >> 56) & 0xff;
    res[8] = tmp[1] & 0xff;
    res[9] = (tmp[1] >> 8) & 0xff;
    res[10] = (tmp[1] >> 16) & 0xff;
    res[11] = (tmp[1] >> 24) & 0xff;
    res[12] = (tmp[1] >> 32) & 0xff;
    res[13] = (tmp[1] >> 40) & 0xff;
    res[14] = (tmp[1] >> 48) & 0xff;
    res[15] = (tmp[1] >> 56) & 0xff;
    res[16] = tmp[2] & 0xff;
    res[17] = (tmp[2] >> 8) & 0xff;
    res[18] = (tmp[2] >> 16) & 0xff;
    res[19] = (tmp[2] >> 24) & 0xff;
    res[20] = (tmp[2] >> 32) & 0xff;
    res[21] = (tmp[2] >> 40) & 0xff;
    res[22] = (tmp[2] >> 48) & 0xff;
    res[23] = (tmp[2] >> 56) & 0xff;
    res[24] = tmp[3] & 0xff;
    res[25] = (tmp[3] >> 8) & 0xff;
    res[26] = (tmp[3] >> 16) & 0xff;
    res[27] = (tmp[3] >> 24) & 0xff;
    res[28] = (tmp[3] >> 32) & 0xff;
    res[29] = (tmp[3] >> 40) & 0xff;
    res[30] = (tmp[3] >> 48) & 0xff;
    res[31] = (tmp[3] >> 56) & 0xff;
}
