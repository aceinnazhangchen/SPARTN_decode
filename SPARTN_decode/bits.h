#ifndef _BITS_H_
#define _BITS_H_

#ifdef __cplusplus
extern "C" {
#endif

void Bp(unsigned char n);
int bitscopy(unsigned char* dest,int dbo/*dest bits offset*/,const unsigned char* src,int sbo/*src bits offset*/,int nbits);
void bits_to_bytes_array(uint8_t *bits_array, uint8_t *bytes_array, uint32_t array_len);

uint32_t getbitu(const uint8_t *buff, int pos, int len);
int32_t getbits(const uint8_t *buff, int pos, int len);
/*--------------------------------------------------------------------*/
#ifdef __cplusplus
}
#endif
#endif
