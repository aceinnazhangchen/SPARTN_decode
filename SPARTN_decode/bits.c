
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "bits.h"

//print bit
void Bp(unsigned char n)  
{  
	int i;
	for (i=7;i>=0;i--)  
	{  
		printf("%u",(n>>i)&1);  
	}
}

//按比特位拷贝
//从src数组首地址跳过sbb个字节，又跳过ssb个比特位，拷贝nbits个比特位的数据到
//dest数组首地址跳过dbb个字节，又跳过dsb个比特位位置
int copybits(const unsigned char* src,int sbb/*source begin byte*/,int ssb/*source skip bit*/,
			unsigned char* dest,int dbb/*dest begin byte*/,int dsb/*dest skip bit*/,int nbits)
{
	// assert(src && dest && sbb>=0 && ssb>=0 && dbb>=0 && dsb>=0 && nbits>=0);
	if(src ==NULL || dest == NULL)return -1;
	if(sbb < 0 || ssb < 0 || dbb < 0 || dsb <0)return -2;
	if(nbits==0)return 0;
	
	if(ssb == dsb){
		//边界对其情况
		//1拷贝对齐部分
		int copybyte=(nbits -(8-ssb))/8;
		memmove(dest+dbb+1,src+sbb+1,copybyte);
		//2拷贝前端不完整字节
		if(ssb != 0){
			unsigned char s=src[sbb];
			s &= 0xFF>>ssb;
			dest[dbb] &= 0xFF<<(8-ssb);
			dest[dbb] |= s;
		}
		//拷贝后端不完整字节
		int endbit=(nbits - (8- ssb))%8;
		if(endbit != 0){
			unsigned char s=src[sbb+1+copybyte];
			s &= 0xFF<<(8-endbit);
			dest[dbb+1 + copybyte] &= 0xFF>>endbit;
			dest[dbb+1 + copybyte] |= s;	
		}
		return (8 - endbit);
	}
	//-------------------------------------------------
	
	int sbgb = sbb*8 + ssb;	//源开始的比特位置
	int dbgb = dbb*8 + dsb;	//目标开始比特位置
	int i,ret;
	int k1,m1,k2,m2;
	unsigned char s;
	if(((dest - src)*8+dbgb) < sbgb  ){
		// 目标开始位置在源开始位置左边
		for(i=0;i<nbits;++i){
			//拷贝某个位
			//1、源位置			目标位置  	 >>3 == /8   &0x7 == %8
			k1=(sbgb+i)>>3;		k2=(dbgb+i)>>3;			//得到字节位
			m1=(sbgb+i)&0x7;	m2=(dbgb+i)&0x7;		//得到比特位 0~7
			s=src[k1];
			s &= 0x80>>m1;	//获取源比特位
			if(m1!=m2){	//偏移位
				s = m1<m2? (s>>(m2-m1)):(s<<(m1-m2));
			}
			dest[k2] &= (~(0x80>>m2));	//目标位置0
			dest[k2] |= s;	//目标位赋值
		}
	}
	else{
		for(i=nbits-1; i >=0 ;--i){
			//拷贝某个位
			//1、源位置			目标位置
			k1=(sbgb+i)>>3;		k2=(dbgb+i)>>3;
			m1=(sbgb+i)&0x7;	m2=(dbgb+i)&0x7;
			s=src[k1];
			s &= 0x80>>m1;	//获取源比特位
			if(m1!=m2){	//偏移位
				s = m1<m2? (s>>(m2-m1)):(s<<(m1-m2));
			}
			dest[k2] &= (~(0x80>>m2));	//目标位置0
			dest[k2] |= s;	//目标位赋值
		}

	}
	return (8 - (dbgb+nbits)%8);
}


 int bitscopy(unsigned char* dest,int dbo/*dest bits offset*/,const unsigned char* src,int sbo/*src bits offset*/,int nbits)
{
	if(src ==NULL || dest == NULL)return -1;
	if(dbo < 0 || dbo < 0)return -2;
	if(nbits<=0)return 0;

	int i;
	int k1,m1,k2,m2;
	unsigned char s;
	//1、源位置		目标位置  	 	
	k1=sbo>>3;		k2=dbo>>3;		//得到字节位  >>3 == /8
	m1=sbo&0x7;		m2=dbo&0x7;		//得到比特位  &0x7 == %8   0~7
	if(m1 == m2)
	{
		int copybyte=(nbits -(8-m1))/8;
		memcpy(dest+k2+1,src+k1+1,copybyte);
		//2拷贝前端不完整字节
		if(m1 != 0){
			unsigned char s=src[k1];
			s &= 0xFF>>m1;
			dest[k2] &= 0xFF<<(8-m1);
			dest[k2] |= s;
		}
		else {
			dest[k2] = src[k1];
		}
		//拷贝后端不完整字节
		int endbit=(nbits - (8- m1))%8;
		if(endbit != 0){
			unsigned char s=src[k1+1+copybyte];
			s &= 0xFF<<(8-endbit);
			dest[k2+1 + copybyte] &= 0xFF>>endbit;
			dest[k2+1 + copybyte] |= s;	
		}
	}
	else
	{
		for(i=0;i<nbits;++i){
			//拷贝某个位
			//1、源位置			目标位置  	 			 
			k1=(sbo+i)>>3;		k2=(dbo+i)>>3;		//得到字节位  >>3 == /8
			m1=(sbo+i)&0x7;		m2=(dbo+i)&0x7;		//得到比特位  &0x7 == %8   0~7
			s=src[k1];
			s &= 0x80>>m1;	//获取源比特位
			if(m1!=m2){	//偏移位
				s = m1<m2? (s>>(m2-m1)):(s<<(m1-m2));
			}
			dest[k2] &= (~(0x80>>m2));	//目标位置0
			dest[k2] |= s;	//目标位赋值
		}
	}
	
	return dbo+nbits;
}

 void bits_to_bytes_array(uint8_t *bits_array, uint8_t *bytes_array, uint32_t array_len) {
	 int i, bi;
	 for (i = 0; i < array_len; i += 8) {
		 uint8_t c = bits_array[i / 8];
		 for (bi = 7; bi >= 0; bi--) {
			 bytes_array[i + 7 - bi] = ((c >> bi) & 1);
		 }
	 }
	 for (i = 0; i < array_len; i++) {
		 printf("%d", bytes_array[i]);
	 }
 }

// int main()
// {
// 	int i;
// 	unsigned char src[10]={0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff};//{0,0x1F,0x0F,0x0F,0x70,0,0,1,1};
// 	unsigned char dst[10]={0,0,0,0,0,0,0,0,0};

// 	//printf("%d\n",copybits(src,1,3,dst,0,5,24));
// 	printf("%d\n",bitscopy(dst,18,src,18,32));
// 	printf("%d\n",bitscopy(dst,50,src,18,8));

// 	for(i=0;i<10;++i){
// 		//printf("\t%2x\t%2x\n",src[i],dst[i]);
// 		Bp(src[i]);
// 		putchar(' ');
// 	}
// 	putchar('\n');
// 	for(i=0;i<10;++i){
// 		//printf("\t%2x\t%2x\n",src[i],dst[i]);
// 		Bp(dst[i]);
// 		putchar(' ');
// 	}
// 	putchar('\n');
// 	return 0;
// }