#include "crc.h"

#define CRC8_POLY 0x07

unsigned char cal_table_high_first(unsigned char value)
{
	unsigned char i, crc;

	crc = value;
	/* 数据往左移了8位，需要计算8次 */
	for (i = 8; i > 0; --i)
	{
		if (crc & 0x80)  /* 判断最高位是否为1 */
		{
			crc = (crc << 1) ^ CRC8_POLY;
		}
		else
		{
			crc = (crc << 1);
		}
	}
	return crc;
}

void  create_crc_table()
{
	unsigned short i;
	unsigned char j;

	for (i = 0; i <= 0xFF; i++)
	{
		if (0 == (i % 8))
			printf("\n");

		j = i & 0xFF;
		printf("0x%.2x, ", cal_table_high_first(j));  /*依次计算每个字节的crc校验值*/
	}
}

uint8_t crc4_itu(uint8_t *data, uint32_t len)
{
	uint8_t i;
	uint8_t crc = 0;                // Initial value
	while (len--)
	{
		crc ^= *data++;                 // crc ^= *data; data++;
		for (i = 0; i < 8; ++i)
		{
			if (crc & 1)
				crc = (crc >> 1) ^ 0x09;// 0x0C = reverse 0x03 | 0x09 = reverse 0x09
			else
				crc = (crc >> 1);
		}
	}
	return crc;
}

uint8_t crc8_ccitt(uint8_t *data, uint32_t len)
{
	uint8_t crc_value = 0;
	while (len--)
	{
		crc_value = CRC8_CCITT_TABLE[crc_value ^ *data++];
	}
	return (crc_value);
}

uint16_t crc16_ccitt(uint8_t *data, uint32_t len)
{
	uint16_t crc_value = 0;
	while (len--)
	{
		crc_value = ((crc_value << 8) ^ CRC16_CCITT_TABLE[((crc_value >> 8) ^ *data++) & 0xFFU]) & 0xFFFFU;
	}
	return crc_value;
}

uint32_t crc24_radix(uint8_t *data, uint32_t len)
{
	uint32_t crc_value = 0;
	int i;
	while (len--)
	{
		crc_value = ((crc_value << 8) & 0xFFFFFFU) ^ CRC24_TABLE[(crc_value >> 16) ^ *data++];
	}
	return crc_value;
}

uint32_t crc32_ccitt(uint8_t *data, uint32_t len)
{
#define INIT  0xFFFFFFFFU
#define XOROT 0xFFFFFFFFU

	uint32_t crc_value = INIT;

	while (len--)
	{
		crc_value = CRC32_CCITT_TABLE[(crc_value ^ *data++) & 0xFFU] ^ (crc_value >> 8);
	}
	/* XOR the output value */
	return crc_value ^ XOROT;
}