#ifndef SPARTN_H
#define SPARTN_H

#include <stdint.h>

typedef struct {
	uint32_t nbyte;
	uint8_t buff[1200];
	uint32_t type;
	uint32_t len;
	uint32_t EAF;
	uint32_t CRC_type;
	uint32_t Frame_CRC;
	uint32_t Subtype;
	uint32_t Time_tag_type;
	uint32_t GNSS_time_type;
	uint32_t Solution_ID;
	uint32_t Solution_processor_ID;
	uint32_t Encryption_ID;
	uint32_t ESN;						//Encryption Sequence Number
	uint32_t AI;						//Authentication Indicator 
	uint32_t EAL;						//Embedded Authentication Length
	uint32_t Payload_offset;
	uint32_t Message_CRC;				//
} spartn_t;

#endif // !CRC_H
