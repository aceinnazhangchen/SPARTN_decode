#include "spartn.h"
#include "log.h"

void decode_GPS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint8_t *satellite_mask_len) {
	int i,offset = *pos;
	int tab = 2;
	uint8_t SF011_Type = getbitu(data, offset, 2);  offset += 2; log(LOG_DEBUG, tab, "SF011_Type = %d", SF011_Type);
	uint8_t SF011_Len = 0;
	switch (SF011_Type) {
	case 0:SF011_Len = 32; break;
	case 1:SF011_Len = 44; break;
	case 2:SF011_Len = 56; break;
	case 3:SF011_Len = 64; break;
	}
	*satellite_mask_len = SF011_Len;
	log(LOG_DEBUG, tab, "SF011_Len = %d", SF011_Len);
	uint8_t SF011[8] = { 0 };
	bitscopy(SF011, 0, data, offset, SF011_Len); offset += SF011_Len;
	bits_to_bytes_array(SF011, satellite_mask, *satellite_mask_len);
	char str_satellite_mask[65] = { 0 };
	for (i = 0; i < *satellite_mask_len; i++) {
		str_satellite_mask[i] = satellite_mask[i] ? '1' : '0';
	}
	log(LOG_DEBUG, tab, "SF011 = %s", str_satellite_mask);
	*pos = offset;
}

void decode_GLONASS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint8_t *satellite_mask_len) {
	int i,offset = *pos;
	int tab = 2;
	uint8_t SF012_Type = getbitu(data, offset, 2);  offset += 2; log(LOG_DEBUG, tab, "SF012_Type = %d", SF012_Type);
	uint8_t SF012_Len = 0;
	switch (SF012_Type) {
	case 0:SF012_Len = 24; break;
	case 1:SF012_Len = 36; break;
	case 2:SF012_Len = 48; break;
	case 3:SF012_Len = 63; break;
	}
	*satellite_mask_len = SF012_Len;
	log(LOG_DEBUG, tab, "SF012_Len = %d", SF012_Len);
	uint8_t SF012[8] = { 0 };
	bitscopy(SF012, 0, data, offset, SF012_Len); offset += SF012_Len;
	bits_to_bytes_array(SF012, satellite_mask, *satellite_mask_len);
	char str_satellite_mask[64] = { 0 };
	for (i = 0; i < *satellite_mask_len; i++) {
		str_satellite_mask[i] = satellite_mask[i] ? '1' : '0';
	}
	log(LOG_DEBUG, tab, "SF012 = %s", str_satellite_mask);
	*pos = offset;
}