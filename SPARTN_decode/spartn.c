#include <stdio.h>
#include <memory.h>
#include <stdarg.h>
#include "crc.h"
#include "bits.h"
#include "spartn.h"

#define SPARTN_PREAMB 0x73 

#define LOG_DEBUG 1
void log(int level,int tab,const char *format, ...) {
	va_list ap;
	int i;
	char buffer[2048] = { 0 };
	char tab_buffer[128] = { 0 };
	if (level >= LOG_DEBUG) {
		va_start(ap, format); vsprintf(buffer, format, ap); va_end(ap);
		for (i = 0; i < tab; i++) {
			tab_buffer[i] = ' ';
		}
		printf("%s%s \n", tab_buffer,buffer);
	}
}

uint32_t getbitu(const uint8_t *buff, int pos, int len){
	uint32_t bits = 0;
	int i;
	for (i = pos; i < pos + len; i++) bits = (bits << 1) + ((buff[i / 8] >> (7 - i % 8)) & 1u);
	return bits;
}

int32_t getbits(const uint8_t *buff, int pos, int len){
	uint32_t bits = getbitu(buff, pos, len);
	if (len <= 0 || 32 <= len || !(bits&(1u << (len - 1)))) return (int32_t)bits;
	return (int32_t)(bits | (~0u << len)); /* extend sign */
}

void decode_GPS_satellite_mask(uint8_t* data, int *pos,uint8_t* satellite_mask, uint32_t *satellite_mask_len) {
	int offset = *pos;
	int tab = 2;
	uint32_t SF011_Type = getbitu(data, offset, 2);  offset += 2; log(LOG_DEBUG, tab, "SF011_Type = %d", SF011_Type);
	uint32_t SF011_Len = 0;
	switch (SF011_Type) {
	case 0:SF011_Len = 32; break;
	case 1:SF011_Len = 44; break;
	case 2:SF011_Len = 56; break;
	case 3:SF011_Len = 64; break;
	}
	*satellite_mask_len = SF011_Len;
	log(LOG_DEBUG, tab, "SF011_Len = %d", SF011_Len);
	printf("  SF011 = ");
	uint8_t SF011[8] = { 0 };
	bitscopy(SF011, 0, data, offset, SF011_Len); offset += SF011_Len;
	bits_to_bytes_array(SF011, satellite_mask, *satellite_mask_len);
	log(LOG_DEBUG, tab, "");
	*pos = offset;
}

void decode_GLONASS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint32_t *satellite_mask_len) {
	int offset = *pos;
	int tab = 2;
	uint32_t SF012_Type = getbitu(data, offset, 2);  offset += 2; log(LOG_DEBUG, tab, "SF012_Type = %d", SF012_Type);
	uint32_t SF012_Len = 0;
	switch (SF012_Type) {
	case 0:SF012_Len = 24; break;
	case 1:SF012_Len = 36; break;
	case 2:SF012_Len = 48; break;
	case 3:SF012_Len = 63; break;
	}
	*satellite_mask_len = SF012_Len;
	log(LOG_DEBUG, tab, "SF012_Len = %d", SF012_Len);
	printf("  SF012 = ");
	uint8_t SF012[8] = { 0 };
	bitscopy(SF012, 0, data, offset, SF012_Len); offset += SF012_Len;
	bits_to_bytes_array(SF012, satellite_mask, *satellite_mask_len);
	log(LOG_DEBUG, tab, "");
	*pos = offset;
}

void decode_bias_mask(uint8_t* data,int *pos, uint32_t *bias_array,uint32_t effective_len,uint32_t subType) {
	int i, offset = *pos;
	int tab = 4;
	uint32_t len_flag = getbitu(data, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "len_flag = %d", len_flag);
	uint32_t max_len = 0;
	if (len_flag == 0) {
		max_len = subType ? 5:6;
	}
	else if (len_flag == 1) {
		max_len = subType ? 9:11;
	}

	for (i = 0; i < max_len; i++) {
		bias_array[i] = getbitu(data, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "bias_array[%d] = %d", i, bias_array[i]);
		if (i == effective_len - 1) { offset += (max_len - effective_len); break; }//just have 3/2 type now
	}
	*pos = offset;
}

void decode_phase_bias_block(uint8_t* data, int *pos, uint32_t *bias_array, uint32_t effective_len) {
	int i, offset = *pos;
	int tab = 4;
	uint32_t SF023_Fix_flag;
	uint32_t SF015_Continuity_indicator;
	uint32_t SF020_Phase_bias_correction;
	for (i = 0; i < effective_len; i++) {
		if (bias_array[i] == 1) {
			SF023_Fix_flag = getbitu(data, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF023_Fix_flag = %d", SF023_Fix_flag);
			SF015_Continuity_indicator = getbitu(data, offset, 3);  offset += 3; log(LOG_DEBUG, tab, "SF015_Continuity_indicator = %d", SF015_Continuity_indicator);
			SF020_Phase_bias_correction = getbitu(data, offset, 14);  offset += 14; log(LOG_DEBUG, tab, "SF020_Phase_bias_correction = %f", SF020_Phase_bias_correction*0.002-16.382);
		}
	}
	*pos = offset;
}

void decode_code_bias_correction(uint8_t* data, int *pos, uint32_t *bias_array, uint32_t effective_len) {
	int i, offset = *pos;
	int tab = 4;
	for (i = 0; i < effective_len; i++) {
		if (bias_array[i] == 1) {
			uint32_t SF029_Code_bias_correction = getbitu(data, offset, 11);  offset += 11; log(LOG_DEBUG, tab, "SF029_Code_bias_correction = %f", SF029_Code_bias_correction*0.02-20.46);
		}
	}
	*pos = offset;
}

void decode_satellite_block(spartn_t* spartn, int *pos,uint32_t SF008_Yaw_present_flag) {
	int i, offset = *pos; 
	int tab = 3;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	//Table 6.4 satellite block
	uint32_t SF013_DNU = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF013_DNU = %d", SF013_DNU);
	uint32_t SF014_Orbit_block_0 = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF014_Orbit_block_0 = %d", SF014_Orbit_block_0);
	uint32_t SF014_Clock_block_1 = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF014_Clock_block_1 = %d", SF014_Clock_block_1);
	uint32_t SF014_Bias_block_2 = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF014_Bias_block_2 = %d", SF014_Bias_block_2);
	uint32_t SF015 = getbitu(payload, offset, 3);  offset += 3; log(LOG_DEBUG, tab, "SF015 = %d", SF015);
	//Table 6.5 orbit block 
	if (SF014_Orbit_block_0) {
		if (spartn->Subtype == 0) {
			uint32_t SF018_IODE = getbitu(payload, offset, 8);  offset += 8; log(LOG_DEBUG, tab, "SF018_IODE = %d", SF018_IODE);
		}
		else if (spartn->Subtype == 1) {
			uint32_t SF019_IODE = getbitu(payload, offset, 7);  offset += 7; log(LOG_DEBUG, tab, "SF019_IODE = %d", SF019_IODE);
		}
		uint32_t SF020_radial = getbitu(payload, offset, 14);  offset += 14; log(LOG_DEBUG, tab, "SF020_radial = %f", SF020_radial*0.002 - 16.382);
		uint32_t SF020_along = getbitu(payload, offset, 14);  offset += 14; log(LOG_DEBUG, tab, "SF020_along = %f", SF020_along*0.002 - 16.382);
		uint32_t SF020_cross = getbitu(payload, offset, 14);  offset += 14; log(LOG_DEBUG, tab, "SF020_cross = %f", SF020_cross*0.002 - 16.382);
		if (SF008_Yaw_present_flag == 1) {
			uint32_t SF021 = getbitu(payload, offset, 6);  offset += 6; log(LOG_DEBUG, tab, "SF021 = %d", SF021*6);
		}
	}
	//Table 6.6 clock block 
	if (SF014_Clock_block_1) {
		uint32_t SF022_IODE_continuity = getbitu(payload, offset, 3);  offset += 3; log(LOG_DEBUG, tab, "SF022_IODE_continuity = %d", SF022_IODE_continuity);
		uint32_t SF020_Clock_correction = getbitu(payload, offset, 14);  offset += 14; log(LOG_DEBUG, tab, "SF020_Clock_correction = %f", SF020_Clock_correction*0.002-16.382);
		uint32_t SF024_User_range_error = getbitu(payload, offset, 3);  offset += 3; log(LOG_DEBUG, tab, "SF024_User_range_error = %d", SF024_User_range_error);
	}
	//bias block
	if (SF014_Bias_block_2) {
		if (spartn->Subtype == 0) {
			//Table 6.7 GPS bias block
			uint32_t SF025_phase_bias[11] = { 0 };
			uint32_t SF025_phase_bias_effective_len = 3;
			decode_bias_mask(payload, &offset, SF025_phase_bias, SF025_phase_bias_effective_len, spartn->Subtype);
			//Table 6.9 Phase bias block (Repeated)
			decode_phase_bias_block(payload, &offset, SF025_phase_bias, SF025_phase_bias_effective_len);
			//SF027
			uint32_t SF027_phase_bias[11] = { 0 };
			uint32_t SF027_phase_bias_effective_len = 3;
			decode_bias_mask(payload, &offset, SF027_phase_bias, SF027_phase_bias_effective_len, spartn->Subtype);
			decode_code_bias_correction(payload, &offset, SF027_phase_bias, SF027_phase_bias_effective_len);
		}
		else if (spartn->Subtype == 1) {
			//Table 6.8 GLONASS bias block
			uint32_t SF026_phase_bias[9] = { 0 };
			uint32_t SF026_phase_bias_effective_len = 2;
			decode_bias_mask(payload, &offset, SF026_phase_bias, SF026_phase_bias_effective_len, spartn->Subtype);
			//Table 6.9 Phase bias block (Repeated)
			decode_phase_bias_block(payload, &offset, SF026_phase_bias, SF026_phase_bias_effective_len);
			//SF028
			uint32_t SF028_phase_bias[9] = { 0 };
			uint32_t SF028_phase_bias_effective_len = 2;
			decode_bias_mask(payload, &offset, SF028_phase_bias, SF028_phase_bias_effective_len, spartn->Subtype);
			decode_code_bias_correction(payload, &offset, SF028_phase_bias, SF028_phase_bias_effective_len);
		}
	}
	*pos = offset;
}

int decode_OCB_message(spartn_t* spartn) {
	int i, bi, tab = 2;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	int offset = 0;
	//Table 6.3 Header block 
	uint32_t SF005_SIOU = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF005_SIOU = %d", SF005_SIOU);
	uint32_t SF010_EOS = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF010_EOS = %d", SF010_EOS);
	uint32_t SF069_Reserved = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF069_Reserved = %d", SF069_Reserved);
	uint32_t SF008_Yaw_present_flag = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF008_Yaw_present_flag = %d", SF008_Yaw_present_flag);
	uint32_t SF009_Satellite_reference_datum = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF009_Satellite_reference_datum = %d", SF009_Satellite_reference_datum);
	uint32_t SF016_017_Ephemeris_type = getbitu(payload, offset, 2);  offset += 2; log(LOG_DEBUG, tab, "SF016_017_Ephemeris_type = %d", SF016_017_Ephemeris_type);
	uint8_t satellite_mask[64] = { 0 };
	uint32_t satellite_mask_len = 0;
	if (spartn->Subtype == 0) {
		decode_GPS_satellite_mask(payload, &offset, satellite_mask,&satellite_mask_len);
	}
	else if (spartn->Subtype == 1) {
		decode_GLONASS_satellite_mask(payload, &offset, satellite_mask, &satellite_mask_len);
	}
	//Table 6.4 Satellite block (Repeated) 
	for (i = 0; i < satellite_mask_len; i++) {
		if (satellite_mask[i]) {
			log(LOG_DEBUG, tab, "PRN_ID = %d", i+1);
			decode_satellite_block(spartn, &offset, SF008_Yaw_present_flag);
		}
	}	
	log(LOG_DEBUG, tab, "offset = %d bits", offset);
	return 1;
}

void decode_ionosphere_satellite_block(spartn_t* spartn, int *pos, uint32_t SF040_Iono, uint32_t SF054_Ionosphere_equation_type, uint32_t SF039_Number_grid_points_present) {
	int i, offset = *pos;
	int tab = 4;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	//Table 6.17 Ionosphere satellite block
	if (SF040_Iono == 1 || SF040_Iono == 2) {
		uint32_t SF055_Ionosphere_quality = getbitu(payload, offset, 4); offset += 4; log(LOG_DEBUG, tab, "SF055_Ionosphere_quality = %d", SF055_Ionosphere_quality);
		uint32_t SF056_Ionosphere_satellite_polynomial_block = getbitu(payload, offset, 1); offset += 1; log(LOG_DEBUG, tab, "SF056_Ionosphere_satellite_polynomial_block = %d", SF056_Ionosphere_satellite_polynomial_block);
		if (SF056_Ionosphere_satellite_polynomial_block) {
			//Table 6.19 
			uint32_t SF060_C00 = getbitu(payload, offset, 14); offset += 14; log(LOG_DEBUG, tab, "SF060_C00 = %f", SF060_C00*0.04 - 327.64);//0,1,2
			if (SF054_Ionosphere_equation_type == 1 || SF054_Ionosphere_equation_type == 2) {
				uint32_t SF061_C01 = getbitu(payload, offset, 14); offset += 14; log(LOG_DEBUG, tab, "SF061_C01 = %f", SF061_C01*0.008 - 65.528);
				uint32_t SF061_C10 = getbitu(payload, offset, 14); offset += 14; log(LOG_DEBUG, tab, "SF061_C10 = %f", SF061_C10*0.008 - 65.528);
			}
			if (SF054_Ionosphere_equation_type == 2) {
				uint32_t SF062_C11 = getbitu(payload, offset, 15); offset += 15; log(LOG_DEBUG, tab, "SF062_C11 = %f", SF062_C11*0.002 - 32.766);
			}
		}
		else {
			//Table 6.18 
			uint32_t SF057_C00 = getbitu(payload, offset, 12); offset += 12; log(LOG_DEBUG, tab, "SF057_C00 = %f", SF057_C00 * 0.04 - 81.88);//0,1,2
			if (SF054_Ionosphere_equation_type == 1 || SF054_Ionosphere_equation_type == 2) {
				uint32_t SF058_C01 = getbitu(payload, offset, 12); offset += 12; log(LOG_DEBUG, tab, "SF058_C01 = %f", SF058_C01 * 0.008 - 16.376);
				uint32_t SF058_C10 = getbitu(payload, offset, 12); offset += 12; log(LOG_DEBUG, tab, "SF058_C10 = %f", SF058_C10 * 0.008 - 16.376);
			}
			if (SF054_Ionosphere_equation_type == 2) {
				uint32_t SF059_C11 = getbitu(payload, offset, 13); offset += 13; log(LOG_DEBUG, tab, "SF059_C11 = %f", SF059_C11 * 0.002 - 8.190);
			}
		}
	}
	if (SF040_Iono == 2) {
		uint32_t SF063_Ionosphere_residual_field_size = getbitu(payload, offset, 2); offset += 2; log(LOG_DEBUG, tab, "SF063_Ionosphere_residual_field_size = %d", SF063_Ionosphere_residual_field_size);
		switch (SF063_Ionosphere_residual_field_size) {
		case 0:
		{
			uint32_t SF064 = 0;
			for (i = 0; i < SF039_Number_grid_points_present; i++) {
				SF064 = getbitu(payload, offset, 4); offset += 4; log(LOG_DEBUG, tab, "SF064 = %f", SF064*0.04 - 0.28);
			}
		}
		break;
		case 1:
		{
			uint32_t SF065 = 0;
			for (i = 0; i < SF039_Number_grid_points_present; i++) {
				SF065 = getbitu(payload, offset, 7); offset += 7; log(LOG_DEBUG, tab, "SF065 = %f", SF065*0.04 - 2.52);
			}
		}
		break;
		case 2:
		{
			uint32_t SF066 = 0;
			for (i = 0; i < SF039_Number_grid_points_present; i++) {
				SF066 = getbitu(payload, offset, 10); offset += 10; log(LOG_DEBUG, tab, "SF066 = %f", SF066*0.04 - 20.44);
			}
		}
		break;
		case 3:
		{
			uint32_t SF067 = 0;
			for (i = 0; i < SF039_Number_grid_points_present; i++) {
				SF067 = getbitu(payload, offset, 14); offset += 14; log(LOG_DEBUG, tab, "SF067 = %f", SF067*0.04 - 327.64);
			}
		}
		break;
		}
	}

	*pos = offset;
}
void decode_Atmosphere_block(spartn_t* spartn , int *pos) {
	int i, offset = *pos;
	int tab = 3;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	//Table 6.12 Area data block 
	uint32_t SF031_Area_ID = getbitu(payload, offset, 8);  offset += 8; log(LOG_DEBUG, tab, "SF031_Area_ID = %d", SF031_Area_ID);
	uint32_t SF039_Number_grid_points_present = getbitu(payload, offset, 7);  offset += 7; log(LOG_DEBUG, tab, "SF039_Number_grid_points_present = %d", SF039_Number_grid_points_present);
	uint32_t SF040_Tropo = getbitu(payload, offset, 2);  offset += 2; log(LOG_DEBUG, tab, "SF040_Tropo = %d", SF040_Tropo);
	uint32_t SF040_Iono = getbitu(payload, offset, 2);  offset += 2; log(LOG_DEBUG, tab, "SF040_Iono = %d", SF040_Iono);
	//Table 6.13 Troposphere data block
	//Troposphere polynomial coefficient block 
	if (SF040_Tropo == 1 || SF040_Tropo == 2) {
		uint32_t SF041_Troposphere_equation_type = getbitu(payload, offset, 3);  offset += 3; log(LOG_DEBUG, tab, "SF041_Troposphere_equation_type = %d", SF041_Troposphere_equation_type);
		uint32_t SF042_Troposphere_quality = getbitu(payload, offset, 3);  offset += 3; log(LOG_DEBUG, tab, "SF042_Troposphere_quality = %d", SF042_Troposphere_quality);
		uint32_t SF043_Area_average_vertical_hydrostatic_delay = getbitu(payload, offset, 8);  offset += 8; log(LOG_DEBUG, tab, "SF043_Area_average_vertical_hydrostatic_delay = %f", SF043_Area_average_vertical_hydrostatic_delay*0.004 - 0.508);
		uint32_t SF044_Troposphere_polynomial_coefficient_size_indicator = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF044_Troposphere_polynomial_coefficient_size_indicator = %d", SF044_Troposphere_polynomial_coefficient_size_indicator);
		if (SF044_Troposphere_polynomial_coefficient_size_indicator) {
			//Table 6.15
			if (SF041_Troposphere_equation_type == 0 || SF041_Troposphere_equation_type == 1 || SF041_Troposphere_equation_type == 2) {
				uint32_t SF048_T00 = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF048_T00 = %f", SF048_T00 * 0.004 - 1.020);//0,1,2
			}
			if (SF041_Troposphere_equation_type == 1 || SF041_Troposphere_equation_type == 2) {
				uint32_t SF049_T01 = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF049_T01 = %f", SF049_T01 * 0.001 - 0.255);
				uint32_t SF049_T10 = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF049_T10 = %f", SF049_T10 * 0.001 - 0.255);
			}
			if (SF041_Troposphere_equation_type == 2) {
				uint32_t SF050_T11 = getbitu(payload, offset, 11); offset += 11; log(LOG_DEBUG, tab, "SF050_T11 = %f", SF050_T11 * 0.0002 - 0.2046);
			}
		}
		else {
			//Table 6.14
			if (SF041_Troposphere_equation_type == 0 || SF041_Troposphere_equation_type == 1 || SF041_Troposphere_equation_type == 2) {
				uint32_t SF045_T00 = getbitu(payload, offset, 7); offset += 7; log(LOG_DEBUG, tab, "SF045_T00 = %f", SF045_T00 * 0.004 - 0.252);//0,1,2
			}
			if (SF041_Troposphere_equation_type == 1 || SF041_Troposphere_equation_type == 2) {
				uint32_t SF046_T01 = getbitu(payload, offset, 7); offset += 7; log(LOG_DEBUG, tab, "SF046_T01 = %f", SF046_T01 * 0.001 - 0.063);
				uint32_t SF046_T10 = getbitu(payload, offset, 7); offset += 7; log(LOG_DEBUG, tab, "SF046_T10 = %f", SF046_T10 * 0.001 - 0.063);
			}
			if (SF041_Troposphere_equation_type == 2) {
				uint32_t SF047_T11 = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF047_T11 = %f", SF047_T11 * 0.0002 - 0.0510);
			}
		}
	}
	//Troposphere grid block 
	if (SF040_Tropo == 2) {
		uint32_t SF051_Troposphere_residual_field_size = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF051_Troposphere_residual_field_size = %d", SF051_Troposphere_residual_field_size);
		if (SF051_Troposphere_residual_field_size) {
			//SF053
			uint32_t SF053 = 0;
			for (i = 0; i < SF039_Number_grid_points_present; i++) {
				SF053 = getbitu(payload, offset, 8); offset += 8; log(LOG_DEBUG, tab, "SF053 = %f", SF053 * 0.004 - 0.508);
			}
		}
		else {
			//SF052
			uint32_t SF052 = 0;
			for (i = 0; i < SF039_Number_grid_points_present; i++) {
				SF052 = getbitu(payload, offset, 6); offset += 6; log(LOG_DEBUG, tab, "SF052 = %f", SF052 * 0.004 - 0.124);
			}
		}
	}
	//Table 6.16 Ionosphere block 
	if (SF040_Iono == 1 || SF040_Iono == 2) {
		uint32_t SF054_Ionosphere_equation_type = getbitu(payload, offset, 3); offset += 3; log(LOG_DEBUG, tab, "SF054_Ionosphere_equation_type = %d", SF054_Ionosphere_equation_type);
		uint8_t satellite_mask[64] = { 0 };
		uint32_t satellite_mask_len = 0;
		if (spartn->Subtype == 0) {
			decode_GPS_satellite_mask(payload, &offset, satellite_mask, &satellite_mask_len);
		}
		else  if (spartn->Subtype == 1) {
			decode_GLONASS_satellite_mask(payload, &offset, satellite_mask, &satellite_mask_len);
		}
		//Table 6.17 Ionosphere satellite block (Repeated)
		for (i = 0; i < satellite_mask_len; i++) {
			if (satellite_mask[i]) {
				log(LOG_DEBUG, tab, "PRN_ID = %d", i + 1);
				decode_ionosphere_satellite_block(spartn, &offset, SF040_Iono, SF054_Ionosphere_equation_type, SF039_Number_grid_points_present);
			}
		}
	}
	*pos = offset;
}
int decode_HPAC_message(spartn_t* spartn) {
	int i, bi, tab = 2;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	int offset = 0;
	//Table 6.10 Header block
	uint32_t SF005_SIOU = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF005_SIOU = %d", SF005_SIOU);
	uint32_t SF069_Reserved = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF069_Reserved = %d", SF069_Reserved);
	uint32_t SF068_AIOU = getbitu(payload, offset, 4);  offset += 4; log(LOG_DEBUG, tab, "SF068_AIOU = %d", SF068_AIOU);
	uint32_t SF030_Area_count = getbitu(payload, offset, 5)+1;  offset += 5; log(LOG_DEBUG, tab, "SF030_Area_count = %d", SF030_Area_count);
	//Table 6.11 Atmosphere block (Repeated)
	for (i = 0; i < SF030_Area_count; i++) {
		decode_Atmosphere_block(spartn, &offset);
		//break;
	}	
	log(LOG_DEBUG, tab, "offset = %d bits", offset);
	return 1;
}

void decode_Area_definition_block(spartn_t* spartn, int *pos) {
	int i, offset = *pos;
	int tab = 3;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	//Table 6.21  Overall layout of SPARTN message 
	uint32_t SF031_Area_ID = getbitu(payload, offset, 8);  offset += 8; log(LOG_DEBUG, tab, "SF031_Area_ID = %d", SF031_Area_ID);
	uint32_t SF032_Area_reference_latitude = getbitu(payload, offset, 11);  offset += 11; log(LOG_DEBUG, tab, "SF032_Area_reference_latitude = %f", SF032_Area_reference_latitude*0.1- 90);
	uint32_t SF033_Area_reference_longitude = getbitu(payload, offset, 12);  offset += 12; log(LOG_DEBUG, tab, "SF033_Area_reference_longitude = %f", SF033_Area_reference_longitude*0.1 - 180);
	uint32_t SF034_Area_latitude_grid_node_count = getbitu(payload, offset, 3);  offset += 3; log(LOG_DEBUG, tab, "SF034_Area_latitude_grid_node_count = %d", SF034_Area_latitude_grid_node_count);
	uint32_t SF035_Area_longitude_grid_node_count = getbitu(payload, offset, 3);  offset += 3; log(LOG_DEBUG, tab, "SF035_Area_longitude_grid_node_count = %d", SF035_Area_longitude_grid_node_count);
	uint32_t SF036_Area_latitude_grid_node_spacing = getbitu(payload, offset, 5);  offset += 5; log(LOG_DEBUG, tab, "SF036_Area_latitude_grid_node_spacing = %f", SF036_Area_latitude_grid_node_spacing*0.1 + 0.1);
	uint32_t SF037_Area_longitude_grid_node_spacing = getbitu(payload, offset, 5);  offset += 5; log(LOG_DEBUG, tab, "SF037_Area_longitude_grid_node_spacing = %f", SF037_Area_longitude_grid_node_spacing*0.1 + 0.1);
	*pos = offset;
}

int decode_GAD_message(spartn_t* spartn) {
	int i, bi, tab = 2;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	int offset = 0;
	//Table 6.20 Header block
	uint32_t SF005_SIOU = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF005_SIOU = %d", SF005_SIOU);
	uint32_t SF069_Reserved = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF069_Reserved = %d", SF069_Reserved);
	uint32_t SF068_AIOU = getbitu(payload, offset, 4);  offset += 4; log(LOG_DEBUG, tab, "SF068_AIOU = %d", SF068_AIOU);
	uint32_t SF030_Area_count = getbitu(payload, offset, 5)+1;  offset += 5; log(LOG_DEBUG, tab, "SF030_Area_count = %d", SF030_Area_count);
	//Table 6.21 (Repeated)
	for (i = 0; i < SF030_Area_count; i++) {
		decode_Area_definition_block(spartn, &offset);
		//break;
	}
	log(LOG_DEBUG, tab, "offset = %d bits", offset);
	return 1;
}

void decode_LPAC_grid_node_VTEC_block(spartn_t* spartn, int *pos) {
	int i, offset = *pos;
	int tab = 4;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	uint32_t SF055_VTEC_quality = getbitu(payload, offset, 4); offset += 4; log(LOG_DEBUG, tab, "SF055_VTEC_quality = %d", SF055_VTEC_quality);
	uint32_t SSF081_VTEC_size_indicator = getbitu(payload, offset, 1); offset += 1; log(LOG_DEBUG, tab, "SSF081_VTEC_size_indicator = %d", SSF081_VTEC_size_indicator);
	if (SSF081_VTEC_size_indicator) {
		//SF083
		uint32_t SF083_VTEC_residual = getbitu(payload, offset, 11); offset += 11; log(LOG_DEBUG, tab, "SF083_VTEC_residual = %f", SF083_VTEC_residual*0.25 - 255.75);
	}
	else {
		//SF082
		uint32_t SF082_VTEC_residual = getbitu(payload, offset, 7); offset += 7; log(LOG_DEBUG, tab, "SF082_VTEC_residual = %f", SF082_VTEC_residual*0.25- 15.75);
	}
	*pos = offset;
}

void decode_LPAC_area_block(spartn_t* spartn, int *pos) {
	int i, offset = *pos;
	int tab = 3;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	//Table 6.23 LPAC area block
	//Table 6.24 LPAC area data block 
	uint32_t SF072_LPAC_area_ID = getbitu(payload, offset, 2); offset += 2; log(LOG_DEBUG, tab, "SF072_LPAC_area_ID = %d", SF072_LPAC_area_ID);
	uint32_t SF073_LPAC_area_reference_latitude = getbitu(payload, offset, 8); offset += 8; log(LOG_DEBUG, tab, "SF073_LPAC_area_reference_latitude = %d", SF073_LPAC_area_reference_latitude-85);
	uint32_t SF074_LPAC_area_reference_longitude = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF074_LPAC_area_reference_longitude = %d", SF074_LPAC_area_reference_longitude - 180);
	uint32_t SF075_LPAC_area_latitude_grid_node_count = getbitu(payload, offset, 4)+1; offset += 4; log(LOG_DEBUG, tab, "SF075_LPAC_area_latitude_grid_node_count = %d", SF075_LPAC_area_latitude_grid_node_count);
	uint32_t SF076_LPAC_area_longitude_grid_node_count = getbitu(payload, offset, 4)+1; offset += 4; log(LOG_DEBUG, tab, "SF076_LPAC_area_longitude_grid_node_count = %d", SF076_LPAC_area_longitude_grid_node_count);
	uint32_t SF077_LPAC_area_latitude_grid_node_spacing = getbitu(payload, offset, 2); offset += 2; log(LOG_DEBUG, tab, "SF077_LPAC_area_latitude_grid_node_spacing = %d", SF077_LPAC_area_latitude_grid_node_spacing);
	uint32_t SF078_LPAC_area_longitude_grid_node_spacing = getbitu(payload, offset, 2); offset += 2; log(LOG_DEBUG, tab, "SF078_LPAC_area_longitude_grid_node_spacing = %d", SF078_LPAC_area_longitude_grid_node_spacing);
	uint32_t SF080_Average_area_VTEC = getbitu(payload, offset, 12); offset += 12; log(LOG_DEBUG, tab, "SF080_Average_area_VTEC = %f", SF080_Average_area_VTEC*0.25- 511.75);
	uint32_t mask_len = SF075_LPAC_area_latitude_grid_node_count * SF076_LPAC_area_longitude_grid_node_count; log(LOG_DEBUG, tab, "mask_len = %d", mask_len);
	uint8_t SF079_Grid_node_present_mask[32] = { 0 };
	uint8_t Grid_node_present_mask[256] = {0};
	bitscopy(SF079_Grid_node_present_mask, 0, payload, offset, mask_len); offset += mask_len;
	bits_to_bytes_array(SF079_Grid_node_present_mask, Grid_node_present_mask, mask_len);
	log(LOG_DEBUG, tab, "");
	//Table 6.25 LPAC grid node VTEC block  (Repeated)
	for (i = 0; i < mask_len; i++) {
		if (Grid_node_present_mask[i]) {
			decode_LPAC_grid_node_VTEC_block(spartn, &offset);
		}
	}
	
	*pos = offset;
}

int decode_LPAC_message(spartn_t* spartn) {
	int i, bi, tab = 2;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	int offset = 0;
	//Table 6.22 Header block
	uint32_t SF005_SIOU = getbitu(payload, offset, 9); offset += 9; log(LOG_DEBUG, tab, "SF005_SIOU = %d", SF005_SIOU);
	uint32_t SF069_Reserved = getbitu(payload, offset, 1);  offset += 1; log(LOG_DEBUG, tab, "SF069_Reserved = %d", SF069_Reserved);
	uint32_t SF070_Ionosphere_shell_height = getbitu(payload, offset, 2);  offset += 2; log(LOG_DEBUG, tab, "SF070_Ionosphere_shell_height = %d", SF070_Ionosphere_shell_height);
	uint32_t SF071_LPAC_area_count = getbitu(payload, offset, 2)+1;  offset += 2; log(LOG_DEBUG, tab, "SF071_LPAC_area_count = %d", SF071_LPAC_area_count);
	//Table 6.23 LPAC area block (Repeated)
	for (i = 0; i < SF071_LPAC_area_count; i++) {
		decode_LPAC_area_block(spartn, &offset);
	}
	log(LOG_DEBUG, tab, "offset = %d bits", offset);
	return 1;
}

int decode_Dynamic_Key(spartn_t* spartn) {
	return 1;
}

int decode_Group_Authentication(spartn_t* spartn) {
	int i, bi, tab = 2;
	uint8_t* payload = spartn->buff + spartn->Payload_offset;
	int offset = 0;
	//Table 6.27 Message ID block 
	uint32_t SF089_Count_of_message_IDs = getbitu(payload, offset, 5);  offset += 5; log(LOG_DEBUG, tab, "SF089_Count_of_message_IDs = %d", SF089_Count_of_message_IDs);
	uint32_t SF090_Group_authentication_type = getbitu(payload, offset, 4);  offset += 4; log(LOG_DEBUG, tab, "SF090_Group_authentication_type = %d", SF090_Group_authentication_type);
	uint32_t SF091_Computed_authentication_data_length = getbitu(payload, offset, 4);  offset += 4; log(LOG_DEBUG, tab, "SF091_Computed_authentication_data_length = %d", SF091_Computed_authentication_data_length);
	//Table 6.28 
	uint32_t TF002_Message_type = getbitu(payload, offset, 7);  offset += 7; log(LOG_DEBUG, tab, "TF002_Message_type = %d", TF002_Message_type);
	uint32_t TF003_Message_sub_type = getbitu(payload, offset, 4);  offset += 4; log(LOG_DEBUG, tab, "TF003_Message_sub_type = %d", TF003_Message_sub_type);
	uint32_t TF014_Encryption_sequence_number = getbitu(payload, offset, 6);  offset += 6; log(LOG_DEBUG, tab, "TF014_Encryption_sequence_number = %d", TF014_Encryption_sequence_number);
	//SF092 
	return 1;
}

int decode_spartn(spartn_t* spartn) {
	switch (spartn->type)
	{
	case 0:
		return decode_OCB_message(spartn);
		break;
	case 1:
		return decode_HPAC_message(spartn);//warning
		break;
	case 2:
		return decode_GAD_message(spartn);
		break;
	case 3:
		return decode_LPAC_message(spartn);
		break;
	case 4:
		//if (spartn->Subtype == 0) {
		//	decode_Dynamic_Key(spartn);
		//}
		//else if (spartn->Subtype == 1) {
		//	decode_Group_Authentication(spartn);
		//}
		break;
	}
	return 0;
}

int input_spartn_data(spartn_t* spartn, uint8_t data) {
	int tab = 1;
	if (spartn->nbyte == 0) {
		if (data == SPARTN_PREAMB) {
			spartn->buff[spartn->nbyte++] = data;
		}
		return 0;
	}
	int CRC_Len = 0;
	if (spartn->nbyte > 4) {
		// type bit-len byte-len
		//	0		8		1
		//	1		16		2
		//	2		24		3
		//	3		32		4
		CRC_Len = spartn->CRC_type + 1;
	}
	int Time_tag_type_len = 0;
	if (spartn->nbyte > 6) {
		if (spartn->Time_tag_type == 0) {
			Time_tag_type_len = 2;
		}
		else if (spartn->Time_tag_type == 1) {
			Time_tag_type_len = 4;
		}
	}
	int EA_Len = 0;
	int EADL = 0;
	if (spartn->EAF == 1) {
		EA_Len = 2;
		if (spartn->AI > 1) {
			switch (spartn->EAL) {
			case 0:
				EADL = 8;//64
				break;
			case 1:
				EADL = 12;//96
				break;
			case 2:
				EADL = 16;//128
				break;
			case 3:
				EADL = 32;//256
				break;
			case 4:
				EADL = 64;//512
				break;
			}
		}
	}

	spartn->buff[spartn->nbyte++] = data;
	if (spartn->nbyte == 2) {//16
		spartn->type = getbitu(spartn->buff, 8, 7); log(LOG_DEBUG, tab, "type = %d", spartn->type);
	}
	else if (spartn->nbyte == 4) {//32
		spartn->len = getbitu(spartn->buff, 15, 10); log(LOG_DEBUG, tab, "len = %d BYTES", spartn->len);
		spartn->EAF = getbitu(spartn->buff, 25, 1); log(LOG_DEBUG, tab, "EAF = %d", spartn->EAF);
		spartn->CRC_type = getbitu(spartn->buff, 26, 2); log(LOG_DEBUG, tab, "CRC_type = %d", spartn->CRC_type);
		spartn->Frame_CRC = getbitu(spartn->buff, 28, 4);
		char Frame_CRC_Buffer[3] = { 0 };
		bitscopy(Frame_CRC_Buffer, 0, spartn->buff+1,0,20);
		uint8_t Result_Frame_CRC = crc4_itu(Frame_CRC_Buffer, 3);
		log(LOG_DEBUG, tab, "Frame_CRC = %d : %d", spartn->Frame_CRC,Result_Frame_CRC);
		if (spartn->Frame_CRC != Result_Frame_CRC) {
			memset(spartn, 0, sizeof(spartn_t));
			return -1;
		}
	}
	else if (spartn->nbyte == 5) {//40
		spartn->Subtype = getbitu(spartn->buff, 32, 4); log(LOG_DEBUG, tab, "Subtype = %d", spartn->Subtype);
		spartn->Time_tag_type = getbitu(spartn->buff, 36, 1); log(LOG_DEBUG, tab, "Time_tag_type = %d", spartn->Time_tag_type);
	}
	else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + Time_tag_type_len) {//48
		spartn->GNSS_time_type = getbitu(spartn->buff, 37, Time_tag_type_len*8); log(LOG_DEBUG, tab, "GNSS_time_type = %d", spartn->GNSS_time_type);
		spartn->Solution_ID = getbitu(spartn->buff, 37 + Time_tag_type_len * 8, 7); log(LOG_DEBUG, tab, "Solution_ID = %d", spartn->Solution_ID);
		spartn->Solution_processor_ID = getbitu(spartn->buff, 44 + Time_tag_type_len * 8, 4); log(LOG_DEBUG, 1, "Solution_processor_ID = %d", spartn->Solution_processor_ID);
	}
	else if (spartn->EAF == 1 && Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len) {//64
		spartn->Encryption_ID = getbitu(spartn->buff, 48 + Time_tag_type_len * 8, 4); log(LOG_DEBUG, tab, "Encryption_ID = %d", spartn->Encryption_ID);
		spartn->ESN = getbitu(spartn->buff, 52 + Time_tag_type_len * 8, 6); log(LOG_DEBUG, tab, "ESN = %d", spartn->ESN);
		spartn->AI = getbitu(spartn->buff, 55 + Time_tag_type_len * 8, 3); log(LOG_DEBUG, tab, "AI = %d", spartn->AI);
		spartn->EAL = getbitu(spartn->buff, 58 + Time_tag_type_len * 8, 3); log(LOG_DEBUG, tab, "EAL = %d", spartn->EAL);
	}
	else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len) {
		spartn->Payload_offset = 6 + EA_Len + Time_tag_type_len;
		log(LOG_DEBUG, tab, "Payload_offset = %d", spartn->Payload_offset);
	}
	else if (spartn->EAF == 1 && spartn->AI > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len + EADL) {

	}
	else if (Time_tag_type_len > 0 && spartn->nbyte == 6 + EA_Len + Time_tag_type_len + spartn->len + EADL + CRC_Len) {
		spartn->Message_CRC = getbitu(spartn->buff, (spartn->nbyte - CRC_Len) * 8, CRC_Len*8);
		log(LOG_DEBUG, tab, "nbyte = %d", spartn->nbyte);

		uint32_t Result_Message_CRC = crc24_radix(spartn->buff + 1, 5 + EA_Len + Time_tag_type_len + spartn->len + EADL);
		log(LOG_DEBUG, tab, "Message_CRC = %d : %d", spartn->Message_CRC,Result_Message_CRC);
		if (spartn->Message_CRC == Result_Message_CRC) {
			log(LOG_DEBUG, tab, "==========");
			decode_spartn(spartn);
			log(LOG_DEBUG, tab, "==========");
		}
		memset(spartn, 0, sizeof(spartn_t));
		return 1;
	}
	return 0;
}

int main() {
	FILE* fp;
	fp = fopen("RawSpartnPreEncrypt.cap","rb");
	if (fp == NULL) {
		return 0;
	}
	char buff = ' ';
	size_t currentCount = 0;
	size_t frameSize = 0;
	size_t frameCount = 0;
	size_t readCount = 0;
	spartn_t spartn;
	memset(&spartn, 0, sizeof(spartn));
	while (!feof(fp)) {
		memset(&buff, 0, sizeof(buff));
		readCount = fread(&buff, sizeof(char), 1, fp);
		if (readCount < 1)
		{
			/* file error or eof of file */
			break;
		}
		currentCount += readCount;
		frameSize += readCount;
		if (buff == SPARTN_PREAMB) {
			log(LOG_DEBUG, 0, "frameSize = %d",frameSize);
		}
		int ret = input_spartn_data(&spartn, buff);
		if (ret == 1) {
			frameSize = 0;
			frameCount++;
			log(LOG_DEBUG, 0, "frame = %d", frameCount);
			log(LOG_DEBUG, 0, "");
			if (frameCount == 4) {
				break;
			}
		}
	}
	fclose(fp);
	system("pause");
	return 0;
}