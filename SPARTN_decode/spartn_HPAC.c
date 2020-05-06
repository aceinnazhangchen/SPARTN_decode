#include "spartn.h"
#include "log.h"
#include "bits.h"

FILE*  hpac_table_file = NULL;

void log_hpac_title_to_table() {
	table_log_ex(hpac_table_file,"%9s,%4s,%3s,%3s,%3s,%7s,%3s,%7s,%7s,%7s", "Time", "area", "NGp", "Tro", "Ion", "delay", "Tgp", "T00", "T01", "T10");
}

void open_hpac_table_file(const char* filename) {
	if (filename) {
		open_table_file_ex(&hpac_table_file, filename);
	}
	else {
		open_table_file_ex(&hpac_table_file, "../HPAC_message.log");
	}
	log_hpac_title_to_table();
}

void close_hpac_table_file() {
	close_table_file_ex(&hpac_table_file);
}

void log_hpac_to_table(raw_spartn_t* spartn, HPAC_t* hpac) {
	int i,j;
	char sys = ' ';
	if (spartn->Subtype == 0) {
		sys = 'G';
	}
	else {
		sys = 'R';
	}
	uint32_t time = spartn->GNSS_time_type;
	for (i = 0; i < hpac->header.SF030_Area_count; i++) {
		HPAC_area_t* area = &hpac->atmosphere[i].area;
		HPAC_troposphere_t* troposphere = &hpac->atmosphere[i].troposphere;

		HPAC_troposphere_large_t* large_coefficient = &(troposphere->large_coefficient);
		table_log_ex(hpac_table_file, "%9d,%4d,%3d,%3d,%3d,%7.3f,%3d,%7.3f,%7.3f,%7.3f", time, area->SF031_Area_ID, area->SF039_Number_grid_points_present, area->SF040_Tropo, area->SF040_Iono, troposphere->SF043_Area_average_vertical_hydrostatic_delay, 
			troposphere->SF044_Troposphere_polynomial_coefficient_size_indicator,large_coefficient->SF048_T00, large_coefficient->SF049_T01, large_coefficient->SF049_T10, large_coefficient->SF050_T11);
		table_log_ex(hpac_table_file, "%30s %3s,%3s,%7s,%7s,%7s", "", "sat", "Igp", "C00", "C01", "C10");
		HPAC_ionosphere_t* ionosphere = &hpac->atmosphere[i].ionosphere;
		for (j = 0; j < ionosphere->ionosphere_satellite_num; j++){
			HPAC_ionosphere_satellite_t* sat = &ionosphere->ionosphere_satellite[j];
			table_log_ex(hpac_table_file, "%30s %c%02d,%3d,%7.3f,%7.3f,%7.3f", "", sys,sat->PRN_ID, sat->SF056_Ionosphere_satellite_polynomial_block, sat->small_coefficient.SF057_C00, sat->small_coefficient.SF058_C01, sat->small_coefficient.SF058_C10);
		}
	}
	table_log_ex(hpac_table_file, "==============================================================");
}
//Table 6.12 Area data block 
void decode_area_data_block(raw_spartn_t* spartn, HPAC_area_t* area, int tab) {
	uint8_t* payload = spartn->payload;
	area->SF031_Area_ID = getbitu(payload, spartn->offset, 8);  spartn->offset += 8; slog(LOG_DEBUG, tab, "SF031_Area_ID = %d", area->SF031_Area_ID);
	area->SF039_Number_grid_points_present = getbitu(payload, spartn->offset, 7);  spartn->offset += 7; slog(LOG_DEBUG, tab, "SF039_Number_grid_points_present = %d", area->SF039_Number_grid_points_present);
	area->SF040_Tropo = getbitu(payload, spartn->offset, 2);  spartn->offset += 2; slog(LOG_DEBUG, tab, "SF040_Tropo = %d", area->SF040_Tropo);
	area->SF040_Iono = getbitu(payload, spartn->offset, 2);  spartn->offset += 2; slog(LOG_DEBUG, tab, "SF040_Iono = %d", area->SF040_Iono);
}
//Table 6.14 troposphere small coefficient block 
void decode_troposphere_small_coefficient_block(raw_spartn_t* spartn, HPAC_troposphere_t* troposphere, HPAC_troposphere_small_t* small_coefficient, int tab) {
	uint8_t* payload = spartn->payload;
	small_coefficient->SF045_T00 = getbitu(payload, spartn->offset, 7) * 0.004 - 0.252 + 0.252; spartn->offset += 7; slog(LOG_DEBUG, tab, "SF045_T00 = %f", small_coefficient->SF045_T00);//0,1,2
	if (troposphere->SF041_Troposphere_equation_type == 1 || troposphere->SF041_Troposphere_equation_type == 2) {
		small_coefficient->SF046_T01 = getbitu(payload, spartn->offset, 7) * 0.001 - 0.063; spartn->offset += 7; slog(LOG_DEBUG, tab, "SF046_T01 = %f", small_coefficient->SF046_T01);
		small_coefficient->SF046_T10 = getbitu(payload, spartn->offset, 7) * 0.001 - 0.063; spartn->offset += 7; slog(LOG_DEBUG, tab, "SF046_T10 = %f", small_coefficient->SF046_T10);
	}
	if (troposphere->SF041_Troposphere_equation_type == 2) {
		small_coefficient->SF047_T11 = getbitu(payload, spartn->offset, 9)* 0.0002 - 0.0510; spartn->offset += 9; slog(LOG_DEBUG, tab, "SF047_T11 = %f", small_coefficient->SF047_T11 );
	}
}
//Table 6.15 troposphere large coefficient block 
void decode_troposphere_large_coefficient_block(raw_spartn_t* spartn, HPAC_troposphere_t* troposphere, HPAC_troposphere_large_t* large_coefficient, int tab) {
	uint8_t* payload = spartn->payload;
	large_coefficient->SF048_T00 = getbitu(payload, spartn->offset, 9) * 0.004 - 1.020 + 0.252; spartn->offset += 9; slog(LOG_DEBUG, tab, "SF048_T00 = %f", large_coefficient->SF048_T00);//0,1,2
	if (troposphere->SF041_Troposphere_equation_type == 1 || troposphere->SF041_Troposphere_equation_type == 2) {
		large_coefficient->SF049_T01 = getbitu(payload, spartn->offset, 9) * 0.001 - 0.255; spartn->offset += 9; slog(LOG_DEBUG, tab, "SF049_T01 = %f", large_coefficient->SF049_T01);
		large_coefficient->SF049_T10 = getbitu(payload, spartn->offset, 9) * 0.001 - 0.255; spartn->offset += 9; slog(LOG_DEBUG, tab, "SF049_T10 = %f", large_coefficient->SF049_T10);
	}
	if (troposphere->SF041_Troposphere_equation_type == 2) {
		large_coefficient->SF050_T11 = getbitu(payload, spartn->offset, 11) * 0.0002 - 0.2046; spartn->offset += 11; slog(LOG_DEBUG, tab, "SF050_T11 = %f", large_coefficient->SF050_T11);
	}
}
//Table 6.13 Troposphere data block
void decode_troposphere_block(raw_spartn_t* spartn, HPAC_area_t* area, HPAC_troposphere_t* troposphere, int tab) {
	int i;
	uint8_t* payload = spartn->payload;
	//Troposphere polynomial coefficient block 
	if (area->SF040_Tropo == 1 || area->SF040_Tropo == 2) {
		troposphere->SF041_Troposphere_equation_type = getbitu(payload, spartn->offset, 3);  spartn->offset += 3;                                    slog(LOG_DEBUG, tab, "SF041_Troposphere_equation_type = %d", troposphere->SF041_Troposphere_equation_type);
		troposphere->SF042_Troposphere_quality = getbitu(payload, spartn->offset, 3);  spartn->offset += 3;                                          slog(LOG_DEBUG, tab, "SF042_Troposphere_quality = %d", troposphere->SF042_Troposphere_quality);
		troposphere->SF043_Area_average_vertical_hydrostatic_delay = getbitu(payload, spartn->offset, 8)*0.004 - 0.508 + 2.30;  spartn->offset += 8; slog(LOG_DEBUG, tab, "SF043_Area_average_vertical_hydrostatic_delay = %f", troposphere->SF043_Area_average_vertical_hydrostatic_delay);
		troposphere->SF044_Troposphere_polynomial_coefficient_size_indicator = getbitu(payload, spartn->offset, 1);  spartn->offset += 1;            slog(LOG_DEBUG, tab, "SF044_Troposphere_polynomial_coefficient_size_indicator = %d", troposphere->SF044_Troposphere_polynomial_coefficient_size_indicator);
		if (troposphere->SF044_Troposphere_polynomial_coefficient_size_indicator) {
			//Table 6.15
			HPAC_troposphere_large_t* large_coefficient = &(troposphere->large_coefficient);
			decode_troposphere_large_coefficient_block(spartn, troposphere, large_coefficient,tab+1);
		}
		else {
			//Table 6.14
			HPAC_troposphere_small_t* small_coefficient = &(troposphere->small_coefficient);
			decode_troposphere_small_coefficient_block(spartn, troposphere, small_coefficient, tab+1);
		}
	}
	//Troposphere grid block 
	if (area->SF040_Tropo == 2) {
		troposphere->SF051_Troposphere_residual_field_size = getbitu(payload, spartn->offset, 1);  spartn->offset += 1; slog(LOG_DEBUG, tab, "SF051_Troposphere_residual_field_size = %d", troposphere->SF051_Troposphere_residual_field_size);
		if (troposphere->SF051_Troposphere_residual_field_size) {
			//SF053
			for (i = 0; i < area->SF039_Number_grid_points_present; i++) {
				troposphere->SF053[i] = getbitu(payload, spartn->offset, 8); spartn->offset += 8; slog(LOG_DEBUG, tab, "SF053[%d] = %f",i, troposphere->SF053[i] * 0.004 - 0.508);
			}
		}
		else {
			//SF052
			for (i = 0; i < area->SF039_Number_grid_points_present; i++) {
				troposphere->SF052[i] = getbitu(payload, spartn->offset, 6); spartn->offset += 6; slog(LOG_DEBUG, tab, "SF052[%d] = %f",i, troposphere->SF052[i] * 0.004 - 0.124);
			}
		}
	}
}
//Table 6.18 ionosphere small coefficient block 
void decode_ionosphere_small_coefficient_block(raw_spartn_t* spartn, HPAC_ionosphere_t* ionosphere, HPAC_ionosphere_small_t* small_coefficient, int tab) {
	uint8_t* payload = spartn->payload;
	small_coefficient->SF057_C00 = getbitu(payload, spartn->offset, 12) * 0.04 - 81.88; spartn->offset += 12; slog(LOG_DEBUG, tab, "SF057_C00 = %f", small_coefficient->SF057_C00);//0,1,2
	if (ionosphere->SF054_Ionosphere_equation_type == 1 || ionosphere->SF054_Ionosphere_equation_type == 2) {
		small_coefficient->SF058_C01 = getbitu(payload, spartn->offset, 12) * 0.008 - 16.376; spartn->offset += 12; slog(LOG_DEBUG, tab, "SF058_C01 = %f", small_coefficient->SF058_C01);
		small_coefficient->SF058_C10 = getbitu(payload, spartn->offset, 12) * 0.008 - 16.376; spartn->offset += 12; slog(LOG_DEBUG, tab, "SF058_C10 = %f", small_coefficient->SF058_C10);
	}
	if (ionosphere->SF054_Ionosphere_equation_type == 2) {
		small_coefficient->SF059_C11 = getbitu(payload, spartn->offset, 13) * 0.002 - 8.190; spartn->offset += 13; slog(LOG_DEBUG, tab, "SF059_C11 = %f", small_coefficient->SF059_C11);
	}
}
//Table 6.19 ionosphere large coefficient block 
void decode_ionosphere_large_coefficient_block(raw_spartn_t* spartn, HPAC_ionosphere_t* ionosphere, HPAC_ionosphere_large_t* large_coefficient, int tab) {
	uint8_t* payload = spartn->payload;
	large_coefficient->SF060_C00 = getbitu(payload, spartn->offset, 14) * 0.04 - 327.64; spartn->offset += 14; slog(LOG_DEBUG, tab, "SF060_C00 = %f", large_coefficient->SF060_C00);//0,1,2
	if (ionosphere->SF054_Ionosphere_equation_type == 1 || ionosphere->SF054_Ionosphere_equation_type == 2) {
		large_coefficient->SF061_C01 = getbitu(payload, spartn->offset, 14) * 0.008 - 65.528; spartn->offset += 14; slog(LOG_DEBUG, tab, "SF061_C01 = %f", large_coefficient->SF061_C01);
		large_coefficient->SF061_C10 = getbitu(payload, spartn->offset, 14) * 0.008 - 65.528; spartn->offset += 14; slog(LOG_DEBUG, tab, "SF061_C10 = %f", large_coefficient->SF061_C10);
	}
	if (ionosphere->SF054_Ionosphere_equation_type == 2) {
		large_coefficient->SF062_C11 = getbitu(payload, spartn->offset, 15) * 0.002 - 32.766; spartn->offset += 15; slog(LOG_DEBUG, tab, "SF062_C11 = %f", large_coefficient->SF062_C11);
	}
}
//Table 6.17 Ionosphere satellite block
void decode_ionosphere_satellite_block(raw_spartn_t* spartn, HPAC_area_t* area, HPAC_ionosphere_t* ionosphere, HPAC_ionosphere_satellite_t* sat,int tab) {
	int i;
	uint8_t* payload = spartn->payload;
	//Table 6.17 Ionosphere satellite block
	slog(LOG_DEBUG, tab, "PRN_ID = %d", sat->PRN_ID);
	if (area->SF040_Iono == 1 || area->SF040_Iono == 2) {
		sat->SF055_Ionosphere_quality = getbitu(payload, spartn->offset, 4); spartn->offset += 4; slog(LOG_DEBUG, tab, "SF055_Ionosphere_quality = %d", sat->SF055_Ionosphere_quality);
		sat->SF056_Ionosphere_satellite_polynomial_block = getbitu(payload, spartn->offset, 1); spartn->offset += 1; slog(LOG_DEBUG, tab, "SF056_Ionosphere_satellite_polynomial_block = %d", sat->SF056_Ionosphere_satellite_polynomial_block);
		if (sat->SF056_Ionosphere_satellite_polynomial_block) {
			//Table 6.19 ionosphere large coefficient block 
			HPAC_ionosphere_large_t* large_coefficient = &(sat->large_coefficient);
			decode_ionosphere_large_coefficient_block(spartn, ionosphere, large_coefficient,tab+1);
		}
		else {
			//Table 6.18 ionosphere small coefficient block 
			HPAC_ionosphere_small_t* small_coefficient = &(sat->small_coefficient);
			decode_ionosphere_small_coefficient_block(spartn, ionosphere, small_coefficient, tab + 1);
		}
	}
	if (area->SF040_Iono == 2) {
		sat->SF063_Ionosphere_residual_field_size = getbitu(payload, spartn->offset, 2); spartn->offset += 2; slog(LOG_DEBUG, tab, "SF063_Ionosphere_residual_field_size = %d", sat->SF063_Ionosphere_residual_field_size);
		switch (sat->SF063_Ionosphere_residual_field_size) {
		case 0:{
			for (i = 0; i < area->SF039_Number_grid_points_present; i++) {
				sat->ionosphere_residual_slant_delay[i] = getbitu(payload, spartn->offset, 4); spartn->offset += 4; slog(LOG_DEBUG, tab, "SF064[%d] = %f", i, sat->ionosphere_residual_slant_delay[i] * 0.04 - 0.28);
			}
		}
		break;
		case 1:{
			for (i = 0; i < area->SF039_Number_grid_points_present; i++) {
				sat->ionosphere_residual_slant_delay[i] = getbitu(payload, spartn->offset, 7); spartn->offset += 7; slog(LOG_DEBUG, tab, "SF065[%d] = %f", i, sat->ionosphere_residual_slant_delay[i] * 0.04 - 2.52);
			}
		}
		break;
		case 2:{
			for (i = 0; i < area->SF039_Number_grid_points_present; i++) {
				sat->ionosphere_residual_slant_delay[i] = getbitu(payload, spartn->offset, 10); spartn->offset += 10; slog(LOG_DEBUG, tab, "SF066[%d] = %f", i, sat->ionosphere_residual_slant_delay[i] * 0.04 - 20.44);
			}
		}
		break;
		case 3:{
			for (i = 0; i < area->SF039_Number_grid_points_present; i++) {
				sat->ionosphere_residual_slant_delay[i] = getbitu(payload, spartn->offset, 14); spartn->offset += 14; slog(LOG_DEBUG, tab, "SF067[%d] = %f", i, sat->ionosphere_residual_slant_delay[i] * 0.04 - 327.64);
			}
		}
		break;
		}
	}
}
//Table 6.16 Ionosphere block 
void decode_ionosphere_block(raw_spartn_t* spartn, HPAC_area_t* area, HPAC_ionosphere_t* ionosphere, int tab) {
	int i;
	uint8_t* payload = spartn->payload;
	if (area->SF040_Iono == 1 || area->SF040_Iono == 2) {
		ionosphere->SF054_Ionosphere_equation_type = getbitu(payload, spartn->offset, 3); spartn->offset += 3; slog(LOG_DEBUG, tab, "SF054_Ionosphere_equation_type = %d", ionosphere->SF054_Ionosphere_equation_type);
		ionosphere->satellite_mask[64];
		ionosphere->satellite_mask_len = 0;
		if (spartn->Subtype == 0) {
			decode_GPS_satellite_mask(payload, &spartn->offset, ionosphere->satellite_mask, &ionosphere->satellite_mask_len);
		}
		else  if (spartn->Subtype == 1) {
			decode_GLONASS_satellite_mask(payload, &spartn->offset, ionosphere->satellite_mask, &ionosphere->satellite_mask_len);
		}
		//Table 6.17 Ionosphere satellite block (Repeated)
		for (i = 0; i < ionosphere->satellite_mask_len; i++) {
			if (ionosphere->satellite_mask[i]) {
				HPAC_ionosphere_satellite_t* sat = &(ionosphere->ionosphere_satellite[ionosphere->ionosphere_satellite_num]);
				sat->PRN_ID = i + 1;
				decode_ionosphere_satellite_block(spartn, area, ionosphere, sat,tab+1);
				ionosphere->ionosphere_satellite_num++;
			}
		}
	}
}
//Table 6.11 Atmosphere block 
void decode_atmosphere_block(raw_spartn_t* spartn, HPAC_atmosphere_t* atmosphere,int tab) {
	uint8_t* payload = spartn->payload;
	//Table 6.12 Area data block 
	HPAC_area_t* area = &(atmosphere->area);
	decode_area_data_block(spartn, area, tab);
	//Table 6.13 Troposphere data block
	HPAC_troposphere_t* troposphere = &(atmosphere->troposphere);
	decode_troposphere_block(spartn, area, troposphere,tab+1);
	//Table 6.16 Ionosphere block 
	HPAC_ionosphere_t* ionosphere = &(atmosphere->ionosphere);
	decode_ionosphere_block(spartn, area, ionosphere, tab + 1);
}
//Table 6.10 Header block
void decode_Header_block(raw_spartn_t* spartn, HPAC_header_t* hearder,int tab) {
	uint8_t* payload = spartn->payload;
	hearder->SF005_SIOU = getbitu(payload, spartn->offset, 9); spartn->offset += 9; slog(LOG_DEBUG, tab, "SF005_SIOU = %d", hearder->SF005_SIOU);
	hearder->SF069_Reserved = getbitu(payload, spartn->offset, 1);  spartn->offset += 1; slog(LOG_DEBUG, tab, "SF069_Reserved = %d", hearder->SF069_Reserved);
	hearder->SF068_AIOU = getbitu(payload, spartn->offset, 4);  spartn->offset += 4; slog(LOG_DEBUG, tab, "SF068_AIOU = %d", hearder->SF068_AIOU);
	hearder->SF030_Area_count = getbitu(payload, spartn->offset, 5) + 1;  spartn->offset += 5; slog(LOG_DEBUG, tab, "SF030_Area_count = %d", hearder->SF030_Area_count);
}
// SM 1-0/1-1  HPAC messages 
extern int decode_HPAC_message(raw_spartn_t* spartn)
{
	if (!spartn) return 0;
	if (!spartn->spartn_out) return 0;
	if (!spartn->spartn_out->hpac) return 0;
	int i, tab = 2;
	spartn->payload = spartn->buff + spartn->Payload_offset;
	spartn->offset = 0;
	HPAC_t* hpac = spartn->spartn_out->hpac;
	memset(hpac, 0, sizeof(HPAC_t));
	HPAC_header_t* hpac_header = &(hpac->header);
	decode_Header_block(spartn, hpac_header, tab);
	//Table 6.11 Atmosphere block (Repeated)
	HPAC_atmosphere_t* atmosphere = NULL;
	for (i = 0; i < hpac_header->SF030_Area_count; i++) {
		atmosphere = &(hpac->atmosphere[i]);
		decode_atmosphere_block(spartn, atmosphere, tab+1);
		//break;
	}
	transform_spartn_ssr(spartn);
	slog(LOG_DEBUG, tab, "offset = %d bits", spartn->offset);
	slog(LOG_DEBUG, tab, "size of HPAC_t = %d ", sizeof(HPAC_t));
	log_hpac_to_table(spartn, hpac);
	return 1;
}