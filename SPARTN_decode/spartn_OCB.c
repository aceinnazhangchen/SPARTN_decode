#include "spartn.h"
#include "log.h"

FILE*  ocb_table_file = NULL;

void open_ocb_table_file(const char* filename) {
	if (filename) {
		open_table_file_ex(&ocb_table_file, filename);
	}
	else {
		open_table_file_ex(&ocb_table_file, "../OCB_message.log");
	}
}

void close_ocb_table_file() {
	close_table_file_ex(&ocb_table_file);
}

void log_ocb_to_table(raw_spartn_t* spartn, OCB_t* ocb) {
	int i;
	uint32_t time = spartn->GNSS_time_type;
	if (spartn->Subtype == 0) {
		table_log_ex(ocb_table_file, "%d%8s, %3s,%3s,%9s,%9s,%9s,%9s,%9s,%9s,%9s,%9s,%9s,%9s", ocb->header.SF010_EOS, "Time", "Sat", "Iod", "Ad", "Cd", "Rd", "Clk", "L1C", "L2W", "L2L", "C1C", "C2W", "C2L");
	}
	else {
		table_log_ex(ocb_table_file, "%d%8s, %3s,%3s,%9s,%9s,%9s,%9s,%9s,%9s,%9s %9s,%9s %9s", ocb->header.SF010_EOS, "Time", "Sat", "Iod", "Ad", "Cd", "Rd", "Clk", "L1C", "L2C", "", "C1C", "C2C", "");
	}
	for (i = 0; i < ocb->satellite_num; i++) {
		if (spartn->Subtype == 0) {
			table_log_ex(ocb_table_file, "%9d, G%02d,%3d,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f", time, ocb->satellite[i].PRN_ID, ocb->satellite[i].orbit.SF018_SF019_IODE,
				ocb->satellite[i].orbit.SF020_along, ocb->satellite[i].orbit.SF020_cross, ocb->satellite[i].orbit.SF020_radial,ocb->satellite[i].clock.SF020_Clock_correction, 
				ocb->satellite[i].GPS_bias.Phase_bias[0].SF020_Phase_bias_correction, ocb->satellite[i].GPS_bias.Phase_bias[1].SF020_Phase_bias_correction, ocb->satellite[i].GPS_bias.Phase_bias[2].SF020_Phase_bias_correction,
				ocb->satellite[i].GPS_bias.SF029_Code_bias_correction[0], ocb->satellite[i].GPS_bias.SF029_Code_bias_correction[1], ocb->satellite[i].GPS_bias.SF029_Code_bias_correction[2]);
		}
		else {
			table_log_ex(ocb_table_file, "%9d, R%02d,%3d,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9s %9.3f,%9.3f %9s", time, ocb->satellite[i].PRN_ID, ocb->satellite[i].orbit.SF018_SF019_IODE,
				ocb->satellite[i].orbit.SF020_along, ocb->satellite[i].orbit.SF020_cross, ocb->satellite[i].orbit.SF020_radial,ocb->satellite[i].clock.SF020_Clock_correction, 
				ocb->satellite[i].GLONASS_bias.Phase_bias[0].SF020_Phase_bias_correction, ocb->satellite[i].GLONASS_bias.Phase_bias[1].SF020_Phase_bias_correction, "",
				ocb->satellite[i].GLONASS_bias.SF029_Code_bias_correction[0], ocb->satellite[i].GLONASS_bias.SF029_Code_bias_correction[1], "");
		}
	}
}

void decode_bias_mask(uint8_t* data, int *pos, uint8_t *mask_array, uint32_t effective_len, uint32_t subType) {
	int i, offset = *pos;
	int tab = 4;
	uint32_t len_flag = getbitu(data, offset, 1);  offset += 1; slog(LOG_DEBUG, tab, "len_flag = %d", len_flag);
	uint32_t max_len = 0;
	if (len_flag == 0) {
		max_len = subType ? 5 : 6;
	}
	else if (len_flag == 1) {
		max_len = subType ? 9 : 11;
	}

	for (i = 0; i < max_len; i++) {
		mask_array[i] = getbitu(data, offset, 1);  offset += 1; slog(LOG_DEBUG, tab, "bias_mask_array[%d] = %d", i, mask_array[i]);
		if (i == effective_len - 1) {
			offset += (max_len - effective_len);
			break; 
		}//just have 3 or 2 type now
	}
	*pos = offset;
}
//Table 6.9 phase bias block 
void decode_phase_bias_block(raw_spartn_t* spartn, uint8_t* mask_array, OCB_Phase_bias_t* bias_array, uint32_t effective_len,int tab) {
	uint8_t* payload = spartn->payload;
	int i;
	for (i = 0; i < effective_len; i++) {
		if (mask_array[i] == 1) {
			bias_array[i].SF023_Fix_flag = getbitu(payload, spartn->offset, 1);  spartn->offset += 1; slog(LOG_DEBUG, tab, "SF023_Fix_flag = %d", bias_array[i].SF023_Fix_flag);
			bias_array[i].SF015_Continuity_indicator = getbitu(payload, spartn->offset, 3);  spartn->offset += 3; slog(LOG_DEBUG, tab, "SF015_Continuity_indicator = %d", bias_array[i].SF015_Continuity_indicator);
			bias_array[i].SF020_Phase_bias_correction = getbitu(payload, spartn->offset, 14)*0.002 - 16.382;  spartn->offset += 14; slog(LOG_DEBUG, tab, "SF020_Phase_bias_correction = %f", bias_array[i].SF020_Phase_bias_correction);
		}
	}
}
//Code bias correction
void decode_code_bias_correction(raw_spartn_t* spartn, uint8_t *mask_array, double *bias_array, uint32_t effective_len, int tab) {
	uint8_t* payload = spartn->payload;
	int i;
	for (i = 0; i < effective_len; i++) {
		if (mask_array[i] == 1) {
			bias_array[i] = getbitu(payload, spartn->offset, 11)*0.02 - 20.46;  spartn->offset += 11; slog(LOG_DEBUG, tab, "SF029_Code_bias_correction = %f", bias_array[i] );
		}
	}
}
//Table 6.5 orbit block 
void decode_orbit_block(raw_spartn_t* spartn, OCB_orbit_t* orbit, uint32_t SF008_Yaw_present_flag, int tab) {
	uint8_t* payload = spartn->payload;
	if (spartn->Subtype == 0) {
		orbit->SF018_SF019_IODE = getbitu(payload, spartn->offset, 8);  spartn->offset += 8; slog(LOG_DEBUG, tab, "SF018_IODE = %d", orbit->SF018_SF019_IODE);
	}
	else if (spartn->Subtype == 1) {
		orbit->SF018_SF019_IODE = getbitu(payload, spartn->offset, 7);  spartn->offset += 7; slog(LOG_DEBUG, tab, "SF019_IODE = %d", orbit->SF018_SF019_IODE);
	}
	orbit->SF020_radial = getbitu(payload, spartn->offset, 14)*0.002 - 16.382;  spartn->offset += 14; slog(LOG_DEBUG, tab, "SF020_radial = %f", orbit->SF020_radial);
	orbit->SF020_along = getbitu(payload, spartn->offset, 14)*0.002 - 16.382;  spartn->offset += 14; slog(LOG_DEBUG, tab, "SF020_along = %f", orbit->SF020_along);
	orbit->SF020_cross = getbitu(payload, spartn->offset, 14)*0.002 - 16.382;  spartn->offset += 14; slog(LOG_DEBUG, tab, "SF020_cross = %f", orbit->SF020_cross);
	if (SF008_Yaw_present_flag == 1) {
		orbit->SF021_Satellite_yaw = getbitu(payload, spartn->offset, 6);  spartn->offset += 6; slog(LOG_DEBUG, tab, "SF021 = %d", orbit->SF021_Satellite_yaw * 6);
	}
}
//Table 6.6 clock block 
void decode_clock_block(raw_spartn_t* spartn, OCB_clock_t* clock, int tab) {
	uint8_t* payload = spartn->payload;
	clock->SF022_IODE_continuity = getbitu(payload, spartn->offset, 3);  spartn->offset += 3; slog(LOG_DEBUG, tab, "SF022_IODE_continuity = %d", clock->SF022_IODE_continuity);
	clock->SF020_Clock_correction = getbitu(payload, spartn->offset, 14)*0.002 - 16.382;  spartn->offset += 14; slog(LOG_DEBUG, tab, "SF020_Clock_correction = %f", clock->SF020_Clock_correction);
	clock->SF024_User_range_error = getbitu(payload, spartn->offset, 3);  spartn->offset += 3; slog(LOG_DEBUG, tab, "SF024_User_range_error = %d", clock->SF024_User_range_error);
}
//Table 6.7 GPS bias block
void decode_GPS_bias_block(raw_spartn_t* spartn, OCB_GPS_bias_t* GPS_bias, int tab) {
	uint8_t* payload = spartn->payload;
	decode_bias_mask(payload, &spartn->offset, GPS_bias->SF025_phase_bias, SF025_Phase_Bias_Effective_Len, spartn->Subtype);
	//Table 6.9 Phase bias block (Repeated)
	decode_phase_bias_block(spartn,GPS_bias->SF025_phase_bias, GPS_bias->Phase_bias, SF025_Phase_Bias_Effective_Len, tab + 1);
	//SF027
	decode_bias_mask(payload, &spartn->offset, GPS_bias->SF027_code_bias, SF027_Phase_Bias_Effective_Len, spartn->Subtype);
	decode_code_bias_correction(spartn, GPS_bias->SF027_code_bias, GPS_bias->SF029_Code_bias_correction, SF027_Phase_Bias_Effective_Len, tab + 1);
}
//Table 6.8 GLONASS bias block
void decode_GLONASS_bias_block(raw_spartn_t* spartn, OCB_GLONASS_bias_t* GLONASS_bias, int tab) {
	uint8_t* payload = spartn->payload;
	decode_bias_mask(payload, &spartn->offset, GLONASS_bias->SF026_phase_bias, SF026_Phase_Bias_Effective_Len, spartn->Subtype);
	//Table 6.9 Phase bias block (Repeated)
	decode_phase_bias_block(spartn, GLONASS_bias->SF026_phase_bias, GLONASS_bias->Phase_bias, SF026_Phase_Bias_Effective_Len, tab + 1);
	//SF028
	decode_bias_mask(payload, &spartn->offset, GLONASS_bias->SF028_code_bias, SF028_Phase_Bias_Effective_Len, spartn->Subtype);
	decode_code_bias_correction(spartn, GLONASS_bias->SF028_code_bias, GLONASS_bias->SF029_Code_bias_correction, SF028_Phase_Bias_Effective_Len, tab + 1);
}
//Table 6.4 satellite block
void decode_satellite_block(raw_spartn_t* spartn, OCB_Satellite_t* sat, uint32_t SF008_Yaw_present_flag, int tab) {
	uint8_t* payload = spartn->payload;
	slog(LOG_DEBUG, tab, "PRN_ID = %d", sat->PRN_ID);
	sat->SF013_DNU = getbitu(payload, spartn->offset, 1);  spartn->offset += 1; slog(LOG_DEBUG, tab, "SF013_DNU = %d", sat->SF013_DNU);
	sat->SF014_Orbit_block_0 = getbitu(payload, spartn->offset, 1);  spartn->offset += 1; slog(LOG_DEBUG, tab, "SF014_Orbit_block_0 = %d", sat->SF014_Orbit_block_0);
	sat->SF014_Clock_block_1 = getbitu(payload, spartn->offset, 1);  spartn->offset += 1; slog(LOG_DEBUG, tab, "SF014_Clock_block_1 = %d", sat->SF014_Clock_block_1);
	sat->SF014_Bias_block_2 = getbitu(payload, spartn->offset, 1);  spartn->offset += 1; slog(LOG_DEBUG, tab, "SF014_Bias_block_2 = %d", sat->SF014_Bias_block_2);
	sat->SF015_Continuity_indicator = getbitu(payload, spartn->offset, 3);  spartn->offset += 3; slog(LOG_DEBUG, tab, "SF015 = %d", sat->SF015_Continuity_indicator);
	//Table 6.5 orbit block 
	if (sat->SF014_Orbit_block_0) {
		decode_orbit_block(spartn, &sat->orbit, SF008_Yaw_present_flag, tab+1);
	}
	//Table 6.6 clock block 
	if (sat->SF014_Clock_block_1) {
		decode_clock_block(spartn, &sat->clock, tab+1);
	}
	//bias block
	if (sat->SF014_Bias_block_2) {
		if (spartn->Subtype == 0) {
			//Table 6.7 GPS bias block
			decode_GPS_bias_block(spartn, &sat->GPS_bias, tab + 1);
		}
		else if (spartn->Subtype == 1) {
			//Table 6.8 GLONASS bias block
			decode_GLONASS_bias_block(spartn, &sat->GLONASS_bias, tab + 1);
		}
	}
}
//Table 6.3 Header block 
void decode_OCB_hearder(raw_spartn_t* spartn, OCB_header_t* ocb_header,int tab) {
	uint8_t* payload = spartn->payload;
	ocb_header->SF005_SIOU = getbitu(payload, spartn->offset, 9); spartn->offset += 9; slog(LOG_DEBUG, tab, "SF005_SIOU = %d", ocb_header->SF005_SIOU);
	ocb_header->SF010_EOS = getbitu(payload, spartn->offset, 1); spartn->offset += 1; slog(LOG_DEBUG, tab, "SF010_EOS = %d", ocb_header->SF010_EOS);
	ocb_header->SF069_Reserved = getbitu(payload, spartn->offset, 1); spartn->offset += 1; slog(LOG_DEBUG, tab, "SF069_Reserved = %d", ocb_header->SF069_Reserved);
	ocb_header->SF008_Yaw_present_flag = getbitu(payload, spartn->offset, 1); spartn->offset += 1; slog(LOG_DEBUG, tab, "SF008_Yaw_present_flag = %d", ocb_header->SF008_Yaw_present_flag);
	ocb_header->SF009_Satellite_reference_datum = getbitu(payload, spartn->offset, 1); spartn->offset += 1; slog(LOG_DEBUG, tab, "SF009_Satellite_reference_datum = %d", ocb_header->SF009_Satellite_reference_datum);
	ocb_header->SF016_SF017_Ephemeris_type = getbitu(payload, spartn->offset, 2); spartn->offset += 2; slog(LOG_DEBUG, tab, "SF016_SF017_Ephemeris_type = %d", ocb_header->SF016_SF017_Ephemeris_type);
	if (spartn->Subtype == 0) {
		decode_GPS_satellite_mask(payload, &spartn->offset, ocb_header->SF011_SF012_satellite_mask, &ocb_header->Satellite_mask_len);
	}
	else if (spartn->Subtype == 1) {
		decode_GLONASS_satellite_mask(payload, &spartn->offset, ocb_header->SF011_SF012_satellite_mask, &ocb_header->Satellite_mask_len);
	}
}
// SM 0-0/0-1  OCB messages 
int decode_OCB_message(raw_spartn_t* spartn) 
{
	if (!spartn) return 0;
	if (!spartn->spartn_out) return 0;
	if (!spartn->spartn_out->ocb) return 0;
	int i,tab = 2;
	spartn->payload = spartn->buff + spartn->Payload_offset;
	spartn->offset = 0;
	OCB_t* ocb= spartn->spartn_out->ocb;
	memset(ocb, 0, sizeof(OCB_t));
	OCB_header_t* ocb_header = &ocb->header;
	//Table 6.3 Header block 
	decode_OCB_hearder(spartn,ocb_header, tab);
	//Table 6.4 Satellite block (Repeated) 
	for (i = 0; i < ocb_header->Satellite_mask_len; i++) {
		if (ocb_header->SF011_SF012_satellite_mask[i]) {
			ocb->satellite[ocb->satellite_num].PRN_ID = i + 1;
			decode_satellite_block(spartn, &ocb->satellite[ocb->satellite_num], ocb_header->SF008_Yaw_present_flag, tab+1);
			ocb->satellite_num++;
		}
	}
	transform_spartn_ssr(spartn);
	slog(LOG_DEBUG, tab, "offset = %d bits", spartn->offset);
	slog(LOG_DEBUG, tab, "size of OCB_t = %d ", sizeof(OCB_t));
	log_ocb_to_table(spartn, ocb);
	return 1;
}