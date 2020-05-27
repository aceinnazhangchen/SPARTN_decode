#include "spartn.h"
#include "log.h"
#include "bits.h"

FILE*  gad_table_file = NULL;

void log_gad_title_to_table() {
	table_log_ex(gad_table_file, "%9s,%5s,%8s,%8s,%5s,%5s,%7s,%7s", "Time", "area", "lat", "lon", "lat_c", "lon_c", "lat_s", "lon_s");
}

void open_gad_table_file(const char* filename) {
	if (filename) {
		open_table_file_ex(&gad_table_file, filename);
	}
	else {
		open_table_file_ex(&gad_table_file, "../GAD_message.log");
	}
	log_gad_title_to_table();
}

void close_gad_table_file() {
	close_table_file_ex(&gad_table_file);
}

void log_gad_to_table(raw_spartn_t* spartn, GAD_t* gad) {
	int i;
	uint32_t time = spartn->GNSS_time_type;
	for (i = 0; i < gad->header.SF030_Area_count; i++) {
		GAD_area_t* area = &gad->areas[i];
		table_log_ex(gad_table_file, "%9d,%5d,%8.3f,%8.3f,%5d,%5d,%7.3f,%7.3f", time, area->SF031_Area_ID, 
			area->SF032_Area_reference_latitude, area->SF033_Area_reference_longitude,
			area->SF034_Area_latitude_grid_node_count, area->SF035_Area_longitude_grid_node_count,
			area->SF036_Area_latitude_grid_node_spacing, area->SF037_Area_longitude_grid_node_spacing);
	}
	table_log_ex(gad_table_file, "");
}
//Table 6.21 Area definition block 
void decode_Area_definition_block(raw_spartn_t* spartn, GAD_area_t* area, int tab) {
	uint8_t* payload = spartn->payload;
	area->SF031_Area_ID = getbitu(payload, spartn->offset, 8);  spartn->offset += 8; slog(LOG_DEBUG, tab, "SF031_Area_ID = %d", area->SF031_Area_ID);
	area->SF032_Area_reference_latitude = getbitu(payload, spartn->offset, 11)*0.1 - 90;  spartn->offset += 11; slog(LOG_DEBUG, tab, "SF032_Area_reference_latitude = %f", area->SF032_Area_reference_latitude);
	area->SF033_Area_reference_longitude = getbitu(payload, spartn->offset, 12)*0.1 - 180;  spartn->offset += 12; slog(LOG_DEBUG, tab, "SF033_Area_reference_longitude = %f", area->SF033_Area_reference_longitude);
	area->SF034_Area_latitude_grid_node_count = getbitu(payload, spartn->offset, 3);  spartn->offset += 3; slog(LOG_DEBUG, tab, "SF034_Area_latitude_grid_node_count = %d", area->SF034_Area_latitude_grid_node_count);
	area->SF035_Area_longitude_grid_node_count = getbitu(payload, spartn->offset, 3);  spartn->offset += 3; slog(LOG_DEBUG, tab, "SF035_Area_longitude_grid_node_count = %d", area->SF035_Area_longitude_grid_node_count);
	area->SF036_Area_latitude_grid_node_spacing = getbitu(payload, spartn->offset, 5)*0.1 + 0.1;  spartn->offset += 5; slog(LOG_DEBUG, tab, "SF036_Area_latitude_grid_node_spacing = %f", area->SF036_Area_latitude_grid_node_spacing);
	area->SF037_Area_longitude_grid_node_spacing = getbitu(payload, spartn->offset, 5)*0.1 + 0.1;  spartn->offset += 5; slog(LOG_DEBUG, tab, "SF037_Area_longitude_grid_node_spacing = %f", area->SF037_Area_longitude_grid_node_spacing);
}
//Table 6.20 Header block
void decode_GAD_header_block(raw_spartn_t* spartn, GAD_header_t* header,int tab) {
	uint8_t* payload = spartn->payload;
	header->SF005_SIOU = getbitu(payload, spartn->offset, 9); spartn->offset += 9; slog(LOG_DEBUG, tab, "SF005_SIOU = %d", header->SF005_SIOU);
	header->SF069_Reserved = getbitu(payload, spartn->offset, 1);  spartn->offset += 1; slog(LOG_DEBUG, tab, "SF069_Reserved = %d", header->SF069_Reserved);
	header->SF068_AIOU = getbitu(payload, spartn->offset, 4);  spartn->offset += 4; slog(LOG_DEBUG, tab, "SF068_AIOU = %d", header->SF068_AIOU);
	header->SF030_Area_count = getbitu(payload, spartn->offset, 5) + 1;  spartn->offset += 5; slog(LOG_DEBUG, tab, "SF030_Area_count = %d", header->SF030_Area_count);
}
// SM 2-0 GAD messages
extern int decode_GAD_message(raw_spartn_t* spartn)
{
	if (!spartn) return 0;
	if (!spartn->spartn_out) return 0;
	int i, tab = 2;
	spartn->payload = spartn->buff + spartn->Payload_offset;
	spartn->offset = 0;
	GAD_t gad_o = { 0 };
	GAD_t* gad = &gad_o;
	//Table 6.20 Header block
	GAD_header_t* header = &(gad->header);
	decode_GAD_header_block(spartn, header, tab);
	//Table 6.21 (Repeated)
	for (i = 0; i < header->SF030_Area_count; i++) {
		GAD_area_t* area = &(gad->areas[i]);
		decode_Area_definition_block(spartn, area, tab+1);
		//break;
	}
	transform_spartn_ssr(spartn, NULL, NULL, gad, NULL);
	slog(LOG_DEBUG, tab, "offset = %d bits", spartn->offset);
	slog(LOG_DEBUG, tab, "size of GAD_t = %d ", sizeof(GAD_t));
	log_gad_to_table(spartn, gad);
	return 1;
}