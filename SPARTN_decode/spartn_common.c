#include "spartn.h"
#include "log.h"
#include "bits.h"
#define Leap_Sec 18.0
#define GLO_GPS_TD  10800

void decode_GPS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint8_t *satellite_mask_len) {
	int i,offset = *pos;
	int tab = 2;
	uint8_t SF011_Type = getbitu(data, offset, 2);  offset += 2; slog(LOG_DEBUG, tab, "SF011_Type = %d", SF011_Type);
	uint8_t SF011_Len = 0;
	switch (SF011_Type) {
	case 0:SF011_Len = 32; break;
	case 1:SF011_Len = 44; break;
	case 2:SF011_Len = 56; break;
	case 3:SF011_Len = 64; break;
	}
	*satellite_mask_len = SF011_Len;
	slog(LOG_DEBUG, tab, "SF011_Len = %d", SF011_Len);
	uint8_t SF011[8] = { 0 };
	bitscopy(SF011, 0, data, offset, SF011_Len); offset += SF011_Len;
	bits_to_bytes_array(SF011, satellite_mask, *satellite_mask_len);
	char str_satellite_mask[65] = { 0 };
	for (i = 0; i < *satellite_mask_len; i++) {
		str_satellite_mask[i] = satellite_mask[i] ? '1' : '0';
	}
	slog(LOG_DEBUG, tab, "SF011 = %s", str_satellite_mask);
	*pos = offset;
}

void decode_GLONASS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint8_t *satellite_mask_len) {
	int i,offset = *pos;
	int tab = 2;
	uint8_t SF012_Type = getbitu(data, offset, 2);  offset += 2; slog(LOG_DEBUG, tab, "SF012_Type = %d", SF012_Type);
	uint8_t SF012_Len = 0;
	switch (SF012_Type) {
	case 0:SF012_Len = 24; break;
	case 1:SF012_Len = 36; break;
	case 2:SF012_Len = 48; break;
	case 3:SF012_Len = 63; break;
	}
	*satellite_mask_len = SF012_Len;
	slog(LOG_DEBUG, tab, "SF012_Len = %d", SF012_Len);
	uint8_t SF012[8] = { 0 };
	bitscopy(SF012, 0, data, offset, SF012_Len); offset += SF012_Len;
	bits_to_bytes_array(SF012, satellite_mask, *satellite_mask_len);
	char str_satellite_mask[64] = { 0 };
	for (i = 0; i < *satellite_mask_len; i++) {
		str_satellite_mask[i] = satellite_mask[i] ? '1' : '0';
	}
	slog(LOG_DEBUG, tab, "SF012 = %s", str_satellite_mask);
	*pos = offset;
}

sap_ssr_t* suitable_ssr(spartn_t* spartn, int prn, int sys) {
    int j = 0, n = 0;
    sap_ssr_t* ssr = NULL;
    for (j = 0; j < SSR_NUM; j++) {
        ssr = &spartn->ssr[j];
        if (ssr->prn == 0) {
            spartn->ssr_offset++;
            break;
        }
        if (ssr->prn == prn && ssr->sys == sys) {
            break;
        }
        ssr = NULL;
    }
    double early_time = 0.0;
    if (ssr == NULL) {
        for (j = 0; j < SSR_NUM; j++) {
            ssr = &spartn->ssr[j];
            if (early_time == 0.0) {
                early_time = ssr->t0[1];
            }
            else {
                if (early_time > ssr->t0[1]) {
                    early_time = ssr->t0[1];
                    n = j;
                }
            }
        }
        ssr = &spartn->ssr[n];
        memset(ssr, 0, sizeof(sap_ssr_t));
    }
    return ssr;
}

void ssr_append_ocb_sat(spartn_t* spartn, OCB_Satellite_t* sat_obc) {
	int j = 0;

	sap_ssr_t* ssr = suitable_ssr(spartn, sat_obc->PRN_ID, spartn->Subtype);
	ssr->prn = sat_obc->PRN_ID;
	ssr->sys = spartn->Subtype;
	ssr->sat = ssr->sys * 40 + ssr->prn;

	if (sat_obc->orbit.SF018_SF019_IODE != 0) {
		ssr->t0[0] = (double)spartn->time + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
		if (ssr->sys == 1 && ssr->t0[0] < 0) ssr->t0[0] += DAY_SECONDS / 2;
		ssr->iod[0] = sat_obc->orbit.SF018_SF019_IODE;
		ssr->deph[0] = sat_obc->orbit.SF020_radial;
		ssr->deph[1] = sat_obc->orbit.SF020_along;
		ssr->deph[2] = sat_obc->orbit.SF020_cross;
		ssr->yaw_ang = sat_obc->orbit.SF021_Satellite_yaw;
	}

	ssr->t0[1] = (double)spartn->time + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
	if (ssr->sys == 1 && ssr->t0[1] < 0)   ssr->t0[1] += DAY_SECONDS / 2;
	ssr->iod[1] = sat_obc->clock.SF022_IODE_continuity;
	ssr->ure = sat_obc->clock.SF024_User_range_error;
	ssr->dclk = sat_obc->clock.SF020_Clock_correction;

	int update_num = 0;
	for (j = 0; j < Bias_Effective_Len; ++j) {
		if (sat_obc->GPS_bias.SF027_code_bias[j] == 1) {
			ssr->cbias[j] = sat_obc->GPS_bias.SF029_Code_bias_correction[j];
			update_num++;
		}
	}
	if (update_num > 0) {
		ssr->t0[2] = (double)spartn->time + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
		if (ssr->sys == 1 && ssr->t0[2] < 0)   ssr->t0[2] += DAY_SECONDS / 2;
	}
	update_num = 0;
	for (j = 0; j < Bias_Effective_Len; ++j) {
		if (sat_obc->GPS_bias.SF025_phase_bias[j] == 1) {
			ssr->iod[j + 2] = sat_obc->GPS_bias.Phase_bias[j].SF015_Continuity_indicator;
			ssr->fix_flag[j] = sat_obc->GPS_bias.Phase_bias[j].SF023_Fix_flag;
			ssr->pbias[j] = sat_obc->GPS_bias.Phase_bias[j].SF020_Phase_bias_correction;
			update_num++;
		}
	}
	if (update_num > 0) {
		ssr->t0[3] = (double)spartn->time + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
		if (ssr->sys == 1 && ssr->t0[3] < 0)   ssr->t0[3] += DAY_SECONDS / 2;
	}
}

void ssr_append_hpac_sat(spartn_t* spartn, HPAC_atmosphere_t* atmosphere) {
	int j = 0, n = 0, m = 0;
	sap_ssr_t* ssr = NULL;
	for (j = 0; j < SSR_NUM; ++j) {
		ssr = &spartn->ssr[j];
		ssr->rap_num = 0;
		if (ssr->prn == 0) break;
		if (ssr->sys != spartn->Subtype) continue;
		ssr->t0[4] = (double)spartn->time + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
		ssr->t0[5] = (double)spartn->time + (ssr->sys ? (-GLO_GPS_TD + Leap_Sec) : 0);
		if (ssr->sys == 1)
		{
			if (ssr->t0[4] < 0) ssr->t0[4] += DAY_SECONDS / 2;
			if (ssr->t0[5] < 0) ssr->t0[5] += DAY_SECONDS / 2;
		}

		for (m = 0; m < RAP_NUM; ++m) {
			if (ssr->areaId[m] == atmosphere->area.SF031_Area_ID) break;
			if (ssr->areaId[m] == 0) break;
		}

		ssr->areaId[m] = atmosphere->area.SF031_Area_ID;
		ssr->ave_htd = atmosphere->troposphere.SF043_Area_average_vertical_hydrostatic_delay;
		ssr->tro_coef[(m + ssr->rap_num) * 3 + 0] = atmosphere->troposphere.small_coefficient.SF045_T00;
		ssr->tro_coef[(m + ssr->rap_num) * 3 + 1] = atmosphere->troposphere.small_coefficient.SF046_T01;
		ssr->tro_coef[(m + ssr->rap_num) * 3 + 2] = atmosphere->troposphere.small_coefficient.SF046_T10;
		for (n = 0; n < atmosphere->ionosphere.ionosphere_satellite_num && n < SAT_MAX; ++n) {
			if (ssr->prn == atmosphere->ionosphere.ionosphere_satellite[n].PRN_ID) {
				ssr->stec_coef[(m + ssr->rap_num) * 3 + 0] = atmosphere->ionosphere.ionosphere_satellite[n].small_coefficient.SF057_C00;
				ssr->stec_coef[(m + ssr->rap_num) * 3 + 1] = atmosphere->ionosphere.ionosphere_satellite[n].small_coefficient.SF058_C01;
				ssr->stec_coef[(m + ssr->rap_num) * 3 + 2] = atmosphere->ionosphere.ionosphere_satellite[n].small_coefficient.SF058_C10;
				break;
			}
		}
	}
}

void ssr_append_gad_sat(spartn_t* spartn, GAD_area_t* area) {
	int j = 0;
	gad_ssr_t* ssr_gad = NULL;
	for (j = 0; j < RAP_NUM; ++j) {
		ssr_gad = &spartn->ssr_gad[j];
		if (ssr_gad->areaId == area->SF031_Area_ID) {
			break;
		}
		if (ssr_gad->areaId == 0) {
			break;
		}
	}
	ssr_gad->areaId = area->SF031_Area_ID;
	ssr_gad->rap_lon = area->SF033_Area_reference_longitude;
	ssr_gad->rap_lat = area->SF032_Area_reference_latitude;
	ssr_gad->nc_lon = area->SF035_Area_longitude_grid_node_count;
	ssr_gad->nc_lat = area->SF034_Area_latitude_grid_node_count;
	ssr_gad->spa_lon = area->SF037_Area_longitude_grid_node_spacing;
	ssr_gad->spa_lat = area->SF036_Area_latitude_grid_node_spacing;
}

void ssr_append_lpac_area(spartn_t * spartn, LPAC_area_t * area)
{
	int j = 0;
	vtec_t* vtec = NULL;
	for (j = 0; j < AREA_NUM; ++j) {
		vtec = &spartn->vtec[j];
		if (vtec->areaId == area->SF072_LPAC_area_ID) {
			break;
		}
		if (vtec->areaId == 0) {
			break;
		}
	}

	vtec->areaId = area->SF072_LPAC_area_ID;
	vtec->rap_lon = area->SF073_LPAC_area_reference_latitude;
	vtec->rap_lat = area->SF074_LPAC_area_reference_longitude;
	vtec->nc_lon = area->SF075_LPAC_area_latitude_grid_node_count;
	vtec->nc_lat = area->SF076_LPAC_area_longitude_grid_node_count;
	vtec->spa_lon = area->SF077_LPAC_area_latitude_grid_node_spacing;
	vtec->spa_lat = area->SF078_LPAC_area_longitude_grid_node_spacing;
	vtec->avg_vtec = area->SF080_Average_area_VTEC;
	
	for (j = 0; j < VTEC_NUM; ++j) {
		vtec->residual[j] = area->VTEC[j].SF082_VTEC_residual;
	}
}

//void ocb_to_ssr(spartn_t* spartn, OCB_t* ocb) {
//	if (!spartn) return;
//	if (!ocb) return;
//	int i = 0;
//	OCB_Satellite_t* sat_obc = NULL;
//	for (i = 0; i < ocb->satellite_num && i < SAT_MAX; ++i) {
//		sat_obc = &ocb->satellite[i];
//		ssr_append_ocb_sat(spartn, sat_obc);
//	}
//	spartn->eos = ocb->header.SF010_EOS;
//}

//void hpac_to_ssr(spartn_t* spartn, HPAC_t* hpac) {
//	if (!spartn) return;
//	if (!hpac) return;
//	int i = 0;
//	sap_ssr_t* ssr = NULL;
//	for (i = 0; i < hpac->header.SF030_Area_count; ++i) {
//		//ssr = &spartn->ssr[j];
//		ssr_append_hpac_sat(spartn, &hpac->atmosphere[i]);
//	}
//	//ssr->rap_num += hpac->header.SF030_Area_count;
//}

//void gad_to_ssr(spartn_t* spartn, GAD_t* gad) {
//	if (!spartn) return;
//	if (!gad) return;
//
//	int i = 0;
//	gad_ssr_t* ssr_gad = NULL;
//	for (i = 0; i < gad->header.SF030_Area_count; ++i) {
//		ssr_append_gad_sat(spartn, &gad->areas[i]);
//	}
//}

//void transform_spartn_ssr(spartn_t* spartn, OCB_t* ocb, HPAC_t* hpac, GAD_t* gad, LPAC_t* lpac)
//{
//	if (!spartn) return;
//    int i = 0, j = 0, n = 0, m=0;
//
//    if (spartn->type == 0) {
//		ocb_to_ssr(spartn, ocb);
//    }
//    else if (spartn->type == 1) {
//		hpac_to_ssr(spartn, hpac);
//    }
//    else if (spartn->type == 2) {
//		gad_to_ssr(spartn, gad);
//    }
//}

void expanded_full_time(raw_spartn_t* raw_spartn) {
    //to full time
    if (raw_spartn->Subtype == 0) {
    	if (raw_spartn->Time_tag_type) {
            raw_spartn->GNSS_time_type = raw_spartn->GNSS_time_type % (DAY_SECONDS);
    	}
    	else {
    		//if (GPS_DAYS > 0) {
    		//	raw_spartn->GNSS_time_type = GPS_DAYS * DAY_SECOND + raw_spartn->GNSS_time_type;
    		//	raw_spartn->Time_tag_type = 1;//modify Time_tag_type
    		//}
    	}
    }
    else if (raw_spartn->Subtype == 1) {//
    	if (raw_spartn->Time_tag_type) {
            raw_spartn->GNSS_time_type = raw_spartn->GNSS_time_type % (DAY_SECONDS);
    	}
    	else {
    		//if (GLONASS_DAYS > 0) {
    		//	raw_spartn->GNSS_time_type = GLONASS_DAYS * DAY_SECOND + raw_spartn->GNSS_time_type;
    		//	raw_spartn->Time_tag_type = 1;//modify Time_tag_type
    		//}
    	}
    }
    //raw_spartn->GNSS_time_type += 1262304000;
}