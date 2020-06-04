#ifndef SPARTN_H
#define SPARTN_H

#include <stdint.h>
#include <memory.h>
#ifdef __cplusplus
extern "C"
{
#endif

#define SPARTN_PREAMB 0x73 
#define RAP_NUM       42
#define SSR_NUM       24
#define AREA_NUM      4
#define VTEC_NUM	  64
#define SAT_MAX		  32
#define HPAC_MAX_ARAE 10

#define SSR_SAP 
#define DAY_SECONDS 86400
//#define USE_Tropo_Iono_2

//=============================
// SM 0-0/0-1  OCB messages 
//=============================
#define SF025_Phase_Bias_Effective_Len 3
#define SF026_Phase_Bias_Effective_Len 2
#define SF027_Phase_Bias_Effective_Len 3
#define SF028_Phase_Bias_Effective_Len 2

#define Bias_Effective_Len 3

typedef struct {
	uint32_t SF005_SIOU;
	uint8_t SF010_EOS;
	uint8_t SF069_Reserved;
	uint8_t SF008_Yaw_present_flag;
	uint8_t SF009_Satellite_reference_datum;
	uint8_t SF016_SF017_Ephemeris_type;
	uint8_t SF011_SF012_satellite_mask[64];
	uint8_t Satellite_mask_len;
} OCB_header_t;

typedef struct {						//Present if SF014_Orbit_block_0 == 1
	uint32_t SF018_SF019_IODE;
	double SF020_radial;
	double SF020_along;
	double SF020_cross;
	uint32_t SF021_Satellite_yaw;		//Present if SF008_Yaw_present_flag == 1
}OCB_orbit_t;

typedef struct {						//Present if SF014_Clock_block_1 == 1
	uint32_t SF022_IODE_continuity;
	double SF020_Clock_correction;
	uint32_t SF024_User_range_error;
}OCB_clock_t;

typedef struct {
	uint32_t SF023_Fix_flag;
	uint32_t SF015_Continuity_indicator;
	double SF020_Phase_bias_correction;
}OCB_Phase_bias_t;

typedef struct {						//Present if SF014_Bias_block_2 == 1
	uint8_t SF025_phase_bias[Bias_Effective_Len];
	OCB_Phase_bias_t Phase_bias[Bias_Effective_Len];
	uint8_t SF027_code_bias[Bias_Effective_Len];
	double SF029_Code_bias_correction[Bias_Effective_Len];
}OCB_GPS_bias_t;

typedef struct {						//Present if SF014_Bias_block_2 == 1
	uint8_t SF026_phase_bias[Bias_Effective_Len];
	OCB_Phase_bias_t Phase_bias[Bias_Effective_Len];
	uint8_t SF028_code_bias[Bias_Effective_Len];
	double SF029_Code_bias_correction[Bias_Effective_Len];
}OCB_GLONASS_bias_t;

typedef struct {						
	uint8_t  PRN_ID;
	uint32_t SF013_DNU;
	uint32_t SF014_Orbit_block_0;
	uint32_t SF014_Clock_block_1;
	uint32_t SF014_Bias_block_2;
	uint32_t SF015_Continuity_indicator;
	OCB_orbit_t orbit;
	OCB_clock_t clock;
	union {
		OCB_GPS_bias_t GPS_bias;
		OCB_GLONASS_bias_t GLONASS_bias;
	};
} OCB_Satellite_t;

//typedef struct {
//	OCB_header_t header;
//	OCB_Satellite_t satellite[SAT_MAX];
//	uint8_t satellite_num;
//}OCB_t;
//=============================
// SM 1-0/1-1  HPAC messages 
//=============================
typedef struct {
	uint16_t SF005_SIOU;
	uint8_t SF069_Reserved;
	uint8_t SF068_AIOU;
	uint8_t SF030_Area_count;
} HPAC_header_t;

typedef struct {
	uint8_t SF031_Area_ID;
	uint8_t SF039_Number_grid_points_present;
	uint8_t SF040_Tropo;
	uint8_t SF040_Iono;
}HPAC_area_t;

typedef struct {
	double SF045_T00;
	double SF046_T01;
	double SF046_T10;
	double SF047_T11;
} HPAC_troposphere_small_t;

typedef struct {
	double SF048_T00;
	double SF049_T01;
	double SF049_T10;
	double SF050_T11;
} HPAC_troposphere_large_t;

typedef struct {
	uint8_t SF041_Troposphere_equation_type;
	uint8_t SF042_Troposphere_quality;
	double SF043_Area_average_vertical_hydrostatic_delay;
	uint8_t SF044_Troposphere_polynomial_coefficient_size_indicator;
	union {
		HPAC_troposphere_small_t small_coefficient;
		HPAC_troposphere_large_t large_coefficient;
	};
	//SF040_Tropo == 2
#ifdef USE_Tropo_Iono_2
	uint8_t SF051_Troposphere_residual_field_size;
	union {
		uint8_t SF052[128];
		uint8_t SF053[128];
	};
#endif // USE_Tropo_Iono_2
}HPAC_troposphere_t;

typedef struct {
	double SF057_C00;
	double SF058_C01;
	double SF058_C10;
	double SF059_C11;
} HPAC_ionosphere_small_t;

typedef struct {
	double SF060_C00;
	double SF061_C01;
	double SF061_C10;
	double SF062_C11;
} HPAC_ionosphere_large_t;

typedef struct {
	uint8_t PRN_ID;
	uint8_t SF055_Ionosphere_quality;
	uint8_t SF056_Ionosphere_satellite_polynomial_block;
	union {
		HPAC_ionosphere_small_t small_coefficient;
		HPAC_ionosphere_large_t large_coefficient;
	};
	//SF040_Iono == 2
#ifdef USE_Tropo_Iono_2
	uint8_t SF063_Ionosphere_residual_field_size;
	uint16_t ionosphere_residual_slant_delay[128];
#endif // USE_Tropo_Iono_2
}HPAC_ionosphere_satellite_t;

typedef struct {
	uint8_t SF054_Ionosphere_equation_type;
	HPAC_ionosphere_satellite_t ionosphere_satellite[SAT_MAX];
	uint8_t ionosphere_satellite_num;
}HPAC_ionosphere_t;

typedef struct {
	HPAC_area_t area;
	HPAC_troposphere_t troposphere;
	HPAC_ionosphere_t ionosphere;
}HPAC_atmosphere_t;

//typedef struct {
//	HPAC_header_t header;
//	HPAC_atmosphere_t atmosphere[HPAC_MAX_ARAE];
//}HPAC_t;
//=============================
// SM 2-0 GAD messages 
//=============================
typedef struct {
	uint8_t SF031_Area_ID;
	double SF032_Area_reference_latitude;
	double SF033_Area_reference_longitude;
	uint8_t SF034_Area_latitude_grid_node_count;
	uint8_t SF035_Area_longitude_grid_node_count;
	double SF036_Area_latitude_grid_node_spacing;
	double SF037_Area_longitude_grid_node_spacing;
}GAD_area_t;

typedef struct {
	uint16_t SF005_SIOU;
	uint8_t SF069_Reserved;
	uint8_t SF068_AIOU;
	uint8_t SF030_Area_count;
}GAD_header_t;

//typedef struct {
//	GAD_header_t header;
//	GAD_area_t areas[32];
//}GAD_t;
//=============================
// SM 3-0 LPAC messages 
//=============================
typedef struct {
	uint8_t SF055_VTEC_quality;
	uint8_t SSF081_VTEC_size_indicator;
	union {
		double SF082_VTEC_residual;
		double SF083_VTEC_residual;
	};
}LPAC_VTEC_t;

typedef struct {
	uint8_t SF072_LPAC_area_ID;
	int16_t SF073_LPAC_area_reference_latitude;
	int16_t SF074_LPAC_area_reference_longitude;
	uint8_t SF075_LPAC_area_latitude_grid_node_count;
	uint8_t SF076_LPAC_area_longitude_grid_node_count;
	uint8_t SF077_LPAC_area_latitude_grid_node_spacing;
	uint8_t SF078_LPAC_area_longitude_grid_node_spacing;
	double SF080_Average_area_VTEC;
	uint8_t  SF079_Grid_node_present_mask[256];
	LPAC_VTEC_t VTEC[64];
}LPAC_area_t;

typedef struct {
	uint16_t SF005_SIOU;
	uint8_t SF069_Reserved;
	uint8_t SF070_Ionosphere_shell_height;
	uint8_t SF071_LPAC_area_count;
}LPAC_header_t;

//typedef struct {
//	LPAC_header_t header;
//	LPAC_area_t areas[4];
//}LPAC_t;

typedef struct {
	int     areaId;
	double  rap_lon;
	double  rap_lat;
	int     nc_lon;
	int     nc_lat;
	double  spa_lon;
	double  spa_lat;
}gad_ssr_t;

typedef struct {                            /* SSR correction type */
    unsigned char prn;
	unsigned char sys;
    unsigned char sat;
    double  t0[6];                         /* epoch time (GPST) {eph,clk,cbias,pbias,tro,ion} */
    int     iod[5];                        /* iod ssr {eph,clk,pbias} */
	double  ure;                           /* user range error*/
	int     fix_flag[3];                   /*0: float, 1:fixed*/
    double  deph[3];                       /* delta orbit {radial,along,cross} (m) */
    double  dclk;                          /* delta clock {c0} (m) */
    double  cbias[3];                      /* code biases (m) */
    double  pbias[3];                      /* phase biases (m) */
    double  yaw_ang;                       /* yaw angle and yaw rate (deg,deg/s) */
	int     areaId[RAP_NUM];
    double  ave_htd;
    double  tro_coef[3 * RAP_NUM];                  /* T00,T01,T10*/
    double  stec_coef[3 * RAP_NUM];                 /* C00,C01.C10*/
    int     rap_num;
} sap_ssr_t;

typedef struct {                            /* SSR correction type */
    double time;
    int    num;
    double corr1[7];                        /*geo,clk,etide,phase windup,shi,trop,sion*/
    double corr2[7];
    double diff [7];
} sap_cor_dif;

typedef struct {
	double time;
	uint8_t  areaId;
	int16_t  rap_lon;
	int16_t  rap_lat;
	uint8_t  nc_lon;
	uint8_t  nc_lat;
	double  spa_lon;
	double  spa_lat;
	double avg_vtec;
	double residual[VTEC_NUM];
} vtec_t;

typedef struct {
	uint32_t type;
	uint32_t Subtype;
	uint32_t len;
	uint32_t time;
	gad_ssr_t ssr_gad[RAP_NUM];
	uint8_t ssr_offset;
    uint8_t eos;
	sap_ssr_t ssr[SSR_NUM];
	vtec_t vtec[AREA_NUM];
	uint16_t is_height;			//Ionosphere shell height
} spartn_t;

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
    uint8_t* payload;
    uint32_t offset;
} raw_spartn_t;

void decode_GPS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint8_t *satellite_mask_len);
void decode_GLONASS_satellite_mask(uint8_t* data, int *pos, uint8_t* satellite_mask, uint8_t *satellite_mask_len);

int decode_OCB_message(raw_spartn_t* spartn, spartn_t* spartn_out);
int decode_HPAC_message(raw_spartn_t* spartn, spartn_t* spartn_out);
int decode_GAD_message(raw_spartn_t* spartn, spartn_t* spartn_out);
int decode_LPAC_message(raw_spartn_t* spartn, spartn_t* spartn_out);

int input_spartn_data(raw_spartn_t* spartn, spartn_t* spartn_out, uint8_t data);

void open_ocb_table_file(const char* filename);
void open_hpac_table_file(const char* filename);
void open_gad_table_file(const char* filename);
void open_lpac_table_file(const char* filename);

void close_ocb_table_file();
void close_hpac_table_file();
void close_gad_table_file();
void close_lpac_table_file();

void ssr_append_ocb_sat(spartn_t* spartn, OCB_Satellite_t* sat_obc);
void ssr_append_hpac_sat(spartn_t * spartn, HPAC_atmosphere_t * atmosphere);
void ssr_append_gad_sat(spartn_t * spartn, GAD_area_t * area);
void ssr_append_lpac_area(spartn_t * spartn, LPAC_area_t * area);

//void transform_spartn_ssr(spartn_t* spartn, OCB_t* ocb, HPAC_t* hpac, GAD_t* gad, LPAC_t* lpac);
void expanded_full_time(raw_spartn_t* raw_spartn);

#ifdef __cplusplus
}
#endif

#endif // !CRC_H
