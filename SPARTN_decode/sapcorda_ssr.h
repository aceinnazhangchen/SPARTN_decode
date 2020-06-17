#pragma once
#include <stdio.h>
#include "spartn.h"
#include <mutex>
#include "rtcm.h"
#include <vector>
#include <map>
using namespace std;

class sapcorda_ssr
{
public:
	static sapcorda_ssr* getInstance();
	static void createInstance();
private:
	static sapcorda_ssr* m_instance;
	static std::once_flag m_flag;
public:
	gnss_rtcm_t m_rtcm;
	obs_t m_obs_vrs;
	spartn_t m_spartn_out;
	map<int, eph_t> m_last_eph_map;
	map<int, geph_t> m_last_geph_map;
	map<int, sap_ssr_t> m_last_ssr_map;
	//vector<eph_t> m_last_eph;
	//vector<geph_t> m_last_geph;
	//vector<sap_ssr_t> m_last_ssr;
private:
	raw_spartn_t m_spartn;
	FILE* m_fLOG;
public:
	sapcorda_ssr();
	~sapcorda_ssr();
	void input_ssr_stream(unsigned char* buffer, uint32_t len);
	void input_eph_stream(unsigned char* buffer, uint32_t len);
	void save_last_eph(nav_t * last_nav, nav_t * nav);
	void save_last_geph(nav_t * last_nav, nav_t * nav);
	void save_last_ssr(sap_ssr_t * last_ssr, uint8_t ssr_offset, spartn_t * spartn);
	unsigned char* merge_ssr_to_obs(double * rovpos, unsigned char*out_buffer, uint32_t *len);
};

