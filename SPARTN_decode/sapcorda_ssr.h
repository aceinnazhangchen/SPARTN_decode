#pragma once
#include <stdio.h>
#include "spartn.h"
#include <mutex>

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
private:
	raw_spartn_t m_spartn;
	OCB_t ocb;
	HPAC_t hpac;
	GAD_t gad;
	LPAC_t lpac;
public:
	sapcorda_ssr();
	~sapcorda_ssr();
	void input_ssr_stream(unsigned char* buffer, uint32_t len);
	void input_eph_stream(unsigned char* buffer, uint32_t len);
};

unsigned char* merge_ssr_to_obs(double * rovpos, unsigned char*out_buffer, uint32_t *len);