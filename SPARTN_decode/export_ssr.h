#pragma once
#include <stdint.h>

void input_ssr(unsigned char* buffer, uint32_t len);
void input_eph(unsigned char* buffer, uint32_t len);
unsigned char* input_gga(char* buffer, uint32_t *len);
unsigned char* merge_ssr_to_obs(double * rovpos, uint32_t *len);