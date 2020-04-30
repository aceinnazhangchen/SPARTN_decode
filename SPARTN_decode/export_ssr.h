#pragma once
#include <stdint.h>

void input_ssr(unsigned char* buffer, uint32_t len);
void input_eph(unsigned char* buffer, uint32_t len);
void input_gga(char* buffer, unsigned char* out_buffer, uint32_t* len);