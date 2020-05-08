#ifndef EXPORT_SSR
#define EXPORT_SSR

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


	void input_ssr(unsigned char* buffer, uint32_t len);
	void input_eph(unsigned char* buffer, uint32_t len);
	void input_gga(char* buffer, unsigned char* out_buffer, uint32_t* len);

	void input_ssr_test(unsigned char* buffer, uint32_t len);
	void input_gga_test(char* buffer, unsigned char* out_buffer, uint32_t* len);

#ifdef __cplusplus
}
#endif

#endif // ! EXPORT_SSR


