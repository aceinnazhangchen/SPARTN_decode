#ifndef _LOG_H_
#define _LOG_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

//#define LOG
//#define TABLE_LOG

#define LOG_DEBUG 1
#define LOG_INFO 2
void slog(int level, int tab, const char *format, ...);
void open_table_file(const char* filename);
void close_table_file();
void table_log(const char *format, ...);

void open_table_file_ex(FILE** file, const char* filename);
void close_table_file_ex(FILE** file);
void table_log_ex(FILE * file, const char * format, ...);

#ifdef __cplusplus
}
#endif
#endif