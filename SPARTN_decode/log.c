#include "log.h"
#include <stdarg.h>

FILE*  table_file = NULL;

void log(int level, int tab, const char *format, ...) {
#ifdef LOG
	if (level >= LOG_DEBUG) {
		va_list ap;
		int i;
		char buffer[2048] = { 0 };
		char tab_buffer[128] = { 0 };
		va_start(ap, format); vsprintf(buffer, format, ap); va_end(ap);
		for (i = 0; i < tab; i++) {
			tab_buffer[i] = ' ';
		}
		printf("%s%s \n", tab_buffer, buffer);
	}
#endif // LOG
}

void open_table_file(const char* filename) {
#ifdef TABLE_LOG
	table_file = fopen(filename, "w");
#endif // TABLE_LOG
}

void close_table_file() {
#ifdef TABLE_LOG
	if (table_file) {
		fclose(table_file);
		table_file = NULL;
	}
#endif // TABLE_LOG
}

void table_log(const char *format, ...) {
#ifdef TABLE_LOG
	va_list ap;
	char buffer[2048] = { 0 };
	va_start(ap, format); vsprintf(buffer, format, ap); va_end(ap);
	if (table_file) {
		fprintf(table_file, "%s \n", buffer);
	}
	else {
		printf("%s \n", buffer);
	}
#endif // TABLE_LOG
}

