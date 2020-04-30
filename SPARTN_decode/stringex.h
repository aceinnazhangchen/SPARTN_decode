#ifndef COMMON_FUNCTION_H
#define COMMON_FUNCTION_H

#include <string>
#include <vector>

using namespace std;

vector<string> split(const string& s, const string& sep);
void get_timediff_fom_gga(const string& gga,int64_t& diff);

#endif