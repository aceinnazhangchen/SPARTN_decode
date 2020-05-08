#include "stringex.h"
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>

time_t mkgmtime_g(struct tm *unixdate)
{
    //assert(unixdate != nullptr);

    time_t     fakeUnixTime = mktime(unixdate); 
    struct tm *fakeDate     = gmtime(&fakeUnixTime);

    int32_t nOffset = fakeDate->tm_hour - unixdate->tm_hour;
    if (nOffset > 12)
    {
        nOffset = 24 - nOffset;
    }
    return fakeUnixTime - nOffset * 3600;
}

vector<string> split(const string& s, const string& sep)
{
    vector<string> v; 
    string::size_type pos1, pos2;
    pos2 = s.find(sep);
    pos1 = 0;
    while(string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2-pos1));
         
        pos1 = pos2 + sep.size();
        pos2 = s.find(sep, pos1);
    }
    if(pos1 != s.length())
        v.push_back(s.substr(pos1));
    return v;
}

 void get_timediff_fom_gga(const string& gga,int64_t& diff)
 {
    std::vector<std::string> vgga = split(gga,",");
    if(vgga.size() == 15)
    {
        std::string str_gga_time = vgga[1];
        if(str_gga_time.size() == 9)
        {
            int _hour = std::stoi(str_gga_time.substr(0,2));
            int _min = std::stoi(str_gga_time.substr(2,2));
            int _sec = std::stoi(str_gga_time.substr(4,2));
            int _milli = std::stoi(str_gga_time.substr(7,2));

            struct timeb tb;
            ftime(&tb);

            struct tm* local_time = localtime(&tb.time);
            time_t n_local_time = mktime(local_time);
            struct tm* gga_time = gmtime(&tb.time);
            time_t n_gm_time = mktime(gga_time);

            gga_time->tm_hour = _hour;
            gga_time->tm_min = _min;
            gga_time->tm_sec = _sec;
            

            time_t gga_t = mktime(gga_time);
            gga_t += n_local_time - n_gm_time;
            int64_t milli_gga_time = gga_t*1000 + _milli;
            //char sz_gga_time[128] = {0};
            //std::sprintf(sz_gga_time,"%d.%d",gga_t,_milli);
            //double double_gga_time = atof(sz_gga_time);

            int64_t milli_current_time = tb.time*1000 + tb.millitm;
            //char sz_current_time[128] = {0};
            //std::sprintf(sz_current_time,"%d.%d",t,tb.millitm);
            //double double_current_time = atof(sz_current_time);

            diff = milli_current_time - milli_gga_time;  
        }
    }
 }
