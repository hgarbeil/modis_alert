#ifndef ALERTDB_H
#define ALERTDB_H 1
#endif
//#include "support/configCosmos.h"
//#include "math/mathlib.h"
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
using std::string;
#include <vector>
using std::vector;
//#include "plotfunc.h"

#define TEMPB21(x)   (3634.22/log(1.2246e11/(1e6*x)+1.))
#define RADB21(x) (1.2246e5/(exp(3634.22/(x))-1.))
#define ALERT_SAT_TERRA			1
#define ALERT_SAT_AQUA			2
const double DTOR=(M_PIl / (double)180.);
#define RADOF(deg)  (double)(DTOR * (deg))
#define RADOF(deg)  (double)(DTOR * (deg))

enum class ByteOrder : std::uint8_t {
    //! Big Endian byte order
    BIGENDIAN=0, // was previouly ORDER_BIGENDIAN, replace by ByteOrder::BIGENDIAN
    //! PowerPC byte order
    PPC=ByteOrder::BIGENDIAN,
    //! Motorola byte order
    MOTOROLA=ByteOrder::BIGENDIAN,
    //! Little Endian byte order
    LITTLEENDIAN=1, // was previouly ORDER_LITTLEENDIAN
    //! Intel byte order
    INTEL=ByteOrder::LITTLEENDIAN,
    //! Network byte order
    NETWORK=ByteOrder::BIGENDIAN
};


struct alert_entry
    {
    double date;
    double longitude;
    double latitude;
    uint16_t sat;
    float modis21;
    float modis22;
    float modis6;
    float modis31;
    float modis32;
    float ratio;
    uint16_t line;
    uint16_t sample;
    float zenith;
    float azimuth;
    float sun_zenith;
    float sun_azimuth;
    float glint_angle;
    int32_t next_entry;
    };

struct gst_handle
{
	FILE *mf;
	FILE *sf;
} ;

const size_t lst_size = 1200;
const double lst_radius = 6371007.181;
const double lst_step = M_PIl / 18.;
const size_t lst_max_quad = 36;

enum class lst_period : size_t
{
    AQUA_NIGHT,
    TERRA_DAY,
    AQUA_DAY,
    TERRA_NIGHT
};

struct lst_coord
{
    size_t vgrid;
    size_t hgrid;
    size_t vindex;
    size_t hindex;
    double lat;
    double lon;
    size_t month;
    lst_period period;
    double mean;
    double std;
};

struct lst_quad
{
    size_t vgrid;
    size_t hgrid;
    size_t month;
    lst_period period;
    double date;
//    double ullat;
//    double ullon;
//    double lrlat;
//    double lrlon;
    FILE *ff;
    vector<vector<double> > mean;
    vector<vector<double> > std;
//    double mean[lst_size][lst_size];
//    double std[lst_size][lst_size];
} ;

struct lst_handle
{
    size_t hgrid;
    size_t vgrid;
    size_t qindex;
    size_t qcount;
    vector<lst_quad> quads;
} ;

struct calstruc
{
    int32_t year;
    int32_t month;
    int32_t dom;
    int32_t doy;
    int32_t hour;
    int32_t minute;
    int32_t second;
    int32_t nsecond;
};



int32_t gst_read(alert_entry alert, double result[]);
int32_t gst_open();
int32_t gst_close();
int32_t lst_read(vector<lst_coord> &coords);
int32_t lst_open();
int32_t lst_load(lst_coord &coord);
int32_t lst_close();
int32_t lst_index(lst_coord &coord);
calstruc mjd2cal(double mjd);
int32_t mjd2ymd(double mjd, int32_t &year, int32_t &month, double &day, double &doy);
double currentmjd(double offset);
double currentmjd();
int16_t isleap(int32_t year);
double unix2utc(struct timeval unixtime);
float floatfrom(uint8_t *pointer, ByteOrder order);
ByteOrder local_byte_order();

//#ifdef __cplusplus
//}
//#endif

