#ifndef m_rs 
#define m_rs

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "modis_hdf.h"
#include "sinuProjection.h"

using namespace std;

class modis_process {
    bool alertFlag ;

	int nsamps, nlines, ns, nl, nx_grid, ny_grid ;
	float startlat, startlon, endlat, endlon, gspace ;
	float *distarr, *bandsarr, *temp, *alerts, *glob_mnstdev, *basedata, *basevals_full ;
	float day_limit, night_limit ;
	modis_hdf *geom, *therm ;
        
	// alert indices
	vector <int> alinds ;


public:
    float gridspace;
    modis_process();
    ~modis_process();
    sinuProjection *sp ;
    void write_output(char *);
    void set_month(int m);
    void set_modis_hdfs(modis_hdf *geom, modis_hdf *therm);
    void set_bounds(float ulc_lat, float ulc_lon, float llc_lat,
            float llc_lon, float gridspace);
    float bb_radtotemp(float, float);
    void process();
    //void procalert();
    void write_header(char *, int);
    void write_alert_textfile(char *);
    void extract_from_baseline_file();
    float calcdist(float, float);
    void calc_alert(char *);
    void zscore () ;
    void get_nearest_pixel(float lat, float lon, int *ind, float *pargs);
    int num_alerts_full, mon ;
                        


} ;

#endif
