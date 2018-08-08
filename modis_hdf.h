#ifndef mhdf
#define mhdf 

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

#include "mfhdf.h"

using namespace std;

class modis_hdf {
    char inhdfname [420];
    float *modbands_em, *modbands_refl;

    int32 sd_id, n_datasets, n_fileattrs;


public:
    bool thermFlag, aquaflag, dayflag;
    uint16 *thermdata, *refldata;
    float *latarr, *lonarr;
    int16 *solzen, *solaz, *senszen, *sensaz;
    float *th_scales_offsets, *refl_scales_offsets;
    float *raddata_cal, *refdata_cal, *solsens;
    modis_hdf(char *);
    ~modis_hdf();
    int open_hdf();
    void init_MOD21();
    void init_MOD03();
    void load_geometry();
    void load_thermal_bands();
    void calib_thermal_bands();
    void load_refSB_bands();
    void calib_refSB_bands();
    void calib_geometry();
    void get_date_period(char *infile, int *datearr);
    void get_file_name(char *, char *);



};




#endif

	

