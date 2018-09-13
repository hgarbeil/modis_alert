/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sinuProjection.cpp
 * Purpose : Various functions for the sinuProjection class. This will relate the latitude and longitude 
 * based to an xy cartesian system based upon the sinusoidal projection. 
 * Author: harold garbeil
 * 
 * Created on August 3, 2018, 9:29 AM
 */

#include "sinuProjection.h"
#include "surftemp.h"

sinuProjection::sinuProjection() {
    Radius = 6371007.181 ;
    cent_merid = 0. ;
    mDegToRad = 0.01745329 ;
    wrld_starty =  10007554.677899 ;
    wrld_startx = -20015109.35496 ;
    gspace = 4633.1271655 ;
}


/**
 * allows the user to set the radius and central meridian if different than the default.
 * @param radius spherical radius of the planet - default is 6371007.181 m
 * @param centmeridian central meridian specifying the x zero location, long to 
 * the west is -ve, to the east is +ve -default is zero
 */
void sinuProjection::set_parameters (float radius, float centmeridian) {
    Radius = radius ;
    cent_merid = centmeridian ;
}

// latlon in degrees, gets converted to radians in function 
void sinuProjection::latlon_to_xy(float lat, float lon, float* x, float* y) {
    double latd = lat ;
    double lond = lon ;
    *x = Radius * (lond-cent_merid) * mDegToRad * cos (latd * mDegToRad) ;
    *y = Radius * latd * mDegToRad ;
}

void sinuProjection::xy_to_latlon (float x, float y, float *lat, float *lon) {
    *lat = y / Radius /mDegToRad ;
    *lon = cent_merid +( x / (Radius * cos(*lat * mDegToRad)))/mDegToRad ;
}

sinuProjection::sinuProjection(const sinuProjection& orig) {
}

sinuProjection::~sinuProjection() {
}
 
/**
 * 
 * @param mon       int - month to fill the grid with, should add something for which satellite, etc
 * @param aq_terra_flag int - aqua (0-night 2-day) terra (3-night 1-day)
 * @param ulc_lat   float - northernmost latitude
 * @param ulc_lon   float - westernmost longitude
 * @param gridspace float - grid spacing in decimal degrees 
 * @param nx        int - number of samples in output grid
 * @param ny        int - number of lines in output grid
 * @param outarr    float * - 8 banded nx * ny float array, memory allocated before calling function
 */
void sinuProjection::fillGrid (int mon, int aq_terra_flag, float ulc_lat, float ulc_lon, float gridspace, int nx, int ny, float *outarr) {
    
    int i, j, ixloc, iyloc, itype, npix_grid, iretn ;
    float lat, lon, xloc, yloc ;
    char bfile [420] ;
    m02ssh_coordvalue lst_val ;
    vector <m02ssh_coordvalue> mpix ;
    
    npix_grid = nx * ny ;
    
    // the ifil will depend upon the band, and the month to be extracted from the ave,stdev sinusoidal files.
    char const *outstr[] = {"rad221", "rad32", "nti"} ;
    int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;

    

    for (itype = 0; itype < 3; itype++) {
        sprintf(basef, "%s_%03d_%02d.bsq", outstr[itype], mon_days[mon], aq_terra_flag) ;
        sprintf(bfile, "/local/worldbase/02ssh/%s_%03d_%02d.bsq\0", outstr[itype], mon_days[mon], aq_terra_flag);
        iretn = m02ssh_open (outstr[itype],mon_days[mon], modis_period(aq_terra_flag)) ;
        if (iretn <0) {
            cout << "Problem opening ... " << bfile << endl ;
        }
        cout << "opening baseline file of : " <<bfile << endl ;
        //FILE *fil = fopen(bfile, "r");
        //fread(aldat, 8, 4320L * 8640, fil);
        //fclose(fil);
        lst_val.jday = mon_days[mon] ;
        for (i = 0; i < ny; i++) {
            lat = ulc_lat - gridspace * i ;
            mpix.clear() ;
            lst_val.lat = RADOF (lat) ;
            
            for (j = 0; j < nx; j++) {
                lon = ulc_lon + gridspace * j  ;
                lst_val.lon = RADOF (lon) ;
                mpix.push_back (lst_val) ;
                
                //outarr [2 * itype * npix_grid + npix_grid + i * nx + j] = aldat [8640L * 4320 + iyloc * 8640L + ixloc];
            }
            m02ssh_read (mpix) ;
            for (j=0; j<nx; j++) {
                outarr [2 * itype * npix_grid + i * nx + j] = mpix[j].mean ;
                outarr [2 * itype * npix_grid + npix_grid + i * nx + j] = mpix[j].std ;
            }
        }
    }
    /*
    FILE *fout = fopen ("/home/harold/aldat", "w") ;
    fwrite (outarr, 8, nx * ny, fout) ;
    fclose (fout) ;
    */
    
}



void sinuProjection::getCorrespondingValues(vector<int> alinds, int mon, int aq_terra_flag,  float *lat, float *lon, float *vals) {
    char bfile [240] ;
    int i, ind, itype, num_alerts, iretn ;
    char const *outstr[] = {"rad221" , "rad32", "nti"} ;
    int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;
    int ixloc, iyloc ;
    float xloc, yloc ;
    m02ssh_coordvalue lst_val ;
    vector <m02ssh_coordvalue> mpix ;
    
    num_alerts = alinds.size() ;
    
    if (num_alerts <= 0) return ;
    for (itype = 0; itype < 3; itype++) {
        mpix.clear() ;
        sprintf(basef, "%s_%03d_%02d.bsq", outstr[itype], mon_days[mon], aq_terra_flag) ;
        sprintf(bfile, "/local/worldbase/02ssh/%s_%03d_%02d.bsq\0", outstr[itype], mon_days[mon], aq_terra_flag);
        iretn = m02ssh_open (outstr[itype],mon_days[mon], modis_period(aq_terra_flag)) ;
        if (iretn < 0)
            cout << "problem opening baseline file of : " <<bfile << endl ;
        
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period (aq_terra_flag) ;
        for (i = 0; i < num_alerts; i++) {
            ind = alinds.at(i);
            this->latlon_to_xy(lat[ind], lon[ind], &xloc, &yloc);
            lst_val.lat = RADOF (lat[ind]) ;
            lst_val.lon = RADOF (lon[ind]) ;
            mpix.push_back (lst_val) ;
            
            

        }
        m02ssh_read (mpix) ;
        for (i = 0; i < num_alerts; i++) {
            vals [i*6+2 * itype] = mpix[i].mean ;
            vals [i*6+2 * itype+1] = mpix[i].std ;
            
        }
        
    }
    
    
}