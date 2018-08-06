/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sinuProjection.cpp
 * Author: harold
 * 
 * Created on August 3, 2018, 9:29 AM
 */

#include "sinuProjection.h"

sinuProjection::sinuProjection() {
    Radius = 6371007.181 ;
    cent_merid = 0. ;
}


void sinuProjection::set_parameters (float radius, float centmeridian) {
    Radius = radius ;
    cent_merid = centmeridian ;
}

// latlon in degrees, gets converted to radians in function 
void sinuProjection::latlon_to_xy(float lat, float lon, float* x, float* y) {
    *x = Radius * (lon-cent_merid) * DTOR * cos (lat * DTOR) ;
    *y = Radius * lat * DTOR ;
}

void sinuProjection::xy_to_latlon (float x, float y, float *lat, float *lon) {
    *lat = y / Radius /DTOR ;
    *lon = cent_merid +( x / (Radius * cos(lon * DTOR)))/DTOR ;
}

sinuProjection::sinuProjection(const sinuProjection& orig) {
}

sinuProjection::~sinuProjection() {
}
 

