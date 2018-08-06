/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sinuProjection.h
 * Author: harold
 *
 * Created on August 3, 2018, 9:29 AM
 */

#ifndef SINUPROJECTION_H
#define SINUPROJECTION_H

#define DTOR  0.017453

class sinuProjection {
public:
    sinuProjection();
    sinuProjection(const sinuProjection& orig);
    virtual ~sinuProjection();
    void set_parameters (float radius, float centmeridian)  ;
    void latlon_to_xy (float lat, float lon, float *x, float *y) ;
    void xy_to_latlon (float x, float y, float *lat, float *lon) ;
    
private:
    float Radius ;
    float cent_merid ;
};

#endif /* SINUPROJECTION_H */

