#include "modis_process.h"
#include <math.h>
#include "surftemp.h"
#include "sinuProjection.h"

modis_process::modis_process() {

    distarr = 0L;
    bandsarr = 0L;
    temp = 0L;
    alerts = 0L;
    nl = 2030;
    ns = 1354;
    gridspace = 0.1;
    day_limit = -.6;
    night_limit = -.8;
    alerts = new float [ns * nl];
    sp = 0L ;
    glob_mnstdev = 0L ;
    basevals_full = 0l ;
    alertFlag = false ;


}

modis_process::~modis_process() {

    delete [] bandsarr;
    delete [] distarr;
    delete [] alerts;
    if (temp) delete [] temp;
    if (basevals_full) delete[] basevals_full ;
    delete [] glob_mnstdev ;
}

void modis_process::set_modis_hdfs(modis_hdf *g, modis_hdf *th) {
    geom = g;
    therm = th;
    alertFlag = false ;
}

void modis_process::set_bounds(float ulc_lat, float ulc_lon, float lrc_lat,
        float lrc_lon, float gspace) {

    int i, npix_grid ;


    gridspace = gspace;

    ny_grid = int (((ulc_lat - lrc_lat) / gridspace) + 1);
    nx_grid = int (((lrc_lon - ulc_lon) / gridspace) + 1);
    npix_grid = nx_grid * ny_grid ;
    startlat = ulc_lat;
    startlon = ulc_lon;
    endlat = lrc_lat;
    endlon = lrc_lon;

    cout << "Number of lines is : " << ny_grid << endl;
    cout << "Number of samples is : " << nx_grid << endl;

    if (distarr) {
        delete [] distarr;
        delete [] bandsarr;
        delete [] glob_mnstdev ;
    }
    distarr = new float [npix_grid];
    bandsarr = new float [7 * npix_grid];
    glob_mnstdev = new float [6 * npix_grid] ;
    // alloc memory for the baseline data b221, b32, may do all b6 and ndti
    
    for (i = 0; i < npix_grid; i++) {
        distarr [i] = 9999.;
    }
    for (i = 0; i < npix_grid * 7; i++) {
        bandsarr [i] = -1.;
    }

}


void modis_process::set_month(int mon) {
    this->mon = mon;
}

void modis_process::get_nearest_pixel(float lat, float lon, int *index, float *pixvals) {
    int i, minind;
    float plat, plon;
    double xdist, ydist, dist, mperdegree, mindist;

    mindist = 1.E9;
    mperdegree = 110000.;

    // go through the data set lat and lon arrays :finding pixel closest to the ref pixel
    // defined by the function arguments (lat, lon)
    // the function will load up the index pointer with the closest pixel index and the
    // pixvals will be loaded with
    // [0] = mindist
    // [1] = b21 radiance
    // [2] = b22 b temp
    // [3] = b22 radiance
    // [4] = b32 radiance
    // [5] = solar zenith
    int npix = 1354L * 2030;
    for (i = 0; i < npix; i++) {
        ydist = (geom->latarr[i] - lat) * mperdegree;
        xdist = (geom->lonarr[i] - lon) * mperdegree;
        dist = sqrt(xdist * xdist + ydist * ydist);
        if (dist < mindist) {
            mindist = dist;
            minind = i;
        }
    }
    *index = minind;
    pixvals[0] = mindist;
    // b21
    pixvals[1] = therm->raddata_cal[minind];
    pixvals[2] = this->bb_radtotemp(3.992, pixvals[1]);
    // b22
    pixvals[3] = therm->raddata_cal[npix + minind];
    // b31
    pixvals[4] = therm->raddata_cal[2 * npix + minind];
    // solar zenith
    pixvals[5] = geom->solzen[minind] / 100.;
}
/*
void modis_process::procalert() {

    float plat, plon, dist_pixel;
    int i;
    size_t period;
    bool aflag, dflag;
    float b21val, b21temp;
    vector<lst_coord> mpix;
    lst_coord tcoord;
    int npix = 2030 * 1354;
    temp = new float [npix];

    // get the period 
    aflag = therm->aquaflag;
    dflag = therm->dayflag;


    // if aqua, daytime is 2, nighttime is 1
    if (aflag) {
        if (dflag)
            period = 2;
        else
            period = 0;
    }        // if terra, daytime is 1, nighttime is 3
    else {
        if (dflag)
            period = 1;
        else period = 3;
    }

    tcoord.period = lst_period(period);



    for (i = 0; i < npix; i++) {
        tcoord.lat = RADOF(geom->latarr[i]);
        tcoord.lon = RADOF(geom->lonarr[i]);
        plat = geom->latarr[i];
        plon = geom->lonarr[i];
        tcoord.month = mon;
        mpix.push_back(tcoord);

        b21val = therm->raddata_cal [i];
        // convert the b21 radiance to temperature
        b21temp = this->bb_radtotemp(3.992, b21val);
        temp[i] = b21temp;


    }


    // now get the mean and std dev
    lst_open();
    lst_read(mpix);
    float zval;
    for (i = 0; i < npix; i++) {
        tcoord = mpix.at(i);
        if (tcoord.mean < 200) {
            continue;
        }
        zval = (temp[i] - tcoord.mean) / tcoord.std;
        if (zval > 1) {
            //cout << "Mean temp is : " << tcoord.mean <<  "  "  << tcoord.std<< "  " << temp[i] << endl ;
        }
    }


    //delete [] temp ;	


}
*/
float modis_process::calcdist(float plat, float plon) {

    double xdist, ydist, tdist;
    float reflat, reflon;
    reflat = 19.33;
    reflon = -155.33;
    xdist = abs(reflat - plat) * 100000.;
    ydist = abs(reflon - plon) * 100000.;
    tdist = sqrt(xdist * xdist + ydist * ydist);
    return (float(tdist));
}

void modis_process::process() {

    int i, j, ib, snum, is, js, ixloc, iyloc, grid_x, grid_y, gridloc, npix_grid;
    float xloc, yloc, xdist, ydist, dist, latval, lonval;
    npix_grid = ny_grid * nx_grid;
    ;

    //if (!alertFlag) this->calc_alert() ;

    for (i = 0; i < 2030; i++) {
        for (j = 0; j < 1354; j++) {
             
            snum = i * 1354 + j;
            latval = geom->latarr[snum];
            lonval = geom->lonarr[snum];
            if (latval > startlat || latval < endlat)
                continue;
            if (lonval < startlon || lonval > endlon)
                continue;

           
            
            xloc = (lonval - startlon) / gridspace;
            yloc = (startlat - latval) / gridspace;
            ixloc = int(xloc + 0.5);
            iyloc = int(yloc + 0.5);

            for (is = -2; is <= 2; is++) {

                grid_y = iyloc + is;
                if (grid_y < 0) continue;
                if (grid_y >= ny_grid) continue;

                ydist = yloc - grid_y;
                for (js = -2; js <= 2; js++) {
                    grid_x = js + ixloc;
                    if (grid_x < 0) continue;
                    if (grid_x >= nx_grid) continue;
                    gridloc = grid_y * nx_grid + grid_x;
                    xdist = xloc - grid_x;
                    dist = sqrt(xdist * xdist + ydist * ydist);
                    if (dist < distarr[gridloc]) {
                        distarr[gridloc] = dist;
                        bandsarr[gridloc] = therm->refdata_cal[snum];
                        // resample the three bands of radiance to bandsarr
                        
                        for (ib = 0; ib < 3; ib++)
                            bandsarr[(ib + 1L) * npix_grid + gridloc] = therm->raddata_cal[ib * 2030L * 1354 + snum];
                        bandsarr[4*npix_grid+gridloc] = *(alerts+snum) ;

                        //bandsarr[4L * npix_grid + gridloc] = geom->solsens[snum];
                        //bandsarr[5L * npix_grid + gridloc] = geom->solsens[2030L * 1354 + snum];
                        //bandsarr[6L * npix_grid + gridloc] = geom->solsens[2 * 2030L * 1354 + snum];
                        //bandsarr[7L * npix_grid + gridloc] = geom->solsens[3 * 2030L * 1354 + snum];

                    }
                }
            }
        }
    }
}

void modis_process::calc_alert(char *alertfile) {
    char hdrfile [420], outstr[240];
    int i, ind, ib, npix, num_alerts_full;
    float zval, val21, val22, val32, val6, alval, alert_limit ;

    npix = 2030 * 1354L;
    alinds.clear();
    alertFlag = true ;

    alert_limit = night_limit ;
    if (geom->dayflag) {
        alert_limit = day_limit ;
    }
    FILE *fout = fopen (alertfile, "w") ;
    for (i = 0; i < npix; i++) {
        // band 6,21,22, 32
        val6 = therm->refdata_cal[i];
        val21 = therm->raddata_cal[i];
        val22 = therm->raddata_cal[npix + i];
        val32 = therm->raddata_cal[2 * npix + i];
        // check if val22 is saturated, then use val21
        if (val22 > 2.0 || val22 < 0.001)
            val22 = val21;
        // also check to make sure 21 and 32 are good 
        if (val21 > 105. || val32 > 40.) {
            *(alerts+i) = -2. ;
            continue ;
        }
        if (val32 > 40.) continue ;
        alval = (val22 - val32) / (val22 + val32);
        if (alval > -0.05) {
            int iline = i / 1354 ;
            int ipix = i - iline * 1354 ;
            cout << "crazy value hit -- line samp  : "  << iline << "  " << ipix << "  " << endl ;
            
        }
        if (alval > alert_limit)
            alinds.push_back(i);
        *(alerts + i) = alval;
       
    }
    num_alerts_full = alinds.size() ;
    if (basevals_full) delete[] basevals_full ;
    basevals_full = new float [6 * num_alerts_full] ;
    cout <<"Number of alerts is modis scene is "  << num_alerts_full << endl ;
    sp->getCorrespondingValues (alinds, mon, geom->aq_terra_flag, geom->latarr, geom->lonarr, basevals_full) ;
    for (i=0; i< num_alerts_full; i++) {
        
        ind = alinds.at(i) ;
        zval = (alerts[ind] - basevals_full[6*i+4])/basevals_full[6*i+5] ;
            sprintf (outstr, "%d\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%3.1f\t%7.3f\t%7.3f\r\n", 
                    ind, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
                    alerts[ind], zval, basevals_full[6*i+4], basevals_full[6*i+5]) ;
        fputs (outstr, fout)  ;  
    }
    fclose (fout) ;
}
      
// write out each alert to an appended file with the following columns
// mod21 filename
// year
// month
// day
// sol zenith
// b21
// b22
// b32
// alert value

void modis_process::write_alert_textfile(char *outfile) {
    FILE *fout = fopen(outfile, "a+");
    int nalerts = alinds.size();

}

void modis_process::write_output(char *outfile) {
    char hdrfile [420];
    int i, npix;
    npix = nx_grid * ny_grid;
    cout << "Number of samples is : " << nx_grid << endl;
    cout << "Number of lines is : " << ny_grid << endl;

    strcpy(hdrfile, outfile);
    strcat(hdrfile, ".hdr");
    write_header(hdrfile, 13);
    cout << "Output file is : " << outfile << endl ;

    FILE *fout = fopen(outfile, "w");
    if (fout == NULL) {
        cout << "Could not open " << outfile << endl;
        return;
    }
    
    
    fwrite((char *) bandsarr, 4, 7 * npix, fout);
    
    
    // band21,22,32,alert, mean stdev
    fwrite((char *) this->glob_mnstdev, 4, 6 * npix, fout) ;
    
    
   // fwrite ((char *) this->bandsarr, 4, 8* npix, fout) ;
    
    fclose(fout);

}

void modis_process::write_header(char *outfile, int nbands) {

    FILE *hdrout = fopen(outfile, "w");
    char bnames [1200];
    strcpy(bnames, "band names = {\nModisBand6,ModisBand21,ModisBand22,ModisBand32,Alert,Alert_z, Alert_rat,10yr_2122,10yr_2221_stdv,10yr_32,10yr_32_stdv, 10yr_nti, 10yr_nti_stdv}\n");

    fprintf(hdrout, "ENVI\ndescription = {\nMOD021KM - MOD03  }\n");
    fprintf(hdrout, "samples    = %5d\n", nx_grid);
    fprintf(hdrout, "lines      = %5d\n", ny_grid);
    fprintf(hdrout, "bands      = %3d\n", nbands);
    fprintf(hdrout, bnames);
    fprintf(hdrout, "header offset = 0 \n");
    fprintf(hdrout, "file type = ENVI Standard \n");
    fprintf(hdrout, "data type = 4 \n");
    fprintf(hdrout, "interleave = bsq \n");
    fprintf(hdrout, "sensor type = MODIS \n");
    fprintf(hdrout, "byte order = 0\n");
    fprintf(hdrout, "map info = {Geographic Lat/Lon, 1.0000, 1.0000, %12.8f, %12.8f,",
            startlon, startlat);
    fprintf(hdrout, "%6.4f, %6.4f, WGS-84, units=Degrees}\n", gridspace, gridspace);

    fclose(hdrout);
}

float modis_process::bb_radtotemp(float wave, float rad) {

    float wave_m, l_m;
    double h, k, c, c1, c2, val0, Temp;

    // convert to m
    l_m = rad * 1.E6;
    wave_m = wave * 1.E-6;
    h = 6.62606755E-34;
    k = 1.380658E-23;
    c = 2.9979246E8;
    c1 = 2. * h * c * c;
    c2 = h * c / k;

    val0 = c2 / (wave_m * log(c1 / (l_m * pow(wave_m, 5)) + 1.));

    return float(val0);
}

void modis_process::zscore () {
    int i, npix ;
    npix = nx_grid * ny_grid ;
    float alval ;
    for (i=0; i<npix; i++) {
        
        alval = (bandsarr[4*npix+i]-glob_mnstdev[4*npix+i]) / glob_mnstdev[5*npix+i] ;
        bandsarr[5*npix+i] = alval ;
        bandsarr[6*npix+i] =   glob_mnstdev[4*npix+i] / bandsarr[4*npix+i];
    }
    
    
}


void modis_process::extract_from_baseline_file () {
    int i ;
    sp = new sinuProjection () ; 
    
    
    //sp->set_projection_parameters ()
    sp->fillGrid (this->mon, this->geom->aq_terra_flag, startlat, startlon, gridspace, nx_grid, ny_grid, glob_mnstdev) ;
    //for (i=0; i<nx_grid*ny_grid; i++) {
    
}



/*
https://cimss.ssec.wisc.edu/dbs/China2011/Day2/Lectures/MODIS_DB_Land_Surface_Temperature_reference.pdf
float modis_process::bb_radtotemp (float wave, float rad) {

        float wave_m, l_m ;
        double h, k, c, K1, K2, val0, Temp ;

        // wave in microns 
        l_m = rad   ;
        wave_m = wave  ;
        K2 = 1.43883E4 / wave;
        K1 = 1.19107E8 ;
        val0 = K1 / (l_m*pow(wave,5)) + 1. ;
        Temp = (K2 / log (val0))  ;
        return float(Temp) ;
}


 */
