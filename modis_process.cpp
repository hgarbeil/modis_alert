#include "modis_process.h"
#include <math.h>
#include "surftemp.h"



modis_process::modis_process () {

	distarr = 0L ;
	bandsarr = 0L ;
	temp = 0L ;
	nl = 2030 ;
	ns = 1354 ;
	gridspace = 0.1 ;


}

modis_process::~modis_process () {

	delete [] bandsarr ;
	delete [] distarr ;
	if (temp) delete [] temp ;
}

void modis_process::set_modis_hdfs (modis_hdf *g, modis_hdf *th) {
	geom = g ;
	therm = th ;
}

void modis_process::set_bounds  (float ulc_lat, float ulc_lon, float lrc_lat, 
        float lrc_lon, float gspace) {

		int i ;
		

		gridspace = gspace ; 

		nl = int (((ulc_lat - lrc_lat) / gridspace )+ 1) ;
		ns = int (((lrc_lon - ulc_lon) / gridspace )+ 1) ;
		startlat = ulc_lat ;
		startlon = ulc_lon ;
		endlat = lrc_lat ;
		endlon = lrc_lon ;

		cout << "Number of lines is : " << nl << endl ;
		cout << "Number of samples is : " << ns << endl ;

		if (distarr) {
			delete [] distarr ;
			delete [] bandsarr ;
		}
		distarr = new float [ns * nl] ;
		bandsarr = new float [8 * ns * nl] ;
		for (i=0; i<ns * nl; i++) {
			distarr [i] = 9999. ;
		}
		for (i=0; i<ns * nl * 8; i++) {
			bandsarr [i] = -1. ;
		}

}

void modis_process::set_month(int mon) {
	this->mon = mon ;
}

void modis_process::get_nearest_pixel (float lat, float lon, int *index, float *pixvals) {
    int i , minind;
    float plat, plon ;
    double xdist, ydist, dist, mperdegree, mindist ;
    
    mindist = 1.E9 ;
    mperdegree = 110000. ;
    
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
    int npix = 1354L * 2030 ;
    for (i=0; i<npix; i++) {
        ydist = (geom->latarr[i]-lat) * mperdegree ;
        xdist = (geom->lonarr[i]-lon) * mperdegree ;
        dist = sqrt (xdist * xdist + ydist * ydist) ;
        if (dist < mindist){
            mindist = dist ;
            minind = i ;
        }
    }
    *index = minind ;
    pixvals[0] = mindist ;
    // b21
    pixvals[1] = therm->raddata_cal[minind] ;
    pixvals[2] = this->bb_radtotemp (3.992, pixvals[1]) ;
    // b22
    pixvals[3] = therm->raddata_cal[npix+minind] ;
    // b31
    pixvals[4] = therm->raddata_cal[2*npix+minind] ;
    // solar zenith
    pixvals[5] = geom->solzen[minind]/100. ;
}

void modis_process::procalert () {
	
    float plat, plon, dist_pixel ;
	int i ;
	size_t period ;
	bool aflag, dflag ; 
	float b21val, b21temp ;
	vector<lst_coord> mpix ;
	lst_coord tcoord ;
	int npix = 2030 * 1354 ;
	temp = new float [npix] ;

	// get the period 
	aflag = therm->aquaflag ;
	dflag = therm->dayflag ;
	

	// if aqua, daytime is 2, nighttime is 1
	if (aflag) {
		if (dflag) 
			period = 2 ; 
		else 
			period = 0 ;
	}
	// if terra, daytime is 1, nighttime is 3
	else {
		if (dflag) 
			period = 1 ;
		else period = 3 ;
	}

	tcoord.period = lst_period(period) ;

	

    for (i=0; i<npix; i++) {
		tcoord.lat = RADOF(geom->latarr[i]) ;
		tcoord.lon = RADOF(geom->lonarr[i]) ;
        plat = geom->latarr[i] ;
        plon = geom->lonarr[i] ;
		tcoord.month = mon ;
		mpix.push_back (tcoord) ;
		
		b21val = therm->raddata_cal [i] ;
		// convert the b21 radiance to temperature
		b21temp = this->bb_radtotemp (3.992, b21val) ;
        temp[i] = b21temp ;
        

	}
  

	// now get the mean and std dev
	lst_open() ;
	lst_read(mpix) ;
	float zval ;	
	for (i=0; i<npix; i++) {
		tcoord = mpix.at(i) ;
		if (tcoord.mean <200) {
			continue ;
		}
		zval = (temp[i]-tcoord.mean) / tcoord.std ;
		if (zval>1) {
			//cout << "Mean temp is : " << tcoord.mean <<  "  "  << tcoord.std<< "  " << temp[i] << endl ;
		}
	}


	//delete [] temp ;	
		
			
}

float modis_process::calcdist (float plat, float plon){
    
	double xdist, ydist, tdist ;
	float reflat, reflon ;
    reflat = 19.33 ;
    reflon = -155.33 ;
    xdist = abs(reflat - plat) * 100000. ;
	ydist = abs(reflon - plon) * 100000. ;
	tdist = sqrt (xdist * xdist + ydist * ydist) ;
	return (float(tdist)) ;
}


void modis_process::process () {

	int i, j, ib, snum, is, js, ixloc, iyloc, grid_x, grid_y, gridloc,npix ;
	float xloc, yloc, xdist, ydist, dist, latval, lonval ;
	npix = nl * ns ;
	

	for (i=0; i<2030; i++) {
		for (j=0; j<1354; j++) {
			snum = i*1354+j ;
			latval = geom->latarr[snum] ;
			lonval = geom->lonarr[snum] ;
			if (latval > startlat || latval < endlat) 
				continue ;
			if (lonval < startlon || lonval > endlon) 
				continue ;

			xloc = (lonval - startlon) / gridspace ;
			yloc = (startlat - latval) / gridspace ;
			ixloc = int(xloc+0.5) ;
			iyloc = int(yloc+0.5) ;

			for (is=-2; is<=2; is++) {
				
				grid_y = iyloc + is ;
				if (grid_y < 0) continue ;
				if (grid_y >= nl) continue ;

				ydist = yloc - grid_y ;
				for (js=-2; js<=2; js++) {
					grid_x = is + ixloc ;
					if (grid_x < 0) continue ;
					if (grid_x >= ns) continue ;
					gridloc = grid_y * ns + grid_x ;
					xdist = xloc - grid_x ;
					dist = sqrt (xdist * xdist + ydist * ydist) ;
					if (dist < distarr[grid_y * ns + grid_x])
					{
						distarr[iyloc*ns+ixloc]= dist ;
						bandsarr[gridloc] = therm->refdata_cal[i*ns+j] ;
						for (ib=0; ib<3; ib++) 
							bandsarr[(ib+1) * npix + gridloc] = therm->raddata_cal[ib*2030*1354+snum] ;
						bandsarr[4*npix+gridloc] = geom->solsens[snum] ; 
						bandsarr[5*npix+gridloc] = geom->solsens[2030L*1354+snum] ;
						bandsarr[6*npix+gridloc] = geom->solsens[2*2030L*1354+snum] ;
						bandsarr[7*npix+gridloc] = geom->solsens[3*2030L*1354+snum] ;
					}
				}
			}
		}
	}
}


void modis_process::calc_alert (char *outfile) {
	char hdrfile [420] ;
	int i, ib, npix ;
	float val21,val22,val32, val6 ;

	npix = 2030 * 1354L ;
	float *alert = new float [npix] ;
	
	if (therm->dayflag) 
	for (i=0; i<npix; i++) {
		// band 6,21,22, 32
		val6  = therm->refdata_cal[i] ;	
		val21  = therm->raddata_cal[i] ;	
		val22  = therm->raddata_cal[npix+i] ;	
		val32  = therm->raddata_cal[2*npix+i] ;	
		if (val21 > val22) 
			val22 = val21 ;
		val22 = val22 - .0426 * val6 ;
		*(alert+i) = (val22 - val32) / (val22 + val32) ;
	}
	else 
	for (i=0; i<npix; i++) {
		// band 6,21,22, 32
		val6  = therm->refdata_cal[i] ;	
		val21  = therm->raddata_cal[i] ;	
		val22  = therm->raddata_cal[npix+i] ;	
		val32  = therm->raddata_cal[2*npix+i] ;	
		// check if val22 is saturated, then use val21
		if (val21 > val22) 
			val22 = val21 ;
		*(alert+i) = (val22 - val32) / (val22 + val32) ;
	}

	FILE *fout = fopen (outfile, "w") ;
	fwrite ((char *)alert, 4, npix, fout) ;
	fclose (fout) ;

	strcpy (hdrfile, outfile) ;
	strcat (hdrfile, ".hdr") ;
	write_header (hdrfile, 1) ;

	delete [] alert ;

}

void modis_process::write_output (char *outfile) {
	char hdrfile [420] ;
	int i, npix ;
	npix = ns * nl ;
	cout << "Number of samples is : " << ns << endl ;
	cout << "Number of lines is : " << nl << endl ;

	strcpy (hdrfile, outfile) ;
	strcat (hdrfile, ".hdr") ;
	write_header (hdrfile, 6) ;

	FILE *fout = fopen (outfile, "w") ;
	if (fout == NULL) {
		cout << "Could not open " << outfile << endl ;
		return ;
	}
	for (i=0; i<3; i++) {
		fwrite ((char *) &therm->raddata_cal[i * npix], 4, npix, fout) ;
	}
	fwrite ((char *) &temp[0], 4, npix, fout) ;
	fwrite ((char *) geom->latarr,4,npix,fout) ;
	fwrite ((char *) geom->lonarr,4,npix,fout) ;
	fclose (fout) ;

}


void modis_process::write_header (char *outfile, int nbands)
{

	FILE *hdrout = fopen (outfile, "w") ;
	char bnames [1200] ;
	strcpy (bnames, "band names = {\nModisBand6,ModisBand21,ModisBand22,ModisBand32,\nSolarAzimuth,SolarZenith,SensorAzimuth,SensorZenith}\n") ;

    fprintf (hdrout, "ENVI\ndescription = {\nMOD021KM - MOD03  }\n") ;
    fprintf (hdrout, "samples    = %5d\n",ns) ;
	fprintf (hdrout, "lines      = %5d\n",nl) ;
	fprintf (hdrout, "bands      = nbands\n") ;
	fprintf (hdrout, bnames) ;
	fprintf (hdrout, "header offset = 0 \n") ;
	fprintf (hdrout, "file type = ENVI Standard \n") ;
	fprintf (hdrout, "data type = 4 \n") ;
	fprintf (hdrout, "interleave = bsq \n") ;
	fprintf (hdrout, "sensor type = MODIS \n") ;
	fprintf (hdrout, "byte order = 0\n") ;
	fprintf (hdrout, "map info = {Geographic Lat/Lon, 1.0000, 1.0000, %12.8f, %12.8f,",
	         startlon, startlat) ;
	fprintf (hdrout, "%6.4f, %6.4f, WGS-84, units=Degrees}\n", gridspace, gridspace) ;

	fclose (hdrout) ;
}

float modis_process::bb_radtotemp (float wave, float rad) {

	float wave_m, l_m ;
	double h, k, c, c1, c2, val0, Temp ;

	// convert to m
	l_m = rad * 1.E6  ;
	wave_m = wave * 1.E-6 ;
	h = 6.62606755E-34 ;
	k = 1.380658E-23 ;
	c = 2.9979246E8 ;
	c1 = 2. * h * c * c ;
	c2 = h * c / k ;

	val0 = c2 / (wave_m * log(c1/ (l_m * pow(wave_m,5))+1.)) ;

	return float(val0)  ;
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



