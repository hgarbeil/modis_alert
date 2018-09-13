#include "modis_hdf.h"
#include "time.h"

modis_hdf::modis_hdf (char *infile) {

	aquaflag = false ;
	thermFlag = false ;
	dayflag = false ;
        aq_terra_flag = 0 ;
	if (strstr (infile, "MYD")) {
		aquaflag = true ;
	}
	// check for MOD
	if (strstr(infile, "021KM")) {
		thermFlag = true ;
	}
	 
	strcpy (inhdfname, infile) ;
	open_hdf () ;
	if (thermFlag) 
		init_MOD21 () ;
	else 
		init_MOD03 () ;
	

}

void modis_hdf::get_date_period (char *infile, int *datearr) {
	
	//int imon, days_in_month[] ={31,28,31,30,31,30,31,31,30,31,30,31} ;
        int imon, days_in_month[] = {32,30,30,31,30,31,30,31,30,30,31,30} ;

	char cfile [420], *tmp ;
	int year, iyear, yday, yday1, mday, mon, hr ;
	strcpy (cfile, infile) ;

	tmp = strtok (cfile, ".") ;
	tmp = strtok (NULL, ".") ;
	year = atoi (&tmp[1]) ;
	iyear = year / 1000 ;
	if (iyear %4 ==0) 
		days_in_month[1]++ ;
	yday = year - (iyear *1000) ;
	yday1 = yday ;

	mday = yday ;
	for (imon=0; imon<12; imon++) {
		if (yday1 <= days_in_month[imon]) 
			break ;
		yday1 -= days_in_month[imon] ;
		mday = yday1 ;
	}
	mon = imon + 1 ;
	printf ("Year is %d  day is %d\r\n", iyear, mday ) ;
	tmp = strtok (NULL, ".") ;
	hr = atoi (tmp) /100 ;

	datearr[0] = iyear ;
	datearr[1] = mon ;
	datearr[2] = mday ;
	datearr[3] = hr ;


	


}

void modis_hdf::init_MOD21() {
    thermdata = new uint16 [16 * 1354 * 2030L];
    th_scales_offsets = new float [2 * 16];
    refldata = new uint16 [5 * 1354 * 2030];
    refl_scales_offsets = new float [2 * 5];
    raddata_cal = new float [3 * 1354 * 2030L];
    refdata_cal = new float [1 * 1354 * 2030L];
    load_refSB_bands();
    calib_refSB_bands();
    load_thermal_bands(); //calib_thermal_bands() ;
    calib_thermal_bands();


}

// this function takes the full path of the hdf file (ifile) and will output only the 
// filename, this is used for output file naming  and writing to a pre-specified directory
void modis_hdf::get_file_name (char *ifile, char *outfname) {
    char workfil[420], fname[420], *tmp ;
    
    strcpy (workfil, ifile) ;
    tmp = strtok (workfil, "\\/") ;
    // search for the last file separator
    while (1){
        tmp = strtok (NULL, "\\/") ;
        if (tmp == NULL) break ;
        strcpy (fname, tmp) ;
    }
    strcpy (outfname, fname) ;
    
}


void modis_hdf::init_MOD03 () {
	latarr = new float [1354 * 2030] ;
	lonarr = new float [1354 * 2030] ;
	solaz = new int16 [1354 * 2030] ;
	solzen = new int16 [1354 * 2030] ;
	sensaz = new int16 [1354 * 2030] ;
	senszen = new int16 [1354 * 2030] ;
	solsens = new float [4 * 1354 * 2030] ;
	load_geometry () ;
	calib_geometry () ;
        get_aq_terra_flag () ;

}




modis_hdf::~modis_hdf() {
	if (thermFlag) {
		delete [] thermdata ;
		delete [] th_scales_offsets ;
		delete [] refldata ;
		delete [] refl_scales_offsets ;
		delete [] raddata_cal ;
		delete [] refdata_cal ;
	}
	else {
		delete [] solzen ;
		delete [] solaz ;
		delete [] senszen ;
		delete [] sensaz ;
		delete [] solsens ;
		delete [] latarr ;
		delete [] lonarr ;
	}
}



int modis_hdf::open_hdf () {
	sd_id = SDstart (inhdfname, DFACC_RDONLY) ;

	if (sd_id <= 0){
		cout << "Problem with : " << inhdfname << endl ;
		return -1 ;
	}
	cout << "opened " << inhdfname << " with id : " << sd_id << endl ;
	n_datasets = 0 ;
	n_fileattrs = 0 ;
	SDfileinfo (sd_id, &n_datasets, &n_fileattrs) ;
	cout << "num datasets " << n_datasets << " file atts : " << 
		n_fileattrs << endl ;

	return (sd_id) ;

}


void modis_hdf::load_geometry () {
	int i, is ;
	float zenval ;
	char name [240], attrb_name[40] ;
	int32 sds_id, rank, dim_sizes[5], attributes, num_type, attnum ;

	int32 startarr [2] = {0,0} ;
	int32 stride [2] = {1,1} ;
	int32 edge [2] = {2030,1354} ;

	for (i=0; i<n_datasets; i++) {
		sds_id = SDselect (sd_id, i) ;
		SDgetinfo (sds_id, name, &rank, dim_sizes, &num_type, &attributes) ;
		//printf ("%d  : %s\r\n", sd_id, name) ;
		if (!strcmp(name, "Latitude")) {
			SDreaddata (sds_id, startarr, stride, edge, latarr) ; 
			SDendaccess (sds_id) ;
		}
		if (!strcmp(name, "Longitude")) {
			SDreaddata (sds_id, startarr, stride, edge, lonarr) ; 
			SDendaccess (sds_id) ;
		}
		if (!strcmp(name, "SolarZenith")) {
			SDreaddata (sds_id, startarr, stride, edge, solzen) ; 
			SDendaccess (sds_id) ;
		}
		if (!strcmp(name, "SolarAzimuth")) {
			SDreaddata (sds_id, startarr, stride, edge, solaz) ; 
			SDendaccess (sds_id) ;
		}
		if (!strcmp(name, "SensorAzimuth")) {
			SDreaddata (sds_id, startarr, stride, edge, sensaz) ; 
			SDendaccess (sds_id) ;
		}
		if (!strcmp(name, "SensorZenith")) {
			SDreaddata (sds_id, startarr, stride, edge, senszen) ; 
			SDendaccess (sds_id) ;
		}
	}

	SDend (sd_id) ;
	cout << "Cent Lat is : " << latarr[1015*1354+670] << endl ;
	cout << "Cent Lon is : " << lonarr[1015*1354+670] << endl ;
	zenval =  solzen[1015*1354+670]/100. ;
	cout << "Cent zenith is : " << zenval << endl ;
	if (zenval > 90) dayflag = false ;
	else dayflag= true ;


}
/**
 * 
 * @return is the aq_terra_flag which distinguished if the file is aqua, terra and nighttime or daytime
 *  aqua (0 night 2 day) terra (3 night 1 day)
 */
void modis_hdf::get_aq_terra_flag () {
    // if daytime, enter this loop
    
    if (dayflag) {
        if (aquaflag) {
            aq_terra_flag = 2 ;
            return ;
        }
        else {
            aq_terra_flag = 1 ;
            return ;
        } 
            
        
    }
    // if still here, this is a nighttime image
    if (aquaflag) {
        aq_terra_flag = 0 ;
        return  ;
    }    
    aq_terra_flag = 3 ;
}
	

void modis_hdf::calib_geometry () {
	int i, npix = 1354L * 2030 ;

	for (i=0; i<npix; i++) {
		solsens[i] = solaz[i] * .01 ;
		solsens[npix + i] = solzen[i] * .01 ;
		solsens[2*npix + i] = sensaz[i] * .01 ;
		solsens[3*npix + i] = senszen[i] * .01 ;
	}
}

void modis_hdf::load_thermal_bands() {
	int i, is ;
	char name [240], attrb_name[40] ;
	int32 sds_id, rank, dim_sizes[5], attributes, num_type, attnum ;

	int32 startarr [3] = {0,0,0} ;
	int32 stride [3] = {1,1,1} ;
	int32 edge [3] = {16,2030,1354} ;

	for (i=0; i<n_datasets; i++) {
		sds_id = SDselect (sd_id, i) ;
		SDgetinfo (sds_id, name, &rank, dim_sizes, &num_type, &attributes) ;
		if (!strcmp (name, "EV_1KM_Emissive")) {
			strcpy (attrb_name, "radiance_scales") ;
			attnum = SDfindattr (sds_id, attrb_name);
			SDreadattr (sds_id, attnum, th_scales_offsets) ;
			strcpy (attrb_name, "radiance_offsets") ;
			attnum = SDfindattr (sds_id, attrb_name);
			SDreadattr (sds_id, attnum, &th_scales_offsets[16]) ;
//			for (is=0; is<16; is++) {
//				cout << "band : "<< is << "  scale : " << th_scales_offsets [is] << 
//					"  offsets : " << th_scales_offsets[16+is] << endl ;
//
//			}
			SDreaddata (sds_id, startarr, stride, edge, thermdata) ;
			SDendaccess (sds_id) ;
			
		}

			
	}	

	//FILE *fout = fopen ("outdat", "w") ;
	//fwrite ((char *)thermdata,2,16*1354l*2030L, fout) ;
	//fclose(fout) ;
			
	SDend (sd_id) ;

}

void modis_hdf::load_refSB_bands() {
	int i, is ;
	char name [240], attrb_name[40] ;
	int32 sds_id, rank, dim_sizes[5], attributes, num_type, attnum ;

	int32 startarr [3] = {0,0,0} ;
	int32 stride [3] = {1,1,1} ;
	int32 edge [3] = {4,2030,1354} ;

	for (i=0; i<n_datasets; i++) {
		sds_id = SDselect (sd_id, i) ;
		SDgetinfo (sds_id, name, &rank, dim_sizes, &num_type, &attributes) ;
		if (!strcmp (name, "EV_500_Aggr1km_RefSB")) {
			strcpy (attrb_name, "radiance_scales") ;
			attnum = SDfindattr (sds_id, attrb_name);
			SDreadattr (sds_id, attnum, refl_scales_offsets) ;
			strcpy (attrb_name, "radiance_offsets") ;
			attnum = SDfindattr (sds_id, attrb_name);
			SDreadattr (sds_id, attnum, &refl_scales_offsets[5]) ;
//			for (is=0; is<5; is++) {
//				cout << "band : "<< is << "  scale : " << refl_scales_offsets [is] << 
//					"  offsets : " << refl_scales_offsets[5+is] << endl ;
//
//			}
			SDreaddata (sds_id, startarr, stride, edge, refldata) ;
			SDendaccess (sds_id) ;
			
		}

			
	}	

			

}


void modis_hdf::calib_thermal_bands () {
	int i, is, bnum ;
	float val, scalev, offv ;
	// radiance = (dn - offset) * scale

	// load up bands 21,22,32 in the calib arrays
	int npix = 2030L * 1354 ;
	int inbands [3] = {1,2,11} ;
	for (i=0 ; i<3 ; i++) {
		bnum = inbands [i] ;
		offv = th_scales_offsets [16+bnum] ;
		scalev =  th_scales_offsets [bnum] ;
		for (is=0; is<npix; is++) {
                       
			val = thermdata[bnum * npix + is] ;
			raddata_cal[i*npix+is]=(val-offv) * scalev ;
		}
	}
	

}

void modis_hdf::calib_refSB_bands () {
	int i, is, bnum ;
	float val, scalev, offv ;
	// radiance = (dn - offset) * scale

	// load up bands 21,22,32 in the calib arrays
	int npix = 2030L * 1354 ;
	// only need band 6
	int inbands [1] = {3} ;
	for (i=0 ; i<1 ; i++) {
		bnum = inbands [i] ;
		offv = refl_scales_offsets [5+bnum] ;
		scalev =  refl_scales_offsets [bnum] ;
		for (is=0; is<npix; is++) {
			val = refldata[bnum * npix + is] ;
			refdata_cal[i*npix+is]=(val-offv) * scalev ;
		}
	}
}



