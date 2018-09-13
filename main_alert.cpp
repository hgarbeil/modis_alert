/**
 * main_alert.cpp
 * Purpose : get the alert values for all pixels in the image
 * can put them out in a huge text file or in a resampled image file
 * @author Harold Garbeil HIGP/SOEST/UHM
 * @version 1.0 6/30/2018
 * 
 **/



#include <cstdlib>
#include "modis_hdf.h"
#include "modis_process.h"
#include "surftemp.h"
#include <boost/filesystem/path.hpp>

 
using namespace std;

/**
 * Main function of the modis_alert program.
 *  
 **/ 
int main(int argc, char** argv) {


    int count ;
    vector<lst_coord> mc;
    lst_coord tcoord;
    char outstr[800],  tfile [420], mfile[420], ofile[420], ofile1[420],flist[420], outpref[420], outputfile[420], alertfile[420], mstdv_file[420] ;
    char proclogfile [420] ;
    char *sufind, *sufind1 ;
    int minind ;
    float testtemp;
    float rlat, rlon, pvals[6];
    rlat = 19.333;
    rlon = -155.333;
    count = 0 ;
    int stuff [4], oldmonth = -1;
    modis_hdf *geom, *therm;
    modis_process *mproc;

    
    
    // get ASCII file with emiss, geom hdf files
    strcpy (flist, *++argv) ;
    cout << "File list is " << flist << endl ;
    // get output path 
    strcpy (outpref, "/home/harold/workdir/data/modis/hawaii") ;
    strcpy (proclogfile, outpref) ;
    strcat (proclogfile, "/proclog.txt\0") ;
    cout << "outpref is " << outpref << endl ;
    // and the log file
    
    
    FILE *fproclog = fopen (proclogfile, "w") ;
    if (fproclog== NULL) {
        cout << "Could not open " << proclogfile << endl ;
        exit (-1) ;
    }
    
    
    cout << "File list is : " << flist << endl ;
    // get MOD021KM file

    //strcpy(tfile, *++argv);
    // get MOD03 file
    //strcpy(mfile, *++argv);
    // get output file
    //strcpy (ofile, *++argv) ;


/**     outpref [strlen(outpref)] = sep ;
        outpref [strlen(outpref)-1] ='\0' ;
*/  
    FILE *fin = fopen(flist, "r");
   
    if (fin == NULL) {
        cout << " Could not open file list :  " << flist << endl;
        exit(-1);
    }
    // controlling process class
    mproc = new modis_process();
    
    
    // go through each pair of MODIS MOD021KM and corresponding MOD03 files
    // calculating NTI values for each pixel in the 1354 x 2030 original array
    // as well as resampling radiance bands and NTI values to the predetermined
    // Hawaii ROI
    while (!feof(fin) && count < 3) {
        fscanf(fin, "%s %s", tfile, mfile);
        

        therm = new modis_hdf(tfile);
        therm->get_file_name (tfile, ofile) ;
        strcpy (ofile1, ofile) ;
        
        sufind = strstr (ofile, ".hdf") ;
        
        
        strncpy (sufind, "_out", 4) ;
        strcpy (outputfile, outpref) ;
        strcat (outputfile, "/") ;
        
        strcat (outputfile, ofile) ;
        cout << "Output file will be " << outputfile << endl ;
        
        // alert text file
        sufind = strstr (ofile1, ".hdf") ;
        strncpy (sufind, "_alert.txt", 10) ;
        strcpy (alertfile, outpref) ;
        strcat (alertfile, "/") ;
        strcat (alertfile, ofile1) ;
         cout << "Alert file will be " << alertfile << endl ;
        
        
         
        
        therm->get_date_period(tfile, stuff);
        cout << "HDF file month is : " << stuff [1] << endl ;
        geom = new modis_hdf(mfile);
        mproc->set_month(stuff[1]-1);
        
        mproc->set_modis_hdfs(geom, therm);
        if (count==0)
        mproc->set_bounds (22.0, -157.8, 18.804, -154.6, .008) ;
        if (geom->dayflag) continue ;
        if (mproc->mon !=oldmonth ) {
            mproc->extract_from_baseline_file() ;
            oldmonth = mproc->mon ;
        }
        // calc alert is the full scene alert calculation- this fills the alinds vector
        mproc->calc_alert(alertfile) ;
        
        
        
        mproc->process() ;
        
        //mproc->procalert() ;
        
        mproc->get_nearest_pixel(rlat, rlon, &minind, pvals);
        cout << "Dayflag is " << geom->dayflag << endl;
        
        //strcpy (mstdv_file, "/local/worldbase/02ssh/rad32_001_00.bsq") ;
        // get the alert from MODIS file
        
        
        mproc->zscore() ;
        mproc->write_output (outputfile) ;
        cout << flush;
        sprintf (outstr, "%s\t%s\t%d", tfile, mproc->sp->basef, mproc->num_alerts_full) ;
        fprintf (fproclog, "%s\n", outstr) ;
        delete geom;
        delete therm; 
        
        count++ ;

    }
    delete mproc;

    fclose(fin);
    fclose (fproclog) ;

   // fclose(fout);
    return 0;
}

