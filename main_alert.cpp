/*
 * main_alert.cpp
 * get the alert values for all pixels in the image
 * can put them out in a huge text file or in a resampled image file
 */

/* 
 * File:   main.cpp
 * Author: harold
 *
 * Created on May 1, 2018, 12:09 PM
 */

#include <cstdlib>
#include "modis_hdf.h"
#include "modis_process.h"
#include "surftemp.h"
#include <boost/filesystem/path.hpp>


using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {


    vector<lst_coord> mc;
    lst_coord tcoord;
    char outstr[800],  tfile [420], mfile[420], ofile[420], flist[420], outpref[420], outputfile[420];
    char *sufind ;
    int minind ;
    float testtemp;
    float rlat, rlon, pvals[6];
    rlat = 19.333;
    rlon = -155.333;

    int stuff [4];
    modis_hdf *geom, *therm;
    modis_process *mproc;

    // get ASCII file with emiss, geom hdf files
    strcpy (flist, *++argv) ;
    // get output path 
    strcpy (outpref, *++argv) ;
    
    
    cout << "File list is : " << flist << endl ;
    // get MOD021KM file

    //strcpy(tfile, *++argv);
    // get MOD03 file
    //strcpy(mfile, *++argv);
    // get output file
    //strcpy (ofile, *++argv) ;


    boost::filesystem::path p("/home/harold/workdir");
    //char sep = p.preferred_separator ;
    char sep = boost::filesystem::path::preferred_separator;
    outpref [strlen(outpref)] = sep ;
    outpref [strlen(outpref)-1] ='\0' ;
    
    cout << "Separator is " << sep << endl;
   
    FILE *fin = fopen(flist, "r");
   
    if (fin == NULL) {
        cout << " Could not open " << flist << endl;
        exit(-1);
    }
    mproc = new modis_process();
    while (!feof(fin)) {
        fscanf(fin, "%s %s", tfile, mfile);


        therm = new modis_hdf(tfile);
        therm->get_file_name (tfile, ofile) ;
        
        sufind = strstr (ofile, ".hdf") ;
        strncpy (sufind, "_out", 4) ;
        strcpy (outputfile, outpref) ;
        strcat (outputfile, ofile) ;
        cout << "Output file will be " << outputfile << endl ;
        therm->get_date_period(tfile, stuff);
        geom = new modis_hdf(mfile);
        mproc->set_modis_hdfs(geom, therm);
        mproc->set_bounds (20.33, -156.2, 18.8, -154.6, .008) ;
        mproc->process() ;
        mproc->set_month(stuff[1]);
        //mproc->procalert() ;
        mproc->write_output (ofile) ;
        mproc->get_nearest_pixel(rlat, rlon, &minind, pvals);
        cout << "Dayflag is " << geom->dayflag << endl;
        if (!geom->dayflag && pvals[0] < 1200.) {
            sprintf(outstr,
                    "%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f",
                    tfile, stuff[0], stuff[1], stuff[2], stuff[3],
                    pvals[0], pvals[5], pvals[1], pvals[2], pvals[3], pvals[4]);
            //fprintf(fout, "%s\r\n", outstr);
            //fflush(fout);
        }
        cout << flush;
        delete geom;
        delete therm;
        break ;

    }
    delete mproc;

    fclose(fin);

   // fclose(fout);
    return 0;
}

