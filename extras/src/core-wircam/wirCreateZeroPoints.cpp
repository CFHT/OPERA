/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirCreateZeroPoints
 Version: 1.0
 Description: Module to Create zero point calibration
 to start up with an OPERA module. 
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jul/2013
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2013  Opera Pipeline team, Canada France Hawaii Telescope
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see:
 http://software.cfht.hawaii.edu/licenses
 -or-
 http://www.gnu.org/licenses/gpl-3.0.html
 ********************************************************************/

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include <stdio.h>
#include <getopt.h>
#include <fstream>

#include "operaError.h"
#include "core-wircam/wirCreateZeroPoints.h"

#include "libraries/operaException.h"
#include "libraries/Polynomial.h"	
#include "libraries/operaWIRCamImage.h"	

/* \file wirCreateZeroPoints.cpp */
/* \package core_wircam */

using namespace std;

/* 
 * wirCreateZeroPoints
 * \author Doug Teeple
 * \brief Module to Create zero point calibration.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */
/*
 */

string default_config = 
"CATALOG_TYPE		ASCII_HEAD\n"					// "ASCII_HEAD","ASCII","FITS_1.0" or "FITS_LDAC"
"DETECT_TYPE		CCD\n"
"DETECT_MINAREA		5.0\n"
"THRESH_TYPE		RELATIVE\n"
"DETECT_THRESH		1.5\n"
"ANALYSIS_THRESH	1.5\n"
"FILTER				Y\n"
"FILTER_NAME		~/opera-1.0/config/default.conv\n"
"FILTER_THRESH		0\n"
"DEBLEND_NTHRESH	32\n"
"DEBLEND_MINCONT	0.05\n"
"CLEAN				Y\n"
"CLEAN_PARAM		1.0\n"
"MASK_TYPE			CORRECT\n"
"WEIGHT_TYPE		NONE\n"
"WEIGHT_IMAGE		weight.fits\n"
"WEIGHT_GAIN		N\n"
"FLAG_IMAGE			flag.fits\n"
"FLAG_TYPE			OR\n"
"PHOT_APERTURES		5\n"
//"PHOT_AUTOPARAMS[0]	2.5\n"
//"PHOT_AUTOPARAMS[1]	3.5\n"
//"PHOT_PETROPARAMS[0] 2.0\n"
//"PHOT_PETROPARAMS[1] 3.5\n"
//"PHOT_AUTOAPERS[0]	0.0\n"
//"PHOT_AUTOAPERS[1]	0.0\n"
"PHOT_FLUXFRAC		0.5\n"
"SATUR_LEVEL		50000.0\n"
"SATUR_KEY			SATURATE\n"
"MAG_ZEROPOINT		0.0\n"
"MAG_GAMMA			4.0\n"
"GAIN				0.0\n"
"GAIN_KEY			GAIN\n"
"PIXEL_SCALE		1.0\n"
"SEEING_FWHM		1.2\n"
"STARNNW_NAME		~/opera-1.0/config/default.nnw\n"
"BACK_TYPE			AUTO\n"
"BACK_VALUE			0.0\n"
"BACK_SIZE			64\n"
"BACK_FILTERSIZE	3\n"
"BACKPHOTO_TYPE		GLOBAL\n"
"BACKPHOTO_THICK	24\n"
"BACK_FILTTHRESH	0.0\n"
"CHECKIMAGE_TYPE	NONE\n"
"CHECKIMAGE_NAME	check.fits\n"
"MEMORY_OBJSTACK	3000\n"
"MEMORY_PIXSTACK	300000\n"
"MEMORY_BUFSIZE		1024\n"
//"ASSOC_NAME			sky.list\n"
//"ASSOC_DATA[0]		2\n"
//"ASSOC_DATA[1]		3\n"
//"ASSOC_DATA[2]		4\n"
//"ASSOC_PARAMS[0]	2\n"
//"ASSOC_PARAMS[1]	3\n"
//"ASSOC_PARAMS[2]	4\n"
"ASSOC_RADIUS		2.0\n"
"ASSOC_TYPE			NEAREST\n"
"ASSOCSELEC_TYPE	MATCHED\n"
"VERBOSE_TYPE		NORMAL\n"
"WRITE_XML			N\n"
"XML_NAME			sex.xml\n"
"XSL_URL			~/opera-1.0/config/sextractor.xsl\n"
"NTHREADS			0\n"
"FITS_UNSIGNED		N\n"
"INTERP_MAXXLAG		16\n"
"INTERP_MAXYLAG		16\n"
"INTERP_TYPE		ALL\n"
;

string default_parameters = 
"X_IMAGE\n"
"Y_IMAGE\n"
"MAG_BEST\n"
;

int main(int argc, char *argv[])
{
	int opt;
	string name_output;
	string name_input;
	string configfilename = "/tmp/config.sex";
	string twomasscatalogname = "/tmp/twomass.cat";
	string sexcatalogname = "/tmp/sex.cat";
	string parametersfilename = "/tmp/sex.param";
	string tmc_path_string = "/data/polena/catalogs/2mass.wcs";
	string param_iiwiversion = "3.0";
	string param_procdate;
	    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string datafilename;	
	string scriptfilename;
	
	struct option longopts[] = {
		{"name_input",			optional_argument, NULL, 'i'},	
		{"name_output",			optional_argument, NULL, 'o'},	
		{"configfilename",		optional_argument, NULL, 'c'},
		{"twomasscatalogname",  optional_argument, NULL, 'w'},
		{"sexcatalogname",		optional_argument, NULL, 'x'},
		{"parametersfilename",	optional_argument, NULL, 'a'},
		{"tmc_path",			optional_argument, NULL, 'm'},
		
		{"plotfilename",		optional_argument, NULL, 'P'},
		{"datafilename",		optional_argument, NULL, 'F'},
		{"scriptfilename",		optional_argument, NULL, 'S'},  
		
		{"verbose",				0, NULL, 'v'},
		{"debug",				0, NULL, 'd'},
		{"trace",				0, NULL, 't'},
		{"help",				0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:B:o:c:w:x:a:m:v::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				name_input = optarg;
				break;
			case 'o':
				name_output = optarg;
				break;
			case 'c':
				configfilename = optarg;
				break;
			case 'w':
				twomasscatalogname = optarg;
				break;
			case 'x':
				sexcatalogname = optarg;
				break;
			case 'm':
				tmc_path_string = optarg;
				break;
			case 'a':
				parametersfilename = optarg;
				break;
				
			case 'v':
				verbose = 1;
				break;
			case 'd':
				debug = 1;
				break;
			case 't':
				trace = 1;
				break;         
			case 'h':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			case '?':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
		}
	}	
	
	try {
		if (name_input.empty()) {
			throw operaException("wirCreateZeroPoints: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cerr << "wirCreateZeroPoints: input image = " << name_input << endl; 
			cerr << "wirCreateZeroPoints: output image = " << name_output << endl; 
			cout << "wirCreateZeroPoints: config file name = " << configfilename << endl;
			cout << "wirCreateZeroPoints: parameters file name = " << parametersfilename << endl;
			cout << "wirCreateZeroPoints: two mass file name = " << twomasscatalogname << endl;
			cout << "wirCreateZeroPoints: sextractor file name = " << sexcatalogname << endl;
			cout << "wirCreateZeroPoints: tmc_path = " << tmc_path_string << endl;			
			if (plot) {
                cerr << "wirCreateZeroPoints: plotfilename = " << plotfilename << endl;
                cerr << "wirCreateZeroPoints: datafilename = " << datafilename << endl;
                cerr << "wirCreateZeroPoints: scriptfilename = " << scriptfilename << endl; 
            }            
		}
		
        ofstream *fdata = NULL;
        
		char *tmc_path = getenv("TMC_PATH");		// The environment variable TMC_PATH must be defined for the 2MASS catalogue to be queried.
		if (tmc_path == NULL) {
			throw operaException("wirCreateZeroPoints: The environment variable TMC_PATH must be defined for the 2MASS catalogue to be queried ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		} else if (tmc_path_string.empty()) {
			tmc_path_string = string(tmc_path);
		}
		
        if (verbose) {
			cout << "wirCreateZeroPoints: writing default configuration file: " << configfilename << endl;
		}
		ofstream fout(configfilename.c_str());
		fout << default_config;
		fout.close();
		
        if (verbose) {
			cout << "wirCreateZeroPoints: writing default parameters file: " << parametersfilename << endl;
		}
		ofstream pout(parametersfilename.c_str());
		pout << default_parameters;
		pout.close();
		
        if (verbose) {
			cout << "wirCreateZeroPoints: Calling sextractor." << endl;
		}
		string basename = name_input;
		if (name_input.find("/") != string::npos) {
			basename = name_input.substr(name_input.find_last_of("/")+1);
		}
		sexcatalogname = basename.substr(0, basename.find(".fits")) + ".cat";
		if (verbose) {
			printf("%s/sex %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s -FILTER_NAME %s -THRESH_TYPE RELATIVE -WEIGHT_TYPE NONE -SATUR_LEVEL %d -CHECKIMAGE_TYPE NONE >/tmp/sex.log\n",
				   "/usr/local/bin/",
				   name_input.c_str(),
				   configfilename.c_str(), 
				   parametersfilename.c_str(), 
				   sexcatalogname.c_str(),
				   "~/opera-1.0/config/default.conv",
				   35000);			
		}
		systemf("%s/sex %s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s -FILTER_NAME %s -THRESH_TYPE RELATIVE -WEIGHT_TYPE NONE -SATUR_LEVEL %d -CHECKIMAGE_TYPE NONE >/tmp/sex.log",
				"/usr/local/bin/",
				name_input.c_str(),
				configfilename.c_str(), 
				parametersfilename.c_str(), 
				sexcatalogname.c_str(),
				"~/opera-1.0/config/default.conv",
				35000);
		
		if (verbose) {
			cout << "wirCreateZeroPoints: sex catalog:" << endl;
			ifstream sexcat(sexcatalogname.c_str());
			if (sexcat.is_open()) {
				string dataline;
				while (sexcat.good()) {
					getline(sexcat, dataline);
					cout << dataline << endl;
				} 
				sexcat.close();
			}
		}
		
		if (verbose) {
			cout << "wirCreateZeroPoints: Calling scat." << endl;
		}
		operaWIRCamImage image(name_input, tfloat, READONLY, cNone, false); 
		float absra_center = image.operaFITSGetFloatHeaderValue("RA_DEG");
		float absdec_center = image.operaFITSGetFloatHeaderValue("DEC_DEG");
		// Apply the ra-dec offset between the science target (RA_DEG,DEC_DEG) and the instrument mosaic center
		float instzra = image.operaFITSGetFloatHeaderValue("INSTZRA"); // Should be 52.0 arcsec	NOTE: the signs are inconsistent!
		float instzdec = image.operaFITSGetFloatHeaderValue("INSTZDEC"); // Should be -52.0 arcsec
		
		absra_center = absra_center - instzra/3600.0;			// The header interprets this as the instrument center
		absdec_center = absdec_center + instzdec/3600.0;		// Yes, the header signs are inconsistent!
		
		if (verbose) {
			printf("export TMC_PATH=%s ; %s/scat -c tmc -j -d -r 700,700 -n 2000 %f %f >%s\n",
				   tmc_path_string.c_str(),
				   "/usr/local/bin/",
				   absra_center,
				   absdec_center, 
				   twomasscatalogname.c_str());				
		}
		systemf("export TMC_PATH=%s ; %s/scat -c tmc -j -d -r 700,700 -n 2000 %f %f >%s",
				tmc_path_string.c_str(),
				"/usr/local/bin/",
				absra_center,
				absdec_center, 
				twomasscatalogname.c_str());
		
		if (verbose) {
			cout << "wirCreateZeroPoints: scat catalog:" << endl;
			ifstream twomass(twomasscatalogname.c_str());
			if (twomass.is_open()) {
				string dataline;
				while (twomass.good()) {
					getline(twomass, dataline);
					cout << dataline << endl;
				}
				twomass.close();
			}
		}

        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }          
		
		
        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateZeroPointPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(), false);
            }
        }                
	}
	catch (operaException e) {
		cerr << "wirCreateZeroPoints: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirCreateZeroPoints: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cerr <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n" 
	"  -I, --interactive=<BOOL>\n\n";		
}

void GenerateZeroPointPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display)
{
    FILE *fgnu;
    remove(gnuScriptFileName); // delete any existing file with the same name
	
    fgnu = fopen(gnuScriptFileName,"w");
    
    fprintf(fgnu,"unset key\n");
    fprintf(fgnu,"set view 0,0\n");
    fprintf(fgnu,"set iso 100\n");
    fprintf(fgnu,"set samples 100\n");
    fprintf(fgnu,"set pm3d at s\n");
    fprintf(fgnu,"set ticslevel 0\n");   
    
    fprintf(fgnu,"set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
    fprintf(fgnu,"set output \"%s\"\n",outputPlotEPSFileName);
	
    if (display) {
		fprintf(fgnu,"set output\n");
		fprintf(fgnu,"set terminal x11\n");
		fprintf(fgnu,"replot\n");       
		fclose(fgnu);   
		systemf("gnuplot -persist %s",gnuScriptFileName);
    } else {
		fclose(fgnu);  
	}
}
