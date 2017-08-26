/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirPhotometry
 Version: 1.0
 Description: Apply Photometry.
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
#include "core-wircam/wirPhotometry.h"

#include "libraries/operaException.h"
#include "libraries/Polynomial.h"	
#include "libraries/operaWIRCamImage.h"	

/* \file wirPhotometry.cpp */
/* \package core_wircam */

using namespace std;

/* 
 * wirPhotometry
 * \author Doug Teeple
 * \brief Module to Apply Photometry.
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
 * This program creates a matched catalogue for every slice of a cube assuming
 * the WCS in the header is good and determines zero point and measure the
 * absorption. This program is used in conjonction with wirAstrometry.pro
 * which does the WCS.
 */
string default_config = 
" -DETECT_MINAREA 3"
" -THRESH_TYPE RELATIVE"
" -DETECT_THRESH 2.0"
" -ANALYSIS_THRESH 2.0"
" -WEIGHT_TYPE NONE"
" -SEEING_FWHM 0.7"
" -CHECKIMAGE_TYPE NONE"
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
	string param_iiwiversion = "3.0";
	string param_procdate;

	string sextractor_configuration_parameters = default_config;
	string matchedcatalogname = "/tmp/sex.matchedcat";
	string configfilename = "/tmp/config.sex";
	string twomasscatalogname = "/tmp/twomass.cat";
	string sexcatalogname = "/tmp/sex.cat";
	string parametersfilename = "/tmp/sex.param";
	string tmc_path_string = "/data/polena/catalogs/2mass.wcs";
	bool linearize = false;
	float zp_base  = 0.0;
	float zp_err   = 0.0;
	float zp_wonlc = 0.0;
	float zp_ext1  = 0.0;
	float zp_ext2  = 0.0;
	float zp_ext3  = 0.0;
	float zp_ext4  = 0.0;
	float detect_threshold = 2.0;
	float detect_minarea = 3.0;
	float analysis_thresh = 2.0;

	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string datafilename;	
	string scriptfilename;
	    
	struct option longopts[] = {
		{"name_input",			1, NULL, 'i'},	
		{"name_output",			1, NULL, 'o'},	
		{"name_sexcat",			1, NULL, 's'},	
		{"configfilename",		1, NULL, 'c'},
		{"twomasscatalogname",  1, NULL, 'w'},
		{"sexcatalogname",		1, NULL, 'x'},
		{"parametersfilename",	1, NULL, 'a'},
		{"tmc_path",			1, NULL, 'm'},
		{"name_matchedcat",		1, NULL, 'e'},	
		{"detect_threshold",	1, NULL, '0'},	
		{"detect_minarea",		1, NULL, '1'},	
		{"analysis_thresh",		1, NULL, '2'},	
		{"zp_base",				1, NULL, '3'},	
		{"zp_err",				1, NULL, '4'},	
		{"zp_wonlc",			1, NULL, '5'},	
		{"zp_ext1",				1, NULL, '6'},	
		{"zp_ext2",				1, NULL, '7'},	
		{"zp_ext3",				1, NULL, '8'},	
		{"zp_ext4",				1, NULL, '9'},	
		{"linearize",			1, NULL, 'l'},	
		
		{"plotfilename",		1, NULL, 'P'},
		{"datafilename",		1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},  
		
		{"verbose",				0, NULL, 'v'},
		{"debug",				0, NULL, 'd'},
		{"trace",				0, NULL, 't'},
		{"help",				0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:s:c:w:x:a:m:e:0:1:2:3:4:5:6:7:8:9:l:v::d::t::p::h", 
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
			case 's':
				sexcatalogname = optarg;
				break;
			case 'e':
				matchedcatalogname = optarg;
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
			case '0':
				detect_threshold = atof(optarg);
				sextractor_configuration_parameters += " -DETECT_THRESH "+string(optarg);
				break;
			case '1':
				detect_minarea = atof(optarg);
				sextractor_configuration_parameters += " -DETECT_MINAREA "+string(optarg);
				break;
			case '2':
				analysis_thresh = atof(optarg);
				sextractor_configuration_parameters += " -ANALYSIS_THRESH "+string(optarg);
				break;
			case '3':
				zp_base = atof(optarg);
				break;
			case '4':
				zp_err = atof(optarg);
				break;
			case '5':
				zp_wonlc = atof(optarg);
				break;
			case '6':
				zp_ext1 = atof(optarg);
				break;
			case '7':
				zp_ext2 = atof(optarg);
				break;
			case '8':
				zp_ext3 = atof(optarg);
				break;
			case '9':
				zp_ext4 = atof(optarg);
				break;
			case 'l':
				linearize = atoi(optarg) == 1;
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
			throw operaException("wirPhotometry: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cerr << "wirPhotometry: input image = " << name_input << endl; 
			cerr << "wirPhotometry: output image = " << name_output << endl; 
			cerr << "wirPhotometry: sexcatalogname = " << sexcatalogname << endl; 
			cerr << "wirPhotometry: matchedcatalogname = " << matchedcatalogname << endl; 
			cout << "wirPhotometry: config file name = " << configfilename << endl;
			cout << "wirPhotometry: parameters file name = " << parametersfilename << endl;
			cout << "wirPhotometry: two mass file name = " << twomasscatalogname << endl;
			cout << "wirPhotometry: sextractor file name = " << sexcatalogname << endl;
			cout << "wirPhotometry: tmc_path = " << tmc_path_string << endl;			
			cerr << "wirPhotometry: detect_threshold = " << detect_threshold << endl; 
			cerr << "wirPhotometry: detect_minarea = " << detect_minarea << endl; 
			cerr << "wirPhotometry: analysis_thresh = " << analysis_thresh << endl; 
			cerr << "wirPhotometry: zp_base = " << zp_base << endl; 
			cerr << "wirPhotometry: zp_err = " << zp_err << endl; 
			cerr << "wirPhotometry: zp_wonlc = " << zp_wonlc << endl; 
			cerr << "wirPhotometry: zp_ext1 = " << zp_ext1 << endl; 
			cerr << "wirPhotometry: zp_ext2 = " << zp_ext2 << endl; 
			cerr << "wirPhotometry: zp_ext3 = " << zp_ext3 << endl; 
			cerr << "wirPhotometry: zp_ext4 = " << zp_ext4 << endl; 
			cerr << "wirPhotometry: linearize = " << linearize << endl; 
            if (plot) {
                cerr << "wirPhotometry: plotfilename = " << plotfilename << endl;
                cerr << "wirPhotometry: datafilename = " << datafilename << endl;
                cerr << "wirPhotometry: scriptfilename = " << scriptfilename << endl; 
            }            
		}
		
        ofstream *fdata = NULL;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }          
		
		char *tmc_path = getenv("TMC_PATH");		// The environment variable TMC_PATH must be defined for the 2MASS catalogue to be queried.
		if (tmc_path == NULL) {
			throw operaException("wirPhotometry: The environment variable TMC_PATH must be defined for the 2MASS catalogue to be queried ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		} else if (tmc_path_string.empty()) {
			tmc_path_string = string(tmc_path);
		}
		
        if (verbose) {
			cout << "wirPhotometry: writing default configuration file: " << configfilename << endl;
		}
		ofstream fout(configfilename.c_str());
		fout << default_config;
		fout.close();
		
        if (verbose) {
			cout << "wirPhotometry: writing default parameters file: " << parametersfilename << endl;
		}
		ofstream pout(parametersfilename.c_str());
		pout << default_parameters;
		pout.close();

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
			cout << "wirPhotometry: sex catalog:" << endl;
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

        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GeneratePhotometryPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(), false);
            }
        }                
	}
	catch (operaException e) {
		cerr << "wirPhotometry: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirPhotometry: " << operaStrError(errno) << endl;
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

void GeneratePhotometryPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display)
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
