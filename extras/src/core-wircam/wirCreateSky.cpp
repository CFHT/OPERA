/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirCreateSky
 Version: 1.0
 Description: Module to use the sky list to create a sky.
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
#include <vector>

#include "operaError.h"
#include "core-wircam/wirCreateSky.h"

#include "libraries/operaException.h"
#include "libraries/Polynomial.h"
#include "libraries/operaWIRCamImage.h"
#include "libraries/operaHelio.h"

/* \file wirCreateSky.cpp */
/* \package core_wircam */

using namespace std;

/*
 * wirCreateSky
 * \author Doug Teeple
 * \brief Module to use the sky list to create a sky.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core_wircam
 */

#define MAXSKIES 10

int main(int argc, char *argv[])
{
	int opt;
	string name_skies[MAXSKIES];
	string name_sky;
	unsigned sky_count = 0;
	string iiwiversion = "3.0";
	string procdate;
	
	vector<string> skylistvector;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;
	string datafilename;
	string scriptfilename;
	
	struct option longopts[] = {
		{"name_sky",			1, NULL, 'o'},
		{"name_skies",			1, NULL, 'y'},
		
		{"plotfilename",		1, NULL, 'P'},
		{"datafilename",		1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},
		
		{"verbose",				0, NULL, 'v'},
		{"debug",				0, NULL, 'd'},
		{"trace",				0, NULL, 't'},
		{"help",				0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "o:y:v::d::t::p::h",
							 longopts, NULL))  != -1) {
		switch(opt) {
            case 'o':
                name_sky = optarg;
            break;
            case 'y':
                name_skies[sky_count++] = optarg;
                if (sky_count > MAXSKIES) {
                    throw operaException("wirCreateSky: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
                }
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
		if (name_sky.empty()) {
			throw operaException("wirCreateSky: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (name_skylist.empty()) {
			throw operaException("wirCreateSky: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (sky_count == 0) {
			throw operaException("wirCreateSky: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
 		if (verbose) {
			for (unsigned i=0; i<sky_count; i++) {
				cout << "wirCreateSky: sky image[" << itos(i) << "] = " << name_skies[i] << endl;
			}
			cout << "wirCreateSky: name_sky = " << name_sky << endl;
            if (plot) {
                cerr << "wirCreateSky: plotfilename = " << plotfilename << endl;
                cerr << "wirCreateSky: datafilename = " << datafilename << endl;
                cerr << "wirCreateSky: scriptfilename = " << scriptfilename << endl;
            }
		}
		
        ofstream *fdata = NULL;
		operaWIRCamImage sky(name_sky, WIRCAM_NAXIS1, WIRCAM_NAXIS2, 1, WIRCAM_EXTENSIONS, tfloat);
        operaWIRCamImage *skies[MAXSKIES];
		for (unsigned n = 0; n < sky_count; n++) {
            skies[i] = new operaWIRCamImage(name_skies[i], tfloat, READONLY);
        }
        
        // Copy headers to output files
		sky.operaMultiExtensionFITSCubeCopyHeader(name_skies[0]);

        sky.createSky(skies, sky_count);

		// Save sky image
		if (verbose) {
			cout << "wirCreateSky: Saving sky as as " << name_sky << endl;
 		}
		sky.operaWIRCamImageSave();
		
        sky.operaFITSImageClose();
		for (unsigned n = 0; n < sky_count; n++) {
            skies[i]->operaFITSImageClose();
        }
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());
        }
		
        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateSkyPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(), false);
            }
        }
	}
	catch (operaException e) {
		cerr << "wirCreateSky: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirCreateSky: " << operaStrError(errno) << endl;
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

void GenerateSkyPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display)
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
