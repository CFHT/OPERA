/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirSubtractSky
 Version: 1.0
 Description: Module to subtract a sky
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
#include "core-wircam/wirSubtractSky.h"

#include "libraries/operaException.h"
#include "libraries/Polynomial.h"	
#include "libraries/operaWIRCamImage.h"	

/* \file wirSubtractSky.cpp */
/* \package core_wircam */

using namespace std;

/* 
 * wirSubtractSky
 * \author Doug Teeple
 * \brief Module to subtract a sky.
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

/*
 * scaletoexptime
 * If set, it uses the exposure time ratio (sky frame / current image) rather than the
 * current image median to subtract the sky. To be used for the nodding strategy where
 * the current image (assumed with an extended source) does not represent the atmopsheric
 * sky intensity. Default is param_scaletoexptime = 0 for DP, STARING. For WDP, it is a
 * grey area but has been turned to 0 for all 2005-2009 processing.
 * Example: if the sky image has exptime of 10sec and sky median of 5000 adu, and if the
 * image to sky subtract has an exptime of 15sec then the sky will be scaled to 7500 adu
 * and subtracted from the image to sky subtract, regardless of the fact that the image
 * may have a sky median different than 7500 adu.
 *
 * skylevel 
 * A matrix passed as input that gives the sky background values for all 4 extensions.
 * matrix_skylvl[n,ext]. where n is the number of slices in the cube. Can only use this 
 * technique one image at at time, i.e. if nimage = 1. It uses the list_sky fits file.
 */
#define MAXSKIES 20

int main(int argc, char *argv[])
{
	int opt;
	string name_input;
	string name_skysubtracted;
	string name_sky;
	string name_badpix;
	string param_iiwiversion = "3.0";
	string param_procdate;
	
	bool scaletoexptime = false;
	float etime = 1.0;
	unsigned int sky_count = 0;

	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string datafilename;	
	string scriptfilename;
	    
	struct option longopts[] = {
		{"name_input",			1, NULL, 'i'},	
		{"name_skysubtracted",	1, NULL, 's'},	
		{"name_badpix",			1, NULL, 'B'},	
		{"name_sky",			1, NULL, 'o'},	
		{"skycount",			1, NULL, 'y'},	
		{"scaletoexptime",		1, NULL, 'c'},	
		{"etime",				1, NULL, 'e'},	
		
		{"plotfilename",		1, NULL, 'P'},
		{"datafilename",		1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},  
		
		{"verbose",				0, NULL, 'v'},
		{"debug",				0, NULL, 'd'},
		{"trace",				0, NULL, 't'},
		{"help",				0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:s:B:c:e:y:v::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				name_input = optarg;
				break;
			case 'o':
				name_sky = optarg;
				break;
			case 's':
				name_skysubtracted = optarg;
				break;
			case 'B':
				name_badpix = optarg;
				break;
			case 'c':
				scaletoexptime = true;
				break;
			case 'e':
				if (scaletoexptime)
					etime = atof(optarg);
				break;
			case 'y':
					sky_count = atoi(optarg);
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
			throw operaException("wirSubtractSky: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (name_sky.empty()) {
			throw operaException("wirSubtractSky: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (name_skysubtracted.empty()) {
			throw operaException("wirSubtractSky: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cerr << "wirSubtractSky: input image = " << name_input << endl; 
			cerr << "wirSubtractSky: skysubtracted image = " << name_skysubtracted << endl; 
			cerr << "wirSubtractSky: sky image = " << name_sky << endl; 
			cerr << "wirSubtractSky: sky_count = " << sky_count << endl; 
			cerr << "wirSubtractSky: scaletoexptime = " << scaletoexptime << endl; 
			cerr << "wirSubtractSky: etime = " << ftos(etime) << endl; 
			
            if (plot) {
                cerr << "wirSubtractSky: plotfilename = " << plotfilename << endl;
                cerr << "wirSubtractSky: datafilename = " << datafilename << endl;
                cerr << "wirSubtractSky: scriptfilename = " << scriptfilename << endl; 
            }            
		}
		
        ofstream *fdata = NULL;
        
		operaWIRCamImage in(name_input, tfloat, READONLY);
		operaWIRCamImage sky(name_sky, tfloat, READONLY);
		operaWIRCamImage out(name_skysubtracted, WIRCAM_NAXIS1, WIRCAM_NAXIS2, 1, WIRCAM_EXTENSIONS, tfloat);

        // Copy headers to output files
		out.operaMultiExtensionFITSCubeCopyHeader(&in);
		
		// Copy raw pixels to out
        out = in;
		out -= out.getChipBias();
		out.skySubtraction(out, sky, etime);
		
		// add new headers
		out.operaFITSSetHeaderValue("N_ODOSKY", itos(sky_count), "Nbr of odometers used in constructing sky");
		if (scaletoexptime) {
			out.operaFITSSetHeaderValue("SKYSBTEC", "SCALE_TO_EXPTIME", "Sky Subtraction Technique");
		}

        // Save detrended image
 		out.operaWIRCamImageSave();
		
        // Close files
		in.operaFITSImageClose();
		sky.operaFITSImageClose();
		out.operaFITSImageClose();
		
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
		cerr << "wirSubtractSky: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirSubtractSky: " << operaStrError(errno) << endl;
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
