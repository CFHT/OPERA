/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirMasterTwilightFlat
 Version: 1.0
 Description: Module to create a master flat
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
#include "core-wircam/wirMasterTwilightFlat.h"

#include "libraries/operaException.h"
#include "libraries/Polynomial.h"	
#include "libraries/operaWIRCamImage.h"	

/* \file wirMasterTwilightFlat.cpp */
/* \package core_wircam */

using namespace std;

/* 
 * wirMasterTwilightFlat
 * \author Doug Teeple
 * \brief Module to create a master twilight flat.
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

#define MAXFLATS 20

int main(int argc, char *argv[])
{
	int opt;
	string name_flats[MAXFLATS];
	string name_mastertwilightflat;
	string name_weightmap;
	string name_badpix;
	string param_iiwiversion = "3.0";
	string param_procdate;
	unsigned flat_count = 0;

	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string datafilename;	
	string scriptfilename;
	    
	struct option longopts[] = {
		{"name_flat",					1, NULL, 'i'},	
		{"name_badpix",					1, NULL, 'B'},	
		{"name_mastertwilightflat",		1, NULL, 'o'},	
		{"name_weightmap",				1, NULL, 'w'},	
		
		{"plotfilename",				1, NULL, 'P'},
		{"datafilename",				1, NULL, 'F'},
		{"scriptfilename",				1, NULL, 'S'},  
		
		{"verbose",						0, NULL, 'v'},
		{"debug",						0, NULL, 'd'},
		{"trace",						0, NULL, 't'},
		{"help",						0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:w:B:v::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				name_flats[flat_count++] = optarg;
				break;
			case 'o':
				name_mastertwilightflat = optarg;
				break;
			case 'w':
				name_weightmap = optarg;
				break;
			case 'B':
				name_badpix = optarg;
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
		if (flat_count == 0) {
			throw operaException("wirMasterTwilightFlat: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (name_mastertwilightflat.empty()) {
			throw operaException("wirMasterTwilightFlat: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (name_weightmap.empty()) {
			throw operaException("wirMasterTwilightFlat: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			for (unsigned i=0; i<flat_count; i++) {
				cerr << "wirMasterTwilightFlat: input flat[" << i << "] = " << name_flats[i] << endl; 
			}
			cerr << "wirMasterTwilightFlat: output image = " << name_mastertwilightflat << endl; 
			cerr << "wirMasterTwilightFlat: weight map image = " << name_weightmap << endl; 
			
            if (plot) {
                cerr << "wirMasterTwilightFlat: plotfilename = " << plotfilename << endl;
                cerr << "wirMasterTwilightFlat: datafilename = " << datafilename << endl;
                cerr << "wirMasterTwilightFlat: scriptfilename = " << scriptfilename << endl; 
            }            
		}
		
        ofstream *fdata = NULL;
        
		if (verbose) {
	        cout << "wirMasterTwilightFlat: Reading files in to memory " << endl;
		}
		operaWIRCamImage *flats[MAXFLATS];
		for (unsigned i=0; i<flat_count; i++) {
			flats[i] =  new operaWIRCamImage(name_flats[i], tfloat, READONLY);
		}
		operaWIRCamImage out(name_mastertwilightflat, WIRCAM_NAXIS1, WIRCAM_NAXIS2, 1, WIRCAM_EXTENSIONS, tfloat);
		operaWIRCamImage weightmap(name_weightmap, tfloat, READONLY);
		
        // Copy headers to output files
		out.operaMultiExtensionFITSCubeCopyHeader(flats[0]);
				
		// Create twilight flat, note that the algorithm is the same as master darks...
		out.medianStackParallel(flats, flat_count); 

        // Save the master twilight flat image
		if (verbose) {
	        cout << "wirMasterTwilightFlat: Saving output" << name_mastertwilightflat << endl;
		}
 		out.operaWIRCamImageSave();
		
        // Close files
		for (unsigned i=0; i<flat_count; i++) {
			flats[i]->operaFITSImageClose();
		}
		out.operaFITSImageClose();
		weightmap.operaFITSImageClose();
		
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }          
		
		
        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateMasterTwilightFlatPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(), false);
            }
        }                
	}
	catch (operaException e) {
		cerr << "wirMasterTwilightFlat: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirMasterTwilightFlat: " << operaStrError(errno) << endl;
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
	;		
}

void GenerateMasterTwilightFlatPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, bool display)
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
