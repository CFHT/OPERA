/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: wirdarkanalysis
 Version: 1.0
 Description: This module performs a statistical analysis of WIRCam darks
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope
 Location: Hawaii USA
 Date: Jan/2011
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
 
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

#include "globaldefines.h"
#include "operaError.h"

#include "libraries/operaException.h"
#include "libraries/operaWIRCamImage.h"

#include "libraries/operaLibCommon.h"
#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"					// for MAXORDERS
#include "libraries/operaFit.h"	

/*! \file wirdarkanalysis.cpp */

using namespace std;

/*! 
 * wirdarkanalysis
 * \author Megan Tannock
 * \brief Performs a statistical analysis of WIRCam darks
 * \arg argc
 * \arg argv
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdthp] <dark filenames>" +
    " --plotfilename=<EPS_FILE>"
	" --datafilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -I, --interactive=<BOOL>\n\n";
}

int main(int argc, char *argv[])
{
	int opt;
	
	string badpixelmaskFile; 
    
    string plotfilename;
	string datafilename;
	string scriptfilename;
    bool interactive = false;
    bool mediancollapse = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
	struct option longopts[] = {
		{"badpixelmask",	1, NULL, 'm'},
		{"mediancollapse",	1, NULL, 'c'},
        
		{"plotfilename",	1, NULL, 'P'},
		{"datafilename",	1, NULL, 'F'},
		{"scriptfilename",	1, NULL, 'S'},
		
		{"interactive",	optional_argument, NULL, 'I'},         
		{"plot",		optional_argument, NULL, 'p'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "m:c:P:F:S:I::p::v::d::t::h",
							 longopts, NULL))  != -1)
	{
		
		switch(opt) 
		{
			case 'm':		// badpixelmask
				badpixelmaskFile = optarg;
				break;
			case 'c':
				mediancollapse = atoi(optarg) == 1;
				break;            
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break;
			case 'F':
				datafilename = optarg;
				break;
			case 'S':
				scriptfilename = optarg;
				break;
			case 'I':		// for interactive plots
				interactive = true;
				break;
				
			case 'p':
				plot = 1; 
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
		
		if (plot && (datafilename.empty() || plotfilename.empty() || scriptfilename.empty())) {
			throw operaException("wirdarkanalysis: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (badpixelmaskFile.empty()) {
			throw operaException("wirdarkanalysis: badpixelmask "+badpixelmaskFile+" ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		
		if (verbose) {
			cout << "wirdarkanalysis: badpixelmaskFile = " << badpixelmaskFile << endl;
			cout << "wirdarkanalysis: mediancollapse = " << mediancollapse << endl;
            if (plot) {
                cout << "wirdarkanalysis: plotfilename = " << plotfilename << endl;
                cout << "wirdarkanalysis: datafilename = " << datafilename << endl;
                cout << "wirdarkanalysis: scriptfilename = " << scriptfilename << endl;
                if(interactive) {
                    cout << "wirdarkanalysis: interactive = YES" << endl;
                } else {
                    cout << "wirdarkanalysis: interactive = NO" << endl;
                }
            }            
		}
        
        ofstream *fdata = NULL;
        
        if (plot) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());
        }
        
        operaMultiExtensionFITSImage badpixelmask(badpixelmaskFile, tfloat, READONLY, cNone, false); // Read in the bad pixel mask
		
		while (optind < argc) {
			string imagefilename = argv[optind++];
			
			string basefilename = imagefilename.substr(imagefilename.find_last_of("/")+1);
			string odometer = basefilename.substr(0, basefilename.find_last_of(".")-1);
			if (verbose) {
                cout << "wirdarkanalysis: darkfilename = " << imagefilename		<< endl;
			}
			unsigned XDimension, YDimension, ZDimension, Extensions;
			edatatype Datatype;
			long Npixels; 
			
			unsigned regionx1 = 48;
			unsigned regionx2 = 2000;
			unsigned regiony1 = 48;
			unsigned regiony2 = 2000; 
			unsigned xsize = regionx2 - regionx1;
			unsigned ysize = regiony2 - regiony1;
			unsigned long myregionsize = xsize*ysize;
			
			float *myregion = new float[myregionsize];
			
			getFITSImageInformation(imagefilename, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
			if (verbose) {
				cout << "wirdarkanalysis: " << imagefilename << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
				cout << "wirdarkanalysis: " << " regionx1= " << regionx1 << " regionx2= " << regionx2  << " regiony1= " << regiony1 << " regiony2= " << regiony2 << " xsize= " << xsize << " ysize= " << ysize << " myregionsize= " << myregionsize << endl;
                cout << "wirdarkanalysis: mediancollapse = " << mediancollapse << endl;
            }
			if (verbose) {
				cout << "wirdarkanalysis: opening " << imagefilename << endl;
			}
			operaWIRCamImage input(imagefilename, tfloat, READONLY, cNone, false);
            float exptime = input.operaFITSGetFloatHeaderValue("EXPTIME");
			if (mediancollapse) {
				if (verbose) {
					cout << "wirdarkanalysis: starting median collapse..." << endl;
                }
				input.medianCollapseParallel(Extensions);
				ZDimension = 1;
				if (verbose) {
					cout << "wirdarkanalysis: median collapse done..." << endl;
                }
			}
			if (verbose) {
				cout << "wirdarkanalysis: subtracting chip bias..." << endl;
			}
			input -= input.getChipBias();
			
			unsigned long counter = 0;
			for (unsigned extension=1; extension<=Extensions; extension++) {
				for (unsigned slice=1; slice<=ZDimension; slice++) {
					unsigned long goodpixelcount = 0;
					operaFITSImage &extensionslice = (operaFITSImage &)(input[extension][slice]);
					float *p = myregion;
					// get the region in an array
					if (verbose) {
						cout << "wirdarkanalysis: extension = " << extension << " slice = " << slice << endl;
					}
					for (unsigned y=regiony1; y<regiony2; y++) {
						for (unsigned x=regionx1; x<regionx2; x++) {
                            if (badpixelmask[extension][y][x] == 1) { // Check if this pixel is flagged as good !CAREFUL! WATCH OUT FOR NANS!
								*p++ = extensionslice[y][x]; 
								goodpixelcount++;
							}
						}
					}
					float stddev = operaArraySigma(goodpixelcount, myregion);
					float regionmedian = operaArrayMedianQuick(goodpixelcount, myregion);
					p = myregion;
					goodpixelcount = 0;
					if (verbose) {
						cout << "wirdarkanalysis: starting good pixel count. stddev = " << stddev << " median = " << regionmedian << endl;
					}
					for (unsigned y=regiony1; y<regiony2; y++) {
						for (unsigned x=regionx1; x<regionx2; x++) {
							if (extensionslice[y][x] > regionmedian - stddev*TwoSigma &&
								extensionslice[y][x] < regionmedian + stddev*TwoSigma) { // Check if this pixel within sigma of the median
								*p++ = extensionslice[y][x]; // If this is not a bad pixel add it to the array
								goodpixelcount++;
							}
						}
					}
					// sigma clip > sigma * sigma coeff. where sigma is the standard deviation
					if (verbose) {
						cout << "wirdarkanalysis: good pixel count = " << goodpixelcount << endl;
					}
					float mean = operaArrayMean(goodpixelcount, myregion);
					float sigma = operaArraySigma(goodpixelcount, myregion);
					float median = operaArrayMedianQuick(goodpixelcount, myregion);
					float sigmaclip = operaArrayAvgSigmaClip(goodpixelcount, myregion, 3);
					float mediansigma = operaArrayMedianSigma(goodpixelcount, myregion, median);
					
					if (verbose) {
						cout << "wirdarkanalysis: " << "ext=" << extension << " slice=" << slice << " mean=" << mean << " median=" << median << " sigma=" << sigma << " 3 sigmaclip=" << sigmaclip << " mediansigma=" << mediansigma << " goodpixelcount=" << goodpixelcount << endl;
					}
					if (plot) {
						counter++;
						*fdata << counter << ' ' << extension << ' ' << slice << ' ' << mean << ' '  << median << ' '  << sigma << ' '  << exptime << ' ' << odometer << ' ' << mediansigma << endl;
					}
				}
				if (plot) {
					*fdata << endl;	// separate extensions by a blank line
				}
			}
			if (plot) {
				*fdata << endl;	// separate files by a blank line
			}
			if (verbose) {
                cout << "wirdarkanalysis: closing " << imagefilename		<< endl;
			}
			input.operaFITSImageClose();
			delete[] myregion;
		}
		if (plot) {
			fdata->close();
		}
	}
	catch (operaException e) {
		cerr << "wirdarkanalysis: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wirdarkanalysis: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

