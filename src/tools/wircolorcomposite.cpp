/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: wircolorcomposite
 Version: 1.0
 Description: Perform various tests on the operaMEFFITSImage and MEFFITSCube classes.
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

#include <string.h>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMultiExtensionFITSImage.h"
#include "libraries/operaMultiExtensionFITSCube.h"
#include "libraries/operaFITSCube.h"
#include "libraries/operaWIRCamImage.h"

/*! \file wircolorcomposite.cpp */

using namespace std;

/*! 
 * wircolorcomposite
 * \author Doug Teeple / Megan Tannock
 * \brief create a color WIRCam image from 3 images.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup tools
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cout << " Usage: wircolorcomposite --R=<filename> --G=<filename> --B=<filename> --extension=[1..4] --badpixelmask==<filename> --output=<file name> --split -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	string r, g, b;
	float RPercent = 100.0;
	float GPercent = 100.0;
	float BPercent = 100.0;
	string output;
	string badpixelmask;
	unsigned extension = 4;
	bool split = false;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"R",				1, NULL, 'R'},	// red input image
		{"G",				1, NULL, 'G'},	// green input image
		{"B",				1, NULL, 'B'},	// blue input image
		{"r",				1, NULL, 'r'},	// red %
		{"g",				1, NULL, 'g'},	// green %
		{"b",				1, NULL, 'b'},	// blue %
		{"extension",		1, NULL, 'e'},	// extension (1 .. 4) where the object is, default 4
		{"badpixelmask",	1, NULL, 'x'},	// badpixelmask
		{"split",			0, NULL, 's'},	// split into 3 images
		{"output",			1, NULL, 'o'},	// output fits file
		
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "R:G:B:r:g:b:e:x:so:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'R':
					r = optarg;
					break;
				case 'G':
					g = optarg;
					break;
				case 'B':
					b = optarg;
					break;
				case 'r':
					RPercent = atof(optarg);
					break;
				case 'g':
					GPercent = atof(optarg);
					break;
				case 'b':
					BPercent = atof(optarg);
					break;
				case 'o':
					output = optarg;
					break;
				case 's':
					split = true;
					break;
				case 'x':
					badpixelmask = optarg;
					break;
				case 'e':
					extension = atoi(optarg);
					break;
					
				case 'v':
					verbose = 1;
					break;
				case 'p':
					plot = 1;
					break;
				case 'd':
					debug = 1;
					break;
				case 't':
					trace = 1; 
					break;         
				case 'h':
					printUsageSyntax();
					exit(EXIT_SUCCESS);
					break;
				default:
					printUsageSyntax();
					exit(EXIT_SUCCESS);
					break;
			}	// switch
		}	// while
		
        unsigned XDimension, YDimension, ZDimension, Extensions;
		edatatype Datatype;
		long Npixels;
        
		if (badpixelmask.empty()) {
			throw operaException("wircolorcomposite: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (output.empty()) {
			throw operaException("wircolorcomposite: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (r.empty()) {
			throw operaException("wircolorcomposite: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (g.empty()) {
			throw operaException("wircolorcomposite: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.empty()) {
			throw operaException("wircolorcomposite: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        operaWIRCamImage red(r, tfloat, READONLY);
        operaWIRCamImage green(g, tfloat, READONLY);
        operaWIRCamImage blue(b, tfloat, READONLY);
        operaWIRCamImage badpix(badpixelmask, tfloat, READONLY);
        operaFITSCube out(output, red.getXDimension(), red.getYDimension(), 3, tfloat, cNone);
		out.operaFITSImageCopyHeader(&red);	

        // Get some information about the image...
		getFITSImageInformation(r, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
        if (verbose) {
			cout << "wircolorcomposite: red: filter:" << red.operaFITSGetHeaderValue("FILTER") << " " << r << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;        
		}
		getFITSImageInformation(g, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
        if (verbose) {
			cout << "wircolorcomposite: green: filter: " << green.operaFITSGetHeaderValue("FILTER") << " " << g << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;        
		}
		getFITSImageInformation(b, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
        if (verbose) {
			cout << "wircolorcomposite: blue: ilter: " << blue.operaFITSGetHeaderValue("FILTER") << " " << b << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;        
		}
		// median collapse to get rid of cosmic rays
		red.medianCollapse();
		green.medianCollapse();
		blue.medianCollapse();

		// create FITS images of the extension
		operaFITSImage redfits(((operaMultiExtensionFITSImage &)red)[extension]);
		operaFITSImage greenfits(((operaMultiExtensionFITSImage &)green)[extension]);
		operaFITSImage bluefits(((operaMultiExtensionFITSImage &)blue)[extension]);
		operaFITSImage badpixfits(((operaMultiExtensionFITSImage &)badpix)[extension]);
		
		/*
		 * 1. remove the chip bias
		 * 2. get the image median
		 * 3. set the inverse of the badpix to image median
		 * 4. mask the guidewindow to image median
		 * 5. adjust the rgb pixel level by the rgb levels on command line
		 */
		redfits -= red.getChipBias();
		float imageMedian = operaArrayMedian(redfits.getnpixels(), (float *)redfits.getpixels());
		redfits[where(!badpixfits)] = imageMedian;
		redfits[new operaImageVector(red.getGuideWindow(extension),1)] = imageMedian;
		redfits *= RPercent / 100.0;
        if (verbose)
			cout << "wircolorcomposite: median red: " << imageMedian << endl;        
		
		greenfits -= green.getChipBias();
		imageMedian = operaArrayMedian(greenfits.getnpixels(), (float *)greenfits.getpixels());
		greenfits[where(!badpixfits)] = imageMedian;
		greenfits[new operaImageVector(green.getGuideWindow(extension),1)] = imageMedian;
		greenfits *= GPercent / 100.0;
        if (verbose)
			cout << "wircolorcomposite: median green: " << imageMedian << endl;        
		
		bluefits -= blue.getChipBias();
		imageMedian = operaArrayMedian(bluefits.getnpixels(), (float *)bluefits.getpixels());
		bluefits[where(!badpixfits)] = imageMedian;
		bluefits[new operaImageVector(blue.getGuideWindow(extension),1)] = imageMedian;
		bluefits *= BPercent / 100.0;
        if (verbose)
			cout << "wircolorcomposite: median blue: " << imageMedian << endl;        
		
		// copy each color into the cube
		out[1] = redfits;
		out[2] = greenfits;
		out[3] = bluefits;
		
		out.operaFITSImageSave();
		out.operaFITSImageClose();
		if (split) {
			operaFITSImage outr(output.substr(0,output.find_last_of("."))+"_r.fits", red.getXDimension(), red.getYDimension(), tfloat, cNone);
			operaFITSImage outg(output.substr(0,output.find_last_of("."))+"_g.fits", green.getXDimension(), green.getYDimension(), tfloat, cNone);
			operaFITSImage outb(output.substr(0,output.find_last_of("."))+"_b.fits", blue.getXDimension(), blue.getYDimension(), tfloat, cNone);
			//outr.operaFITSImageCopyHeader(&red);	
			//outg.operaFITSImageCopyHeader(&green);	
			//outb.operaFITSImageCopyHeader(&blue);	
			outr = redfits;
			outg = greenfits;
			outb = bluefits;
			outr.operaFITSImageSave();
			outg.operaFITSImageSave();
			outb.operaFITSImageSave();
			outr.operaFITSImageClose();
			outg.operaFITSImageClose();
			outb.operaFITSImageClose();
		}
		red.operaFITSImageClose();
		green.operaFITSImageClose();
		blue.operaFITSImageClose();
	}
	catch (operaException e) {
		cerr << "wircolorcomposite: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "wircolorcomposite: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

