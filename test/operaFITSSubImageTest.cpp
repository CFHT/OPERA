/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFITSImageTest
 Version: 1.0
 Description: Perform various tests on the operaFITSImage class.
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
#include <libgen.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "operaFITSImageTest.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"

/*! \file operaFITSImageTest.cpp */

using namespace std;

/*! 
 * operaFITSSubImageTest
 * \author Doug Teeple
 * \brief Perform various tests on the operaFITSImage class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

#define MAXIMAGES 1000

int main(int argc, char *argv[])
{
	int opt;
	string image;
	string output;
	string bias;
	string flat1;
	string flat2;
	string badpix;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"bias",			1, NULL, 'i'},	// bias image
		{"badpix",			1, NULL, 'x'},	// badpix image
		{"flat1",			1, NULL, '1'},	// flat image
		{"flat2",			1, NULL, '2'},	// flat image
		{"output",			1, NULL, 'o'},	//  output fits file
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:o:x:1:2:vdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// bias
					bias = optarg;
					break;
				case '1':		// flat1
					flat1 = optarg;
					break;					
				case '2':		// flat2
					flat2 = optarg;
					break;					
				case 'x':		// flat2
					badpix = optarg;
					break;					
				case 'o':		// output
					output = optarg;
					break;
					
				case 'v':
					verbose = true;
					break;
				case 'd':
					debug = true;
					break;
				case 't':
					trace = true; 
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
		
		/* 
		 * Read image path from input file and append to list of input images	
		 */ 
		if (output.empty()) {
			throw operaException("operaFITSSubImageTest: need an output file name ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (bias.empty()) {
			throw operaException("operaFITSSubImageTest: need a bias file name ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (flat1.empty()) {
			throw operaException("operaFITSSubImageTest: need a flat1 file name ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (flat2.empty()) {
			throw operaException("operaFITSSubImageTest: need a flat2 file name ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (badpix.empty()) {
			throw operaException("operaFITSSubImageTest: need a badpix file name ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (verbose) {
			cout << " bias = " << bias << " flat1 = " << flat1 << " flat2 = " << flat2 << " badpix = " << badpix << endl;
		} 
		operaFITSImage    biasFITSImage(bias, tfloat, READONLY);
		operaFITSImage    flat1FITSImage(flat1, tfloat, READONLY);
		operaFITSImage    flat2FITSImage(flat2, tfloat, READONLY);
		operaFITSImage    badpixFITSImage(badpix, tfloat, READONLY);
		operaFITSImage    outputFITSImage(output, biasFITSImage.getnaxis1(), biasFITSImage.getnaxis2(), tfloat, 0);
	
		outputFITSImage.operaFITSImageCopyHeader(&biasFITSImage);

		operaFITSSubImage biasSubImage(biasFITSImage, 0, 0, 1024, 2048);
		operaFITSSubImage flat1SubImage(flat1FITSImage, 0, 0, 1024, 2048);
		operaFITSSubImage flat2SubImage(flat2FITSImage, 0, 0, 1024, 2048);
		operaFITSSubImage badpixSubImage(badpixFITSImage, 0, 0, 1024, 2048);
		operaFITSSubImage outputSubImage(outputFITSImage, 0, 0, 1024, 2048);
		
		flat1SubImage -= biasSubImage;		// remove the bias from the flats
		flat2SubImage -= biasSubImage;		// remove the bias from the flats
		outputSubImage = flat1SubImage - flat2SubImage;		// subtract the flats
		outputSubImage *= badpixSubImage;	// zero the badpixels
		outputSubImage += 100.0;			// add in a bias
		
		cout << "The output image set from the subimage" << endl;
		for (unsigned y=0; y<outputSubImage.getny(); y++) {
			for (unsigned x=0; x<outputSubImage.getnx(); x++) {
				outputFITSImage[y][x] = outputSubImage[y][x];
				if (x >= 0 && x <= 7 && y >= 0 && y <= 7) {
					cout << setw(8) << setprecision(4) << outputSubImage[y][x] << ' ' ;	// print the image flux value
				}
			}
			if (y >= 0 && y <= 7) 
				cout << endl;
		}
		cout << endl;
		// transpose rows for cols, be sure to get the new matrix...
		outputSubImage.transpose();

		outputFITSImage.operaFITSImageSave();
		outputFITSImage.operaFITSImageClose();
		
		// here we are re-dimensioning the subimage as a vector, in two different ways...
		//operaFITSSubImage subImage1(flat1SubImage, where(outputSubImage > 200.0 && outputSubImage <= 400.0));
		//operaFITSSubImage *subImage2 = flat1SubImage[where(outputSubImage > 200.0 && outputSubImage <= 400.0)];
		//delete subImage2;
	}
	catch (operaException e) {
		cerr << "operaTestFITSImage: " << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaTestFITSImage: " << operaStrError(errno) << "\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cerr << " Usage: operaFITSSubImageTest --bias=<filename> --flat1=<filename> --flat2=<filename> --badpix=<filename> --output=<file name> -[dvth]\n";
}	
