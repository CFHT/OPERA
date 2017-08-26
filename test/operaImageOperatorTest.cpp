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
#include "operaImageOperatorTest.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaImageVector.h"
#include "libraries/operaFITSSubImage.h"

#include "libraries/operaImage.h"

/*! \file operaFITSImageTest.cpp */

using namespace std;

/*! 
 * operaImageOperatortTest
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
			throw operaException("operaImageOperatorTest: need an output file name ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (bias.empty()) {
			throw operaException("operaImageOperatorTest: need a bias file name ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (flat1.empty()) {
			throw operaException("operaImageOperatorTest: need a flat1 file name ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (flat2.empty()) {
			throw operaException("operaImageOperatorTest: need a flat2 file name ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (badpix.empty()) {
			throw operaException("operaImageOperatorTest: need a badpix file name ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (verbose) {
			cout << " bias = " << bias << " flat1 = " << flat1 << " flat2 = " << flat2 << " badpix = " << badpix << endl;
		}
		
		operaFITSImage biasImage(bias, tfloat, READONLY);
		operaFITSImage flat1Image(flat1, tfloat, READONLY);
		operaFITSImage flat2Image(flat2, tfloat, READONLY);
		operaFITSImage badpixImage(badpix, tfloat, READONLY);
		operaFITSImage outputImage(output, biasImage.getnaxis1(), biasImage.getnaxis2(), tfloat, 0);
		outputImage.operaFITSImageCopyHeader(&biasImage);

		outputImage = 0.0;							// zero the whole output image
		flat1Image -= biasImage;					// remove the bias from the flats
		flat2Image -= biasImage;					// remove the bias from the flats
		flat1Image *= flat1Image < 35000.0;			// remove saturated pixels
		flat2Image *= flat2Image < 35000.0;			// remove saturated pixels
		outputImage = flat1Image - flat2Image;		// subtract the flats
		outputImage /= flat1Image;					// divide by flat1
		outputImage *= outputImage > 0.0;			// remove any values less than zero
		outputImage *= badpixImage;					// zero the badpixels
		outputImage += 100.0;						// add in a bias
		outputImage[where(badpixImage == 1.0)] = 0.0;// mask bad pixels
		// print some pixel values from the center of the image
		for (unsigned y=2048; y<2056; y++) {
			for (unsigned x=1024; x<1032; x++) {
				cout << setw(8) << setprecision(4) << (float)outputImage[y][x] << ' ';	// output the pixel value
			}
			cout << endl;
		}
		cout << endl;
		float *rowpointer = outputImage[1204];
		cout << "row 1204 pixel 5 = " << rowpointer[5] << endl;
		float rowmean = 0.0;
		for (unsigned col=0; col<outputImage.getXDimension(); col++) {
			rowmean += rowpointer[col];
		}
		rowmean /= outputImage.getXDimension();
		cout << "row 1204 mean = " << rowmean << endl;

		// Note also that the class is indexed [y][x]
		for (unsigned y=2048; y<2056; y++) {
			for (unsigned x=1024; x<1032; x++) {
				outputImage[y][x] = outputImage[x][y];							// swap the pixels around
				cout << setw(8) << setprecision(4) << (float)outputImage[y][x] << ' ';	// output the pixel value
				outputImage[y][x] = 0.0;										// zero the pixel value
			}
			cout << endl;
		}

		outputImage.operaFITSImageSave();
		outputImage.operaFITSImageClose();

	}
	catch (operaException e) {
		cerr << "operaImageOperatortest: " << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaImageOperatortest: " << operaStrError(errno) << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cerr << " Usage: operaImageOperatortest --bias=<filename> --flat1=<filename> --flat2=<filename> --badpix=<filename> --output=<file name> -[dvth]\n";
}	
