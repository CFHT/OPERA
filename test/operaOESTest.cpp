/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaOESTest
 Version: 1.0
 Description: Perform various tests on an OES Images.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA 
 Date: Jan/2011
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

#include <string.h>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSCube.h"
#include "libraries/operaLibCommon.h"

/*! \file operaOESTest.cpp */

using namespace std;

/*! 
 * operaOESTest
 * \author Doug Teeple / Megan Tannock
 * \brief Perform various tests on the operaCube classes.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	cerr << " Usage: operaOESTest -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
    string inFileName = "/data/oes/c201307040002.fit";   
    string inFileName2 = "/data/oes/c201307040003.fit";
    string inFileName3 = "/data/oes/c201307040004.fit"; 
    string inFileName4 = "/data/oes/c201307040014.fit";
	
	struct option longopts[] = {
		{"inFileName",			1, NULL, 'i'},
		{"inFileName2",			1, NULL, '2'},
		{"inFileName3",			1, NULL, '3'},
		{"inFileName4",			1, NULL, '4'},
		
		{"plot",				optional_argument, NULL, 'p'},
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
        while ((opt = getopt_long(argc, argv, "i:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
                    
				case 'i':
					inFileName = optarg;
					break;
				case '2':
					inFileName2 = optarg;
					break;
				case '3':
					inFileName3 = optarg;
					break;
				case '4':
					inFileName4 = optarg;
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
        
		if (inFileName.empty()) {
			throw operaException("operaOESTest: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
        // Get some information about the image.
        // getFITSImageInformation assigns values from inFileName in to variables XDimension, YDimension,
        // ZDimension, Extensions, Datatype. Will re-assign if called again with different file.
		getFITSImageInformation(inFileName, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels);
        cout << "operaOESTest: IMAGE: " << inFileName << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
        
		cout << "operaOESTest: Calling constructors" << endl;
		operaFITSImage inImage(inFileName, tfloat, READONLY, cNone);
		operaFITSImage outImage("outImage.fits", YDimension, XDimension, tfloat, cNone);
		
		cout << "operaOESTest: Copying headers from inImage to outImage..." << endl;
		outImage.operaFITSImageCopyHeader(&inImage);
		cout << "operaOESTest: Copying pixels from inImage to outImage..." << endl;
		//outImage = inImage;
		outImage.transpose(inImage);
		
		// Testing operaFITSCubeSave
		cout << "operaOESTest: Saving outImage as outImage.fits..." << endl;
		outImage.operaFITSImageSave(); 
		
		// Print some pixels to test they were copied correctly... index by [y][x]
		cout << "operaOESTest: inImage[1632][1464]=" << inImage[1632][1464] << " outImage[1464][1632]=" << outImage[1464][1632] << endl;
		inImage.operaFITSImageClose();
		outImage.operaFITSImageClose();
	}
	catch (operaException e) {
		cerr << "operaOESTest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaOESTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
    
    cout << "operaOESTest: TEST COMPLETE" << endl;
}

