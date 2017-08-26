
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaNIFSImageTest
 Version: 1.0
 Description: Perform various tests on an NIFS image.
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

#include <getopt.h>
#include <iostream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMultiExtensionFITSImage.h"

/*! \file operaNIFSImageTest.cpp */

using namespace std;

/*!
 * operaNIFSImageTest
 * \author Eder Marioli
 * \brief Perform various tests on  an NIFS image.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: operaNIFSImageTest -[dvth]\n";
}

int main(int argc, char *argv[])
{
	int opt;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
    
    string inputImage;
    string outputImage;
    
	struct option longopts[] = {
		{"inputImage",			1, NULL, 'i'},
		{"outputImage",			1, NULL, 'o'},
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:o:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':
					inputImage = optarg;
					break;
				case 'o':
					outputImage = optarg;
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
        
		if (inputImage.empty()) {
			throw operaException("operaNIFSImageTest: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		
        getFITSImageInformation(inputImage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels);
        cout << "operaNIFSImageTest: " << inputImage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
		
        operaMultiExtensionFITSImage geminiNIFS(inputImage, tfloat, false);
        unsigned npixels = geminiNIFS.getnpixels();
        
		cout << "npixels="<<npixels<<" geminiNIFS[1][1024][1024] = " << geminiNIFS[1][1024][1024] << endl;
        
		if (!outputImage.empty()) {
			operaMultiExtensionFITSImage *output = new operaMultiExtensionFITSImage(outputImage, XDimension, YDimension, ZDimension, tfloat, cNone);
			
			for(unsigned y=0;y<geminiNIFS.getnaxis2();y++) {
				for(unsigned x=0;x<geminiNIFS.getnaxis1();x++) {
					(*output)[1][y][x] = geminiNIFS[1][y][x];
				}
			}
			//*output = geminiNIFS; // is faster
			
			output->operaMultiExtensionFITSImageCopyHeader(&geminiNIFS);
			output->operaFITSImageSave();
			output->operaFITSImageClose();			
		}
        geminiNIFS.operaFITSImageClose();
    }
	catch (operaException e) {
		cerr << "operaNIFSImageTest: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaNIFSImageTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

