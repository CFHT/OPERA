
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaNICIImageTest 
 Version: 1.0
 Description: Perform various tests on an NICI image.
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

/*! \file operaNICIImageTest.cpp */

using namespace std;

/*! 
 * operaNICIImageTest
 * \author Edr martioli
 * \brief Perform various tests on an NICI image.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: operaNICIImageTest -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	   
    string inputImage;
    string outputImage1,outputImage2;
    
	struct option longopts[] = {
		{"inputImage",			1, NULL, 'i'},
		{"outputImage1",			1, NULL, '1'},
		{"outputImage2",			1, NULL, '2'},
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:1:2:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':
					inputImage = optarg;
					break;
				case '1':
					outputImage1 = optarg;
					break;
				case '2':
					outputImage2 = optarg;
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
			throw operaException("operaNICIImageTest: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}

		getFITSImageInformation(inputImage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
        cout << "operaNICIImageTest: " << inputImage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
		operaMultiExtensionFITSImage NICI(inputImage, tfloat, READONLY, cNone, false);
		cout << "NICI[1][512][512] = " << NICI[1][512][512] << endl;
		cout << "NICI[2][512][512] = " << NICI[2][512][512] << endl;
        
		if (!outputImage1.empty()) {
			operaMultiExtensionFITSImage *output1 = new operaMultiExtensionFITSImage(outputImage1, NICI.getnaxis1(), NICI.getnaxis2(),1, tfloat, cNone);
			
			// Note: These two loops provide identical results....
			
			//float *outputData = (float *)output->getpixels();
			
			//for(unsigned pixIndex=0;pixIndex<npixels;pixIndex++) {
			//	outputData[pixIndex] = inputData[pixIndex];
			//}
			for(unsigned y=0;y<NICI.getnaxis2();y++) {
				for(unsigned x=0;x<NICI.getnaxis1();x++) {
					(*output1)[1][y][x] = NICI[1][y][x];
				}
			}

			output1->operaMultiExtensionFITSImageCopyHeader(&NICI);
			output1->operaFITSImageSave();
			output1->operaFITSImageClose();
		}
        
		if (!outputImage2.empty()) {
			operaFITSImage output2(outputImage2, NICI.getnaxis1(), NICI.getnaxis2(), tfloat, cNone);
			
			// Note: These two loops provide identical results....
			
			//float *outputData = (float *)output->getpixels();
			
			//for(unsigned pixIndex=0;pixIndex<npixels;pixIndex++) {
			//	outputData[pixIndex] = inputData[pixIndex];
			//}
			for(unsigned y=0;y<NICI.getnaxis2();y++) {
				for(unsigned x=0;x<NICI.getnaxis1();x++) {
					output2[y][x] = NICI[1][y][x];
				}
			}
			output2.operaFITSImageCopyHeader(&NICI[1]);
			output2.operaFITSImageSave();
			output2.operaFITSImageClose();
		}
        
        NICI.operaFITSImageClose();
    } 
	catch (operaException e) {
		cerr << "operaNICIImageTest: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaNICIImageTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

