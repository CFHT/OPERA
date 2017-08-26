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

#define MAXIMAGES 100

/*! \file operaFITSImageTest.cpp */

using namespace std;

/*! 
 * operaFITSImageTest
 * \author Doug Teeple
 * \brief Perform various tests on the operaFITSImage class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{
	int opt;
	string image;
	string outputs[MAXIMAGES];
	string images[MAXIMAGES];
	operaFITSImage *Ins[MAXIMAGES];
	string listofimages;
	unsigned imageIndex = 0;
	unsigned outImageIndex = 0;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"images",			1, NULL, 'i'},	// input images
		{"list",			1, NULL, 'l'},	// list of input flats		
		{"output",			1, NULL, 'o'},	// output fits file
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:o:vdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// image
					images[imageIndex++] = optarg;
					break;
				case 'l':		// list of images
					listofimages = optarg;
					break;					
				case 'o':		// output
					outputs[outImageIndex++] = optarg;
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
		if (!listofimages.empty()) {
			ifstream flist(listofimages.c_str());		
			if (flist.is_open())
			{
				while (flist.good()) {
					getline (flist,images[imageIndex++]);
					if (images[imageIndex-1].size() == 0 || images[imageIndex-1][0] == '#')
						imageIndex--;					
				}	
				flist.close();
			}
		}
		if (debug) {
			for(unsigned i=0;i<imageIndex;i++)
				cout << images[i] << " " << i << endl;
		}
		if (imageIndex == 0) {
			throw operaException("operaFITSImageTest: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (outImageIndex == 0) {
			throw operaException("operaFITSImageTest: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}

		// open the input file
		unsigned i;
		for (i=0; i<imageIndex; i++) {
			cerr << "Opening " << images[i] << '\n';
			Ins[i] = new operaFITSImage(images[i], tfloat, READONLY);
			cerr << "DETECTOR= " << Ins[i]->operaFITSGetHeaderValue("DETECTOR") << '\n';
			cerr << "AMPLIST= " << Ins[i]->operaFITSGetHeaderValue("AMPLIST") << '\n';
			cerr << "INSTMODE= " << Ins[i]->operaFITSGetHeaderValue("INSTMODE") << '\n';
			cerr << "EREADSPD= " << Ins[i]->operaFITSGetHeaderValue("EREADSPD") << '\n';
			cerr << "CMMTSEQ= " << Ins[i]->operaFITSGetHeaderValue("CMMTSEQ") << '\n';
			if (i == (imageIndex-1)) {
				// clone to the output file
				cerr << "Creating " << outputs[0] << ' ' << Ins[i]->getnaxis1() << 'x' << Ins[i]->getnaxis2() << '\n';
				operaFITSImage *Out = new operaFITSImage(outputs[0], Ins[i]->getnaxis1(), Ins[i]->getnaxis2(), tushort, 0);
				cerr << "Copying header from " << images[i] << " to "<< outputs[0] << '\n';
				Out->operaFITSImageCopyHeader(Ins[i]);
				Out->operaFITSImageSetData(Ins[i]->operaFITSImageClonePixelsUSHORT());
				cerr << "DETECTOR= " << Out->operaFITSGetHeaderValue("DETECTOR") << '\n';
				cerr << "AMPLIST= " << Out->operaFITSGetHeaderValue("AMPLIST") << '\n';
				cerr << "INSTMODE= " << Out->operaFITSGetHeaderValue("INSTMODE") << '\n';
				cerr << "EREADSPD= " << Out->operaFITSGetHeaderValue("EREADSPD") << '\n';
				cerr << "CMMTSEQ= " << Out->operaFITSGetHeaderValue("CMMTSEQ") << '\n';
				cerr << "Saving " << outputs[0] << '\n';
				Out->operaFITSImageSave();
				cerr << "Closing " << outputs[0] << '\n';
				Out->operaFITSImageClose();
				// convert input to float
				cerr << "Converting " << images[i] << " to float " << '\n';
				Ins[i]->operaFITSImageConvertImage(tfloat);
				if (outImageIndex == 2) {
					cerr << "Saving " << images[i] << " as " << outputs[1] << '\n';
					Ins[i]->operaFITSImageSaveAs(outputs[1]);
				}
				cerr << "Closing " << images[i] << '\n';
				Ins[i]->operaFITSImageClose();
				cerr << "Done " << '\n';
			}
		}
	}
	catch (operaException e) {
		cerr << "operaTestFITSImage: " << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaTestFITSImage: " << operaStrError(errno) << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cerr << " Usage: operaFITSImageTest [--image=<filename>]+ --output=<file name> -[dvth]\n";
}	
