
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: sitelletest.cpp
 Version: 1.0
 Description: Perform various tests on sitelle images.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA 
 Date: Dec/2012
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2012  Opera Pipeline team, Canada France Hawaii Telescope
 
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
#include "libraries/operaStats.h"						// for operaArrayMedianQuick

/*! \file sitelletest.cpp */

using namespace std;

/*! 
 * sitelletest
 * \author Foug Teeple
 * \brief Perform various tests on sitelle images.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: sitelletest [--first=[1]] --last=[1-144]] -[pdvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		
		{"first",			optional_argument, NULL, 'f'},
		{"last",			optional_argument, NULL, 'l'},
		
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	unsigned first = 1;;
	unsigned last = 1;//144;

	try  {
		while ((opt = getopt_long(argc, argv, "f:l:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'f':
					first = atoi(optarg);
					break;
				case 'l':
					last = atoi(optarg);
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

		for (unsigned index=1; index<=last; index++) {
#ifdef __APPLE__
			string basepath = "/data/";
#else
			string basepath = "/data/uhane5/opera/";
#endif
			string fitsimage = basepath + "sitelle/test"+itos(index)+"o.fits";
			string outimage = "out"+itos(index)+"o.fits";
			string outasimage = "outas"+itos(index)+"o.fits";
			
            cout << "sitelletest: IMAGE " << fitsimage << endl;
            getFITSImageInformation(fitsimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "sitelletest: fitsimage " << fitsimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
            
			operaFITSImage infitsimage(fitsimage, tfloat, READONLY, false);
			operaFITSImage outfitsimage(outimage, infitsimage.getXDimension(), infitsimage.getYDimension(), tfloat, cNone, false);
			
			cout << "sitelletest: Copying infitsimage header to outfitsimage..." << endl;
			outfitsimage.operaFITSImageCopyHeader(&infitsimage);
			cout << "sitelletest: Copying infitsimage to " << outimage << endl;
			outfitsimage = infitsimage;		// copy entire image pixels
			
			float median = operaArrayMedianQuick(infitsimage.getnpixels(), (float *)infitsimage.getpixels());
            cout << "sitelletest: fitsimage " << fitsimage << " median= " << median << endl;
			
			cout << "sitelletest: Saving fitsimage.fits..." << endl;
			outfitsimage.operaFITSImageSave();
			
			cout << "sitelletest: Saving as " << outasimage << endl;
			outfitsimage.operaFITSImageSaveAs(outasimage);
			
			cout << "sitelletest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			cout << "sitelletest: fitsimage test complete." << endl;
		}
	}
	catch (operaException e) {
		cerr << "sitelletest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "sitelletest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

