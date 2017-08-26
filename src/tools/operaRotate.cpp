
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaRotate.cpp
 Version: 1.0
 Description: operaRotate.
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
#include "libraries/operaFITSCube.h"
#include "libraries/operaMultiExtensionFITSImage.h"
#include "libraries/operaMultiExtensionFITSCube.h"

/*! \file operaRotate.cpp */

using namespace std;

/*! 
 * operaRotate
 * \author Doug Teeple
 * \brief Rotate an image 90 degrees.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: operaRotate -[dvth] [--replace=<text>] <filelist>\n";
	cerr << " if --replace=<text> is given then that text is added to the basename of the image file\n";
	cerr << " i.e. --replace=_90_ 1000o.fits -> _90_1000o.fits\n";
	cerr << " if --replace is not given then the image file is overwritten\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	string inputfilename;
	string outputfilename;
	string replacetemplate;
	
	struct option longopts[] = {
		
		{"input",			1, NULL, 'i'},
		{"output",			1, NULL, 'o'},
		{"replace",			1, NULL, 'r'},
		
		{"plot",			optional_argument, NULL, 'p'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:o:r:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':
					inputfilename = optarg;
					break;
				case 'o':
					outputfilename = optarg;
					break;
				case 'r':
					replacetemplate = optarg;
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
		
		while (optind < argc) {
			string imagefilename = argv[optind++];
			cout << "operaRotate: " << imagefilename << endl;
			
			string basefilename = imagefilename.substr(imagefilename.find_last_of("/")+1);
			string directory = imagefilename.substr(0, imagefilename.find_last_of("/")+1);
			string newfilename = "/tmp/" + basefilename;
			if (!replacetemplate.empty()) {
				newfilename = directory + replacetemplate + basefilename;
			}
			operaFITSImage infitsimage(imagefilename, tfloat, READONLY, false);
			operaFITSImage outfitsimageushort(newfilename, infitsimage.getXDimension(), infitsimage.getYDimension(), tushort, cNone, false);
			
			outfitsimageushort.operaFITSImageCopyHeader(&infitsimage);
			
			for (unsigned y=0; y<infitsimage.getYDimension(); y++) {
				for (unsigned x=0; x<infitsimage.getXDimension(); x++) {
					outfitsimageushort.setpixel((unsigned short)(infitsimage[y][x]), x, y);
				}
			}
			cout << "operaRotate: Rotating..." << endl;
			outfitsimageushort.rotate90();		//  rotate 90 degrees
			
			cout << "operaRotate: Saving " << newfilename << endl;
			outfitsimageushort.operaFITSImageSave();
			
			cout << "operaRotate: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimageushort.operaFITSImageClose();
			if (replacetemplate.empty()) {
				remove(imagefilename.c_str());
				rename(newfilename.c_str(), imagefilename.c_str());
			}
			cout << "operaRotate: " << imagefilename << " complete." << endl;
		}
	}
	catch (operaException e) {
		cerr << "operaRotate: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaRotate: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

