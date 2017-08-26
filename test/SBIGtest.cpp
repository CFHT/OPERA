
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: SBIGtest.cpp
 Version: 1.0
 Description: Perform various tests on the operaFITSImage class for SBIG images.
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

/*! \file SBIGtest.cpp */

using namespace std;

/*! 
 * SBIGtest
 * \author Doug Teeple
 * \brief Perform various tests on the operaFITSImage class for SBIG images.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: SBIGtest -[dvth] <filelist>\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	string inputfilename;
	string outputfilename;
	unsigned test = 1;
	
	struct option longopts[] = {
		
		{"input",			1, NULL, 'i'},
		{"output",			1, NULL, 'o'},
		{"test",			1, NULL, 'e'},
		
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:o:e:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':
					inputfilename = optarg;
					break;
				case 'o':
					outputfilename = optarg;
					break;
				case 'e':
					test = atoi(optarg);
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
		

		if (test == 1) {
			unsigned XDimension, YDimension, ZDimension, Extensions;
			edatatype Datatype;
			long Npixels;
			
#ifdef __APPLE__
			string basepath = "/data/";
#else
			string basepath = "/data/uhane5/opera/";
#endif
			string fitsimage = basepath + "SBIG/4slice/InFocusFlat007.FIT";
			
            cout << "SBIGtest: SIMPLE FITS IMAGE" << endl;
            getFITSImageInformation(fitsimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "SBIGtest: fitsimage " << fitsimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
            
			operaFITSImage infitsimage(fitsimage, tfloat, READONLY, false);
			operaFITSImage outfitsimage("fitsimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), tfloat, cNone, false);
			operaFITSImage outfitsimageushort("fitsimagetushort.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), tushort, cNone, false);
			
			cout << "SBIGtest: Copying infitsimage to outfitsimage..." << endl;
			for (unsigned y=0; y<infitsimage.getYDimension(); y++) {
				for (unsigned x=0; x<infitsimage.getXDimension(); x++) {
					outfitsimage.setpixel(infitsimage[y][x], x, y);
				}
			}
			cout << "SBIGtest: Rotating outfitsimage..." << endl;
			outfitsimage.rotate90();		// rotate 90 degrees
			
			for (unsigned y=0; y<infitsimage.getYDimension(); y++) {
				for (unsigned x=0; x<infitsimage.getXDimension(); x++) {
					outfitsimageushort.setpixel((unsigned short)(infitsimage[y][x]), x, y);
				}
			}
			cout << "SBIGtest: Rotating fitsimagetushort..." << endl;
			outfitsimageushort.rotate90();		//  rotate 90 degrees
			
			cout << "SBIGtest: Saving fitsimage.fits..." << endl;
			outfitsimage.operaFITSImageSave();
			
			cout << "SBIGtest: Saving fitsimagetushort.fits..." << endl;
			outfitsimageushort.operaFITSImageSave();
			
			cout << "SBIGtest: Saving as fitsimagesaveas.fits..." << endl;
			outfitsimage.operaFITSImageSaveAs("fitsimagesaveas.fits");
			
			cout << "SBIGtest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			outfitsimageushort.operaFITSImageClose();
			cout << "SBIGtest: fitsimage test complete." << endl;
		}
		if (test == 2) {
			while (optind < argc) {
				string imagefilename = argv[optind++];
				cout << "SBIGtest: SBIG IMAGE " << imagefilename << endl;
				
				string basefilename = imagefilename.substr(imagefilename.find_last_of("/")+1);
				string directory = imagefilename.substr(0, imagefilename.find_last_of("/")+1);
				string newfilename = directory + basefilename.substr(0, basefilename.find(".FIT")) + ".fits";
				
				operaFITSImage infitsimage(imagefilename, tfloat, READONLY, false);
				operaFITSImage outfitsimageushort(newfilename, infitsimage.getXDimension(), infitsimage.getYDimension(), tushort, cNone, false);
				
				outfitsimageushort.operaFITSImageCopyHeader(&infitsimage);

				for (unsigned y=0; y<infitsimage.getYDimension(); y++) {
					for (unsigned x=0; x<infitsimage.getXDimension(); x++) {
						outfitsimageushort.setpixel((unsigned short)(infitsimage[y][x]), x, y);
					}
				}
				cout << "SBIGtest: Rotating..." << endl;
				outfitsimageushort.rotate90();		//  rotate 90 degrees
				
				cout << "SBIGtest: Saving " << newfilename << endl;
				outfitsimageushort.operaFITSImageSave();
				
				cout << "SBIGtest: Closing file(s)..." << endl;
				infitsimage.operaFITSImageClose();
				outfitsimageushort.operaFITSImageClose();
				cout << "SBIGtest: " << imagefilename << " complete." << endl;
			}
		}
	}
	catch (operaException e) {
		cerr << "SBIGtest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "SBIGtest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

