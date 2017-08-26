
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: basicFITSImageTest.cpp
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

#include <getopt.h>
#include <iostream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSCube.h"
#include "libraries/operaMultiExtensionFITSImage.h"
#include "libraries/operaMultiExtensionFITSCube.h"

/*! \file basicFITSImageTest.cpp */

using namespace std;

/*! 
 * basicFITSImageTest
 * \author Megan Tannock
 * \brief Perform various tests on the operaFITSImage class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: basicFITSImageTest -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt, step=-1;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		
		{"step",			optional_argument, NULL, 's'},
		
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "s:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 's':
					step = atoi(optarg);
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
		
#ifdef __APPLE__
		string basepath = "/data/";
        string fitsimage = basepath + "espadons/11AQ14-Jul08/1315157o.fits";
#else
		string basepath = "/data/uhane5/opera/";
        string fitsimage = basepath + "11AQ14-Jul08/1315157o.fits";
#endif
		string mefimage = basepath + "/wircam/12AQ03-Mar01/1529925o.fits";
		string mefcubeimage = basepath + "/wircam/12AQ03-Mar01/1530725d.fits";
		string cubeimage = basepath + "/Cube/Cube.fits";

		/*
		 * Step 1 : Simple FITSImage
		 */
		if (step < 0 || step == 1) {

            cout << "basicFITSImageTest: SIMPLE FITS IMAGE" << endl;
            getFITSImageInformation(fitsimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "basicFITSImageTest: fitsimage " << fitsimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
            
			operaFITSImage infitsimage(fitsimage, tfloat, READONLY, false);
			operaFITSImage outfitsimage("fitsimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), tfloat, cNone, false);
			operaFITSImage outfitsimageushort("fitsimagetushort.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), tushort, cNone, false);
			
			cout << "basicFITSImageTest: Copying infitsimage header to outfitsimage..." << endl;
			outfitsimage.operaFITSImageCopyHeader(&infitsimage);
			cout << "basicFITSImageTest: Copying infitsimage to outfitsimage..." << endl;
			outfitsimage = infitsimage;		// copy entire image pixels
			
			for (unsigned y=0; y<infitsimage.getYDimension(); y++) {
				for (unsigned x=0; x<infitsimage.getXDimension(); x++) {
					outfitsimageushort.setpixel((unsigned short)infitsimage[y][x], x, y);
				}
			}
			cout << "basicFITSImageTest: Saving fitsimage.fits..." << endl;
			outfitsimage.operaFITSImageSave();
			
			cout << "basicFITSImageTest: Saving fitsimagetushort.fits..." << endl;
			outfitsimageushort.operaFITSImageSave();
			
			cout << "basicFITSImageTest: Saving as fitsimagesaveas.fits..." << endl;
			outfitsimage.operaFITSImageSaveAs("fitsimagesaveas.fits");
			
			cout << "basicFITSImageTest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			outfitsimageushort.operaFITSImageClose();
			cout << "basicFITSImageTest: fitsimage test complete." << endl;
		}
		
		/*
		 * Step 2 : MEF FITSImage
		 */
		if (step < 0 || step == 2) {
            cout << "basicFITSImageTest: MEF FITS IMAGE" << endl;
            getFITSImageInformation(mefimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "basicFITSImageTest: mefimage " << mefimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;

			operaMultiExtensionFITSImage infitsimage(mefimage, tfloat, READONLY, cNone, false);
			operaMultiExtensionFITSImage outfitsimage("mefimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), infitsimage.getNExtensions(), tfloat, cNone, false);
			
			cout << "basicFITSImageTest: Copying mefimage header to outfitsimage..." << endl;
			outfitsimage.operaMultiExtensionFITSImageCopyHeader(&infitsimage);
			cout << "basicFITSImageTest: Copying mefimage to outfitsimage..." << endl;
			outfitsimage = infitsimage;		// copy entire image pixels
			
			cout << "basicFITSImageTest: Saving mefimage.fits..." << endl;
			outfitsimage.operaMultiExtensionFITSImageSave();
			 
			cout << "basicFITSImageTest: Saving as mefsimagesaveas.fits..." << endl;
			outfitsimage.operaMultiExtensionFITSImageSaveAs("mefimagesaveas.fits");
			
			cout << "basicFITSImageTest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			cout << "basicFITSImageTest: mefimage test complete." << endl;
		}
        
        /*
         * Step 3 : Simple MEF
         */
        if (step < 0 || step == 3) {
            cout << "basicFITSImageTest: SIMPLE MEF FITS IMAGE" << endl;            
			operaMultiExtensionFITSImage infitsimage(mefimage, tfloat, READONLY, cNone, false);
			operaMultiExtensionFITSImage outfitsimage("simplemefimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), infitsimage.getNExtensions(), tfloat, cNone, false);
			
			cout << "basicFITSImageTest: Saving mefimage.fits..." << endl;
			outfitsimage.operaMultiExtensionFITSImageSave();
			cout << "basicFITSImageTest: Saving as mefsimagesaveas.fits..." << endl;
			outfitsimage.operaMultiExtensionFITSImageSaveAs("simplemefimagesaveas.fits");
			
			cout << "basicFITSImageTest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			cout << "basicFITSImageTest: simple mefimage test complete." << endl;
        }
		
		/*
		 * Step 4 : Cube FITSImage
		 */ 
		if (step < 0 || step == 4) {
            cout << "basicFITSImageTest: CUBE FITS IMAGE" << endl;
            getFITSImageInformation(cubeimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "basicFITSImageTest: cubeimage " << cubeimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;

			operaFITSCube infitsimage(cubeimage, tfloat, READONLY, cNone, false);
			operaFITSCube outfitsimage("cubeimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), infitsimage.getZDimension(), tfloat, cNone, false);
		 	
	 		cout << "basicFITSImageTest: Copying cubeimage header to outfitsimage..." << endl;
			outfitsimage.operaFITSImageCopyHeader(&infitsimage);
			cout << "basicFITSImageTest: Copying cubeimage to outfitsimage..." << endl;
			outfitsimage = infitsimage;		// copy entire image pixels
			
			cout << "basicFITSImageTest: Saving cubeimage.fits..." << endl;
			outfitsimage.operaFITSCubeSave();
			
			cout << "basicFITSImageTest: Saving as cubeimagesaveas.fits..." << endl;
			outfitsimage.operaFITSCubeSaveAs("cubeimagesaveas.fits");
			
			cout << "basicFITSImageTest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			cout << "basicFITSImageTest: cubeimage test complete." << endl;
		}
		
		/*
		 * Step 5 : MEFCubeImage
		 */
		if (step < 0 || step == 5) {
            cout << "basicFITSImageTest: MEF CUBE IMAGE" << endl;
            getFITSImageInformation(mefcubeimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "basicFITSImageTest: mefcubeimage " << mefcubeimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
            
			operaMultiExtensionFITSCube infitsimage(mefcubeimage, tfloat, READONLY, cNone, false);
			operaMultiExtensionFITSCube outfitsimage("mefcubeimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), infitsimage.getZDimension(), infitsimage.getNExtensions(), tfloat, cNone, false);
			
			cout << "basicFITSImageTest: Copying operaMultiExtensionFITSCube header to outfitsimage..." << endl;
			outfitsimage.operaMultiExtensionFITSCubeCopyHeader(&infitsimage);
			cout << "basicFITSImageTest: Copying operaMultiExtensionFITSCube to outfitsimage..." << endl;
			outfitsimage = infitsimage;		// copy entire image pixels
			
			cout << "basicFITSImageTest: Saving mefcubeimage.fits..." << endl;
			outfitsimage.operaMultiExtensionFITSCubeSave();
			
			cout << "basicFITSImageTest: Saving as mefcubeimagesaveas.fits..." << endl;
			outfitsimage.operaMultiExtensionFITSCubeSaveAs("mefcubeimagesaveas.fits");
			
			cout << "basicFITSImageTest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			cout << "basicFITSImageTest: mefcubeimage test complete." << endl;
		}

		/*
		 * Step 6 : Test Cloning
		 */
		if (step < 0 || step == 6) {
			
            cout << "basicFITSImageTest: Clone FITS IMAGE" << endl;
            getFITSImageInformation(fitsimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "basicFITSImageTest: fitsimage " << fitsimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
            
			operaFITSImage infitsimage(fitsimage, tfloat, READONLY, false);
			operaFITSImage outfitsimage("fitsimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), tfloat, cNone, false);
			operaFITSImage outfitsimageushort("fitsimagetushort.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), tushort, cNone, false);
			
			cout << "basicFITSImageTest: Copying infitsimage header to outfitsimage..." << endl;
			outfitsimage.operaFITSImageCopyHeader(&infitsimage);
			cout << "basicFITSImageTest: Copying infitsimage header to outfitsimageushort..." << endl;
			outfitsimageushort.operaFITSImageCopyHeader(&infitsimage);
			
			cout << "basicFITSImageTest: Copying infitsimage to outfitsimage..." << endl;
			outfitsimage = infitsimage;		// copy entire image pixels
			
			for (unsigned y=0; y<infitsimage.getYDimension(); y++) {
				for (unsigned x=0; x<infitsimage.getXDimension(); x++) {
					outfitsimageushort.setpixel((unsigned short)infitsimage[y][x], x, y);
				}
			}
			
			unsigned x = 1024;
			unsigned y = 0;	// must be 0 for this test as the clone has only one dimension...
			
			cout << "basicFITSImageTest: Cloning infitsimage to floatclone..." << endl;
			float *floatclone = outfitsimage.operaFITSImageClonePixels();
			cout << "x= " << x << " y= " << y << " original= " << outfitsimage[y][x] <<  " address= " << floatclone << " clone= " << floatclone[x] << endl;

			x = 1036;
			cout << "basicFITSImageTest: Cloning outfitsimageushort to ushortclone..." << endl;
			unsigned short *ushortclone = outfitsimageushort.operaFITSImageClonePixelsUSHORT();
			cout << "x= " << x << " y= " << y << " original= "  << outfitsimageushort.getpixelUSHORT(x, y) <<  " address= " << ushortclone << " clone= " << ushortclone[x] << endl;

			cout << "basicFITSImageTest: Saving fitsimage.fits..." << endl;
			outfitsimage.operaFITSImageSave();
			
			cout << "basicFITSImageTest: Saving fitsimagetushort.fits..." << endl;
			outfitsimageushort.operaFITSImageSave();
			
			cout << "basicFITSImageTest: Saving as fitsimagesaveas.fits..." << endl;
			outfitsimage.operaFITSImageSaveAs("fitsimagesaveas.fits");
			
			cout << "basicFITSImageTest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			outfitsimageushort.operaFITSImageClose();
			cout << "basicFITSImageTest: fitsimage test complete." << endl;
		}
		
	}
	catch (operaException e) {
		cerr << "basicFITSImageTest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "basicFITSImageTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

