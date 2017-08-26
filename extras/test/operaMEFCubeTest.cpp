/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMEFCubeTest
 Version: 1.0 
 Description: Perform various tests on the operaMEFFITSImage and MEFFITSCube classes.
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
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMultiExtensionFITSImage.h"
#include "libraries/operaMultiExtensionFITSCube.h"
#include "libraries/operaFITSCube.h"

/*! \file operaMEFCubeTest.cpp */

using namespace std;

/*! 
 * operaMEFCubeTest
 * \author Doug Teeple / Megan Tannock
 * \brief Perform various tests on the operaMEFFITSImage and MEFFITSCube classes.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: operaMEFCubeTest -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;
	string image;
	string output;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		
		{"plot",			optional_argument, NULL, 'v'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {					
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
#else
		string basepath = "/data/uhane5/opera/";
#endif
        string mefcubeimage = basepath + "/wircam/12AQ03-Mar01/1530166o.fits";
        
        /*
         * STEP ONE: Read entire image in to memory, basic tests
         */
        if (0) {
			// Get some information about the image...
			cout << "operaMEFCubeTest: Calling constructors" << endl;
			operaMultiExtensionFITSCube inmefcubeimage(mefcubeimage, tfloat, READONLY, cNone, false);
			operaMultiExtensionFITSCube outmefcubeimage("mefcubeimage.fits", inmefcubeimage.getXDimension(), inmefcubeimage.getYDimension(), inmefcubeimage.getZDimension(), inmefcubeimage.getNExtensions(), tfloat, cNone, false);
			
			cout << "operaMEFCubeTest: Testing if file is a MEF Cube" << endl;
			if (!inmefcubeimage.isMEFCube()) {
				throw operaException("operaMEFCubeTest: "+mefcubeimage+" is not a MEF Cube. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
			}
			
			getFITSImageInformation(mefcubeimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
			cout << "operaMEFCubeTest: " << mefcubeimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
			
            cout << "operaMEFCubeTest: STEP 1: Entire image in memory tests... " << endl;
            
            // Print some information about the images...
            cout << "operaMEFCubeTest: inmefcubeimage: XDimension = " << inmefcubeimage.getXDimension() << " YDimension = " << inmefcubeimage.getYDimension() << " ZDimension = " <<  inmefcubeimage.getZDimension() << " NExtensions = " << inmefcubeimage.getNExtensions() << " npixels = " << inmefcubeimage.getnpixels() << endl;  
          
            if (1) {
				cout << "operaMEFCubeTest: outmefcubeimage: XDimension = " << outmefcubeimage.getXDimension() << " YDimension = " << outmefcubeimage.getYDimension() << " ZDimension = " <<  outmefcubeimage.getZDimension() << " NExtensions = " << outmefcubeimage.getNExtensions() << " npixels = " << outmefcubeimage.getnpixels() << endl;
				// Copy headers and pixels from inmefcubeimage to outmefcubeimage
				cout << "\noperaMEFCubeTest: Copying inmefcubeimage headers to outmefcubeimage..." << endl;   // takes a long time...
				outmefcubeimage.operaMultiExtensionFITSCubeCopyHeader(&inmefcubeimage);
				cout << "operaMEFCubeTest: Copying inmefcubeimage pixels to outfitsimage..." << endl;
				outmefcubeimage = inmefcubeimage;		// copy image pixels
				cout << "operaMEFFITSImageTest: Copied inmefcubeimage headers and pixels to outmefcubeimage..." << endl;

				// Print some pixels to test they were copied correctly... index by [extension][slice][y][x]
				cout << "operaMEFCubeTest: Printing pixels for extension 1 (first 4 slices)" << endl;
				for (unsigned slicecounter=1; slicecounter<=4; slicecounter++) {
					cout << "operaMEFCubeTest: outmefcubeimage[1][" << slicecounter << "][1024][1024]=" << outmefcubeimage[1][slicecounter][1024][1024] << " inmefcubeimage[1][" << slicecounter << "][1024][1024]=" << inmefcubeimage[1][slicecounter][1024][1024] << endl;
				}           
				// Testing basic operaMultiExtensionFITSCubeSave, operaMultiExtensionFITSCubeSaveAs
				cout << "operaMEFCubeTest: Saving mefcubeimage.fits... (operaMultiExtensionFITSCubeSave)" << endl;
				outmefcubeimage.operaMultiExtensionFITSCubeSave();
				cout << "operaMEFCubeTest: Saving as mefcubeimageSaveAs.fits... (operaMultiExtensionFITSCubeSaveAs)" << endl;
				outmefcubeimage.operaMultiExtensionFITSCubeSaveAs("mefcubeimageSaveAs.fits");
				
				// Subtract the bias from the first slice of each extension
				cout << "\noperaMEFCubeTest: Subtracting bias (7000) from first slice of each extension... " << endl;
				for (unsigned ext=1; ext<=Extensions; ext++) {
					outmefcubeimage[ext][1] -= 7000; // extension ext, slice 1
					cout << "operaMEFCubeTest: outmefcubeimage[" << ext << "][1][1024][1024]=" << outmefcubeimage[ext][1][1024][1024] << " inmefcubeimage[" << ext << "][1][1024][1024]=" << inmefcubeimage[ext][1][1024][1024] << endl;            
				}
				
				cout << "operaMEFCubeTest: Saving outmefcubeimage as mefcubeimageSaveAsBiasRemoved.fits..." << endl;
				outmefcubeimage.operaMultiExtensionFITSCubeSaveAs("mefcubeimageSaveAsBiasRemoved.fits");
				
				cout << "\noperaMEFCubeTest: Copying slice 3 of first extension from inmefcubeimage to slice 2 of first extension of outmefcubeimage..." << endl;
				outmefcubeimage[1][2] = inmefcubeimage[1][3];	// copy one slice
				cout << "operaMEFCubeTest: slice switch outmefcubeimage[1][2][1024][1024]=" << outmefcubeimage[1][2][1024][1024] << " inmefcubeimage[1][3][1024][1024]=" << inmefcubeimage[1][3][1024][1024] << endl;
				cout << "operaMEFCubeTest: Saving outmefcubeimage as mefcubeimageSaveAsSlicesSwitched.fits..." << endl;
				outmefcubeimage.operaMultiExtensionFITSCubeSaveAs("mefcubeimageSaveAsSlicesSwitched.fits");
            }
			if (1) {
				cout << "\noperaMEFCubeTest: medianCollapse..." << endl;
				//inmefcubeimage.medianCollapseSingleExtension(1);
				inmefcubeimage.medianCollapse();
				cout << "operaMEFCubeTest: medianCollapse complete" << endl;
				cout << "operaMEFCubeTest: creating MEF" << endl;
				operaMultiExtensionFITSImage collapsed(inmefcubeimage, 1 ,true);
				cout << "operaMEFCubeTest: Saving as mefimageCollapse.fits..." << endl;
				collapsed.operaMultiExtensionFITSImageSaveAs("mefimageCollapse.fits");
				collapsed.operaFITSImageClose();				
			}
            if (1) {
				// REQUIRES TESTING:
				// downsize outmefcubeimage to a MEF FITS Image
				cout << "\noperaMEFCubeTest: downsizing extension 1 of outmefcubeimage to a FITS Cube Image..." << endl;
				operaFITSCube downSize(outmefcubeimage[1], false, true);   // is this correct?
				cout << "operaMEFCubeTest: Copying outmefcubeimage header to downSize..." << endl;
				// Note that "outmefcubeimage" is a MEF Cube but "downSize" is just a regular FITS Cube...
				// Print centre pixel of slice 1. For a FITS Cube: [slice][y][x]
				cout << "operaMEFCubeTest: downSize[1][1024][1024]=" << downSize[1][1024][1024] << endl;
				// Save downSize as MEF FITS image...
				cout << "operaMEFCubeTest: Saving downSize as outmefcubeimageSaveAsDownSize.fits... (operaFITSCubeSaveAs)" << endl;
				downSize.operaFITSImageClose();
				downSize.operaFITSCubeSaveAs("outmefcubeimageSaveAsDownSize.fits");
            }
            // Close the files
            cout << "\noperaMEFCubeTest: Closing files..." << endl;
            inmefcubeimage.operaFITSImageClose();
            outmefcubeimage.operaFITSImageClose();
        }
        
		/*  
		 * Step 2 extension save/saveas -- LAZY READ
		 */
        // REQUIRES TESTING
        if (1) {
            cout << "operaMEFCubeTest: STEP 2: Lazy read tests... " << endl;
            
            cout << "operaMEFCubeTest: Calling constructors" << endl;
            operaMultiExtensionFITSCube inmefcubeimage(mefcubeimage, tfloat, READONLY, cNone, true);
            operaMultiExtensionFITSCube outmefcubeimage("lazymefcubeimage.fits", inmefcubeimage.getXDimension(), inmefcubeimage.getYDimension(), inmefcubeimage.getZDimension(), inmefcubeimage.getNExtensions(), tfloat, cNone);
            
            cout << "operaMEFCubeTest: Testing if file is a MEF Cube" << endl;
            if (!inmefcubeimage.isMEFCube()) {
                throw operaException("operaMEFCubeTest: "+mefcubeimage+" is not a MEF Cube. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
            }
            
            cout << "operaMEFCubeTest: inmefcubeimage npixels= " << inmefcubeimage.getnpixels() << endl;
            cout << "operaMEFCubeTest: outmefcubeimage npixels= " << outmefcubeimage.getnpixels() << endl;
            
            // Copy headers
            cout << "operaMEFCubeTest: Copying headers..." << endl;
            //outmefcubeimage.operaMultiExtensionFITSCubeCopyHeader(&inmefcubeimage);            
            
            // Loop through and copy each extension one by one.
            cout << "operaMEFCubeTest: Copying extensions..." << endl;
            unsigned lastextension = inmefcubeimage.getNExtensions();
            for (unsigned extension=1; extension<=lastextension; extension++) {
                outmefcubeimage[extension] = inmefcubeimage[extension]; // Exercises readExtension, saveExtension inside of []
                cout << "operaMEFCubeTest: outmefcubeimage["<<extension<<"][1][1024][1024]=" << outmefcubeimage[extension][1][1024][1024] << " inmefcubeimage["<<extension<<"][1][1024][1024]=" << inmefcubeimage[extension][1][1024][1024] << endl;
            }
            // Save file
            cout << "operaMEFCubeTest: Saving outmefcubeimage..." << endl;
            outmefcubeimage.operaMultiExtensionFITSCubeSave();
            // now copy whole image
			outmefcubeimage = inmefcubeimage;
            cout << "operaMEFCubeTest: SavingAS outmefcubeimageAll..." << endl;
            outmefcubeimage.operaMultiExtensionFITSCubeSaveAs("lazymefcubeimageAll.fits");
            // close the files
            cout << "operaMEFCubeTest: Closing files..." << endl;
            inmefcubeimage.operaFITSImageClose();
            outmefcubeimage.operaFITSImageClose();
        }
        
	}
	catch (operaException e) {
		cerr << "operaMEFCubeTest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaMEFCubeTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;

    cout << "operaCubeTest: TEST COMPLETE" << endl;
}

