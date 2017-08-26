/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaAOBImageTest
 Version: 1.0
 Description: Perform various tests on an AOB Image.
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

/*! \file operaAOBImageTest.cpp */

using namespace std;

/*! 
 * operaAOBImageTest
 * \author Doug Teeple / Megan Tannock
 * \brief Perform various tests on the operaCube classes.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	cerr << " Usage: operaAOBImageTest -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
    string inputImage;
	
	struct option longopts[] = {
		{"inputImage",			1, NULL, 'i'},
		{"plot",			optional_argument, NULL, 'p'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
        while ((opt = getopt_long(argc, argv, "i:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
                    
				case 'i':
					inputImage = optarg;
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
			throw operaException("operaAOBImageTest: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
        // Get some information about the image.
        // getFITSImageInformation assigns values from inputImage in to variables XDimension, YDimension,
        // ZDimension, Extensions, Datatype. Will re-assign if called again with different file.
		getFITSImageInformation(inputImage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels);
        cout << "operaAOBImageTest: IMAGE: " << inputImage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
        
        /*
         * Read entire image in to memory, basic tests
         */
        {
            cout << "operaAOBImageTest: Entire image in memory tests... " << endl;
            
            cout << "operaAOBImageTest: Calling constructors" << endl;
            operaFITSCube inCube(inputImage, tfloat, READONLY, cNone);
            operaFITSCube outCube("outCube.fits", inCube.getXDimension(), inCube.getYDimension(), inCube.getZDimension(), tfloat, cNone);
            operaFITSImage outImage("outImage.fits", inCube.getXDimension(), inCube.getYDimension(), tfloat, cNone);
            
            cout << "operaAOBImageTest: Testing if file is a FITS Cube" << endl;
            if (!(inCube.isCube())) {
                throw operaException("operaAOBImageTest: "+inputImage+" is not a Cube. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
            }          
            
            cout << "operaAOBImageTest: Copying headers from inCube to outImage and outCube..." << endl;
            outImage.operaFITSImageCopyHeader(&inCube);
            outCube.operaFITSImageCopyHeader(&inCube); 
            cout << "operaAOBImageTest: Copying pixels from inCube to outCube..." << endl;
            outCube = inCube;	

            // Testing operaFITSCubeSave, operaFITSCubeSaveAs
            cout << "operaAOBImageTest: Saving outCube as outCube.fits... (operaFITSCubeSave)" << endl;
            outCube.operaFITSCubeSave(); 
            cout << "operaAOBImageTest: Saving outCube as outCubeSaveAs.fits... (operaFITSCubeSaveAs)" << endl;
            outCube.operaFITSCubeSaveAs("outCubeSaveAs.fits"); 
            
            // Print some pixels to test they were copied correctly... index by [slice][y][x]
            cout << "operaAOBImageTest: inCube[1][1632][2464]=" << inCube[1][1632][2464] << " outCube[1][1632][2464]=" << outCube[1][1632][2464] << endl;
            cout << "operaAOBImageTest: inCube[2][1632][2464]=" << inCube[2][1632][2464] << " outCube[2][1632][2464]=" << outCube[2][1632][2464] << endl;
            cout << "operaAOBImageTest: inCube[3][1632][2464]=" << inCube[3][1632][2464] << " outCube[3][1632][2464]=" << outCube[3][1632][2464] << endl;
            
            // Do a median collapse
            cout << "operaAOBImageTest: Doing median collapse of all " << outCube.getslices() << " slices..." << endl;
            outCube.medianCollapse();  // Commented because it takes a very long time to run through
            outImage = outCube[1];
            cout << "operaAOBImageTest: Saving outImage (median collapse of outCube) as outCubeMedian.fits... (operaFTSImageSaveAs)" << endl;
            outImage.operaFITSImageSaveAs("outCubeMedian.fits");     
            
            // Take the mean.... First need to re-copy pixels (scrambled by medianCollapse)
            cout << "operaAOBImageTest: Re-copying pixels from inCube to outCube..." << endl;
            outCube = inCube;
            cout << "operaAOBImageTest: Taking the mean of all " << outCube.getslices() << " slices..." << endl;
            for (unsigned slice=2; slice<=outCube.getslices(); slice++) {
                outCube[1] += outCube[slice];   // add all slices together
            }
            outCube[1] /= outCube.getslices();  // divide by number of slices
            outImage = outCube[1];  // [] exercises setSlice method, gets just the first slice
            cout << "operaAOBImageTest: Saving outImage as outCubeMean.fits... (operaFTSImageSaveAs)" << endl;
            outImage.operaFITSImageSaveAs("outCubeMean.fits");
           
            // Now, create an operaFITSImage from a single slice (downsize from a FITS Cube to a FITS Image)
            cout << "operaAOBImageTest: Downsizing slice 3 to a FITS Image..." << endl;
            operaFITSImage downSize(outCube[3], false, true);
            cout << "operaAOBImageTest: downSize[1024][1024]=" << downSize[1024][1024] << endl;
            cout << "operaAOBImageTest: Saving downsize as outCubeDownSizeSlice3.fits..." << endl;
            downSize.operaFITSImageSaveAs("outCubeDownSizeSlice3.fits");
            
            cout << "operaAOBImageTest: Closing file(s)..." << endl;
            outCube.operaFITSImageClose();
            inCube.operaFITSImageClose();
            outImage.operaFITSImageClose();
            downSize.operaFITSImageClose();
        }
    }
	catch (operaException e) {
		cerr << "operaAOBImageTest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaAOBImageTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
    
    cout << "operaAOBImageTest: TEST COMPLETE" << endl;
}

