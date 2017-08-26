/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaCubeTest
 Version: 1.0
 Description: Perform various tests on the operaCube and operaCube classes.
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
#include "libraries/operaFITSCube.h"
#include "libraries/operaLibCommon.h"

/*! \file operaCubeTest.cpp */

using namespace std;

/*! 
 * operaCubeTest
 * \author Doug Teeple / Megan Tannock
 * \brief Perform various tests on the operaCube classes.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	cerr << " Usage: basicCubeTest -[dvth]\n";
}	

int main(int argc, char *argv[])
{
	int opt;	
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
		
#ifdef __APPLE__
		string basepath = "/data/";
#else
		string basepath = "/data/uhane5/opera/";
#endif
        string cubeimage = "Cube/Cube.fits";
        unsigned XDimension, YDimension, ZDimension, Extensions;
		edatatype Datatype;
		long Npixels;
        
        // Get some information about the image.
        // getFITSImageInformation assigns values from cubeimage in to variables XDimension, YDimension,
        // ZDimension, Extensions, Datatype. Will re-assign if called again with different file.
		getFITSImageInformation(cubeimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
        cout << "operaCubeTest: IMAGE: " << cubeimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;  
        
        /*
         * Read entire image in to memory, basic tests
         */
        {
            cout << "operaCubeTest: Entire image in memory tests... " << endl;
            
            cout << "operaCubeTest: Calling constructors" << endl;
            operaFITSCube inCube(cubeimage, tfloat, READONLY, cNone);
            operaFITSCube outCube("outCube.fits", inCube.getXDimension(), inCube.getYDimension(), inCube.getZDimension(), tfloat, cNone);
            operaFITSImage outImage("outImage.fits", inCube.getXDimension(), inCube.getYDimension(), tfloat, cNone);
            
            cout << "operaCubeTest: Testing if file is a FITS Cube" << endl;
            if (!(inCube.isCube())) {
                throw operaException("operaCubeTest: "+cubeimage+" is not a Cube. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
            }          
            
            cout << "operaCubeTest: Copying headers from inCube to outImage and outCube..." << endl;
            outImage.operaFITSImageCopyHeader(&inCube);
            outCube.operaFITSImageCopyHeader(&inCube); 
            cout << "operaCubeTest: Copying pixels from inCube to outCube..." << endl;
            outCube = inCube;	

            // Testing operaFITSCubeSave, operaFITSCubeSaveAs
            cout << "operaCubeTest: Saving outCube as outCube.fits... (operaFITSCubeSave)" << endl;
            outCube.operaFITSCubeSave(); 
            cout << "operaCubeTest: Saving outCube as outCubeSaveAs.fits... (operaFITSCubeSaveAs)" << endl;
            outCube.operaFITSCubeSaveAs("outCubeSaveAs.fits"); 
            
            // Print some pixels to test they were copied correctly... index by [slice][y][x]
            cout << "operaCubeTest: inCube[1][1632][2464]=" << inCube[1][1632][2464] << " outCube[1][1632][2464]=" << outCube[1][1632][2464] << endl;
            cout << "operaCubeTest: inCube[2][1632][2464]=" << inCube[2][1632][2464] << " outCube[2][1632][2464]=" << outCube[2][1632][2464] << endl;
            cout << "operaCubeTest: inCube[3][1632][2464]=" << inCube[3][1632][2464] << " outCube[3][1632][2464]=" << outCube[3][1632][2464] << endl;
            
            // Do a median collapse
            cout << "operaCubeTest: Doing median collapse of all " << outCube.getslices() << " slices..." << endl;
            outCube.medianCollapse();  // Commented because it takes a very long time to run through
            outImage = outCube[1];
            cout << "operaCubeTest: Saving outImage (median collapse of outCube) as outCubeMedian.fits... (operaFTSImageSaveAs)" << endl;
            outImage.operaFITSImageSaveAs("outCubeMedian.fits");     
            
            // Take the mean.... First need to re-copy pixels (scrambled by medianCollapse)
            cout << "operaCubeTest: Re-copying pixels from inCube to outCube..." << endl;
            outCube = inCube;
            cout << "operaCubeTest: Taking the mean of all " << outCube.getslices() << " slices..." << endl;
            for (unsigned slice=2; slice<=outCube.getslices(); slice++) {
                outCube[1] += outCube[slice];   // add all slices together
            }
            outCube[1] /= outCube.getslices();  // divide by number of slices
            outImage = outCube[1];  // [] exercises setSlice method, gets just the first slice
            cout << "operaCubeTest: Saving outImage as outCubeMean.fits... (operaFTSImageSaveAs)" << endl;
            outImage.operaFITSImageSaveAs("outCubeMean.fits");
           
            // Now, create an operaFITSImage from a single slice (downsize from a FITS Cube to a FITS Image)
            cout << "operaCubeTest: Downsizing slice 3 to a FITS Image..." << endl;
            operaFITSImage downSize(outCube[3], false, true);
            cout << "operaCubeTest: downSize[1024][1024]=" << downSize[1024][1024] << endl;
            cout << "operaCubeTest: Saving downsize as outCubeDownSizeSlice3.fits..." << endl;
            downSize.operaFITSImageSaveAs("outCubeDownSizeSlice3.fits");
            
            cout << "operaCubeTest: Closing file(s)..." << endl;
            outCube.operaFITSImageClose();
            inCube.operaFITSImageClose();
            outImage.operaFITSImageClose();
            downSize.operaFITSImageClose();
        }
    }
	catch (operaException e) {
		cerr << "operaCubeTest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaCubeTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
    
    cout << "operaCubeTest: TEST COMPLETE" << endl;
}

