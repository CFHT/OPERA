
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: FITSImageVectorTest.cpp
 Version: 1.0
 Description: Perform various tests on the operaFITSImageVector class.
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
#include "libraries/operaFITSImageVector.h"

/*! \file FITSImageVectorTest.cpp */

using namespace std;

/*! 
 * FITSImageVectorTest
 * \author Doug Teeple
 * \brief Perform various tests on the operaFITSImageVector class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: FITSImageVectorTest -[dvth]\n";
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
        string fitsimages[] = {
			"/data/espadons/11AQ14-Jul08/1315157o.fits",
			"/data/espadons/11AQ14-Jul08/1315158o.fits",
			"/data/espadons/11AQ14-Jul08/1315159o.fits",
			"/data/espadons/11AQ14-Jul08/1315160o.fits",
			"/data/espadons/11AQ14-Jul08/1315161o.fits",
			""
		};
#else
		string basepath = "/data/uhane5/opera/";
        string fitsimages[] = {
			"/data/espadons/11AQ14-Jul08/1315157o.fits",
			"/data/espadons/11AQ14-Jul08/1315158o.fits",
			"/data/espadons/11AQ14-Jul08/1315159o.fits",
			"/data/espadons/11AQ14-Jul08/1315160o.fits",
			"/data/espadons/11AQ14-Jul08/1315161o.fits",
			""
		};
#endif
        string fitsimage = basepath + "espadons/11AQ14-Jul08/1315157o.fits";
		string mefimage = basepath + "/wircam/12AQ03-Mar01/1529925o.fits";
		string mefcubeimage = basepath + "/wircam/12AQ03-Mar01/1530725d.fits";
		string cubeimage = basepath + "/Cube/Cube.fits";

		/*
		 * Step 1 : Simple FITSImage
		 */
		if (step < 0 || step == 1) {

            const unsigned saturation = 50000;
            const unsigned bias = 300;
			
			cout << "FITSImageVectorTest: SIMPLE FITS IMAGE LIST" << endl;
			operaFITSImageVector<operaFITSImage>infitsimages(fitsimages);
			infitsimages.addImage("/data/espadons/11AQ14-Jul08/1315162o.fits");
			
			cout << "FITSImageVectorTest: removing bias (" << bias << ")..." << endl;
			infitsimages.removebias(bias);
			
			operaFITSImage *stack = infitsimages.stack();
			cout << "FITSImageVectorTest: stack[1024][1024]=" << (*stack)[1024][1024] << endl;
			
			operaFITSImage *median = infitsimages.median();
			cout << "FITSImageVectorTest: median[1024][1024]=" << (*median)[1024][1024] << endl;
			
			operaFITSImage *mean = infitsimages.mean();
			cout << "FITSImageVectorTest: mean[1024][1024]=" << (*mean)[1024][1024] << endl;
			
			operaFITSImage *meanunsaturated = infitsimages.meanunsaturated(saturation);
			cout << "FITSImageVectorTest: < " << saturation << " meanunsaturated[1024][1024]=" << (*meanunsaturated)[1024][1024] << endl;
			
			operaFITSImage *meanunsaturatedbyetime = infitsimages.meanunsaturatedbyetime(saturation);
			cout << "FITSImageVectorTest: < " << saturation << " meanunsaturatedbyetime[1024][1024]=" << (*meanunsaturatedbyetime)[1024][1024] << endl;
			
			cout << "FITSImageVectorTest: Saving stack as stack.fits..." << endl;
			stack->operaFITSImageSaveAs("stack.fits");
			
			cout << "FITSImageVectorTest: Saving median as median.fits..." << endl;
			median->operaFITSImageSaveAs("median.fits");
			
			cout << "FITSImageVectorTest: Saving mean as mean.fits..." << endl;
			mean->operaFITSImageSaveAs("mean.fits");
			
			cout << "FITSImageVectorTest: Saving meanunsaturated as meanunsaturated.fits..." << endl;
			meanunsaturated->operaFITSImageSaveAs("meanunsaturated.fits");
			
			cout << "FITSImageVectorTest: Saving meanunsaturatedbyetime as meanunsaturatedbyetime.fits mean etime=" << infitsimages.getMeanEtime() << endl;
			meanunsaturatedbyetime->operaFITSImageSaveAs("meanunsaturatedbyetime.fits");
			
			cout << "FITSImageVectorTest: Closing files..." << endl;
			meanunsaturatedbyetime->operaFITSImageClose();
			meanunsaturated->operaFITSImageClose();
			mean->operaFITSImageClose();
			median->operaFITSImageClose();
			stack->operaFITSImageClose();
			infitsimages.close();
			
			cout << "FITSImageVectorTest: fitsimage test complete." << endl;
		}
		
		/*
		 * Step 2 : MEF FITSImage
		 */
		if (step < 0 || step == 2) {
            cout << "FITSImageVectorTest: MEF FITS IMAGE" << endl;
            getFITSImageInformation(mefimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "FITSImageVectorTest: mefimage " << mefimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;

			operaMultiExtensionFITSImage infitsimage(mefimage, tfloat, READONLY, cNone, false);
			operaMultiExtensionFITSImage outfitsimage("mefimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), infitsimage.getNExtensions(), tfloat, cNone, false);
			
			cout << "FITSImageVectorTest: mefimage test complete." << endl;
		}
        
        /*
         * Step 3 : Simple MEF
         */
        if (step < 0 || step == 3) {
            cout << "FITSImageVectorTest: SIMPLE MEF FITS IMAGE" << endl;            
			operaMultiExtensionFITSImage infitsimage(mefimage, tfloat, READONLY, cNone, false);
			operaMultiExtensionFITSImage outfitsimage("simplemefimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), infitsimage.getNExtensions(), tfloat, cNone, false);
			
			outfitsimage.operaFITSImageClose();
			cout << "FITSImageVectorTest: simple mefimage test complete." << endl;
        }
		
		/*
		 * Step 4 : Cube FITSImage
		 */ 
		if (step < 0 || step == 4) {
            cout << "FITSImageVectorTest: CUBE FITS IMAGE" << endl;
            getFITSImageInformation(cubeimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "FITSImageVectorTest: cubeimage " << cubeimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;

			operaFITSCube infitsimage(cubeimage, tfloat, READONLY, cNone, false);
			operaFITSCube outfitsimage("cubeimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), infitsimage.getZDimension(), tfloat, cNone, false);
		 	
			cout << "FITSImageVectorTest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			cout << "FITSImageVectorTest: cubeimage test complete." << endl;
		}
		
		/*
		 * Step 5 : MEFCubeImage
		 */
		if (step < 0 || step == 5) {
            cout << "FITSImageVectorTest: MEF CUBE IMAGE" << endl;
            getFITSImageInformation(mefcubeimage, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
            cout << "FITSImageVectorTest: mefcubeimage " << mefcubeimage << " X= " << XDimension << " Y= " << YDimension << " Z= " << ZDimension << " ext= " << Extensions << " Datatype= " << Datatype << " Npixels= " << Npixels << endl;
            
			operaMultiExtensionFITSCube infitsimage(mefcubeimage, tfloat, READONLY, cNone, false);
			operaMultiExtensionFITSCube outfitsimage("mefcubeimage.fits", infitsimage.getXDimension(), infitsimage.getYDimension(), infitsimage.getZDimension(), infitsimage.getNExtensions(), tfloat, cNone, false);
			
			cout << "FITSImageVectorTest: Closing file(s)..." << endl;
			infitsimage.operaFITSImageClose();
			outfitsimage.operaFITSImageClose();
			cout << "FITSImageVectorTest: mefcubeimage test complete." << endl;
		}
		
	}
	catch (operaException e) {
		cerr << "FITSImageVectorTest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "FITSImageVectorTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

