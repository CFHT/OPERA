/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaWIRCamImageTest
 Version: 1.0
 Description: Perform various tests on the operaWIRCamImageTest classes.
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
#include "libraries/Polynomial.h"	
#include "libraries/operaException.h"
#include "libraries/operaGeometricShapes.h"				// for Box
#include "libraries/operaStats.h"						// for operaArrayMedianQuick
#include "libraries/operaWIRCamImage.h"

/*! \file operaWIRCamImageTest.cpp */

using namespace std;

/*! 
 * operaWIRCamImageTest
 * \author Doug Teeple / Megan Tannock
 * \brief Perform various tests on the operaWIRCamImage classes.
 * \arg argc
 * \arg argv
 * \return EXIT_STATUS
 * \ingroup test
 */

/* Print out the proper program usage syntax */
void printUsageSyntax() {
	
	cerr << " Usage: basicFITSImageTest -[dvth]\n";
    
}	

int main(int argc, char *argv[])
{
	int opt;
	
	unsigned debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		
		{"plot",			optional_argument, NULL, 'p'},
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
        
        cout << "operaWIRCamImageTest: ---------- DETREND TEST ----------" << endl;
        // This will be a test of the steps for the detrend module.
        // Uses AGGRESSIVE read only, lazy will need to be tested
        
        unsigned XDimension, YDimension, ZDimension, Extensions;
		edatatype Datatype;
		long Npixels;        
		string basepath = "/data/WIRCam/12BF12/";
		// all files also available at: /h/tannock/OPERAreduction-1579272/

        // Define the paths to the images
        string inimagename = basepath + "1579272o.fits";      // raw image
        string outimagename = "1579272s-opera.fits";     // detrended image - to be created
		string weightmapname = basepath + "1579272w.fits";    // weight map - to be created
		string mastertwilightflatname = basepath + "mastertwilightflat_H_12Bw02_v100.fits";  // mastertwilight flat
		string masterdarkname = basepath + "masterdark_010s_12Bw02_v200.fits"; // master dark
		string badpixelmaskname = basepath + "badpix16_20121004HST055237_v100.fits"; // bad pixel map
		string saturationmapname = "1579272saturation.fits"; // saturation map - to be created
		string referencepixelsmapname = "1579272referencepixels.fits"; // referencepixels map - to be created
		
        // Check that the raw image, master twilight flat, master dark, and bad pixel mask exist
		if (inimagename.empty()) {
			throw operaException("operaWIRCamImageTest: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (mastertwilightflatname.empty()) {
			throw operaException("operaWIRCamImageTest: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (masterdarkname.empty()) {
			throw operaException("operaWIRCamImageTest: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (badpixelmaskname.empty()) {
			throw operaException("operaWIRCamImageTest: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
        
        // Print some information about the raw image (inImage).
        getFITSImageInformation(inimagename, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &Npixels); 
        cout << "operaWIRCamImageTest: inImage = " << inimagename << "   X= " << XDimension << "   Y= " << YDimension << "   Z= " << ZDimension << "   Extensions= " << Extensions << "   Datatype= " << Datatype << "   Npixels= " << Npixels << endl;
		
		// Read in and collapse the raw science image
        cout << "operaWIRCamImageTest: Reading raw file in to memory: " << inimagename << endl;
		operaWIRCamImage inImage(inimagename, tfloat, READONLY, cNone, false);
        cout << "operaWIRCamImageTest: Median collapsing the in image...." << endl;
		inImage % ZDimension;	// median combine all slices
        cout << "operaWIRCamImageTest: inImage = " << inimagename << "   X= " << (XDimension=inImage.getXDimension()) << "   Y= " << (YDimension=inImage.getYDimension()) << "   Z= " << (ZDimension=inImage.getZDimension()) << "   Extensions= " << (Extensions=inImage.getNExtensions()) << "   Datatype= " << Datatype << "   Npixels= " << (Npixels=inImage.getnpixels()) << endl;
 
		// Create output and reference pixel images
		cout << "operaWIRCamImageTest: Creating output image: " << outimagename << endl;
		operaWIRCamImage outImage(outimagename, inImage.getXDimension(), inImage.getYDimension(), inImage.getZDimension(), inImage.getNExtensions(), tfloat, cNone, false);
        cout << "operaWIRCamImageTest: Creating referencepixels map: " << referencepixelsmapname << endl;
        operaWIRCamImage referencepixels(referencepixelsmapname, inImage.getXDimension(), inImage.getYDimension(), inImage.getZDimension(), inImage.getNExtensions(), tfloat, cNone, false);
		// Read in master twilight flat, master dark, and bad pixel map
		cout << "operaWIRCamImageTest: Reading master twilight flat in to memory: " << mastertwilightflatname << endl;
		operaWIRCamImage masterTwilightFlat(mastertwilightflatname, tfloat, READONLY, cNone, false);
        cout << "operaWIRCamImageTest: Reading master dark in to memory: " << masterdarkname << endl;
		operaWIRCamImage masterDark(masterdarkname, tfloat, READONLY, cNone, false);
        cout << "operaWIRCamImageTest: Reading bad pixel mask in to memory: " << badpixelmaskname << endl;
		operaWIRCamImage badPixelMask(badpixelmaskname, tfloat, READONLY, cNone, false);
   
        cout << "operaWIRCamImageTest: Copying inImage header to outImage..." << endl;
        outImage.operaMultiExtensionFITSCubeCopyHeader(&inImage); // Copy image header
        cout << "operaWIRCamImageTest: Copying inImage pixels to outImage...." << endl;
        outImage = inImage;		// copy entire image pixels
        
        // Detrend Steps
        cout << "operaWIRCamImageTest: Starting detrend steps." << endl;
        
        // DETREND STEPS:
        // 1. Reference pixels correction
        // 2. Dark subtraction
        // 3. Apply non-linearity correction for each chip
        // 4. Flat field correction
        // 5. Bad pixel mask
                     
        // 1. referencepixels correction (removes reference pixels)
        cout << "operaWIRCamImageTest: Copying inImage header to referencepixels..." << endl;
        //referencepixels.operaMultiExtensionFITSCubeCopyHeader(&inImage);
        cout << "operaWIRCamImageTest: Creating referencepixels map." << endl;
        referencepixels.createReferencePixelImage(inImage);
        cout << "operaWIRCamImageTest: Saving referencepixels.fits..." << endl;
        referencepixels.operaWIRCamImageSave();
        cout << "operaWIRCamImageTest: Applying referencepixels correction (removing reference pixels)." << endl;
        outImage -= referencepixels;
		outImage.operaWIRCamImageSaveAs("opera-1-1579272-refpixremoved.fits"); // Save the intermediate file
               
        // 2. Dark subtraction
        cout << "operaWIRCamImageTest: Applying dark subtraction." << endl;
        outImage -= masterDark;
		outImage.operaWIRCamImageSaveAs("opera-2-1579272-darkremoved.fits"); // Save the intermediate file
        
        // 3. Apply non-linearity correction for each chip
        cout << "operaWIRCamImageTest: Applying non linearity correction." << endl;    
        // Declare four 2nd order polynomial for each chip (non linearity correction)
        cout << "operaWIRCamImageTest: Creating polynomials for each chip." << endl;
        Polynomial chip1(3);
        Polynomial chip2(3);
        Polynomial chip3(3);
        Polynomial chip4(3);
        
        cout << "operaWIRCamImageTest: Assigning polynomial coefficients for each chip." << endl;
//        chip1.setCoefficient(0, 0.991276); 	// April 5 2008
//        chip1.setCoefficient(1, 1.72141e-06); // April 5 2008
//        chip1.setCoefficient(2, 7.57266e-12); // April 5 2008
//        chip2.setCoefficient(0, 0.993746); 	// April 5 2008
//        chip2.setCoefficient(1, 1.82717e-06); // April 5 2008
//        chip2.setCoefficient(2, 1.93950e-11); // April 5 2008
//        chip3.setCoefficient(0, 0.994254); 	// April 5 2008
//        chip3.setCoefficient(1, 1.92776e-06); // April 5 2008
//        chip3.setCoefficient(2, 1.69437e-11); // April 5 2008
//        chip4.setCoefficient(0, 0.994799); 	// April 5 2008
//        chip4.setCoefficient(1, 1.49037e-06); // April 5 2008
//        chip4.setCoefficient(2, 2.85603e-11); // April 5 2008
		chip1.setCoefficient(0, 0.998810);		// July 16 2007
        chip1.setCoefficient(1, 7.92059e-07);	// July 16 2007
        chip1.setCoefficient(2, 5.39334e-11);	// July 16 2007
        chip2.setCoefficient(0, 0.996922); 		// July 16 2007
        chip2.setCoefficient(1, 1.28674e-06);	// July 16 2007
        chip2.setCoefficient(2, 4.18800e-11);	// July 16 2007
        chip3.setCoefficient(0, 0.998038); 		// July 16 2007
        chip3.setCoefficient(1, 1.15189e-06);	// July 16 2007
        chip3.setCoefficient(2, 4.63836e-11);	// July 16 2007
        chip4.setCoefficient(0, 0.996537);		// July 16 2007
        chip4.setCoefficient(1, 9.29331e-07);	// July 16 2007
        chip4.setCoefficient(2, 5.25666e-11);	// July 16 2007
        
        cout << "operaWIRCamImageTest: Applying non-linear correction " << endl;
		outImage.applyNonLinearCorrection(1, chip1);
        cout << "operaWIRCamImageTest: chip1 non-linear correction complete " << endl;
		outImage.applyNonLinearCorrection(2, chip2);
        cout << "operaWIRCamImageTest: chip2 non-linear correction complete " << endl;
		outImage.applyNonLinearCorrection(3, chip3);
        cout << "operaWIRCamImageTest: chip3 non-linear correction complete " << endl;
		outImage.applyNonLinearCorrection(4, chip4);
        cout << "operaWIRCamImageTest: chip4 non-linear correction complete " << endl;
        
        // Update Headers
        cout << "operaWIRCamImageTest: Updating header values. " << endl;
        outImage.operaFITSSetHeaderValue("NLCORR", "yes", "Non-linearity correction applied?");
		outImage.operaFITSSetHeaderValue("NLC_FUNC","xc/xm=a0+a1*xm+a2*xm^2","Non-linearity function");
		outImage.operaFITSSetHeaderValue("NLC_A0",chip1.Get(0),"Non-linearity function parameters", 1);
		outImage.operaFITSSetHeaderValue("NLC_A1",chip1.Get(1),"Non-linearity function parameters", 1);
		outImage.operaFITSSetHeaderValue("NLC_A2",chip1.Get(2),"Non-linearity function parameters", 1);
		outImage.operaFITSSetHeaderValue("NLC_A0",chip2.Get(0),"Non-linearity function parameters", 2);
		outImage.operaFITSSetHeaderValue("NLC_A1",chip2.Get(1),"Non-linearity function parameters", 2);
		outImage.operaFITSSetHeaderValue("NLC_A2",chip2.Get(2),"Non-linearity function parameters", 2);
		outImage.operaFITSSetHeaderValue("NLC_A0",chip3.Get(0),"Non-linearity function parameters", 3);
		outImage.operaFITSSetHeaderValue("NLC_A1",chip3.Get(1),"Non-linearity function parameters", 3);
		outImage.operaFITSSetHeaderValue("NLC_A2",chip3.Get(2),"Non-linearity function parameters", 3);
		outImage.operaFITSSetHeaderValue("NLC_A0",chip4.Get(0),"Non-linearity function parameters", 4);
		outImage.operaFITSSetHeaderValue("NLC_A1",chip4.Get(1),"Non-linearity function parameters", 4);
		outImage.operaFITSSetHeaderValue("NLC_A2",chip4.Get(2),"Non-linearity function parameters", 4);

		outImage.operaWIRCamImageSaveAs("opera-3-1579272-nonlinearitycorrected.fits"); // Save the intermediate file
        
        // 4. Flat field correction
        cout << "operaWIRCamImageTest: Applying flat field correction." << endl;
        outImage /= masterTwilightFlat;
		outImage.operaWIRCamImageSaveAs("opera-4-1579272-flatfieldcorrected.fits"); // Save the intermediate file
		
        // 5. Bad pixel mask
        cout << "operaWIRCamImageTest: Setting bad pixels in outImage to 0" << endl;
		outImage *= badPixelMask; // Set bad pixel values in outImage to 0
		outImage.operaWIRCamImageSaveAs("opera-5-1579272-badpixmask.fits"); // Save the intermediate file

        // Detrend complete
        // Save Detrended Image
        cout << "operaWIRCamImageTest: Saving outImage.fits..." << endl;
        outImage.operaWIRCamImageSave();
        cout << "operaWIRCamImageTest: Detrend steps complete." << endl;
        
        // Close files
        cout << "operaWIRCamImageTest: Closing files..." << endl;
        inImage.operaFITSImageClose();
		outImage.operaFITSImageClose();
		masterTwilightFlat.operaFITSImageClose();
		masterDark.operaFITSImageClose();
		badPixelMask.operaFITSImageClose();
		referencepixels.operaFITSImageClose();
        
        cout << "operaWIRCamImageTest: Test complete." << endl;
	}
	catch (operaException e) {
		cerr << "operaWIRCamImageTest: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaWIRCamImageTest: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

