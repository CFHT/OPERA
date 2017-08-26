/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaPixelSensitivityMap
 Version: 1.0
 Description: Create a normalized map of pixel-by-pixel sensitivity variations  
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

#include "libraries/operaFITSImage.h"
#include "libraries/operaFit.h"	
#include "libraries/operaArgumentHandler.h"

/*! \file operaPixelSensitivityMap.cpp */

using namespace std;

/*! 
 * operaPixelSensitivityMap
 * \author Eder Martioli
 * \brief Create a normalized map of pixel-by-pixel sensitivity variations.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
	
	string outputfilename;
	string masterbias; 
	string masterflat; 
	string badpixelmask;
	unsigned compression_val = cNone;
    
    args.AddRequiredArgument("outputfilename", outputfilename, "Ouput pixel-by-pixel sensitivity map FITS image");
    args.AddRequiredArgument("masterbias", masterbias, "Input Master Bias FITS image");
    args.AddRequiredArgument("masterflat", masterflat, "Input Master Flat-Field FITS image");
    args.AddRequiredArgument("badpixelmask", badpixelmask, "FITS image with badpixel mask");
    args.AddOptionalArgument("compressiontype", compression_val, cNone, "Compression type");
	
	try {
		args.Parse(argc, argv);
		
		// we need a masterflat...
		if (masterflat.empty()) {
			throw operaException("operaPixelSensitivityMap: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}			
        // we need an output...
		if (outputfilename.empty()) {
			throw operaException("operaPixelSensitivityMap: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		
		eCompression compression = eCompression(compression_val);
		
		if (args.verbose) {
			cout << "operaPixelSensitivityMap: outputfilename = " << outputfilename << endl; 
			cout << "operaPixelSensitivityMap: masterflat = " << masterflat << endl; 	
			cout << "operaPixelSensitivityMap: masterbias = " << masterbias << endl; 		
			cout << "operaPixelSensitivityMap: badpixelmask = " << badpixelmask << endl;
            cout << "operaPixelSensitivityMap: compression = " << compression << endl;
		}
        
		operaFITSImage flat(masterflat, tfloat, READONLY);
		
        operaFITSImage *bias = NULL;
        operaFITSImage *badpix = NULL;
        
		if (!masterbias.empty()){
			bias = new operaFITSImage(masterbias, tfloat, READONLY);
		} else {
            bias = new operaFITSImage(flat.getnaxis1(),flat.getnaxis2(),tfloat);
            *bias = 0.0;
        }
        
		if (!badpixelmask.empty()){
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(flat.getnaxis1(),flat.getnaxis2(),tfloat);
            *badpix = 1.0;
        }

		flat -= *bias;			// remove bias from masterflat
		
        long npixels = flat.getnpixels();
        
        float *flatData = (float *)flat.getpixels();
        float *badpixData = (float *)badpix->getpixels();
        
        float maxvalue = -BIG;
        
        for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
            if(maxvalue < flatData[pixIndex] && badpixData[pixIndex]) {
                maxvalue = flatData[pixIndex];
            }
        }
        
        operaFITSImage outputImage(outputfilename, flat.getnaxis1(), flat.getnaxis2(), tfloat, compression);
		outputImage.operaFITSImageCopyHeader(&flat);
        
		for (unsigned y=0; y<flat.getnaxis2(); y++) {
			for (unsigned x=0; x<flat.getnaxis1(); x++) {
				outputImage[y][x] = (flat[y][x] + (*bias)[y][x]) / maxvalue;
			}
		}
        
        outputImage.operaFITSImageSave();
		outputImage.operaFITSImageClose();

        flat.operaFITSImageClose();
		        
        if(bias) delete bias;
        if(badpix) delete badpix;
	}
	catch (operaException e) {
		cerr << "operaPixelSensitivityMap: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaPixelSensitivityMap: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
