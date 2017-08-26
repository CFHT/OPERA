/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaMasterCalibration
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: May/2015
 
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

#include "libraries/operaFITSImage.h"
#include "libraries/operaImage.h"
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"
#include "core-espadons/operaMasterCalibration.h"

#define MAXIMAGES 1000

using namespace std;

/*! 
 * operaMasterCalibration
 * \author Doug Teeple / Christopher Usher
 * \brief Encapsulates the functionality of modules that create master calibration images.
 * \file operaMasterCalibration.cpp
 * \ingroup libraries
 */

int MasterCalibrationCreation(int argc, char *argv[], const string moduleName, const string imagetype)
{
	operaArgumentHandler args;
	
	string listofimages;
	string imagelistfile;
	string output;
	unsigned pick;
	bool rotate;
	unsigned compressionval;
	string version;
	string date;
	
	args.AddOptionalArgument("images", listofimages, "", "List of input " + imagetype + " fits files, seperated by spaces");
	args.AddOptionalArgument("imagelistfile", imagelistfile, "", "File containing list of input " + imagetype + " fits files");
	args.AddRequiredArgument("output", output, "Master " + imagetype + " output file");
	args.AddOptionalArgument("pick", pick, 0, "Index of a specific " + imagetype + " to use (one-based, not zero-based)");
	args.AddSwitch("rotate", rotate, "Rotate output by 90 degrees");
	args.AddOptionalArgument("compressiontype", compressionval, cNone, "Ouput compression type");
	args.AddOptionalArgument("version", version, "OPERA-1.0", "The version of OPERA, to be inserted into the output fits header");
	args.AddOptionalArgument("date", date, "", "The reduction date, to be inserted into the output fits header");
	
	try  {
		args.Parse(argc, argv);
		eCompression compression = (eCompression)compressionval;
		
		string images[MAXIMAGES];
		unsigned imageIndex = 0;
		SplitStringIntoArray(listofimages, images, imageIndex, MAXIMAGES); // Split list of images into array
		ReadStringsFromFileIntoArray(imagelistfile, images, imageIndex, MAXIMAGES); // Read list of images from file
		
		if (imageIndex == 0) {
			throw operaException(moduleName + ": ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (output.empty()) {
			throw operaException(moduleName + ": ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (args.debug) for(unsigned i=0;i<imageIndex;i++) cout << images[i] << " " << i << endl;
		if (args.verbose) {
			cout << moduleName << ": output= " << output << endl;
			cout << moduleName << ": imageIndex= " << imageIndex << endl;
			cout << moduleName << ": pick= " << pick << endl;
			cout << moduleName << ": OPERA version= " << version << endl;
			cout << moduleName << ": Reduction date= " << date << endl;
		}
		// if there are not enough images to median combine, just pick one
		if ((imageIndex < 3 && pick == 0) || (imageIndex < pick)) {
			pick = 1;	
			if (args.verbose) cout << moduleName << ": too few images (" << imageIndex << "), picking " << pick << endl;
		}
		
		operaFITSImage *masterImg = NULL;
		if (pick) { // choose single input image
			if (pick > imageIndex) throw (operaErrorPickOutofRange);
			if (args.verbose) cout << moduleName << ": picking " << images[pick-1] << endl;
			operaFITSImage *imgIn = new operaFITSImage(images[pick-1], READONLY);
			masterImg = new operaFITSImage(output, imgIn->getnaxis1(), imgIn->getnaxis2(), tushort, compression);
			masterImg->operaFITSImageCopyHeader(imgIn);
			*masterImg = *imgIn;	// copy the pixels
			imgIn->operaFITSImageClose();
			delete imgIn;
		} else { // otherwise median combine stack
			if (args.verbose) cout << moduleName << ": median of " << imageIndex << " images." << endl;
			unsigned short *masterData = NULL;
			unsigned short *imgVals[MAXIMAGES];
			long npixels = 0;
			for (unsigned i=0; i<imageIndex; i++) {
				operaFITSImage *imgIn = new operaFITSImage(images[i], READONLY);
				if (i == 0) {
					masterImg = new operaFITSImage(output, imgIn->getnaxis1(), imgIn->getnaxis2(), tushort, compression);
					masterImg->operaFITSImageCopyHeader(imgIn);
					masterData = (unsigned short *)masterImg->getpixels();
					npixels = masterImg->getnpixels();
				}
				imgVals[i] = imgIn->operaFITSImageClonePixelsUSHORT();
				imgIn->operaFITSImageClose();
				delete imgIn;
			}
			imgVals[imageIndex] = NULL;
			masterData = operaArrayMedianCombineUSHORT(imageIndex, npixels, masterData, imgVals);
		}
		masterImg->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
		masterImg->operaFITSAddComment(version);
		if(pick) masterImg->operaFITSAddComment("Picking a single "+imagetype+" "+images[pick-1]);
		else {
			masterImg->operaFITSAddComment("A median of "+itos(imageIndex)+" images.");
			for (unsigned i=0; i<imageIndex; i++) {
				masterImg->operaFITSAddComment("Using "+imagetype+" image "+images[i]);
			}
		}
		masterImg->operaFITSSetHeaderValue("FILENAME", output, "Filename");
		masterImg->operaFITSDeleteHeaderKey("EXPNUM");
		masterImg->operaFITSDeleteHeaderKey("OBSID");
		if (rotate) masterImg->rotate90();
		masterImg->operaFITSImageSave();
		masterImg->operaFITSImageClose();
		delete masterImg;
	}
	catch (operaException e) {
		cerr << moduleName << ": " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << moduleName << ": " << operaStrError(errno) << endl;
	}
	return EXIT_SUCCESS;
}
