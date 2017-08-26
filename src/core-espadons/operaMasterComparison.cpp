/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMasterComparison
 Version: 1.0
 Description: Create a master comparison.
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
#include "libraries/operaImage.h"
#include "libraries/operaStats.h" 
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

#define MAXIMAGES 1000
#ifndef SATURATIONLIMIT
#define SATURATIONLIMIT (unsigned short)65535  // this should be retrieved from the config/param file
#endif

/*! \file operaMasterComparison.cpp */

using namespace std;

/*! 
 * operaMasterComparison
 * \author Doug Teeple
 * \brief Creates a master Comparison FITS image from a list of input Comparison FITS file names.
 * \arg argc
 * \arg argv
 * \note [--images=...]* --output=... [ --pick=\<posint\>0\> ]
 * \note Pick one or median stack images.
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \ingroup core
 * \return EXIT_STATUS
 */
 
int main(int argc, char *argv[])
{
	operaArgumentHandler args;
	
	string listofimages;
	string imagelistfile;
	string output;
	string badpixelmask;
	string masterbias;
	unsigned pick = 0;
	unsigned compressionval;
	string expTimeFITSKeyword = "EXPTIME";
    unsigned combineMethod = 1;
    /*
     * combineMethod 1: Median combine stack.
     *
     * combineMethod 2: Stack and save total exposure time per pixel. Note that
     * a given pixel may be saturated in some images and not saturated in others
     * so the only way to bring all pixel counts to the same units is
     * by knowing the total exposure time accounted for each pixel. The other trick
     * here is that the masterData should be initiallized as biasData to avoid negative
     * assignments to the unsigned variable.
     *
     * combineMethod 3: A simple pixel-by-pixel stack ignoring saturated pixels.
     */
	
	unsigned saturationLimit = SATURATIONLIMIT;
    double outputExposureTime = 60;    
    bool biasConstant = true;
    bool truncateOuputFluxToSaturation = true;
	string version;
	string date;
	
	args.AddOptionalArgument("images", listofimages, "", "List of input comparison fits files, seperated by spaces");
	args.AddOptionalArgument("imagelistfile", imagelistfile, "", "File containing list of input comparison fits files");
	args.AddRequiredArgument("output", output, "Master comparison output file");
	args.AddRequiredArgument("badpixelmask", badpixelmask, "Bad pixel mask fits image of pixels to ignore");
	args.AddOptionalArgument("masterbias", masterbias, "", "Master bias fits file");
	args.AddOptionalArgument("pick", pick, 0, "Index of a specific comparison to use (one-based, not zero-based)");
	args.AddOptionalArgument("compressiontype", compressionval, cNone, "Ouput compression type");
	args.AddRequiredArgument("expTimeFITSKeyword", expTimeFITSKeyword, "String to identify exposure time header keyword");
	args.AddRequiredArgument("combineMethod", combineMethod, "Method for combining images: 1 = Median, 2 = Mean weighted by exposure time, 3 = Sum");
	args.AddRequiredArgument("saturationLimit", saturationLimit, "Define saturation or linearity limit");
	args.AddRequiredArgument("outputExposureTime", outputExposureTime, "Exposure time to reset output image (only used in combineMethod=2)");
	args.AddRequiredArgument("biasConstant", biasConstant, "Use median bias constant to add to output image. Useful to avoid negative numbers or when bias has low gradient.");
	args.AddRequiredArgument("truncateOuputFluxToSaturation", truncateOuputFluxToSaturation, "Limit maximum possible flux value to saturation");
	args.AddOptionalArgument("version", version, "OPERA-1.0", "The version of OPERA, to be inserted into the output fits header");
	args.AddOptionalArgument("date", date, "", "The reduction date, to be inserted into the output fits header");
	
	try {
		args.Parse(argc, argv);
		eCompression compression = (eCompression)compressionval;
		
		string images[MAXIMAGES];
		unsigned imageIndex = 0;
		SplitStringIntoArray(listofimages, images, imageIndex, MAXIMAGES); // Split list of images into array
		ReadStringsFromFileIntoArray(imagelistfile, images, imageIndex, MAXIMAGES); // Read list of images from file
        
		if (args.verbose) {
            cout << "operaMasterComparison: imagelistfile = " << imagelistfile << endl;
            cout << "operaMasterComparison: output = " << output << endl;
            cout << "operaMasterComparison: compression = " << compression << endl;
            cout << "operaMasterComparison: pick = " << pick << endl;
            cout << "operaMasterComparison: badpixelmask = " << badpixelmask << endl;
            cout << "operaMasterComparison: masterbias = " << masterbias << endl;
 			cout << "operaMasterComparison: combineMethod = " << combineMethod << endl;
 			cout << "operaMasterComparison: saturationLimit = " << saturationLimit << endl;
 			cout << "operaMasterComparison: outputExposureTime = " << outputExposureTime << endl;
 			cout << "operaMasterComparison: biasConstant = " << biasConstant << endl;
 			cout << "operaMasterComparison: truncateOuputFluxToSaturation = " << truncateOuputFluxToSaturation << endl;
			cout << "operaMasterComparison: expTimeFITSKeyword = " << expTimeFITSKeyword << endl;
			cout << "operaMasterComparison: OPERA version= " << version << endl;
			cout << "operaMasterComparison: Reduction date= " << date << endl;
			for(unsigned i=0;i<imageIndex;i++) {
                cout << "operaMasterComparison: image[" << i << "] = " << images[i] << endl;
            }
        }
        
		// we need a masterbias for methods 2 or 3
		if (masterbias.empty() && (combineMethod == 2 || combineMethod == 3)) {
			throw operaException("operaMasterComparison: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
		if (imageIndex == 0) {
			throw operaException("operaMasterComparison: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (output.empty()) {
			throw operaException("operaMasterComparison: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		// if there are not enough images to median combine, just pick one
		if ((imageIndex < 3 && pick == 0) || (imageIndex < pick)) {
			pick = 1;	
			if (args.verbose) cout << "operaMasterComparison: too few images (" << imageIndex << "), picking " << pick << endl;
		}		
		
		// choose one input image -- one-based indexing
		if (pick) {
			if (pick > imageIndex) {
				throw (operaErrorPickOutofRange);
			}
			if (args.verbose) cout << "operaMasterComparison: picking " << images[pick-1] << endl;
			operaFITSImage *compIn = new operaFITSImage(images[pick-1], READONLY);
			operaFITSImage *masterComparison = new operaFITSImage(output, compIn->getnaxis1(), compIn->getnaxis2(), tushort, compression);
			masterComparison->operaFITSImageCopyHeader(compIn);
			*masterComparison = *compIn;	// copy the pixels
			masterComparison->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
			masterComparison->operaFITSAddComment(version);
			masterComparison->operaFITSAddComment("Picking a single comparison "+images[pick-1]);
			masterComparison->operaFITSSetHeaderValue("FILENAME", output, "Filename");
			masterComparison->operaFITSDeleteHeaderKey("EXPNUM");
			masterComparison->operaFITSDeleteHeaderKey("OBSID");
			masterComparison->operaFITSImageSave();
			masterComparison->operaFITSImageClose();
			delete masterComparison;
			delete compIn;
			return EXIT_SUCCESS;
		}
		
		if (args.verbose) {
			switch (combineMethod) {
				case 1:
					cout << "operaMasterComparison: median of " << imageIndex << " images." << endl;
					break;
				case 2:
					cout << "operaMasterComparison: etime-based combination of " << imageIndex << " images." << endl;
					break;
				case 3:
					cout << "operaMasterComparison: sum ignoring saturated pixels of " << imageIndex << " images." << endl;
					break;
				default:
					break;
			}
		}

        string exposuretimes[MAXIMAGES];
        operaFITSImage *masterComparison = NULL;
        operaFITSImage *badpix = NULL;
        operaFITSImage *bias = NULL;
        unsigned short *masterData = NULL;
        unsigned short *comparisons[MAXIMAGES];        
        unsigned short *biasData = NULL;        
        unsigned short *badpixData = NULL;
        float *exptimePerPixel = NULL;
        unsigned short biasConstantValue = 0;
        long npixels = 0;
        
        for (unsigned i=0; i<imageIndex; i++) {
            operaFITSImage *comparisonIn = new operaFITSImage(images[i], READONLY);

            if (i == 0) {
                masterComparison = new operaFITSImage(output, comparisonIn->getnaxis1(), comparisonIn->getnaxis2(), tushort, compression);
                masterComparison->operaFITSImageCopyHeader(comparisonIn);
                masterData = (unsigned short *)masterComparison->getpixels();
                npixels = masterComparison->getnpixels();
                if(args.verbose) cout << "operaMasterComparison: base image =" << images[i]  << " NAXIS1=" << comparisonIn->getnaxis1() << " NAXIS2=" << comparisonIn->getnaxis2()<< " npixels=" << npixels << endl;

                if (!masterbias.empty()) bias = new operaFITSImage(masterbias, tushort, READONLY);
                else bias = new operaFITSImage(comparisonIn->getnaxis1(),comparisonIn->getnaxis2(),tushort);
                biasData = (unsigned short *)bias->getpixels();
                if(biasConstant) biasConstantValue = operaArrayMedianQuickUSHORT(bias->getnpixels(),bias->operaFITSImageClonePixelsUSHORT());
                if(args.debug) cerr << "operaMasterComparison: biasConstantValue=" << biasConstantValue << endl;
                
                if (!badpixelmask.empty()) badpix = new operaFITSImage(badpixelmask, tushort, READONLY);
                else badpix = new operaFITSImage(comparisonIn->getnaxis1(),comparisonIn->getnaxis2(),tushort);
                badpixData = (unsigned short *)badpix->getpixels();
                
                if(combineMethod==2 || combineMethod==3) {
                    exptimePerPixel = new float[npixels];
                    for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
                        if(badpixelmask.empty()) badpixData[pixIndex] = 1;
                        exptimePerPixel[pixIndex] = 0.0;
                        if(biasConstant) masterData[pixIndex] = biasConstantValue;
                        else masterData[pixIndex] = biasData[pixIndex];
                    }
                }
            }
            comparisons[i] = comparisonIn->operaFITSImageClonePixelsUSHORT();
            
            if(combineMethod==2 || combineMethod==3) {                
                exposuretimes[i] = comparisonIn->operaFITSGetHeaderValue(expTimeFITSKeyword.c_str());
                float exptime = atof(exposuretimes[i].c_str());
                if(args.verbose) cout << "operaMasterComparison: image= " << images[i] << " -> " << expTimeFITSKeyword << "= " << exptime << endl;
                
                for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
                    if(args.debug) cerr << "operaMasterComparison: badpixData[pixIndex]=" << badpixData[pixIndex] << " comparisons[i][pixIndex]=" << comparisons[i][pixIndex] << " saturationLimit=" << saturationLimit << " npixels=" << npixels << " exptimePerPixel[" <<pixIndex<<"]=" << exptimePerPixel[pixIndex] << " exptime=" << exptime << " masterData[pixIndex]=" << masterData[pixIndex] << endl;
                    
                    if(comparisons[i][pixIndex] < saturationLimit && badpixData[pixIndex] == 1) {
                        masterData[pixIndex] += int(comparisons[i][pixIndex]) - int(biasData[pixIndex]); // WARNING: if *master is not initialized as *bias or biasConstant this could be negative
                        exptimePerPixel[pixIndex] += exptime;
                    }
                }
            }
            comparisonIn->operaFITSImageClose();
            delete comparisonIn;
        }
        comparisons[imageIndex] = NULL;

		masterComparison->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
		masterComparison->operaFITSAddComment(version);
        if(combineMethod==1) {
            masterData = operaArrayMedianCombineUSHORT(imageIndex, npixels, masterData, comparisons);
            masterComparison->operaFITSAddComment("A median of "+itos(imageIndex)+" images.");
            for (unsigned i=0; i<imageIndex; i++) {
                masterComparison->operaFITSAddComment("Using comparison image "+images[i]);
            }
        } else if(combineMethod==2 || combineMethod==3) {
            for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
				if(args.debug && combineMethod==2) cerr << "operaMasterComparison: exptimePerPixel[pixIndex]=" << exptimePerPixel[pixIndex] << " masterData[pixIndex]=" << masterData[pixIndex] << endl;
                
                if(exptimePerPixel[pixIndex] > 0) {
					if(combineMethod==2) {
						float biasval = biasConstant ? float(biasConstantValue) : float(biasData[pixIndex]);
						masterData[pixIndex] = (unsigned)((float(masterData[pixIndex]) - biasval) * float(outputExposureTime)/exptimePerPixel[pixIndex] + biasval);
					}
                    if(masterData[pixIndex] > saturationLimit && truncateOuputFluxToSaturation) masterData[pixIndex] = saturationLimit;
                } else {
                    masterData[pixIndex] = saturationLimit;
                }
            }
            masterComparison->operaFITSAddComment("A stack of "+itos(imageIndex)+" images.");
            if (combineMethod==2) masterComparison->operaFITSAddComment("Combined flux is equivalent to an exptime of "+ftos(float(outputExposureTime))+" s");
            if (combineMethod==3) masterComparison->operaFITSAddComment("Combined flux by the sum of all images");
            for (unsigned i=0; i<imageIndex; i++) {
                masterComparison->operaFITSAddComment("Using comparison image "+images[i]+" exptime="+exposuretimes[i]+" s");
            }
        } else {
            throw operaException("operaMasterComparison: combineMethod not supported. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        masterComparison->operaFITSSetHeaderValue("FILENAME", output, "Filename");
        masterComparison->operaFITSDeleteHeaderKey("EXPNUM");
        masterComparison->operaFITSDeleteHeaderKey("OBSID");
        masterComparison->operaFITSImageSave();
        masterComparison->operaFITSImageClose();

        delete masterComparison;
        if(badpix) delete badpix;
        if(bias) delete bias;
	}
	catch (operaException e) {
		cerr << "operaMasterComparison: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMasterComparison: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
