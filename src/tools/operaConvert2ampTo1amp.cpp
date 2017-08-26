/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaConvert2ampTo1amp
 Version: 1.0
 Description: THis module converts 2amp images into 1amp
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

/*
 * Algorithm
 *
 *
 */

#include <sys/stat.h>					// mkdir
#include <getopt.h>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "tools/operaConvert2ampTo1amp.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"

/*! \file operaConvert2ampTo1amp.cpp */

using namespace std;

/*! 
 * operaConvert2ampTo1amp
 * \author Eder Martioli / Doug Teeple
 * \brief Convert images taken with 2 amplifiers into the same format as images taken with only 1 amplifier.
 * \arg argc
 * \arg argv
 * \note [--images=...]*
 * \note 
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOutput
 * \ingroup tools
 * \return EXIT_STATUS
 */
/*
 * Here are the header differences...
 *  OLAPA-a
 *  NAXIS   =                    2 / Number of axes
 *  NAXIS1  =                 2080 / Number of pixel columns
 *  NAXIS2  =                 4640 / Number of pixel rows
 *  DETSIZE = '[1:2048,1:4608]'    / Total data pixels in full mosaic
 *  CCDSIZE = '[1:2048,1:4608]'    / Detector imaging area size
 *  CCDSEC  = '[1:2048,1:4608]'    / Read out area of the detector (unbinned)
 *  DATASEC = '[1:2048,1:4608]'    / Imaging area of the detector
 *  BIASSEC = '[2049:2080,1:4608]' / Overscan (bias) area of the detector
 *  
 *  OLAPA-ab
 *  NAXIS   =                    2 / Number of axes
 *  NAXIS1  =                 2088 / Number of pixel columns
 *  NAXIS2  =                 4608 / Number of pixel rows
 *  DETSIZE = '[1:2048,1:4608]'    / Total data pixels in full mosaic
 *  CCDSIZE = '[1:2048,1:4608]'    / Detector imaging area size
 *  CCDSEC  = '[21:2068,1:4608]'   / Read out area of the detector (unbinned)
 *  DATASEC = '[21:2068,1:4608]'   / Imaging area of the detector
 *  TRIMSEC = '[21:2068,4:4605]'   / Useful imaging area of the detector
 *  BSECA   = '[1:20,1:4608]'      / Overscan (bias) area from Amp A
 *  BSECB   = '[2069:2088,1:4608]' / Overscan (bias) area from Amp B
 *  CSECA   = '[21:1044,1:4608]'   / Section in full CCD for DSECA
 *  CSECB   = '[1045:2068,1:4608]' / Section in full CCD for DSECB
 *  DSECA   = '[21:1044,1:4608]'   / Imaging area from Amp A
 *  DSECB   = '[1045:2068,1:4608]' / Imaging area from Amp B
 *  TSECA   = '[21:1044,4:4605]'   / Trim section for Amp A
 *  TSECB   = '[1045:2068,4:4605]' / Trim section for Amp B
 *
 */

int main(int argc, char *argv[])
{
	int opt;
	string images[MAXIMAGES];
	string listofimages;
	string outputdir;
	unsigned imageIndex = 0;
	unsigned OLAPA_ab_overscan = 12;
	const unsigned AmpBStartCol = 1045;
	const unsigned NAXIS1 = 2080;
	const unsigned NAXIS2 = 4640;
	const unsigned CCDX = 2048;
	const unsigned CCDY = 4608;
	float gainA = 1.0;
	float gainB = 1.0;
	unsigned biasA = 0;
	unsigned biasB = 0;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"images",			1, NULL, 'i'},	// series of input flats
		{"list",			1, NULL, 'l'},	// list of input flats		
		{"outputdir",		1, NULL, 'o'},	// directory for output		
		{"gaina",			1, NULL, 'A'},	// gain amp A		
		{"gainb",			1, NULL, 'B'},	// gain amp B		
		{"biasa",			1, NULL, 'a'},	// bias amp A		
		{"biasb",			1, NULL, 'b'},	// bias amp B		
		{"noshift",			0, NULL, 's'},	// subtract bias don't shift		
		
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:l:o:a:b:svdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// images
					images[imageIndex++] = optarg;
					break;
				case 'l':		// list of images
					listofimages = optarg;
					break;	
				case 'o':		// directory for output
					outputdir = optarg;
					break;						
				case 'A':		// gain amp A
					gainA = atof(optarg);
					break;						
				case 'B':		// gain amp A
					gainB = atof(optarg);
					break;						
				case 'a':		// bias amp A
					biasA = atoi(optarg);
					break;						
				case 'b':		// bias amp A
					biasB = atoi(optarg);
					break;						
				case 's':		// no pixel shift
					OLAPA_ab_overscan = 0;
					break;						
					
				case 'v':
					verbose = true;
					break;
				case 'd':
					debug = true;
					break;
				case 't':
					trace = true; 
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
		
		/* 
		 * Read image path from input file and append to list of input images	
		 */ 
		if (!listofimages.empty()) {
			ifstream flist(listofimages.c_str());		
			if (flist.is_open())
			{
				while (flist.good()) {
					getline (flist,images[imageIndex++]);
					if (images[imageIndex-1].size() == 0 || images[imageIndex-1][0] == '#')
						imageIndex--;					
				}	
				flist.close();
			}
		}
		if (debug) {
			for(unsigned i=0;i<imageIndex;i++)
				cout << "Image #" << i << ": " << images[i]  << endl;
		}
		
		/*
		 * end of reading list of images
		 */
		if (imageIndex == 0) {
			throw operaException("operaConvert2ampTo1amp: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (outputdir.empty()) {
			throw operaException("operaConvert2ampTo1amp: please specify an output directory ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (gainB == 0.0 || gainA == 0.0) {
			throw operaException("operaConvert2ampTo1amp: neither gainB nor gainA may be zero ", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		// make sure output directory exists
		struct stat St;
		if (stat(outputdir.c_str(), &St ) != 0) {
			mkdir(outputdir.c_str(), 0777);
		}
		
		bool biasAIsHigher = (biasA > biasB);
		bool biasBIsHigher = (biasA < biasB);
		unsigned short biasDiff = 0;
		if (biasAIsHigher) {
			biasDiff = biasA - biasB;
		} else {
			biasDiff = biasB - biasA;
		}
		float maxGain = (gainA>gainB?gainA:gainB);
		float minGain = (gainA<gainB?gainA:gainB);
		
		
		for (unsigned i=0; i<imageIndex; i++) {
			
			string basefilename = images[i].substr(images[i].find_last_of('/'));
			string output = outputdir + basefilename;			
			if(verbose){
				cout << "operaConvert2ampTo1amp: Converting image: " << images[i] << " --> " << output << endl;					
			}
			
			// remove existing file - cfitsio returns an error if it exists...
			remove(output.c_str());			
			
			operaFITSImage outputImage(output, NAXIS1, NAXIS2, tushort, 0);
			operaFITSImage inputImage(images[i], READONLY);
			
			if (verbose) {
				if (gainA > gainB) {	// gainA > gainB, lower amp A
					cout << "operaConvert2ampTo1amp: gain of amp A is greater than amp B. ampA=" << gainA << " ampB=" << gainB  << endl;					
				} else if (gainA < gainB){
					cout << "operaConvert2ampTo1amp: gain of amp B is greater than amp A. ampA=" << gainA << " ampB=" << gainB  << endl;					
				}
				if (OLAPA_ab_overscan == 0) {
					cout << "operaConvert2ampTo1amp: no pixel shift applied." << endl;					
				} else {
					cout << "operaConvert2ampTo1amp: pixel shift of " << OLAPA_ab_overscan << " applied." << endl;					
				}
			}
			
			for (unsigned x=OLAPA_ab_overscan; x<CCDX+OLAPA_ab_overscan; x++) {
				if (gainA > gainB) {	// gainA > gainB, lower amp A
					if (x >= AmpBStartCol) {	// ampB
						for (unsigned y=0; y<CCDY; y++) {
							unsigned short pixelvalue = inputImage.getpixelUSHORT(x, y);
							if (biasBIsHigher) {
								if (pixelvalue != 65535) {	// retain saturation
									pixelvalue = pixelvalue - biasDiff;
								}
							}
							outputImage.setpixel(pixelvalue, x-OLAPA_ab_overscan, y);
						}
					} else {					// ampA
						for (unsigned y=0; y<CCDY; y++) {
							unsigned short pixelvalue = inputImage.getpixelUSHORT(x, y);
							if (pixelvalue != 65535) {	// retain saturation
								if (biasAIsHigher) {
									pixelvalue = pixelvalue - biasDiff;
								}
								pixelvalue = (unsigned short)((float)pixelvalue*minGain/maxGain);
							}
							outputImage.setpixel(pixelvalue, x-OLAPA_ab_overscan, y);
						}
						
					}
				} else if (gainA < gainB) {	// gainB > gainA, lower amp B
					if (x < AmpBStartCol) {	// ampA
						for (unsigned y=0; y<CCDY; y++) {
							unsigned short pixelvalue = inputImage.getpixelUSHORT(x, y);
							if (biasAIsHigher) {
								if (pixelvalue != 65535) {	// retain saturation
									pixelvalue = pixelvalue - biasDiff;
								}
							}
							outputImage.setpixel(pixelvalue, x-OLAPA_ab_overscan, y);
						}
					} else {				// ampB
						for (unsigned y=0; y<CCDY; y++) {
							unsigned short pixelvalue = inputImage.getpixelUSHORT(x, y);
							if (pixelvalue != 65535) {	// retain saturation
								if (biasBIsHigher) {
									pixelvalue = pixelvalue - biasDiff;
								}
								pixelvalue = (unsigned short)((float)pixelvalue*minGain/maxGain);
							}
							outputImage.setpixel(pixelvalue, x-OLAPA_ab_overscan, y);
						}
						
					}
				} else {	// gains are the same
					if (x < AmpBStartCol) {	// ampA
						for (unsigned y=0; y<CCDY; y++) {
							unsigned short pixelvalue = inputImage.getpixelUSHORT(x, y);
							if (biasAIsHigher) {
								if (pixelvalue != 65535) {	// retain saturation
									pixelvalue = pixelvalue - biasDiff;
								}
							}
							outputImage.setpixel(pixelvalue, x-OLAPA_ab_overscan, y);
						}
					} else {	// ampB
						for (unsigned y=0; y<CCDY; y++) {
							unsigned short pixelvalue = inputImage.getpixelUSHORT(x, y);
							if (biasBIsHigher) {
								if (pixelvalue != 65535) {	// retain saturation
									pixelvalue = pixelvalue - biasDiff;
								}
							}
							outputImage.setpixel(pixelvalue, x-OLAPA_ab_overscan, y);
						}
					}
				}
			}
			outputImage.operaFITSImageCopyHeader(&inputImage);
			// setheader overwrites the existing keyword
			outputImage.operaFITSSetHeaderValue("DETSIZE", "[1:2048,1:4608]", "Total data pixels in full mosaic");
			outputImage.operaFITSSetHeaderValue("CCDSIZE", "[1:2048,1:4608]", "Detector imaging area size");
			outputImage.operaFITSSetHeaderValue("CCDSEC",  "[1:2048,1:4608]", "Read out area of the detector (unbinned)");
			outputImage.operaFITSSetHeaderValue("DATASEC", "[1:2048,1:4608]", "Imaging area of the detector");
			outputImage.operaFITSSetHeaderValue("BIASSEC", "[2049:2080,1:4608]", "Overscan (bias) area of the detector");
			outputImage.operaFITSDeleteHeaderKey("TRIMSEC");	// This got added to OLAPA-ab...
			outputImage.operaFITSAddComment("NOTE: This image has been modified by a preprocessing step.");
			outputImage.operaFITSAddComment("NOTE: The image has been shifted 12 pixels to the left");
			outputImage.operaFITSAddComment("NOTE: to remove the OLAPA-ab overscan.");
			outputImage.operaFITSImageSave();
			outputImage.operaFITSImageClose();
			inputImage.operaFITSImageClose();				
			
		}
		
	}
	catch (operaException e) {
		cerr << "operaConvert2ampTo1amp: " << e.getFormattedMessage() << '\n';
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaConvert2ampTo1amp: " << operaStrError(errno) << "\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
	
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout <<
	"\n"
	" Usage: operaConvert2ampTo1amp [--images=<filenames>] [--list=<list of files>] -[dvth]\n";
}	

