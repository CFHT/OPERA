/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMedianCombine
 Version: 1.0
 Description: This module median combines a list of images
 Author(s): CFHT OPERA team / Eder Martioli
 Affiliation: Canada France Hawaii Telescope / Laboratorio Nacional de Astrofisica
 Location: Hawaii USA / Itajuba-MG Brazil
 Date: Jan/2014
 Contact: opera@cfht.hawaii.edu / emartioli@lna.br
 
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
 * Pick one or median stack images.
 *
 */

#include <getopt.h>
#include <iostream>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "tools/operaMedianCombine.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaLib.h"     // for itos

#define MAXDIRNAMESIZE 1000

#include "libraries/operaStats.h"
#include "libraries/operaImage.h"

/*! \file operaMedianCombine.cpp */

using namespace std;

/*! 
 * operaMedianCombine
 * \author Eder Martioli
 * \brief This module median combines a list of images
 * \arg argc
 * \arg argv
 * \note operaMedianCombine [--images=...]* --output=...[ --pick=\<posint\>0\> ]
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \ingroup core
 * \return EXIT_STATUS
 */
int main(int argc, char *argv[])
{
	int opt;
	string badpixelmask;
	string images[MAXIMAGES];
	string listofimages;
	string output;
	string version = "OPERA-1.0";
	string date = "";
    
    bool normalize = false;
    
	unsigned imageIndex = 0;
	eCompression compression = cNone;
	
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"images",			1, NULL, 'i'},	// series of input imgs
		{"list",			1, NULL, 'l'},	// list of input imgs
		{"output",			1, NULL, 'o'},	// a single master img output fits file
		{"badpixelmask",	1, NULL, 'm'},	// bad pixel mask fits file to ignore pixels
		{"compressiontype", 1, NULL, 'C'},
		{"version",			1, NULL, 'V'},
		{"date",			1, NULL, 'a'},
		{"normalize",		1, NULL, 'n'},
		
		{"plot",			optional_argument, NULL, 'p'},
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:l:o:m:C:V:a:n:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// images
					images[imageIndex++] = optarg;
					break;
				case 'l':		// list of images
					listofimages = optarg;
					break;
				case 'o':		// output
					output = optarg;
					break;
				case 'm':		// badpixelmask
					badpixelmask = optarg;
					break;
				case 'C':
					compression = (eCompression)atoi(optarg);
					break;
				case 'V':
					version = optarg;
					break;
				case 'a':
					date = optarg;
					break;
				case 'n':
					normalize = (atoi(optarg)?true:false);
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
				cout << images[i] << " " << i << endl;
		}
		/*
		 * end of reading list of images
		 */
		if (imageIndex == 0) {
			throw operaException("operaMedianCombine: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (output.empty()) {
			throw operaException("operaMedianCombine: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		// if there are not enough images to median combine, just pick one
		if (imageIndex < 1) {
			if (verbose)
				cout << "operaMedianCombine: too few images (" << imageIndex << ")" <<  endl;
		}
		if (verbose) {
			cout << "operaMedianCombine: output= " << output << endl;
			cout << "operaMedianCombine: imageIndex= " << imageIndex << endl;
		}
        
		/*
		 * else do a median
		 */
		long npixels = 0;
		unsigned i;
		operaFITSImage *masterFlat = NULL;
        float *flats[MAXIMAGES];
		float *masterData = NULL;
		for (i=0; i<imageIndex; i++) {
            operaFITSImage inputImage(images[i], tfloat, READONLY, cNone, false);
			if (i == 0) {
				masterFlat = new operaFITSImage(output, inputImage.getnaxis1(), inputImage.getnaxis2(), tfloat, compression);
				masterFlat->operaFITSImageCopyHeader(&inputImage);
				masterData = (float *)masterFlat->getpixels();
				npixels = masterFlat->getnpixels();
			}
			flats[i] = inputImage.operaFITSImageClonePixels();
		}
		flats[i] = NULL;
		masterData = medianCombineFloat(imageIndex, npixels, masterData, flats);
		masterFlat->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
		masterFlat->operaFITSAddComment(version);
		masterFlat->operaFITSAddComment("A median of "+itos(imageIndex)+" images.");
		for (i=0; i<imageIndex; i++) {
			masterFlat->operaFITSAddComment("Using image "+images[i]);
		}
		masterFlat->operaFITSSetHeaderValue("FILENAME", output, "Filename");
        
		masterFlat->operaFITSImageSave();
        
        if(normalize) {
            float *pixdata = new float[masterFlat->getnaxis1()*masterFlat->getnaxis2()];
            unsigned np = 0;
            for(unsigned x=0; x<masterFlat->getnaxis1(); x++) {
                for(unsigned y=0; y<masterFlat->getnaxis2(); y++) {
                    if((*masterFlat)[y][x] && !isnan((*masterFlat)[y][x])) {
                        pixdata[np++] = (*masterFlat)[y][x];
                    }
                }
            }
            float medianpix = operaArrayMedian(np,pixdata);
            for(unsigned x=0; x<masterFlat->getnaxis1(); x++) {
                for(unsigned y=0; y<masterFlat->getnaxis2(); y++) {
                    (*masterFlat)[y][x] = (*masterFlat)[y][x]/medianpix;
                }
            }
            masterFlat->operaFITSImageSave();
            delete[] pixdata;
        }
        
		masterFlat->operaFITSImageClose();
		delete masterFlat;
	}
	catch (operaException e) {
		cerr << "operaMedianCombine: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMedianCombine: " << operaStrError(errno) << endl;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout <<
	"\n"
	" Usage: operaMedianCombine [--images=<flat filename>]+ --output=<master flat file name> [badpixlemask=<bad pixel mask file name>] [--pick=<n>] -[dvth]\n";

}	
