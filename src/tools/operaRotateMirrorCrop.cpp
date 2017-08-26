/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaRotateMirrorCrop
 Version: 1.0
 Description: Module to rotate, mirror and crop a list of images
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
 * Rotate, mirror or/and crop a list of images and save them
 *
 */

#include <getopt.h>
#include <iostream>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "tools/operaRotateMirrorCrop.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaLib.h"     // for itos

#define MAXDIRNAMESIZE 1000

#include "libraries/operaStats.h"
#include "libraries/operaImage.h"

/*! \file operaRotateMirrorCrop.cpp */

using namespace std;

/*! 
 * operaRotateMirrorCrop
 * \author Eder Martioli
 * \brief Module to rotate, mirror and crop a list of images
 * \arg argc
 * \arg argv
 * \note operaRotateMirrorCrop [--images=...]*
 * \note Pick one or median stack images.
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \ingroup core
 * \return EXIT_STATUS
 */
int main(int argc, char *argv[])
{
	int opt;
	string images[MAXIMAGES];
	string listofimages;
	string outputdir;
	string sufix;
	string version = "OPERAOES-1.0";
	string date = "";
	unsigned imageIndex = 0;
	eCompression compression = cNone;
	
    bool rotate = false;
	
    bool mirrorcols = false;
    bool mirrorrows = false;
    bool crop = false;
    
    struct subwindow {
		unsigned x0,xf;
		unsigned y0,yf;
	} cropsubwindow = {0,0,0,0};
    
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"images",			1, NULL, 'i'},	// series of input flats
		{"listofimages",    1, NULL, 'l'},	// list of input flats		
		{"outputdir",       1, NULL, 'o'},	// output directory
		{"sufix",           1, NULL, 'S'},	// sufix to add at the end of file names
		{"rotate",          1, NULL, 'R'},	// rotate output by 90 degrees
		{"mirrorcols",		1, NULL, 'c'},	// mirror columns (x-axis)
		{"mirrorrows",		1, NULL, 'r'},	// mirror rows (y-axis)
		{"cropsubwindow",	1, NULL, 's'},	// crop subwindow
		{"compressiontype", 1, NULL, 'C'},
		{"version",			1, NULL, 'V'},
		{"date",			1, NULL, 'a'},
		
		{"plot",			optional_argument, NULL, 'p'},       
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:l:o:k:m:r:C:V:a:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// images
					images[imageIndex++] = optarg;
					break;
				case 'l':		// list of images
					listofimages = optarg;
					break;					
				case 'o':		// output
					outputdir = optarg;
					break;
				case 'S':		// output
					sufix = optarg;
					break;
				case 'R':		// rotate by 90 degrees
					rotate = (atoi(optarg)?true:false);
					break;
				case 'c':		
					mirrorcols = (atoi(optarg)?true:false);
					break;
				case 'r':	
					mirrorrows = (atoi(optarg)?true:false);
					break;
				case 's':		
                    if (strlen(optarg))
                        sscanf(optarg, "%u %u %u %u", &cropsubwindow.x0, &cropsubwindow.xf, &cropsubwindow.y0, &cropsubwindow.yf);
                    crop = true;
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
					if(images[imageIndex-1].size() == 0 || images[imageIndex-1][0] == '#')
						imageIndex--;					
				}	
				flist.close();
			}
		}
        
        
        string outputImages[MAXIMAGES];
        
        for(unsigned i=0;i<imageIndex;i++) {
            
            string imageSuffix =  images[i].substr(images[i].find_last_of("."));
                        
            if(!outputdir.empty()) {
                string basefilename = images[i];
                if (images[i].find_last_of("/") != string::npos) {
                    basefilename = images[i].substr(images[i].find_last_of("/")+1);
                }
                
                if (imageSuffix.compare(sufix) != 0) {
                    outputImages[i] = outputdir + basefilename + sufix;
                } else {
                    outputImages[i] = outputdir + basefilename;
                }
                
            } else {
                if (imageSuffix.compare(sufix) != 0) {
                    outputImages[i] = images[i];
                } else {
                    outputImages[i] = images[i] + sufix;
                }
            }
        }

		/*
		 * end of reading list of images
		 */
		if (imageIndex == 0) {
			throw operaException("operaRotateMirrorCrop: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}

		if (verbose) {
			cout << "operaRotateMirrorCrop: imageIndex= " << imageIndex << endl;
			cout << "operaRotateMirrorCrop: outputdir= " << outputdir << endl;
            cout << "operaRotateMirrorCrop: sufix= " << sufix << endl;
            cout << "operaRotateMirrorCrop: rotate= " << rotate << endl;
			cout << "operaRotateMirrorCrop: mirrorcols= " << mirrorcols << endl;
			cout << "operaRotateMirrorCrop: mirrorrows= " << mirrorrows << endl;
            if(crop) {
                cout << "operaRotateMirrorCrop: cropsubwindow = \"" << cropsubwindow.x0 << " " << cropsubwindow.xf << " " << cropsubwindow.y0 << " " << cropsubwindow.yf << "\"" << endl;
            }
			cout << "operaRotateMirrorCrop: OPERA version= " << version << endl;
			cout << "operaRotateMirrorCrop: Reduction date= " << date << endl;
		}
        
        for (unsigned i=0; i<imageIndex; i++) {
            operaFITSImage inputImage(images[i],tfloat,READONLY);
            if(verbose){
                cout << "operaRotateMirrorCrop: processing image " << i + 1 << " of " << imageIndex << " " << images[i] << " -> " << outputImages[i]<< endl;
            }
            
            operaFITSImage outputImage(outputImages[i],inputImage.getnaxis1(),inputImage.getnaxis2(),inputImage.getdatatype(),compression,false);
			outputImage.operaFITSImageCopyHeader(&inputImage);
            
            outputImage = inputImage;
            
            if (rotate) {
                if(verbose){
                    cout << "operaRotateMirrorCrop: rotating by 90 degrees ... " << endl;
                }
				outputImage.rotate90();
			}

            if (mirrorcols) {
                if(verbose){
                    cout << "operaRotateMirrorCrop: mirroring columns ... " << endl;
                }
                outputImage.mirrorColumns();
            }

            if (mirrorrows) {
                if(verbose){
                    cout << "operaRotateMirrorCrop: mirroring rows ... " << endl;
                }
                outputImage.mirrorRows();
            }
            
            if(crop) {
                if(verbose){
                    cout << "operaRotateMirrorCrop: cropping image ... " << endl;
                }
                
                if((cropsubwindow.x0 == 0 && cropsubwindow.xf == 0 && cropsubwindow.y0 == 0 && cropsubwindow.yf == 0) ||
                   cropsubwindow.x0 >= cropsubwindow.xf  || cropsubwindow.y0 >= cropsubwindow.yf ||
                   cropsubwindow.xf > outputImage.getnaxis1() || cropsubwindow.yf > outputImage.getnaxis2()) {
                    
                    cerr << "operaRotateMirrorCrop: invalid crop region. setting crop region to full image." << endl;
                    
                    cropsubwindow.x0 = 0;
                    cropsubwindow.xf = outputImage.getnaxis1();
                    cropsubwindow.y0 = 0;
                    cropsubwindow.yf = outputImage.getnaxis2();
                }
                if(cropsubwindow.xf > outputImage.getnaxis1() ) {
                    cropsubwindow.xf = outputImage.getnaxis1();
                }
                if(cropsubwindow.yf > outputImage.getnaxis2()) {
                    cropsubwindow.yf = outputImage.getnaxis2();
                }
                
                outputImage.resize(cropsubwindow.x0,cropsubwindow.xf,cropsubwindow.y0,cropsubwindow.yf);
            }
                        
            outputImage.operaFITSImageSave();
            
            outputImage.operaFITSImageClose();
            
            inputImage.operaFITSImageClose();
        }
	}
	catch (operaException e) {
		cerr << "operaRotateMirrorCrop: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaRotateMirrorCrop: " << operaStrError(errno) << endl;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout <<
	"\n"
	" Usage: operaRotateMirrorCrop [--images=<flat filename>]+ --output=<master flat file name> [badpixlemask=<bad pixel mask file name>] [--pick=<n>] -[dvth]\n";

}	
