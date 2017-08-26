/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaExtractionApertureTest
 Version: 1.0
 Description: Perform various tests on the operaGeometricShapes class.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
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
#include <libgen.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaExtractionAperture.h"
#include "libraries/operaGeometricShapes.h"
#include "libraries/operaFITSImage.h"

static void printUsageSyntax();

/*! \file operaExtractionApertureTest.cpp */

using namespace std;

/*! 
 * operaExtractionApertureTest
 * \author Eder Martioli
 * \brief Perform various tests on the operaExtractionAperture class.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(int argc, char *argv[])
{	
    int opt;
	string input;
	string output;
	string hiresoutput;
	
	bool debug=false, verbose=false, trace=false;
	
	struct option longopts[] = {
		{"input",			1, NULL, 'i'},	// input image	
		{"output",			1, NULL, 'o'},	//  output image
		{"hiresoutput",		1, NULL, 'r'},	//  hi-resolution output image        
		{"verbose",			0, NULL, 'v'},
		{"debug",			0, NULL, 'd'},
		{"trace",			0, NULL, 't'},
		{"help",			0, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:o:r:vdth", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// image
					input = optarg;
					break;				
				case 'o':		// output
					output = optarg;
					break;
				case 'r':		// output
					hiresoutput = optarg;
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
			}
		}
		
		/* 
		 * Read image path from input file and append to list of input images	
		 */ 
		if (input.empty()){ 
			throw operaException("operaExtractionAperture: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}	        
		operaFITSImage inImage(input, tfloat, READONLY);					
		operaFITSImage outImage(output, inImage.getnaxis1(), inImage.getnaxis2(), tfloat, 0);		
        
		outImage.operaFITSImageCopyHeader(&inImage);		
        
        
        cout << "Test Module for operaExtractionAperture" << endl << endl;        
        
        /*
         * Test aperture with Circle
         */         
        operaPoint circleCenter(950,2300);
        float radius = 300;
        Circle *mycircle = new Circle(radius, circleCenter);
        
        /*
         * Test aperture with Rectangle
         */   
        operaPoint rectangleCenter(575,3900);
        float width = 60, height = 800, angle = 30;        
        Rectangle *myrectangle = new Rectangle(width,height,angle,rectangleCenter);
        
        /*
         * Test aperture with Polygon
         */        
        unsigned nsides = 4;
        operaPoint polygoncorners[4];
        polygoncorners[0].setPoint(430,243);
        polygoncorners[1].setPoint(341,342);
        polygoncorners[2].setPoint(333,1000);
        polygoncorners[3].setPoint(1111,1053);

        Polygon *myPolygon  = new Polygon(nsides,polygoncorners); 

        /*
         * Test aperture with Lines
         */  
        float lineSlope = -0.54;
        float lineWidth = 100;
        float lineLength = 900;
        operaPoint lineMidPoint(1000,1000);
        Line *myline = new Line(lineSlope,lineWidth,lineLength,lineMidPoint);  
    
        /*
         * Tests for operaExtractionAperture starts below
         */
        
        operaExtractionAperture *circleAperture = new operaExtractionAperture(mycircle,1,1,inImage);
        operaExtractionAperture *rectangleAperture = new operaExtractionAperture(myrectangle,1,1,inImage);
        operaExtractionAperture *polygonAperture = new operaExtractionAperture(myPolygon,1,1,inImage);
        operaExtractionAperture *lineAperture = new operaExtractionAperture(myline,1,1,inImage);        
        
        /*
         * Grab pixel set and writes out only pixels within extraction aperture
         */        
        PixelSet *circpixels = circleAperture->getSubpixels();
        for(unsigned j=0;j<5;j++){
            circleAperture->shiftAperture(50*float(j),80*float(j),inImage);  // apply shifts to aperture
            circpixels = circleAperture->getSubpixels();              
            for(unsigned i=0; i<circpixels->getNPixels(); i++){
                outImage[circpixels->getjIndex(i)][circpixels->getiIndex(i)] = circpixels->getPixelValue(i);
            }        
        }
        
        PixelSet *recpixels = rectangleAperture->getSubpixels();
        for(unsigned i=0; i<recpixels->getNPixels(); i++){
            outImage[recpixels->getjIndex(i)][recpixels->getiIndex(i)] = recpixels->getPixelValue(i);
        }
        
        operaPoint newCenter(275,3900);
        rectangleAperture->recenterAperture(newCenter,inImage);   // recenter aperture
        recpixels = rectangleAperture->getSubpixels();
        for(unsigned i=0; i<recpixels->getNPixels(); i++){
            outImage[recpixels->getjIndex(i)][recpixels->getiIndex(i)] = recpixels->getPixelValue(i);
        }
        
        PixelSet *polypixels = polygonAperture->getSubpixels();        
        for(unsigned i=0; i<polypixels->getNPixels(); i++){
            outImage[polypixels->getjIndex(i)][polypixels->getiIndex(i)] = polypixels->getPixelValue(i);
        }
        
        PixelSet *linepixels = lineAperture->getSubpixels();        
        for(unsigned i=0; i<linepixels->getNPixels(); i++){
            outImage[linepixels->getjIndex(i)][linepixels->getiIndex(i)] = linepixels->getPixelValue(i);
        }        
        /*
         * Tests for high-resolution image starts below
         */
/*        
        unsigned xsampling = 5;
        unsigned ysampling = 5;
        
        operaExtractionAperture *hrcircleAperture = new operaExtractionAperture(mycircle,xsampling,ysampling,inImage);
        
        unsigned naxis1 = ceil(hrcircleAperture->getBoundingBox()->getWidth())*xsampling;
        unsigned naxis2 = ceil(hrcircleAperture->getBoundingBox()->getHeight())*ysampling;
        
		operaFITSImage hiresoutImage(hiresoutput,naxis1,naxis2,tfloat,0);		
//        hiresoutImage.operaFITSImageCopyHeader(&inImage); 
        
        PixelSet *hrcircpixels = hrcircleAperture->getSubpixels();
            for(unsigned i=0; i<hrcircpixels->getNPixels(); i++){
                hiresoutImage[hrcircpixels->getjIndex(i)][hrcircpixels->getiIndex(i)] = hrcircpixels->getPixelValue(i);
            }        
        
        delete hrcircleAperture;
 
        hiresoutImage.operaFITSImageSave();
        hiresoutImage.operaFITSImageClose();   
  */   
        
        delete circleAperture;
        delete rectangleAperture;
        delete polygonAperture;
        delete lineAperture;
        
        delete mycircle;
        delete myrectangle;
        delete myPolygon;
        delete myline;        
        
		
        outImage.operaFITSImageSave();
		outImage.operaFITSImageClose();      
        
		inImage.operaFITSImageClose();
 	}
	catch (operaException e) {
		cerr << "operaExtractionAperture: " << e.getFormattedMessage() << '\n';
	}
	catch (...) {
		cerr << "operaExtractionAperture: " << operaStrError(errno) << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cerr << " Usage: operaExtractionApertureTest [--input=<FITS file name>] --output=<file name> -[dvth]\n";
}	
