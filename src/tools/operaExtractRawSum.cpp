/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaExtractRawSum
 Version: 1.0
 Description: Extract raw spectrum
 to start up with an OPERA module.
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

#include <stdio.h>
#include <getopt.h>

#include "fitsio.h"
#include "globaldefines.h"
#include "operaError.h"
#include "tools/operaExtractRawSum.h"

#include "libraries/operaException.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaIOFormats.h"
#include "libraries/operaInstrumentProfile.h"		// for operaInstrumentProfile
#include "libraries/operaFITSImage.h"

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"					// for MAXORDERS
#include "libraries/operaFit.h"	

#define NOTPROVIDED -999

/*! \file operaExtractRawSum.cpp */

using namespace std;

/*! 
 * operaExtractRawSum
 * \author Eder Martioli
 * \brief Tool to extract raw spectra.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup tools
 */

int main(int argc, char *argv[])
{
	int opt;
	string inputImageName; 
	string inputGeometryName; 	
	string inputWavelengthName; 	
	string inputBiasName; 	
	string inputgain;
	string outputSpectraFile; 
	string badpixelmask;	
	
	int ordernumber = NOTPROVIDED;
	float aperture = 20;	
	float spectralElementHeight = 1.0;
	//float defaultEspadonsBias = 400.0;
	
	int debug=0, verbose=0, trace=0, plot=0;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;    
    bool maxorderprovided = false; 
	
    unsigned binsize = 100;    
    bool usePolynomial = FALSE;
    unsigned orderOfPolynomial = 5;
	bool normalize = false;

	struct option longopts[] = {
		{"inputImage",				1, NULL, 'i'},
		{"inputGeometryFile",		1, NULL, 'g'},		
		{"inputBiasFile",			1, NULL, 'B'},
		{"outputSpectraFile",		1, NULL, 's'},
		{"inputgain",				1, NULL, 'G'},
		{"ordernumber",				1, NULL, 'o'},
		{"minorder",				1, NULL, 'O'},
		{"maxorder",				1, NULL, 'X'},          
		{"badpixelmask",			1, NULL, 'm'},
		{"aperture",				1, NULL, 'A'},	
		{"spectralElementHeight",	1, NULL, 'w'},
		{"wave",					1, NULL, 'W'},
		{"normalize",				1, NULL, 'N'},
		{"binsize",					1, NULL, 'b'},    
		{"usePolynomial",			1, NULL, 'l'},
		{"orderOfPolynomial",		1, NULL, 'r'},
		
		{"plot",					0, NULL, 'p'},
		{"verbose",					0, NULL, 'v'},
		{"debug",					0, NULL, 'd'},
		{"trace",					0, NULL, 't'},
		{"help",					0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:g:B:s:G:o:b:l:r:N:O:X:m:A:W:w:pvdth", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputImageName = optarg;
				break;   
			case 'g':
				inputGeometryName = optarg;
				break;
			case 'W':
				inputWavelengthName = optarg;
				break;
			case 'B':
				inputBiasName = optarg;
				break;
			case 's':
				outputSpectraFile = optarg;
				break;   
			case 'G':
				inputgain = optarg;
				break;   
			case 'o':
				ordernumber = atoi(optarg);
				break;   
			case 'O':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;  
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
				break; 				
			case 'm':		// badpixelmask
				badpixelmask = optarg;
				break;            					
			case 'A':		// aperture in pixels
				aperture = atof(optarg);
				break; 
			case 'w':		// element height in pixels
				spectralElementHeight = atof(optarg);
				break;
			case 'N':		// for normalization
				normalize = atoi(optarg)==1;
				break;
			case 'b':		// binsize
				binsize = atoi(optarg);
				break;
			case 'l':		
				usePolynomial = true;
				break;
			case 'r':
				orderOfPolynomial = atoi(optarg);
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
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			case '?':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
		}
	}	
	
	try {
		// we need an image...
		if (inputImageName.empty()) {
			throw operaException("operaExtractRawSum: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an output...
		if (outputSpectraFile.empty()) {
			throw operaException("operaExtractRawSum: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a geometry file...
		if (inputGeometryName.empty()) {
			throw operaException("operaExtractRawSum: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
		// we need a gain file...
		if (inputgain.empty()) {
			throw operaException("operaExtractRawSum: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
		
		if (verbose) {
			cout << "operaExtractRawSum: input image = " << inputImageName << endl; 
			cout << "operaExtractRawSum: input geometry file = " << inputGeometryName << endl; 			
			cout << "operaExtractRawSum: input wavelength file = " << inputWavelengthName << endl; 			
			cout << "operaExtractRawSum: input gain file = " << inputgain << endl;            
			cout << "operaExtractRawSum: output spectra file = " << outputSpectraFile << endl;
			cout << "operaExtractRawSum: badpixelmask = " << badpixelmask << endl; 				
			cout << "operaExtractRawSum: X-aperture = " << aperture << endl; 
			cout << "operaExtractRawSum: Y-height = " << spectralElementHeight << endl; 
			cout << "operaExtractRawSum: inputBiasName = " << inputBiasName << endl; 
		}
		
		//operaFITSImage *badpix = NULL;
		//if (!badpixelmask.empty()){ 
		//	badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);					
		//}		
		
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputGeometryName);
		if (!inputWavelengthName.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWavelengthName);
		}
		
		unsigned amp = 0;	// for gain/noise
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputgain);
		double gain = spectralOrders.getGainBiasNoise()->getGain(amp);
		double noise = spectralOrders.getGainBiasNoise()->getNoise(amp);
		
		if (verbose)
			cout << "operaExtractRawSum: gain="<< gain << " noise=" << noise << endl;

		unsigned minorder = spectralOrders.getMinorder();
		unsigned maxorder = spectralOrders.getMaxorder();
		
        if(!minorderprovided) {
            minorder = spectralOrders.getMinorder();
        }
        if(!maxorderprovided) {
            maxorder = spectralOrders.getMaxorder();            
        }
		
		if(ordernumber != NOTPROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}        
        
		operaFITSImage inputImage(inputImageName, tfloat, READONLY);	
		operaFITSImage *biasImage = NULL;	
		if (!inputBiasName.empty()) {
			biasImage = new operaFITSImage(inputBiasName, tfloat, READONLY);
			inputImage -= *biasImage;
		} else {
			inputImage -= 400.0;
		}

		/*
		 * Add the order data to the raw spectrum product
		 */
		float spectralElementHeight = 1.0;
		unsigned xsampling = 5; // 5 was found to be the best sampling for ESPaDOnS images.
		unsigned extraAperturePixels = 2; // those are extra-pixels to extend the aperture for the fitting but not for the normalization
		
		for (unsigned order=minorder; order<=maxorder; order++) {
			
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			if (spectralOrder->gethasGeometry()) {
				if ((!inputWavelengthName.empty() && spectralOrder->gethasWavelength()) || inputWavelengthName.empty()) {
					operaGeometry *geometry = spectralOrder->getGeometry();
					spectralOrder->setInstrumentProfileVector((unsigned)aperture + 2*extraAperturePixels + 1, xsampling, 1, 1, 1);
					geometry->setapertureWidth(aperture);	
					geometry->CalculateAndSetOrderLength();			
					spectralOrder->setSpectralElementsByHeight(spectralElementHeight);
					spectralOrder->extractRawSum(inputImage, noise, gain);
					if (!inputWavelengthName.empty()) {
						spectralOrder->CalculateWavelengthSolution();
					}
					if (normalize) {
						spectralOrder->applyNormalization(binsize,orderOfPolynomial,usePolynomial,true,false);
					}
				}
			}
		}
		if (inputWavelengthName.empty()) {
			operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputSpectraFile, RawSpectrum);
		} else {
			operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputSpectraFile, CalibratedRawSpectrum);
		}
		inputImage.operaFITSImageClose();
		
	}
	catch (operaException e) {
		cerr << "operaExtractRawSum: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaExtractRawSum: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	" --inputImage=<FITS_IMAGE>"
	" --inputGeometryFile=<GEOM_FILE>"
	" --outputSpectraFile=<FILE_NAME>"
	" --ordernumber=<INT_VALUE>"
	" --minorder=<INT_VALUE>"
	" --maxorder=<INT_VALUE>"       
	" --badpixelmask=<FITS_IMAGE>"
	" --aperture=<FLT_VALUE>"	
	" --spectralElementHeight=<FLT_VALUE>  \n\n"
	" Example: "+string(modulename)+" -g OLAPAa_sp2_Normal.geom -A 26 -w 1 -i masterfabperot_OLAPAa_sp2_Normal.fits -o 40 -v \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -i, --inputImage=<FITS_IMAGE>, Input FITS image to extract spectrum\n"
	"  -g, --inputGeometryNameetryFile=<GEOM_FILE>, Input geometry file\n"
	"  -s, --outputSpectraFile=<FILE_NAME>, Output file name\n"
	"  -o, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all)\n"
	" -N, --minorder=<INT_VALUE>, Define minimum order number\n"
	" -X, --maxorder=<INT_VALUE>, Define maximum order number\n"     
	"  -m, --badpixelmask=<FITS_IMAGE>, FITS image with badpixel mask\n"
	"  -A, --aperture=<FLT_VALUE>, Aperture for extraction in X-direction\n"
	"  -w, --spectralElementHeight=<FLT_VALUE>, Width of spectral element in Y-direction \n\n";	
}

