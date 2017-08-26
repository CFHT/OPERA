/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaExtractSpectralLines
 Version: 1.0
 Description: Extract spectral lines. 
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
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaExtractSpectralLines.h"

#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaEspadonsImage.h"			// for imtype_t
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaInstrumentProfile.h"		// for operaInstrumentProfile
#include "libraries/operaExtractionAperture.h"

#include "libraries/operaLibCommon.h"
#include "libraries/operaLib.h"						// systemf
#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaFit.h"	
#include "libraries/operaFFT.h"	

#define NOTPROVIDED -999

/*! \file operaExtractSpectralLines.cpp */

using namespace std;

/*! 
 * operaExtractSpectralLines
 * \author Doug Teeple
 * \brief extract spectral lines.
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
	int opt;
	
	string inputgeom;
	string inputprof; 
    
	string outputfile;
	string masterbias;    
	string mastercomparison;	    
	string badpixelmask;
    
	float spectralElementHeight = 1.0;  
	float gain = 1.12;
    float noise = 3.5;
	
    float referenceLineWidth = 2.5;
    
	float LocalMaxFilterWidth = 0;
	float DetectionThreshold = 0;         
    float MinPeakDepth = 0;
	
    int ordernumber = NOTPROVIDED;
    
	string plotfilename;	
	string linesdatafilename;
	string specdatafilename;
	string scriptfilename;
    string gainfilename;   

	bool interactive = false;
    
	int plot=0, verbose=0, debug=0,  trace=0;
	
	struct option longopts[] = {
		{"inputGeometryFile",			1, NULL, 'g'},
		{"inputInstrumentProfileFile",	1, NULL, 'i'},  
		{"outputSpectraFile",			1, NULL, 'o'},
		{"masterbias",					1, NULL, 'b'},        
		{"mastercomparison",			1, NULL, 'c'},
		{"gain",						1, NULL, 'G'},	
		{"badpixelmask",				1, NULL, 'm'},
		{"spectralElementHeight",		1, NULL, 'H'},
		{"referenceLineWidth",			1, NULL, 'W'},
		{"LocalMaxFilterWidth",			1, NULL, 'M'},
		{"DetectionThreshold",			1, NULL, 'T'},        
		{"MinPeakDepth",				1, NULL, 'D'},     
        {"ordernumber",					1, NULL, 'O'},  	
		{"plotfilename",				1, NULL, 'P'},
		{"linesdatafilename",			1, NULL, 'F'},
		{"specdatafilename",			1, NULL, 'C'},        
		{"scriptfilename",				1, NULL, 'S'},
		{"interactive",					0, NULL, 'I'},
		
		{"plot",						optional_argument, NULL, 'p'},       
		{"verbose",						optional_argument, NULL, 'v'},
		{"debug",						optional_argument, NULL, 'd'},
		{"trace",						optional_argument, NULL, 't'},
		{"help",						no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "g:i:o:b:c:m:H:W:M:T:D:O:P:F:C:S:G:Iv::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'g':
				inputgeom = optarg;
				break;    
			case 'i':
				inputprof = optarg;
				break;                 
			case 'o':
				outputfile = optarg;
				break; 
			case 'b':
				masterbias = optarg;
				break;                
			case 'c':
				mastercomparison = optarg;
				break;				
			case 'm':		// badpixelmask
				badpixelmask = optarg;
				break;            
			case 'H':		// aperture in pixels
				spectralElementHeight = atof(optarg);
				break; 	
			case 'W':		// line width for reference (in pixel units)
				referenceLineWidth = atof(optarg);
				break;  
			case 'M':
				LocalMaxFilterWidth = atof(optarg);
				break; 
			case 'T':
				DetectionThreshold = atof(optarg);
				break; 
			case 'D':
				MinPeakDepth = atof(optarg);
				break;                 
			case 'O':
				ordernumber = atoi(optarg);
				break;	                
			case 'G':		// gain
				gainfilename = optarg;
				break;
			case 'P':		// plot please
				plotfilename = optarg;
				plot = 1;
				break; 		                
			case 'F':
				linesdatafilename = optarg;
				break; 	
			case 'C':
				specdatafilename = optarg;
				break;                 
			case 'S':
				scriptfilename = optarg;
				break; 	                				
			case 'I':
				interactive = true;
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
        
		// we need a *.geom...
		if (inputgeom.empty()) {
			throw operaException("operaExtractSpectralLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a *.prof...        
		if (inputprof.empty()) {
			throw operaException("operaExtractSpectralLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a masterbias...
		if (masterbias.empty()) {
			throw operaException("operaExtractSpectralLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}        
		// we need a comparison...
		if (mastercomparison.empty()) {
			throw operaException("operaExtractSpectralLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}							
		
		if (verbose) {
			cout << "operaExtractSpectralLines: inputgeom = " << inputgeom << endl;
			cout << "operaExtractSpectralLines: inputprof = " << inputprof << endl;            
			cout << "operaExtractSpectralLines: outputfile = " << outputfile << endl;
			cout << "operaExtractSpectralLines: masterbias = " << masterbias << endl;                        
			cout << "operaExtractSpectralLines: mastercomparison = " << mastercomparison << endl; 				
			cout << "operaExtractSpectralLines: badpixelmask = " << badpixelmask << endl; 				
		}
        
		operaFITSImage bias(masterbias, tfloat, READONLY);			
		operaFITSImage comp(mastercomparison, tfloat, READONLY);

		operaFITSImage *badpix = NULL;
        
		if (!badpixelmask.empty()){              
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(comp.getnaxis1(),comp.getnaxis2(),tfloat);
            *badpix = 1.0;
        }
		
		operaSpectralOrderVector spectralOrders(inputgeom);
        spectralOrders.ReadIntoSpectralOrders(inputprof);
		spectralOrders.ReadIntoSpectralOrders(gainfilename);
		
		unsigned amp = 0;
		gain = spectralOrders.getGainBiasNoise()->getGain(amp);
		noise = spectralOrders.getGainBiasNoise()->getNoise(amp);
        
        if(LocalMaxFilterWidth==0)
            LocalMaxFilterWidth = 2*referenceLineWidth;
        if(DetectionThreshold==0)
            DetectionThreshold = 0.2;
        if(MinPeakDepth==0)
            MinPeakDepth = 1.5*noise;               
        
		unsigned minorder = spectralOrders.getMinorder();
		unsigned maxorder = spectralOrders.getMaxorder();
        
		if(ordernumber != NOTPROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}
        
        for (unsigned order=minorder; order<=maxorder; order++) {
			if (verbose)
				cout << "operaExtractSpectralLines: #Order number: "<< order << endl;
			// create pointer to current spectral order:
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);	
			if (spectralOrder->gethasGeometry()/* && spectralOrder->gethasInstrumentprofile()*/) { // FIXME instrumentprofile is not used???
																								   // create pointer to geometry for current spectral order:
				// create pointer to spectral lines for current spectral order:
				operaSpectralLines *spectralLines = NULL;
				
				// create a set of spectral elements based on given element height:
				spectralOrder->setSpectralElementsByHeight(spectralElementHeight);
				
				try {
                    spectralOrder->calculateXCorrBetweenIPandImage(comp, *badpix, NULL);  
					spectralOrder->setSpectralLines(comp, *badpix, bias, noise, gain, referenceLineWidth, DetectionThreshold, LocalMaxFilterWidth, MinPeakDepth);           
					spectralOrder->sethasSpectralLines(true);
					spectralLines = spectralOrder->getSpectralLines();
					if (verbose)
						cout << "operaExtractSpectralLines: " << spectralLines->getnLines() << " lines found in order " << order<< endl;
				} catch (operaException e) {
					if (verbose)
						cout << "operaExtractSpectralLines: No lines found, skipping order... " << order << endl;
				}
				
				// setSpectralLines creates  new class instance...
				spectralLines = spectralOrder->getSpectralLines();
				
				if (spectralOrder->gethasSpectralLines() && spectralLines->getnLines() > 0) {
					if (!plotfilename.empty()) {
						if (!linesdatafilename.empty() && !specdatafilename.empty() && !scriptfilename.empty()) {
							
							ofstream flinesdata;
							flinesdata.open(linesdatafilename.c_str());                    
							spectralLines->printLines(&flinesdata);
							flinesdata.close();
							
							ofstream fspecdata;
							fspecdata.open(specdatafilename.c_str());                      
							spectralLines->printReferenceSpectrum(&fspecdata);
							fspecdata.close(); 
							
							GenerateSpectralLinesPlot(scriptfilename.c_str(),plotfilename.c_str(),linesdatafilename.c_str(),specdatafilename.c_str(), interactive);
						}
						
					} else if (debug) {
						spectralLines->printLines(&cout);
					}
				}
			}
		}
		spectralOrders.WriteSpectralOrders(outputfile, Lines);

		bias.operaFITSImageClose();
		comp.operaFITSImageClose();
        if(badpix) {
            delete badpix;
        }
	}
	catch (operaException e) {
		cerr << "operaExtractSpectralLines: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaExtractSpectralLines: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}     

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	" --inputGeometryFile=<GEOM_FILE>"
	" --inputInstrumentProfileFile=<PROF_FILE>"
	" --outputSpectraFile=<LINES_FILE>"    	
	" --masterbias=<FITS_IMAGE>"    
	" --inputMasterComparison=<FITS_IMAGE>"	
	" --badpixelmask=<FITS_IMAGE>"	
	" --spectralElementHeight=<FLT_VALUE>"
	" --referenceLineWidth=<FLT_VALUE>"
	" --LocalMaxFilterWidth=<FLT_VALUE>"
	" --DetectionThreshold=<FLT_VALUE>"        
	" --MinPeakDepth=<FLT_VALUE>"    
	" --ordernumber=<INT_VALUE>"  	
	" --plotfilename=<EPS_FILE>"
	" --linesdatafilename=<DAT_FILE>"
	" --specdatafilename=<DAT_FILE>"    
	" --scriptfilename=<GNU_FILE>  \n\n"    
	
	" Example: "+string(modulename)+" -g OLAPAa_sp2_Normal.geom -i OLAPAa_sp2_Normal.prof -o OLAPAa_sp2_Normal.l -c mastercomparison_OLAPAa_sp2_Normal.fits -m badpixelmask.fits -W 3 -H 1 -O 40 -v \n\n"
	" -h, --help  display help message\n"
	" -v, --verbose,  Turn on message sending\n"
	" -d, --debug,  Turn on debug messages\n"
	" -t, --trace,  Turn on trace messages\n"
	" -g, --inputGeometryFile=<GEOM_FILE>, Input geometry file\n"
	" -i, --inputInstrumentProfileFile=<PROF_FILE>, Input Instrument Profile file\n"
	" -o, --outputSpectraFile=<LINES_FILE>, Output spectral lines file\n"    	
	" -b, --masterbias=<FITS_IMAGE>, Input Master Bias FITS image\n"	    
	" -c, --inputMasterComparison=<FITS_IMAGE>, Input Master Comparison (ThAr) FITS image\n"	
	" -m, --badpixelmask=<FITS_IMAGE>, FITS image with badpixel mask\n"
	" -H, --spectralElementHeight=<FLT_VALUE>, Width of spectral element in Y-direction \n"
	" -W, --referenceLineWidth=<FLT_VALUE>, Reference width of spectral line gaussian profile \n"    
	" -M, --LocalMaxFilterWidth=<FLT_VALUE>, Width within which to apply a local maximum filter\n"
	" -T, --DetectionThreshold=<FLT_VALUE>, Peak detection probability threshold \n"
	" -D, --MinPeakDepth=<FLT_VALUE>, Minimum peak depth cut-off for line detection \n"
	" -P, --plotfilename=<EPS_FILE>, Output plot eps file name \n"
	" -F, --linesdatafilename=<DAT_FILE>, Output spectral lines data file name \n"
	" -C, --specdatafilename=<DAT_FILE>, Output comparison spectrum data file name \n"    
	" -S, --scriptfilename=<GNU_FILE>, Output gnuplot script file name \n\n";
    
}		


void GenerateSpectralLinesPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *linesdataFileName,const char *specdataFileName, bool display)
{
    FILE *fgnu;
    remove(gnuScriptFileName); // delete any existing file with the same name
	
    fgnu = fopen(gnuScriptFileName,"w");
    
    fprintf(fgnu,"unset key\n");
    fprintf(fgnu,"set xlabel \"distance (pixels)\"\n");   
    fprintf(fgnu,"set ylabel \"flux (counts)\"\n");
    
    fprintf(fgnu,"set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
    fprintf(fgnu,"set output \"%s\"\n",outputPlotEPSFileName);
    
    fprintf(fgnu,"plot \"%s\" u 1:2 w l, \"%s\" u 8:($12/2+$4+$6*$8):10 with xerr\n", specdataFileName, linesdataFileName);
    
    if (display) {
		fprintf(fgnu,"set output\n");
		fprintf(fgnu,"set terminal x11\n");
        fprintf(fgnu,"replot\n");
        fclose(fgnu);        
        systemf("gnuplot -persist %s", gnuScriptFileName);
    } else { 
        fclose(fgnu);   
    }
}
