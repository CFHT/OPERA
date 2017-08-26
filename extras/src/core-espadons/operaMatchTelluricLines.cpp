/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMatchTelluricLines
 Version: 1.0
 Description: Wavelength Calibration 
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
#include <stdarg.h>
#include <getopt.h>
#include <math.h>									// for pow

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaWavelengthCalibration.h"

#include "libraries/operaException.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaSpectralLines.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/GainBiasNoise.h"
#include "libraries/Gaussian.h"
#include "libraries/Polynomial.h"						// for Polynomial

#include "libraries/operaLibCommon.h"					// for doubleValue_t
#include "libraries/operaLib.h"							// for itos
#include "libraries/operaMath.h"						// for LengthofPolynomial
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFFT.h"							// for operaXCorrelation

#define DEBUG false

/*! \brief wavelength calibration. */
/*! \file operaMatchTelluricLines.cpp */
/*! \ingroup core */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*
 * the reference atlas spectrum
 */
static unsigned nPointsInAtlasSpectrum = 0;
static double atlasWavelength[MAXFULLREFWAVELENGTHS];
static double atlasIntensity[MAXFULLREFWAVELENGTHS];   
static double atlasVariance[MAXFULLREFWAVELENGTHS];   
 
/*
 * the reference atlas spectral lines
 */
static unsigned thatlaslines = 0;
static double thAtlasWavelength[MAXREFWAVELENGTHS];
static double thAtlasIntensity[MAXREFWAVELENGTHS]; 

/*
 
 Below it follows in a few words a 1st-pass for the wavelength calibration
 algorithm.
 
 For each spectral order do the following steps:
 
 1. Read ThAr raw spectrum; intensity versus distance in pixel units:
 I(d) vs. d
 
 2. Measure total distance "D" (in pixel units) covered by the order. This
 is given by the line integral of the polynomial that describes the center
 of the order.
 
 3. Read wavelength range covered by the order: wl0,wlf
 
 4. Calculate first order solution:
 wl = f(d), where f(d) as first order is given by
 
 f(d) = wl0 + ((wlf - wl0)/D)*d
 
 assuming f(d=0) = wl0.
 
 5. Read ThAr atlas of spectral lines within the range covered by the order
 [wl0:wlf]. The atlas consists of line wavelength (l_wl), error (l_wlerr),
 and line relative intensity (l_i).
 
 7. Once we have a table of l_i, l_d, and l_derr, then we can calculate the
 maximum cross-correlation between this and the atlas data to identify the
 lines. The identification usually doesn't go one-by-one, so we will end up
 having to do some cleaning for either the undetected or over-detected
 lines.
 
 8. Now one can use the table (l_d+/-l_derr) versus (l_wl+/-l_wlerr) to
 find the wavelength solution by fitting a polynomial to these data.
 
 note that the polynomial should be an update to the first-order solution.
 
 The update is intended for two reasons:
 
 First because I have experienced before that the higher order terms are so
 small when compared to the first order that the fitting routine can get in
 trouble to find good solutions.
 
 Another reason is that we want to use our first solution to exclude
 outliers and then run the fitting again as many times as it gets to give
 the best solution. So, we will have to update our solution as we get our
 dataset cleaner or as we gather more information.
 
 */

/* 
 * operaMatchTelluricLines
 * \author Eder Martioli
 * \brief wavelength calibration.
 * \arg argc
 * \arg argv
 * \note --outputWave=...
 * \note --atlas_lines=...
 * \note --thcal=...
 * \note --geom=...
 * \note --wlcal_initialguess=...
 * \note --binsize=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

#define NOTPROVIDED -999

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth] \n\n"    
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -p, --plot,   Plot output \n"
	"  -o, --outputWave=\n"    
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n" 
	"  -I, --interactive=<BOOL>\n"    
	"\n";
}

int main(int argc, char *argv[])
{
	int opt;
	
	string outputWave; 
    
	string atlas_lines; 
    string atlas_spectrum;
	
    string uncalibrated_lines;
    string uncalibrated_spectrum;
    
    double uncalibrated_linewidth = 1.5; // default is 1.5 pixels
    
	string geometryfilename; 
	string wlcal_initialguess; 
    
    string plotfilename;	
	string datafilename;	
	string scriptfilename;	
    
	bool interactive = false;    

    unsigned minorder = 22;
    bool minorderprovided = false;
    unsigned maxorder = 61;    
    bool maxorderprovided = false; 
    
    int ordernumber = NOTPROVIDED;	
	
    /*
     * The parameters below we don't know yet whether they will be input
     */
    double DetectionThreshold = 0.05;    // threshold to regulate the sensitivity of line detection. Must be between 0 and 1.
    double LocalMaxFilterWidth = 3.0;     // parameter to set a window filter to guarantee a line is not detected twice. It's in units of line width
    double MinPeakDepth = 0.3;          // limit that also regulates the sensitity of line detection in units of noise.
    
    double ParRangeSizeInPerCent = 1.0; // define the range within which a coefficient will be changed to calcuate the x-correlation 
    unsigned NpointsPerPar = 100;        // define the number of times a coefficient is changed 
    unsigned maxpolyorder = 2;          // define the coefficients to be changed in order to search for maximum x-correlation 
    unsigned nIterforXcorr = 2;         // number of iterations to search for maximum x-correlation
    
    unsigned maxNIter = 30;                 // maximum number of iterations for shrinking the acceptable mismatch
    unsigned minNumberOfLines = 30;         // minimum number of lines to stop shrinking the acceptable mismatch difference between atlas and comparison
    unsigned maxorderofpolynomial = 5;      // maximum degree of polynomial for wavelength solution
    double dampingFactor = 0.85;            // Damping factor to shrink the size of the quantity acceptableMismatch on each iteration.  This factor may be set between 0 to 1. 
    double initialAcceptableMismatch = 6.0; // initial acceptable mismatch difference between atlas and comparison lines. Used for identification. In units of line width.
    
    unsigned nsig = 3;                      // in units of rms for clipping.
    
	struct option longopts[] = {
		{"atlas_lines",1, NULL, 'a'},
		{"atlas_spectrum",1, NULL, 's'},        
		{"uncalibrated_lines",1, NULL, 'l'},					// operaExtractSpactralLines does this for us
		{"uncalibrated_spectrum",1, NULL, 'u'},					// one may use the spectrum     
		{"uncalibrated_linewidth",1, NULL, 'w'},				// need this if using input spectrum  
		{"geom",		1, NULL, 'g'},	
		{"wlcal_initialguess",		1, NULL, 'r'},				// stores polynomials for each order to get a first guess as to the wavelengths
		{"outputWave",		1, NULL, 'o'},
        {"ordernumber",1, NULL, 'O'}, 
		{"minorder",1, NULL, 'N'},
		{"maxorder",1, NULL, 'X'},
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},
		{"interactive",0, NULL, 'I'}, 		
		{"plot",		optional_argument, NULL, 'p'},       
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "a:s:l:u:w:g:r:o:O:N:X:P:F:S:I:p::v::d::t::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'a':		// atlas_lines
				atlas_lines = optarg;
				break;
			case 's':		// atlas full spectrum
				atlas_spectrum = optarg;
				break;                
			case 'l':		// uncalibrated_lines
				uncalibrated_lines = optarg;
				break;    
			case 'u':		// uncalibrated_spectrum
				uncalibrated_spectrum = optarg;
				break;          
			case 'w':		// uncalibrated_linewidth: necessary only if input spectrum
				uncalibrated_linewidth = atof(optarg);
				break;                 
			case 'g':		// geometryfilename
				geometryfilename = optarg;
				break; 
			case 'r':		// wavelngth calibration reference polynomials
				wlcal_initialguess = optarg;
				break;  
			case 'o':		// output
				outputWave = optarg;
				break;
			case 'O':
				ordernumber = atoi(optarg);
				break;	
			case 'N':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;  
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
				break; 
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break; 		                
			case 'F':
				datafilename = optarg;
				break; 	
			case 'S':
				scriptfilename = optarg;
				break;
			case 'I':		// for interactive plots
				interactive = true;
				break; 
			case 'p':
				plot = 1;
				break;                 
			case 'v':
				verbose = 1;
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
		// we need a atlas_lines lines or spectrum...
		if (atlas_lines.empty() && atlas_spectrum.empty()) {
			throw operaException("operaMatchTelluricLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}

		// we need EITHER thorium uncalibrated lines or spectrum...
		if (uncalibrated_lines.empty() && uncalibrated_spectrum.empty()) {
			throw operaException("operaMatchTelluricLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
		// we need a geometryfilename...
		if (geometryfilename.empty()) {
			throw operaException("operaMatchTelluricLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a wlcal_initialguess, the initial guess at a polynomial...
		if (wlcal_initialguess.empty()) {
			throw operaException("operaMatchTelluricLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
		if (verbose) {
			cout << "operaMatchTelluricLines: atlas_lines = " << atlas_lines << endl;
			cout << "operaMatchTelluricLines: atlas_spectrum = " << atlas_spectrum << endl;            
			cout << "operaMatchTelluricLines: uncalibrated_lines = " << uncalibrated_lines << endl;
			cout << "operaMatchTelluricLines: uncalibrated_spectrum = " << uncalibrated_spectrum << endl;            
			cout << "operaMatchTelluricLines: geometryfilename = " << geometryfilename << endl;            
			cout << "operaMatchTelluricLines: wlcal_initialguess = " << wlcal_initialguess << endl; 
			cout << "operaMatchTelluricLines: outputWave = " << outputWave << endl;            
            if(ordernumber != NOTPROVIDED) {
                cout << "operaMatchTelluricLines: ordernumber = " << ordernumber << endl;            
            }
            if(plot) {
                cout << "operaMatchTelluricLines: plotfilename = " << plotfilename << endl;               
            }            
		}

		ofstream *fdata = NULL;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }
        
		operaSpectralOrderVector spectralOrders(geometryfilename);	// get the geometry

		spectralOrders.ReadIntoSpectralOrders(wlcal_initialguess); // read wavelength calibration reference first guess
        
        
		if (!uncalibrated_lines.empty()) {     
            spectralOrders.ReadIntoSpectralOrders(uncalibrated_lines); // This merges in the uncalibrated lines information
        }

        if (!uncalibrated_spectrum.empty()) { 
            spectralOrders.ReadIntoSpectralOrders(uncalibrated_spectrum); // This merges in the uncalibrated spectrum information
		}
        
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

		if (verbose) {
			cout << "operaMatchTelluricLines: minorder = " << minorder << " maxorder = " << maxorder << endl;            
		}        
        
		/*
		 * Read wavelength calibration reference first guess covered by the order: wl0,wlf
		 */        
		/*
		 * Read ThAr atlas spectrum
		 *		lambda vs. intensity, intensityVariance
		 */
		if (!atlas_spectrum.empty()) {      
            nPointsInAtlasSpectrum = readAtlasSpectrum(atlas_spectrum, atlasWavelength, atlasIntensity, atlasVariance);
        }
		/*
		 * Read ThAr atlas lines
		 *		lambda vs. intensity
		 */        
		if (!atlas_lines.empty()) {         
            thatlaslines = readThoriumArgonAtlas(atlas_lines, thAtlasWavelength, thAtlasIntensity);        
        }
        
        /*
         * vectors for uncalibrated spectral lines
         */
        unsigned rawlinesinorder=0;
        double rawlinecenter[MAXREFWAVELENGTHSPERORDER];
        double rawlinecenterError[MAXREFWAVELENGTHSPERORDER];        
        double rawlineflux[MAXREFWAVELENGTHSPERORDER];
        double rawlinesigma[MAXREFWAVELENGTHSPERORDER];
        float frawlinesigma[MAXREFWAVELENGTHSPERORDER];
        
        /*
         * vectors for atlas spectral lines
         */        
        unsigned atlaslinesinorder=0;
        double atlasLineswl[MAXREFWAVELENGTHSPERORDER]; 
        double atlasLineswlError[MAXREFWAVELENGTHSPERORDER];         
        double atlasLinesflux[MAXREFWAVELENGTHSPERORDER];        
             
        /*
         * cross-correlation
         */          
        double atlasXcorr[MAXPOINTSINSIMULATEDSPECTRUM];        
        double convolvedAtlas[MAXPOINTSINSIMULATEDSPECTRUM];        

        /*
         * Initialize linewidth with input uncalibrated_linewidth plus-minus 20% error. 
         */           
        double rawlinewidth;
        double rawlinewidth_err;
        
        for (unsigned order=minorder; order<=maxorder; order++) {
            unsigned bestnpar = maxorderofpolynomial; 
            
            rawlinewidth = uncalibrated_linewidth;
            rawlinewidth_err = uncalibrated_linewidth*0.2;
            
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
			if (spectralOrder->gethasGeometry() && spectralOrder->gethasWavelength()) {
				operaGeometry *geometry = spectralOrder->getGeometry();
				operaWavelength *wavelength = spectralOrder->getWavelength(); 
                
                operaSpectralElements *spectralElements = NULL;
                if(spectralOrder->gethasSpectralElements()) {
                    spectralElements = spectralOrder->getSpectralElements();
                    for(unsigned elemIndex=0;elemIndex<spectralElements->getnSpectralElements();elemIndex++) {
                        double elmewavelength = wavelength->evaluateWavelength(spectralElements->getdistd(elemIndex));
                        spectralElements->setwavelength(elmewavelength,elemIndex);
                        if(debug) // at this point one can check the initial wavelength calibration and the input comparison flux
                            cout << spectralElements->getwavelength(elemIndex) << " " << spectralElements->getFlux(elemIndex) << endl;
                    }
                    spectralElements->setHasWavelength(TRUE);
                }

                /*
                 * opera spectroscopic classes for the atlas
                 */          
                //operaFluxVector *atlasfluxvector = NULL;
                //operaSpectralElements *atlasSpectrum = NULL;
                //operaSpectralLines *atlasLines = NULL;                  
                
                double dmin = 0.0;
                double dmax = (double)geometry->calculateLength(geometry->getYmin(),geometry->getYmax());
                wavelength->setDmin(dmin);                
                wavelength->setDmax(dmax);
                
                if (verbose) {
					printf("operaMatchTelluricLines: Order %d: [geom] ymin = %.2f ymax = %.2f dmin = %.2f dmax = %.2f \n", spectralOrder->getorder(), geometry->getYmin(), geometry->getYmax(), dmin, dmax); 
				}   

                /*
                 * Below it calculates the initial and final wavelength based on the geometry calibration
                 * Note: this will be used to select the atlas range.  
                 */                     
                double wl0 = wavelength->getinitialWavelength();
                double wlf = wavelength->getfinalWavelength();
                wl0 -= (wavelength->getcentralWavelength())*(ParRangeSizeInPerCent/100);
                wlf += (wavelength->getcentralWavelength())*(ParRangeSizeInPerCent/100);                    
                
				if (verbose) {
					printf("operaMatchTelluricLines: Order %d: [wave] wavelength selected range: wl0 = %.2f wlc = %.2f wlf = %.2f\n",  spectralOrder->getorder(), wl0, wavelength->getcentralWavelength(),wlf);
				}
                                
                /*
                 * Read raw distance and flux from uncalibrated_lines file
                 */
                if (!uncalibrated_lines.empty()) {             
                    rawlinesinorder = getRawLines(spectralOrder, rawlinecenter, rawlinecenterError, rawlineflux, rawlinesigma); 
                    for (unsigned i=0; i<rawlinesinorder; i++)
                        frawlinesigma[i] = (float)rawlinesigma[i];	
                    rawlinewidth = (double)operaArrayMedian(rawlinesinorder,frawlinesigma);
                    rawlinewidth_err = (double)operaArrayMedSig(rawlinesinorder,frawlinesigma,(float)rawlinewidth);                
                }
                /*
                 * Below it reads raw distance and flux of lines detected from comparison input spectrum
                 * Note: the step below overrides the lines obtained from uncalibrated_lines above if both are provided. 
                 */            
                if (!uncalibrated_spectrum.empty() && spectralOrder->gethasSpectralElements()) {   
                    operaSpectralElements *compSpectrum = spectralOrder->getSpectralElements();
                    
                    operaSpectralLines compLines(compSpectrum,rawlinewidth,distance_disp);
                                     
                    operaFluxVector *compfluxvector = compSpectrum->getFluxVector();
                    /*
                    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
                        cout << compSpectrum->getdistd(i) << " " << compfluxvector->getflux(i) << " " << compSpectrum->getXCorrelation(i) << endl;
                    }
                    */
                    //compSpectrum->setHasXCorrelation(false);
                    if(!compSpectrum->getHasXCorrelation()){
						double *compSpectrumdistd = new double[MAXPOINTSINSIMULATEDSPECTRUM];         
						double *compSpectrumflux = new double[MAXPOINTSINSIMULATEDSPECTRUM]; 
                        double *compXcorr = new double[MAXPOINTSINSIMULATEDSPECTRUM]; 
                        
						for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
							compSpectrumdistd[i] = compSpectrum->getdistd(i);
							compSpectrumflux[i] = compfluxvector->getflux(i);
                           // cout << compSpectrumdistd[i] << " " << compSpectrumflux[i] << endl;
						}
                        calculateXCorrWithGaussian(compSpectrum->getnSpectralElements(), compSpectrumdistd,compSpectrumflux,compXcorr,rawlinewidth);                
						for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
							compSpectrum->setXCorrelation(compXcorr[i], i);
                        }
						compSpectrum->setHasXCorrelation(true); 
                        delete[] compSpectrumdistd;
                        delete[] compSpectrumflux;
                        delete[] compXcorr;                        
                    }
                    double CompLocalMaxFilterWidth = LocalMaxFilterWidth*rawlinewidth;

                    double meanVariance = 0;
                    unsigned nvarpoints = 0;
                    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
                        if(!isnan(compfluxvector->getvariance(i))) {
                            meanVariance += compfluxvector->getvariance(i);
                            nvarpoints++;
                        }
                    } 
                    double CompMinPeakDepth = MinPeakDepth*sqrt(meanVariance/(double)nvarpoints);

                    compLines.detectSpectralFeatures(DetectionThreshold,CompLocalMaxFilterWidth,CompMinPeakDepth);    
                                       
                    unsigned line = 0;

                    for(unsigned feature=0;feature<compLines.getNFeatures();feature++) {
                        operaSpectralFeature *currentFeature = compLines.getSpectralFeature(feature);
                        double *center = currentFeature->getGaussianFit()->getCenterVector();
                        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
                        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
                        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector(); 
                        for(unsigned l=0; l<currentFeature->getnLines(); l++) { 
                            rawlinecenter[line] = center[l];
                            rawlinecenterError[line] = centerError[l];
                            rawlineflux[line] = amplitude[l];
                            rawlinesigma[line] = sigma[l];
                            frawlinesigma[line] = (float)rawlinesigma[line];                            
                            if(debug)
                                cout << center[l] <<  " " << centerError[l] << " " << amplitude[l] << " " << sigma[l] << endl;
                            line++;
                        }
                    }

					rawlinesinorder = line;
                    
                    rawlinewidth = (double)operaArrayMedian(rawlinesinorder,frawlinesigma);
                    rawlinewidth_err = (double)operaArrayMedSig(rawlinesinorder,frawlinesigma,(float)rawlinewidth);   
//                    rawlinewidth = (double)operaArrayMean(rawlinesinorder,frawlinesigma);                    
//                    rawlinewidth_err = (double)operaArraySig(rawlinesinorder,frawlinesigma);
                    
                    compSpectrum->setHasWavelength(true);
                } 
                if (rawlinesinorder == 0) {
                    printf("operaMatchTelluricLines: Warning: Order %d: [Comparison] No lines detected from input comparison. Skipping calibration.\n", order); 
                    continue;
                } else {
                    /*
                     * Below it normalizes the raw lines
                     */                                      
                    if (verbose) {
                        printf("operaMatchTelluricLines: Order %d: [Comparison] %d lines in comparison between wl0 = %.2f and wlf = %.2f.\n", order, rawlinesinorder,wavelength->evaluateWavelength(rawlinecenter[0]),wavelength->evaluateWavelength(rawlinecenter[rawlinesinorder-1])); 
                    }
                    wavelength->createComparisonDataVectors(rawlinesinorder,rawlinecenter,rawlinecenterError,rawlineflux);
                }
       
                /*
                 * Read wavelength and flux from atlas_lines file
                 */                                         
                                
                if (!atlas_lines.empty()) {   
                    double *thwl = NULL, *thintensity = NULL;
                    atlaslinesinorder = getThoriumArgonAtlasRange(wl0, wlf, &thwl, &thintensity);  
                    for(unsigned l=0; l<atlaslinesinorder; l++) {
                        atlasLineswl[l] = *thwl++;
                        atlasLinesflux[l] = *thintensity++;
                        atlasLineswlError[l] = wavelength->convertPixelToWavelength(rawlinewidth);
                    }
                } 

                /*
                 * Below it reads wavelength and flux from the atlas spectrum
                 */    
                if (!atlas_spectrum.empty()) {    
                    double *thwl = NULL, *thintensity = NULL, *thvar = NULL;
                    unsigned npatlasspecinorder = getAtlasSpectrumRange(wl0, wlf, &thwl, &thintensity, &thvar);                
                    
                    /*
                     * Below it calculates the cross-correlation between the atlas spectrum and a gaussian function. 
                     */                                   
                    calculateXCorrWithGaussian(npatlasspecinorder,thwl,thintensity,atlasXcorr,wavelength->convertPixelToWavelength(rawlinewidth));                
                    
                    /*
                     * Below it degrades the resolution of the atlas to the resolution of raw lines.
                     * The degradation is done by convolving the spectrum with a gaussian. 
                     */                      
                    
                    convolveSpectrumWithGaussian(npatlasspecinorder,thwl,thintensity,convolvedAtlas,wavelength->convertPixelToWavelength(rawlinewidth));                
 
                    /*
                     * Below it reads the atlas spectrum into an operaSpectralElements class
                     */                     
                    operaSpectralElements atlasSpectrum(npatlasspecinorder);
					for (unsigned i=0; i<npatlasspecinorder; i++) {
						atlasSpectrum.setXCorrelation(atlasXcorr[i], i);
					}
                    atlasSpectrum.setHasXCorrelation(true); 
                    //atlasSpectrum.setwavelengthVector(thwl);
					for (unsigned i=0; i<npatlasspecinorder; i++) {
						atlasSpectrum.setwavelength(thwl[i], i);
					}
                    atlasSpectrum.setHasWavelength(true);
                    operaFluxVector atlasfluxvector(convolvedAtlas,thvar,npatlasspecinorder);
                    atlasSpectrum.setFluxVector(&atlasfluxvector);	// copies, does not set local stack address...
                    atlasSpectrum.setHasRawSpectrum(true);
 
                    /*
                     * Below it creates an operaSpectralLines class for the atlas lines
                     */                     
                    operaSpectralLines atlasLines(&atlasSpectrum, wavelength->convertPixelToWavelength(rawlinewidth), wavelength_disp);
                    if(debug) {
                        for (unsigned i=0; i<npatlasspecinorder; i++) {
                            cout << atlasSpectrum.getwavelength(i)  << " " << atlasSpectrum.getFlux(i) << " " << atlasSpectrum.getFluxVariance(i) << " " << atlasSpectrum.getXCorrelation(i) << endl;
                        }
                    }

                    /*
                     * Below it set detection thresholds and run the algorithm to detect spectral lines in the atlas
                     */  

                    double AtlasLocalMaxFilterWidth = LocalMaxFilterWidth*wavelength->convertPixelToWavelength(uncalibrated_linewidth);
                    double AtlasMinPeakDepth = MinPeakDepth*sqrt(operaArrayMean_d(npatlasspecinorder,thvar));
                    atlasLines.detectSpectralFeatures(DetectionThreshold,AtlasLocalMaxFilterWidth,AtlasMinPeakDepth);
                    
                    /*
                     * Below it reads the atlas lines information
                     */                      
                    atlaslinesinorder = atlasLines.getnLines();
                    
                    unsigned line = 0;
                    
                    for(unsigned feature=0;feature<atlasLines.getNFeatures();feature++) {
                        operaSpectralFeature *currentFeature = atlasLines.getSpectralFeature(feature);
                        double *center = currentFeature->getGaussianFit()->getCenterVector();
                        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();                        
                        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
                        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector(); 
                        for(unsigned l=0; l<currentFeature->getnLines(); l++) { 
                            atlasLineswl[line] = center[l];
                            atlasLineswlError[line] = centerError[l];                            
                            atlasLinesflux[line] = amplitude[l];
                            if(debug)
                                cout << center[l] <<  " " << amplitude[l] << " " << sigma[l] << endl;
                            line++;
                        }
                    }
                }   
                                                  
                if (atlaslinesinorder == 0) {
                    printf("operaMatchTelluricLines: Warning:  Order %d: [Atlas] No lines detected from input atlas. Skipping calibration.\n", order);
                    continue;
                } else {                             
                    if (verbose) {
                        printf("operaMatchTelluricLines: Order %d: [Atlas] %d lines detected in input atlas between wl0 = %.2f and wlf = %.2f .\n", order, atlaslinesinorder, wl0, wlf);
                    }
                    wavelength->createAtlasDataVectors(atlaslinesinorder,atlasLineswl, atlasLineswlError,atlasLinesflux);
                }

                /*
                 * At this point all possible lines either in comparison or in the Atlas have been read.
                 * So, it starts identification of lines and refining wavelength solution.
                 */                    
                wavelength->createDataVectors(MAXREFWAVELENGTHS);
                doubleValue_t ResolutionElementInPixels = {2*rawlinewidth, rawlinewidth_err};
                /*
                 * The loop below adjust the wavelength solution to the maximum correlation between the two sets of lines. 
                 * The adjustment is done by varying each coefficient of the polynomial solution and finding the maximum
                 * x-correlation.
                 * maxpolyorder is the highest coefficient that is changed in the search
                 * NpointsPerPar is the number of trial solutions for a given coefficient
                 * ParRangeSizeInPerCent is a quantity that defined the range for the search. The range is set 
                 * relative to the value of the coefficient.
                 */  

                for(unsigned iter = 0; iter < nIterforXcorr; iter++) {
                    wavelength->calculateSpectralResolution(ResolutionElementInPixels);                
                    wavelength->refineWavelengthSolutionByXCorrelation(NpointsPerPar, ParRangeSizeInPerCent, maxpolyorder);
                }
                double acceptableMismatch = initialAcceptableMismatch; // in units of sigma                 
                
                wavelength->calculateSpectralResolution(ResolutionElementInPixels);
                wavelength->matchAtlaswithComparisonLines(acceptableMismatch);

                if(wavelength->getnDataPoints() <= bestnpar) {
                    bestnpar = wavelength->getnDataPoints()-1;
                }                
                wavelength->CalculateWavelengthSolution(bestnpar,false);

                Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
                double *bestpar = (double *)wavelengthPolynomial->getVector();                
                bestnpar = wavelengthPolynomial->getOrderOfPolynomial();
                
                double minchisqr = wavelengthPolynomial->getChisqr();
                unsigned nochangeinChisqr = 0;
                
                if (verbose) {                  
                    double ComparisonMatchPercentage =  wavelength->getPerCentageOfComparisonMatch();    
                    double AtlasMatchPercentage =  wavelength->getPerCentageOfAtlasMatch();                     
                    printf("\noperaMatchTelluricLines: Order %d: Initial Solution:\n", order);
                    printf("operaMatchTelluricLines: Order %d: ", order);
                    for(unsigned k=0;k<bestnpar;k++) {
                        cout << "par[" << k << "]=" << bestpar[k] << " "; 
                    }                    
                    cout <<  " chisqr=" << minchisqr << endl; 
                    printf("operaMatchTelluricLines: Order %d: %u lines matched between wl0 = %.2f  wlf = %.2f.\n", order,wavelength->getnDataPoints(),wavelength->getinitialWavelength(),wavelength->getfinalWavelength());
                    printf("operaMatchTelluricLines: Order %d: [Atlas]    matched %.2f %% of detected lines.\n", order,AtlasMatchPercentage);
                    printf("operaMatchTelluricLines: Order %d: [Comparison] matched %.2f %% of detected lines.\n", order,ComparisonMatchPercentage);                        
                }                
                
                /*
                 * Below it starts the iterations to improve the polynomial that gives the pixel-to-wavelength solution
                 */  
                for(unsigned iter=0; iter < maxNIter; iter++) {	
    
                    /*** NOTE ***/
                    // the acceptable mismatch can start considerably big and then shrink down as the calibration gets better.
                    // But it will only shrink to a minimum value or minimum number of lines. 
                    if(wavelength->getnDataPoints() > minNumberOfLines) {                                          
                        acceptableMismatch *= dampingFactor; 
                    }                    
                    wavelength->calculateSpectralResolution(ResolutionElementInPixels);                    
                    wavelength->matchAtlaswithComparisonLines(acceptableMismatch);
                    wavelength->filterDataPointsBySigmaClip(nsig);                    
                    
                    wavelength->RefineWavelengthSolution(bestnpar,false);
                             
                    if(wavelength->getnDataPoints() <= bestnpar) {
                            bestnpar = wavelength->getnDataPoints()-1;
                    }
                    
                    wavelength->CalculateWavelengthSolution(bestnpar,false);

                    wavelengthPolynomial = wavelength->getWavelengthPolynomial();
                    bestpar = (double *)wavelengthPolynomial->getVector();
                    bestnpar = wavelengthPolynomial->getOrderOfPolynomial();
                    
                    if(wavelengthPolynomial->getChisqr() < minchisqr) {
                        minchisqr = wavelengthPolynomial->getChisqr();
                        nochangeinChisqr = 0;
                    } else if(wavelengthPolynomial->getChisqr() == minchisqr) {
                        nochangeinChisqr++;
                    }             
                  /*  
                    if(verbose) {
                        double ComparisonMatchPercentage =  wavelength->getPerCentageOfComparisonMatch();    
                        double AtlasMatchPercentage =  wavelength->getPerCentageOfAtlasMatch();                         
                        double *wl = wavelength->getWavelengthDataVector();
                        printf("operaMatchTelluricLines: Order %d: [Iteration] iter= %u. \n", order, iter);                        
                        for(unsigned k=0;k<bestnpar;k++) {
                            cout << " par[" << k << "]=" << bestpar[k] << " "; 
                        }                    
                        cout <<  ". Chisqr=" << minchisqr << endl;    
                        printf("operaMatchTelluricLines: Order %d: %u lines matched between wl0 = %.2f  wlf = %.2f.\n", order,wavelength->getnDataPoints(),wl[0],wl[wavelength->getnDataPoints()-1]);
                        printf("operaMatchTelluricLines: Order %d: [Atlas]    matched %.2f %% of detected lines.\n", order,AtlasMatchPercentage);
                        printf("operaMatchTelluricLines: Order %d: [Comparison] matched %.2f %% of detected lines.\n", order,ComparisonMatchPercentage);                         
                    }
                    */                
                    if(nochangeinChisqr > 3) {
                        break;
                    }  
                } // for(unsigned iter=0; iter < nIter; iter++)    
                
                if(debug) {
                    for(unsigned l=0 ; l<wavelength->getnDataPoints(); l++) {
                        cout << order << " " << wavelength->getDistance(l) << " " << wavelength->getWavelength(l) << " " << wavelength->evaluateWavelength(wavelength->getDistance(l)) << " " << wavelength->getWavelength(l) - wavelength->evaluateWavelength(wavelength->getDistance(l)) << " " << wavelength->getWavelengthError(l) << endl;
                    }
                }
                if (verbose) {
                    double ComparisonMatchPercentage =  wavelength->getPerCentageOfComparisonMatch();    
                    double AtlasMatchPercentage =  wavelength->getPerCentageOfAtlasMatch();                     
                    wavelengthPolynomial = wavelength->getWavelengthPolynomial();
                    bestpar = (double *)wavelengthPolynomial->getVector();
                    bestnpar = wavelengthPolynomial->getOrderOfPolynomial();                    
                    cout << "operaMatchTelluricLines: Order z" << order << ": Final Solution:" << endl;  
                    cout << "operaMatchTelluricLines: Order " << order << ": ";
                    for(unsigned k=0;k<bestnpar;k++) {
                        cout << "par[" << k << "]=" << bestpar[k] << " "; 
                    }                    
                    cout <<  " chisqr=" << minchisqr << endl;                     
                    printf("operaMatchTelluricLines: Order %d: %u lines matched between wl0 = %.2f  wlf = %.2f.\n", order,wavelength->getnDataPoints(),wavelength->getinitialWavelength(),wavelength->getfinalWavelength());
                    printf("operaMatchTelluricLines: Order %d: [Atlas]    matched %.2f %% of detected lines.\n", order,AtlasMatchPercentage);
                    printf("operaMatchTelluricLines: Order %d: [Comparison] matched %.2f %% of detected lines.\n", order,ComparisonMatchPercentage);                        
                }
                
                /*
                 ***** at this point there should be a good wavelength solution *****
                 */ 
                
                /*
                 * Below it calculates the radial velocity precision
                 */                 
                wavelength->calculateRadialVelocityPrecision();
                
                /*
                 * Below it calculates the spectral resolution
                 */                  
                if (verbose){
                    printf("operaMatchTelluricLines: -----------------------------------------------------------------\n");
                    printf("operaMatchTelluricLines: Order %d: Radial velocity precision = %.2f m/s. Using %d spectral lines.\n", order, wavelength->getRadialVelocityPrecision(),wavelength->getnDataPoints());
                }
                if (verbose)
                    printf("operaMatchTelluricLines: Order %d: Wavelength RMS precision = %.10f nm  Median Precision = %.10f nm.\n", order,wavelength->calculateWavelengthRMSPrecision(),wavelength->calculateWavelengthMedianPrecision());
                
                if (verbose)
                    printf("operaMatchTelluricLines: Order %d: [Comparison Lines] median sigma = %.2f +/- %.2f.\n", order, rawlinewidth, rawlinewidth_err);
                
                wavelength->calculateSpectralResolution(ResolutionElementInPixels);
                
                if (debug)
                    printf("%d %.2f %.2f %.2f %.2f %.2f\n", order, wavelength->getcentralWavelength(),ResolutionElementInPixels.value, ResolutionElementInPixels.error, wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
                
                if (verbose) {
                    printf("operaMatchTelluricLines: Order %d: Spectral Resolution = %.2f +/- %.2f.\n", order, wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
                    printf("operaMatchTelluricLines: -----------------------------------------------------------------\n\n");                
                }
               
                if (debug) {                
                    cout << order << " " 
                     << wavelength->getinitialWavelength() <<  " "
                     << wavelength->getcentralWavelength() << " " 
                     << wavelength->getfinalWavelength() << " "
                     << wavelength->calculateWavelengthRMSPrecision() << " "
                     << wavelength->calculateWavelengthMedianPrecision() << " "
                     << wavelength->getRadialVelocityPrecision() << " "
                     << ResolutionElementInPixels.value << " "
                     << ResolutionElementInPixels.error << " "
                     << wavelength->getSpectralResolution().value << " "
                     << wavelength->getSpectralResolution().error << endl;
                }
                
            } else if (!spectralOrder->gethasWavelength()) {
                if (verbose) {
                    printf("operaMatchTelluricLines: Order %d: has no associated wavelength reference calibration data. Wavelength calibration not possible.\n",  spectralOrder->getorder());
                }
            } else if (!spectralOrder->gethasGeometry()) {
                if (verbose) {
                    printf("operaMatchTelluricLines: Order %d: has no associated geometry data. Wavelength calibration not possible.\n",  spectralOrder->getorder());
                }
                spectralOrder->sethasWavelength(false);
            } else {
                if (verbose) {
                    printf("operaMatchTelluricLines: Order %d: has neither geometry nor wavelength reference calibration data. Wavelength calibration not possible.\n",  spectralOrder->getorder());
                }
            }            
		} // for (unsigned order=minorder; order<=maxorder; order++)
                
        if (fdata != NULL) {
            fdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                //GenerateWavelengthSolutionPlot(scriptfilename.c_str(),plotfilename.c_str(),datafilename.c_str(), interactive);
            }
        }       

        spectralOrders.WriteSpectralOrders(outputWave, Wave);  
        
        //delete[] convolvedAtlas;
        exit(EXIT_SUCCESS);
	}
	catch (operaException e) {
		cerr << "operaMatchTelluricLines: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMatchTelluricLines: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
}


/*
 * get the raw lines, between wl0 and wlf
 */
unsigned getRawLines(operaSpectralOrder *spectralOrder, double *rawlinecenter, double *rawlinecenterError, double *rawlineflux, double *rawlinesigma) {
	unsigned rawlinesinorder = 0;
	
	if (spectralOrder->gethasSpectralLines()) {
		operaSpectralLines *spectralLines = spectralOrder->getSpectralLines();
		unsigned nFeatures = spectralLines->getNFeatures();
		for (unsigned featurenumber=0;featurenumber<nFeatures;featurenumber++) {
			operaSpectralFeature *spectralFeature = spectralLines->getSpectralFeature(featurenumber);
			double *center = spectralFeature->getGaussianFit()->getCenterVector();
			double *sigma = spectralFeature->getGaussianFit()->getSigmaVector();
			double *amplitude = spectralFeature->getGaussianFit()->getAmplitudeVector();        
			double *centerError = spectralFeature->getGaussianFit()->getCenterErrorVector();
			//double *sigmaError = spectralFeature->getGaussianFit()->getSigmaErrorVector();
			//double *amplitudeError = spectralFeature->getGaussianFit()->getAmplitudeErrorVector(); 
			
			for (unsigned line=0; line<spectralFeature->getnLines(); line++) {
				rawlinecenter[rawlinesinorder] = *center++;
                rawlinecenterError[rawlinesinorder] = *centerError++;
				rawlineflux[rawlinesinorder] = *amplitude++;
                rawlinesigma[rawlinesinorder] = *sigma++;
				rawlinesinorder++;
			}
		}
	}
	return rawlinesinorder;
}
/*
 * get a subset of the thatlas lines for this order only, between wl0 and wlf
*/ 
unsigned getThoriumArgonAtlasRange(double wl0, double wlf, double **thwl, double **thi) {
	unsigned firstline = 0;
	unsigned line = 0;

	for (line=0; line<thatlaslines; line++) {
		if (thAtlasWavelength[line] >= wl0) {
			if (firstline == 0) {
				*thi = &thAtlasIntensity[line];
				*thwl = &thAtlasWavelength[line];
				firstline = line;
			}
			if (thAtlasWavelength[line] > wlf)
				break;
		}
	}
	if (line)
		line--;
	if (line > firstline) {
		return (line-firstline);
	} else {
		return 0;
	}
}

/*
 * Read the entire thorium argon atlas
 * and normalize the results
*/

unsigned readThoriumArgonAtlas(string atlas_lines, double *thAtlasWavelength, double *thAtlasIntensity) {
	ifstream astream;
	string dataline;
	double tmp = -1.0; 
	double tmpi = -1.0; 
	unsigned line = 0;
	
	astream.open(atlas_lines.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					char buff[8];
					sscanf(dataline.c_str(), "%*f %lf %lf %s", &tmp, &tmpi, buff);
                    
					if (strcmp(buff, "Th") && strcmp(buff, "Ar")) {
						if (DEBUG)
							printf("Skipping non-thorium-argon atlas entry %s.\n", buff); 
					} else {
						tmp *= 0.1;
						tmpi = pow(10, tmpi);
                        
                        thAtlasWavelength[line] = tmp;
                        thAtlasIntensity[line] = tmpi;
                        line++;  
                    }
 				}	// skip comments
			}	// if strlen
		} // while (astream.good())
		line--;
		if (line > 0) {
			if (verbose) {
				printf("          [Atlas] %d lines found wl0=%.2f wlc=%.2f wlf=%.2f\n", line, thAtlasWavelength[0], thAtlasWavelength[line/2], thAtlasWavelength[line-1]);
			}
		} else {
			printf("          [Atlas] no lines found in atlas.\n");
		}
		astream.close();
	}	// if (astream.open()
	return line;
}

/*
 * Read the the full atlas spectrum
 */
unsigned readAtlasSpectrum(string atlas_spectrum, double *atlasWavelength, double *atlasIntensity, double *atlasVariance) {
	ifstream astream;
	string dataline;
    
	double tmpwl = -1.0; 
	double tmpi = -1.0; 
	double tmp1 = -1.0; 
	double tmp2 = -1.0; 
	double tmpvar = -1.0; 
	unsigned np = 0;
	
	astream.open(atlas_spectrum.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf %lf %lf %lf", &tmpwl, &tmpi, &tmp1, &tmp2, &tmpvar);
                    
                    atlasWavelength[np] = 0.1*tmpwl;
                    atlasIntensity[np] = tmpi;
                    atlasVariance[np] = tmpvar;
                    np++;  
                }	// skip comments
            }
		} // while (astream.good())

		if (np > 0) {
			if (verbose) {
				printf("          [Atlas] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", np, atlasWavelength[0], atlasWavelength[np/2], atlasWavelength[np-1]);
			}
		} else {
			printf("          [Atlas] no points found in atlas.\n");
		}
		astream.close();
	}	// if (astream.open()
	return np;
}

/*
 * get a subset of the atlas spectrum for this order only, between wl0 and wlf
 */
unsigned getAtlasSpectrumRange(double wl0, double wlf, double **thwl, double **thi, double **thvar) {
	unsigned firstline = 0;
	unsigned np = 0;
    
	for (np=0; np<nPointsInAtlasSpectrum; np++) {
		if (atlasWavelength[np] >= wl0) {
			if (firstline == 0) {
				*thi = &atlasIntensity[np];
				*thwl = &atlasWavelength[np];
                *thvar = &atlasVariance[np];
				firstline = np;
			}
			if (atlasWavelength[np] > wlf)
				break;
		}
	}
	if (np)
		np--;
	if (np > firstline) {
		return (np-firstline);
	} else {
		return 0;
	}
}
/*
 * Normalize spectrum - version 1
 */
void normalizeSpectrum(unsigned nLines, double *lineflux) {
	double maxflux = 0;
    
    for (unsigned i=0; i<nLines; i++) {
		if(lineflux[i] > maxflux) {
            maxflux = lineflux[i];
        }
	}    
	for (unsigned i=0; i<nLines; i++) {
		lineflux[i] /= maxflux/100;	
	}
}
/*
 * Normalize spectrum - version 2, include variance
 */
void normalizeSpectrum(unsigned nLines, double *lineflux, double *linevariance) {
	double maxflux = 0;
    
    for (unsigned i=0; i<nLines; i++) {
		if(lineflux[i] > maxflux) {
            maxflux = lineflux[i];
        }
	}    
	for (unsigned i=0; i<nLines; i++) {
		lineflux[i] /= maxflux/100;
        linevariance[i] /= (maxflux/100)*(maxflux/100);
	}
}

/*
 * Convolve spectrum with a Gaussian function 
 */                
void convolveSpectrumWithGaussian(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double sigma) {
    // first figure out how many points of the spectrum is covered by 4*sigma
    double wlstep = fabs(wavelength[np-1] - wavelength[0])/(double)np; 
    unsigned window = (unsigned)ceil(4*sigma/wlstep);
   
    memset(convolvedSpectrum, 0, sizeof(double)*np);    
    
    for(unsigned i=window/2;i<(np-window/2);i++) {
        double weighSum = 0;
        for(unsigned j=0;j<window;j++) {
            convolvedSpectrum[i] += flux[i-window/2+j]*exp(-((wavelength[i-window/2+j] - wavelength[i])*(wavelength[i-window/2+j] - wavelength[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma); 
            weighSum += exp(-((wavelength[i-window/2+j] - wavelength[i])*(wavelength[i-window/2+j] - wavelength[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma); 
        }
        convolvedSpectrum[i] /= weighSum;
    }
}

/*
 * Calculate cross-correlation between spectrum and a normalized Gaussian function 
 */  
void calculateXCorrWithGaussian(unsigned np, double *wavelength, double *flux, double *outputXcorr, double sigma) {
    
    // first figure out how many points of the spectrum is covered by 4*sigma
    double wlstep = fabs(wavelength[np-1] - wavelength[0])/(double)np; 
    unsigned window = (unsigned)ceil(4*sigma/wlstep);
    
    //set up window function
    double *windowFunc = (double *) malloc (window * sizeof(double));
    double *mainFunc = (double *) malloc (window * sizeof(double));
    memset(windowFunc, 0, sizeof(double)*window);   
    
    //calculate window function    
    for(unsigned j=0;j<window;j++) {
         windowFunc[j] = exp(-double((j - window/2)*(j - window/2))/(2.0*double(window*window)/(4.0*4.0)))/(sqrt(2.0*M_PI)*double(window/4)); 
    }    
    
    //calculate xcorrelation function    
    for(unsigned i=window/2;i<(np-window/2);i++) {
        for(unsigned j=0;j<window;j++) {
            mainFunc[j] = flux[i-window/2+j];
        }
        outputXcorr[i] = operaCrossCorrelation(window,windowFunc,mainFunc);
    }
}




