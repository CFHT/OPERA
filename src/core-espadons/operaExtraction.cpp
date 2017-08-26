/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaExtraction
 Version: 1.0
 Description: Extract spectrum using various alogorithms.
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

#include <fstream>
#include <pthread.h>
#include "libraries/operaIOFormats.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

/*! \file operaExtraction.cpp */

using namespace std;

/*! 
 * operaExtraction
 * \author Eder Martioli
 * \brief Module to extract spectra using various alogorithms.
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

void GenerateExtraction3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, unsigned numberOfBeams, bool display);

operaArgumentHandler args;

string outputSpectraFile;
string inputImage;
string inputgain;
string inputgeom;
string inputprof;
string inputaper;
string masterbias;
string masterflat;
string normalizedflatfile;
string badpixelmask;
int ordernumber = NOTPROVIDED;
int minorder = NOTPROVIDED;
int maxorder = NOTPROVIDED;
double biasConstantToAdd = 0;
double effectiveApertureFraction = 0.99;
unsigned backgroundBinsize = 1;
double minsigmaclip = 25;
double sigmacliprange;
unsigned iterations = 3;
bool onTargetProfile = false;
bool usePolynomialFit = false;
bool removeBackground = false;
bool starplusskymode = false;
bool starplusskyInvertSkyFiber = false;
double skyOverStarFiberAreaRatio = 1.0;
bool noCrossCorrelation = false;
bool rejectBadpixInRawExtraction = false;
unsigned spectralOrderType_val = RawBeamSpectrum;
string spectrumtypename;
unsigned maxthreads = 1;
string plotfilename;
string datafilename;
string scriptfilename;
bool interactive = false;

operaSpectralOrderVector spectralOrders;
operaSpectralOrder_t spectralOrderType;
GainBiasNoise *gainBiasNoise = NULL;

operaFITSImage *bias = NULL;
operaFITSImage *normalizedflat = NULL;
operaFITSImage *badpix = NULL;
operaFITSImage *flat = NULL;
operaFITSImage *object = NULL;

/*
 * Thread Support to process all orders in parallel
 */

typedef struct thread_args {
	int order;
} thread_args_t;

pthread_t *threads = NULL;
thread_args_t *thread_args = NULL;

void *processOrder(void *argument) {
	thread_args_t *thread_args_s = (thread_args_t *)argument;
	int order = thread_args_s->order;
    
    operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
    
    if (args.verbose) {
		if (!spectralOrder->gethasGeometry()) cout << "operaExtraction: Skipping order number: "<< order << " no geometry." << endl;
		if (!spectralOrder->gethasInstrumentProfile()) cout << "operaExtraction: Skipping order number: "<< order << " no instrument profile." << endl;
		if (!spectralOrder->gethasExtractionApertures()) cout << "operaExtraction: Skipping order number: "<< order << " no extraction aperture." << endl;
	}

    if (spectralOrder->gethasGeometry() &&
        spectralOrder->gethasInstrumentProfile() &&
        spectralOrder->gethasExtractionApertures()) {
        
        if (args.verbose) cout << "operaExtraction: Processing order number: "<< order << endl;
        
        switch (spectralOrderType) {
            case RawBeamSpectrum:
                if(rejectBadpixInRawExtraction) {
                    spectralOrder->extractRawSpectrumRejectingBadpixels(*object, *flat, *normalizedflat, *bias, *badpix, *gainBiasNoise, effectiveApertureFraction, backgroundBinsize, minsigmaclip, iterations, onTargetProfile, usePolynomialFit, removeBackground, args.verbose, !noCrossCorrelation, NULL);
                } else {
                    spectralOrder->extractRawSpectrum(*object, *normalizedflat, *bias, *badpix, *gainBiasNoise, effectiveApertureFraction);
                    if(!noCrossCorrelation) spectralOrder->calculateXCorrBetweenIPandImage(*object,*badpix,NULL);
                }
                break;
            case StandardBeamSpectrum:
                spectralOrder->extractStandardSpectrum(*object, *normalizedflat, *bias, *badpix, *gainBiasNoise, effectiveApertureFraction, backgroundBinsize);
                if(!noCrossCorrelation) spectralOrder->calculateXCorrBetweenIPandImage(*object,*badpix,NULL);
                break;
            case OptimalBeamSpectrum:
				spectralOrder->extractOptimalSpectrum(*object, *flat, *normalizedflat, *bias, *badpix, *gainBiasNoise, effectiveApertureFraction, backgroundBinsize, minsigmaclip, sigmacliprange, iterations, onTargetProfile, usePolynomialFit, removeBackground, args.verbose, !noCrossCorrelation, NULL);
                break;
            default:
                break;
        }
        spectralOrder->normalizeFluxToAperture();
        if(starplusskymode) {
            spectralOrder->calculateStarAndSkyElements(starplusskyInvertSkyFiber, skyOverStarFiberAreaRatio);
        }
    } else {
        spectralOrder->sethasSpectralElements(false);
    }
	return NULL;
}

static void processSingleOrder(int order) {
    thread_args[0].order = order;
    processOrder((void *) &thread_args[0]);
}

static bool spawnthreads(int order, int maxorder, int count) {
	int j = 0;
    for (int i=order; i<=maxorder; i++) {
		thread_args[i].order = i;
		if (pthread_create(&threads[i], NULL, processOrder, (void *) &thread_args[i]) != 0)
			return false;
        if (++j >= count)
            break;
	}
    return true;
}

static bool waitthreads(int order, int maxorder, int count) {
	int j = 0;
	for (int i=order; i<=maxorder; i++) {
		if (pthread_join(threads[i], NULL) != 0)
			return false;
        if (++j >= count)
            break;
	}
    return true;
}

static bool processOrders(int minorder, int maxorder) {
	for (int order=minorder; order<=maxorder; order+=maxthreads) {
        spawnthreads(order, maxorder, maxthreads);
        waitthreads(order, maxorder, maxthreads);
	}
	return true;
}

int main(int argc, char *argv[])
{
	args.AddRequiredArgument("outputSpectraFile", outputSpectraFile, "Output file name");
	args.AddRequiredArgument("inputImage", inputImage, "Input FITS image to extract spectrum");
	args.AddRequiredArgument("inputGainFile", inputgain, "Input noise/gain file");
	args.AddRequiredArgument("inputGeometryFile", inputgeom, "Input geometry file");
	args.AddRequiredArgument("inputInstrumentProfileFile", inputprof, "Input instrument profile file");
	args.AddRequiredArgument("inputApertureFile", inputaper, "Input extraction aperture file");
	args.AddRequiredArgument("masterbias", masterbias, "FITS image with masterbias");
    args.AddRequiredArgument("masterflat", masterflat, "FITS image with masterflat");
	args.AddOptionalArgument("normalizedflat", normalizedflatfile, "", "FITS image with normalized flat-field");
	args.AddRequiredArgument("badpixelmask", badpixelmask, "FITS image with badpixel mask");
    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    args.AddOptionalArgument("effectiveApertureFraction", effectiveApertureFraction, 0.99, "Fraction of aperture width to set spectral element size");
    args.AddOptionalArgument("biasConstantToAdd", biasConstantToAdd, 0.0, "Bias constant to add -- to avoid negative numbers in extracted spectrum");
    args.AddOptionalArgument("backgroundBinsize", backgroundBinsize, 1, "Number of points to bin for IP measurements");
    args.AddOptionalArgument("minsigmaclip", minsigmaclip, 25, "Minimum variance threshold for optimal extraction");
    args.AddOptionalArgument("sigmacliprange", sigmacliprange, 6, "Maximum range of points to filter in single iteration");
    args.AddOptionalArgument("iterations", iterations, 3, "Number iterations for optimal extraction");
    args.AddOptionalArgument("onTargetProfile", onTargetProfile, false, "Measure spatial profile on-target instead of using flat-field");
    args.AddOptionalArgument("usePolynomialFit", usePolynomialFit, false, "Use polynomial instead of median for first measurement of profile");
    args.AddOptionalArgument("removeBackground", removeBackground, false, "Remove background. May be turned Off if target is bright.");
    args.AddOptionalArgument("starplusskymode", starplusskymode, false, "Star+sky: main flux is the sum of right beams minus sum of left beams.");
    args.AddOptionalArgument("starplusskyInvertSkyFiber", starplusskyInvertSkyFiber, false, "Star+sky: invert sky fiber (default is beam[0]=star and beam[1]=sky).");
	args.AddOptionalArgument("skyOverStarFiberAreaRatio", skyOverStarFiberAreaRatio, 1.0, "Area of the sky fiber over the area of the star fiber.");
    args.AddOptionalArgument("noCrossCorrelation", noCrossCorrelation, false, "Turn-off the cross-correlation calculation");
    args.AddOptionalArgument("rejectBadpixInRawExtraction", rejectBadpixInRawExtraction, false, "Apply raw extraction but reject badpixels");
    args.AddRequiredArgument("spectrumtype", spectralOrderType_val, "Method for extraction: 5 = Raw Flux Sum (default), 6 = Standard Flux, 7 = Optimal Extraction, 8 = OPERA Optimal Extraction");
    args.AddRequiredArgument("spectrumtypename", spectrumtypename, "Spectrum type name for verbose mode");
    args.AddRequiredArgument("maxthreads", maxthreads, "Maximum number of threads");
    args.AddPlotFileArguments(plotfilename, datafilename, scriptfilename, interactive);
	
	try {
		args.Parse(argc, argv);
		
		spectralOrderType = operaSpectralOrder_t(spectralOrderType_val);
		
		// we need an image...
		if (inputImage.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a gain file...
		if (inputgain.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}        
		// we need a geometry file...
		if (inputgeom.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
		// we need an instrument profile file...
		if (inputprof.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}	        
		// we need a aperture file...
		if (inputaper.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}	
		// we need a master flat...
		//if (masterflat.empty()) {
		//	throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		//}
		
		if (args.verbose) {
			cout << "operaExtraction: inputImage = " << inputImage << endl; 
			cout << "operaExtraction: outputSpectraFile = " << outputSpectraFile << endl;
			cout << "operaExtraction: inputGainFile = " << inputgain << endl;            
			cout << "operaExtraction: inputGeometryFile = " << inputgeom << endl; 			
			cout << "operaExtraction: inputInstrumentProfileFile = " << inputprof << endl;
			cout << "operaExtraction: inputApertureFile = " << inputaper << endl;            
			cout << "operaExtraction: masterbias = " << masterbias << endl; 	
			cout << "operaExtraction: masterflat = " << masterflat << endl; 	
			cout << "operaExtraction: normalizedflatfile = " << normalizedflatfile << endl; 				
			cout << "operaExtraction: badpixelmask = " << badpixelmask << endl;
			cout << "operaExtraction: spectrumtype = " << spectralOrderType << endl;
			cout << "operaExtraction: spectrumtypename = " << spectrumtypename << endl;
            cout << "operaExtraction: biasConstantToAdd = " << biasConstantToAdd << endl;
            cout << "operaExtraction: effectiveApertureFraction = " << effectiveApertureFraction << endl;
            cout << "operaExtraction: backgroundBinsize = " << backgroundBinsize << endl;
            cout << "operaExtraction: starplusskymode = " << starplusskymode << endl;
            cout << "operaExtraction: starplusskyInvertSkyFiber = " << starplusskyInvertSkyFiber << endl;
			cout << "operaExtraction: skyOverStarFiberAreaRatio = " << skyOverStarFiberAreaRatio << endl;
            cout << "operaExtraction: minsigmaclip = " << minsigmaclip << endl;
            cout << "operaExtraction: sigmacliprange = " << sigmacliprange << endl;
			cout << "operaExtraction: iterations = " << iterations << endl;
			cout << "operaExtraction: onTargetProfile = " << onTargetProfile << endl;
			cout << "operaExtraction: usePolynomialFit = " << usePolynomialFit << endl;
			cout << "operaExtraction: removeBackground = " << removeBackground << endl;
			cout << "operaExtraction: noCrossCorrelation = " << noCrossCorrelation << endl;
            cout << "operaExtraction: rejectBadpixInRawExtraction = " << rejectBadpixInRawExtraction << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaExtraction: ordernumber = " << ordernumber << endl;
            if(args.plot) {
                cout << "operaExtraction: plotfilename = " << plotfilename << endl;
                cout << "operaExtraction: datafilename = " << datafilename << endl;
                cout << "operaExtraction: scriptfilename = " << scriptfilename << endl; 
                cout << "operaExtraction: interactive = " << (interactive ? "YES" : "NO") << endl; 
            }            
		}
        
        ofstream fdata;
        if (!datafilename.empty()) fdata.open(datafilename.c_str());
        
        object = new operaFITSImage(inputImage, tfloat, READONLY);
		
        if (!masterflat.empty()){
            flat = new operaFITSImage(masterflat, tfloat, READONLY);
        } else {
            flat = new operaFITSImage(object->getnaxis1(),object->getnaxis2(),tfloat);
            *flat = 1.0;
        }
        
		if (!masterbias.empty()){
			bias = new operaFITSImage(masterbias, tfloat, READONLY);
            *bias = *bias - (float)biasConstantToAdd;
		} else {
            bias = new operaFITSImage(object->getnaxis1(),object->getnaxis2(),tfloat);
            *bias = 0.0;
        }
        
		if (!badpixelmask.empty()){
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(object->getnaxis1(),object->getnaxis2(),tfloat);
            *badpix = 1.0;
        }
        
		if (!normalizedflatfile.empty()){
			normalizedflat = new operaFITSImage(normalizedflatfile, tfloat, READONLY);
		} else {
            normalizedflat = new operaFITSImage(object->getnaxis1(),object->getnaxis2(),tfloat);
            *normalizedflat = 1.0;
        }

		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputgeom);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputaper);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputprof);

        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
		if (args.verbose) cout << "operaExtraction: minorder ="<< minorder << " maxorder=" << maxorder << endl;
		/*
		 * Add the order data to the raw spectrum product
		 */
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputgain);
        gainBiasNoise = spectralOrders.getGainBiasNoise();
  
        /* 
         * Uncomment below to introduce an artificial gain/noise
         *   to the second amplifier region. This may be wanted to
         *   simulate the effects of gain/noise in the reduction
         */
        // gainBiasNoise->setGain(1,gainBiasNoise->getGain(0)*50);
        // gainBiasNoise->setNoise(1,gainBiasNoise->getNoise(0));
        
        if (args.verbose) {
            for(unsigned amp=0;amp<gainBiasNoise->getAmps();amp++) {
                DATASEC_t datasec;
                gainBiasNoise->getDatasec(amp, datasec);
                
                double gain = gainBiasNoise->getGain(amp);
                double noise = gainBiasNoise->getNoise(amp);
                
                cout << "operaExtraction: "
                << " namps=" << gainBiasNoise->getAmps()
                << " amp="   << amp
                << " gain="  << gain
                << " noise=" << noise
                << " DataSec=[" << datasec.x1 << ":" << datasec.x2 <<","<< datasec.y1 << ":" << datasec.y2 << "]"
                << endl;
            }
		}
        
        unsigned long nthreads = maxorder+1;
        threads = (pthread_t *)calloc(nthreads, sizeof(pthread_t*));
        thread_args = (thread_args_t *)calloc(nthreads, sizeof(thread_args_t));
        
        if (maxthreads > 1) {
            processOrders(minorder, maxorder);
        } else {
            for (int order=minorder; order<=maxorder; order++) {
                processSingleOrder(order);
            }
        }
        
        unsigned NumberofBeams = spectralOrders.GetSpectralOrder(minorder)->getnumberOfBeams(); // for plotting
        
        for (int order=minorder; order<=maxorder; order++) {
           if (fdata.is_open()) {
               operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
               for(unsigned slitview=0; slitview<2; slitview++) {
                    spectralOrder->printBeamSpectrum(itos(slitview),&fdata);
                }
                fdata << endl;
            }
        }
		// output a spectrum...
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputSpectraFile, spectralOrderType);
        
		object->operaFITSImageClose();
		
        flat->operaFITSImageClose();
        
        if(bias) delete bias;
        if(badpix) delete badpix;
        if(normalizedflat) delete normalizedflat;
        if(object) delete object;
        if(flat) delete flat;
        
        if (fdata.is_open()) {
            fdata.close();
            if (!scriptfilename.empty()) {
                GenerateExtraction3DSpecPlot(scriptfilename,plotfilename,datafilename, NumberofBeams, interactive);
            }
        }             
	}
	catch (operaException e) {
		cerr << "operaExtraction: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaExtraction: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaExtraction: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/*
 * Generate 3D plot for spectra of atlas + comparison + identified lines
 */
void GenerateExtraction3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, unsigned numberOfBeams, bool display)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "set view 0,0" << endl;
    
    fgnu << "set palette gray" << endl;
    fgnu << "set palette gamma 2.0" << endl;
    fgnu << "set pm3d map" << endl;
    fgnu << "unset ztics" << endl;
    fgnu << "set cblabel \"flux\"" << endl;

    fgnu << "set xrange[-200:*]" << endl;
   
    fgnu << "\nset xlabel \"distance (pixels)\"" << endl;
    fgnu << "set ylabel \"order number\"" << endl;
    
    unsigned fluxColumnForFirstBeam = 13;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << endl;
        
        fgnu << "splot \"" << datafilename << "\" u 6:($2 + 0.4*$1 - 0.8 + 0.25 + 0.125):" << fluxColumnForFirstBeam <<" w pm3d";
        for(unsigned beam=1; beam<numberOfBeams; beam++) {
            unsigned fluxcol = fluxColumnForFirstBeam + 4*beam;
            fgnu << ",\"\" u 6:($2 + 0.4*$1 - 0.35 + 0.25 + 0.125):" << fluxcol << " w pm3d";
        }
        fgnu << endl;
    
        if (display) {
            fgnu << "\nset terminal x11" << endl;
            fgnu << "set output" << endl;
            fgnu << "replot" << endl;
        } else {
            fgnu << "\n#set terminal x11" << endl;
            fgnu << "#set output" << endl;
            fgnu << "#replot" << endl;
        }
    } else {
        fgnu << endl;
        
        fgnu << "splot \"" << datafilename << "\" u 6:($2 + 0.4*$1 - 0.8 + 0.25 + 0.125):" << fluxColumnForFirstBeam <<" w pm3d";
        for(unsigned beam=1; beam<numberOfBeams; beam++) {
            unsigned fluxcol = fluxColumnForFirstBeam + 4*beam;
            fgnu << ",\"\" u 6:($2 + 0.4*$1 - 0.35 + 0.25 + 0.125):" << fluxcol << " w pm3d";
        }
        fgnu << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}
