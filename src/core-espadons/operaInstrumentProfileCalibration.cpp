/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaInstrumentProfileCalibration
 Version: 1.0
 Description: Create the Instrument Profile 
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

#include <pthread.h>
#include <fstream>
#include "libraries/operaIOFormats.h"
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

/*! \file operaInstrumentProfileCalibration.cpp */

using namespace std;

/*!
 * operaInstrumentProfileCalibration
 * \author Eder Martioli
 * \brief Create the Instrument Profile.
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

void GenerateInstrumentProfile3DPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned minorderWithIP, unsigned maxorderWithIP, unsigned IPxsize, unsigned IPysize, bool display);

// These variables are global for thread support
operaArgumentHandler args;

// Input parameters
string outputprof;
string geometryfilename;
string masterbias;
string masterflat;
string mastercomparison;
string masterfabperot;
string badpixelmask;
string gainfilename;

int method = 1;
double spectralElementHeight = 1.0;
double referenceLineWidth = 2.0;
double LocalMaxFilterWidth = 2.5;
double DetectionThreshold = 0.2;
double MinPeakDepth = 1.5;
unsigned minimumLinesForIPMeasurements = 20;
unsigned binsize = 80;
double tilt = -3.0;
unsigned IPxsize = 0;
unsigned IPysize = 0;
unsigned IPxsampling = 0;
unsigned IPysampling = 0;
unsigned maxthreads = 1;
int ordernumber = NOTPROVIDED;
int minorder = 22;
int maxorder = 62;
string plotfilename;
string datafilename;
string scriptfilename;
bool interactive = false;

// Other variables used across threads
double LocalMaxFilterWidthInPixels;
double MinPeakDepthInElectronUnits;
double gain = 1.12;
double noise = 3.5;
bool fabperot;
operaFITSImage *badpix = NULL;
operaFITSImage *bias = NULL;
operaFITSImage *flat = NULL;
operaFITSImage *comp = NULL;

operaSpectralOrderVector spectralOrders;

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
    
    if (args.verbose) {
        cout << "operaInstrumentProfileCalibration: Processing order = " << order << endl;
    }
    // create pointer to current spectral order:
    operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
    
    // create a set of spectral elements based on given element height:
    spectralOrder->setSpectralElementsByHeight(spectralElementHeight);
    
    // create vector to store information of instrument profile for current spectral order:
    unsigned MaxNDataPoints = spectralOrder->getSpectralElements()->getnSpectralElements();
    
    unsigned sampleElementForPlot = spectralOrder->getSpectralElements()->getnSpectralElements()/2;
    
    spectralOrder->setInstrumentProfileVector(IPxsize,IPxsampling,IPysize,IPysampling,MaxNDataPoints);
    
    spectralOrder->measureInstrumentProfileAlongRowsInto2DWithGaussian(*flat,*badpix,binsize,referenceLineWidth,tilt,false,sampleElementForPlot,NULL,minimumLinesForIPMeasurements);
    
    if(spectralOrder->gethasInstrumentProfile()) {
        string methodName;
        double MaxContamination;
        double amplitudeCutOff;
        unsigned nSigCut;
        if (fabperot) {
            methodName = "Fabry-Perot";
            MaxContamination = 1.0;
            amplitudeCutOff = 3*noise;
            nSigCut = 3;
		} else {
            methodName = "Comparison";
            MaxContamination = 0.05;      // accept up to 1% flux contamination from neighbor line
            amplitudeCutOff = 2*noise;  // limit lines with amplitude greater than 10 x CCD noise
            nSigCut = 2;               // limit lines with width up to twice larger than the median deviation
        }    
		operaSpectralLines *spectralLines = NULL;
		try {
			spectralOrder->calculateXCorrBetweenIPandImage(*comp, *badpix, NULL);
			spectralOrder->setSpectralLines(*comp, *badpix, *bias, noise, gain, referenceLineWidth, DetectionThreshold, LocalMaxFilterWidthInPixels, MinPeakDepthInElectronUnits);
			spectralOrder->sethasSpectralLines(true);
			spectralLines = spectralOrder->getSpectralLines();
			if (args.verbose) cout << "operaInstrumentProfileCalibration: " << spectralLines->getnLines() << " lines found in order " << order << " of " << methodName << "." << endl;
		} catch (operaException e) {
			if (args.verbose) cout << "operaInstrumentProfileCalibration: No lines found in order " << order << " of " << methodName << "." << endl;
			spectralOrder->sethasInstrumentProfile(false);
		}
		try {
			if (spectralOrder->gethasSpectralLines() && spectralLines && spectralLines->getnLines() > 0) {
				if(method == 1) {
					spectralOrder->measureInstrumentProfileUsingWeightedMean(*comp, *badpix, MaxContamination, amplitudeCutOff, nSigCut, sampleElementForPlot, NULL);
				} else if (method == 2) {
					spectralOrder->measureInstrumentProfileUsingMedian(*comp, *badpix, MaxContamination, amplitudeCutOff, nSigCut, sampleElementForPlot, NULL,minimumLinesForIPMeasurements);
				} else if (method == 3) {
					spectralOrder->measureInstrumentProfile(*comp, *badpix, MaxContamination, amplitudeCutOff, nSigCut, sampleElementForPlot, NULL,minimumLinesForIPMeasurements);
				} else if (method == 4) {
					spectralOrder->measureInstrumentProfileWithBinning(*comp, *badpix, binsize, MaxContamination, amplitudeCutOff, nSigCut, sampleElementForPlot, NULL,minimumLinesForIPMeasurements);
				}
				spectralOrder->recenterOrderPosition();
			}
		} catch (operaException e) {
			cout << "operaInstrumentProfileCalibration: " << "non-fatal error in order " << order << endl;
			cout << e.getFormattedMessage() << endl;
		}
    }
    // clean up so we don't run out of memory
    spectralOrder->getInstrumentProfile()->deleteDataCubes();
    if (spectralOrder->gethasSpectralLines()) {
        delete spectralOrder->getSpectralLines();
        spectralOrder->sethasSpectralLines(false);
    }
    if (spectralOrder->gethasSpectralElements()) {
        delete spectralOrder->getSpectralElements();
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
	args.AddRequiredArgument("outputProf", outputprof, "Output instrument profile file");
	args.AddRequiredArgument("geometryfilename", geometryfilename, "Input geometry file");
	args.AddRequiredArgument("masterbias", masterbias, "Input master bias FITS image");
	args.AddRequiredArgument("masterflat", masterflat, "Input master flat-field FITS image");
	args.AddOptionalArgument("mastercomparison", mastercomparison, "", "Input master comparison (ThAr) FITS image (use this or masterfabperot but not both)");
	args.AddOptionalArgument("masterfabperot", masterfabperot, "", "Input master Fabry-Perot FITS image (use this or mastercomparison but not both)");
	args.AddOptionalArgument("badpixelmask", badpixelmask, "", "FITS image for badpixel mask");
	args.AddRequiredArgument("gainfilename", gainfilename, "Input gain/noise file");
	
	args.AddOptionalArgument("method", method, 1, "Method to combine IP measurements from spectral lines: 1 = weighted mean, 2 = median combine, 3 = polynomial fit, 4 = polynomial fit with median binning (use binsize provided)");
	args.AddOptionalArgument("spectralElementHeight", spectralElementHeight, 1.0, "Height of spectral element in Y-direction in pixel units");
	args.AddOptionalArgument("referenceLineWidth", referenceLineWidth, 2.5, "Spectral line width for reference in pixel units");
	args.AddOptionalArgument("LocalMaxFilterWidth", LocalMaxFilterWidth, 6.25, "In units of line width");
	args.AddOptionalArgument("DetectionThreshold", DetectionThreshold, 0.2, "Between 0 and 1");
	args.AddOptionalArgument("MinPeakDepth", MinPeakDepth, 5.25, "In units of noise");
	args.AddOptionalArgument("minimumlines", minimumLinesForIPMeasurements, 20, "Minimum number of lines for IP measurements");
	args.AddOptionalArgument("binsize", binsize, 80, "Number of points to bin for IP measurements");
	args.AddOptionalArgument("tilt", tilt, -3.0, "Tilt in degrees");
	args.AddOptionalArgument("xSize", IPxsize, 30, "IP Dimensions");
	args.AddOptionalArgument("ySize", IPysize, 6, "IP Dimensions");
	args.AddOptionalArgument("xSampling", IPxsampling, 5, "IP Dimensions");
	args.AddOptionalArgument("ySampling", IPysampling, 5, "IP Dimensions");
	args.AddOptionalArgument("maxthreads", maxthreads, 1, "Maximum number of threads");
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
	args.AddPlotFileArguments(plotfilename, datafilename, scriptfilename, interactive);
	
	try {
		args.Parse(argc, argv);
		
        if(method!=1 && method !=2 && method != 3 && method != 4) {
            method = 1;
            cout <<  "operaInstrumentProfileCalibration: Warning, undefined method. Using method = 1 (weighted mean)" << endl;	
        }
		// we need a geometry file...
		if (geometryfilename.empty()) {
			throw operaException("operaInstrumentProfileCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a masterbias...
		if (masterbias.empty()) {
			throw operaException("operaInstrumentProfileCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a masterflat...
		if (masterflat.empty()) {
			throw operaException("operaInstrumentProfileCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a comparison or fabperot, but not both...
		if ((mastercomparison.empty() && masterfabperot.empty()) || (!mastercomparison.empty() && !masterfabperot.empty())) {
			throw operaException("operaInstrumentProfileCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}								
		
		if (args.verbose) {
			cout << "operaInstrumentProfileCalibration: geometryfilename = " << geometryfilename << endl; 
			cout << "operaInstrumentProfileCalibration: outputProf = " << outputprof << endl;
			cout << "operaInstrumentProfileCalibration: masterbias = " << masterbias << endl;            
			cout << "operaInstrumentProfileCalibration: masterflat = " << masterflat << endl; 	
			cout << "operaInstrumentProfileCalibration: mastercomparison = " << mastercomparison << endl; 	
			cout << "operaInstrumentProfileCalibration: badpixelmask = " << badpixelmask << endl;     
			cout << "operaInstrumentProfileCalibration: method = " << method << endl;                 
			cout << "operaInstrumentProfileCalibration: masterfabperot = " << masterfabperot << endl; 			
			cout << "operaInstrumentProfileCalibration: spectralElementHeight = " << spectralElementHeight << endl;            
			cout << "operaInstrumentProfileCalibration: referenceLineWidth = " << referenceLineWidth << endl;
			cout << "operaInstrumentProfileCalibration: MinPeakDepth = " << MinPeakDepth << endl;
			cout << "operaInstrumentProfileCalibration: LocalMaxFilterWidth = " << LocalMaxFilterWidth << endl;
			cout << "operaInstrumentProfileCalibration: DetectionThreshold = " << DetectionThreshold << endl;
			cout << "operaInstrumentProfileCalibration: maxthreads = " << maxthreads << endl;
			if(ordernumber != NOTPROVIDED) cout << "operaInstrumentProfileCalibration: ordernumber = " << ordernumber << endl;            
            if(args.plot) {
                cout << "operaInstrumentProfileCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaInstrumentProfileCalibration: datafilename = " << datafilename << endl;
                cout << "operaInstrumentProfileCalibration: scriptfilename = " << scriptfilename << endl;                
            }            
		}
        
        ofstream fdata;
        if (!datafilename.empty()) fdata.open(datafilename.c_str());
        
        fabperot = !masterfabperot.empty();
		if (fabperot) comp = new operaFITSImage(masterfabperot, tfloat, READONLY);
		else comp = new operaFITSImage(mastercomparison, tfloat, READONLY);
        bias = new operaFITSImage(masterbias, tfloat, READONLY);
		flat = new operaFITSImage(masterflat, tfloat, READONLY);
		
		//flat -= bias; // remove bias from masterflat
		
		if (!badpixelmask.empty()){              
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(flat->getnaxis1(),flat->getnaxis2(),tfloat);
            *badpix = 1.0;
        }
        
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, geometryfilename);
		UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
		
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, gainfilename);
        unsigned amp = 0;
		gain = spectralOrders.getGainBiasNoise()->getGain(amp);
		noise = spectralOrders.getGainBiasNoise()->getNoise(amp);

        MinPeakDepthInElectronUnits = MinPeakDepth*noise;
        LocalMaxFilterWidthInPixels = LocalMaxFilterWidth*referenceLineWidth;
		
		if (args.verbose) {
			cout << "operaInstrumentProfileCalibration: binsize = " << binsize << endl;            
			cout << "operaInstrumentProfileCalibration: ipDimensions = [" << IPxsize << " pxl (" << IPxsampling << " ppp)," << IPysize <<" pxl ("<< IPysampling << " ppp)]" << endl; 
			cout << "operaInstrumentProfileCalibration: minorder = " << minorder << " maxorder = " << maxorder << endl;            
		}

        unsigned long nthreads = maxorder+1;
        threads = (pthread_t *)calloc(nthreads, sizeof(pthread_t*));
        thread_args = (thread_args_t *)calloc(nthreads, sizeof(thread_args_t));

        if (maxthreads > 1) processOrders(minorder, maxorder);
        else for (int order=minorder; order<=maxorder; order++) processSingleOrder(order);
        
        if (fdata.is_open()) {
            int minorderWithIP = NOTPROVIDED;
            int maxorderWithIP = NOTPROVIDED;
            
            for (int order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                if(spectralOrder->gethasGeometry() && spectralOrder->gethasInstrumentProfile()) {
                    operaGeometry *Geometry = spectralOrder->getGeometry();
                    operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
                    float distd = (float)Geometry->CalculateDistance(Geometry->getYmin(), (Geometry->getYmax() - Geometry->getYmin())/2);
                    instrumentProfile->printModel(distd,order,&fdata);
                    if(minorderWithIP == NOTPROVIDED){
                        minorderWithIP = (int)order;
                    }
                    maxorderWithIP = (int)order;
                }
            }
            
            fdata.close();
            if (!scriptfilename.empty()) {
                GenerateInstrumentProfile3DPlot(scriptfilename,plotfilename,datafilename,minorderWithIP,maxorderWithIP, IPxsize, IPysize, interactive);
            }
        }
		
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputprof, Prof);

		bias->operaFITSImageClose();
		flat->operaFITSImageClose();
		comp->operaFITSImageClose();
        
        if (badpix) delete badpix;
		if (bias) delete bias;
		if (flat) delete flat;
		if (comp) delete comp;
	}
	catch (operaException e) {
		cerr << "operaInstrumentProfileCalibration: "  << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaInstrumentProfileCalibration: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaInstrumentProfileCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

void GenerateInstrumentProfile3DPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned minorderWithIP, unsigned maxorderWithIP, unsigned IPxsize, unsigned IPysize, bool display)
{
	if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "set view 0,0" << endl;
    
    fgnu << "set palette color" << endl;
    fgnu << "set palette gamma 2.5" << endl;
    fgnu << "set pm3d map" << endl;
    fgnu << "unset ztics" << endl;
    fgnu << "set cblabel \"flux fraction\"" << endl;
    fgnu << endl;
    fgnu << "set xrange[-5:(5*" << IPxsize << "+5)]" << endl;
    fgnu << "set yrange[-5:(("<< (maxorderWithIP - minorderWithIP) <<"*" << IPysize << "/5)+5)]" << endl;
    
    for(unsigned order=minorderWithIP;order<=maxorderWithIP;order++) {
        double xlabelpos = ((double)order-((double)minorderWithIP+(double)(5*floor((order-minorderWithIP)/5))))*(double)IPxsize+1;
        double ylabelpos = (double)floor((order-minorderWithIP)/5)*(double)IPysize+1;
        fgnu << "set label \"" << order << "\" at " << xlabelpos << "," << ylabelpos << " front font \"Helvetica,7\"" << endl;
    }

    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nsplot \"" << dataFileName << "\" u (($1-("<<minorderWithIP <<"+5*int(($1-"<< minorderWithIP <<")/5)))*" << IPxsize << "+" << (float)IPxsize/2 << "+$4):(int(($1-"<<minorderWithIP <<")/5)*" << IPysize << "+(" << (float)IPysize/2 << ")+$5):7 with pm3d" << endl;

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
        fgnu << "\nsplot \"" << dataFileName << "\" u (($1-("<<minorderWithIP <<"+5*int(($1-"<< minorderWithIP <<")/5)))*" << IPxsize << "+" << (float)IPxsize/2 << "+$4):(int(($1-"<<minorderWithIP <<")/5)*" << IPysize << "+(" << (float)IPysize/2 << ")+$5):7 with pm3d" << endl;
        
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
