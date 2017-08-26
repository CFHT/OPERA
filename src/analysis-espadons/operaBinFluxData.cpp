/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaBinFluxData
 Version: 1.0
 Description: Bin flux data
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
#include "libraries/operaException.h"
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaLibCommon.h"
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/operaStats.h"                   // for operaArrayIndexSort
#include "libraries/Polynomial.h"
#include "libraries/operaFFT.h"

#include "analysis-espadons/operaBinFluxData.h"

#ifndef MAXPOINTSINOUTPUTSPECTRUM
#define MAXPOINTSINOUTPUTSPECTRUM 500000
#endif

#define NOTPROVIDED -999

/*! \brief operaBinFluxData */
/*! \file operaBinFluxData.cpp */
/*! \package operaBinFluxData */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*! 
 * operaBinFluxData
 * \author Eder Martioli
 * \brief Bin polar data.
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

static void printUsageSyntax(char * modulename);

void GenerateBinFluxPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string binnedDataFilename, bool display);

int main(int argc, char *argv[])
{
	int opt;
	
	string inputSpectrum;
	string inputWaveFile;
	string inputGeomFile;
	string outputBinnedSpectrum;
    
    operaSpectralOrder_t spectralOrderType = CalibratedRawSpectrum;
	string spectrumtypename;
	
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;    
    bool maxorderprovided = false;            
	
    unsigned binsize = 100;
    
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string spectrumDataFilename;
	string scriptfilename;	
	
	struct option longopts[] = {
		{"inputSpectrum",	1, NULL, 'i'},      // polarimetry spectrum *p.s  
		{"inputWaveFile",1, NULL, 'w'},             // wavelength calibration file (.wcal)
		{"inputGeomFile",		1, NULL, 'g'},      // geometry calibration file (.wcal)   
		{"outputBinnedSpectrum",	1, NULL, 'o'},   
		{"spectrumtype",		1, NULL, 'T'},	
		{"spectrumtypename",	1, NULL, 'N'},        
		{"ordernumber",			1, NULL, 'O'},	
		{"minorder",			1, NULL, 'M'},
		{"maxorder",			1, NULL, 'X'}, 
		{"binsize",				1, NULL, 'b'},    
		{"plotfilename",		1, NULL, 'P'},
		{"spectrumDataFilename",	1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},  
		{"interactive",			0, NULL, 'I'},
		{"plot",				optional_argument, NULL, 'p'},       
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:w:g:o:T:N:O:M:X:b:P:F:B:S:I:p::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{    
			case 'i':
				inputSpectrum = optarg;	
				break;
			case 'w':
				inputWaveFile = optarg;
				break;
			case 'g':
				inputGeomFile = optarg;
				break;
			case 'o':		// output
				outputBinnedSpectrum = optarg;
				break;
			case 'T':		// spectrum type
				spectralOrderType = (operaSpectralOrder_t)atoi(optarg);
				break;
			case 'N':		// spectrum type name for verbose mode
				spectrumtypename = optarg;
				break;                		                
			case 'O':
				ordernumber = atoi(optarg);
				break;				
			case 'M':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;  
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
				break;     
			case 'b':		// binsize
				binsize = atoi(optarg);
				break;                 
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break; 		                
			case 'F':
				spectrumDataFilename = optarg;
				break; 	                   
			case 'S':
				scriptfilename = optarg;
				break;  
			case 'I':		// for interactive plots
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
		// we need an input uncalibrated spectrum...
		if (inputSpectrum.empty()) {
			throw operaException("operaBinFluxData: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}		// we need an input uncalibrated spectrum...
		if (inputWaveFile.empty()) {
			throw operaException("operaBinFluxData: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (inputGeomFile.empty()) {
			throw operaException("operaBinFluxData: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output calibrated file name...
		if (outputBinnedSpectrum.empty()) {
			throw operaException("operaBinFluxData: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cout << "operaBinFluxData: inputSpectrum = " << inputSpectrum << endl;
			cout << "operaBinFluxData: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaBinFluxData: inputGeomFile = " << inputGeomFile << endl;
			cout << "operaBinFluxData: output Continuum Spectrum = " << outputBinnedSpectrum << endl;
			cout << "operaBinFluxData: spectrum type = " << spectralOrderType << endl;		
			cout << "operaBinFluxData: spectrum type name= " << spectrumtypename << endl;	            
            if(ordernumber != NOTPROVIDED) {
                cout << "operaBinFluxData: ordernumber = " << ordernumber << endl;            
            }   
            cout << "operaBinFluxData: binsize = " << binsize << endl;              
            if(plot) {
                cout << "operaBinFluxData: plotfilename = " << plotfilename << endl;
                cout << "operaBinFluxData: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaBinFluxData: scriptfilename = " << scriptfilename << endl; 
                if(interactive) {
                    cout << "operaBinFluxData: interactive = YES" << endl; 
                } else {
                    cout << "operaBinFluxData: interactive = NO" << endl; 
                }
            }            
            
		}
        
        ofstream *fspecdata = NULL;
        ofstream *fbinneddata = NULL;
        
        if (!spectrumDataFilename.empty()) {
            fspecdata = new ofstream();
            fspecdata->open(spectrumDataFilename.c_str());  
        }    
        
        if (!outputBinnedSpectrum.empty()) {
            fbinneddata = new ofstream();
            fbinneddata->open(outputBinnedSpectrum.c_str());  
        }              
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputSpectrum);
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputGeomFile); // read geometry calibration
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile); // read wavelength calibration
		
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
        
		if (verbose)
			cout << "operaBinFluxData: minorder ="<< minorder << " maxorder=" << maxorder << endl;        
        
        //unsigned NumberofBeams = 0; // for plotting
        float *meanwl = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianflux = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianfluxsig = new float[MAXPOINTSINOUTPUTSPECTRUM];
        int *originalorder = new int[MAXPOINTSINOUTPUTSPECTRUM];
        unsigned nTotalPoints = 0;
        
        float *flux_tmp = new float[binsize + 2];
        float *fluxvar_tmp = new float[binsize + 2];
        float *wavelength_tmp = new float[binsize + 2];
        
		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            
			if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
                // Below it populates SpectralElements wavelengths from calibration
                SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());

				unsigned length = SpectralElements->getnSpectralElements();
				
				if (fspecdata != NULL) {
					for(unsigned index=0;index<length;index++) {
						fspecdata->precision(6);
						*fspecdata << fixed;
						*fspecdata << order << "\t" << index << "\t"
						<< SpectralElements->getFlux(index) << "\t"
						<< SpectralElements->getFluxVariance(index) << endl;
					}
				}
				/*
				 *  Below is for binning
				 */
				unsigned NumberOfSamples = (unsigned)ceil((float)length/(float)binsize);
                
				for(unsigned k=0;k<NumberOfSamples;k++){
					unsigned firstPoint = k*binsize;
					unsigned lastPoint = (k+1)*binsize;
					if(lastPoint>length){
						lastPoint=length;
					}
					
					unsigned np=0;
					for(unsigned i=firstPoint;i<lastPoint;i++)
					{
						wavelength_tmp[np] = (float)SpectralElements->getwavelength(i);
						flux_tmp[np] = (float)SpectralElements->getFlux(i);
						fluxvar_tmp[np] = (float)SpectralElements->getFluxVariance(i);
						np++;
					}
					if(np > 3) {
                        meanwl[nTotalPoints] = operaArrayMean(np,wavelength_tmp);
                        medianflux[nTotalPoints] = operaArrayMedian(np,flux_tmp);
                        medianfluxsig[nTotalPoints] = operaArrayMedianSigma(np,flux_tmp,medianflux[nTotalPoints]);
                        originalorder[nTotalPoints] = order;
                        nTotalPoints++;
                    }
				}
			}
        }
        
        delete[] flux_tmp;
        delete[] fluxvar_tmp;
        delete[] wavelength_tmp;
        
        int *sindex = new int[nTotalPoints];
        
        operaArrayIndexSort((int)nTotalPoints,meanwl,sindex);
        
        for(unsigned index=0; index<nTotalPoints; index++) {
            if (fbinneddata != NULL) {
                fbinneddata->precision(6);
                *fbinneddata << fixed;
                *fbinneddata << originalorder[sindex[index]] << "\t"
                << sindex[index] << "\t"
                << index << "\t"
                << meanwl[sindex[index]] << "\t"
                << medianflux[sindex[index]] << "\t"
                << medianfluxsig[sindex[index]] << "\t"
                << endl;
            }
        }
        delete[] meanwl;
        delete[] medianflux;
        delete[] medianfluxsig;
        delete[] sindex;
        
        if (fbinneddata != NULL && fspecdata != NULL) {
            fbinneddata->close();
            fspecdata->close();
            if (!scriptfilename.empty()) {
                GenerateBinFluxPlot(scriptfilename,plotfilename,spectrumDataFilename,outputBinnedSpectrum,interactive);
            }
        }
        
	}
	catch (operaException e) {
		cerr << "operaBinFluxData: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaBinFluxData: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth] --param=<PARAMETER_VALUE_1> --param=<PARAMETER_VALUE_2> ... --output=<PRODUCT_FILE_NAME> --input=<INPUT_FILE_1> --input=<INPUT_FILE_2> ... \n\n"
	" Example: "+string(modulename)+"  -v -p 10 -p 2 --output=o.fits -i 001.fits -i 002.fits -i bad_pix.dat \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -p, --param=<PARAMETER_VALUE>, Input parameters \n"
	"  -o, --output=<PRODUCT_FILE_NAME>, Output product file \n"
	"  -i, --input=<INPUT_FILE_NAME>, Input files  \n\n";
}

void GenerateBinFluxPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string binnedDataFilename, bool display)
{
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str());  // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;
    *fgnu << "\nset xlabel \"wavelength (nm)\"" << endl;
    *fgnu << "set ylabel \"degree of polarization (%)\"" << endl;
    
    *fgnu << "set pointsize 0.5" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nplot \"" << spectrumDataFilename << "\" u (7575*$1 + $2):($3*100) w d" <<
        ",\"" << binnedDataFilename << "\" u (7575*$1 + $3):($4*100):($5*100) w yerr pt 7" <<
        ",\"" << binnedDataFilename << "\" u (7575*$1 + $3):($4*100) w l" <<
        ",\"" << binnedDataFilename << "\" u (7575*$1 + $3):($6*100) w l" <<
        ",\"" << binnedDataFilename << "\" u (7575*$1 + $3):($6*100):($7*100) w yerr pt 7" << endl;
        
        if (display) {
            *fgnu << "\nset terminal x11" << endl;
            *fgnu << "set output" << endl;
            *fgnu << "replot" << endl;
        } else {
            *fgnu << "\n#set terminal x11" << endl;
            *fgnu << "#set output" << endl;
            *fgnu << "#replot" << endl;
        }
    } else {
        *fgnu << "\nplot \"" << spectrumDataFilename << "\" u (7575*$1 + $2):($3*100) w d" <<
        ",\"" << binnedDataFilename << "\" u (7575*$1 + $3):($4*100):($5*100) w yerr pt 7" <<
        ",\"" << binnedDataFilename << "\" u (7575*$1 + $3):($4*100) w l" <<
        ",\"" << binnedDataFilename << "\" u (7575*$1 + $3):($6*100) w l" <<
        ",\"" << binnedDataFilename << "\" u (7575*$1 + $3):($6*100):($7*100) w yerr pt 7" << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());
    }
}

