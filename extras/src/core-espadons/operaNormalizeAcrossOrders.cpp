/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaNormalizeAcrossOrders
 Version: 1.0
 Description: Flux Calibration with Standard source 
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
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaLibCommon.h"
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/Polynomial.h"
#include "libraries/operaFFT.h"
#include "libraries/operaCCD.h"						// for MAXORDERS

#include "core-espadons/operaNormalizeAcrossOrders.h"

#define NOTPROVIDED -999
#define MAXFLUXREFERENCELENGTH 20000

/*! \file operaNormalizeAcrossOrders.cpp */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*! 
 * operaNormalizeAcrossOrders
 * \author Eder Martioli
 * \brief Normalize spectrum using inter-order information
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
	
	string inputUncalibratedSpectrum;
	string outputNormalizedSpectrum;
	string inputWaveFile;

    string spectrumtypename;
	operaSpectralOrder_t spectralOrderType = CalibratedOperaOptimalBeamSpectrum;

    double wavelength = 0;
    bool pickOnlyBrightestOrder = true;
    unsigned binsize = 100;
    int orderBin = 2;
    unsigned nsigcut = 3;
    
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;    
    bool maxorderprovided = false;            
        
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string spectrumDataFilename;
	string continuumDataFilename;	
	string scriptfilename;	
	
	struct option longopts[] = {
		{"inputUncalibratedSpectrum",	1, NULL, 'i'},
		{"outputNormalizedSpectrum",	1, NULL, 'o'},
		{"spectrumtype",		1, NULL, 'T'},
		{"spectrumtypename",	1, NULL, 'N'},        
		{"inputWaveFile",       1, NULL, 'w'},
		{"wavelength",		1, NULL, 'L'},
		{"binsize",				1, NULL, 'b'},
		{"orderBin",				1, NULL, 'B'},
		{"nsigcut",				1, NULL, 'c'},
		{"ordernumber",			1, NULL, 'O'},	
		{"minorder",			1, NULL, 'M'},
		{"maxorder",			1, NULL, 'X'},               
		{"plotfilename",		1, NULL, 'P'},
		{"spectrumDataFilename",		1, NULL, 'F'},
		{"continuumDataFilename",		1, NULL, 'C'},        
		{"scriptfilename",		1, NULL, 'S'},  
		{"interactive",			0, NULL, 'I'},
		{"plot",				optional_argument, NULL, 'p'},       
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:T:N:w:L:b:B:c:O:M:X:P:F:C:S:I:p::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputUncalibratedSpectrum = optarg;	
				break;
			case 'o':		// output
				outputNormalizedSpectrum = optarg;
				break;
			case 'T':		// spectrum type
				spectralOrderType = (operaSpectralOrder_t)atoi(optarg);
				break;
			case 'N':		// spectrum type name for verbose mode
				spectrumtypename = optarg;
				break;                
			case 'w':
				inputWaveFile = optarg;
				break;                
			case 'L':
				wavelength = atof(optarg);
				break;
			case 'b':	
				binsize = atoi(optarg);
				break;
			case 'B':	
				orderBin = atoi(optarg);
				break;
 			case 'c':		
				nsigcut = atoi(optarg);
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
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break; 		                
			case 'F':
				spectrumDataFilename = optarg;
				break; 	
			case 'C':
				continuumDataFilename = optarg;
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
		if (inputUncalibratedSpectrum.empty()) {
			throw operaException("operaNormalizeAcrossOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an output...
		if (outputNormalizedSpectrum.empty()) {
			throw operaException("operaNormalizeAcrossOrders: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaNormalizeAcrossOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

		if (verbose) {
			cout << "operaNormalizeAcrossOrders: input uncalibrated pectrum file = " << inputUncalibratedSpectrum << endl; 
            cout << "operaNormalizeAcrossOrders: outputNormalizedSpectrum = " << outputNormalizedSpectrum << endl;
			cerr << "operaNormalizeAcrossOrders: spectrumtype = " << spectralOrderType << endl;
			cerr << "operaNormalizeAcrossOrders: spectrumtypename = " << spectrumtypename << endl;
			cout << "operaNormalizeAcrossOrders: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaNormalizeAcrossOrders: wavelength = " << wavelength << " nm" << endl;
            cout << "operaNormalizeAcrossOrders: binsize = " << binsize << endl;
            cout << "operaNormalizeAcrossOrders: orderBin = " << orderBin << endl;
            cout << "operaNormalizeAcrossOrders: nsigcut = " << nsigcut << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaNormalizeAcrossOrders: ordernumber = " << ordernumber << endl;            
            }   
            if(plot) {
                cout << "operaNormalizeAcrossOrders: plotfilename = " << plotfilename << endl;
                cout << "operaNormalizeAcrossOrders: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaNormalizeAcrossOrders: continuumDataFilename = " << continuumDataFilename << endl;
                cout << "operaNormalizeAcrossOrders: scriptfilename = " << scriptfilename << endl; 
                if(interactive) {
                    cout << "operaNormalizeAcrossOrders: interactive = YES" << endl; 
                } else {
                    cout << "operaNormalizeAcrossOrders: interactive = NO" << endl; 
                }
            }            
            
		}
        ofstream *fspecdata = NULL;
        ofstream *fcontinuumdata = NULL;
        
        if (!spectrumDataFilename.empty()) {
            fspecdata = new ofstream();
            fspecdata->open(spectrumDataFilename.c_str());
        }    
        
        if (!continuumDataFilename.empty()) {
            fcontinuumdata = new ofstream();
            fcontinuumdata->open(continuumDataFilename.c_str());  
        }          
        
        unsigned numberOfprintouts = 1; // for 2D plotting
//        unsigned numberOfprintouts = 2; // for 3D plotting
        
		operaSpectralOrderVector spectralOrders(inputUncalibratedSpectrum);
        spectralOrders.ReadIntoSpectralOrders(inputWaveFile);

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
        
        int actualMinOrder = minorder;
        int actualMaxOrder = maxorder;
        
		if (verbose)
			cout << "operaNormalizeAcrossOrders: For continuum: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
            
        unsigned NumberofBeams = 0;
        for (int order=minorder; order<=maxorder; order++) {
            if (spectralOrders.GetSpectralOrder(minorder)->gethasSpectralElements()) {
                NumberofBeams = spectralOrders.GetSpectralOrder(minorder)->getnumberOfBeams();
                break;
            }
        }
        if (verbose) {
			cout << "operaNormalizeAcrossOrders: NumberofBeams = " << NumberofBeams << endl;
        }
        if(NumberofBeams == 0) {
            throw operaException("operaNormalizeAcrossOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
                
        /*
         *   E. Martioli Jul 2 2013.
         *   The function below measures the continuum using the normalization algorithm
         *   plus an additional fit to a robust line for the continuum across
         *   orders. The fit takes only a local sample of orders:  currentOrder +/- orderBin.
         *   The continuum data are populated in the spectralEnergyDistribution class. 
         *   This function also interpolates points to create the uncalibrated elements in 
         *   this class, so they won't need to be created again in the flux calibration routine
         */
        spectralOrders.measureContinuumAcrossOrders(binsize,orderBin,nsigcut);
        
                
        /*
         *  E. Martioli Jul 2 2013.
         *  The function below returns the number of orders (and size of output vectors) and 
         *  two vectors, one with order numbers and the other with element indexes, which
         *  are used as addresses to locate a given input wavelength in the extracted spectra. 
         *  The function parses all orders and all spectral elements using wavelength calibration to 
         *  search the right position in the image of a given wavelength. The reason it returns a 
         *  vector is because in echelle spectra more than one order may contain the same wavelength.
         */

        int *orderWithReferenceFluxForNormalization = new int[MAXORDERS];
        unsigned *elemIndexWithReferenceFluxForNormalization = new unsigned[MAXORDERS];
        unsigned nOrdersPicked = 0;
        int orderpicked = 0;

        
        if(wavelength) {
             nOrdersPicked = spectralOrders.getElemIndexAndOrdersByWavelength(orderWithReferenceFluxForNormalization,elemIndexWithReferenceFluxForNormalization,wavelength);
            
            for(unsigned i=0; i< nOrdersPicked; i++) {
                cout << "operaNormalizeAcrossOrders: wavelength=" << wavelength << " i=" << i << " order=" << orderWithReferenceFluxForNormalization[i] << " elemIndexWithReferenceFluxForNormalization = " <<  elemIndexWithReferenceFluxForNormalization[i] << endl;
            }
            
            
            if(pickOnlyBrightestOrder) {
                double uncalibratedContinuumFluxForNormalization = -BIG;
                for(unsigned i=0;i<nOrdersPicked;i++) {
                    double tmp_uncalibratedContinuumFluxForNormalization = spectralOrders.GetSpectralOrder(orderWithReferenceFluxForNormalization[i])->getSpectralEnergyDistribution()->getUncalibratedFluxElements()->getFlux(elemIndexWithReferenceFluxForNormalization[i]);
                    
                    if(uncalibratedContinuumFluxForNormalization < tmp_uncalibratedContinuumFluxForNormalization) {
                        uncalibratedContinuumFluxForNormalization = tmp_uncalibratedContinuumFluxForNormalization;
                        orderpicked = orderWithReferenceFluxForNormalization[i];
                    }
                }
                nOrdersPicked = 1;
                actualMinOrder = orderpicked;
                actualMaxOrder = orderpicked;
            } else {
                actualMinOrder = orderWithReferenceFluxForNormalization[0];
                actualMaxOrder = orderWithReferenceFluxForNormalization[nOrdersPicked-1];
            }
            
        }
        
		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            //operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            
            if(wavelength != 0 && (order < actualMinOrder || order > actualMaxOrder)) {
                spectralOrder->sethasSpectralElements(false);
            } else {
                if (spectralOrder->gethasSpectralEnergyDistribution()) {
                    if (verbose)
                        cout << "operaNormalizeAcrossOrders: Normalization for order number: "<< order << endl;
                    
                    spectralOrder->applyNormalizationFromExistingContinuum(fspecdata, fcontinuumdata, TRUE, TRUE, numberOfprintouts);
                }
                if (!spectralOrder->gethasWavelength()) {
                    if (verbose)
                        cout << "operaNormalizeAcrossOrders: Skipping order number: "<< order << " no wavelength calibration." << endl;
                }
            }
        }
		/*
		 * and write it out
		 */
		spectralOrders.WriteSpectralOrders(outputNormalizedSpectrum, spectralOrderType);
		
        
        if (fspecdata != NULL && fcontinuumdata != NULL) {
			fspecdata->close();
            fcontinuumdata->close();
            
            if (!scriptfilename.empty()) {
                GenerateNormalizeAcrossOrdersPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(),continuumDataFilename.c_str(), NumberofBeams, interactive);
            }
        }
        
	}
	catch (operaException e) {
		cerr << "operaNormalizeAcrossOrders: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaNormalizeAcrossOrders: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	" --inputUncalibratedSpectrum=<SPEC_FILE>"
	" --outputNormalizedSpectrum=<SPEC_FILE>"
	" --spectrumtype=<UNS_VALUE>"
	" --spectrumtypename=<STRING>"    
	" --inputWaveFile=<WAVE_FILE>"
	" --wavelength=<DOUBLE_VALUE>"
	" --binsize=<UNS_VALUE>"
	" --orderBin=<INT_VALUE>"
	" --nsigcut=<UNS_VALUE>"
	" --ordernumber=<UNS_VALUE>"
	" --minorder=<UNS_VALUE>"
	" --maxorder=<UNS_VALUE>"
	" --usePolynomial=<BOOL>"
	" --orderOfPolynomial=<UNS_VALUE>"
	" --generate3DPlot=<BOOL>"
	" --generateBeamPlot=<BOOL>"
	" --plotContinuum=<BOOL>"
	" --plotfilename=<EPS_FILE>"
	" --spectrumDataFilename=<DATA_FILE>"
	" --continuumDataFilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputUncalibratedSpectrum=1515007.e.gz --inputWaveFile=/Users/edermartioli/opera/calibrations/GalileanMoons/OLAPAa_pol_Normal.wcar.gz --outputNormalizedSpectrum=1515007n.s.gz --binsize=170 --orderBin=2 --nsigcut=3 --wavelength=656.28 -v \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -i, --inputUncalibratedSpectrum=<SPEC_FILE>,  Spectrophotometric standard extracted uncalibrated spectrum \n"
	"  -o, --outputFluxCalibrationFile=<SPEC_FILE>,  Output flux calibration conversion file \n"
	"  -T, --spectrumtype=<UNS_VALUE>, Option for spectrum type \n"
	"  -N, --spectrumtypename=<STRING>, Option for spectrum type \n"    
	"  -w, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file\n"
	"  -L, --wavelength=<DOUBLE_VALUE>, Wavelength (nm) to output order\n"
	"  -b, --binsize=<UNS_VALUE>, Number of points to bin for continuum estimate \n"
	"  -B, --orderBin=<INT_VALUE>, Number of neighbor order to perform local fit \n"
	"  -C, --nsigcut=<UNS_VALUE>, Cut off threshold for normalization algorithm in units of absdev \n"
	"  -O, --ordernumber=<UNS_VALUE>, Absolute order number to extract (default=all)\n"
	"  -M, --minorder=<UNS_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<UNS_VALUE>, Define maximum order number\n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -E, --generate3DPlot=<BOOL>, Switch to generate 3D or 2D plot spectra\n"
	"  -B, --generateBeamPlot=<BOOL>, Switch to generate plot of beams or full slit spectra\n"
	"  -c, --plotContinuum=<BOOL>, Switch to generate plot of continuum or normalized line spectra\n"
	"  -F, --spectrumDataFilename=<DATA_FILE>\n"
	"  -C, --continuumDataFilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateNormalizeAcrossOrdersPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string continuumDataFilename, unsigned NumberofBeams, bool display)
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
    *fgnu << "set ylabel \"flux\"" << endl;
    
    *fgnu << "set pointsize 1.0" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nplot \"" << spectrumDataFilename << "\" u 6:7 w d" <<
        ",\"" << continuumDataFilename << "\" u 4:5 w linespoint lw 2.5" << endl;
        
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
        *fgnu << "\nplot \"" << spectrumDataFilename << "\" u 6:7 w d" <<
        ",\"" << continuumDataFilename << "\" u 4:5 w linespoint lw 2.5" << endl;
        
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

