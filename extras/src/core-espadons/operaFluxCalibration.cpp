/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaFluxCalibration
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

#include "core-espadons/operaFluxCalibration.h"

#define NOTPROVIDED -999

#ifndef MAXPOINTSINOUTPUTSPECTRUM
#define MAXPOINTSINOUTPUTSPECTRUM 500000
#endif


/*! \file operaFluxCalibration.cpp */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;
        
/*! 
 * operaFluxCalibration
 * \author Eder Martioli
 * \brief Flux Calibration with Standard source.
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
	string outputCalibratedSpectrum;
	string inputWaveFile;    
	string inputfcal;

    operaSpectralOrder_t spectralOrderType = CalibratedRawBeamSpectrum;
	string spectrumtypename;
    
    double exposureTime = 1;   
    
    bool AbsoluteCalibration = false;
    int orderBin = 2;
    unsigned binsize = 210;
    
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;    
    bool maxorderprovided = false;            
    
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string spectrumDataFilename;
	string scriptfilename;	
	
	struct option longopts[] = {
		{"inputfcal",	1, NULL, 'i'},        
		{"inputUncalibratedSpectrum",	1, NULL, 'u'},
		{"inputWaveFile",       1, NULL, 'w'},        
		{"outputCalibratedSpectrum",	1, NULL, 'o'},   
		{"spectrumtype",		1, NULL, 'T'},	
		{"spectrumtypename",	1, NULL, 'N'},        
		{"exposureTime",	1, NULL, 'E'},
		{"binsize",	1, NULL, 'B'},
		{"orderBin",	1, NULL, 'b'},
		{"AbsoluteCalibration",	1, NULL, 'A'},
		{"ordernumber",			1, NULL, 'O'},
		{"minorder",			1, NULL, 'M'},
		{"maxorder",			1, NULL, 'X'},               
		{"plotfilename",		1, NULL, 'P'},
		{"spectrumDataFilename",		1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},  
		{"interactive",			0, NULL, 'I'},
		{"plot",				optional_argument, NULL, 'p'},       
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:u:w:o:T:N:E:B:b:A:O:M:X:P:F:S:I:p::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputfcal = optarg;	
				break;    
			case 'u':
				inputUncalibratedSpectrum = optarg;	
				break;
			case 'w':
				inputWaveFile = optarg;
				break;
			case 'o':		// output
				outputCalibratedSpectrum = optarg;
				break;
			case 'T':		// spectrum type
				spectralOrderType = (operaSpectralOrder_t)atoi(optarg);
				break;
			case 'N':		// spectrum type name for verbose mode
				spectrumtypename = optarg;
				break;                
			case 'E':
				exposureTime = atof(optarg);
				break;
			case 'B':
				binsize = atoi(optarg);
				break;
			case 'b':
				orderBin = atoi(optarg);
				break;
			case 'A':
				AbsoluteCalibration  = (atoi(optarg)?true:false);
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
		// we need an input opera flux calibration file...
		if (inputfcal.empty()) {
			throw operaException("operaFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an input uncalibrated spectrum...
		if (inputUncalibratedSpectrum.empty()) {
			throw operaException("operaFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}        
		// we need an output calibrated file name...
		if (outputCalibratedSpectrum.empty()) {
			throw operaException("operaFluxCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if (verbose) {
			cout << "operaFluxCalibration: input Flux Calibration .fcal file = " << inputfcal << endl; 
			cout << "operaFluxCalibration: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaFluxCalibration: input Uncalibrated Spectrum = " << inputUncalibratedSpectrum << endl;
			cout << "operaFluxCalibration: output Calibrated Spectrum file = " << outputCalibratedSpectrum << endl;
			cout << "operaCreateFluxCalibration: spectrum type = " << spectralOrderType << endl;		
			cout << "operaCreateFluxCalibration: spectrum type name = " << spectrumtypename << endl;	            
			cout << "operaFluxCalibration: exposure time = " << exposureTime << endl;
			cout << "operaFluxCalibration: binsize = " << binsize << endl;
			cout << "operaFluxCalibration: orderBin = " << orderBin << endl;
			cout << "operaFluxCalibration: AbsoluteCalibration = " << AbsoluteCalibration << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaFluxCalibration: ordernumber = " << ordernumber << endl;            
            }   
            if(plot) {
                cout << "operaFluxCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaFluxCalibration: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaFluxCalibration: scriptfilename = " << scriptfilename << endl; 
                if(interactive) {
                    cout << "operaFluxCalibration: interactive = YES" << endl; 
                } else {
                    cout << "operaFluxCalibration: interactive = NO" << endl; 
                }
            }            
            
		}
        ofstream *fspecdata = NULL;
        
        if (!spectrumDataFilename.empty()) {
            fspecdata = new ofstream();
            fspecdata->open(spectrumDataFilename.c_str());  
        }         
        
		operaSpectralOrderVector spectralOrders(inputUncalibratedSpectrum);
        spectralOrders.ReadIntoSpectralOrders(inputWaveFile);
        spectralOrders.ReadIntoSpectralOrders(inputfcal);
        
        unsigned NumberofBeams = 0;
        for (int order=minorder; order<=maxorder; order++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            if (spectralOrder->gethasSpectralElements()) {
                NumberofBeams = spectralOrder->getnumberOfBeams();
                break;
            }
        }
        
        unsigned nsigcut = 3;

        double uncalibratedContinuumFluxForNormalization = 0;
        double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS];
        
        spectralOrders.getContinuumFluxesForNormalization(&uncalibratedContinuumFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization,binsize, orderBin, nsigcut);

        double spectralBinConstant = exposureTime;
        double BeamSpectralBinConstant[MAXNUMBEROFBEAMS];
        for(unsigned beam=0; beam < NumberofBeams; beam++) {
            BeamSpectralBinConstant[beam] = exposureTime;
        }

        spectralOrders.ReadIntoSpectralOrders(inputfcal);
        
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
			cout << "operaFluxCalibration: minorder ="<< minorder << " maxorder=" << maxorder << endl;        
        
        float *wl = NULL;
        float *flux = NULL;
        float *variance = NULL;
        int *originalorder = NULL;
        int *sindex = NULL;
        float *beamFlux[MAXNUMBEROFBEAMS],*beamVariance[MAXNUMBEROFBEAMS];
        
        if (fspecdata != NULL) {
            wl = new float[MAXPOINTSINOUTPUTSPECTRUM];
            flux = new float[MAXPOINTSINOUTPUTSPECTRUM];
            variance = new float[MAXPOINTSINOUTPUTSPECTRUM];
            originalorder = new int[MAXPOINTSINOUTPUTSPECTRUM];
            for(unsigned beam = 0; beam < NumberofBeams; beam++) {
                beamFlux[beam] = new float[MAXPOINTSINOUTPUTSPECTRUM];
                beamVariance[beam] = new float[MAXPOINTSINOUTPUTSPECTRUM];                
            }
        }
        unsigned nTotalPoints = 0;
        
		for (int order=minorder; order<=maxorder; order++) {			
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            
            if (spectralOrder->gethasSpectralEnergyDistribution() && spectralOrder->gethasWavelength()) {

                // Below it populates SpectralElements wavelengths from calibration
                SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());

                if (verbose)
					cout << "operaFluxCalibration: Calibrating order number: "<< order << endl;
                
                if(debug) {                
                    // loop below print out the uncalibrated spectrum
                    for(unsigned i=0; i<SpectralElements->getnSpectralElements(); i++) {
                        cout << SpectralElements->getwavelength(i) << " " << SpectralElements->getFlux(i) << " " << SpectralElements->getFluxVariance(i) << endl;
                    } 
                }                

                spectralOrder->applyFluxCalibration(spectralBinConstant, BeamSpectralBinConstant, uncalibratedContinuumFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization, AbsoluteCalibration, NULL);
                
                if (fspecdata != NULL) {
                    for(unsigned indexElem=0; indexElem<SpectralElements->getnSpectralElements(); indexElem++) {
                        originalorder[nTotalPoints] = order;
                        wl[nTotalPoints] = SpectralElements->getwavelength(indexElem);
                        flux[nTotalPoints] = SpectralElements->getFlux(indexElem);
                        variance[nTotalPoints] = SpectralElements->getFluxVariance(indexElem);
                        for(unsigned beam = 0; beam < NumberofBeams; beam++) {
                            beamFlux[beam][nTotalPoints] = spectralOrder->getBeamElements(beam)->getFlux(indexElem);
                            beamVariance[beam][nTotalPoints] = spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem);
                        }
                        nTotalPoints++;
                    }
                }
			}
            
            if (!spectralOrder->gethasSpectralEnergyDistribution()) {
				if (verbose)
					cout << "operaFluxCalibration: Skipping order number: "<< order << " no flux calibration." << endl;
			}            
            if (!spectralOrder->gethasWavelength()) {
				if (verbose)
					cout << "operaFluxCalibration: Skipping order number: "<< order << " no wavelength calibration." << endl;
			}
        }
        
        /*
		 * and write it out
		 */
		spectralOrders.WriteSpectralOrders(outputCalibratedSpectrum, spectralOrderType);
        
        if (fspecdata != NULL) {
            sindex = new int[nTotalPoints];
            
            operaArrayIndexSort((int)nTotalPoints,wl,sindex);
            
            for(int index=0; index<(int)nTotalPoints; index++) {
                fspecdata->precision(6);
                *fspecdata << fixed;
                *fspecdata << originalorder[sindex[index]] << "\t"
                << sindex[index] << "\t"
                << index << "\t"
                << wl[sindex[index]] << "\t"
                << flux[sindex[index]] << "\t"
                << variance[sindex[index]];
                for(unsigned beam = 0; beam < NumberofBeams; beam++) {
                    *fspecdata << "\t" << beam
                               << "\t" << beamFlux[beam][sindex[index]]
                               << "\t" << beamVariance[beam][sindex[index]];
                }
                *fspecdata << endl;
            }

			fspecdata->close();
            if (!scriptfilename.empty()) {
                GenerateFluxCalibrationPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(), NumberofBeams, interactive);
            }
            delete[] wl;
            delete[] flux;
            delete[] variance;
            delete[] originalorder;
            delete[] sindex;
            for(unsigned beam = 0; beam < NumberofBeams; beam++) {
                delete[] beamFlux[beam];
                delete[] beamVariance[beam];
            }
        }
	}
	catch (operaException e) {
		cerr << "operaFluxCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaFluxCalibration: " << operaStrError(errno) << endl;
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
	" --inputfcal=<SPEC_FILE>"
    " --inputWaveFile=<WAVE_FILE>"
	" --outputCalibratedSpectrum=<SPEC_FILE>"
	" --spectrumtype=<UNS_VALUE>"
	" --spectrumtypename=<STRING>"
	" --exposureTime=<FLOAT_VALUE>"
	" --binsize=<UNS_VALUE>"
	" --AbsoluteCalibration=<BOOL>"
	" --ordernumber=<UNS_VALUE>"
	" --minorder=<UNS_VALUE>"
	" --maxorder=<UNS_VALUE>"
	" --plotfilename=<EPS_FILE>"
	" --generate3DPlot=<BOOL>"
    " --generateBeamPlot=<BOOL>"
	" --plotContinuum=<BOOL>"
	" --spectrumDataFilename=<DATA_FILE>"
	" --continuumDataFilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputfcal=master.fcal.gz --inputUncalibratedSpectrum=/Users/edermartioli/opera/spectra/GalileanMoons/1515004.e.gz --inputWaveFile=/Users/edermartioli/opera/calibrations/GalileanMoons/OLAPAa_pol_Normal.wcar.gz --outputCalibratedSpectrum=1515004.s.gz --spectrumDataFilename=1515004fcalib.spec --scriptfilename=1515004fcalib.gnu --exposureTime=1 -v \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -u, --inputUncalibratedSpectrum=<SPEC_FILE>,  Object uncalibrated spectrum \n"
	"  -i, --inputfcal=<SPEC_FILE>,  Input flux calibration conversion file \n"
    "  -w, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file\n"
	"  -o, --outputCalibratedSpectrum=<SPEC_FILE>,  Output calibrated object spectrum \n"
	"  -T, --spectrumtype=<UNS_VALUE>, Option for spectrum type \n"
	"  -N, --spectrumtypename=<STRING>, Option for spectrum type \n"
	"  -E, --exposureTime=<FLT_VALUE>, Exposure time of input uncalibrated spectrum\n"
	"  -B, --binsize=<UNS_VALUE>, Bin size for continuum flux calculation\n"
	"  -A, --AbsoluteCalibration=<BOOL>, Perform absolute flux calibration\n"
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


void GenerateFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, unsigned NumberofBeams, bool display)
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
    
    *fgnu << "set pointsize 0.5" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nplot \"" << spectrumDataFilename << "\" u 4:5 w l" << endl;
        
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
        *fgnu << "\nplot \"" << spectrumDataFilename << "\" u 4:5 w l" << endl;
        
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
