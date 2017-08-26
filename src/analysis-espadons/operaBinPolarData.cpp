/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaBinPolarData
 Version: 1.0
 Description: Bin polar data
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
#include "libraries/Polynomial.h"
#include "libraries/operaFFT.h"
#include "libraries/gzstream.h"

#define NOTPROVIDED -999

#ifndef MAXPOINTSINOUTPUTSPECTRUM
#define MAXPOINTSINOUTPUTSPECTRUM 500000
#endif

/*! \brief operaBinPolarData */
/*! \file operaBinPolarData.cpp */
/*! \package operaBinPolarData */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*! 
 * operaBinPolarData
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

void GenerateBinPolarPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string binnedDataFilename, bool display);

int main(int argc, char *argv[])
{
	int opt;
	
	string inputPolarSpectrum;
	string outputBinnedSpectrum;
    
	string inputWaveFile;    
    
    stokes_parameter_t StokesParameter = StokesV;
	
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
		{"inputPolarSpectrum",	1, NULL, 'i'},      // polarimetry spectrum *p.s  
		{"outputBinnedSpectrum",1, NULL, 'o'},   
		{"inputWaveFile",		1, NULL, 'w'}, 
        {"StokesParameter",		1, NULL, 's'},
		{"ordernumber",			1, NULL, 'O'},	
		{"minorder",			1, NULL, 'M'},
		{"maxorder",			1, NULL, 'X'}, 
		{"binsize",				1, NULL, 'b'},    
		{"plotfilename",		1, NULL, 'P'},
		{"spectrumDataFilename",1, NULL, 'F'},
		{"scriptfilename",		1, NULL, 'S'},  
		{"interactive",			0, NULL, 'I'},
		{"plot",				optional_argument, NULL, 'p'},       
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:w:s:O:M:X:b:P:F:B:S:I:p::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{    
			case 'i':
				inputPolarSpectrum = optarg;	
				break;
			case 'o':		// output
				outputBinnedSpectrum = optarg;
				break;
			case 'w':		// wavelength
				inputWaveFile = optarg;
				break;
			case 's':		// Stokes parameter
				StokesParameter = (stokes_parameter_t)atoi(optarg);
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
		if (inputPolarSpectrum.empty()) {
			throw operaException("operaBinPolarData: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}		// we need an input uncalibrated spectrum...
		// we need an output calibrated file name...
		if (outputBinnedSpectrum.empty()) {
			throw operaException("operaBinPolarData: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
        if (inputWaveFile.empty()) {
			throw operaException("operaBinPolarData: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (verbose) {
			cout << "operaBinPolarData: input polar spectrum = " << inputPolarSpectrum << endl;
			cout << "operaBinPolarData: output Continuum Spectrum = " << outputBinnedSpectrum << endl;
			cout << "operaBinPolarData: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaBinPolarData: StokesParameter = " << StokesParameter << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaBinPolarData: ordernumber = " << ordernumber << endl;            
            }   
            cout << "operaBinPolarData: binsize = " << binsize << endl;
            cout << "operaBinPolarData: plotfilename = " << plotfilename << endl;
            cout << "operaBinPolarData: spectrumDataFilename = " << spectrumDataFilename << endl;
            cout << "operaBinPolarData: scriptfilename = " << scriptfilename << endl;
            if(interactive) {
                cout << "operaBinPolarData: interactive = YES" << endl;
            } else {
                cout << "operaBinPolarData: interactive = NO" << endl;
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
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputPolarSpectrum);
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
			cout << "operaBinPolarData: minorder ="<< minorder << " maxorder=" << maxorder << endl;        
        
        //unsigned NumberofBeams = 0; // for plotting
		
        float *meanwl = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *meandistd = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianfluxStokesI = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianfluxStokesIsig = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianfluxStokesQUV = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianfluxStokesQUVsig = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianDegreeOfPolarization = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianDegreeOfPolarizationsig = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianNulldiff = new float[MAXPOINTSINOUTPUTSPECTRUM];
        float *medianNullratio = new float[MAXPOINTSINOUTPUTSPECTRUM];
        int *originalorder = new int[MAXPOINTSINOUTPUTSPECTRUM];
       
        unsigned nTotalPoints = 0;
        
        float *wavelength_tmp = new float[binsize + 2];
        float *distd_tmp = new float[binsize + 2];
        float *fluxStokesI_tmp = new float[binsize + 2];
        float *fluxStokesQUV_tmp = new float[binsize + 2];
        float *degreeOfPolarization_tmp = new float[binsize + 2];
        float *nulldiff_tmp = new float[binsize + 2];
        float *nullratio_tmp = new float[binsize + 2];
        
		for (int order=minorder; order<=maxorder; order++) {

			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                        
            if (spectralOrder->gethasPolarimetry() &&
                spectralOrder->gethasSpectralElements() &&
                spectralOrder->gethasWavelength()) {
                                
				operaPolarimetry *polarimetry = spectralOrder->getPolarimetry();
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                unsigned length = polarimetry->getLength();
                                
                SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());
                
				if (fspecdata != NULL) {
					for(unsigned index=0;index<length;index++) {
						fspecdata->precision(6);
						*fspecdata << fixed;
						*fspecdata << order << '\t'
                        << SpectralElements->getdistd(index) << '\t'
                        << SpectralElements->getwavelength(index) << '\t'
                        << polarimetry->getStokesParameterFlux(StokesI, index) << '\t';
                        
                        if (polarimetry->getHasStokesQ() && StokesParameter == StokesQ) {
                            *fspecdata << polarimetry->getStokesParameterFlux(StokesQ, index) << '\t'
                                       << polarimetry->getDegreeOfPolarizationFlux(StokesQ,index) << '\t';
                            if (polarimetry->getHasFirstNullPolarization())
                                *fspecdata << polarimetry->getFirstNullPolarizationFlux(StokesQ,index) << '\t';
                            if (polarimetry->getHasSecondNullPolarization())
                                *fspecdata << polarimetry->getSecondNullPolarizationFlux(StokesQ,index) << '\t';
                        } else if (polarimetry->getHasStokesU() && StokesParameter == StokesU) {
                            *fspecdata << polarimetry->getStokesParameterFlux(StokesU, index) << '\t'
                            << polarimetry->getDegreeOfPolarizationFlux(StokesU,index) << '\t';
                            if (polarimetry->getHasFirstNullPolarization())
                                *fspecdata << polarimetry->getFirstNullPolarizationFlux(StokesU,index) << '\t';
                            if (polarimetry->getHasSecondNullPolarization())
                                *fspecdata << polarimetry->getSecondNullPolarizationFlux(StokesU,index) << '\t';
                        } else if (polarimetry->getHasStokesV() && StokesParameter == StokesV) {
                            *fspecdata << polarimetry->getStokesParameterFlux(StokesV, index) << '\t'
                            << polarimetry->getDegreeOfPolarizationFlux(StokesV,index) << '\t';
                            if (polarimetry->getHasFirstNullPolarization())
                                *fspecdata << polarimetry->getFirstNullPolarizationFlux(StokesV,index) << '\t';
                            if (polarimetry->getHasSecondNullPolarization())
                                *fspecdata << polarimetry->getSecondNullPolarizationFlux(StokesV,index) << '\t';
                        }
                        *fspecdata << endl;
					}
                    *fspecdata << endl;
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
                        distd_tmp[np] = (float)SpectralElements->getdistd(i);
                        fluxStokesI_tmp[np] = (float)polarimetry->getStokesParameterFlux(StokesI, i);

                        if (polarimetry->getHasStokesQ() && StokesParameter == StokesQ) {
                            fluxStokesQUV_tmp[np] = (float)polarimetry->getStokesParameterFlux(StokesQ,i);
                            degreeOfPolarization_tmp[np] = (float)polarimetry->getDegreeOfPolarizationFlux(StokesQ,i);
                            if (polarimetry->getHasFirstNullPolarization())
                                nulldiff_tmp[np] = (float)polarimetry->getFirstNullPolarizationFlux(StokesQ,i);
                            if (polarimetry->getHasSecondNullPolarization())
                                nullratio_tmp[np] = (float)polarimetry->getSecondNullPolarizationFlux(StokesQ,i);
                        } else if (polarimetry->getHasStokesU() && StokesParameter == StokesU) {
                            fluxStokesQUV_tmp[np] = (float)polarimetry->getStokesParameterFlux(StokesU,i);
                            degreeOfPolarization_tmp[np] = (float)polarimetry->getDegreeOfPolarizationFlux(StokesU,i);
                            if (polarimetry->getHasFirstNullPolarization())
                                nulldiff_tmp[np] = (float)polarimetry->getFirstNullPolarizationFlux(StokesU,i);
                            if (polarimetry->getHasSecondNullPolarization())
                                nullratio_tmp[np] = (float)polarimetry->getSecondNullPolarizationFlux(StokesU,i);
                        } else if (polarimetry->getHasStokesV() && StokesParameter == StokesV) {
                            fluxStokesQUV_tmp[np] = (float)polarimetry->getStokesParameterFlux(StokesV,i);
                            degreeOfPolarization_tmp[np] = (float)polarimetry->getDegreeOfPolarizationFlux(StokesV,i);
                            if (polarimetry->getHasFirstNullPolarization())
                                nulldiff_tmp[np] = (float)polarimetry->getFirstNullPolarizationFlux(StokesV,i);
                            if (polarimetry->getHasSecondNullPolarization())
                                nullratio_tmp[np] = (float)polarimetry->getSecondNullPolarizationFlux(StokesV,i);
                        } else {
                            cout << "operaBinPolarData: WARNING: Stokes parameter not found. StokesParameter="<< StokesParameter << endl;
                        }
                        
						np++;
					}
                    
                    if(np > 3) {
                        meanwl[nTotalPoints] = operaArrayMean(np,wavelength_tmp);
                        meandistd[nTotalPoints] = operaArrayMean(np,distd_tmp);
                        medianfluxStokesI[nTotalPoints] = operaArrayMedian(np,fluxStokesI_tmp);
                        medianfluxStokesIsig[nTotalPoints] = operaArrayMedianSigma(np,fluxStokesI_tmp,medianfluxStokesI[nTotalPoints]);
                        medianfluxStokesQUV[nTotalPoints] = operaArrayMedian(np,fluxStokesQUV_tmp);
                        medianfluxStokesQUVsig[nTotalPoints] = operaArrayMedianSigma(np,fluxStokesQUV_tmp,medianfluxStokesQUV[nTotalPoints]);
                        medianDegreeOfPolarization[nTotalPoints] = operaArrayMedian(np,degreeOfPolarization_tmp);
                        medianDegreeOfPolarizationsig[nTotalPoints] = operaArrayMedianSigma(np,degreeOfPolarization_tmp,medianDegreeOfPolarization[nTotalPoints]);
                        medianNulldiff[nTotalPoints] = operaArrayMedian(np,nulldiff_tmp);
                        medianNullratio[nTotalPoints] = operaArrayMedian(np,nullratio_tmp);
                        originalorder[nTotalPoints] = order;
                        
                        nTotalPoints++;
                    }
				}
			}
        }
        
        int *sindex = new int[nTotalPoints];
        
        operaArrayIndexSort((int)nTotalPoints,meanwl,sindex);
                
        if (fbinneddata != NULL) {
            for(unsigned index=0; index<nTotalPoints; index++) {
                fbinneddata->precision(6);
                *fbinneddata << fixed;
                *fbinneddata << originalorder[sindex[index]] << '\t'
                << sindex[index] << '\t'
                << index << '\t'
                << meanwl[sindex[index]] << '\t'
                << meandistd[sindex[index]] << '\t'
                << medianfluxStokesI[sindex[index]] << '\t'
                << medianfluxStokesIsig[sindex[index]] << '\t'
                << medianfluxStokesQUV[sindex[index]] << '\t'
                << medianfluxStokesQUVsig[sindex[index]] << '\t'
                << medianDegreeOfPolarization[sindex[index]] << '\t'
                << medianDegreeOfPolarizationsig[sindex[index]] << '\t'
                << medianNulldiff[sindex[index]] << '\t'
                << medianNullratio[sindex[index]] << '\t'
                << endl;
            }
        }
        
        delete[] wavelength_tmp;
        delete[] distd_tmp;
        delete[] fluxStokesI_tmp;
        delete[] fluxStokesQUV_tmp;
        delete[] degreeOfPolarization_tmp;
        delete[] nulldiff_tmp;
        delete[] nullratio_tmp;
       
        delete[] meanwl;
        delete[] meandistd;
        delete[] medianfluxStokesI;
        delete[] medianfluxStokesIsig;
        delete[] medianfluxStokesQUV;
        delete[] medianfluxStokesQUVsig;
        delete[] medianDegreeOfPolarization;
        delete[] medianDegreeOfPolarizationsig;
        delete[] medianNulldiff;
        delete[] medianNullratio;
        delete[] originalorder;

        delete[] sindex;
        
        if (fbinneddata != NULL) {
            fbinneddata->close();
        }
        if (fspecdata != NULL) {
            fspecdata->close();
        }
        if (fbinneddata != NULL && fspecdata != NULL) {
            if (!scriptfilename.empty()) {
                GenerateBinPolarPlot(scriptfilename,plotfilename,spectrumDataFilename,outputBinnedSpectrum,interactive);
            }
        }
	}
	catch (operaException e) {
		cerr << "operaBinPolarData: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaBinPolarData: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth] --param=<PARAMETER_VALUE_1> --param=<PARAMETER_VALUE_2> ... --output=<PRODUCT_FILE_NAME> --input=<INPUT_FILE_1> --input=<INPUT_FILE_2> ... \n\n"
	" Example: "+string(modulename)+"  --inputPolarSpectrum=1599045i-new.p.gz --outputBinnedSpectrum=1599045i-outputbin.dat --binsize=4000 --spectrumDataFilename=1599045i-databin.dat --scriptfilename=1599045i-bin.gnu --StokesParameter=3 -v \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -p, --param=<PARAMETER_VALUE>, Input parameters \n"
	"  -o, --output=<PRODUCT_FILE_NAME>, Output product file \n"
	"  -i, --input=<INPUT_FILE_NAME>, Input files  \n\n";
}

void GenerateBinPolarPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string binnedDataFilename, bool display)
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
        
        *fgnu << "\nplot \"" << binnedDataFilename << "\" u 4:10 w l, \"\" u 4:12 w l, \"\" u 4:13 w l" << endl;
        //*fgnu << "\nreplot \"" << spectrumDataFilename << "\" u 3:6 w d" << endl;
        
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
        *fgnu << "\nplot \"" << binnedDataFilename << "\" u 4:10 w l, \"\" u 4:12 w l, \"\" u 4:13 w l" << endl;
        //*fgnu << "\nreplot \"" << spectrumDataFilename << "\" u 3:6 w d" << endl;

        
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

