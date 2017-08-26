/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCreateFlatFieldFluxCalibration
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

#include <fstream>
#include "libraries/operaIOFormats.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

/*! \file operaCreateFlatFieldFluxCalibration.cpp */

using namespace std;

void GenerateCreateFlatFieldFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string continuumDataFilename, unsigned NumberofBeams, bool display);

double getMedianValueInRange(operaSpectralElements *inputElements, unsigned centerElement, unsigned binsize);

/*! 
 * operaCreateFlatFieldFluxCalibration
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
	operaArgumentHandler args;
	
	string inputMasterFlatSpectrum;     // extracted masterflat spectrum (*.e.gz file)
	string outputFluxCalibrationFile;   // fluxCalibration (*.fcal.gz file)
	string wavelengthCalibration;
    double wavelengthForNormalization = 0.0;
    unsigned binsize = 30;
    int ordernumber = NOTPROVIDED;
    int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;
    string plotfilename;
	string spectrumDataFilename;
	string continuumDataFilename;
	string scriptfilename;
	bool interactive = false;
	
	args.AddRequiredArgument("inputMasterFlatSpectrum", inputMasterFlatSpectrum, "Spectrophotometric standard extracted uncalibrated spectrum");
	args.AddRequiredArgument("outputFluxCalibrationFile", outputFluxCalibrationFile, "Output flux calibration conversion file");
	args.AddRequiredArgument("wavelengthCalibration", wavelengthCalibration, "Input wavelength calibration file");
	args.AddRequiredArgument("wavelengthForNormalization", wavelengthForNormalization, "Choose a wavelength for which the spectrum should be normalized");
	args.AddRequiredArgument("binsize", binsize, "Number of points to bin for continuum estimate");
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
	args.AddOptionalArgument("plotfilename", plotfilename, "", "Output plot eps file name");
	args.AddOptionalArgument("spectrumDataFilename", spectrumDataFilename, "", "Output data file name");
	args.AddOptionalArgument("continuumDataFilename", continuumDataFilename, "", "Output data file name");
	args.AddOptionalArgument("scriptfilename", scriptfilename, "", "Output gnuplot script file name");
	args.AddSwitch("interactive", interactive, "For interactive plots");
	
	try {
		args.Parse(argc, argv);
		// we need an input uncalibrated spectrum...
		if (inputMasterFlatSpectrum.empty()) {
			throw operaException("operaCreateFlatFieldFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an output...
		if (outputFluxCalibrationFile.empty()) {
			throw operaException("operaCreateFlatFieldFluxCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (wavelengthCalibration.empty()) {
			throw operaException("operaCreateFlatFieldFluxCalibration: wcal: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

		if (args.verbose) {
			cout << "operaCreateFlatFieldFluxCalibration: inputMasterFlatSpectrum = " << inputMasterFlatSpectrum << endl;
            cout << "operaCreateFlatFieldFluxCalibration: outputFluxCalibrationFile = " << outputFluxCalibrationFile << endl;
            cout << "operaCreateFlatFieldFluxCalibration: wavelengthCalibration = " << wavelengthCalibration << endl;
            cout << "operaCreateFlatFieldFluxCalibration: wavelengthForNormalization = " << wavelengthForNormalization << endl;
            cout << "operaCreateFlatFieldFluxCalibration: binsize = " << binsize << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaCreateFlatFieldFluxCalibration: ordernumber = " << ordernumber << endl;
            if(args.plot) {
                cout << "operaCreateFlatFieldFluxCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaCreateFlatFieldFluxCalibration: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaCreateFlatFieldFluxCalibration: continuumDataFilename = " << continuumDataFilename << endl;
                cout << "operaCreateFlatFieldFluxCalibration: scriptfilename = " << scriptfilename << endl;
				cout << "operaCreateFlatFieldFluxCalibration: interactive = " << (interactive ? "YES" : "NO")  << endl;
            }            
            
		}
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputMasterFlatSpectrum);
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, wavelengthCalibration);

        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
		if (args.verbose) cout << "operaCreateFlatFieldFluxCalibration: minorder ="<< minorder << " maxorder=" << maxorder << endl;
  
        double maxFluxForNormalization = 0.0; //  DT May 08 2014 
        double maxBeamFluxForNormalization[MAXNUMBEROFBEAMS];
        unsigned numberOfBeams = 0;
        unsigned elemIndexPicked = 0;
        int refOrder = -999;
        
        // First loop over orders to search for order with maximum number of elements
        // --> maxNElements
        spectralOrders.setWavelengthsFromCalibration(minorder, maxorder);
        unsigned maxNElements = spectralOrders.getMaxNumberOfElementsInOrder(minorder, maxorder);
        
        // Calculate reference order and maximum flux for normalization
        // --> refOrder, maxFluxForNormalization
        if(wavelengthForNormalization) {
            // If wavelength for normalization is given, then use max flux at wavelength
            // --> refOrder, maxFluxForNormalization
            int *orderWithReferenceFluxForNormalization = new int[MAXORDERS];
            unsigned *elemIndexWithReferenceFluxForNormalization = new unsigned[MAXORDERS];
            unsigned nOrdersPicked = spectralOrders.getElemIndexAndOrdersByWavelength(orderWithReferenceFluxForNormalization,elemIndexWithReferenceFluxForNormalization,wavelengthForNormalization);

            for(unsigned i=0;i<nOrdersPicked;i++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(orderWithReferenceFluxForNormalization[i]);

                if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                    
                    operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                    double tmp_FluxForNormalization = SpectralElements->getFlux(elemIndexWithReferenceFluxForNormalization[i]);
                    
                    if(maxFluxForNormalization < tmp_FluxForNormalization) {
                        maxFluxForNormalization = tmp_FluxForNormalization;
                        refOrder = orderWithReferenceFluxForNormalization[i];
                        elemIndexPicked = elemIndexWithReferenceFluxForNormalization[i];
                        numberOfBeams = spectralOrder->getnumberOfBeams();
                        
                        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                            operaSpectralElements *BeamElements = spectralOrder->getBeamElements(beam);
                            maxBeamFluxForNormalization[beam] = BeamElements->getFlux(elemIndexPicked);
                        }   
                    }
                }
            }
            
        } else {
            float *flux_tmp = new float[maxNElements];

            for (int order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                
                if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                    operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();

                    unsigned NumberofElementsToBin = binsize; 
                    unsigned NumberOfElementSamples = (unsigned)ceil((float)SpectralElements->getnSpectralElements()/(float)NumberofElementsToBin);
                    
                    for(unsigned k=0;k<NumberOfElementSamples;k++){
                        unsigned firstElement = NumberofElementsToBin*(k);
                        unsigned lastElement =  NumberofElementsToBin*(k+1);
                        if (lastElement > SpectralElements->getnSpectralElements()){
                            lastElement = SpectralElements->getnSpectralElements();   
                        }
                        unsigned np=0;
                        for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
                            flux_tmp[np++] = (float)SpectralElements->getFlux(indexElem);
                        }
                        double medianFlux = (double)operaArrayMedian(np,flux_tmp);
                        
                        if(maxFluxForNormalization < medianFlux) {
                            maxFluxForNormalization = medianFlux;
                            refOrder = order;
                            elemIndexPicked = firstElement + (lastElement - firstElement)/2;
                            numberOfBeams = spectralOrder->getnumberOfBeams();
                            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                                operaSpectralElements *BeamElements = spectralOrder->getBeamElements(beam);
                                np=0;
                                for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
                                    flux_tmp[np++] = (float)BeamElements->getFlux(indexElem);
                                }
                                maxBeamFluxForNormalization[beam] = (double)operaArrayMedian(np,flux_tmp);
                            }
                        }
                    }
                }
            }
        }
 
        if(args.debug) {
            cout << "refOrder=" << refOrder << endl;
            cout << "maxFluxForNormalization=" << maxFluxForNormalization << endl;
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                cout << "maxBeamFluxForNormalization[" << beam << "]=" << maxBeamFluxForNormalization[beam] << endl;
            }
        }
       
        // Normalize flat-field spectrum by maxFluxForNormalization
        // to generate flat-field flux calibration spectrum
        // --> spectralEnergyDistribution
		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);

            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                if (args.verbose) {
					cout << "operaCreateFlatFieldFluxCalibration: Calibrating order number: "<< order << endl;
                }
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                operaWavelength *wavelength = spectralOrder->getWavelength();

                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();                
                spectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralElements *BeamElements = spectralOrder->getBeamElements(beam);
                    operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);
                    BeamSED->setwavelengthForNormalization(wavelengthForNormalization);
                }
                
                // normalize spectral Elements
				operaVector maxFluxesForNormalization;
				maxFluxesForNormalization.insert(maxFluxForNormalization);
				for (unsigned b = 0; b < numberOfBeams; b++) maxFluxesForNormalization.insert(maxBeamFluxForNormalization[b]);
				spectralOrder->normalizeSpectralElementsByConstant(maxFluxesForNormalization);

                for(unsigned b=0; b < spectralOrder->MainAndBeamCount(); b++) {
					operaSpectralEnergyDistribution& sed = spectralOrder->MainAndBeamSED(b);
					sed.setFluxCalibration(spectralOrder->MainAndBeamElements(b).getFluxVector().getflux());
					sed.setThroughput(spectralOrder->MainAndBeamElements(b).getFluxVector().getflux());
					sed.setCalibrationWavelength(spectralOrder->MainAndBeamElements(b).getWavelength());
					sed.setHasFluxCalibration(true);
					sed.setHasInstrumentThroughput(true);
				}
                spectralOrder->sethasSpectralEnergyDistribution(true);
                
                if(args.debug) {
                    for(unsigned i=0; i<SpectralElements->getnSpectralElements(); i++) {
                        cout << order << " " << i << " "
                        << SpectralElements->getwavelength(i) << " "
                        << SpectralElements->getdistd(i) << " "
                        << SpectralElements->getFlux(i) << " "
                        << spectralEnergyDistribution->getFluxCalibration().getflux(i) << " "
                        << spectralEnergyDistribution->getThroughput().getflux(i) << endl;
                    }
               }

			} //if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            
            if (!spectralOrder->gethasSpectralElements() ||  !spectralOrder->gethasWavelength()) {
				if (args.verbose) cout << "operaCreateFlatFieldFluxCalibration: Skipping order number: "<< order << " no spectral elements." << endl;
			}
        }

		// write output
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputFluxCalibrationFile, Fcal);
		
        if (!spectrumDataFilename.empty() && !continuumDataFilename.empty() && !scriptfilename.empty()) {
			GenerateCreateFlatFieldFluxCalibrationPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(),continuumDataFilename.c_str(), 2, interactive);
        }
        
	}
	catch (operaException e) {
		cerr << "operaCreateFlatFieldFluxCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCreateFlatFieldFluxCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
} 

void GenerateCreateFlatFieldFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string continuumDataFilename, unsigned NumberofBeams, bool display)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str());  // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "\nset xlabel \"wavelength (nm)\"" << endl;
    fgnu << "set ylabel \"flux\"" << endl;
    
    fgnu << "set pointsize 1.0" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nplot \"" << spectrumDataFilename << "\" u 6:7 w d" <<
        ",\"" << continuumDataFilename << "\" u 4:5 w linespoint lw 2.5" << endl;
        
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
        fgnu << "\nplot \"" << spectrumDataFilename << "\" u 6:7 w d" <<
        ",\"" << continuumDataFilename << "\" u 4:5 w linespoint lw 2.5" << endl;
        
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

double getMedianValueInRange(operaSpectralElements *inputElements, unsigned centerElement, unsigned binsize)
{
    unsigned nElements = inputElements->getnSpectralElements();
    
    int i1 = (int)centerElement - (int)binsize/2;
    int i2 = (int)centerElement + (int)binsize/2;
    
    if(i1 < 0) {
        i1 = 0;
        i2 = i1 + (int)binsize;
    }
    if((unsigned)i2 > nElements) {
        i2 = (int)nElements;
        i1 = (int)nElements - (int)binsize;
        if(i1 < 0) i1 = 0;
    }

    float *fluxInRange = new float[nElements];
    unsigned count = 0;
    for(unsigned i=(unsigned)i1;i<(unsigned)i2;i++) {
        fluxInRange[count++] = (float)inputElements->getFlux(i);
    }
    
    double medianFlux = (double)operaArrayMedian(count,fluxInRange);
    delete[] fluxInRange;
    
    return medianFlux;
}
