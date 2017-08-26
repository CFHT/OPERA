/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMasterFluxCalibration
 Version: 1.0
 Description: Create a Master Flux Calibration out of many calibrations. 
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
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

#define MAXNUMBEROFFLUXCALIBRATIONFILES 500
#define WAVELENGTH_PRECISION 0.01
#define MINELEMENTS 20

/*! \file operaMasterFluxCalibration.cpp */

using namespace std;

void GenerateMasterFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string outputDataFilename, unsigned NumberofBeams, bool display);
        
/*! 
 * operaMasterFluxCalibration
 * \author Eder Martioli
 * \brief Create a Master Flux Calibration out of many calibrations.
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
	
	string listofinputfcal;
	double inputconstant = 0.0;
    string outputfcal;
	string inputReferenceSpectrum;      // this is a reference opera spectrum that defines the format for the output
    string inputWaveFile;
    unsigned combineMethod = 1; // Not implemented yet. Eder - Jan/09/2013.
	int ordernumber = NOTPROVIDED;
    int minorder = 22;
    int maxorder = 62;    
    string plotfilename;	
	string spectrumDataFilename;
	string outputDataFilename;
	string scriptfilename;
	bool interactive = false;
	
	args.AddOptionalArgument("inputfcal", listofinputfcal, "", "List of input flux calibration files, separated by spaces");
	args.AddOptionalArgument("inputconstant", inputconstant, 0.0, "This value should be input when there are no fcal files available");
	args.AddRequiredArgument("outputfcal", outputfcal, "Output master flux calibration file");
	args.AddRequiredArgument("inputReferenceSpectrum", inputReferenceSpectrum, "Input reference spectrum file");
	args.AddRequiredArgument("inputWaveFile", inputWaveFile, "Input wavelength calibration file");
	args.AddOptionalArgument("combineMethod", combineMethod, 1, "Method for combining images: 1 = Median, 2 = Mean, 3 = Weighted mean -- WARNING: NOT CURRENTLY IMPLEMENTED"); //NOT IMPLEMENTED YET
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
	args.AddOptionalArgument("plotfilename", plotfilename, "", "Output plot eps file name");
	args.AddOptionalArgument("spectrumDataFilename", spectrumDataFilename, "", "Output spectrum data file name");
	args.AddOptionalArgument("outputDataFilename", outputDataFilename, "", "Output data file name");
	args.AddOptionalArgument("scriptfilename", scriptfilename, "", "Output gnuplot script file name");
	args.AddSwitch("interactive", interactive, "For interactive plots");
	
	try {
        args.Parse(argc, argv);
        
        string inputfcal[MAXNUMBEROFFLUXCALIBRATIONFILES];
		unsigned inputFcalIndex = 0;
		SplitStringIntoArray(listofinputfcal, inputfcal, inputFcalIndex, MAXNUMBEROFFLUXCALIBRATIONFILES); // Split list of images into array
        
		// we need at least one input opera flux calibration file..
        if(inputFcalIndex == 0 && inputconstant == 0.0) {
			throw operaException("operaMasterFluxCalibration: no input fcal or constant. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        } else if (inputFcalIndex > MAXNUMBEROFFLUXCALIBRATIONFILES) {
            throw operaException("operaMasterFluxCalibration: number of input fcal exceeds limit. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
		// we need an output calibration file name...
		if (outputfcal.empty()) {
			throw operaException("operaMasterFluxCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input reference spectrum...
		if (inputReferenceSpectrum.empty()) {
			throw operaException("operaMasterFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (args.verbose) {
            for(unsigned index=0; index<inputFcalIndex; index++) {
                cout << "operaMasterFluxCalibration: input flux calibration inputfcal["<<index<<"] = " << inputfcal[index] << endl;
            }
            cout << "operaMasterFluxCalibration: inputconstant = " << inputconstant << endl;
			cout << "operaMasterFluxCalibration: output Flux Calibration .fcal file = " << outputfcal << endl;
			cout << "operaMasterFluxCalibration: inputReferenceSpectrum = " << inputReferenceSpectrum << endl;
			cout << "operaMasterFluxCalibration: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaMasterFluxCalibration: combineMethod = " << combineMethod << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaMasterFluxCalibration: ordernumber = " << ordernumber << endl;
            if(args.plot) {
                cout << "operaMasterFluxCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaMasterFluxCalibration: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaMasterFluxCalibration: outputDataFilename = " << spectrumDataFilename << endl;
                cout << "operaMasterFluxCalibration: scriptfilename = " << scriptfilename << endl;
			    cout << "operaMasterFluxCalibration: interactive = " << (interactive ? "YES" : "NO") << endl;
            }
		}
        
        ofstream fspecdata;
        ofstream foutdata;
        if (!spectrumDataFilename.empty()) fspecdata.open(spectrumDataFilename.c_str());
        if (!outputDataFilename.empty()) foutdata.open(outputDataFilename.c_str());
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputReferenceSpectrum);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile);
        
        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
		if (args.verbose) cout << "operaMasterFluxCalibration: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        unsigned NumberofBeams = 0;
        double wavelengthForNormalization = 0;
        
        spectralOrders.setWavelengthsFromCalibration(minorder, maxorder);
        
        // initialize output vectors
        for (int order=minorder; order<=maxorder; order++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
			if (spectralOrder->gethasSpectralElements()) {
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                unsigned nElements = SpectralElements->getnSpectralElements();
                
                if(nElements==0) {
                    cerr << "operaMasterFluxCalibration: skipping order " << order << " -> reference spectrum has nElements=0" << endl;
                    continue;
                }
                
                if(NumberofBeams == 0) {
                    NumberofBeams = spectralOrder->getnumberOfBeams();
                }
                
                // create output energy distribution elements for output fcals
                spectralOrder->getSpectralEnergyDistribution()->setCalibrationWavelength(SpectralElements->getWavelength());
                for(unsigned b = 0; b < spectralOrder->MainAndBeamCount(); b++) {
					operaSpectralEnergyDistribution& spectralEnergyDistribution = spectralOrder->MainAndBeamSED(b);
					operaFluxVector fcal(nElements);
					if(inputFcalIndex == 0) fcal = inputconstant;
					spectralEnergyDistribution.setFluxCalibration(fcal);
					spectralEnergyDistribution.setThroughput(fcal);
					spectralEnergyDistribution.setHasFluxCalibration(true);
					spectralEnergyDistribution.setHasInstrumentThroughput(true);
				}
				
                spectralOrder->sethasSpectralEnergyDistribution(true);
            }
        }
        
        for(unsigned index=0; index<inputFcalIndex; index++) {
            // read input flux calibration
            operaSpectralOrderVector inputSpectralOrders;
            operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputfcal[index]);
            
            // get wavelengthForNormalization
            if(index==0) {
                for (int order=minorder; order<=maxorder; order++) {
                    operaSpectralOrder *inputSpectralOrder = inputSpectralOrders.GetSpectralOrder(order);
                    if (inputSpectralOrder->gethasSpectralEnergyDistribution()) {
                        operaSpectralEnergyDistribution *inputSpectralEnergyDistribution = inputSpectralOrder->getSpectralEnergyDistribution();
                        if(inputSpectralEnergyDistribution->getwavelengthForNormalization() != 0) {
                            wavelengthForNormalization = inputSpectralEnergyDistribution->getwavelengthForNormalization();
                            break;
                        }
                    }
                }
            }
            
            for (int order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                
                if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralEnergyDistribution()) {
                    
                    if (args.debug) cout << "operaMasterFluxCalibration: processing fcal index=" << index << " order=" << order << endl;
                    
                    operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                    operaSpectralEnergyDistribution *SpectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                    
                    SpectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
                    
                    unsigned nElements = SpectralElements->getnSpectralElements();
                    
                    // DT May 20 2014 you can't do spline fitting without some reasobale number of elements...
					if(nElements < MINELEMENTS) {
                        cerr << "operaMasterFluxCalibration: skipping order " << order << " -> reference spectrum has nElements < MINELEMENTS" << endl;
                        continue;
                    }
                    
                    const operaVector& referenceWavelength = SpectralElements->getWavelength();
                    
                    // grab order of input fcal
                    operaSpectralOrder *inputSpectralOrder = inputSpectralOrders.GetSpectralOrder(order);

                    if (inputSpectralOrder->gethasSpectralEnergyDistribution()) {
                        const operaVector& inputWavelength = inputSpectralOrder->getSpectralEnergyDistribution()->getCalibrationWavelength();
                        unsigned inputnElements = inputWavelength.size();
                        if(args.verbose) cout << "operaMasterFluxCalibration: processing fcal index=" << index << " order=" << order << " refnElements=" << nElements << " inputnElements=" << inputnElements << endl;
                        
                        bool interpolate = false;
                        if(inputnElements != nElements) {
                            interpolate = true;
                            if(args.verbose) cout << "operaMasterFluxCalibration: using interpolation since inputnElements=" << inputnElements << " != nElements=" << nElements << endl;
                        }
                        
                        // collect data from input fcal
                        if(!interpolate) {
							double maxdist = Max(Abs(inputWavelength - referenceWavelength));
							if(maxdist > WAVELENGTH_PRECISION) {
								interpolate = true;
								if(args.verbose) cout << "operaMasterFluxCalibration: using interpolation since |inputWavelength - referenceWavelength| = " << maxdist << " > " << WAVELENGTH_PRECISION << endl;
							}
						}
						
						// below is for plotting
						if(fspecdata.is_open()) {
                            for(unsigned indexElem=0; indexElem<inputnElements; indexElem++) {
                                fspecdata << index << ' ' << order << ' ' << indexElem << ' ';
                                for(unsigned b = 0; b < spectralOrder->MainAndBeamCount(); b++) {
                                    operaSpectralEnergyDistribution& inputSED = inputSpectralOrder->MainAndBeamSED(b);
                                    if(b == 0) fspecdata << inputSED.getCalibrationWavelength()[indexElem] << ' ';
                                    else fspecdata << b-1 << ' ';
                                    fspecdata << inputSED.getFluxCalibration().getflux(indexElem) << ' '
                                    << inputSED.getFluxCalibration().getvariance(indexElem) << ' '
                                    << inputSED.getThroughput().getflux(indexElem) << ' '
                                    << inputSED.getThroughput().getvariance(indexElem) << ' ';
                                }
                                fspecdata << endl;
                            }
                            fspecdata << endl;
                        }
                        
                        if(interpolate && args.verbose) cout << "operaMasterFluxCalibration: Starting interpolations" << endl;
                        
                        for(unsigned b = 0; b < spectralOrder->MainAndBeamCount(); b++) {
                            operaSpectralEnergyDistribution& inputSED = inputSpectralOrder->MainAndBeamSED(b);
							
                            operaFluxVector referenceInputFluxCal;
                            operaFluxVector referenceInputThroughput;
                            if(interpolate) {
								referenceInputFluxCal = fitSpectrum(inputWavelength, inputSED.getFluxCalibration().getflux(), referenceWavelength);
								referenceInputThroughput = fitSpectrum(inputWavelength, inputSED.getThroughput().getflux(), referenceWavelength);
							} else {
								referenceInputFluxCal = inputSED.getFluxCalibration();
								referenceInputThroughput = inputSED.getThroughput();
							}
							
							operaSpectralEnergyDistribution& outputSED = spectralOrder->MainAndBeamSED(b);
							outputSED.setFluxCalibration(outputSED.getFluxCalibration() + referenceInputFluxCal/(double)inputFcalIndex);
							outputSED.setThroughput(outputSED.getThroughput() + referenceInputThroughput/(double)inputFcalIndex);
						}
                        
                        // below is for plotting
                        if(foutdata.is_open() && index==inputFcalIndex-1) {
                            for(unsigned indexElem=0; indexElem<nElements; indexElem++) {
                                foutdata << index << ' ' << order << ' ' << indexElem << ' '
                                << SpectralElements->getwavelength(indexElem) << ' '
                                << SpectralElements->getFlux(indexElem) << ' '
                                << SpectralElements->getFluxVariance(indexElem) << ' '
                                << SpectralEnergyDistribution->getFluxCalibration().getflux(indexElem) << ' '
                                << SpectralEnergyDistribution->getFluxCalibration().getvariance(indexElem) << ' ';
                                for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                                    foutdata << beam << ' '
                                    << spectralOrder->getBeamSED(beam)->getFluxCalibration().getflux(indexElem) << ' '
                                    << spectralOrder->getBeamSED(beam)->getFluxCalibration().getvariance(indexElem) << ' '
                                    << spectralOrder->getBeamSED(beam)->getThroughput().getflux(indexElem) << ' '
                                    << spectralOrder->getBeamSED(beam)->getThroughput().getvariance(indexElem) << ' ';
                                }
                                foutdata << endl;
                            }
                            foutdata << endl;
                        }
                        
                    } else if (!inputSpectralOrder->gethasSpectralEnergyDistribution()) {
                        if (args.verbose) cout << "operaMasterFluxCalibration: Skipping order number: "<< order << "for index=" << index << " no input SpectralEnergyDistribution." << endl;
                    }
                } // if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralEnergyDistribution()) {
            } //for (int order=minorder; order<=maxorder; order++) {
        } //for(unsigned index=0; index<inputFcalIndex; index++) {
        
        /*
         * and write output
         */
        operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputfcal, Fcal);
        
        if (fspecdata.is_open() && foutdata.is_open()) {
            fspecdata.close();
            foutdata.close();
            
            if (!scriptfilename.empty()) {
                GenerateMasterFluxCalibrationPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(),outputDataFilename.c_str(), NumberofBeams, interactive);
            }
        }
	}
	catch (operaException e) {
		cerr << "operaMasterFluxCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMasterFluxCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

void GenerateMasterFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string outputDataFilename, unsigned NumberofBeams, bool display)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str());  // delete any existing file with the same name
    ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "\nset xlabel \"wavelength (nm)\"" << endl;
    fgnu << "set ylabel \"flux calibration\"" << endl;
    
    fgnu << "set pointsize 0.5" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nplot \"" << spectrumDataFilename << "\" u 4:5 w d" <<
        ",\"" << outputDataFilename << "\" u 4:7 w l lw 3" << endl;
        
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
        fgnu << "\nplot \"" << spectrumDataFilename << "\" u 4:5 w d" <<
        ",\"" << outputDataFilename << "\" u 4:7 w l lw 3" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}
