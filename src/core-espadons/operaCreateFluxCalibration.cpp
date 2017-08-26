/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCreateFluxCalibration
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
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/operaSpectralTools.h"			// void calculateUniformSample, getFluxAtWavelength
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

#define MAXFLUXREFERENCELENGTH 20000
#define MAXNUMBEROFREFWLRANGES 1000

/*! \file operaCreateFluxCalibration.cpp */

using namespace std;

operaArgumentHandler args;

/*
 * the reference spectrum
 */
static unsigned nPointsInReferenceSpectrum = 0;
static double referenceWavelength[MAXFLUXREFERENCELENGTH];
static double referenceIntensity[MAXFLUXREFERENCELENGTH];
static double referenceNormIntensity[MAXFLUXREFERENCELENGTH];
static double referenceVariance[MAXFLUXREFERENCELENGTH];  

unsigned getReferenceSpectrumRange(double wl0, double wlf, double **wl, double **flux, double **normflux, double **fluxvar);
unsigned readReferenceSpectrum(string reference_spectrum, double *referenceWavelength, double *referenceIntensity, double *referenceVariance);

void normalizeIntensityByMaximum(unsigned np, double *intensity, double *variance);
void normalizeIntensityByReferenceWavelength(unsigned np, double *intensity, double *wavelength, double *outputNormIntensity, double refWavelength);
double getReferenceFlux(unsigned np, double *intensity, double *wavelength, double refWavelength);

double operaArrayMaxValue_d(unsigned np, const double *xarray, const double *yarray, double *maxx);
unsigned getContinuumFromInputReferenceSpectrum(string inputWavelengthMaskForRefContinuum, float *refContinuumwl,float *refContinuumflux,float *refContinuumNormflux);
unsigned getReferenceSpectrumRange(unsigned nRefContinuum,double *refContinuumwl,double *refContinuumflux,double *refContinuumNormflux,double wl0,double wlf, double **wl, double **flux, double **normflux);

void GenerateCreateFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string continuumDataFilename, unsigned NumberofBeams, bool display);

/*! 
 * operaCreateFluxCalibration
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
	string inputUncalibratedSpectrum;
	string inputCalibratedSpectrum;
	string inputWavelengthMaskForUncalContinuum;
    string inputWavelengthMaskForRefContinuum;
    string inputFlatFluxCalibration;
	string inputWaveFile;
	string inputaper;
	string outputFluxCalibrationFile;
	
    double wavelengthForNormalization = 548;
    double exposureTime = 1;   // in seconds
    unsigned numberOfPointsInUniformSample = 200;
    unsigned numberOfPointsInUniformRefSample = 70;
    unsigned binsize = 100;
    
	int ordernumber = NOTPROVIDED;
    int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;    
	string plotfilename;	
	string spectrumDataFilename;
	string continuumDataFilename;	
	string scriptfilename;	
	bool interactive = false;
	
	args.AddRequiredArgument("inputUncalibratedSpectrum", inputUncalibratedSpectrum, "Spectrophotometric standard extracted uncalibrated spectrum");
	args.AddRequiredArgument("inputCalibratedSpectrum", inputCalibratedSpectrum, "Spectrophotometric standard template calibrated spectrum");
	args.AddRequiredArgument("inputWavelengthMaskForUncalContinuum", inputWavelengthMaskForUncalContinuum, "");  //???
	args.AddRequiredArgument("inputWavelengthMaskForRefContinuum", inputWavelengthMaskForRefContinuum, "");      //???
	args.AddRequiredArgument("inputFlatFluxCalibration", inputFlatFluxCalibration, "Input flat flux calibration file");
	args.AddRequiredArgument("inputWaveFile", inputWaveFile, "Input wavelength calibration file");
	args.AddRequiredArgument("inputApertureFile", inputaper, "Input aperture calibration file");
	args.AddRequiredArgument("outputFluxCalibrationFile", outputFluxCalibrationFile, "Output flux calibration conversion file");
	
	args.AddRequiredArgument("wavelengthForNormalization", wavelengthForNormalization, "Wavelength (nm) for normalization of reference spectrum");
	args.AddRequiredArgument("exposureTime", exposureTime, "Exposure time of input uncalibrated spectrum");
	args.AddRequiredArgument("numberOfPointsInUniformSample", numberOfPointsInUniformSample, "Define lowest order to consider in the fit across orders");
	args.AddRequiredArgument("numberOfPointsInUniformRefSample", numberOfPointsInUniformRefSample, "Define highest order to consider in the fit across orders");
	args.AddRequiredArgument("binsize", binsize, "Number of points to bin for continuum estimate");
	
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
	args.AddOptionalArgument("plotfilename", plotfilename, "", "Output plot eps file name");
	args.AddOptionalArgument("spectrumDataFilename", spectrumDataFilename, "", "Output spectrum data file name");
	args.AddOptionalArgument("continuumDataFilename", continuumDataFilename, "", "Output continuum data file name");
	args.AddOptionalArgument("scriptfilename", scriptfilename, "", "Output gnuplot script file name");
	args.AddSwitch("interactive", interactive, "For interactive plots");
	
	try {
		args.Parse(argc, argv);
		
		// we need an input uncalibrated spectrum...
		if (inputUncalibratedSpectrum.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an input template spectrum...
		if (inputCalibratedSpectrum.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an input flat flux calibration spectrum...
		if (inputFlatFluxCalibration.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output...
		if (outputFluxCalibrationFile.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength mask file for ref...
		if (inputWavelengthMaskForRefContinuum.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength mask file for uncal...
		if (inputWavelengthMaskForUncalContinuum.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a aperture file...
		if (inputaper.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
		if (args.verbose) {
			cout << "operaCreateFluxCalibration: input uncalibrated spectrum file = " << inputUncalibratedSpectrum << endl; 
			cout << "operaCreateFluxCalibration: input calibrated spectrum file = " << inputCalibratedSpectrum << endl;
			cout << "operaCreateFluxCalibration: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
			cout << "operaCreateFluxCalibration: inputWavelengthMaskForRefContinuum = " << inputWavelengthMaskForRefContinuum << endl;
			cout << "operaCreateFluxCalibration: inputWavelengthMaskForUncalContinuum = " << inputWavelengthMaskForUncalContinuum << endl;
			cout << "operaCreateFluxCalibration: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaCreateFluxCalibration: output flux calibration file = " << outputFluxCalibrationFile << endl;
			cout << "operaCreateFluxCalibration: input aperture file = " << inputaper << endl;                        
            cout << "operaCreateFluxCalibration: wavelengthForNormalization= " << wavelengthForNormalization << " nm" << endl;
            cout << "operaCreateFluxCalibration: exposure time= " << exposureTime << " seconds" << endl;
			cout << "operaCreateFluxCalibration: numberOfPointsInUniformSample = " << numberOfPointsInUniformSample << endl;
			cout << "operaCreateFluxCalibration: numberOfPointsInUniformRefSample = " << numberOfPointsInUniformRefSample << endl;
			if(ordernumber != NOTPROVIDED) cout << "operaCreateFluxCalibration: ordernumber = " << ordernumber << endl;            
            cout << "operaCreateFluxCalibration: binsize = " << binsize << endl;  
            if(args.plot) {
                cout << "operaCreateFluxCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaCreateFluxCalibration: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaCreateFluxCalibration: continuumDataFilename = " << continuumDataFilename << endl;
                cout << "operaCreateFluxCalibration: scriptfilename = " << scriptfilename << endl; 
				cout << "operaCreateFluxCalibration: interactive = " << (interactive ? "YES" : "NO") << endl; 
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
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputUncalibratedSpectrum);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputaper);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile);

        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
		if (args.verbose) cout << "operaCreateFluxCalibration: minorder ="<< minorder << " maxorder=" << maxorder << endl;        

		/*
		 * Flux calibration reference file:
         *
		 * Read reference calibrated spectrum
		 *		lambda vs. intensity, intensityVariance (optional)
		 */
        nPointsInReferenceSpectrum = readReferenceSpectrum(inputCalibratedSpectrum, referenceWavelength, referenceIntensity, referenceVariance);        
        normalizeIntensityByReferenceWavelength(nPointsInReferenceSpectrum,referenceIntensity,referenceWavelength,referenceNormIntensity,wavelengthForNormalization);
        
        //---------------------------------
        // Loop over orders to set maximum number of elements, set wavelength and the number of beams
        // --> maxNElements & NumberofBeams
        unsigned NumberofBeams = spectralOrders.getNumberOfBeams(minorder, maxorder);
        unsigned maxNElements = spectralOrders.getMaxNumberOfElementsInOrder(minorder, maxorder);
        spectralOrders.setWavelengthsFromCalibration(minorder, maxorder);
        if (args.verbose) cout << "operaCreateFluxCalibration: NumberofBeams = " << NumberofBeams << endl;

        if(NumberofBeams == 0) {
            throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }

        //---------------------------------
        // below one can define constants for flux calibration
        double spectralBinConstant = exposureTime;
        
        //---------------------------------
        // Correct flat-field
        if (!inputFlatFluxCalibration.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputFlatFluxCalibration);
            spectralOrders.correctFlatField(minorder, maxorder, false);
        }

        //---------------------------------
        // Calculate a clean sample of the continuum from the ref spectrum
        float *refContinuumwl = new float[MAXNUMBEROFREFWLRANGES];
        float *refContinuumflux = new float[MAXNUMBEROFREFWLRANGES];
        float *refContinuumNormflux = new float[MAXNUMBEROFREFWLRANGES];
        unsigned nRefContinuum = getContinuumFromInputReferenceSpectrum(inputWavelengthMaskForRefContinuum,refContinuumwl,refContinuumflux,refContinuumNormflux);

        double refFluxForNormalization = (double)getFluxAtWavelength(nRefContinuum,refContinuumwl,refContinuumflux,wavelengthForNormalization);
        
        operaSpectrum refContinuum;
        for(unsigned i = 0; i < nRefContinuum; i++) refContinuum.insert(refContinuumwl[i], refContinuumflux[i]);
        
        operaSpectrum uniformRef = calculateUniformSample(refContinuum, numberOfPointsInUniformRefSample);
        
        if(args.debug) {
            // original sample
            for(unsigned i=0;i<nRefContinuum;i++) {
                cout << refContinuumwl[i] << " "
                << refContinuumflux[i] << " "
                << refContinuumNormflux[i] << endl;
            }
            // uniform sample
            for (unsigned i=0; i<uniformRef.size(); i++) {
                cout << uniformRef.getwavelength(i) << " " << uniformRef.getflux(i) << endl;
            }
        }

        const double delta_wl = 1.0; // wavelength (in nm) range for stiching non-overlapping orders
		operaSpectrum uniform_spectrum = spectralOrders.calculateCleanUniformSampleOfContinuum(minorder, maxorder, binsize, delta_wl, inputWavelengthMaskForUncalContinuum, numberOfPointsInUniformSample, true);
        
        operaVector uncalFluxForNormalization = getFluxesAtWavelength(uniform_spectrum, wavelengthForNormalization);
        
        //---------------------------------
        // Use clean sample of the continuum from both the reference and uncalibrated spectrum
        // to calculate flux calibration and instrument throughput for all orders
        // Results are saved in the SED classes which will be written out to fcal
        
        for (int order=minorder; order<=maxorder; order++) {
            
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            if (spectralOrder->gethasWavelength() &&
                spectralOrder->gethasSpectralElements() &&
                spectralOrder->gethasSpectralEnergyDistribution()) {
				            
				operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                spectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
                
                unsigned nElements = spectralElements->getnSpectralElements();

                if(!continuumDataFilename.empty()) {
                    for(unsigned i=0;i<uniform_spectrum.size();i++) *fcontinuumdata << uniform_spectrum.getwavelength(i) << ' ' << uniform_spectrum.getflux(i) << endl;
                }
                
                operaVector CalibratedModelFlux = fitSpectrum(uniformRef.wavelengthvector(), uniformRef.fluxvector(), spectralElements->getWavelength());
                spectralEnergyDistribution->setCalibratedFlux(CalibratedModelFlux/refFluxForNormalization);
                for(unsigned b=0;b < spectralOrder->MainAndBeamCount(); b++) {
					operaVector UncalibratedModelFlux = fitSpectrum(uniform_spectrum.wavelengthvector(), uniform_spectrum.fluxvector(b), spectralElements->getWavelength());
					spectralOrder->MainAndBeamSED(b).setFluxCalibration(CalibratedModelFlux/(UncalibratedModelFlux/spectralBinConstant));
                    spectralOrder->MainAndBeamSED(b).setThroughput((CalibratedModelFlux/UncalibratedModelFlux)*(uncalFluxForNormalization[b]/refFluxForNormalization));
				}
                
                spectralEnergyDistribution->setHasCalibratedFlux(true);
                spectralEnergyDistribution->setHasUncalibratedFlux(true);
                operaVector UncalibratedModelFluxAgain = fitSpectrum(uniform_spectrum.wavelengthvector(), uniform_spectrum.fluxvector(), spectralElements->getWavelength());
                if(!spectrumDataFilename.empty()) {
                    for(unsigned i=0; i<spectralElements->getnSpectralElements(); i++) {
                        *fspecdata << order << " " << i << " "
                        << spectralElements->getwavelength(i) << " "
                        << spectralElements->getFlux(i) << " "
                        << UncalibratedModelFluxAgain[i] << " "
                        << CalibratedModelFlux[i] << " "
                        << spectralEnergyDistribution->getUncalibratedFlux().getflux(i) << " "
                        << spectralEnergyDistribution->getCalibratedFlux().getflux(i) << " "
                        << spectralEnergyDistribution->getFluxCalibration().getflux(i) << " "
                        << spectralEnergyDistribution->getThroughput().getflux(i) << " ";
                        for(unsigned beam=0;beam<NumberofBeams;beam++) {
                            *fspecdata << spectralOrder->getBeamSED(beam)->getUncalibratedFlux().getflux(i) << " "
                                << spectralOrder->getBeamSED(beam)->getFluxCalibration().getflux(i) << " "
                                << spectralOrder->getBeamSED(beam)->getThroughput().getflux(i) << " ";
                        }
                        *fspecdata << endl;
                    }
                }
            } else {
                spectralOrder->sethasSpectralEnergyDistribution(FALSE);
            }
        }
        
		/*
		 * and write out fcal
		 */
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputFluxCalibrationFile, Fcal);
		
        if (fspecdata != NULL && fcontinuumdata != NULL) {
			fspecdata->close();
            fcontinuumdata->close();
            
            if (!scriptfilename.empty()) {
                GenerateCreateFluxCalibrationPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(),continuumDataFilename.c_str(), NumberofBeams, interactive);
            }
        }
        
	}
	catch (operaException e) {
		cerr << "operaCreateFluxCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCreateFluxCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

void GenerateCreateFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string continuumDataFilename, unsigned NumberofBeams, bool display)
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

/*
 * Read the the full reference spectrum
 */
unsigned readReferenceSpectrum(string reference_spectrum, double *referenceWavelength, double *referenceIntensity, double *referenceVariance) {
	ifstream astream;
	string dataline;
    
	double tmpwl = -1.0; 
	double tmpi = -1.0; 
	unsigned np = 0;
	
	astream.open(reference_spectrum.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl, &tmpi);
                    referenceWavelength[np] = tmpwl;
                    referenceIntensity[np] = tmpi;
                    referenceVariance[np] = tmpi;
                    np++;  
                }	// skip comments
            }
		} // while (astream.good())
		if (np > 0 && args.verbose) printf("          [Reference] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", np, referenceWavelength[0], referenceWavelength[np/2], referenceWavelength[np-1]);
		else if (args.verbose) printf("          [Reference] no points found in flux reference file.\n");
		astream.close();
	}	// if (astream.open())
	return np;
}

/*
 * get a subset of the reference spectrum for this order only, between wl0 and wlf
 */
unsigned getReferenceSpectrumRange(double wl0, double wlf, double **wl, double **flux, double **normflux, double **fluxvar) {
	unsigned firstline = 0;
	unsigned np = 0;
    
	for (np=0; np<nPointsInReferenceSpectrum; np++) {
		if (referenceWavelength[np] >= wl0) {
			if (firstline == 0 && np) {
                    *flux = &referenceIntensity[np-1];
                    *normflux = &referenceNormIntensity[np-1];
                    *wl = &referenceWavelength[np-1];
                    *fluxvar = &referenceVariance[np-1];
                    firstline = np-1;
			}
			if (referenceWavelength[np] > wlf) {
                np++;
				break;
            }
		}
	}
	if (np == nPointsInReferenceSpectrum)
		np--;
	if (np > firstline) {
		return (np-firstline);
	} else {
		return 0;
	}
}

void normalizeIntensityByReferenceWavelength(unsigned np, double *intensity, double *wavelength, double *outputNormIntensity, double refWavelength) {
    double referenceFlux = getReferenceFlux(np,intensity,wavelength,refWavelength);
    
	for(unsigned i=0;i<np;i++) {
        outputNormIntensity[i] = intensity[i]/referenceFlux;
	}
}

double getReferenceFlux(unsigned np, double *intensity, double *wavelength, double refWavelength) {
    
    float *wavelengthData_f = new float[np];
    float *fluxData_f = new float[np];
    
    for(unsigned i=0;i<np;i++) {
        wavelengthData_f[i] = (float)wavelength[i];
        fluxData_f[i] =  (float)intensity[i];
    }
    
    unsigned nElements = 1;
    
    float *referenceFlux = new float[nElements];
    float *referencewl = new float[nElements];
    
    for(unsigned i=0;i<nElements;i++) {
        referencewl[i] = (float)refWavelength;
    }
    operaFitSpline(np,wavelengthData_f,fluxData_f,nElements,referencewl,referenceFlux);
    
    double outputflux = (double)(referenceFlux[0]);
    
    delete[] wavelengthData_f;
    delete[] fluxData_f;
    delete[] referenceFlux;
    delete[] referencewl;
    
	return outputflux;
}

void normalizeIntensityByMaximum(unsigned np, double *intensity, double *variance) {
	double maxIntensity = -3.4e+38;
	
	for(unsigned i=0;i<np;i++) {
		if(intensity[i] > maxIntensity)
			maxIntensity = intensity[i];
	}
	
	for(unsigned i=0;i<np;i++) {
        intensity[i] /= maxIntensity;
        variance[i] /= maxIntensity;
	}
}


double operaArrayMaxValue_d(unsigned np, const double *xarray, const double *yarray, double *maxx) {
	double ymax = -3.4e+38;
	double xmax = 0;
    
	while (np--) {
		if(*yarray > ymax) {
			ymax = *yarray;
            xmax = *xarray;
        }
		yarray++;
        xarray++;
	}
    *maxx = xmax;
	return ymax;
}


unsigned getContinuumFromInputReferenceSpectrum(string inputWavelengthMaskForRefContinuum, float *refContinuumwl,float *refContinuumflux,float *refContinuumNormflux) {
    
    double *wl0_vector = new double[MAXNUMBEROFREFWLRANGES];
    double *wlf_vector = new double[MAXNUMBEROFREFWLRANGES];
    
    unsigned nRangesInWLMask = readContinuumWavelengthMask(inputWavelengthMaskForRefContinuum,wl0_vector,wlf_vector);
    
    unsigned nTotalPoints = 0;
    for(unsigned k=0;k<nRangesInWLMask; k++){
        
        double *refwl = NULL, *refflux = NULL, *refnormflux = NULL, *reffluxvar = NULL;
        unsigned nPointsInReference = getReferenceSpectrumRange(wl0_vector[k],wlf_vector[k],&refwl,&refflux,&refnormflux,&reffluxvar);
        
        double ref_wl=0;
        double ref_maxFlux = 0;
        double ref_maxNormFlux = 0;
        
        ref_maxFlux = operaArrayMaxValue_d(nPointsInReference,refwl,refflux,&ref_wl);
        ref_maxNormFlux = operaArrayMaxValue_d(nPointsInReference,refwl,refnormflux, &ref_wl);
        
        if(args.debug) {
            cout << k << " "
            << nPointsInReference << " "
            << wl0_vector[k] << " "
            << wlf_vector[k] << " "
            << ref_wl << " "
            << ref_maxFlux << " "
            << ref_maxNormFlux << endl;
        }
        
        if(ref_wl && ref_maxFlux && ref_maxNormFlux) {
            refContinuumwl[nTotalPoints] = (float)ref_wl;
            refContinuumflux[nTotalPoints] = (float)ref_maxFlux;
            refContinuumNormflux[nTotalPoints] = (float)ref_maxNormFlux;
            nTotalPoints++;
        }
    }
    
    delete[] wl0_vector;
    delete[] wlf_vector;
    
    return nTotalPoints;
}
