/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCreateFlatResponse
 Version: 1.0
 Description: Flat Response Flux Calibration with Standard or Moon spectrum
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Feb/2015
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

/*! \file operaCreateFlatResponse.cpp */

using namespace std;

operaArgumentHandler args;

operaSpectrum readReferenceSpectrum(string reference_spectrum);
double getIntensityAtReferenceWavelength(const operaSpectrum& spectrum, double refWavelength);
operaSpectrum getContinuumFromInputReferenceSpectrum(string inputWavelengthMaskForRefContinuum, const operaSpectrum& referenceSpectrum, double referenceFluxForNormalization);

/*! 
 * operaCreateFlatResponse
 * \author Eder Martioli
 * \brief Flat Response Flux Calibration with Standard or Moon spectrum
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
	const double DELTA_WL = 1.0; // wavelength (in nm) range for stiching non-overlapping orders
    
	string inputUncalibratedSpectrum;
    string inputSpectrumFITSImage;
	string inputCalibratedSpectrum;
    string inputFlatFluxCalibration;
	string inputWaveFile;
    string inputWavelengthMaskForRefContinuum;
    string inputWavelengthMaskForUncalContinuum;
	string outputFlatResponseFile;
    double wavelengthForNormalization = 548;
    int ordernumber = NOTPROVIDED;
    int minorder = 22;
    int maxorder = 62;    
	unsigned numberOfPointsInUniformSample = 200;
    unsigned numberOfPointsInUniformRefSample = 70;
    unsigned binsize = 100;
    bool outputFITS = false;
    
    args.AddRequiredArgument("inputUncalibratedSpectrum", inputUncalibratedSpectrum, "Spectrophotometric standard extracted uncalibrated spectrum");
    args.AddRequiredArgument("inputSpectrumFITSImage", inputSpectrumFITSImage, "Raw FITS image for input uncalibrated spectrum");
	args.AddRequiredArgument("inputCalibratedSpectrum", inputCalibratedSpectrum, "Spectrophotometric standard template calibrated spectrum");
	args.AddRequiredArgument("inputFlatFluxCalibration", inputFlatFluxCalibration, "Flat flux calibration data file (.fcal.gz)");
	args.AddRequiredArgument("inputWaveFile", inputWaveFile, "Input wavelength calibration file");
	args.AddRequiredArgument("inputWavelengthMaskForRefContinuum", inputWavelengthMaskForRefContinuum, "Wavelength mask to detect continuum in reference spectrum");
	args.AddRequiredArgument("inputWavelengthMaskForUncalContinuum", inputWavelengthMaskForUncalContinuum, "Wavelength mask to detect continuum in uncalibrated spectrum");
	args.AddRequiredArgument("outputFlatResponseFile", outputFlatResponseFile, "Output flat response data file");
	args.AddOptionalArgument("wavelengthForNormalization", wavelengthForNormalization, 548, "Wavelength (nm) for normalization of reference spectrum");
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
	args.AddOptionalArgument("numberOfPointsInUniformSample", numberOfPointsInUniformSample, 200, "Define number of points in output data file");
	args.AddOptionalArgument("numberOfPointsInUniformRefSample", numberOfPointsInUniformRefSample, 70, "Define number of poins in reference sample");
	args.AddOptionalArgument("binsize", binsize, 100, "Number of points to bin for continuum estimate");
    args.AddSwitch("outputFITS", outputFITS, "output data as FITS file? otherwise output is in ASCII LE format");
    
	try {
		args.Parse(argc, argv);
		
		// we need an input uncalibrated spectrum...
		if (inputUncalibratedSpectrum.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an input template spectrum...
		if (inputCalibratedSpectrum.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
        if(outputFITS) {
            // if output is FITS then we need the input FITS image of the spectrum to grab header info...
            if (inputSpectrumFITSImage.empty()) {
                throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
            }
        }
		// we need an input flat flux calibration spectrum...
		if (inputFlatFluxCalibration.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength mask file for ref...
		if (inputWavelengthMaskForRefContinuum.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength mask file for uncal...
		if (inputWavelengthMaskForUncalContinuum.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output...
		if (outputFlatResponseFile.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
        }
        
		if (args.verbose) {
			cout << "operaCreateFlatResponse: input uncalibrated spectrum file = " << inputUncalibratedSpectrum << endl;
			cout << "operaCreateFlatResponse: input calibrated spectrum file = " << inputCalibratedSpectrum << endl;
			cout << "operaCreateFlatResponse: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
			cout << "operaCreateFlatResponse: inputWavelengthMaskForRefContinuum = " << inputWavelengthMaskForRefContinuum << endl;
			cout << "operaCreateFlatResponse: inputWavelengthMaskForUncalContinuum = " << inputWavelengthMaskForUncalContinuum << endl;
			cout << "operaCreateFlatResponse: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaCreateFlatResponse: output flux calibration file = " << outputFlatResponseFile << endl;
            cout << "operaCreateFlatResponse: wavelengthForNormalization= " << wavelengthForNormalization << " nm" << endl;
			cout << "operaCreateFlatResponse: numberOfPointsInUniformSample = " << numberOfPointsInUniformSample << endl;
			cout << "operaCreateFlatResponse: numberOfPointsInUniformRefSample = " << numberOfPointsInUniformRefSample << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaCreateFlatResponse: ordernumber = " << ordernumber << endl;
            cout << "operaCreateFlatResponse: binsize = " << binsize << endl;
		}
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputUncalibratedSpectrum);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile);

        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        if (args.verbose) cout << "operaCreateFlatResponse: minorder ="<< minorder << " maxorder=" << maxorder << endl;        

		/*
		 * Flux calibration reference file:
         *
		 * Read reference calibrated spectrum
		 *		lambda vs. intensity, intensityVariance (optional)
		 */
        operaSpectrum referenceSpectrum = readReferenceSpectrum(inputCalibratedSpectrum);
        double referenceFluxForNormalization = getIntensityAtReferenceWavelength(referenceSpectrum, wavelengthForNormalization);
        
        //---------------------------------
        // Loop over orders to set maximum number of elements, set wavelength and the number of beams
        // --> maxNElements & NumberofBeams
        unsigned NumberofBeams = spectralOrders.getNumberOfBeams(minorder, maxorder);
        spectralOrders.setWavelengthsFromCalibration(minorder, maxorder);
        if (args.verbose) cout << "operaCreateFlatResponse: NumberofBeams = " << NumberofBeams << endl;
        if(NumberofBeams == 0) throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);

        
        //---------------------------------
        // Correct flat-field
        if (!inputFlatFluxCalibration.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputFlatFluxCalibration);
            spectralOrders.correctFlatField(minorder, maxorder, false);
        }
        
        //---------------------------------
        // Calculate a clean sample of the continuum from the ref spectrum
        operaSpectrum refContinuum = getContinuumFromInputReferenceSpectrum(inputWavelengthMaskForRefContinuum, referenceSpectrum, referenceFluxForNormalization);

        operaSpectrum uniformRef = calculateUniformSample(refContinuum, numberOfPointsInUniformRefSample);

        operaSpectrum uniform = spectralOrders.calculateCleanUniformSampleOfContinuum(minorder, maxorder, binsize, DELTA_WL, inputWavelengthMaskForUncalContinuum, numberOfPointsInUniformSample, true);

        const operaVector& uniform_wl = uniform.wavelengthvector();
        operaVector calibratedModelFlux = fitSpectrum(uniformRef.wavelengthvector(), uniformRef.fluxvector(), uniform_wl);

        operaVector flatResp = uniform.fluxvector() / calibratedModelFlux;

        double flatRespForNormalization = getFluxAtWavelength(uniform_wl, flatResp, wavelengthForNormalization);

        if(outputFITS) {
            /*
             * and write out flatresponse (FITS *.fits.gz)
             */
            operaFITSImage inImage(inputSpectrumFITSImage, tfloat, READONLY);
            operaFITSImage outFlatResp(outputFlatResponseFile, numberOfPointsInUniformSample, 2, tfloat);
            outFlatResp.operaFITSImageCopyHeader(&inImage);
            for(unsigned i=0;i<numberOfPointsInUniformSample;i++) {
                flatResp[i] /= flatRespForNormalization;
                outFlatResp.setpixel(float(uniform_wl[i]),i,0);
                outFlatResp.setpixel(float(flatResp[i]),i,1);
            }
            
            outFlatResp.operaFITSAddComment("Created by the OPERA 1.0");
            outFlatResp.operaFITSAddComment("Flat response calibration file");
            outFlatResp.operaFITSAddComment("1st row: Wavelength [nm]");
            outFlatResp.operaFITSAddComment("2nd row: Relative flat response [abitrary units]");
            
            outFlatResp.operaFITSImageSave();
            outFlatResp.operaFITSImageClose();
            
        } else {
            /*
             * and write out flatresponse (LE *.s)
             */
            ofstream frespoutput(outputFlatResponseFile.c_str());
            frespoutput << "***" << endl;
            frespoutput << numberOfPointsInUniformSample << " 1" << endl;
            for(unsigned i=0;i<numberOfPointsInUniformSample;i++) {
                flatResp[i] /= flatRespForNormalization;
                frespoutput << uniform_wl[i] << ' ' << flatResp[i] << endl;
            }
            frespoutput.close();
        }
        
	}
	catch (operaException e) {
		cerr << "operaCreateFlatResponse: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCreateFlatResponse: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
} 

operaSpectrum readReferenceSpectrum(string reference_spectrum) {
	operaSpectrum referenceSpectrum;
	ifstream astream(reference_spectrum.c_str());
	if (astream.is_open()) {
		string dataline;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') {
				double tmpwl, tmpi;
				sscanf(dataline.c_str(), "%lf %lf", &tmpwl, &tmpi);
                referenceSpectrum.insert(tmpwl, tmpi, tmpi);
            }
		}
        
        if (args.verbose) {
			if (referenceSpectrum.size() > 0) printf("          [Reference] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", referenceSpectrum.size(), referenceSpectrum.firstwl(), referenceSpectrum.midwl(), referenceSpectrum.lastwl());
			else printf("          [Reference] no points found in flux reference file.\n");
		}
		astream.close();
	}
	return referenceSpectrum;
}

double getIntensityAtReferenceWavelength(const operaSpectrum& spectrum, double refWavelength) {
    double refFlux;
    operaFitSplineDouble(spectrum.size(), spectrum.wavelength_ptr(), spectrum.flux_ptr(), 1, &refWavelength, &refFlux);
	return refFlux;
}

operaSpectrum getContinuumFromInputReferenceSpectrum(string inputWavelengthMaskForRefContinuum, const operaSpectrum& referenceSpectrum, double referenceFluxForNormalization) {
    operaWavelengthRanges wlranges = readContinuumWavelengthMask(inputWavelengthMaskForRefContinuum);
    operaSpectrum refContinuum;
    operaVector refContinuumNormflux; //Not used.
    for(unsigned k=0;k<wlranges.size(); k++){
        operaSpectrum subSpectrum = getSpectrumWithinRange(wlranges.getrange(k), referenceSpectrum);
        if(subSpectrum.size() == 0) continue;
        
        unsigned maxindex = MaxIndex(subSpectrum.fluxvector());
        double ref_wl = subSpectrum.getwavelength(maxindex);
        double ref_maxFlux = subSpectrum.getflux(maxindex);
        double ref_maxNormFlux = ref_maxFlux / referenceFluxForNormalization;
        
        if(args.debug) {
            cout << k << " " << subSpectrum.size() << " " << wlranges.wl0(k) << " " << wlranges.wlf(k) << " " << ref_wl << " " << ref_maxFlux << " " << ref_maxNormFlux << endl;
        }
        
        if(ref_wl && ref_maxFlux && ref_maxNormFlux) {
            refContinuum.insert(ref_wl, ref_maxFlux);
            refContinuumNormflux.insert(ref_maxNormFlux);
        }
    }
    return refContinuum;
}
