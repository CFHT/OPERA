/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaTelluricWavelengthCorrection
 Version: 1.0
 Description: Apply wavelength correction based on telluric lines
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2015
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
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaSpectralTools.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFit.h"							// for operaFitSplineDouble
#include "libraries/gzstream.h"							// for gzstream - read compressed reference spectra
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"
#include "core-espadons/operaTelluricWavelengthCorrection.h"

/*! \file operaTelluricWavelengthCorrection.cpp */

using namespace std;

operaArgumentHandler args;

/*! 
 * operaTelluricWavelengthCorrection
 * \author Eder Martioli
 * \brief Calculate and apply wavelength correction based on telluric lines.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \note --wave=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{    
	string inputWaveFile;
    string inputObjectSpectrum;
    string inputFlatFluxCalibration;
	string outputWaveFile;
    string telluric_lines; // HITRAN Library
    string inputWavelengthMaskForTelluric;
    bool hitranformat;

	int ordernumber = NOTPROVIDED;
    int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;
    bool StarPlusSky = false;
    bool starplusskyInvertSkyFiber = false;
	double skyOverStarFiberAreaRatio = 1.0;
    
    // The parameters below we don't know yet whether they would be useful if used as input
    double spectralResolution = 80000;
    double radialVelocityRange = 10;
    double radialVelocityStep = 0.05;
    double XCorrelationThreshold = 0.05;
    unsigned normalizationBinsize = 110;
    bool useFitToFindMaximum = false;
    
	double LocalMaxFilterWidth=4.0;
	double MinPeakDepth=3;
	double DetectionThreshold=0.05;
	double nsigclip=3.0;
	bool emissionSpectrum = false;
	unsigned minNumberOfMatchedLines = 10;
	double duplicateLineThreshold = 0.001;

    unsigned RVCorrectionMethod = 1; // 1. line matching; 2. x-correlation;
    
    string rvcorrsplotfilename;
    string specplotfilename;
    string rvcorrscriptfilename;
    string specscriptfilename;
    string rvcorrdatafilename;
    string specdatafilename;
    string rvcorrfitdatafilename;
    string linesdatafilename;
    
    args.AddRequiredArgument("inputWaveFile", inputWaveFile, "input wavelength calibration file (.wcal)");
	args.AddRequiredArgument("inputObjectSpectrum", inputObjectSpectrum, "input object spectrum file (.e or .p)");
    args.AddOptionalArgument("inputFlatFluxCalibration", inputFlatFluxCalibration, "", "flat field spectrum ff_");
    args.AddRequiredArgument("outputWaveFile", outputWaveFile, "output telluric wavelength calibration file (.tell)");
    args.AddRequiredArgument("telluric_lines", telluric_lines, "atlas of telluric lines");
    args.AddRequiredArgument("inputWavelengthMaskForTelluric", inputWavelengthMaskForTelluric, "telluric wavelength mask");
    args.AddOptionalArgument("hitranformat", hitranformat, true, "telluric line atlas is in HITRAN format");
    
    args.AddRequiredArgument("spectralResolution", spectralResolution, "input spectral resolution (wl/dwl) as reference for line detection");
    args.AddRequiredArgument("radialVelocityRange", radialVelocityRange, "radial Velocity Range (in km/s) to scan for first order correction");
    args.AddRequiredArgument("radialVelocityStep", radialVelocityStep, "radial Velocity Step step (in km/s) to scan for first order correction");
    args.AddRequiredArgument("normalizationBinsize", normalizationBinsize, "binsize to normalize input object spectrum");
    
    args.AddRequiredArgument("RVCorrectionMethod", RVCorrectionMethod, "Method for measuring RV correction: 1 = Line matching (default), 2 = Spectral cross-correlation");
    
    args.AddOptionalArgument("LocalMaxFilterWidth", LocalMaxFilterWidth, 0.0, "Method 1: To set a window filter to guarantee a line is not detected twice, in units of line width");
    args.AddOptionalArgument("DetectionThreshold", DetectionThreshold, 0.0, "Method 1: Threshold to regulate the sensitivity of line detection, between 0 and 1");
    args.AddOptionalArgument("MinPeakDepth", MinPeakDepth, 0.0, "Method 1: Limit that regulates the sensitity of line detection, in units of noise");
    args.AddOptionalArgument("nsigclip", nsigclip, 0.0, "Method 1: Threshold to filter detected lines by line width compared to median line width, in units of line width error");
    args.AddSwitch("emissionSpectrum", emissionSpectrum, "Method 1: Using emission spectrum for object");
    args.AddOptionalArgument("minNumberOfMatchedLines", minNumberOfMatchedLines, 0, "Method 1: Arbitrary threshold to avoid small number statistics");
    args.AddOptionalArgument("duplicateLineThreshold", duplicateLineThreshold, 0, "Method 1: Minimum distance between two object lines");
    
    args.AddOptionalArgument("XCorrelationThreshold", XCorrelationThreshold, 0.0, "Method 2: X-correlation lower threshold to consider a match between telluric and object spectra");
    args.AddSwitch("useFitToFindMaximum", useFitToFindMaximum, "Method 2: use gaussian fit to find RV correction");
    
    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    args.AddSwitch("StarPlusSky", StarPlusSky, "star plus sky mode");
    args.AddSwitch("starplusskyInvertSkyFiber", starplusskyInvertSkyFiber, "Invert sky fiber (default is beam[0]=star, beam[1]=sky)");
	args.AddOptionalArgument("skyOverStarFiberAreaRatio", skyOverStarFiberAreaRatio, 1.0, "Area of the sky fiber over the area of the star fiber.");
    
    args.AddOptionalArgument("rvcorrsplotfilename", rvcorrsplotfilename, "", "Output RV correction plot eps file name");
    args.AddOptionalArgument("specplotfilename", specplotfilename, "", "Output spectrum plot eps file name");
    args.AddOptionalArgument("rvcorrscriptfilename", rvcorrscriptfilename, "", "Output cross-correlation gnuplot script file name");
    args.AddOptionalArgument("specscriptfilename", specscriptfilename, "", "Output spectrum gnuplot script file name");
    args.AddOptionalArgument("rvcorrdatafilename", rvcorrdatafilename, "", "Output RV correction data file name");
    args.AddOptionalArgument("specdatafilename", specdatafilename, "", "Output spectrum data file name");
	args.AddOptionalArgument("rvcorrfitdatafilename", rvcorrfitdatafilename, "", "Output RV correction fit data file name");
	args.AddOptionalArgument("linesdatafilename", linesdatafilename, "", "Output object lines data file name");
	
	try {
		args.Parse(argc, argv);
		
		// we need an input wavelength calibration file ...
		if (inputWaveFile.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input object spectrum file ...        
		if (inputObjectSpectrum.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output wavelength calibration file ...
		if (outputWaveFile.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input atlas of telluric lines ...
		if (telluric_lines.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        if (inputWavelengthMaskForTelluric.empty()) {
            throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        if(RVCorrectionMethod == 1) {
			if(!LocalMaxFilterWidth || !MinPeakDepth || !DetectionThreshold || !nsigclip || !minNumberOfMatchedLines || !duplicateLineThreshold)
				throw operaException("operaTelluricWavelengthCorrection: missing parameters required for method 1", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
		}
		if(RVCorrectionMethod == 2) {
			if(!XCorrelationThreshold)
				throw operaException("operaTelluricWavelengthCorrection: missing parameters required for method 2", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (args.verbose) {
			cout << "operaTelluricWavelengthCorrection: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaTelluricWavelengthCorrection: inputObjectSpectrum = " << inputObjectSpectrum << endl;
            cout << "operaTelluricWavelengthCorrection: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
			cout << "operaTelluricWavelengthCorrection: outputWaveFile = " << outputWaveFile << endl;
			cout << "operaTelluricWavelengthCorrection: telluric_lines =" << telluric_lines << endl;
			cout << "operaTelluricWavelengthCorrection: hitranformat =" << hitranformat << endl;
			cout << "operaTelluricWavelengthCorrection: spectralResolution =" << spectralResolution << endl;
			cout << "operaTelluricWavelengthCorrection: radialVelocityRange =" << radialVelocityRange << endl;
			cout << "operaTelluricWavelengthCorrection: radialVelocityStep =" << radialVelocityStep << endl;
			cout << "operaTelluricWavelengthCorrection: normalizationBinsize =" << normalizationBinsize << endl;
            cout << "operaTelluricWavelengthCorrection: StarPlusSky = " << StarPlusSky << endl;
            cout << "operaTelluricWavelengthCorrection: starplusskyInvertSkyFiber = " << starplusskyInvertSkyFiber << endl;
			cout << "operaTelluricWavelengthCorrection: skyOverStarFiberAreaRatio = " << skyOverStarFiberAreaRatio << endl;
            cout << "operaTelluricWavelengthCorrection: inputWavelengthMaskForTelluric = " << inputWavelengthMaskForTelluric << endl;
            cout << "operaTelluricWavelengthCorrection: RVCorrectionMethod = " << RVCorrectionMethod << endl;
            if(RVCorrectionMethod == 1) {
				cout << "operaTelluricWavelengthCorrection: LocalMaxFilterWidth =" << LocalMaxFilterWidth << endl;
				cout << "operaTelluricWavelengthCorrection: MinPeakDepth =" << MinPeakDepth << endl;
				cout << "operaTelluricWavelengthCorrection: DetectionThreshold =" << DetectionThreshold << endl;
				cout << "operaTelluricWavelengthCorrection: nsigclip =" << nsigclip << endl;
				cout << "operaTelluricWavelengthCorrection: emissionSpectrum =" << emissionSpectrum << endl;
				cout << "operaTelluricWavelengthCorrection: minNumberOfMatchedLines =" << minNumberOfMatchedLines << endl;
				cout << "operaTelluricWavelengthCorrection: duplicateLineThreshold =" << duplicateLineThreshold << endl;
			} else if(RVCorrectionMethod == 2) {
				cout << "operaTelluricWavelengthCorrection: XCorrelationThreshold =" << XCorrelationThreshold << endl;
				cout << "operaTelluricWavelengthCorrection: useFitToFindMaximum = " << useFitToFindMaximum << endl;
			}
            if(ordernumber != NOTPROVIDED) cout << "operaTelluricWavelengthCorrection: ordernumber = " << ordernumber << endl;
            if(args.plot) {
                cout << "operaTelluricWavelengthCorrection: rvcorrsplotfilename = " << rvcorrsplotfilename << endl;
                cout << "operaTelluricWavelengthCorrection: specplotfilename = " << specplotfilename << endl;
                cout << "operaTelluricWavelengthCorrection: rvcorrscriptfilename = " << rvcorrscriptfilename << endl;
                cout << "operaTelluricWavelengthCorrection: specscriptfilename = " << specscriptfilename << endl;
                cout << "operaTelluricWavelengthCorrection: rvcorrdatafilename = " << rvcorrdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: rvcorrfitdatafilename = " << rvcorrfitdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: specdatafilename = " << specdatafilename << endl;
            }
		}
		
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputObjectSpectrum);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile);
        
        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        if (args.verbose) cout << "operaTelluricWavelengthCorrection: minorder = "<< minorder << " maxorder = " << maxorder << endl;
        
        spectralOrders.setWavelengthsFromCalibration(minorder, maxorder);
        
		// Read telluric lines database lambda vs. intensity
		if (args.debug) cout << "operaTelluricWavelengthCorrection: reading telluric lines database " << telluric_lines << endl;
		operaSpectrum telluricLines;
		if (hitranformat) telluricLines = readTelluricLinesHITRAN(telluric_lines);
		else telluricLines = readTelluricLinesRaw(telluric_lines);
        
        // Ignore orders outside of our telluric mask
        spectralOrders.getOrdersByWavelengthRange(operaWavelengthRange(telluricLines.firstwl(), telluricLines.lastwl()), minorder, maxorder);
        if (args.verbose) cout << "operaTelluricWavelengthCorrection: useful minorder = " << minorder << " maxorder = " << maxorder << endl;
        
        // Correct for flat-field
        if (!inputFlatFluxCalibration.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputFlatFluxCalibration);
            spectralOrders.correctFlatField(minorder, maxorder, StarPlusSky, starplusskyInvertSkyFiber, skyOverStarFiberAreaRatio);
        }
       
        // Initialize rvshift to zero, so if telluric SNR is low the RV shift will not do anything.
        double rvshift = 0;
        double rvshifterror = 0;
        
        /* 
         *   Below one can use one of the following methods to measure the RV shift :
         * (RVCorrectionMethod=1): detect lines in observed spectrum and match them with telluric HITRAN lines
         * (RVCorrectionMethod=2): cross-correlation between a synthetic telluric and the observed spectra
         */
        
        switch (RVCorrectionMethod) {
                
            case 1: {
                // Detect absorption lines in the observed spectrum within telluric regions defined in inputWavelengthMaskForTelluric
				operaSpectralLineList objectLines = spectralOrders.detectSpectralLinesWithinWavelengthMask(inputWavelengthMaskForTelluric, minorder, maxorder,true, normalizationBinsize, spectralResolution,emissionSpectrum,LocalMaxFilterWidth,MinPeakDepth,DetectionThreshold,nsigclip);
				if (args.verbose) cout << "operaTelluricWavelengthCorrection: " << objectLines.size() << " lines detected in object spectrum" << endl;
				
				if(args.debug){
					for (unsigned l=0; l<objectLines.size(); l++) {
						if(l == 0 || objectLines.getwavelength(l) - objectLines.getwavelength(l-1) > duplicateLineThreshold) cout << objectLines.getwavelength(l) << " " << objectLines.getflux(l) << " " << objectLines.getsigma(l) << endl;
						else cout << "skipping line, too close to previous line: " << objectLines.getwavelength(l) << " " << objectLines.getflux(l) << " " << objectLines.getsigma(l) << endl;
					}
				}
				
				if (args.verbose) cout << "operaTelluricWavelengthCorrection: calculating RV shift by matching lines for radialVelocityRange=" << radialVelocityRange << " km/s and radialVelocityStep=" << radialVelocityStep << " km/s" << endl;
                
                // Note: spectralResolution can be updated with measurements from Object spectrum! -- could use median of all orders
                operaVector telluricMatchedWavelengths;
                operaSpectrum objectMatchedLines;
                operaVector radialVelocities;
                
                matchTelluricLines(telluricLines, objectLines, telluricMatchedWavelengths, objectMatchedLines, radialVelocities, spectralResolution, radialVelocityRange, duplicateLineThreshold);
                
                if(telluricMatchedWavelengths.size() > minNumberOfMatchedLines) {
                    rvshift = Median(radialVelocities);
                    rvshifterror = MedianStdDev(radialVelocities, rvshift);
                }
                
                if (args.verbose) cout << "operaTelluricWavelengthCorrection: (Line Match Method) Radial Velocity correction = " << rvshift << " +/- " << rvshifterror << " km/s" << endl;
                
                if (!specdatafilename.empty()) {
					ofstream fspecdata(specdatafilename.c_str());
					for(unsigned i=0; i<telluricLines.size(); i++) fspecdata << telluricLines.getwavelength(i) << " " << telluricLines.getflux(i) << endl;
				}
				
				if (!linesdatafilename.empty()) {
					ofstream flinesdata(linesdatafilename.c_str());
					for(unsigned i=0; i<objectLines.size(); i++) flinesdata << objectLines.getwavelength(i) << " " << objectLines.getflux(i) << endl;
				}
				
                if (!rvcorrdatafilename.empty()) {
					ofstream frvcorrdata(rvcorrdatafilename.c_str());
                    for (unsigned l=0; l<telluricMatchedWavelengths.size(); l++) {
                        double linewidth = telluricMatchedWavelengths[l]/spectralResolution;
                        frvcorrdata << l << setprecision(5) << fixed << " " << telluricMatchedWavelengths[l] << " " <<  objectMatchedLines.getwavelength(l) << " " << objectMatchedLines.getflux(l) << " " << linewidth << " " << radialVelocities[l] << endl;
                    }
                }
                
                if (!rvcorrfitdatafilename.empty()) {
					operaVector rvVector;
					operaVector probDensity;
					operaVector histWaveVector;
					
					generateHistogramData(telluricMatchedWavelengths, radialVelocities, radialVelocityRange, radialVelocityStep, rvVector, probDensity, histWaveVector);
					ofstream frvcorrfitdata(rvcorrfitdatafilename.c_str());
                    for (unsigned i=0; i<histWaveVector.size(); i++) {
                        frvcorrfitdata << i << " " << histWaveVector[i] << " " <<  rvVector[i] << " " << probDensity[i] << endl;
                    }
                }
                
                if (!specdatafilename.empty() && !rvcorrdatafilename.empty() && !specscriptfilename.empty()) {
					GenerateTelluricLineMatchPlot(specscriptfilename, specplotfilename, specdatafilename, rvcorrdatafilename);
                }
                
                if (!rvcorrdatafilename.empty() && !rvcorrfitdatafilename.empty() && !rvcorrscriptfilename.empty()) {
					GenerateTelluricRVCorrPlot(rvcorrscriptfilename, rvcorrsplotfilename, rvcorrdatafilename, rvcorrfitdatafilename, rvshift, rvshifterror);
                }
                
                break;
            }
            case 2: {
                /* -- E. Martioli Aug 17 2015 -- below is the old way of calculating RVshift for telluric wavelength
                 correction, where it uses a cross-correlation between observed and telluric spectra rather than
                 individual line positions. This method seems to be biased by~300m/s because the line profiles are
                 not symmetric.
                 */
                
                // Get the object spectrum within telluric regions defined in inputWavelengthMaskForTelluric (used for method 2)
				operaSpectrum objectSpectrum = spectralOrders.getSpectrumWithinTelluricMask(inputWavelengthMaskForTelluric, minorder, maxorder, true, normalizationBinsize);
				if(args.debug){
					for (unsigned l=0; l<objectSpectrum.size(); l++) {
						cout << objectSpectrum.getwavelength(l) << " " << objectSpectrum.getflux(l) << " " << objectSpectrum.getvariance(l) << endl;
					}
				}
        
				// Spectrum plot: plot observed and reference telluric spectra.
				if (!specdatafilename.empty()) {
					//Use HITRAN lines to generate synthetic spectrum sampled to the same points as the object spectrum
					operaVector hitranTelluricSpectrum = generateSyntheticTelluricSpectrumUsingLineProfile(telluricLines, objectSpectrum.wavelengthvector(), spectralResolution, GAUSSIAN);
					
					ofstream fspecdata(specdatafilename.c_str());
					for(unsigned i=0; i<objectSpectrum.size(); i++) fspecdata << objectSpectrum.getwavelength(i) << " " << objectSpectrum.getflux(i) << " " << hitranTelluricSpectrum[i] << endl;
				}
                
                ofstream frvcorrdata;
				ofstream frvcorrfitdata;
				if (!rvcorrdatafilename.empty()) frvcorrdata.open(rvcorrdatafilename.c_str());
				if (!rvcorrfitdatafilename.empty()) frvcorrfitdata.open(rvcorrfitdatafilename.c_str());
                
                if (args.verbose) cout << "operaTelluricWavelengthCorrection: calculating cross-correlation for radialVelocityRange=" << radialVelocityRange << " km/s and radialVelocityStep=" << radialVelocityStep << " km/s" << endl;
                double maxcorr=-BIG, chisqr=0;
                
                bool validXCorrelation = calculateRVShiftByXCorr(telluricLines, objectSpectrum, radialVelocityRange, radialVelocityStep, XCorrelationThreshold, rvshift, rvshifterror, maxcorr, frvcorrdata, frvcorrfitdata, spectralResolution, useFitToFindMaximum, chisqr);
                
                if(!validXCorrelation) {
                    rvshift = 0;
                    rvshifterror = 0;
                }
                
                if(args.verbose) cout << "operaTelluricWavelengthCorrection: (X-Corr Method) Radial Velocity correction = " << rvshift << " +/- " << rvshifterror << " km/s, maxXCorr=" << maxcorr << ", chisqr=" << chisqr << "\n" << endl;
                
                // Telluric spectrum plot:
                if (!specdatafilename.empty() && !specscriptfilename.empty()) {
                    GenerateTelluricSpecPlot(specscriptfilename, specplotfilename, specdatafilename);
                }
                
                // Telluric wavelength correction plot:
                if (frvcorrdata.is_open() && frvcorrfitdata.is_open()) {
                    frvcorrdata.close();
                    frvcorrfitdata.close();
                    if (!rvcorrscriptfilename.empty()) {
                        GenerateTelluricXCorrelationPlot(rvcorrscriptfilename, rvcorrsplotfilename, rvcorrdatafilename, rvcorrfitdatafilename);
                    }
                }
                
                break;
            }
                
            default:
                break;
        }
        
        FormatHeader outputheader("Telluric Radial Velocity Correction");
        outputheader << "radialvelocity (km/s)" << "radialvelocityerror (km/s)" << newline;
        FormatData outputdata;
        outputdata << rvshift << rvshifterror << endl;
        operaIOFormats::WriteCustomFormat("tell", outputheader, outputdata, outputWaveFile);

    }
	catch (operaException e) {
		cerr << "operaTelluricWavelengthCorrection: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaTelluricWavelengthCorrection: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

// Generate multiple plots containing statistical info about telluric wavelength correction
void GenerateTelluricXCorrelationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string cleanDataFileName)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
	
    fgnu << "reset" << endl;

    fgnu << "\nset xlabel \"Radial Velocity (km/s)\"" << endl;
    fgnu << "set ylabel \"cross-correlation\"" << endl;
    fgnu << "set pointsize 1.5" << endl;

    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced mono solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << endl;
        
        fgnu << "plot \"" << dataFileName << "\" u 1:2 t \"gaussian fit\" w l lt 3 lw 2, ";
        fgnu << "\"" << dataFileName << "\" u 1:3:4 t \"XCorr data\" w yerr pt 6, ";
        fgnu << "\"" << cleanDataFileName << "\" u 1:4:2:5 t \"fit data\" w xyerr pt 7" << endl;
        
        fgnu << "\n#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    } else {
        fgnu << endl;
        
        fgnu << "plot \"" << dataFileName << "\" u 1:2 t \"gaussian fit\" w l lt 3 lw 2, ";
        fgnu << "\"" << dataFileName << "\" u 1:3:4 t \"XCorr data\" w yerr pt 6, ";
        fgnu << "\"" << cleanDataFileName << "\" u 1:4:2:5 t \"fit data\" w xyerr pt 7" << endl;
        
        fgnu << endl;
        
        fgnu << "\n#set terminal postscript enhanced mono solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();

    if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

// Generate multiple plots containing statistical info about telluric wavelength correction
void GenerateTelluricRVCorrPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string histDataFileName, float rvshift, float rvshifterror)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
    remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
    ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "\nreset" << endl;

    if(!outputPlotEPSFileName.empty()) {
        
        fgnu << "\nNX=1; NY=2" << endl;
        fgnu << "\nDX=0.1; DY=0.1; SX=0.98; SY=0.42" << endl;
        fgnu << "\nset bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY" << endl;
        
        fgnu << "\nset size SX*NX+DX*2,SY*NY+DY*2" << endl;
        
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nmedianRV = " << rvshift << endl;
        fgnu << "\nRVerr = " << rvshifterror << endl;
        
        fgnu << "\nset multiplot" << endl;
        
        fgnu << "\nset size 0.87*SX,1.35*SY" << endl;
        fgnu << "set origin DX,DY+2*SY/3" << endl;
        fgnu << "set xlabel \"Radial velocity (km/s)\"" << endl;
        fgnu << "set ylabel \"Probability Density\"" << endl;
        
        fgnu << "\nset label \"{/Symbol D}RV = " << fixed << setprecision(3) <<  rvshift << "+/-" << rvshifterror << "\" at " << rvshift + 2*rvshifterror << ",0.125" << endl;
        
        fgnu << "\nset arrow from medianRV,0.1 to medianRV,0 lt 3 lw 2" << endl;
        fgnu << "set arrow from medianRV-RVerr,0 to medianRV-RVerr,0.1 nohead lt 2 lw 1" << endl;
        fgnu << "set arrow from medianRV+RVerr,0 to medianRV+RVerr,0.1 nohead lt 2 lw 1" << endl;
        
        fgnu << "\nplot \"" << histDataFileName << "\" u 3:4 t \"RV shift PDF\" w boxes" << endl;
        
        fgnu << "\nset size 0.87*SX,0.45*SY" << endl;
        fgnu << "set origin DX,DY" << endl;
        fgnu << "set xlabel \"{/Symbol l} (nm)\"" << endl;
        fgnu << "set ylabel \"Radial velocity (km/s)\"" << endl;
        
        fgnu << "\nf(x) = medianRV" << endl;
        fgnu << "fl(x) = medianRV - RVerr" << endl;
        fgnu << "fh(x) = medianRV + RVerr" << endl;
        
        fgnu << "\nset yrange[medianRV - 8*RVerr:medianRV + 8*RVerr]" << endl;
        
        fgnu << "\nplot \"" << dataFileName << "\" u 2:6 notitle w p pt 7, f(x) notitle w l lt 3 lw 2, fl(x) notitle w l lt 2 lw 1, fh(x) notitle w l lt 2 lw 1" << endl;
        
        fgnu << "\nunset multiplot" << endl;
        
        fgnu << "\n#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    
    } else {
        
        fgnu << "\nNX=1; NY=2" << endl;
        fgnu << "\nDX=0.1; DY=0.1; SX=0.98; SY=0.42" << endl;
        fgnu << "\nset bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY" << endl;
        
        fgnu << "\nset size SX*NX+DX*2,SY*NY+DY*2" << endl;
        
        fgnu << "\nmedianRV = " << rvshift << endl;
        fgnu << "\nRVerr = " << rvshifterror << endl;
        
        fgnu << "\nset multiplot" << endl;
        
        fgnu << "\nset size 0.87*SX,1.35*SY" << endl;
        fgnu << "set origin DX,DY+2*SY/3" << endl;
        fgnu << "set xlabel \"Radial velocity (km/s)\"" << endl;
        fgnu << "set ylabel \"Probability Density\"" << endl;
        
        fgnu << "\nset label \"{/Symbol D}RV = " << fixed << setprecision(3) << rvshift << "+/-" << rvshifterror << "\" at " << rvshift + 2*rvshifterror << ",0.125" << endl;
        
        fgnu << "\nset arrow from medianRV,0.1 to medianRV,0 lt 3 lw 2" << endl;
        fgnu << "set arrow from medianRV-RVerr,0 to medianRV-RVerr,0.1 nohead lt 2 lw 1" << endl;
        fgnu << "set arrow from medianRV+RVerr,0 to medianRV+RVerr,0.1 nohead lt 2 lw 1" << endl;
        
        fgnu << "\nplot \"" << histDataFileName << "\" u 3:4 t \"RV shift PDF\" w boxes" << endl;
        
        fgnu << "\nset size 0.87*SX,0.45*SY" << endl;
        fgnu << "set origin DX,DY" << endl;
        fgnu << "set xlabel \"{/Symbol l} (nm)\"" << endl;
        fgnu << "set ylabel \"Radial velocity (km/s)\"" << endl;
        
        fgnu << "\nf(x) = medianRV" << endl;
        fgnu << "fl(x) = medianRV - RVerr" << endl;
        fgnu << "fh(x) = medianRV + RVerr" << endl;
        
        fgnu << "\nset yrange[medianRV - 8*RVerr:medianRV + 8*RVerr]" << endl;
        
        fgnu << "\nplot \"" << dataFileName << "\" u 2:6 notitle w p pt 7, f(x) notitle w l lt 3 lw 2, fl(x) notitle w l lt 2 lw 1, fh(x) notitle w l lt 2 lw 1" << endl;
        
        fgnu << "\nunset multiplot" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    }
    
    fgnu.close();
    
    if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

// Generate 2D plot for spectra of atlas + identified lines
void GenerateTelluricLineMatchPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename, string matchdatafilename)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
    remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
    ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "#unset key" << endl;
    
    fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;
    
    fgnu << "set ylabel \"norm flux\"" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << endl;
        
        fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Telluric Reference (HITRAN)\" w l lt 3, ";
        fgnu << "\"" << matchdatafilename << "\" u 3:(1 - (1-$4)/2):5 t \"Matched lines\" w xerr pt 7 lt 2" << endl;

        fgnu << endl;
        
        fgnu << "\n#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    } else {
        fgnu << endl;
        
        fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Telluric Reference (HITRAN)\" w l lt 3, ";
        fgnu << "\"" << matchdatafilename << "\" u 3:(1 - (1-$4)/2):5 t \"Matched lines\" w xerr pt 7 lt 2" << endl;
        
        fgnu << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

// Generate 2D plot for spectra of atlas + comparison
void GenerateTelluricSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename)
{    
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "#unset key" << endl;

    fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;

    fgnu << "set ylabel \"norm flux\"" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << endl;
        
        fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Object Spectrum\" w l lt 4, ";
        fgnu << "\"" << specdatafilename << "\" u 1:3 t \"Telluric Reference (HITRAN)\" w l lt 3";
        fgnu << endl;
        
        fgnu << "\n#set terminal x11" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#replot" << endl;
    } else {
        fgnu << endl;
        
        fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Object Spectrum\" w l lt 4, ";
        fgnu << "\"" << specdatafilename << "\" u 1:3 t \"Telluric Reference (HITRAN)\" w l lt 3";
        fgnu << endl;
        fgnu << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

operaSpectrum readTelluricLinesHITRAN(string telluric_database_file)
{
	const double N_OVER_V = TYPICAL_PRESSURE_AT_MAUNAKEA/(TYPICAL_TEMPERATURE_AT_MAUNAKEA*k_BOLTZMANN_CONSTANT);
	
    operaSpectrum telluricLines;
    igzstream astream(telluric_database_file.c_str());
	if (astream.is_open()) {
		string dataline;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') { // skip blank lines and comments
				double wave_number;
				float intensity;
				if(!sscanf(dataline.c_str(), "%*d %lf %G %*[^\n]", &wave_number, &intensity)) continue; //skip over bad line
				//if (intensity < 1.15e-26) continue; //skip over lines below a certain threshold
				double wavelength_in_nm = 1e7/wave_number;
				telluricLines.insert(convertVacuumToAirWavelength(wavelength_in_nm*10)/10, ((double)intensity/(N_OVER_V*1e-6))/TYPICAL_ATMOSPHERE_PATH_LENGTH);
			}
		}
		astream.close();
	}
	telluricLines.reverse();
    if (args.verbose) {
		if (telluricLines.empty()) printf("          [Telluric] no lines found in telluric database.\n");
		else printf("          [Telluric] %d lines found wl0=%.2f wlc=%.2f wlf=%.2f\n", telluricLines.size(), telluricLines.firstwl(), telluricLines.midwl(), telluricLines.lastwl());
	}
	return telluricLines;
}

// Read the entire set of telluric lines in HITRAN database
operaSpectrum readTelluricLinesRaw(string telluric_database_file)
{
	operaSpectrum telluricLines;
    igzstream astream(telluric_database_file.c_str());
	if (astream.is_open()) {
		string dataline;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') { // skip blank lines and comments
				double wavelength, intensity;
				istringstream ss(dataline);
				ss >> wavelength >> intensity;
				telluricLines.insert(wavelength, intensity);
			}
		}
		astream.close();
	}
    if (args.verbose) {
		if (telluricLines.empty()) printf("          [Telluric] no lines found in telluric database.\n");
		else printf("          [Telluric] %d lines found wl0=%.2f wlc=%.2f wlf=%.2f\n", telluricLines.size(), telluricLines.firstwl(), telluricLines.midwl(), telluricLines.lastwl());
	}
	return telluricLines;
}

// Finds the first and last index between wl0 and wlf. Assumes wavelength vector is in increasing order.
void getWavelengthSubrange(const operaVector& wavelength, double wl0, double wlf, unsigned& startindex, unsigned& endindex)
{
    endindex = startindex = 0;
    operaWavelengthRange wlrange(wl0, wlf);
	for (unsigned i = 0; i < wavelength.size(); i++) {
		if (wlrange.contains(wavelength[i])) {
			if (endindex == 0) endindex = startindex = i;
			endindex++;
		}
	}
}

// Generates a spectrum along the points in wavelengthVector by using a Gaussian profile to fit telluricLines.
operaVector generateSyntheticTelluricSpectrumUsingLineProfile(const operaSpectrum& telluricLines, const operaVector& wavelengthVector, double resolution, ProfileMethod profile)
{
	operaVector outputSpectrum(wavelengthVector.size());
	outputSpectrum.fill(1.0); //Initialize outputSpectrum to uniform 1.0
	for(unsigned i=0; i<wavelengthVector.size(); i++) {
		double gaussianWidth = (wavelengthVector[i]/resolution);
		double wl0 = wavelengthVector[i] - 5*gaussianWidth;
		double wlf = wavelengthVector[i] + 5*gaussianWidth;
		unsigned startindex, endindex;
		getWavelengthSubrange(telluricLines.wavelengthvector(), wl0, wlf, startindex, endindex);
		if(args.debug) cout << "operaTelluricWavelengthCorrection: " << i << " gaussianwidth=" << gaussianWidth << " wl0=" << wl0 << " wlf=" << wlf << " nlinesInRange=" << endindex-startindex << endl;

		//Uncomment the following line for old functionality (recreate bug?)
		//if(endindex > 0) endindex--;
		for(unsigned j=startindex; j<endindex; j++) {
			double gamma = 1.0; // set to one for testing
			double opticaldepth = telluricLines.getflux(j);
			if (profile == VOIGT || profile == GAUSSIAN) opticaldepth *= exp(-((telluricLines.getwavelength(j) - wavelengthVector[i])*(telluricLines.getwavelength(j) - wavelengthVector[i])/(2*gaussianWidth*gaussianWidth)))/(sqrt(2*M_PI)*gaussianWidth);
			if (profile == VOIGT || profile == LORENTZ) opticaldepth *= (1/(M_PI*gamma))*(gamma*gamma)/(gamma*gamma + (telluricLines.getwavelength(j) - wavelengthVector[i])*(telluricLines.getwavelength(j) - wavelengthVector[i]));
			outputSpectrum[i] *= exp(-opticaldepth);
		}
	}
	return outputSpectrum;
}

bool calculateRVShiftByXCorr(const operaSpectrum& telluricLines, const operaSpectrum& objectSpectrum, double radialVelocityRange, double radialVelocityStep, double threshold, double& maxRV, double& sigRV, double& maxcorr, ofstream& frvcorrdata, ofstream& frvcorrfitdata, double spectralResolution, bool useFitToFindMaximum, double& chisqr)
{
    int jmax = -1;
	maxcorr = 0;
	maxRV = 0;
	sigRV = 0;
	chisqr = 0;
    
    operaVector crosscorrelation;
    operaVector crosscorrerror;
    operaVector dRV;
    
    double xcorrerror = 2e-04; //why this value in particular?
    
    for(double deltaRV = -radialVelocityRange/2.0; deltaRV <= radialVelocityRange/2.0; deltaRV+=radialVelocityStep) {
        operaVector telluricSpectrumWavelength;
        // Initalize telluricSpectrum wavelength with wavelength of objectSpectrum shifted by deltaRV
        for (unsigned i=0; i<objectSpectrum.size(); i++) {
            double DWavelength = deltaRV * objectSpectrum.getwavelength(i) / SPEED_OF_LIGHT_KMS;
            telluricSpectrumWavelength.insert(objectSpectrum.getwavelength(i) + DWavelength);
        }
        // Generate a spectrum in telluricSpectrum along points in wavelength vector using the provided telluricLines
        operaVector telluricSpectrumFlux = generateSyntheticTelluricSpectrumUsingLineProfile(telluricLines, telluricSpectrumWavelength, spectralResolution, GAUSSIAN);
        
        // Calculate the x-corr between the generated shifted telluric spectrum and the object spectrum
        double xcorr = operaCrossCorrelation(telluricSpectrumFlux.size(), objectSpectrum.flux_ptr(), telluricSpectrumFlux.datapointer());
        if(args.debug) cout << deltaRV << " " << xcorr << endl;
        
        // Check if this is the highest x-corr we have found so far, but filter out values under threshold
        if(xcorr > threshold && (jmax < 0 || xcorr > maxcorr)) {
            maxcorr = xcorr;
            maxRV = deltaRV;
            sigRV = radialVelocityStep;
            jmax = crosscorrelation.size();
        }
        
        crosscorrelation.insert(xcorr);
        crosscorrerror.insert(xcorrerror);
        dRV.insert(deltaRV);
    }
    
    if (jmax < 0) return false; // Didn't find any x-corr values above threshold
    
    if(useFitToFindMaximum) {
		// Set initial values for our Gaussian using the maximum x-corr.
        double a = crosscorrelation[jmax]; //Initial amplitude
        double x0 = dRV[jmax]; //Initial center
        double sig = radialVelocityRange/4.0; //Initial sigma
        double ea;
        double ex0;
        double esig;
        double fitchisqr;
        
        // Update the initial values and get errors for each along with the fit chi-squared.
        operaMPFitGaussian(crosscorrelation.size(), dRV.datapointer(), crosscorrelation.datapointer(), crosscorrerror.datapointer(), &a, &ea, &x0, &ex0, &sig, &esig, &fitchisqr);
        //operaLMFitGaussian(np, peakXdata, peakYdata, &a, &x0, &sig, &chisqr);
        
        if(args.debug) {
            cout << a << "+/-" << ea << endl;
            cout << x0 << "+/-" << ex0 << endl;
            cout << sig << "+/-" << esig <<  " fitchisqr=" << fitchisqr << endl;
        }
        
        // For plotting
        if(frvcorrdata.is_open()) {
            for(unsigned j=0; j<crosscorrelation.size(); j++) {
                double x = (double)dRV[j];
                double gaussfunc = a*exp(-(x-x0)*(x-x0)/(2*sig*sig));
                frvcorrdata << dRV[j] << " " <<  gaussfunc << " " <<  crosscorrelation[j] << " " <<  crosscorrerror[j] << " " << crosscorrelation[j] - gaussfunc << endl;
            }
            frvcorrdata << endl;
        }
        if(frvcorrfitdata.is_open()) {
            frvcorrfitdata  << x0 << " " << ex0 << " " <<  a << " " <<  maxcorr  <<  " " <<  crosscorrerror[jmax]  << " " << maxcorr - a << endl;
        }

        maxcorr = a;
        maxRV = x0;
        sigRV = ex0;
        chisqr = fitchisqr;
    } else {
		// For plotting
        if(frvcorrdata.is_open()) {
            for(unsigned j=0; j<crosscorrelation.size(); j++) {
                frvcorrdata  << dRV[j] << " " <<  crosscorrelation[j] << " " <<  crosscorrelation[j] <<  " " << crosscorrerror[j] << " " << 0.0 << endl;
            }
            frvcorrdata << endl;
        }
        if(frvcorrfitdata.is_open()) {
            frvcorrfitdata  << maxRV << " " <<  sigRV << " " << maxcorr << " " <<  maxcorr  <<  " " <<  crosscorrerror[jmax] <<  " " << 0.0 << endl;
        }
    }
    return true;
}

// Function to match telluric lines
void matchTelluricLines(const operaSpectrum& atlasLines, const operaSpectralLineList& objectLines, operaVector& atlasMatchedWavelengths, operaSpectrum& objectMatchedLines, operaVector& radialVelocities, double spectralResolution, double radialVelocityRange, double duplicateLineThreshold)
{
    unsigned match_index = 0; // Keep track of the last matched line in the atlas so we can skip ahead to there
    for (unsigned j=0; j<objectLines.size(); j++) {
		double objectwl = objectLines.getwavelength(j);
		if(j > 0 && objectwl - objectLines.getwavelength(j-1) < duplicateLineThreshold) continue; // Skip over object lines which are most likely duplicates
        
        const operaVector& atlaswlvector = atlasLines.wavelengthvector();
        match_index = findClosestInSortedRange(atlaswlvector, objectwl, match_index, atlaswlvector.size()); // Search atlaswlvector in range [match_index, size) for the closest wavelength to objectwl
        
        if(match_index < atlaswlvector.size()) {
			double atlaswl = atlaswlvector[match_index];
			double deltarv = calculateDeltaRadialVelocityInKPS(atlaswl, objectwl);
			
			// Make sure that this atlas line (and only this atlas line) is within the rv range of the object line
			if(fabs(deltarv) >= radialVelocityRange/2.0) continue;
			if(match_index > 0 && fabs(calculateDeltaRadialVelocityInKPS(atlaswlvector[match_index-1], objectwl)) < radialVelocityRange/2.0) continue;
			if(match_index + 1 < atlaswlvector.size() && fabs(calculateDeltaRadialVelocityInKPS(atlaswlvector[match_index+1], objectwl)) < radialVelocityRange/2.0) continue;

			double linewidth = objectwl/spectralResolution;
			if(fabs(atlaswl - objectwl) < linewidth) {					
				double objectflux = objectLines.getflux(j);
				atlasMatchedWavelengths.insert(atlaswl);
				objectMatchedLines.insert(objectwl, objectflux);
				radialVelocities.insert(deltarv);
			}
		}
		else break; // No atlas lines remaining to match
	}
}

void generateHistogramData(const operaVector& telluricMatchedWavelengths, const operaVector& radialVelocities, double radialVelocityRange, double radialVelocityStep, operaVector& rvVector, operaVector& probDensity, operaVector& wavelengthVector) {
    unsigned nTotal = 0;
    for(double deltaRV = -radialVelocityRange/2.0; deltaRV <= radialVelocityRange/2.0; deltaRV+=radialVelocityStep) {
        for (unsigned l=0; l<telluricMatchedWavelengths.size(); l++) {
            if(radialVelocities[l] > deltaRV - radialVelocityStep/2 && radialVelocities[l] <= deltaRV + radialVelocityStep/2) {
                nTotal++;
            }
        }
    }
    
    for(double deltaRV = -radialVelocityRange/2.0; deltaRV <= radialVelocityRange/2.0; deltaRV+=radialVelocityStep) {
        unsigned npbin = 0;
        double meanwl = 0;
        
        for (unsigned l=0; l<telluricMatchedWavelengths.size(); l++) {
            if(radialVelocities[l] > deltaRV - radialVelocityStep/2 && radialVelocities[l] <= deltaRV + radialVelocityStep/2) {
                meanwl += telluricMatchedWavelengths[l];
                npbin++;
            }
        }
        
        rvVector.insert(deltaRV);
        probDensity.insert(static_cast<double>(npbin)/nTotal);
        wavelengthVector.insert(meanwl/npbin);
    }
}
