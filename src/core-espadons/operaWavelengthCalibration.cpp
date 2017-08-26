/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaWavelengthCalibration
 Version: 1.0
 Description: Wavelength Calibration 
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

#include <iomanip>
#include "libraries/operaIOFormats.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaSpectralTools.h"
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/gzstream.h"							// for gzstream - read compressed reference spectra
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"
#include "core-espadons/operaWavelengthCalibration.h"

/*! \file operaWavelengthCalibration.cpp */
/*! \ingroup core */

using namespace std;

/*
 
 Below it follows in a few words a 1st-pass for the wavelength calibration
 algorithm.
 
 For each spectral order do the following steps:
 
 1. Read ThAr raw spectrum; intensity versus distance in pixel units:
 I(d) vs. d
 
 2. Measure total distance "D" (in pixel units) covered by the order. This
 is given by the line integral of the polynomial that describes the center
 of the order.
 
 3. Read wavelength range covered by the order: wl0,wlf
 
 4. Calculate first order solution:
 wl = f(d), where f(d) as first order is given by
 
 f(d) = wl0 + ((wlf - wl0)/D)*d
 
 assuming f(d=0) = wl0.
 
 5. Read ThAr atlas of spectral lines within the range covered by the order
 [wl0:wlf]. The atlas consists of line wavelength (l_wl), error (l_wlerr),
 and line relative intensity (l_i).
 
 7. Once we have a table of l_i, l_d, and l_derr, then we can calculate the
 maximum cross-correlation between this and the atlas data to identify the
 lines. The identification usually doesn't go one-by-one, so we will end up
 having to do some cleaning for either the undetected or over-detected
 lines.
 
 8. Now one can use the table (l_d+/-l_derr) versus (l_wl+/-l_wlerr) to
 find the wavelength solution by fitting a polynomial to these data.
 
 note that the polynomial should be an update to the first-order solution.
 
 The update is intended for two reasons:
 
 First because I have experienced before that the higher order terms are so
 small when compared to the first order that the fitting routine can get in
 trouble to find good solutions.
 
 Another reason is that we want to use our first solution to exclude
 outliers and then run the fitting again as many times as it gets to give
 the best solution. So, we will have to update our solution as we get our
 dataset cleaner or as we gather more information.
 
 */

/* 
 * operaWavelengthCalibration
 * \author Doug Teeple
 * \brief wavelength calibration.
 * \arg argc
 * \arg argv
 * \note --outputWave=...
 * \note --atlas_lines=...
 * \note --thcal=...
 * \note --geom=...
 * \note --wlcal_initialguess=...
 * \note --binsize=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

operaArgumentHandler args;

int main(int argc, char *argv[])
{
	string outputWave;
	string outputResolution;
	string atlas_lines;
    string atlas_spectrum;
    string uncalibrated_lines;
    string uncalibrated_spectrum;
    double uncalibrated_linewidth = 1.5;
	string geometryfilename;
	string wlcal_initialguess;
    string inputLineSetFilename;
    
	bool parseSolution = false;
    bool normalizeUncalibratedSpectrum = false;
    unsigned normalizationBinSize = 150;
    
    int ordernumber = NOTPROVIDED;
    int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;
    
    string ordersplotfilename;
    string specplotfilename;
    string atlasdatafilename;
	string compdatafilename;
	string linesdatafilename;
	string ordersdatafilename;
	string ordersscriptfilename;
    string specscriptfilename;
    bool generate3DPlot = false;
    bool subtractCentralWavelength = true;
    bool interactive = false;
    
    /*
     * The parameters below we don't know yet whether they will be input
     */
    double DetectionThreshold = 0.05;    // threshold to regulate the sensitivity of line detection. Must be between 0 and 1.
    double LocalMaxFilterWidth = 3.0;    // parameter to set a window filter to guarantee a line is not detected twice. It's in units of line width
    double MinPeakDepth = 0.25;           // limit that also regulates the sensitity of line detection in units of noise.
    
    double ParRangeSizeInPerCent = 0.1;  // define the range within which a coefficient will be changed to calcuate the x-correlation
    unsigned NpointsPerPar = 300;        // define the number of times a coefficient is changed
    
    unsigned maxNIter = 30;                 // maximum number of iterations for shrinking the acceptable mismatch
    unsigned minNumberOfLines = 40;         // minimum number of lines to stop shrinking the acceptable mismatch difference between atlas and comparison
    unsigned maxorderofpolynomial = 4;      // maximum degree of polynomial for wavelength solution
    double dampingFactor = 0.90;            // Damping factor to shrink the size of the quantity acceptableMismatch on each iteration.  This factor may be set between 0 to 1.
    double initialAcceptableMismatch = 1.0; // initial acceptable mismatch difference between atlas and comparison lines. Used for identification. In units of line width.
    double nsigclip = 3.0;						// Threshold (in units of rms) for clipping lines.
    
    int nOrdersToSearchAround = 2;
    int referenceOrder = 0;
    
    args.AddRequiredArgument("outputWaveFile", outputWave, "Output wavelength calibration file to store final solution");
    args.AddOptionalArgument("outputResolutionFile", outputResolution, "", "Output spectral resolution file to store additional details");
    args.AddRequiredArgument("inputGeomFile", geometryfilename, "Input geometry calibration file");
    args.AddOptionalArgument("inputWaveFile", wlcal_initialguess, "", "Input wavelength calibration file (initial guess)");
    args.AddOptionalArgument("inputLineSetFilename", inputLineSetFilename, "", "Input wavelength and distance list of manually identified lines");
    args.AddOptionalArgument("atlas_lines", atlas_lines, "", "File containing the atlas of reference lines");
    args.AddOptionalArgument("atlas_spectrum", atlas_spectrum, "", "File containing the spectrum of reference atlas, overrides atlas_lines");
    args.AddOptionalArgument("uncalibrated_lines", uncalibrated_lines, "", "File containing the uncalibrated raw lines"); // operaExtractSpactralLines does this for us
    args.AddOptionalArgument("uncalibrated_spectrum", uncalibrated_spectrum, "", "File containing the uncalibrated raw spectrum, overrides uncalibrated_lines");
    args.AddOptionalArgument("uncalibrated_linewidth", uncalibrated_linewidth, 1.5, "Line width in pixels, used for detecting lines in uncalibrated_spectrum");
    
    args.AddOptionalArgument("normalizeUncalibratedSpectrum", normalizeUncalibratedSpectrum, false, "Normalize uncalibrated_spectrum");
    args.AddOptionalArgument("normalizationBinSize", normalizationBinSize, 150, "Binsize to be used for normalization");
    args.AddRequiredArgument("LocalMaxFilterWidth", LocalMaxFilterWidth, "To set a window filter to guarantee a line is not detected twice, in units of line width");
    args.AddRequiredArgument("DetectionThreshold", DetectionThreshold, "Threshold to regulate the sensitivity of line detection, between 0 and 1");
    args.AddRequiredArgument("MinPeakDepth", MinPeakDepth, "Limit that regulates the sensitity of line detection, in units of noise");
    args.AddOptionalArgument("parseSolution", parseSolution, false, "Vary all coefficients to refine initial solution - use only when initial guess is poor");
    args.AddRequiredArgument("ParRangeSizeInPerCent", ParRangeSizeInPerCent, "The percent range within which a coefficient will be varied to find a solution");
    args.AddRequiredArgument("NpointsPerPar", NpointsPerPar, "The number of times a coefficient will be changed to find a solution");
    args.AddRequiredArgument("initialAcceptableMismatch", initialAcceptableMismatch, "Initial acceptable mismatch difference between atlas and comparison lines, in units of line width");
    args.AddRequiredArgument("maxNIter", maxNIter, "Maximum number of iterations for shrinking the acceptable mismatch");
    args.AddRequiredArgument("dampingFactor", dampingFactor, "Damping factor between 0 and 1 to multiply the acceptable mismatch by after each iteration");
    args.AddRequiredArgument("minNumberOfLines", minNumberOfLines, "Minimum number of lines to keep when shrinking the acceptable mismatch");
    args.AddRequiredArgument("maxorderofpolynomial", maxorderofpolynomial, "Maximum degree of polynomial for wavelength solution");
    args.AddRequiredArgument("nsigclip", nsigclip, "Threshold in units of RMS used to clip matched lines which are too far apart, and to filter detected lines that deviate from the median line width");
    args.AddOptionalArgument("nOrdersToSearchAround", nOrdersToSearchAround, 0, "Number of spectral orders to search around for order shifts");
    args.AddOptionalArgument("referenceOrder", referenceOrder, 0, "Order number to be used as referece to look for order shifts");
    
    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    args.AddOptionalArgument("ordersplotfilename", ordersplotfilename, "", "Output orders plot eps file name");
    args.AddOptionalArgument("specplotfilename", specplotfilename, "", "Output spectrum plot eps file name");
    args.AddOptionalArgument("ordersdatafilename", ordersdatafilename, "", "Output orders data file name");
	args.AddOptionalArgument("atlasdatafilename", atlasdatafilename, "", "Output atlas data file name");
	args.AddOptionalArgument("linesdatafilename", linesdatafilename, "", "Output lines data file name");
	args.AddOptionalArgument("compdatafilename", compdatafilename, "", "Output comparison data file name");
	args.AddOptionalArgument("ordersscriptfilename", ordersscriptfilename, "", "Output orders gnuplot script file name");
	args.AddOptionalArgument("specscriptfilename", specscriptfilename, "", "Output spectrum gnuplot script file name");
	args.AddOptionalArgument("generate3DPlot", generate3DPlot, false, "Choose a 3D plot of the spectra instead of 2D");
	args.AddOptionalArgument("subtractCentralWavelength", subtractCentralWavelength, true, "Choose to subtract order central wavelength for plot");
	args.AddSwitch("interactive", interactive, "For interactive plots");
	
	try {
		args.Parse(argc, argv);
		
		// we need a atlas_lines lines or spectrum...
		if (atlas_lines.empty() && atlas_spectrum.empty()) {
			throw operaException("operaWavelengthCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need EITHER thorium uncalibrated lines or spectrum...
		if (uncalibrated_lines.empty() && uncalibrated_spectrum.empty()) {
			throw operaException("operaWavelengthCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a geometryfilename...
		if (geometryfilename.empty()) {
			throw operaException("operaWavelengthCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        // we need either an initial guess at a polynomial in wlcal_initialguess, or an input set of lines in file inputLineSetFilename ...
        if (wlcal_initialguess.empty() && inputLineSetFilename.empty()) {
            throw operaException("operaWavelengthCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
		if (args.verbose) {
			cout << "operaWavelengthCalibration: atlas_lines = " << atlas_lines << endl;
			cout << "operaWavelengthCalibration: atlas_spectrum = " << atlas_spectrum << endl;            
			cout << "operaWavelengthCalibration: uncalibrated_lines = " << uncalibrated_lines << endl;
            cout << "operaWavelengthCalibration: uncalibrated_spectrum = " << uncalibrated_spectrum << endl;
            cout << "operaWavelengthCalibration: uncalibrated_linewidth = " << uncalibrated_linewidth << endl;
            cout << "operaWavelengthCalibration: inputLineSetFilename = " << inputLineSetFilename << endl;
			cout << "operaWavelengthCalibration: geometryfilename = " << geometryfilename << endl;            
			cout << "operaWavelengthCalibration: wlcal_initialguess = " << wlcal_initialguess << endl; 
			cout << "operaWavelengthCalibration: outputWave = " << outputWave << endl;
			cout << "operaWavelengthCalibration: outputResolution = " << outputResolution << endl;
			cout << "operaWavelengthCalibration: parseSolution = " << parseSolution << endl;
			cout << "operaWavelengthCalibration: normalizeUncalibratedSpectrum = " << normalizeUncalibratedSpectrum << endl;
			cout << "operaWavelengthCalibration: normalizationBinSize = " << normalizationBinSize << endl;
			cout << "operaWavelengthCalibration: LocalMaxFilterWidth = " << LocalMaxFilterWidth << endl;
			cout << "operaWavelengthCalibration: ParRangeSizeInPerCent = " << ParRangeSizeInPerCent << endl;
			cout << "operaWavelengthCalibration: NpointsPerPar = " << NpointsPerPar << endl;
			cout << "operaWavelengthCalibration: maxNIter = " << maxNIter << endl;
			cout << "operaWavelengthCalibration: minNumberOfLines = " << minNumberOfLines << endl;
			cout << "operaWavelengthCalibration: maxorderofpolynomial = " << maxorderofpolynomial << endl;
			cout << "operaWavelengthCalibration: dampingFactor = " << dampingFactor << endl;
			cout << "operaWavelengthCalibration: initialAcceptableMismatch = " << initialAcceptableMismatch << endl;
			cout << "operaWavelengthCalibration: nsigclip = " << nsigclip << endl;
			cout << "operaWavelengthCalibration: nOrdersToSearchAround = " << nOrdersToSearchAround << endl;
			cout << "operaWavelengthCalibration: referenceOrder = " << referenceOrder << endl;
			cout << "operaWavelengthCalibration: DetectionThreshold = " << DetectionThreshold << endl;
			cout << "operaWavelengthCalibration: MinPeakDepth = " << MinPeakDepth << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaWavelengthCalibration: ordernumber = " << ordernumber << endl;            
            }
            if(args.plot) {                
                cout << "operaWavelengthCalibration: ordersplotfilename = " << ordersplotfilename << endl;
                cout << "operaWavelengthCalibration: specplotfilename = " << specplotfilename << endl;
                cout << "operaWavelengthCalibration: ordersscriptfilename = " << ordersscriptfilename << endl;
                cout << "operaWavelengthCalibration: specscriptfilename = " << specscriptfilename << endl;
                cout << "operaWavelengthCalibration: ordersdatafilename = " << ordersdatafilename << endl;
                cout << "operaWavelengthCalibration: atlasdatafilename = " << atlasdatafilename << endl;
                cout << "operaWavelengthCalibration: compdatafilename = " << compdatafilename << endl;
                cout << "operaWavelengthCalibration: linesdatafilename = " << linesdatafilename << endl;
                cout << "operaWavelengthCalibration: generate3DPlot = " << generate3DPlot << endl;
                cout << "operaWavelengthCalibration: subtractCentralWavelength = " << subtractCentralWavelength << endl;
            }            
		}
		
		WavelengthCalibration::generate3DPlot = generate3DPlot;

		ofstream& fatlasdata = WavelengthCalibration::fatlasdata;
		ofstream& fcompdata = WavelengthCalibration::fcompdata;
		ofstream& flinesdata = WavelengthCalibration::flinesdata;
		ofstream& fordersdata = WavelengthCalibration::fordersdata;
		
		if (!atlasdatafilename.empty()) fatlasdata.open(atlasdatafilename.c_str());
		if (!compdatafilename.empty()) fcompdata.open(compdatafilename.c_str());
		if (!linesdatafilename.empty()) flinesdata.open(linesdatafilename.c_str());
		if (!ordersdatafilename.empty()) fordersdata.open(ordersdatafilename.c_str());
        
        FormatHeader specresheader("Spectral Resolution Data");
        specresheader << "number of orders" << newline << "order number" << "mean spectral resolution" << "spectral resolution dispersion" << "rv precision" << "spectral lines used" << newline;
        FormatData specresdata;
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, geometryfilename); // Read in the geometry
		
        // Read wavelength calibration initial guess
        if (!wlcal_initialguess.empty()) operaIOFormats::ReadIntoSpectralOrders(spectralOrders, wlcal_initialguess); // Read wavelength calibration reference first guess
        if (!uncalibrated_lines.empty()) operaIOFormats::ReadIntoSpectralOrders(spectralOrders, uncalibrated_lines); // Read in uncalibrated lines information
        if (!uncalibrated_spectrum.empty()) operaIOFormats::ReadIntoSpectralOrders(spectralOrders, uncalibrated_spectrum); // Read in uncalibrated spectrum information

        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
		if (args.verbose) cout << "operaWavelengthCalibration: minorder = " << minorder << " maxorder = " << maxorder << endl;
        
        // Read ThAr atlas spectrum - lambda vs. intensity, intensityVariance
		operaSpectrum atlasSpectrum;
		if (!atlas_spectrum.empty()) {      
			if (args.verbose) cout << "operaWavelengthCalibration: reading atlas spectrum " << atlas_spectrum << endl;
            atlasSpectrum = readAtlasSpectrum(atlas_spectrum);
        }
		
		// Read ThAr atlas lines - lambda vs. intensity
		operaSpectralLineList atlasLines;
		if (!atlas_lines.empty()) {         
			if (args.verbose) cout << "operaWavelengthCalibration: reading atlas lines " << atlas_lines << endl;
			atlasLines = readThoriumArgonAtlas(atlas_lines);
        }
        
		// Normalize uncalibrated spectrum
		if (normalizeUncalibratedSpectrum) {
			for (int order = minorder; order <= maxorder; order++) {
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements()) {
					spectralOrder->applyNormalizationForEmissionSpectrum(normalizationBinSize, 0, false, true, false);
				}
			}
		}

        // Initialize detection linewidth with uncalibrated_linewidth plus-minus 20% error.
        DetectionParameters detectionParameters = {uncalibrated_linewidth, uncalibrated_linewidth*0.2, LocalMaxFilterWidth, MinPeakDepth, DetectionThreshold};
        
        // Save original input polynomial solutions in case there has been any order shift
		WavelengthSolutions initialSolutions;
		initialSolutions.SetFromSpectralOrders(spectralOrders, minorder-nOrdersToSearchAround, maxorder+nOrdersToSearchAround);
        if (!inputLineSetFilename.empty()) {
            initialSolutions.CalculateFromLineSet(inputLineSetFilename, minorder-nOrdersToSearchAround, maxorder+nOrdersToSearchAround, maxorderofpolynomial);
        }
        
		// Investigate the first reference order with an existing initial solution to determine if there has been any order shift
		int referenceMinOrder = referenceOrder ? referenceOrder : minorder;
		int referenceMaxOrder = referenceOrder ? referenceOrder : maxorder;
		int selectedOrderShift = DetermineOrderShift(spectralOrders, referenceMinOrder, referenceMaxOrder, initialSolutions, nOrdersToSearchAround, atlasSpectrum, atlasLines, ParRangeSizeInPerCent, detectionParameters, nsigclip, NpointsPerPar, initialAcceptableMismatch, dampingFactor, minNumberOfLines, maxNIter);
        
        unsigned validorders = 0;
        ostringstream outputResolutionData;
        
        // Start the wavelength calibration performed on all orders.
        for (int order=minorder; order<=maxorder; order++) {
			
			DetectionParameters orderDetectionParameters = detectionParameters;
            
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			WavelengthCalibration orderCalibration(spectralOrder, atlasSpectrum, atlasLines);
			orderCalibration.SetFromInitialSolution(initialSolutions, order-selectedOrderShift);
			if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasGeometry() && spectralOrder->gethasWavelength()) {
				orderCalibration.CalculateWavelengthSolution(ParRangeSizeInPerCent, detectionParameters, nsigclip, parseSolution, NpointsPerPar, initialAcceptableMismatch, dampingFactor, minNumberOfLines, maxNIter);
                
                //CU Jul 16, 2015 - Commented out since this seems to throw away all the calculations that were already done...
                /*wavelength->createComparisonDataVectors(rawlinesinorder,rawlinecenter,rawlinecenterError,rawlineflux);
                wavelength->createAtlasDataVectors(atlaslinesinorder,atlasLineswl, atlasLineswlError,atlasLinesflux);
                wavelength->calculateSpectralResolution(ResolutionElementInPixels);
                wavelength->matchAtlaswithComparisonLines(ResolutionElementInPixels.value/2);*/
                operaWavelength* wavelength = spectralOrder->getWavelength();
                if (wavelength->getnDataPoints() > 1) {
					unsigned bestnpar = maxorderofpolynomial;
					if(wavelength->getnDataPoints() <= bestnpar) {
						bestnpar = wavelength->getnDataPoints()-1;
					}
                    //wavelength->filterDataPointsBySigmaClip((double)nsigclip/2); //CU Jul 16, 2015 - We no longer need to filter points out at this stage
					wavelength->CalculateWavelengthSolution(bestnpar,false);
                    
                    if (args.verbose) {
						printf("operaWavelengthCalibration: *****************************************************************\n");
						printf("operaWavelengthCalibration: Order %d: Final Solution:\n", order);
						orderCalibration.PrintSolution();
						printf("operaWavelengthCalibration: *****************************************************************\n");
					}
                    
                    orderCalibration.WritePlotOrdersData();
                    
					wavelength->calculateRadialVelocityPrecision();
					wavelength->calculateSpectralResolution();

					if (args.verbose){
						orderCalibration.PrintPrecisionAndResolution();
					}
                    
					if (args.debug)
						printf("%d %.2f %.2f %.2f %.2f %.2f\n", order, wavelength->getcentralWavelength(),wavelength->getResolutionElementInPixels().value, wavelength->getResolutionElementInPixels().error, wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
                    
                    if (!outputResolution.empty()) {
                        validorders++;
						outputResolutionData << order << " " << wavelength->getSpectralResolution().value << " " << wavelength->getSpectralResolution().error/2 << " " << wavelength->getRadialVelocityPrecision() << " " << wavelength->getnDataPoints() << endl;
					}
                }
                
				orderCalibration.WritePlotComparisonData(orderDetectionParameters.RawLineWidth);
                
                orderCalibration.WritePlotLinesData();
                
            } else if (!spectralOrder->gethasWavelength()) {
                if (args.verbose) {
                    printf("operaWavelengthCalibration: Order %d: has no associated wavelength reference calibration data. Wavelength calibration not possible.\n",  spectralOrder->getorder());
                }
            } else if (!spectralOrder->gethasGeometry()) {
                if (args.verbose) {
                    printf("operaWavelengthCalibration: Order %d: has no associated geometry data. Wavelength calibration not possible.\n",  spectralOrder->getorder());
                }
                spectralOrder->sethasWavelength(false);
            } else if (!spectralOrder->gethasSpectralElements()) {
                if (args.verbose) {
                    printf("operaWavelengthCalibration: Order %d: has no associated spectral elements data. Wavelength calibration not possible.\n",  spectralOrder->getorder());
                }
                spectralOrder->sethasWavelength(false);
            } else {
                if (args.verbose) {
                    printf("operaWavelengthCalibration: Order %d: has neither geometry nor wavelength reference calibration data. Wavelength calibration not possible.\n",  spectralOrder->getorder());
                }
                spectralOrder->sethasWavelength(false);
            }
		}
        
        // Wavelength Orders Info Plot: plot spectral resolution, rms in nm, radial velocity precision, data and polynomial solution
        if (fordersdata.is_open()) {
            fordersdata.close();
            if (!ordersscriptfilename.empty()) {
               GenerateWavelengthOrdersPlot(ordersscriptfilename, ordersplotfilename, ordersdatafilename, interactive);
            }
        }

        // Wavelength Spectrum Plot: plot atlas and comparison spectra and final set of matched lines.
        if(generate3DPlot) {
			if (fatlasdata.is_open() && fcompdata.is_open()) {
				fatlasdata.close();
				fcompdata.close();
                if (!specscriptfilename.empty()) {
                    GenerateWavelength3DSpecPlot(specscriptfilename, specplotfilename, atlasdatafilename, compdatafilename, subtractCentralWavelength, interactive);
                }
            }
        } else {
			if (fatlasdata.is_open() && fcompdata.is_open() && flinesdata.is_open()) {
				fatlasdata.close();
				fcompdata.close();
				flinesdata.close();
                
                if (!specscriptfilename.empty()) {
                    GenerateWavelengthSpecPlot(specscriptfilename, specplotfilename, atlasdatafilename, compdatafilename, linesdatafilename, subtractCentralWavelength, interactive);
                }
            }
        }
        // Write number of orders to spectral resolution file
        if (!outputResolution.empty()) {
            specresdata << validorders << endl << outputResolutionData.str();
            operaIOFormats::WriteCustomFormat("spectralres", specresheader, specresdata, outputResolution);
        }

		spectralOrders.setMinorder(minorder);
		spectralOrders.setMaxorder(maxorder);
		operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputWave, Wave);
	}
	catch (operaException e) {
		cerr << "operaWavelengthCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaWavelengthCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

// Read the entire thorium argon atlas and normalize the results
operaSpectralLineList readThoriumArgonAtlas(string atlas_lines) {
	operaSpectralLineList tharAtlas;
	igzstream astream(atlas_lines.c_str());
	if (astream.is_open()) {
		string dataline;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') {
				istringstream ss(dataline);
				double wn, wl, intensity;
				string marker;
				ss >> wn >> wl >> intensity >> marker;
				if (marker == "Th" || marker == "Ar") {
					tharAtlas.center.insert(wl * 0.1);
					tharAtlas.amplitude.insert(pow(10, intensity));
				} else {
					if (args.debug) printf("Skipping non-thorium-argon atlas entry %s.\n", marker.c_str());
				}
			}
		}
		unsigned lines = tharAtlas.size();
		if (args.verbose) {
			if (lines > 0) printf("          [Atlas] %d lines found wl0=%.2f wlc=%.2f wlf=%.2f\n", lines, tharAtlas.center[0], tharAtlas.center[lines / 2], tharAtlas.center[lines - 1]);
			else printf("          [Atlas] no lines found in atlas.\n");
		}
		astream.close();
	}
	return tharAtlas;
}

// Read the the full atlas spectrum
operaSpectrum readAtlasSpectrum(string atlas_spectrum) {
	operaSpectrum atlasSpectrum;
	igzstream astream(atlas_spectrum.c_str());
	if (astream.is_open()) {
		string dataline;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') {
				double tmpwl, tmpi, tmpvar, tmp1, tmp2;
				sscanf(dataline.c_str(), "%lf %lf %lf %lf %lf", &tmpwl, &tmpi, &tmp1, &tmp2, &tmpvar);
				atlasSpectrum.insert(0.1*tmpwl, tmpi, tmpvar);
            }
		}
		if (args.verbose) {
			if(atlasSpectrum.size() > 0) printf("          [Atlas] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", atlasSpectrum.size(), atlasSpectrum.firstwl(), atlasSpectrum.midwl(), atlasSpectrum.lastwl());
			else printf("          [Atlas] no points found in atlas.\n");
		}
		astream.close();
	}
	return atlasSpectrum;
}

int DetermineOrderShift(operaSpectralOrderVector& spectralOrders, int referenceMinOrder, int referenceMaxOrder, const WavelengthSolutions& initialSolutions, int nOrdersToSearchAround, const operaSpectrum& atlasSpectrum, const operaSpectralLineList& atlasLines, double ParRangeSizeInPerCent, DetectionParameters detectionParams, double nsigclip, unsigned NpointsPerPar, double initialAcceptableMismatch, double dampingFactor, unsigned minNumberOfLines, unsigned maxNIter) {
	if (nOrdersToSearchAround == 0) return 0;

	WavelengthCalibration::skipPlots = true;
	int selectedOrderShift = 0;
	double minRV = +BIG;
    
    // Find the first available reference order with initial solution, and is valid when shifted, to avoid looping over all orders
    int referenceOrderUsed = 0;
    for (int order = referenceMinOrder; order <= referenceMaxOrder; order++) {
        if (initialSolutions.HasSolutionForOrder(order)) {
			referenceOrderUsed = order;
            // Make sure if each shift in search range is valid
            for (int o = -nOrdersToSearchAround; o <= nOrdersToSearchAround; o++) {
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(referenceOrderUsed + o);
                if ((!spectralOrder->gethasSpectralElements() && !spectralOrder->gethasSpectralLines()) || !spectralOrder->gethasGeometry()) {
                    referenceOrderUsed = 0;
                    break;
                }
            }
            if (referenceOrderUsed) break;
        }
    }
    
	if(referenceOrderUsed) {
	    // Test each order shift to find the lowest RV error
        for (int o = -nOrdersToSearchAround; o <= nOrdersToSearchAround; o++) {
            // Calculate a wavelength solution for the shifted reference order and determine if the RV error is lower with this shift
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(referenceOrderUsed + o);
            WavelengthCalibration ordercal(spectralOrder, atlasSpectrum, atlasLines);
	        ordercal.SetFromInitialSolution(initialSolutions, referenceOrderUsed);
	        
	        ofstream closedstream;
	        ordercal.CalculateWavelengthSolution(ParRangeSizeInPerCent, detectionParams, nsigclip, false, NpointsPerPar, initialAcceptableMismatch, dampingFactor, minNumberOfLines, maxNIter);
                
            operaWavelength* wavelength = spectralOrder->getWavelength();
            if(wavelength->getnDataPoints() > 0) {
                if(wavelength->getRadialVelocityPrecision() < minRV) {
                    minRV = wavelength->getRadialVelocityPrecision();
                    selectedOrderShift = o;
                }
                if(args.verbose) cout << "operaWavelengthCalibration: Order " << referenceOrderUsed << " minRV=" << minRV << " currRV=" << wavelength->getRadialVelocityPrecision() << " ordershift=" << o << endl;
            }
        }
    }
    if(args.verbose) cout << "operaWavelengthCalibration: Determined ordershift=" << selectedOrderShift << " using referenceOrder=" << referenceOrderUsed << endl;
	WavelengthCalibration::skipPlots = false;
    return selectedOrderShift;
}

void WavelengthSolutions::SetFromSpectralOrders(const operaSpectralOrderVector& spectralOrders, int minorder, int maxorder) {
	for (int order = minorder; order <= maxorder; order++) {
		const operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
		if (spectralOrder->gethasWavelength()) {
			polynomials[order] = *spectralOrder->getWavelength()->getWavelengthPolynomial();
		}
	}
}

void WavelengthSolutions::CalculateFromLineSet(string inputLineSetFilename, int minorder, int maxorder, unsigned maxorderofpolynomial) {
	for (int order = minorder; order <= maxorder; order++) {
        operaVector wavelengthData, distanceData;
		readLineSet(inputLineSetFilename, order, wavelengthData, distanceData);
		operaVector wavelengthErrors(wavelengthData.size());

		if (wavelengthData.size() > 0) {
			operaWavelength wavelength(maxorderofpolynomial);
			wavelength.createDataVectors(wavelengthData, wavelengthErrors, distanceData);
			wavelength.CalculateWavelengthSolution(maxorderofpolynomial, false);

			if (args.debug) {
				wavelength.getWavelengthPolynomial()->printEquation(&cout);
				cout << "order " << order << " chisqr=" << wavelength.getWavelengthPolynomial()->getChisqr() << endl;
			}
			polynomials[order] = *wavelength.getWavelengthPolynomial();
		}
	}
}

bool WavelengthSolutions::HasSolutionForOrder(int order) const {
	return polynomials.count(order) > 0;
}

const Polynomial& WavelengthSolutions::GetSolutionForOrder(int order) const {
	return polynomials.find(order)->second;
}

void WavelengthSolutions::readLineSet(string inputLineSetFilename, int order, operaVector& wavelength, operaVector& distance) const {
	igzstream astream(inputLineSetFilename.c_str());
	if (astream.is_open()) {
		string dataline;
		unsigned np = 0;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') { // skip blank lines and comments
				int tmpo = 0;
				double tmpwl = -1.0;
				double tmpdist = -1.0;
				sscanf(dataline.c_str(), "%d %lf %lf", &tmpo, &tmpwl, &tmpdist);
				if (tmpo == order) {
					wavelength.insert(tmpwl);
					distance.insert(tmpdist);
					np++;
				}
				else if (tmpo != order && np) {
					break;
				}
			}
		}
		astream.close();
	}
}


WavelengthCalibration::WavelengthCalibration(operaSpectralOrder* spectralOrder, const operaSpectrum& atlasSpectrum, const operaSpectralLineList& atlasLines)
: spectralOrder(spectralOrder), atlasSpectrumFull(atlasSpectrum), atlasLinesFull(atlasLines) {
	order = spectralOrder->getorder();
	if(!spectralOrder->gethasWavelength()) spectralOrder->createWavelength(0);
	wavelength = spectralOrder->getWavelength();
}

void WavelengthCalibration::SetFromInitialSolution(const WavelengthSolutions& initialSolutions, int order) {
	if (initialSolutions.HasSolutionForOrder(order)) {
		Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
		*wavelengthPolynomial = initialSolutions.GetSolutionForOrder(order);
		spectralOrder->sethasWavelength(true);
	}
	else {
		spectralOrder->sethasWavelength(false);
	}
}

void WavelengthCalibration::CalculateWavelengthSolution(double ParRangeSizeInPerCent, DetectionParameters detectionParams, double nsigclip, bool parseSolution, unsigned NpointsPerPar, double initialAcceptableMismatch, double dampingFactor, unsigned minNumberOfLines, unsigned maxNIter) {
	InitializeDistanceLimitsFromGeometry();

	FindAndSetComparisonAndAtlasLines(detectionParams, nsigclip);
	if (!hasComparisonLines || !hasAtlasLines) return;
	
	// At this point all lines in both comparison lines (distance) and the Atlas lines (wavelength) have been set. Start identifying matching lines and refining wavelength solution.
	doubleValue_t ResolutionElementInPixels = {2.0*comparisonLineWidth, comparisonLineWidthErr};
	wavelength->setResolutionElementInPixels(ResolutionElementInPixels);
	wavelength->calculateSpectralResolution();
	if (parseSolution) {
		// Vary all 3 coefficients to search for the 2nd degree (parabola) wavelength solution with maximum correlation between simulated spectra for the comparison and atlas.
		wavelength->refineWavelengthSolutionOfSecondOrderByXCorrelation(NpointsPerPar, ParRangeSizeInPerCent);
	} else {
		// Vary the zeroth-order coefficient and check for matching lines in order to find the highest matching rate.
		// This allows the spectral lines to be identified when there is a shift in the detector position of the observed spectrum.
		wavelength->refineWavelengthSolutionByFindingMaxMatching(NpointsPerPar, ParRangeSizeInPerCent, initialAcceptableMismatch);
	}

	// Find matching atlas and comparison lines, and perform a polynomial fit using these points
	FitSolutionPolynomialUsingMatchingLines(initialAcceptableMismatch);

	// Iterate to improve the polynomial that gives the pixel-to-wavelength solution
	RefineSolutionPolynomialFit(initialAcceptableMismatch, dampingFactor, minNumberOfLines, nsigclip, maxNIter);

	// At this point there should be a good wavelength solution
	FinishedShrinkingPolynomialFit();
}

void WavelengthCalibration::InitializeDistanceLimitsFromGeometry() {
	operaGeometry *geometry = spectralOrder->getGeometry();

	double dmin = 0.0;
	double dmax = geometry->CalculateAndSetOrderLength();
	wavelength->setDmin(dmin);
	wavelength->setDmax(dmax);

	if (args.verbose) {
		printf("operaWavelengthCalibration: Order %d: [geom] ymin = %.2f ymax = %.2f dmin = %.2f dmax = %.2f \n", order, geometry->getYmin(), geometry->getYmax(), dmin, dmax);
	}
}

void WavelengthCalibration::FindAndSetComparisonAndAtlasLines(DetectionParameters detectionParams, double nsigclip) {
	// Calculate the initial and final wavelength based on the geometry calibration. This will be used to select the atlas range.
	double wl0 = wavelength->getinitialWavelength();
	double wlf = wavelength->getfinalWavelength();
	wl_central = wavelength->getcentralWavelength();
	operaWavelengthRange wlrange(wl0, wlf);
	if (args.verbose) {
		printf("operaWavelengthCalibration: Order %d: [wave] wavelength selected range: wl0 = %.2f wlc = %.2f wlf = %.2f\n", order, wl0, wl_central, wlf);
	}

	// Load spectral lines for our comparison. If the order has a spectrum, we use it to detect lines. Otherwise we try using the existing lines in the order.
	if (HasUncalSpectrum()) {
		SetComparisonLines(GetLinesFromUncalSpectrum(detectionParams, nsigclip));
	}
	else if (HasUncalLines()) {
		SetComparisonLines(GetUncalLines());
	}
	if (!hasComparisonLines) return;
	
    // Load spectral lines for our Atlas. If the atlasSpectrum is loaded, we use it to detect lines, otherwise we try using the atlasLines inside wlrange.
	if (atlasSpectrumFull.size() > 0) {
		detectionParams.RawLineWidth = comparisonLineWidth; // Use the width of the comparison lines for detecting atlas lines.
        SetAtlasLines(GetLinesFromAtlasSpectrum(wlrange, detectionParams));
	}
	else if (atlasLinesFull.size() > 0) {
		SetAtlasLines(GetAtlasLinesInRange(wlrange, comparisonLineWidth));
	}
}

bool WavelengthCalibration::HasUncalLines() const {
	return spectralOrder->gethasSpectralLines();
}

bool WavelengthCalibration::HasUncalSpectrum() const {
	return spectralOrder->gethasSpectralElements();
}

operaSpectralLineList WavelengthCalibration::GetUncalLines() {
	if (args.verbose) cout << "operaWavelengthCalibration: using uncalibrated lines" << endl;
	operaSpectralLineList rawLines = getAllSpectralLines(*spectralOrder->getSpectralLines());
	return rawLines;
}

operaSpectralLineList WavelengthCalibration::GetLinesFromUncalSpectrum(DetectionParameters detection, double nsigclip) {
	if (args.verbose) cout << "operaWavelengthCalibration: detecting lines in uncalibrated spectrum" << endl;
	operaSpectralLineList rawLines = spectralOrder->getRawLinesFromUncalSpectrum(detection.RawLineWidth, detection.LocalMaxFilterWidth, detection.MinPeakDepth, detection.DetectionThreshold, nsigclip, distance_disp);
	return rawLines;
}

operaSpectralLineList WavelengthCalibration::GetAtlasLinesInRange(operaWavelengthRange wlrange, double rawLineWidth) {
	if (args.verbose) cout << "operaWavelengthCalibration: using atlas lines " << endl;
	operaSpectralLineList atlasLines = getSpectrumWithinRange(wlrange, atlasLinesFull);
	atlasLines.centerError = wavelength->convertPixelToWavelength(rawLineWidth);
	return atlasLines;
}

operaSpectralLineList WavelengthCalibration::GetLinesFromAtlasSpectrum(operaWavelengthRange wlrange, DetectionParameters detection) {
	if (args.verbose) cout << "operaWavelengthCalibration: detecting lines in atlas spectrum " << endl;
	detection.RawLineWidth = wavelength->convertPixelToWavelength(detection.RawLineWidth);
	operaSpectralLineList atlasLines = DetectAtlasLines(wlrange, detection);
	return atlasLines;
}

void WavelengthCalibration::SetComparisonLines(const operaSpectralLineList& comparisonLines) {
	if (comparisonLines.size() > 0) {
		if (args.verbose) {
			printf("operaWavelengthCalibration: Order %d: [Comparison] %d lines in comparison between wl0 = %.2f and wlf = %.2f.\n", order, comparisonLines.size(), wavelength->evaluateWavelength(comparisonLines.center[0]), wavelength->evaluateWavelength(comparisonLines.center[comparisonLines.size() - 1]));
		}
		wavelength->setComparisonDataVectors(comparisonLines);
		comparisonLineWidth = comparisonLines.medianWidth;
		comparisonLineWidthErr = comparisonLines.medianWidthError;
		hasComparisonLines = true;
	}
	else {
		if (args.verbose) {
			printf("operaWavelengthCalibration: Warning: Order %d: [Comparison] No lines detected from input comparison. Skipping calibration.\n", order);
		}
		hasComparisonLines = false;
	}
}

void WavelengthCalibration::SetAtlasLines(const operaSpectralLineList& atlasLines) {
	if (atlasLines.size() > 0) {
		if (args.verbose) {
			printf("operaWavelengthCalibration: Order %d: [Atlas] %d lines detected in input atlas between wl0 = %.2f and wlf = %.2f .\n", order, atlasLines.size(), atlasLines.center[0], atlasLines.center[atlasLines.size()-1]);
		}
		wavelength->setAtlasDataVectors(atlasLines);
		hasAtlasLines = true;
	}
	else {
		printf("operaWavelengthCalibration: Warning:  Order %d: [Atlas] No lines detected from input atlas. Skipping calibration.\n", order);
		hasAtlasLines = false;
	}
}

operaSpectralLineList WavelengthCalibration::DetectAtlasLines(operaWavelengthRange wlrange, DetectionParameters detection) {
	operaSpectrum atlasRegion = getSpectrumWithinRange(wlrange, atlasSpectrumFull);

	WritePlotAtlasData(atlasRegion);

	// Calculate the cross-correlation between the atlas spectrum and a gaussian function.
	if (args.verbose) cout << "operaWavelengthCalibration: Calculating cross correlation with Gaussian.." << endl;
	operaVector atlasXcorr = calculateXCorrWithGaussian(atlasRegion.wavelengthvector(), atlasRegion.fluxvector(), detection.RawLineWidth);

	// Degrade the resolution of the atlas to the resolution of raw lines. The degradation is done by convolving the spectrum with a gaussian.
	if (args.verbose) cout << "operaWavelengthCalibration: Convolving specctrum with Gaussian.." << endl;
	operaVector convolvedAtlas = convolveSpectrumWithGaussian(atlasRegion.wavelengthvector(), atlasRegion.fluxvector(), detection.RawLineWidth);

	// Below it reads the atlas spectrum into an operaSpectralElements class
	operaSpectralElements atlasElements(atlasRegion.size());
	atlasElements.setXCorrelation(atlasXcorr);
	atlasElements.setWavelength(atlasRegion.wavelengthvector());
	atlasElements.setFluxVector(operaFluxVector(convolvedAtlas, atlasRegion.variancevector()));
	atlasElements.setHasXCorrelation(true);
	atlasElements.setHasWavelength(true);
	atlasElements.setHasRawSpectrum(true);

	if (args.debug) {
		for (unsigned i = 0; i<atlasElements.getnSpectralElements(); i++) {
			cout << atlasElements.getwavelength(i) << " " << atlasElements.getFlux(i) << " " << atlasElements.getFluxVariance(i) << " " << atlasElements.getXCorrelation(i) << endl;
		}
	}

	// Detect operaSpectralLines in the atlas
	operaSpectralLines atlasLines = DetectSpectralLines(&atlasElements, detection.RawLineWidth, detection.LocalMaxFilterWidth, detection.MinPeakDepth, detection.DetectionThreshold, wavelength_disp);

	// Get the individual atlas lines in the wavelength range
	return getSpectralLinesInWavelengthRange(atlasLines, wlrange);
}

void WavelengthCalibration::FitSolutionPolynomialUsingMatchingLines(double acceptableMismatch) {
	// Recalculate the matching elements using the current wavelength solution.
	wavelength->calculateSpectralResolution();
	wavelength->matchAtlaswithComparisonLines(acceptableMismatch);

	// List wavelength and distance for matching lines
	if (args.debug) {
		for (unsigned index = 0; index<wavelength->getnDataPoints(); index++) {
			cout << order << " " << wavelength->getWavelength(index) << " " << wavelength->getDistance(index) << endl;
		}
	}

	Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
	unsigned bestnpar = wavelengthPolynomial->getOrderOfPolynomial();

	if (wavelength->getnDataPoints() > 0) {
		if (wavelength->getnDataPoints() <= bestnpar) {
			bestnpar = wavelength->getnDataPoints() - 1;
		}
		wavelength->CalculateWavelengthSolution(bestnpar, false);

		if (args.verbose) {
			printf("\noperaWavelengthCalibration: Order %d: Initial Solution:\n", order);
			PrintSolution();
		}
	}
	else {
		printf("operaWavelengthCalibration: Order %d: ZERO points to calculate wavelength solution, skipping order...\n", order);
	}
}

void WavelengthCalibration::RefineSolutionPolynomialFit(double acceptableMismatch, double dampingFactor, unsigned minNumberOfLines, double nsigclip, unsigned maxNIter) {
	if (wavelength->getnDataPoints() == 0) return;
	double minchisqr = wavelength->getWavelengthPolynomial()->getChisqr();
	unsigned nochangeinChisqr = 0;
	unsigned bestnpar;
	for (unsigned iter = 0; iter < maxNIter; iter++) {
		if (args.debug) printf("operaWavelengthCalibration: Refinement iteration %d\n", iter);

		bestnpar = wavelength->getWavelengthPolynomial()->getOrderOfPolynomial();

		// The acceptable mismatch can start considerably big and then shrink down as calibration gets better.
		// However, it will only shrink either to a minimum value or minimum number of lines.
		if (wavelength->getnDataPoints() > minNumberOfLines) acceptableMismatch *= dampingFactor;

		// Recalculate the spectral resolution so we know the correct line width.
		wavelength->calculateSpectralResolution();

		// Set the wavelengths data points to matching atlas and comparison lines.
		wavelength->matchAtlaswithComparisonLines(acceptableMismatch);

		// Filter out points where |wl - p(dist)| >= nsigclip * rms(wl - p(dist)
		wavelength->filterDataPointsBySigmaClip(nsigclip);
		if (wavelength->getnDataPoints() == 0) {
			if (args.verbose) cout << "operaWavelengthCalibration: No data points left, stopping refinement" << endl;
			break;
		}

		// Perform a new polynomial fit (of order bestnpar) and keep it only if the chi-square is lowered. Don't use polynomial errors.
		wavelength->RefineWavelengthSolution(bestnpar, false);

		// It's pointless to use a polynomial of higher order than the number of remaining points.
		if (wavelength->getnDataPoints() <= bestnpar) {
			bestnpar = wavelength->getnDataPoints() - 1;
		}

		// Try to fit polynomials of all orders up to bestnpar, and keep the one with the lowest chi-square. Don't use polynomial errors.
		wavelength->CalculateWavelengthSolution(bestnpar, false);

		if (args.debug) {
			printf("\noperaWavelengthCalibration: Order %d: Revised Solution:\n", order);
			PrintSolution();
		}
		
		// Exit the loop if the chi-sqaure doesn't drop for 3 consecutive iterations.
		const Polynomial* wavelengthPolynomial = wavelength->getWavelengthPolynomial();
		if (wavelengthPolynomial->getChisqr() < minchisqr) {
			minchisqr = wavelengthPolynomial->getChisqr();
			nochangeinChisqr = 0;
		}
		else if (wavelengthPolynomial->getChisqr() == minchisqr) {
			nochangeinChisqr++;
			// The line below is commented out to reproduce a bug in a previous version, where this condition was never met, which gave better results.
			// The intended result might be able to be acheieved (without greatly reducing precision) by raising the value compared against from 3.
			//if (nochangeinChisqr > 3) break;
		}
	}
}

void WavelengthCalibration::FinishedShrinkingPolynomialFit() {
	if (wavelength->getnDataPoints() > 0) {
		if (args.debug) {
			for (unsigned l = 0; l<wavelength->getnDataPoints(); l++) {
				cout << order << " " << wavelength->getDistance(l) << " " << wavelength->getWavelength(l) << " " << wavelength->evaluateWavelength(wavelength->getDistance(l)) << " " << wavelength->getWavelength(l) - wavelength->evaluateWavelength(wavelength->getDistance(l)) << " " << wavelength->getWavelengthError(l) << endl;
			}
		}

		if (args.verbose) {
			printf("\noperaWavelengthCalibration: Order %d: Wavelength solution after done shrinking:\n", order);
			PrintSolution();
		}

		// Calculate the radial velocity precision and spectral resolution
		wavelength->calculateRadialVelocityPrecision();
		wavelength->calculateSpectralResolution();

		if (args.verbose){
			PrintPrecisionAndResolution();
		}

		if (args.debug) {
			cout << order << " "
				<< wavelength->getnDataPoints() << " "
				<< wavelength->getinitialWavelength() << " "
				<< wavelength->getcentralWavelength() << " "
				<< wavelength->getfinalWavelength() << " "
				<< wavelength->calculateWavelengthRMSPrecision() << " "
				<< wavelength->calculateWavelengthMedianPrecision() << " "
				<< wavelength->getRadialVelocityPrecision() << " "
				<< wavelength->getResolutionElementInPixels().value << " "
				<< wavelength->getResolutionElementInPixels().error << " "
				<< wavelength->getSpectralResolution().value << " "
				<< wavelength->getSpectralResolution().error << endl;
		}
	}
}

void WavelengthCalibration::PrintSolution() {
	double ComparisonMatchPercentage = wavelength->getPerCentageOfComparisonMatch();
	double AtlasMatchPercentage = wavelength->getPerCentageOfAtlasMatch();
	printf("operaWavelengthCalibration: Order %d: ", order);
	wavelength->getWavelengthPolynomial()->printEquation(&cout);
	cout << " chisqr=" << wavelength->getWavelengthPolynomial()->getChisqr() << endl;
	printf("operaWavelengthCalibration: Order %d: %u lines matched between wl0 = %.2f  wlf = %.2f.\n", order, wavelength->getnDataPoints(), wavelength->getinitialWavelength(), wavelength->getfinalWavelength());
	printf("operaWavelengthCalibration: Order %d: [Atlas]    matched %.2f %% of detected lines.\n", order, AtlasMatchPercentage);
	printf("operaWavelengthCalibration: Order %d: [Comparison] matched %.2f %% of detected lines.\n", order, ComparisonMatchPercentage);
}

void WavelengthCalibration::PrintPrecisionAndResolution() {
	printf("operaWavelengthCalibration: -----------------------------------------------------------------\n");
	printf("operaWavelengthCalibration: Order %d: Radial velocity precision = %.2f m/s. Using %d spectral lines.\n", order, wavelength->getRadialVelocityPrecision(), wavelength->getnDataPoints());
	printf("operaWavelengthCalibration: Order %d: Wavelength RMS precision = %.10f nm  Median Precision = %.10f nm.\n", order, wavelength->calculateWavelengthRMSPrecision(), wavelength->calculateWavelengthMedianPrecision());
	printf("operaWavelengthCalibration: Order %d: [Comparison Lines] median sigma = %.2f +/- %.2f.\n", order, comparisonLineWidth, comparisonLineWidthErr);
	printf("operaWavelengthCalibration: Order %d: Spectral Resolution = %.2f +/- %.2f.\n", order, wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
	printf("operaWavelengthCalibration: -----------------------------------------------------------------\n\n");
}

void WavelengthCalibration::WritePlotAtlasData(const operaSpectrum& atlasRegion) {
	if (fatlasdata.is_open() && !skipPlots) {
		double maxatlasflux = Max(atlasRegion.fluxvector());
		if (generate3DPlot) {
			// Below it produces data for a 3D plot of spectrum.
			for (unsigned slitview = 0; slitview < 2; slitview++){
				for (unsigned i = 0; i < atlasRegion.size(); i++) {
					fatlasdata << order << " " << atlasRegion.getwavelength(i) << " " << atlasRegion.getflux(i) / maxatlasflux << " " << slitview << " " << wl_central << endl;
				}
				fatlasdata << endl;
			}
			fatlasdata << endl;
		}
		else {
			// Below it produces data for a 2D plot of spectrum.
			for (unsigned i = 0; i < atlasRegion.size(); i++) {
				fatlasdata << order << " " << atlasRegion.getwavelength(i) << " " << atlasRegion.getflux(i) / maxatlasflux << " " << wl_central << endl;
			}
			fatlasdata << endl;
		}
	}
}

void WavelengthCalibration::WritePlotComparisonData(double rawlinewidth) {
	if (fcompdata.is_open() && !skipPlots) {
		operaSpectralElements* spectralElements = spectralOrder->getSpectralElements();
		double maxflux = Max(spectralElements->getFluxVector().getflux());
		if (generate3DPlot) {
			for (unsigned slitview = 0; slitview < 2; slitview++){
				unsigned lastline = 0;
				for (unsigned i = 0; i < spectralElements->getnSpectralElements(); i++) {
					double dist = spectralElements->getdistd(i);
					double wl = wavelength->evaluateWavelength(dist);
					double flux = spectralElements->getFlux(i);
					double matchlinesflux = 0;

					for (unsigned lines = lastline; lines<wavelength->getnDataPoints(); lines++){
						if (dist > wavelength->getDistance(lines) - rawlinewidth && dist < wavelength->getDistance(lines) + rawlinewidth) {
							matchlinesflux = 1;
							lastline = lines;
							break;
						}
						else if (dist > wavelength->getDistance(lines) + rawlinewidth) {
							lastline++;
							break;
						}
					}
					fcompdata << order << setprecision(8) << " " << dist << " " << wl << " " << flux / maxflux << " " << matchlinesflux << " " << slitview << " " << endl;
				}
				fcompdata << endl;
			}
		}
		else {
			fcompdata << "order dist wl flux centralwl" << endl;
			for (unsigned i = 0; i < spectralElements->getnSpectralElements(); i++) {
				double dist = spectralElements->getdistd(i);
				double wl = wavelength->evaluateWavelength(dist);
				double flux = spectralElements->getFlux(i);
				fcompdata << order << setprecision(8) << " " << fixed << dist << " " << fixed << wl << " " << scientific << flux / maxflux << " " << fixed << wl_central << endl;
			}
		}
		fcompdata << endl;
	}
}

void WavelengthCalibration::WritePlotLinesData() {
	if (flinesdata.is_open() && !skipPlots) {
		operaSpectralElements* spectralElements = spectralOrder->getSpectralElements();
		double maxflux = Max(spectralElements->getFluxVector().getflux());
		flinesdata << "order dist wavelength comparisonflux resolution centralwl comparisonwl atlaswl atlasflux" << endl;
		for (unsigned lines = 0; lines < wavelength->getnDataPoints(); lines++){
			flinesdata << order << " " << setprecision(8);
			flinesdata << fixed << wavelength->getDistance(lines) << " " << wavelength->getWavelength(lines) << " ";
			flinesdata << scientific << (wavelength->getcomparisonLinesflux(wavelength->getMatchComparisonIndex(lines)) / maxflux) / 2.0 << " ";
			flinesdata << fixed << wavelength->convertPixelToWavelength(wavelength->getResolutionElementInPixels().value) << " " << wl_central << " " << wavelength->getcomparisonLineswl(wavelength->getMatchComparisonIndex(lines)) << " " << wavelength->getatlasLineswl(wavelength->getMatchAtlasIndex(lines)) << " ";
			flinesdata << scientific << (wavelength->getatlasLinesflux(wavelength->getMatchAtlasIndex(lines)) / maxflux) / 2.0 << endl;
		}
		flinesdata << endl;
	}
}

void WavelengthCalibration::WritePlotOrdersData() {
	if (fordersdata.is_open() && !skipPlots) {
		double oldprecision = wavelength->getRadialVelocityPrecision();
		doubleValue_t oldresolution = wavelength->getSpectralResolution();
		wavelength->calculateRadialVelocityPrecision();
		wavelength->calculateSpectralResolution();
		fordersdata << order << " "
			<< wavelength->getinitialWavelength() << " "
			<< wavelength->getcentralWavelength() << " "
			<< wavelength->getfinalWavelength() << " "
			<< wavelength->calculateWavelengthRMSPrecision() << " "
			<< wavelength->calculateWavelengthMedianPrecision() << " "
			<< wavelength->getResolutionElementInPixels().value << " "
			<< wavelength->getResolutionElementInPixels().error << " "
			<< wavelength->getnDataPoints() << " "
			<< wavelength->getPerCentageOfComparisonMatch() << " "
			<< wavelength->getPerCentageOfAtlasMatch() << " "
			<< oldprecision << " "
			<< wavelength->getRadialVelocityPrecision() << " "
			<< wavelength->getSpectralResolution().value << " "
			<< wavelength->getSpectralResolution().error / 2 << endl;
		wavelength->setRadialVelocityPrecision(oldprecision);
		wavelength->setSpectralResolution(oldresolution);
	}
}

// Generates multiple plots containing statistical info about wavelength calibration
void GenerateWavelengthOrdersPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display)
{
	if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());

	fgnu << "reset" << endl;
	fgnu << "unset key\n" << endl;
	fgnu << "NX=2; NY=2" << endl;
	fgnu << "DX=0.1; DY=0.1; SX=0.42; SY=0.42" << endl;
	fgnu << "set bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY" << endl;
	fgnu << "set size SX*NX+DX*2,SY*NY+DY*2" << endl;

	if (!outputPlotEPSFileName.empty()) {
		fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
		fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;

		fgnu << "set multiplot" << endl;

		fgnu << "set size 0.9*SX,0.9*SY" << endl;
		fgnu << "unset xlabel" << endl;
		fgnu << "unset y2tics; set ytics" << endl;
		fgnu << "unset y2tics; set ylabel \"{/Symbol l} precision (nm)\"" << endl;
		fgnu << "set origin 0.75*DX,DY+SY" << endl;
		fgnu << "set key" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:5 t \"rms of residuals\" w linespoint lw 1.5, \"\" u 1:6  t \"median of residuals\" w linespoint lw 1.5" << endl;

		fgnu << "unset key" << endl;
		fgnu << "unset xlabel" << endl;
		fgnu << "unset ytics; set y2tics mirror" << endl;
		fgnu << "unset ylabel; set y2label \"spectral resolution {/Symbol l}/{/Symbol Dl}\"" << endl;
		fgnu << "set origin DX+SX,DY+SY" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:14:15 w yerr pt 7 lw 1.5 lt 1,\"\" u 1:14 w linespoint lt 1" << endl;

		fgnu << "set key" << endl;
		fgnu << "unset x2label; set xlabel \"order number\"" << endl;
		fgnu << "unset y2tics; set ytics" << endl;
		fgnu << "unset y2label; set ylabel \"Radial velocity precision (m/s)\"" << endl;
		fgnu << "set origin 0.75*DX,DY" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:13 t \"full set of lines\" w linespoint pt 7 lt 3 lw 2, \"\" u 1:12 t \"clean sample\" w linespoint pt 7 lt 4 lw 2" << endl;
		fgnu << "unset key" << endl;

		fgnu << "unset x2label; set xlabel \"order number\"" << endl;
		fgnu << "unset y2tics; unset y2label" << endl;
		fgnu << "set origin DX+SX,DY" << endl;
		fgnu << "set ytics nomirror" << endl;
		fgnu << "set ylabel \"% matching lines\" offset +1.5,0l" << endl;
		fgnu << "set key bottom" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:10 t \"% of comparison lines\" w linespoint pt 7 lt 4 lw 1, \"\" u 1:11 t \"% of atlas lines\" w linespoint pt 6 lt 4 lw 1" << endl;

		fgnu << "unset key" << endl;
		fgnu << "unset ytics; set y2tics" << endl;
		fgnu << "unset ylabel; set y2label \"Number of matching lines\"" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:9 w linespoint lt 2 lw 3" << endl;

		fgnu << "unset multiplot" << endl;

		if (display) {
			fgnu << "\nset terminal x11" << endl;
			fgnu << "set output" << endl;
			fgnu << "replot" << endl;
		}
		else {
			fgnu << "\n#set terminal x11" << endl;
			fgnu << "#set output" << endl;
			fgnu << "#replot" << endl;
		}
	}
	else {
		fgnu << "set multiplot" << endl;

		fgnu << "set size 0.9*SX,0.9*SY" << endl;
		fgnu << "unset xlabel" << endl;
		fgnu << "unset y2tics; set ytics" << endl;
		fgnu << "unset y2tics; set ylabel \"{/Symbol l} precision (nm)\"" << endl;
		fgnu << "set origin 0.75*DX,DY+SY" << endl;
		fgnu << "set key" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:5 t \"rms of residuals\" w linespoint lw 1.5, \"\" u 1:6  t \"median of residuals\" w linespoint lw 1.5" << endl;

		fgnu << "unset key" << endl;
		fgnu << "unset xlabel" << endl;
		fgnu << "unset ytics; set y2tics mirror" << endl;
		fgnu << "unset ylabel; set y2label \"spectral resolution {/Symbol l}/{/Symbol Dl}\"" << endl;
		fgnu << "set origin DX+SX,DY+SY" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:14:15 w yerr pt 7 lw 1.5 lt 1,\"\" u 1:14 w linespoint lt 1" << endl;

		fgnu << "set key" << endl;
		fgnu << "unset x2label; set xlabel \"order number\"" << endl;
		fgnu << "unset y2tics; set ytics" << endl;
		fgnu << "unset y2label; set ylabel \"Radial velocity precision (m/s)\"" << endl;
		fgnu << "set origin 0.75*DX,DY" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:13 t \"full set of lines\" w linespoint pt 7 lt 3 lw 2, \"\" u 1:12 t \"clean sample\" w linespoint pt 7 lt 4 lw 2" << endl;
		fgnu << "unset key" << endl;

		fgnu << "unset x2label; set xlabel \"order number\"" << endl;
		fgnu << "unset y2tics; unset y2label" << endl;
		fgnu << "set origin DX+SX,DY" << endl;
		fgnu << "set ytics nomirror" << endl;
		fgnu << "set ylabel \"% matching lines\" offset +1.5,0l" << endl;
		fgnu << "set key bottom" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:10 t \"% of comparison lines\" w linespoint pt 7 lt 4 lw 1, \"\" u 1:11 t \"% of atlas lines\" w linespoint pt 6 lt 4 lw 1" << endl;

		fgnu << "unset key" << endl;
		fgnu << "unset ytics; set y2tics" << endl;
		fgnu << "unset ylabel; set y2label \"Number of matching lines\"" << endl;
		fgnu << "plot \"" << dataFileName << "\" u 1:9 w linespoint lt 2 lw 3" << endl;

		fgnu << endl;

		fgnu << "unset multiplot" << endl;

		fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
		fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
		fgnu << "#replot" << endl;
		fgnu << "#set terminal x11" << endl;
		fgnu << "#set output" << endl;
	}

	fgnu.close();

	if (display) systemf("gnuplot -persist %s", gnuScriptFileName.c_str());
	else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s", gnuScriptFileName.c_str());
}

// Generates 2D plot for spectra of atlas + comparison + identified lines
void GenerateWavelengthSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, string linesdatafilename, bool subtractCentralWavelength, bool display)
{
	if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());

	fgnu << "reset" << endl;
	fgnu << "unset key" << endl;
	if (subtractCentralWavelength) {
		fgnu << "\nset xlabel \"{/Symbol l} - {/Symbol l}_c (nm)\"" << endl;
	}
	else {
		fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;
	}
	fgnu << "set ylabel \"order number + norm flux\"" << endl;

	if (!outputPlotEPSFileName.empty()) {
		fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
		fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
		fgnu << endl;

		if (subtractCentralWavelength) {
			fgnu << "plot \"" << atlasdatafilename << "\" u ($2-$4):($1 + $3) w l lt 4, ";
			fgnu << "\"" << compdatafilename << "\" u ($3-$5):($1 + $4) w l lt 3, ";
			fgnu << "\"" << linesdatafilename << "\" u ($3-$6):($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
		}
		else {
			fgnu << "plot \"" << atlasdatafilename << "\" u 2:($1 + $3) w l lt 4, ";
			fgnu << "\"" << compdatafilename << "\" u 3:($1 + $4) w l lt 3, ";
			fgnu << "\"" << linesdatafilename << "\" u 3:($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
		}

		if (display) {
			fgnu << "\nset terminal x11" << endl;
			fgnu << "set output" << endl;
			fgnu << "replot" << endl;
		}
		else {
			fgnu << "\n#set terminal x11" << endl;
			fgnu << "#set output" << endl;
			fgnu << "#replot" << endl;
		}
	}
	else {
		fgnu << endl;

		if (subtractCentralWavelength) {
			fgnu << "plot \"" << atlasdatafilename << "\" u ($2-$4):($1 + $3) w l lt 4, ";
			fgnu << "\"" << compdatafilename << "\" u ($3-$5):($1 + $4) w l lt 3, ";
			fgnu << "\"" << linesdatafilename << "\" u ($3-$6):($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
		}
		else {
			fgnu << "plot \"" << atlasdatafilename << "\" u 2:($1 + $3) w l lt 4, ";
			fgnu << "\"" << compdatafilename << "\" u 3:($1 + $4) w l lt 3, ";
			fgnu << "\"" << linesdatafilename << "\" u 3:($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
		}

		fgnu << endl;

		fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
		fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
		fgnu << "#replot" << endl;
		fgnu << "#set terminal x11" << endl;
		fgnu << "#set output" << endl;
	}

	fgnu.close();

	if (display) systemf("gnuplot -persist %s", gnuScriptFileName.c_str());
	else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s", gnuScriptFileName.c_str());
}

// Generate 3D plot for spectra of atlas + comparison + identified lines 
void GenerateWavelength3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, bool subtractCentralWavelength, bool display)
{
	if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());

	fgnu << "reset" << endl;
	fgnu << "unset key" << endl;
	fgnu << "set view 0,0" << endl;

	fgnu << "set palette gray" << endl;
	fgnu << "set palette gamma 2.0" << endl;
	fgnu << "set pm3d map" << endl;
	fgnu << "unset ztics" << endl;
	fgnu << "set cblabel \"normalized flux\"" << endl;

	if (subtractCentralWavelength) {
		fgnu << "\nset xlabel \"{/Symbol l} - {/Symbol l}_c (nm)\"" << endl;
	}
	else {
		fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;
	}
	fgnu << "set ylabel \"order number\"" << endl;

	fgnu << endl;

	if (!outputPlotEPSFileName.empty()) {
		fgnu << "\nset terminal postscript enhanced mono solid lw 1.5 \"Helvetica\" 14" << endl;
		fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
		if (subtractCentralWavelength) {
			fgnu << "\nsplot \"" << atlasdatafilename << "\" u ($2-$5):($1 + $4/4 - 0.5 + 0.25 + 0.25 + 0.175):3 w pm3d, ";
			fgnu << "\"" << compdatafilename << "\" u ($3-$7):($1 + $6/4 - 0.5 + 0.25 + 0.125):(1-$5) w pm3d, ";
			fgnu << "\"" << compdatafilename << "\" u ($3-$7):($1 + $6/4 - 0.5 + 0.075):4 w pm3d" << endl;
		}
		else {
			fgnu << "\nsplot \"" << atlasdatafilename << "\" u ($2):($1 + $4/4 - 0.5 + 0.25 + 0.25 + 0.175):3 w pm3d, ";
			fgnu << "\"" << compdatafilename << "\" u ($3):($1 + $6/4 - 0.5 + 0.25 + 0.125):(1-$5) w pm3d, ";
			fgnu << "\"" << compdatafilename << "\" u ($3):($1 + $6/4 - 0.5 + 0.075):4 w pm3d" << endl;

		}
		if (display) {
			fgnu << "\nset terminal x11" << endl;
			fgnu << "set output" << endl;
			fgnu << "replot" << endl;
		}
		else {
			fgnu << "\n#set terminal x11" << endl;
			fgnu << "#set output" << endl;
			fgnu << "#replot" << endl;
		}
	}
	else {

		if (subtractCentralWavelength) {
			fgnu << "\nsplot \"" << atlasdatafilename << "\" u ($2-$5):($1 + $4/4 - 0.5 + 0.25 + 0.25 + 0.175):3 w pm3d, ";
			fgnu << "\"" << compdatafilename << "\" u ($3-$7):($1 + $6/4 - 0.5 + 0.25 + 0.125):(1-$5) w pm3d, ";
			fgnu << "\"" << compdatafilename << "\" u ($3-$7):($1 + $6/4 - 0.5 + 0.075):4 w pm3d" << endl;
		}
		else {
			fgnu << "\nsplot \"" << atlasdatafilename << "\" u ($2):($1 + $4/4 - 0.5 + 0.25 + 0.25 + 0.175):3 w pm3d, ";
			fgnu << "\"" << compdatafilename << "\" u ($3):($1 + $6/4 - 0.5 + 0.25 + 0.125):(1-$5) w pm3d, ";
			fgnu << "\"" << compdatafilename << "\" u ($3):($1 + $6/4 - 0.5 + 0.075):4 w pm3d" << endl;

		}
		fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
		fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
		fgnu << "#replot" << endl;
		fgnu << "#set terminal x11" << endl;
		fgnu << "#set output" << endl;
	}

	fgnu.close();

	if (display) systemf("gnuplot -persist %s", gnuScriptFileName.c_str());
	else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s", gnuScriptFileName.c_str());
}

// Generate plot for wavelength solution - NOT IMPLEMENTED YET
void GenerateWavelengthSolutionPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned npolynomials, Polynomial *polynomials[], bool display)
{
	if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name;
	ofstream fgnu(gnuScriptFileName.c_str());

	fgnu << "reset" << endl;
	fgnu << "unset key" << endl;
	fgnu << "\nset xlabel \"image rows (pixels)\"" << endl;
	fgnu << "set ylabel \"image cols (pixels)\"" << endl;

	fgnu << "set pointsize 0.5" << endl;

	//    fgnu << "set xrange[" << row0 << ":" << rowf << "]" << endl;
	//    fgnu << "set yrange[" << col0 << ":" << colf << "]" << endl;

	for (unsigned k = 0; k<npolynomials; k++) {
		fgnu << "poly" << k;
		polynomials[k]->printEquation(&fgnu);
	}

	if (!outputPlotEPSFileName.empty()) {
		fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
		fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;

		fgnu << "\nplot \"" << dataFileName << "\" u 4:3 w p pt 7";

		for (unsigned k = 0; k<npolynomials; k++) {
			fgnu << ", poly" << k << "f(x)";
		}
		fgnu << endl;

		if (display) {
			fgnu << "\nset terminal x11" << endl;
			fgnu << "set output" << endl;
			fgnu << "replot" << endl;
		}
		else {
			fgnu << "\n#set terminal x11" << endl;
			fgnu << "#set output" << endl;
			fgnu << "#replot" << endl;
		}
	}
	else {
		fgnu << "\nplot \"" << dataFileName << "\" u 4:3 w p pt 7";

		for (unsigned k = 0; k<npolynomials; k++) {
			fgnu << ", poly" << k << "f(x)";
		}
		fgnu << endl;

		fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
		fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
		fgnu << "#replot" << endl;
		fgnu << "#set terminal x11" << endl;
		fgnu << "#set output" << endl;
	}

	fgnu.close();

	if (display) systemf("gnuplot -persist %s", gnuScriptFileName.c_str());
	else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s", gnuScriptFileName.c_str());
}

ofstream WavelengthCalibration::fatlasdata;
ofstream WavelengthCalibration::fcompdata;
ofstream WavelengthCalibration::flinesdata;
ofstream WavelengthCalibration::fordersdata;
bool WavelengthCalibration::generate3DPlot;
bool WavelengthCalibration::skipPlots;
