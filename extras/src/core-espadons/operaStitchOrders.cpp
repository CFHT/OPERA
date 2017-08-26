/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaStitchOrders
 Version: 1.0
 Description: This module stitches orders together
 to start up with an OPERA module.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope
 Location: Hawaii USA
 Date: Jan/2014
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
#include "libraries/operaSpectralTools.h"
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFit.h"							// for operaFitSplineDouble
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

bool calculateWavelengthShiftByXCorrInRange(operaSpectralElements *RefSpectrum, operaSpectralElements *compSpectrum, double wl0, double wlf, double DWavelengthStep, double DWavelengthRange,  float nsigcut, double threshold, double &maxDWavelength, double &maxcorr);
bool calculateRadialVelocityShiftByXCorrInRange(operaSpectralElements *RefSpectrum, operaSpectralElements *compSpectrum, double wl0, double wlf, double RVStep, double RVRange,  float nsigcut, double threshold, double &maxRVShift, double &maxcorr);
void GenerateStitchOrdersPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display, int refOrder);

/*! \brief stitch orders together. */
/*! \file operaStitchOrders.cpp */
/*! \package operaStitchOrders */

using namespace std;

/*!
 * operaStitchOrders
 * \author Eder Martioli
 * \brief Stitch orders.
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

operaArgumentHandler args;

int main(int argc, char *argv[])
{
    string inputSpectrum;
    string inputWaveFile;
    string outputWaveFile;
    int orderOfReference = 51;
    
    double RVStep = 0.1;
    double RVRange = 2;
    double XCorrelationThreshold = 0.05;
    double sigmaThreshold = 1.0;
    
    string plotfilename;
    string datafilename;
    string scriptfilename;
    bool interactive = false;
    
    args.AddRequiredArgument("inputSpectrum", inputSpectrum, "Input object spectrum file (.s)");
	args.AddRequiredArgument("outputWaveFile", outputWaveFile, "Output wavelength calibration file to store final solution");
	args.AddRequiredArgument("inputWaveFile", inputWaveFile, "Input wavelength calibration file");
    args.AddRequiredArgument("orderOfReference", orderOfReference, "Order taken as reference (fixed) for stitching");
    args.AddRequiredArgument("RVRange", RVRange, "Radial velocity (km/s) range to search for shift");
    args.AddRequiredArgument("RVStep", RVStep, "Radial velocity precision (km/s) to search for shift");
    args.AddRequiredArgument("XCorrelationThreshold", XCorrelationThreshold, "XCorrelation minimum threshold to accept shift");
    args.AddRequiredArgument("sigmaThreshold", sigmaThreshold, "Threshold in units of sigma to find peak correlation");
    args.AddPlotFileArguments(plotfilename, datafilename, scriptfilename, interactive);

	try {
		args.Parse(argc, argv);
		
		// we need an input spectrum...
		if (inputSpectrum.empty()) {
			throw operaException("operaStitchOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        if (outputWaveFile.empty()) {
			throw operaException("operaStitchOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

        if (args.verbose) {
			cout << "operaStitchOrders: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaStitchOrders: inputSpectrum = " << inputSpectrum << endl;            
			cout << "operaStitchOrders: outputWaveFile = " << outputWaveFile << endl;
			cout << "operaStitchOrders: orderOfReference = " << orderOfReference << endl;
			cout << "operaStitchOrders: RVRange = " << RVRange << endl;
			cout << "operaStitchOrders: RVStep = " << RVStep << endl;
			cout << "operaStitchOrders: XCorrelationThreshold = " << XCorrelationThreshold << endl;
			cout << "operaStitchOrders: sigmaThreshold = " << sigmaThreshold << endl;
            if(args.plot) {
                cout << "operaExtractionApertureCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaExtractionApertureCalibration: datafilename = " << datafilename << endl;
                cout << "operaExtractionApertureCalibration: scriptfilename = " << scriptfilename << endl;
                cout << "operaExtractionApertureCalibration: interactive = " << (interactive ? "YES" : "NO") << endl;
            }
		}
        
        
        ofstream fdata;
        if (!datafilename.empty()) fdata.open(datafilename.c_str());

 		operaSpectralOrderVector spectralOrders;
 		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputSpectrum);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile); // This merges in the wavelength calibration information
        
        int ordernumber = NOTPROVIDED, minorder = NOTPROVIDED, maxorder = NOTPROVIDED;
        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        if (args.verbose) cout << "operaStitchOrders: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        for(unsigned order=(unsigned)(orderOfReference)+1; order<=(unsigned)maxorder; order++) {
            
            if (args.verbose) cout << "operaStitchOrders: calibrating higher orders: order=" << order << " up to " << maxorder << endl;

            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                
                operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
                operaWavelength *wavelength = spectralOrder->getWavelength();
                spectralElements->setwavelengthsFromCalibration(wavelength);
                unsigned neighborOrder = order-1;
                
                operaSpectralOrder *neighborSpectralOrder = spectralOrders.GetSpectralOrder(neighborOrder);
                
                if (neighborSpectralOrder->gethasSpectralElements() && neighborSpectralOrder->gethasWavelength()) {
                    
                    operaSpectralElements *neighborElements = neighborSpectralOrder->getSpectralElements();
                    operaWavelength *neighborWavelength = neighborSpectralOrder->getWavelength();
                    neighborElements->setwavelengthsFromCalibration(neighborWavelength);
                    
                    double wl0=0, wlf=0;
                    bool overlap = getOverlappingWLRange(neighborElements,spectralElements,wl0,wlf);
                    if(overlap) {
                        double rvshift = 0;
                        double maxcorr = 0;
                        
                        // find rv shift that best match neighbor order
                        bool xcorrstatus = calculateRadialVelocityShiftByXCorrInRange(neighborElements,spectralElements,wl0,wlf, RVStep, RVRange,sigmaThreshold,XCorrelationThreshold,rvshift,maxcorr);

                        if(xcorrstatus) {
                            
                            // for plotting
                            if (fdata.is_open()) {
                                double wlc = (wl0 + (wlf-wl0)/2.0);
                                double wlShift = rvshift * wlc / SPEED_OF_LIGHT_KMS;
                                fdata << order << " " << wlc << " " << rvshift << " " << wlShift << " " << maxcorr << endl;
                            }
        
                            // update wavelength solution with new rv shift
                            wavelength->applyRadialVelocityCorrection(rvshift);
                        }
                        
                    } else {
                        if (args.verbose) cout << "operaStitchOrders: breaking at order " << order << " -> no overlapping with neighbor order " << neighborOrder << endl;
                        break;
                    }
                } else {
                    if (args.verbose) cout << "operaStitchOrders: breaking at order " << order << " -> no elements/wavelength in neighbor order " << neighborOrder << endl;
                    break;
                }
            } else {
                if (args.verbose) cout << "operaStitchOrders: breaking at order " << order << " -> no elements/wavelength" << endl;
                break;
            }
        }
        
        // for plotting
        if (fdata.is_open()) {
            fdata << endl;
            fdata << orderOfReference << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 1.0 << endl;
            fdata << endl;
        }
        
        for(unsigned order=(unsigned)(orderOfReference)-1; order>=(unsigned)minorder; order--) {
            
            if (args.verbose) cout << "operaStitchOrders: calibrating lower orders: order=" << order << " down to " << minorder << endl;
            
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                
                operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
                operaWavelength *wavelength = spectralOrder->getWavelength();
                spectralElements->setwavelengthsFromCalibration(wavelength);
                unsigned neighborOrder = order+1;
                
                operaSpectralOrder *neighborSpectralOrder = spectralOrders.GetSpectralOrder(neighborOrder);
                
                if (neighborSpectralOrder->gethasSpectralElements() && neighborSpectralOrder->gethasWavelength()) {
                    
                    operaSpectralElements *neighborElements = neighborSpectralOrder->getSpectralElements();
                    operaWavelength *neighborWavelength = neighborSpectralOrder->getWavelength();
                    neighborElements->setwavelengthsFromCalibration(neighborWavelength);
                    
                    double wl0=0, wlf=0;
                    bool overlap = getOverlappingWLRange(neighborElements,spectralElements,wl0,wlf);
                    if(overlap) {
                        double rvshift = 0;
                        double maxcorr = 0;
                        // find rv shift that best match neighbor order
                        bool xcorrstatus = calculateRadialVelocityShiftByXCorrInRange(neighborElements,spectralElements,wl0,wlf, RVStep, RVRange,sigmaThreshold,XCorrelationThreshold,rvshift,maxcorr);
                        
                        if(xcorrstatus) {
                            
                            // for plotting
                            if (fdata.is_open()) {
                                double wlc = (wl0 + (wlf-wl0)/2.0);
                                double wlShift = rvshift * wlc / SPEED_OF_LIGHT_KMS;
                                fdata << order << " " << wlc << " " << rvshift << " " << wlShift << " " << maxcorr << endl;
                            }

                            // update wavelength solution with new rv shift
                            wavelength->applyRadialVelocityCorrection(rvshift);
                        }
                        
                    } else {
                        if (args.verbose) cout << "operaStitchOrders: breaking at order " << order << " -> no overlapping with neighbor order " << neighborOrder << endl;
                        break;
                    }
                } else {
                    if (args.verbose) cout << "operaStitchOrders: breaking at order " << order << " -> no elements/wavelength in neighbor order " << neighborOrder << endl;
                    break;
                }
            } else {
                if (args.verbose) cout << "operaStitchOrders: breaking at order " << order << " -> no elements/wavelength" << endl;
                break;
            }
        }

        operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputWaveFile, Wave);
        
        if (fdata.is_open()) {
            fdata.close();
            if (!scriptfilename.empty()) {
                GenerateStitchOrdersPlot(scriptfilename,plotfilename,datafilename,interactive, orderOfReference);
            }
        }

	}
	catch (operaException e) {
		cerr << "operaStitchOrders: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaStitchOrders: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


void GenerateStitchOrdersPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display, int refOrder)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
    remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
    ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    
    fgnu << "set xlabel \"order number\"" << endl;
    fgnu << "set ylabel \"Radial velocity shift (km/s)\"" << endl;

    fgnu << "unset key" << endl;

    fgnu << "set arrow from " <<  refOrder << ",0.95 to " << refOrder << ",0.2 lw 2" << endl;
    fgnu << "set label \"Ref order\" at " << float(refOrder) - 0.5 << ",1.1" << endl;
    
    fgnu << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nplot \"" << dataFileName << "\" u 1:3 w linespoint pt 7 ps 2";
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
        fgnu << "\nplot \"" << dataFileName << "\" u 1:3 w linespoint pt 7 ps 2";

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


bool calculateRadialVelocityShiftByXCorrInRange(operaSpectralElements *RefSpectrum, operaSpectralElements *compSpectrum, double wl0, double wlf, double RVStep, double RVRange,  float nsigcut, double threshold, double &maxRVShift, double &maxcorr)
{
    bool status = false;
    
    double *refSpectrumFlux = new double[RefSpectrum->getnSpectralElements()];
    double *refSpectrumWavelength = new double[RefSpectrum->getnSpectralElements()];
    unsigned nref = getSpectrumWithinWLRange(RefSpectrum,wl0,wlf,refSpectrumFlux,refSpectrumWavelength);
    for (unsigned i=0; i<nref; i++) {
        if(isnan(refSpectrumFlux[i])) refSpectrumFlux[i] = 0;
    }
    
    double wlc = wl0 + (wlf - wl0)/2.0;
    double DWavelengthRange = RVRange*wlc/SPEED_OF_LIGHT_KMS;
    
    double *compSpectrumFlux = new double[compSpectrum->getnSpectralElements()];
    double *compSpectrumWavelength = new double[compSpectrum->getnSpectralElements()];
    
    unsigned ncomp = getSpectrumWithinWLRange(compSpectrum,wl0-DWavelengthRange/2.0,wlf+DWavelengthRange/2.0,compSpectrumFlux,compSpectrumWavelength);
    for (unsigned i=0; i<ncomp; i++) {
        if(isnan(compSpectrumFlux[i])) compSpectrumFlux[i] = 0;
    }
    
    double *compSpectrumWavelength_mod = new double[compSpectrum->getnSpectralElements()];
    double *compSpectrumFlux_mod = new double[RefSpectrum->getnSpectralElements()];
    
    unsigned nDataPoints = (unsigned)ceil( (RVRange / RVStep) + 1.0);
    double rvshift = -RVRange/2.0;
    
    unsigned jmax = 0;
    maxcorr = -BIG;
    maxRVShift = 0;
    
    if (nref == 0 || ncomp ==0) {
        maxcorr = NAN;
        maxRVShift = NAN;
        return false;
    }
    
    float *crosscorrelation = new float [nDataPoints];
    float *drv = new float [nDataPoints];
    
    for (unsigned j=0; j<nDataPoints;j++) {
        for (unsigned i=0; i<ncomp; i++) {
            double DWavelength = rvshift * compSpectrumWavelength[i] / SPEED_OF_LIGHT_KMS;
            compSpectrumWavelength_mod[i] = compSpectrumWavelength[i] + DWavelength;
        }
        
        operaFitSplineDouble(ncomp,compSpectrumWavelength_mod,compSpectrumFlux,nref,refSpectrumWavelength,compSpectrumFlux_mod);
        
        crosscorrelation[j] = (float)(operaCrossCorrelation(nref,refSpectrumFlux,compSpectrumFlux_mod));
        drv[j] = (float)rvshift;
        
        if((double)(crosscorrelation[j]) > maxcorr && (double)(crosscorrelation[j]) > threshold) {
            maxcorr = (double)(crosscorrelation[j]);
            maxRVShift = rvshift;
            jmax = j;
            status = true;
        }
        
        if(args.debug) cout << drv[j] << " "  << crosscorrelation[j] << endl;
        
        rvshift+=RVStep;
    }
    
    double firstmaxRVShift = maxRVShift;
    double firstMaxcorr = maxcorr;
    
    if (firstmaxRVShift == -RVRange/2.0) {
        maxcorr = NAN;
        maxRVShift = NAN;
        status = false;
    }
    
    if(args.debug) cout << "calculateWavelengthShiftByXCorrInRange: status=" << status << " maxRVShift=" << maxRVShift  <<  " maxcorr=" << maxcorr << endl;
    
    /*
     * Two tests will be performed:
     *  I) Whether it is an isolated maximum point, i.e. three points before and three points after all must be lower than max
     *
     *  II) Whether it's a significant maximum, i.e. maxcorrelation > median + 3*sigma
     */
    
    unsigned halfslitsize = 5;
    
    double *peakXdata = new double[2*halfslitsize+2];
    double *peakYdata = new double[2*halfslitsize+2];
    unsigned np = 0;
    // Test (I)
    
    if((int)jmax - (int)halfslitsize < 0 || jmax + halfslitsize >= nDataPoints) {
        maxcorr = NAN;
        maxRVShift = NAN;
        status = false;
    } else {
        float goingup = -BIG;
        for(unsigned j=jmax-halfslitsize;j<jmax;j++) {
            peakXdata[np] = (double)drv[j];
            peakYdata[np] = (double)crosscorrelation[j];
            np++;
            if(crosscorrelation[j] > goingup) {
                goingup = crosscorrelation[j];
            } else {
                maxcorr = NAN;
                maxRVShift = NAN;
                np=0;
                status = false;
                break;
            }
        }
        
        peakXdata[np] = (double)drv[jmax];
        peakYdata[np] = (double)crosscorrelation[jmax];
        np++;
        
        float goingdown = BIG;
        for(unsigned j=jmax+1;j<=jmax+halfslitsize;j++) {
            peakXdata[np] = (double)drv[j];
            peakYdata[np] = (double)crosscorrelation[j];
            np++;
            if(crosscorrelation[j] < goingdown) {
                goingdown = crosscorrelation[j];
            } else {
                maxcorr = NAN;
                maxRVShift = NAN;
                np=0;
                status = false;
                break;
            }
        }
    }
    if(args.debug)
        cout << "calculateWavelengthShiftByXCorr: Test I status=" << status << endl;
    
    // Test (II)
    
    float medianXcorr = operaArrayMedian(nDataPoints,crosscorrelation);
    float medsigXcorr = operaArrayMedianSigma(nDataPoints,crosscorrelation,medianXcorr);
    
    if(args.debug)
        cout << "calculateWavelengthShiftByXCorr:maxcorr=" << crosscorrelation[jmax] << "  medianXcorr=" << medianXcorr << "  medsigXcorr=" << medsigXcorr << "  median+n*sig=" << medianXcorr + nsigcut*medsigXcorr << endl;
    
    if(crosscorrelation[jmax] < medianXcorr + nsigcut*medsigXcorr) {
        maxcorr = NAN;
        maxRVShift = NAN;
        status = false;
    }
    
    if(args.debug) cout << "calculateWavelengthShiftByXCorr: Test II status=" << status << endl;
    
    if(status) {
        double a=(double)crosscorrelation[jmax];
        double x0=(double)drv[jmax];
        double sig=(double)(fabs(drv[jmax+halfslitsize/2]-drv[jmax-halfslitsize/2]))/2.0;
        double chisqr;
        
        operaLMFitGaussian(np, peakXdata, peakYdata, &a, &x0, &sig, &chisqr);
        maxcorr = a;
        maxRVShift = x0;
    } else if(firstMaxcorr > threshold) {
        maxRVShift = firstmaxRVShift;
        maxcorr = firstMaxcorr;
        status = true;
    }
    
    if(args.debug) cout << "calculateWavelengthShiftByXCorr: maxRVShift=" << maxRVShift << endl;
    
    delete[] peakXdata;
    delete[] peakYdata;
    
    delete[] crosscorrelation;
    delete[] drv;
    
    delete[] refSpectrumFlux;
    delete[] refSpectrumWavelength;
    delete[] compSpectrumFlux;
    delete[] compSpectrumWavelength_mod;
    delete[] compSpectrumFlux_mod;
    
    return status;
}




bool calculateWavelengthShiftByXCorrInRange(operaSpectralElements *RefSpectrum, operaSpectralElements *compSpectrum, double wl0, double wlf, double DWavelengthStep, double DWavelengthRange,  float nsigcut, double threshold, double &maxDWavelength, double &maxcorr)
{
    bool status = false;
    
    double *refSpectrumFlux = new double[RefSpectrum->getnSpectralElements()];
    double *refSpectrumWavelength = new double[RefSpectrum->getnSpectralElements()];
    unsigned nref = getSpectrumWithinWLRange(RefSpectrum,wl0,wlf,refSpectrumFlux,refSpectrumWavelength);
    for (unsigned i=0; i<nref; i++) {
        if(isnan(refSpectrumFlux[i])) refSpectrumFlux[i] = 0;
    }
    
    double *compSpectrumFlux = new double[compSpectrum->getnSpectralElements()];
    double *compSpectrumWavelength = new double[compSpectrum->getnSpectralElements()];
    unsigned ncomp = getSpectrumWithinWLRange(compSpectrum,wl0-DWavelengthRange/2.0,wlf+DWavelengthRange/2.0,compSpectrumFlux,compSpectrumWavelength);
    for (unsigned i=0; i<ncomp; i++) {
        if(isnan(compSpectrumFlux[i])) compSpectrumFlux[i] = 0;
    }
    
    double *compSpectrumWavelength_mod = new double[compSpectrum->getnSpectralElements()];
    double *compSpectrumFlux_mod = new double[RefSpectrum->getnSpectralElements()];

    unsigned nDataPoints = (unsigned)ceil( (DWavelengthRange / DWavelengthStep) + 1.0);
    double DWavelength = -DWavelengthRange/2.0;
    
    unsigned jmax = 0;
    maxcorr = -BIG;
    maxDWavelength = 0;
	
	if (nref == 0 || ncomp ==0) {
		maxcorr = NAN;
		maxDWavelength = NAN;
		return false;
	}

    float *crosscorrelation = new float [nDataPoints];
    float *dwl = new float [nDataPoints];

    for (unsigned j=0; j<nDataPoints;j++) {
        for (unsigned i=0; i<ncomp; i++) {
            compSpectrumWavelength_mod[i] = compSpectrumWavelength[i] + DWavelength;
        }
        
        operaFitSplineDouble(ncomp,compSpectrumWavelength_mod,compSpectrumFlux,nref,refSpectrumWavelength,compSpectrumFlux_mod);
        
        crosscorrelation[j] = (float)(operaCrossCorrelation(nref,refSpectrumFlux,compSpectrumFlux_mod));
        dwl[j] = (float)DWavelength;
                
        if((double)(crosscorrelation[j]) > maxcorr && (double)(crosscorrelation[j]) > threshold) {
            maxcorr = (double)(crosscorrelation[j]);
            maxDWavelength = DWavelength;
            jmax = j;
            status = true;
        }

        if(args.debug) cout << dwl[j] << " "  << crosscorrelation[j] << endl;
        
        DWavelength+=DWavelengthStep;
    }
    
    double firstMaxDWavelength = maxDWavelength;
    double firstMaxcorr = maxDWavelength;    
  
    if (maxDWavelength == -DWavelengthRange/2.0) {
		maxcorr = NAN;
		maxDWavelength = NAN;
        status = false;
	}
    
    if(args.debug) cout << "calculateWavelengthShiftByXCorrInRange: status=" << status << " wlshift=" << maxDWavelength  <<  " maxcorr=" << maxcorr << endl;

     /*
     * Two tests will be performed:
     *  I) Whether it is an isolated maximum point, i.e. three points before and three points after all must be lower than max
     *
     *  II) Whether it's a significant maximum, i.e. maxcorrelation > median + 3*sigma
     */
    
    unsigned halfslitsize = 5;
    
    double *peakXdata = new double[2*halfslitsize+2];
    double *peakYdata = new double[2*halfslitsize+2];
    unsigned np = 0;
    // Test (I)
    
    if((int)jmax - (int)halfslitsize < 0 || jmax + halfslitsize >= nDataPoints) {
        maxcorr = NAN;
		maxDWavelength = NAN;
        status = false;
    } else {
        float goingup = -BIG;
        for(unsigned j=jmax-halfslitsize;j<jmax;j++) {
            peakXdata[np] = (double)dwl[j];
            peakYdata[np] = (double)crosscorrelation[j];
            np++;
            if(crosscorrelation[j] > goingup) {
                goingup = crosscorrelation[j];
            } else {
                maxcorr = NAN;
                maxDWavelength = NAN;
                np=0;
                status = false;
                break;
            }
        }
        
        peakXdata[np] = (double)dwl[jmax];
        peakYdata[np] = (double)crosscorrelation[jmax];
        np++;
        
        float goingdown = BIG;
        for(unsigned j=jmax+1;j<=jmax+halfslitsize;j++) {
            peakXdata[np] = (double)dwl[j];
            peakYdata[np] = (double)crosscorrelation[j];
            np++;
            if(crosscorrelation[j] < goingdown) {
                goingdown = crosscorrelation[j];
            } else {
                maxcorr = NAN;
                maxDWavelength = NAN;
                np=0;
                status = false;
                break;
            }
        }
    }
    if(args.debug)
        cout << "calculateWavelengthShiftByXCorr: Test I status=" << status << endl;
    
    // Test (II)
    
    float medianXcorr = operaArrayMedian(nDataPoints,crosscorrelation);
    float medsigXcorr = operaArrayMedianSigma(nDataPoints,crosscorrelation,medianXcorr);
    
    if(args.debug)
        cout << "calculateWavelengthShiftByXCorr:maxcorr=" << crosscorrelation[jmax] << "  medianXcorr=" << medianXcorr << "  medsigXcorr=" << medsigXcorr << "  median+n*sig=" << medianXcorr + nsigcut*medsigXcorr << endl;
    
    if(crosscorrelation[jmax] < medianXcorr + nsigcut*medsigXcorr) {
        maxcorr = NAN;
        maxDWavelength = NAN;
        status = false;
    }
    
    if(args.debug) cout << "calculateWavelengthShiftByXCorr: Test II status=" << status << endl;
    
    if(status) {
        double a=(double)crosscorrelation[jmax];
        double x0=(double)dwl[jmax];
        double sig=(double)(fabs(dwl[jmax+halfslitsize/2]-dwl[jmax-halfslitsize/2]))/2.0;
        double chisqr;
        
        operaLMFitGaussian(np, peakXdata, peakYdata, &a, &x0, &sig, &chisqr);
        maxcorr = a;
        maxDWavelength = x0;
    } else if(firstMaxcorr > threshold) {
		maxDWavelength = firstMaxDWavelength;
		maxcorr = firstMaxcorr;
		status = true;
    }
    
    if(args.debug) cout << "calculateWavelengthShiftByXCorr: maxDWavelength=" << maxDWavelength << endl;
    
    delete[] peakXdata;
    delete[] peakYdata;

    delete[] crosscorrelation;
    delete[] dwl;
    
    delete[] refSpectrumFlux;
    delete[] refSpectrumWavelength;
    delete[] compSpectrumFlux;
    delete[] compSpectrumWavelength_mod;
    delete[] compSpectrumFlux_mod;

	return status;
}
