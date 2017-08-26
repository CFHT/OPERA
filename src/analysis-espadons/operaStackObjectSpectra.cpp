/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaStackObjectSpectra
 Version: 1.0
 Description: Stack several object spectra into a master higher SNR spectrum
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2016
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
#include <iomanip>
#include "libraries/operaIOFormats.h"
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"
#include "libraries/operaSpectralTools.h"

#define MAXNUMBEROFSPECTRUMFILES 1000

/*! \file operaStackObjectSpectra.cpp */

using namespace std;

/*! 
 * operaStackObjectSpectra
 * \author Eder Martioli
 * \brief Stack several object spectra into a master higher SNR spectrum
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
	
    string listofinputspectra;
    string inputspectrum;
    
    string listofinputRVfiles;
    string inputRVfile;
    
    string object;

    unsigned combineMethod = 0; // Sum=1, Mean=2, Median=3, WeightedMean=0 (default).
    unsigned spectrumTypeToExtract = 0; // 1=Normalized, 2=FluxCalibrated, 0=Default=Raw
    
    double snrClip = 10;
    unsigned numberOfPointsToCutInOrderEnds = 400;
    
    double RadialVelocityBin = 2.0; // km/s
    
    double firstWavelength = 0;
    double lastWavelength = 0;

    bool applyTelluricWaveCorrection = TRUE;
    bool applyHeliocentricRVCorrection = TRUE;

    string outputspectrum;
    
    int ordernumber = NOTPROVIDED;
    int minorder = 22;
    int maxorder = 62;
    
    args.AddOptionalArgument("inputspectrum", inputspectrum, "", "input spectrum filename");
    args.AddOptionalArgument("listofinputspectra", listofinputspectra, "", "List of input spectra");
    args.AddOptionalArgument("inputRVfile", inputRVfile, "", "input Radial Velocity filename");
    args.AddOptionalArgument("listofinputRVfiles", listofinputRVfiles, "", "List of input RV files");
    args.AddOptionalArgument("object", object, "UNKNOWN_OBJECT", "Object name");
    
    args.AddOptionalArgument("combineMethod", combineMethod, 0, "Method for combining spectra:  0 (default) = Weighted mean, 1 = Sum, 2 = Mean, 3 = Median");
    args.AddOptionalArgument("spectrumTypeToExtract", spectrumTypeToExtract, 0, "Spectrum Type to extract: 0 (Default) = Raw, 1 = Normalized, 2 = FluxCalibrated");

    args.AddOptionalArgument("RadialVelocityBin", RadialVelocityBin, 2.0, "Radial velocity bin to define sampling of output spectrum");
    args.AddOptionalArgument("firstWavelength", firstWavelength, 0, "first wavelength in output spectrum");
    args.AddOptionalArgument("lastWavelength", lastWavelength, 0, "last wavelength in output spectrum");

    args.AddOptionalArgument("snrClip", snrClip, 1.0, "Mininum SNR threshold to accept points in spectrum");
    args.AddOptionalArgument("numberOfPointsToCutInOrderEnds", numberOfPointsToCutInOrderEnds, 0, "Number of points to exclude in order edges");
    
    args.AddOptionalArgument("applyTelluricWaveCorrection", applyTelluricWaveCorrection, true, "Use telluric wavelength correction");
    args.AddOptionalArgument("applyHeliocentricRVCorrection", applyHeliocentricRVCorrection, true, "Use heliocentric radial velocity correction");

    args.AddRequiredArgument("outputspectrum", outputspectrum, "output stacked spectrum file");
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    
	try {
        args.Parse(argc, argv);
        
        /*
         * Get all input spectrum files both from --inputspectrum=.. and from list
         */
        string inputspectra[MAXNUMBEROFSPECTRUMFILES];
		unsigned inputIndex = 0;
		SplitStringIntoArray(inputspectrum, inputspectra, inputIndex, MAXNUMBEROFSPECTRUMFILES); // Split list of images into array
        ReadStringsFromFileIntoArray(listofinputspectra, inputspectra, inputIndex, MAXNUMBEROFSPECTRUMFILES); // Read list of images from file

        if(inputIndex == 0) {
			throw operaException("operaStackObjectSpectra: no input spectrum. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        } else if (inputIndex > MAXNUMBEROFSPECTRUMFILES) {
            throw operaException("operaStackObjectSpectra: number of input spectra exceeds limit. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        //---

        /*
         * Get all input RV files both from --inputRVfile=.. and from list
         */
        string inputrv[MAXNUMBEROFSPECTRUMFILES];
        unsigned inputRVIndex = 0;
        SplitStringIntoArray(inputRVfile, inputrv, inputRVIndex, MAXNUMBEROFSPECTRUMFILES); // Split list of rv files into array
        ReadStringsFromFileIntoArray(listofinputRVfiles, inputrv, inputRVIndex, MAXNUMBEROFSPECTRUMFILES); // Read list of RV files

        bool hasAllRVFiles = true;
        if (inputRVIndex != inputIndex) {
            hasAllRVFiles = false;
            //throw operaException("operaStackObjectSpectra: the number of RV files is different than the number of input spectra. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        //---
        
        // We need an output spectrum file name
		if (outputspectrum.empty()) {
			throw operaException("operaStackObjectSpectra: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (args.verbose) {
            for(unsigned index=0; index<inputIndex; index++) {
                cout << "operaStackObjectSpectra: inputspectra["<<index<<"] = " << inputspectra[index] << endl;
                if (hasAllRVFiles) {
                    cout << "operaStackObjectSpectra: inputrv["<<index<<"] = " << inputrv[index] << endl;
                }
            }
            if (!hasAllRVFiles) {
                cout << "operaStackObjectSpectra: RV files do not match spectral files -> using RV=0 for source." << endl;
            }
			cout << "operaStackObjectSpectra: output spectrum file = " << outputspectrum << endl;
            cout << "operaStackObjectSpectra: combineMethod = " << combineMethod << endl;
            cout << "operaStackObjectSpectra: spectrumTypeToExtract = " << spectrumTypeToExtract << endl;

            cout << "operaStackObjectSpectra: RadialVelocityBin = " << RadialVelocityBin << endl;
            cout << "operaStackObjectSpectra: firstWavelength = " << firstWavelength << endl;
            cout << "operaStackObjectSpectra: lastWavelength = " << lastWavelength << endl;
            cout << "operaStackObjectSpectra: snrClip = " << snrClip << endl;
            cout << "operaStackObjectSpectra: numberOfPointsToCutInOrderEnds = " << numberOfPointsToCutInOrderEnds << endl;

            cout << "operaStackObjectSpectra: applyTelluricWaveCorrection = " << applyTelluricWaveCorrection << endl;
            cout << "operaStackObjectSpectra: applyHeliocentricRVCorrection = " << applyHeliocentricRVCorrection << endl;

            if(ordernumber != NOTPROVIDED) cout << "operaStackObjectSpectra: ordernumber = " << ordernumber << endl;
		}
        
        operaSpectrum spectrum[MAXNUMBEROFSPECTRUMFILES];
        
        double avg_mjdate = 0;
        double avg_HJD_UTC = 0;
        double avg_HJD_TT = 0;
        double np = 0;
        
        for (unsigned index=0; index<inputIndex; index++) {
            
            double RV_KPS = 0.0;

            if (hasAllRVFiles) {
                if (!inputrv[index].empty()) {
                    FormatData rvdata;
                    operaIOFormats::ReadCustomFormat("rv", rvdata, inputrv[index]);
                    
                    avg_mjdate += rvdata.extract<double>();
                    avg_HJD_UTC += rvdata.extract<double>();
                    avg_HJD_TT += rvdata.extract<double>();
                    //double wl0tmp = rvdata.extract<double>();
                    //double wlftmp = rvdata.extract<double>();
                    RV_KPS = rvdata.extract<double>()/1000.0;
                    double RVErr_mps = rvdata.extract<double>();
                    np++;
                }
            }

            if (args.verbose) {
                cout << "operaStackObjectSpectra: index=" << index << " applying source RV correction of " << RV_KPS << " km/s" << endl;
            }
            operaSpectralOrderVector spectralOrderVector;
            operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, inputspectra[index]);
            UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrderVector);
            
            spectrum[index] = spectralOrderVector.getExtendedSpectrum(minorder,maxorder,spectrumTypeToExtract,applyTelluricWaveCorrection,applyHeliocentricRVCorrection,snrClip,numberOfPointsToCutInOrderEnds,RV_KPS);
        }
        
        avg_mjdate /= (double)np; // average MJD
        avg_HJD_UTC /= (double)np; // average HJD_UTC
        avg_HJD_TT /= (double)np; // average HJD_TT
        
        /*
         *  Below we use the first and last wavelength in base spectrum if
         *  either one of these quantities are not provided
         */
        if (firstWavelength==0) {
            firstWavelength = spectrum[0].firstwl();
        }
        if (lastWavelength==0) {
            lastWavelength = spectrum[0].lastwl();
        }
        
        double wl = firstWavelength;
        
        unsigned firsti[MAXNUMBEROFSPECTRUMFILES];
        for (unsigned index=0; index<inputIndex; index++) {
            firsti[index] = 0;
        }
        
        float *fluxdata = new float[MAX_ESPADONS_VECTOR_LENGTH];
        
        operaSpectrum outspectrum;
        
        while (wl<lastWavelength) {
            double dwl = (RadialVelocityBin * wl)/SPEED_OF_LIGHT_KMS;
            double wlend = wl + dwl;
            
            operaWavelengthRange wlrange(wl,wlend);
            
            double flux = 0;
            double variance = 0;
            unsigned nin = 0;

            switch (combineMethod) {
                case 1: // Sum
                    
                    for (unsigned index=0; index<inputIndex; index++) {
                        for (unsigned i=firsti[index]; i<spectrum[index].size(); i++) {
                            if (spectrum[index].getwavelength(i) < wl) {
                                continue;
                            } else {
                                
                                if (wlrange.contains(spectrum[index].getwavelength(i))) {
                                    if (!isnan(spectrum[index].getflux(i))) {
                                        variance += spectrum[index].getvariance(i);
                                        flux += spectrum[index].getflux(i);
                                    }
                                } else {
                                    firsti[index] = i;
                                    break;
                                }
                            }
                        }
                    }
                    if (flux) {
                        double wlc = wl + dwl/2.0;
                        outspectrum.insert(wlc, flux, variance);
                    }

                    break;  // end of case 1
                    
                case 2: // Mean
                    
                    for (unsigned index=0; index<inputIndex; index++) {
                        for (unsigned i=firsti[index]; i<spectrum[index].size(); i++) {
                            if (spectrum[index].getwavelength(i) < wl) {
                                continue;
                            } else {
                                if (wlrange.contains(spectrum[index].getwavelength(i))) {
                                    if (spectrum[index].getflux(i) && !isnan(spectrum[index].getflux(i))) {
                                        variance += spectrum[index].getvariance(i);
                                        flux += spectrum[index].getflux(i);
                                        nin++;
                                    }
                                } else {
                                    firsti[index] = i;
                                    break;
                                }
                            }
                        }
                    }
                    if (nin) {
                        flux /= (double)nin;
                        variance /= (double)nin;
                        double wlc = wl + dwl/2.0;
                        outspectrum.insert(wlc, flux, variance);
                    }

                    break;  // end of case 2
                    
                case 3: // Median
                    
                    for (unsigned index=0; index<inputIndex; index++) {
                        for (unsigned i=firsti[index]; i<spectrum[index].size(); i++) {
                            if (spectrum[index].getwavelength(i) < wl) {
                                continue;
                            } else {
                                if (wlrange.contains(spectrum[index].getwavelength(i))) {
                                    if (spectrum[index].getflux(i) && !isnan(spectrum[index].getflux(i))) {
                                        fluxdata[nin] = (float)spectrum[index].getflux(i);
                                        nin++;
                                    }
                                } else {
                                    firsti[index] = i;
                                    break;
                                }
                            }
                        }
                    }
                    if (nin) {
                        
                        if (nin > 2) {
                            flux = (double)operaArrayMedian(nin,fluxdata);
                            float flux_err = operaArrayMedianSigma(nin, fluxdata, flux);
                            variance = (double)(flux_err*flux_err);
                        } else {
                            flux = (double)operaArrayMean(nin,fluxdata);
                            float flux_err = operaArraySigma(nin, fluxdata);
                            variance = (double)(flux_err*flux_err);
                        }

                        double wlc = wl + dwl/2.0;
                        outspectrum.insert(wlc, flux, variance);
                    }

                    
                    break;  // end of case 3
                    
                default: // Weighted Mean
                    
                    double weightsum = 0;
                    for (unsigned index=0; index<inputIndex; index++) {
                        
                        for (unsigned i=firsti[index]; i<spectrum[index].size(); i++) {
                            if (spectrum[index].getwavelength(i) < wl) {
                                continue;
                            } else {
                                if (wlrange.contains(spectrum[index].getwavelength(i))) {
                                    double weight = 1.0;
                                    if (spectrum[index].getflux(i) && !isnan(spectrum[index].getflux(i))) {
                                        if (spectrum[index].getvariance(i)) {
                                            weight = 1.0/(spectrum[index].getvariance(i))*(spectrum[index].getvariance(i));
                                            variance += spectrum[index].getvariance(i)*weight;
                                        }
                                        flux += spectrum[index].getflux(i)*weight;
                                        weightsum += weight;
                                    }
                                } else {
                                    firsti[index] = i;
                                    break;
                                }
                            }
                        }
                    }
                    if (weightsum) {
                        flux /= weightsum;
                        variance /= weightsum;
                        double wlc = wl + dwl/2.0;
                        outspectrum.insert(wlc, flux, variance);
                    }
                    
                    break; // end of Default case
            } // end of switch
    
            wl = wlend;
        }
        
        /*
         * write output
         */
        if (!outputspectrum.empty()) {
            ofstream fout;
            fout.open(outputspectrum.c_str());
            
            fout << "# MJD  =  " << fixed << setprecision(10) <<  avg_mjdate << endl;
            fout << "# HJD_UTC  =  " << fixed << setprecision(10) << avg_HJD_UTC << endl;
            fout << "# HJD_TT  =  " << fixed << setprecision(10) << avg_HJD_TT << endl;
            fout << "# Object  =  " << object << endl;
            fout << "# FirstWavelength_nm  =  " << firstWavelength << endl;
            fout << "# LastWavelength_nm  =  " << lastWavelength << endl;
            fout << "# RadialVelocityBin_KPS  =  " << RadialVelocityBin << endl;
            fout << "# SNRCut  =  " << snrClip << endl;
            fout << "# numberOfPointsToCutInOrderEnds  =  " << numberOfPointsToCutInOrderEnds << endl;
            
            switch (spectrumTypeToExtract) {
                case 1:
                    fout << "# spectrumTypeToExtract  =  Normalized" << endl;
                    break;
                case 2:
                    fout << "# spectrumTypeToExtract  =  FluxCalibrated" << endl;
                    break;
                default:
                    fout << "# spectrumTypeToExtract  =  Raw" << endl;
                    break;
            }

            switch (combineMethod) {
                case 1:
                    fout << "# combineMethod  =  Sum " <<  endl;
                    break;
                case 2:
                    fout << "# combineMethod  =  Mean" << endl;
                    break;
                case 3:
                    fout << "# combineMethod  =  Median" << endl;
                    break;
                default:
                    fout << "# combineMethod  =  WeightedMean" << endl;
                    break;
            }
            
            if (hasAllRVFiles) {
                fout << "# ApplySourceRVCorrection =  YES " << endl;
            } else {
                fout << "# ApplySourceRVCorrection =  NO " << endl;
            }
            
            if (applyTelluricWaveCorrection) {
                fout << "# ApplyTelluricWaveCorrection =  YES " << endl;
            } else {
                fout << "# ApplyTelluricWaveCorrection =  NO " << endl;
            }
            
            if (applyHeliocentricRVCorrection) {
                fout << "# ApplyHeliocentricRVCorrection =  YES " << endl;
            } else {
                fout << "# ApplyHeliocentricRVCorrection =  NO " << endl;
            }
            
            fout << "# Stack of " << inputIndex << " spectra:" << endl;

            for (unsigned index=0; index<inputIndex; index++) {
                fout << "# " << index  << " " << inputspectra[index] <<  endl;
            }
            fout << "# NumberOfDataPoints = " << outspectrum.size() << endl;

            fout << "# Data format is: " << endl;
            fout << "# <wavelength(nm)> <flux> <variance>" << endl;
            
            for (unsigned i=0; i<outspectrum.size(); i++) {
                fout << outspectrum.getwavelength(i) << " " << outspectrum.getflux(i) << " " << outspectrum.getvariance(i) << endl;
            }
            fout.close();
        }
        

    }
	catch (operaException e) {
		cerr << "operaStackObjectSpectra: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaStackObjectSpectra: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 
