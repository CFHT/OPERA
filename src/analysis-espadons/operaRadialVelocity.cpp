/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaRadialVelocity
 Version: 1.0
 Description: Measure radial velocity shift from source
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope
 Location: Hawaii USA
 Date: Feb/2016
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
#include "analysis-espadons/operaRadialVelocity.h"


#define MAXNUMBEROFLINESWITHINCHUNK 1000
#define MINSIZEOFSPECTRALCHUNK 400

/*! \file operaRadialVelocity.cpp */

using namespace std;

operaVector generateSyntheticTelluricSpectrumUsingLineProfile(const operaSpectrum& telluricLines, const operaVector& wavelengthVector, double resolution, ProfileMethod profile);
bool calculateRVShiftByXCorr(const operaSpectrum& objectSpectrum, const operaSpectrum& templateSpectrum, const operaSpectrum& telluricLines, double radialVelocityRange, double radialVelocityStep, double threshold, double& maxRV, double& sigRV, double& maxcorr, ofstream& frvcorrdata, ofstream& frvcorrfitdata, double spectralResolution, bool useFitToFindMaximum, double& chisqr, double heliocentricRV_mps);
void getWavelengthSubrange(const operaVector& wavelength, double wl0, double wlf, unsigned& startindex, unsigned& endindex);

void matchTelluricLines(const operaSpectrum& telluricLinesFromAtlas, const operaSpectrum& telluricLinesFromObject, operaVector& telluricMatchedWavelengths, operaSpectrum& objectMatchedLines, operaVector& radialVelocities, double spectralResolution, double radialVelocityRange);
operaVector selectWavelengthsInRange(operaWavelengthRange wlrange, operaVector wavelengths);

operaSpectrum readTemplateSpectrum(string inputTemplateSpectrum);
operaSpectrum readSpectrumFromSPCFile(string inputSpectrumFileName, int ordernumber, int minorder, int maxorder, bool applyTelluricWaveCorrection, bool applyHeliocentricRVCorrection);

operaArgumentHandler args;

/*!
 * operaRadialVelocity
 * \author Eder Martioli
 * \brief Measure radial velocity shift from source.
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
    
    string inputObjectSpectrum;

    string outputRVFile;
    
    string telluric_lines; // HITRAN Library
    
    string template_spectrum; // Template spectrum of the source
    
    string inputWavelengthRangesForRVMeasurements; // Wavelength ranges for RV measurements
    string inputHeliocentricCorrection;

    bool printIndividualMeasurements = false;
    
    int ordernumber = NOTPROVIDED;
    int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;

    double mjdate = 0.0;
    double spectralResolution = 80000;

    double radialVelocityRange = 150;
    double radialVelocityStep = 0.3;
    double threshold = 0.05;
    bool useFitToFindMaximum = false;
    
    
    args.AddRequiredArgument("inputObjectSpectrum", inputObjectSpectrum, "input object spectrum file (.e or .p)");
    args.AddRequiredArgument("outputRVFile", outputRVFile, "output RV file (.rv)");
    args.AddRequiredArgument("telluric_lines", telluric_lines, "atlas of telluric lines (HITRAN)");
    args.AddOptionalArgument("template_spectrum", template_spectrum, "", "template spectrum of the source");
    args.AddRequiredArgument("inputWavelengthRangesForRVMeasurements", inputWavelengthRangesForRVMeasurements, "wavelength ranges for RV measurements");
    args.AddRequiredArgument("inputHeliocentricCorrection", inputHeliocentricCorrection, "radial velocity heliocentric corrections");
    args.AddSwitch("printIndividualMeasurements", printIndividualMeasurements,"print out individual measurements otherwise print collapsed RV");
    args.AddRequiredArgument("mjdate", mjdate, "input Modified Julian Date (MJD) at start of exposure");
    args.AddRequiredArgument("spectralResolution", spectralResolution, "input spectral resolution (wl/dwl) as reference for line detection");
    
    args.AddOptionalArgument("radialVelocityRange", radialVelocityRange, 150, "Radial velocity search range in km/s");
    args.AddOptionalArgument("radialVelocityStep", radialVelocityStep, 0.3, "Radial velocity search step in km/s");
    args.AddOptionalArgument("threshold", threshold, 0.05, "Cross-correlation threshold (must lie between -1.0 and 1.0)");
    args.AddSwitch("useFitToFindMaximum", useFitToFindMaximum,"Activate if want to use gaussian fit to measure maximum x-correlation");

    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    
    try {
        args.Parse(argc, argv);
        
        // we need an input object spectrum file ...
        if (inputObjectSpectrum.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        // we need an output wavelength calibration file ...
        if (outputRVFile.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        // we need an input atlas of telluric lines ...
        if (telluric_lines.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        // we need a template spectrum of the source ...
        if (template_spectrum.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        if (inputWavelengthRangesForRVMeasurements.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        if (args.verbose) {
            cout << "operaRadialVelocity: inputObjectSpectrum = " << inputObjectSpectrum << endl;
            cout << "operaRadialVelocity: outputRVFile = " << outputRVFile << endl;
            cout << "operaRadialVelocity: telluric_lines =" << telluric_lines << endl;
            cout << "operaRadialVelocity: template_spectrum =" << template_spectrum  << endl;
            cout << "operaRadialVelocity: mjdate =" << mjdate << endl;
            cout << "operaRadialVelocity: spectralResolution =" << spectralResolution << endl;
            cout << "operaRadialVelocity: inputWavelengthRangesForRVMeasurements = " << inputWavelengthRangesForRVMeasurements << endl;
            cout << "operaRadialVelocity: inputHeliocentricCorrection = " << inputHeliocentricCorrection << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaRadialVelocity: ordernumber = " << ordernumber << endl;
        }
        
        operaSpectrum inputSpectrum = readSpectrumFromSPCFile(inputObjectSpectrum,ordernumber,minorder,maxorder, true, false);

        double heliocentricRV_mps = 0.0;
        double HJD_UTC = 0.0;
        double HJD_TT = 0.0;
        /*
         * Load Heliocentric RV wavelength correction
         */
        if (!inputHeliocentricCorrection.empty()) {
            
            FormatData rveldata;
            operaIOFormats::ReadCustomFormat("rvel", rveldata, inputHeliocentricCorrection);

            double heliocentricRV_kps = rveldata.extract<double>();
            double lunar_rvel = rveldata.extract<double>();
            double orbital_rvel = rveldata.extract<double>();
            double diurnal_rvel = rveldata.extract<double>();
            HJD_UTC = rveldata.extract<double>();
            HJD_TT = rveldata.extract<double>();
            
            heliocentricRV_mps = heliocentricRV_kps*1000;
        }
        
        /*
         * Load telluric lines database lambda vs. intensity
         */
        operaSpectrum telluricLines = readTelluricLines(telluric_lines);
        //---

        /*
         * Load template spectrum: lambda vs. intensity
         */
        operaSpectrum templateSpectrum = readTemplateSpectrum(template_spectrum);


        if (args.debug) {
            for (unsigned i=0; i<templateSpectrum.size(); i++) {
                cout << templateSpectrum.getwavelength(i) << " " << templateSpectrum.getflux(i) << " " << templateSpectrum.getvariance(i) << endl;
            }
        }
        
        /*
         * Load selected spectral ranges within which will be performed RV measurements
         */
        operaWavelengthRanges wlranges = readContinuumWavelengthMask(inputWavelengthRangesForRVMeasurements);
        //---
        
        float *radialvelocities = new float[wlranges.size()];
        float *radialvelocitiesErr = new float[wlranges.size()];
        double *wl0 = new double[wlranges.size()];
        double *wlf = new double[wlranges.size()];
        unsigned nRVs = 0;
        
        /*
         *   Below it performs RV measurements for each spectral chunk separately
         */
        for (unsigned chunk=0; chunk<wlranges.size(); chunk++) {
            
            double maxcorr = 0;
            double chisqr = 0;
            ofstream frvcorrdata;
            ofstream frvcorrfitdata;
            
            double rvshift = 0;
            double rvshifterror = 0;
            
            if (args.verbose) cout << "operaRadialVelocity: Processing chunk #" << chunk << " wl0=" << wlranges.wl0(chunk) << " wlf=" << wlranges.wlf(chunk) << endl;
            
            /*
             *  Collect observed spectrum within chunk
             */
            operaSpectrum spectrumChunk = getSpectrumWithinRange(wlranges.getrange(chunk),inputSpectrum);

            if (spectrumChunk.size() < MINSIZEOFSPECTRALCHUNK) {
                if (args.verbose) { cout << "operaRadialVelocity: skipping chunk #" << chunk << "  chunk size = " << spectrumChunk.size() << " < " << MINSIZEOFSPECTRALCHUNK << endl; }
                continue;
            }
            
            if(args.debug){
                for (unsigned l=0; l<spectrumChunk.size(); l++) {
                    cout << spectrumChunk.getwavelength(l) << " " << spectrumChunk.getflux(l) << " " << spectrumChunk.getvariance(l) << endl;
                }
            }

            /*
             *  Collect template spectrum within chunk
             *  Get wider range to allow RV search
             */
            double dlambda = (radialVelocityRange/2.0)*spectrumChunk.midwl()/SPEED_OF_LIGHT_KMS;
            operaWavelengthRange widerRange(spectrumChunk.firstwl()-dlambda,spectrumChunk.lastwl()+dlambda);
            
            operaSpectrum templateChunk = getSpectrumWithinRange(widerRange,templateSpectrum);

            /*
             *  Collect telluric lines within chunk
             */
            operaSpectrum telluricChunk = getSpectrumWithinRange(wlranges.getrange(chunk),telluricLines);

            bool xcorrect = calculateRVShiftByXCorr(spectrumChunk,templateChunk,telluricChunk,radialVelocityRange,radialVelocityStep,threshold,rvshift,rvshifterror,maxcorr,frvcorrdata,frvcorrfitdata,spectralResolution,useFitToFindMaximum,chisqr,heliocentricRV_mps);

            if(args.debug)
                cout << chunk << " " << spectrumChunk.firstwl() << " " << spectrumChunk.lastwl() << " " << rvshift << " " << rvshifterror << " " << maxcorr << endl;
            
            // The output radial velocity is given in m/s
            radialvelocities[nRVs] = (float)rvshift*1000;
            radialvelocitiesErr[nRVs] = (float)rvshifterror*1000;
            wl0[nRVs] = spectrumChunk.firstwl();
            wlf[nRVs] = spectrumChunk.lastwl() ;
            nRVs++;
        }
        
        FormatHeader outputheader("Source Radial Velocity");
        
        if (printIndividualMeasurements) {
            outputheader << "MJD" << "HJD_UTC" << "HJD_TT" << "wl0(nm)" << "wlf(nm)" <<"radialvelocity (m/s)" << "radialvelocityerror (m/s)" << newline;
            FormatData outputdata;

            for (unsigned i=0; i<nRVs; i++) {
                outputdata << fixed << setprecision(10) << mjdate << HJD_UTC << HJD_TT << fixed << setprecision(5) << wl0[i] << wlf[i] << radialvelocities[i] << radialvelocitiesErr[i] << endl;
            }
            operaIOFormats::WriteCustomFormat("rv", outputheader, outputdata, outputRVFile);
        } else {
            outputheader << "MJD" << "HJD_UTC" << "HJD_TT" <<"radialvelocity (m/s)" << "radialvelocityerror (m/s)" << newline;
            FormatData outputdata;

            float medianRV = operaArrayMedian(nRVs,radialvelocities);
            float medianSigmaRV = operaArrayMedianSigma(nRVs,radialvelocities,medianRV);
            
            outputdata << fixed << setprecision(10) << mjdate << HJD_UTC << HJD_TT << fixed << setprecision(5) << medianRV << medianSigmaRV << endl;
            operaIOFormats::WriteCustomFormat("rv", outputheader, outputdata, outputRVFile);
        }
        
        
    }
    catch (operaException e) {
        cerr << "operaRadialVelocity: " << e.getFormattedMessage() << endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        cerr << "operaRadialVelocity: " << operaStrError(errno) << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

// Read the entire set of telluric lines in HITRAN database
operaSpectrum readTelluricLines(string telluric_database_file)
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

// Function to match telluric lines
void matchTelluricLines(const operaSpectrum& telluricLinesFromAtlas, const operaSpectrum& telluricLinesFromObject, operaVector& telluricMatchedWavelengths, operaSpectrum& objectMatchedLines, operaVector& radialVelocities, double spectralResolution, double radialVelocityRange)
{
    unsigned initiaAtlasIndex = 0;
    
    for (unsigned j=0; j<telluricLinesFromObject.size(); j++) {
        double linewidth = telluricLinesFromObject.getwavelength(j)/spectralResolution;
        
        double diffmin = BIG;
        bool thereIsAMatch = false;
        unsigned matched_jindex = 0;
        unsigned matched_iindex = 0;
        
        for (unsigned i=initiaAtlasIndex; i<telluricLinesFromAtlas.size(); i++) {
            double wl_diff = fabs(telluricLinesFromObject.getwavelength(j) - telluricLinesFromAtlas.getwavelength(i));
            
            if(wl_diff < linewidth && wl_diff < diffmin) {
                thereIsAMatch = true;
                matched_jindex = j;
                matched_iindex = i;
                diffmin = wl_diff;
            } else {
                if (telluricLinesFromAtlas.getwavelength(i) > telluricLinesFromObject.getwavelength(j) + 2*linewidth) {
                    if (thereIsAMatch) {
                        double atlaswl = telluricLinesFromAtlas.getwavelength(matched_iindex);
                        double objectwl = telluricLinesFromObject.getwavelength(matched_jindex);
                        double objectflux = telluricLinesFromObject.getflux(matched_jindex);
                        double deltarv = calculateDeltaRadialVelocityInKPS(atlaswl, objectwl);
                        if(fabs(deltarv) < radialVelocityRange/2.0) {
                            telluricMatchedWavelengths.insert(atlaswl);
                            objectMatchedLines.insert(objectwl, objectflux);
                            radialVelocities.insert(deltarv);
                        }
                        // cout << nmatches << " " << telluricLinesFromObject.wavelength[matched_jindex] << " " <<  telluricLinesFromAtlas.wavelength[matched_iindex] << " " << radialVelocities[nmatches] << endl;
                    }
                    initiaAtlasIndex = i+1;
                    break;
                }
            }
        }
    }
}

operaVector selectWavelengthsInRange(operaWavelengthRange wlrange, operaVector wavelengths) {
    operaVector selectedWavelengths;
    for (unsigned l=0; l<wavelengths.size();l++) {
        if (wlrange.contains(wavelengths[l])) {
            selectedWavelengths.insert(wavelengths[l]);
        }
    }
    return selectedWavelengths;
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

/*
 * Read template spectrum from file. Data format must be lambda, flux, fluxvar
 */
operaSpectrum readTemplateSpectrum(string inputTemplateSpectrum) {

    operaSpectrum templateSpectrum;

    double wl = -1.0;
    double flux = -1.0;
    double fluxvar = -1.0;

    operaistream fspec(inputTemplateSpectrum.c_str());
    
    if (fspec.is_open()) {
        string dataline;
        
        while (fspec.good()) {
            getline(fspec, dataline);
            if (strlen(dataline.c_str())) {
                stringstream ss (stringstream::in | stringstream::out);
                ss << dataline.c_str();
                
                if (dataline.c_str()[0] == '#') { 	// skip comments
                    
                } else {
                    ss >> wl;
                    ss >> flux;
                    ss >> fluxvar;
                    
                    templateSpectrum.insert(wl,flux,fluxvar);
                }

            }
        }
        fspec.close();
    }
    return templateSpectrum;
}

/*
 * Read spectrum from OPERA spc product.
 */
operaSpectrum readSpectrumFromSPCFile(string inputSpectrumFileName, int ordernumber, int minorder, int maxorder, bool applyTelluricWaveCorrection, bool applyHeliocentricRVCorrection) {
    
    operaSpectrum outputSpectrum;
    
    operaSpectralOrderVector spectralOrderVector;
    operaIOFormats::ReadIntoSpectralOrders(spectralOrderVector, inputSpectrumFileName);
    UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrderVector);

    unsigned spectrumTypeToExtract = 1; // 1=Normalized, 2=FluxCalibrated, 0=Default=Raw
    double snrClip = 1; // signal-to-noise threshold
    unsigned numberOfPointsToCutInOrderEnds = 0;
    double RV_KPS = 0.0;  // source radial velocity to subtract
    
    outputSpectrum = spectralOrderVector.getExtendedSpectrum(minorder,maxorder,spectrumTypeToExtract,applyTelluricWaveCorrection,applyHeliocentricRVCorrection,snrClip,numberOfPointsToCutInOrderEnds,RV_KPS);
    
    return outputSpectrum;
}


bool calculateRVShiftByXCorr(const operaSpectrum& objectSpectrum, const operaSpectrum& templateSpectrum, const operaSpectrum& telluricLines, double radialVelocityRange, double radialVelocityStep, double threshold, double& maxRV, double& sigRV, double& maxcorr, ofstream& frvcorrdata, ofstream& frvcorrfitdata, double spectralResolution, bool useFitToFindMaximum, double& chisqr, double heliocentricRV_mps)
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

    operaVector templateIntensityVector = convolveSpectrum(templateSpectrum, spectralResolution);
    //operaVector templateIntensityVector = (templateSpectrum.getintensity()).getflux();
    
    for(double deltaRV = -radialVelocityRange/2.0; deltaRV <= radialVelocityRange/2.0; deltaRV+=radialVelocityStep) {
        operaVector syntheticSpectrumWavelength;
        // Initalize synthetic wavelength with wavelength of objectSpectrum shifted by deltaRV
        for (unsigned i=0; i<objectSpectrum.size(); i++) {
            double DWavelength = (deltaRV + heliocentricRV_mps/1000) * objectSpectrum.getwavelength(i) / SPEED_OF_LIGHT_KMS;
            syntheticSpectrumWavelength.insert(objectSpectrum.getwavelength(i) - DWavelength);
        }
        
        // Generate a spectrum in telluricSpectrumFlux along points in wavelength vector using the provided telluricLines
        operaVector telluricSpectrumFlux = generateSyntheticTelluricSpectrumUsingLineProfile(telluricLines, objectSpectrum.wavelengthvector(), spectralResolution, GAUSSIAN);
        
        // Generate a spectrum in templateSpectrumFlux along points in wavelength vector using the provided templateSpectrum
        operaVector templateSpectrumFlux = fitSpectrum(templateSpectrum.wavelengthvector(), templateIntensityVector, syntheticSpectrumWavelength);
        
        // Generate a sythetic spectrum consisting of telluricSpectrumFlux * templateSpectrumFlux
        operaVector syntheticFlux(syntheticSpectrumWavelength.size());
        
        for(unsigned i=0; i<syntheticSpectrumWavelength.size(); i++) {
            syntheticFlux[i] = telluricSpectrumFlux[i] * templateSpectrumFlux[i];
            if(args.debug)
                cout << syntheticSpectrumWavelength[i] << " " << objectSpectrum.getflux(i) << " " << telluricSpectrumFlux[i] << " " << templateSpectrumFlux[i] << " " << syntheticFlux[i] << endl;
        }

        // Calculate the x-corr between the generated shifted synthetic spectrum and the object spectrum
        double xcorr = operaCrossCorrelation(syntheticSpectrumWavelength.size(), objectSpectrum.flux_ptr(), syntheticFlux.datapointer());

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
        unsigned HalfNumberOfPointsToUseInFit = 5;

        operaVector crosscorrelation_s;
        operaVector crosscorrerror_s;
        operaVector dRV_s;
        
        unsigned firstj = 0;
        unsigned lastj = 0;
        if ((int)jmax - (int)HalfNumberOfPointsToUseInFit < 0) {
            firstj = 0;
        } else {
            firstj = jmax - HalfNumberOfPointsToUseInFit;
        }
        
        if ((int)jmax + (int)HalfNumberOfPointsToUseInFit >= crosscorrelation.size()) {
            lastj = crosscorrelation.size();
        } else {
            lastj = jmax + HalfNumberOfPointsToUseInFit;
        }
        
        for (unsigned j=firstj; j<lastj; j++) {
            crosscorrelation_s.insert(crosscorrelation[j]);
            crosscorrerror_s.insert(crosscorrerror[j]);
            dRV_s.insert(dRV[j]);
        }
        
        // Set initial values for our Gaussian using the maximum x-corr.
        double a = crosscorrelation[jmax]; //Initial amplitude
        double x0 = dRV[jmax]; //Initial center
        double sig = radialVelocityRange/4.0; //Initial sigma
        double ea;
        double ex0;
        double esig;
        double fitchisqr;
        
        // Update the initial values and get errors for each along with the fit chi-squared.
        operaMPFitGaussian(crosscorrelation_s.size(), dRV_s.datapointer(), crosscorrelation_s.datapointer(), crosscorrerror_s.datapointer(), &a, &ea, &x0, &ex0, &sig, &esig, &fitchisqr);
        //operaLMFitGaussian(np, peakXdata, peakYdata, &a, &x0, &sig, &chisqr);
        
        if(args.debug) {
            cout << a << "+/-" << ea << endl;
            cout << x0 << "+/-" << ex0 << endl;
            cout << sig << "+/-" << esig <<  " fitchisqr=" << fitchisqr << endl;
        }
        
        // For plotting
        if(frvcorrdata.is_open()) {
            for(unsigned j=0; j<crosscorrelation_s.size(); j++) {
                double x = (double)dRV[j];
                double gaussfunc = a*exp(-(x-x0)*(x-x0)/(2*sig*sig));
                frvcorrdata << dRV_s[j] << " " <<  gaussfunc << " " <<  crosscorrelation_s[j] << " " <<  crosscorrerror_s[j] << " " << crosscorrelation_s[j] - gaussfunc << endl;
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

