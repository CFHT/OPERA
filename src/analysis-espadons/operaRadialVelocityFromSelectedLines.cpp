/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaRadialVelocityFromSelectedLines
 Version: 1.0
 Description: Measure radial velocity shift from source
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

/*! \file operaRadialVelocityFromSelectedLines.cpp */

using namespace std;

void writegaussfitenv(string envfile, string datafile, string paramfile, string rvparamfile, bool robust);
void writegaussfitdata(string datafile, operaWavelengthRanges linesranges, operaSpectrum spectrum, operaSpectrum sourceLines);
void writegaussfitmodel(string modelfile);
void writegaussfitpars(string paramfile, operaSpectrum lines);
void writegaussfitrvpars(string rvparamfile, double telluricRV, double helioRV, double sourceRV);
void readgaussfitrvpars(string rvparamfile, double *rv, double *rverr);
void GenerateSelectedLinesPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafile, operaSpectrum sourceLines);

operaArgumentHandler args;

/*! 
 * operaRadialVelocityFromSelectedLines
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
	string inputWaveFile;
    string inputObjectSpectrum;
    string inputFlatFluxCalibration;
	string outputRVFile;
    string telluric_lines; // HITRAN Library
    string source_lines; // User Library of stellar lines
    string inputWavelengthMaskForTelluric;
    string inputHeliocentricCorrection;
    string inputTelluricCorrection;
    
	int ordernumber = NOTPROVIDED;
    int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;
    bool StarPlusSky = false;

    double mjdate = 0.0;

    // The parameters below we don't know yet whether they would be useful if used as input
    double spectralResolution = 80000;
    unsigned normalizationBinsize = 110;
    
	double LocalMaxFilterWidth=4.0;
	double MinPeakDepth=3;
	double DetectionThreshold=0.05;
	double nsigclip=3.0;
	bool emissionSpectrum = false;
	unsigned minNumberOfMatchedLines = 10;
    
    double initialRVguess=0.0;

    bool robustFit = false;

    double sourceLineWidthResolution = 2500;
    
    string outputPlotEPSFileName;
    string gnuScriptFileName;

    args.AddRequiredArgument("inputWaveFile", inputWaveFile, "input wavelength calibration file (.wcal)");
	args.AddRequiredArgument("inputObjectSpectrum", inputObjectSpectrum, "input object spectrum file (.e or .p)");
    args.AddOptionalArgument("inputFlatFluxCalibration", inputFlatFluxCalibration, "", "flat field spectrum ff_");
    args.AddRequiredArgument("outputRVFile", outputRVFile, "output telluric wavelength calibration file (.tell)");
    args.AddRequiredArgument("telluric_lines", telluric_lines, "atlas of telluric lines (HITRAN)");
    args.AddRequiredArgument("source_lines", source_lines, "library of source (star/nebular/HII) lines");
    args.AddRequiredArgument("inputWavelengthMaskForTelluric", inputWavelengthMaskForTelluric, "telluric wavelength mask");
    args.AddRequiredArgument("inputHeliocentricCorrection", inputHeliocentricCorrection, "radial velocity heliocentric corrections");
    args.AddRequiredArgument("inputTelluricCorrection", inputTelluricCorrection, "radial velocity telluric correction");
    
    args.AddRequiredArgument("mjdate", mjdate, "input Modified Julian Date (MJD) at start of exposure");
    args.AddRequiredArgument("spectralResolution", spectralResolution, "input spectral resolution (wl/dwl) as reference for line detection");
    args.AddRequiredArgument("normalizationBinsize", normalizationBinsize, "binsize to normalize input object spectrum");
    
    args.AddOptionalArgument("LocalMaxFilterWidth", LocalMaxFilterWidth, 4.0, "Used for line matching");
    args.AddOptionalArgument("MinPeakDepth", MinPeakDepth, 3, "Used for line matching");
    args.AddOptionalArgument("DetectionThreshold", DetectionThreshold, 0.05, "Used for line matching");
    args.AddOptionalArgument("nsigclip", nsigclip, 3.0, "Used for line matching");
    args.AddSwitch("emissionSpectrum", emissionSpectrum, "Used for line matching");

    args.AddOptionalArgument("minNumberOfMatchedLines", minNumberOfMatchedLines, 10, "Arbitrary threshold to avoid small number statistics");
    
    args.AddOptionalArgument("initialRVguess", initialRVguess, 0.0, "An initial guess for the source radial velocity in km/s");

    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    args.AddSwitch("StarPlusSky", StarPlusSky, "star plus sky mode");
    args.AddSwitch("robustFit", robustFit, "Used for robust fit");

    args.AddOptionalArgument("sourceLineWidthResolution", sourceLineWidthResolution, 2500, "Source line resolution -- usually much smaller than instrument resolution");

    
    args.AddOptionalArgument("outputPlotEPSFileName", outputPlotEPSFileName, "", "Output plot eps file name");
    args.AddOptionalArgument("gnuScriptFileName", gnuScriptFileName, "", "Output gnuplot script file name");

	try {
		args.Parse(argc, argv);
		
		// we need an input wavelength calibration file ...
		if (inputWaveFile.empty()) {
			throw operaException("operaRadialVelocityFromSelectedLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input object spectrum file ...        
		if (inputObjectSpectrum.empty()) {
			throw operaException("operaRadialVelocityFromSelectedLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output wavelength calibration file ...
		if (outputRVFile.empty()) {
			throw operaException("operaRadialVelocityFromSelectedLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        // we need an input atlas of telluric lines ...
        if (telluric_lines.empty()) {
            throw operaException("operaRadialVelocityFromSelectedLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        // we need an input atlas of stellar lines ...
        if (source_lines.empty()) {
            throw operaException("operaRadialVelocityFromSelectedLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        if (inputWavelengthMaskForTelluric.empty()) {
            throw operaException("operaRadialVelocityFromSelectedLines: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
		if (args.verbose) {
			cout << "operaRadialVelocityFromSelectedLines: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaRadialVelocityFromSelectedLines: inputObjectSpectrum = " << inputObjectSpectrum << endl;
            cout << "operaRadialVelocityFromSelectedLines: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
			cout << "operaRadialVelocityFromSelectedLines: outputRVFile = " << outputRVFile << endl;
            cout << "operaRadialVelocityFromSelectedLines: telluric_lines =" << telluric_lines << endl;
            cout << "operaRadialVelocityFromSelectedLines: source_lines =" << source_lines << endl;
            cout << "operaRadialVelocityFromSelectedLines: mjdate =" << mjdate << endl;
            cout << "operaRadialVelocityFromSelectedLines: spectralResolution =" << spectralResolution << endl;
            cout << "operaRadialVelocityFromSelectedLines: normalizationBinsize =" << normalizationBinsize << endl;
            cout << "operaRadialVelocityFromSelectedLines: StarPlusSky = " << StarPlusSky << endl;
            cout << "operaRadialVelocityFromSelectedLines: inputWavelengthMaskForTelluric = " << inputWavelengthMaskForTelluric << endl;
            cout << "operaRadialVelocityFromSelectedLines: inputHeliocentricCorrection = " << inputHeliocentricCorrection << endl;
            cout << "operaRadialVelocityFromSelectedLines: inputTelluricCorrection = " << inputTelluricCorrection << endl;
            cout << "operaRadialVelocityFromSelectedLines: initialRVguess = " << initialRVguess << endl;
            cout << "operaRadialVelocityFromSelectedLines: robustFit = " << robustFit << endl;
            cout << "operaRadialVelocityFromSelectedLines: sourceLineWidthResolution = " << sourceLineWidthResolution << endl;            
            if(ordernumber != NOTPROVIDED) cout << "operaRadialVelocityFromSelectedLines: ordernumber = " << ordernumber << endl;
		}
		
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputObjectSpectrum);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputWaveFile);

        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        if (args.verbose) cout << "operaRadialVelocityFromSelectedLines: minorder ="<< minorder << " maxorder=" << maxorder << endl;
    
        // Correct for flat-field
        if (!inputFlatFluxCalibration.empty()) {
			operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputFlatFluxCalibration);
            bool starplusskyInvertSkyFiber = false;
            spectralOrders.correctFlatField(minorder, maxorder, StarPlusSky, starplusskyInvertSkyFiber);
        }
        
        double heliocentricRV_mps = 0.0;
        double HJD_UTC = 0.0;
        double HJD_TT = 0.0;

        // Load Heliocentric RV wavelength correction
        if (!inputHeliocentricCorrection.empty()) {
            FormatData rveldata;
            operaIOFormats::ReadCustomFormat("rvel", rveldata, inputHeliocentricCorrection);
            double heliocentricRV = rveldata.extract<double>();
            double lunar_rvel = rveldata.extract<double>();
            double orbital_rvel = rveldata.extract<double>();
            double diurnal_rvel = rveldata.extract<double>();
            HJD_UTC = rveldata.extract<double>();
            HJD_TT = rveldata.extract<double>();
            heliocentricRV_mps = heliocentricRV*1000;
        }
        
        double telluricRV_mps = 0.0;
        // Load telluric correction for wavelength calibration
        if (!inputTelluricCorrection.empty()) {
            FormatData telldata;
            operaIOFormats::ReadCustomFormat("tell", telldata, inputTelluricCorrection);
            double telluricRV = telldata.extract<double>();
            telluricRV_mps = telluricRV*1000;
        }

        // Read telluric lines database lambda vs. intensity
		if (args.debug) cout << "operaRadialVelocityFromSelectedLines: reading telluric lines database " << telluric_lines << endl;
		operaSpectrum telluricLines = readTelluricLines(telluric_lines);
        
        
        // Read stellar lines database lambda vs. intensity
        if (args.debug) cout << "operaRadialVelocityFromSelectedLines: reading source lines database " << source_lines << endl;
        operaSpectrum sourceLines = readInputListOfLines(source_lines);
        
        // Initialize rvshift to zero, so if telluric SNR is low the RV shift will not do anything.
        double rvshift = 0;
        double rvshifterror = 0;
        
        /* 
         *   Below we measure the RV shift :
         */
        
        // Detect absorption lines in the observed spectrum within telluric regions defined in inputWavelengthMaskForTelluric
        operaSpectralLineList objectLines = spectralOrders.detectSpectralLinesWithinWavelengthMask(inputWavelengthMaskForTelluric, minorder, maxorder,true, normalizationBinsize,spectralResolution,emissionSpectrum,LocalMaxFilterWidth,MinPeakDepth,DetectionThreshold,nsigclip);
        if(args.debug){
            for (unsigned l=0; l<objectLines.size(); l++) {
                cout << objectLines.getwavelength(l) << " " << objectLines.getflux(l) << " " << objectLines.getsigma(l) << endl;
            }
        }
        
        // Detect absorption source lines in the observed spectrum
        double snrClip = 1;
        double numberOfPointsToCutInOrderEnds = 400;
        double nsigwidth = 5.0;  // extract +- nsig around line
        
        operaWavelengthRanges linesranges = getWavelengthMaskAroundLines(sourceLines,sourceLineWidthResolution,nsigwidth);

        operaSpectrum spectrumAroundLines = spectralOrders.getSpectrumAroundLines(sourceLines, minorder, maxorder, true, normalizationBinsize, sourceLineWidthResolution, nsigwidth, snrClip, numberOfPointsToCutInOrderEnds);

        if(args.debug){
            for (unsigned i=0; i<spectrumAroundLines.size(); i++) {
                cout << spectrumAroundLines.getflux(i) << " " << spectrumAroundLines.getvariance(i) << endl;
            }
        }
    
        /*
         * GaussFit files:
         */
        string envfile = "envfile";
        string datafile = "datafile";
        string modelfile = "modelfile";
        string paramfile = "paramfile";
        string rvparamfile = "rvparamfile";
        
        // -> Create command line to run gaussfit
        string command = "/Users/edermartioli/Local/gaussfit/gaussfit " + modelfile + " " + envfile;
        // -> Print command line on screen
        if(args.verbose) cout << endl << command << endl << endl;

        /*
         * First iteration to fit calibration parameters
         */
        bool robustFitForFirstPass = false;
        
        // -> Write gaussfit environment file
        writegaussfitenv(envfile, datafile, paramfile, rvparamfile, robustFitForFirstPass);
        
        // -> Write gaussfit data file
        writegaussfitdata(datafile, linesranges, spectrumAroundLines, sourceLines);
        
        // -> Write gaussfit model file
        writegaussfitmodel(modelfile);
        
        // -> Write gaussfit parameter files
        writegaussfitpars(paramfile, sourceLines);
        
        double initialRVguess_mps = initialRVguess*1000;
        writegaussfitrvpars(rvparamfile, telluricRV_mps, heliocentricRV_mps, initialRVguess_mps);
        
        //***********************
        // -> Run gaussfit
        system(command.c_str());
        //***********************
    
        if (robustFit) {
            writegaussfitenv(envfile, datafile, paramfile, rvparamfile, robustFit);
            
            //***********************
            // -> Run gaussfit
            system(command.c_str());
            //***********************
        }
        readgaussfitrvpars(rvparamfile, &rvshift, &rvshifterror);
        
        FormatHeader outputheader("Source Radial Velocity");
        outputheader << "MJD" << "HJD_UTC" << "HJD_TT" << "radialvelocity (m/s)" << "radialvelocityerror (m/s)" << newline;
        FormatData outputdata;
        outputdata << fixed << setprecision(10) << mjdate << HJD_UTC << HJD_TT << fixed << setprecision(5)<< rvshift << rvshifterror << endl;
        operaIOFormats::WriteCustomFormat("rv", outputheader, outputdata, outputRVFile);
        
        if (!gnuScriptFileName.empty()) {
            GenerateSelectedLinesPlot(gnuScriptFileName,outputPlotEPSFileName,datafile,sourceLines);
        }
        
    }
	catch (operaException e) {
		cerr << "operaRadialVelocityFromSelectedLines: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaRadialVelocityFromSelectedLines: " << operaStrError(errno) << endl;
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

operaSpectrum readInputListOfLines(string inputListOfLines) {
    
    operaistream file(inputListOfLines.c_str());
    
    operaSpectrum spectralLines;
    
    if (file.is_open()) {
        unsigned line = 0;
        string dataline;
        
        double wavelength_in_nm;
        double intensity;
        
        while (file.good()) {
            getline(file, dataline);
            if (strlen(dataline.c_str())) {
                stringstream ss (stringstream::in | stringstream::out);
                ss << dataline.c_str();
        
                if (dataline[0] == '#') {
                    line++;
                } else {
                    
                    ss >> wavelength_in_nm;
                    ss >> intensity;
                    
                    spectralLines.insert(wavelength_in_nm,intensity);
//                    cout << wavelength_in_nm << " " << intensity << endl;
                    
                    line++;
                }
            }
        }
        file.close();
        return spectralLines;
    } else {
        return EXIT_FAILURE;
    }
}


void writegaussfitenv(string envfile, string datafile, string paramfile, string rvparamfile, bool robust) {
    /*
     * Create env file
     */
    ofstream *fenv = NULL;
    
    if (!envfile.empty()) {
        fenv = new ofstream();
        fenv->open(envfile.c_str());
        
        *fenv << "data1 = '" << datafile << "'" << endl;
        *fenv << "param1 = '" << rvparamfile << "'" << endl;
        *fenv << "param2 = '" << paramfile << "'" << endl;
        
        *fenv <<"iters   =	 30.0" << endl;
        *fenv <<"prmat   =	 0.0" << endl;
        *fenv <<"prvar   =	 1.0" << endl;
        *fenv <<"results =	'res'" << endl;
        *fenv <<"tol     =	 0.0001" << endl;
        *fenv <<"triang  =	 0.0" << endl;
        
        if(robust == true)
            *fenv << "fair   =	 0.9" << endl;
        
        *fenv <<"END" << endl;
        
        fenv->close();
    }
}


void writegaussfitdata(string datafile, operaWavelengthRanges linesranges, operaSpectrum spectrum, operaSpectrum sourceLines) {
    
    ofstream *fdata = NULL;
    
    if (!datafile.empty()) {
        fdata = new ofstream();
        fdata->open(datafile.c_str());
        
        // write data header containing data symbols
        *fdata << "line wlc wavelength flux flux_flux" << endl;
        //----
        
        // write data type (double, float, int, etc)
        *fdata << "double double double double double" << endl;
        //----
        
        // write data
        for (unsigned l=0; l<linesranges.size(); l++) {
            for (unsigned i=0; i<spectrum.size(); i++) {
                if (linesranges.contains(spectrum.getwavelength(i),l) ) {
                    *fdata << l << " " << sourceLines.getwavelength(l) << " " << spectrum.getwavelength(i) << " " << spectrum.getflux(i) << " " << spectrum.getvariance(i) << endl;
                }
            }
        }
        
        fdata->close();
    }
    
}

void writegaussfitmodel(string modelfile) {
    
    ofstream *fmodel = NULL;
    
    if (!modelfile.empty()) {
        fmodel = new ofstream();
        fmodel->open(modelfile.c_str());
        
        *fmodel << "parameter rv;" << endl;
        *fmodel << "constant telluricRV;" << endl;
        *fmodel << "constant helioRV;" << endl;
        *fmodel << "parameter level[line];" << endl;
        *fmodel << "parameter a[line];" << endl;
        *fmodel << "parameter sig[line];" << endl;
        *fmodel << "constant wl0[line];" << endl;
        
        *fmodel << "data wavelength;" << endl;
        *fmodel << "observation flux;" << endl;
        
        *fmodel << "variable c = 299792458;" << endl;
        *fmodel << "variable f;" << endl;
        
        *fmodel << endl;
        
        *fmodel << "main()"  << endl;
        *fmodel << "{"  << endl;
        *fmodel << " while(import())"  << endl;
        *fmodel << " {"  << endl;
        *fmodel << "    gaussianModel();"  << endl;
        *fmodel << " }"  << endl;
        *fmodel << "}"  << endl;
        
        *fmodel << endl;
        
        *fmodel << "gaussianModel() {" << endl;
        *fmodel << "    f = level + (1.0 - a*exp(-(wavelength - wl0 - (rv + telluricRV + helioRV)*wavelength/c)*(wavelength - wl0 - rv*wavelength/c)/(2*sig*sig)));" << endl;
        *fmodel << "    export(flux - f);" << endl;
        *fmodel << "}" << endl;
        
        fmodel->close();
    }
}


void writegaussfitpars(string paramfile, operaSpectrum lines) {
    ofstream *fpars = NULL;
    
    if (!paramfile.empty()) {
        fpars = new ofstream();
        fpars->open(paramfile.c_str());
        
        *fpars << "line wl0 level a sig" << endl;
        *fpars << "double double double double double" << endl;
        for (unsigned l=0; l<lines.size(); l++) {
            *fpars << l << " " << lines.getwavelength(l) << " 0.1 0.2 0.2" << endl;
        }
        
        fpars->close();
    }
}


void writegaussfitrvpars(string rvparamfile, double telluricRV, double helioRV, double sourceRV) {
    ofstream *fpars = NULL;
    
    if (!rvparamfile.empty()) {
        fpars = new ofstream();
        fpars->open(rvparamfile.c_str());
        
        *fpars << "rv telluricRV helioRV" << endl;
        *fpars << "double double double" << endl;
        *fpars << sourceRV << " " << telluricRV << " " << helioRV << endl;
        
        fpars->close();
    }
}


void readgaussfitrvpars(string rvparamfile, double *rv, double *rverr) {
    
    operaistream fpars(rvparamfile.c_str());
    
    if (fpars.is_open()) {
        unsigned line = 0;
        string dataline;
        double buff=0;
        
        while (fpars.good()) {
            getline(fpars, dataline);
            if (strlen(dataline.c_str())) {
                stringstream ss (stringstream::in | stringstream::out);
                ss << dataline.c_str();
                
                if (line==0 || line==1) {
                    line++;
                } else {
                    ss >> *rv;
                    ss >> buff;
                    ss >> buff;

                    ss >> buff;
                    ss >> *rverr;
                    
                    line++;
                }
            }
        }
        fpars.close();
    }
}



/*
 * Generate plot
 */
void GenerateSelectedLinesPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafile, operaSpectrum sourceLines)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
    remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
    ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "\nset xlabel \"{/Symbol l} - {/Symbol l}_c (nm)\"" << endl;
    fgnu << "set ylabel \"relative flux + shift\"" << endl;
    
    unsigned nl = sourceLines.size();

    fgnu << "set xrange[-2.2:1.5]" << endl;
    for (unsigned l=0; l<nl; l++) {
        fgnu << "set label \"#" << l <<" {/Symbol l}_c=" << sourceLines.getwavelength(l) << "\"  at " << "-2.0," << l+0.9 << endl;
    }
    
    double lineAmplificationForBetterView = 5.0;
    
    if(!outputPlotEPSFileName.empty()) {
        
        fgnu << "plot \"" << datafile << "\" u ($3-$2):(($4-1)*" << lineAmplificationForBetterView << "+1+$1):(sqrt($5)) w yerr pt 7 ps 0.5" << endl;
        fgnu << "replot \"" << datafile << "\" u ($3-$2):((($4-$6)-1)*" << lineAmplificationForBetterView << "+1+$1) w l lw 2" << endl;
        
        fgnu << endl;

        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        fgnu << "replot" << endl;
        fgnu << "set output" << endl;
        fgnu << "set terminal x11" << endl;

    } else {
        
        fgnu << "plot \"" << datafile << "\" u ($3-$2):(($4-1)*" << lineAmplificationForBetterView << "+1+$1):(sqrt($5)) w yerr pt 7 ps 0.5" << endl;
        fgnu << "replot \"" << datafile << "\" u ($3-$2):((($4-$6)-1)*" << lineAmplificationForBetterView << "+1+$1) w l lw 2" << endl;

        fgnu << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set output" << endl;
        fgnu << "#set terminal x11" << endl;
    }
    
    fgnu.close();
    
    if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}


