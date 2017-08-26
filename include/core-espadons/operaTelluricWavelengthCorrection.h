#ifndef OPERATELLURICWAVELENGTHCORRECTION_H
#define OPERATELLURICWAVELENGTHCORRECTION_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaTelluricWavelengthCorrection
 Version: 1.0
 Description: Apply wavelength correction based on telluric lines
 to start up with an OPERA module.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: SEP/2012
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

/*! \brief Apply wavelength correction based on telluric lines. */
/*! \file operaTelluricWavelengthCorrection.h */
/*! \ingroup core */

#define MAXNUMBEROFPOINTSINTELLURICSPECTRUM 2400000
#define MAXNUMBEROFLINESINTELLURICDATABASE 100000
#define MAXLENGTHOFLINEINTELLURICDATABASE 160
#define TYPICAL_PRESSURE_AT_MAUNAKEA 61000  // Pa
#define TYPICAL_TEMPERATURE_AT_MAUNAKEA 273 // K
#define k_BOLTZMANN_CONSTANT 1.3806503e23 // m2 kg s-2 K-1
#define TYPICAL_ATMOSPHERE_PATH_LENGTH  843500 // cm

enum ProfileMethod { GAUSSIAN, LORENTZ, VOIGT };

operaSpectrum readTelluricLinesHITRAN(string telluric_database_file);

operaSpectrum readTelluricLinesRaw(string telluric_database_file);

operaVector generateSyntheticTelluricSpectrumUsingLineProfile(const operaSpectrum& telluricLines, const operaVector& wavelengthVector, double resolution, ProfileMethod profile);

void getWavelengthSubrange(const operaVector& wavelength, double wl0, double wlf, unsigned& startindex, unsigned& endindex);

bool calculateRVShiftByXCorr(const operaSpectrum& telluricLines, const operaSpectrum& objectSpectrum, double radialVelocityRange, double radialVelocityStep, double threshold, double& maxRV, double& sigRV, double& maxcorr, ofstream& fxcorrdata, ofstream& fxcorrfitdata, double spectralResolution, bool useFitToFindMaximum, double& chisqr);

void GenerateTelluricXCorrelationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string cleanDataFileName);

void GenerateTelluricRVCorrPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string histDataFileName, float rvshift, float rvshifterror);

void GenerateTelluricLineMatchPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename, string matchdatafilename);

void GenerateTelluricSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename);

void matchTelluricLines(const operaSpectrum& atlasLines, const operaSpectralLineList& objectLines, operaVector& atlasMatchedWavelengths, operaSpectrum& objectMatchedLines, operaVector& radialVelocities, double spectralResolution, double radialVelocityRange, double duplicateLineDist);

void generateHistogramData(const operaVector& telluricMatchedWavelengths, const operaVector& radialVelocities, double radialVelocityRange, double radialVelocityStep, operaVector& rvVector, operaVector& probDensity, operaVector& wavelengthVector);

#endif
