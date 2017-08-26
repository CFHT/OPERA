#ifndef OPERAWAVELENGTHCALIBRATION_H
#define OPERAWAVELENGTHCALIBRATION_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operWavelengthCalibration
 Version: 1.0
 Description: Perform wavelength calibration 
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

/*! \brief Perform wavelength calibration. */
/*! \file operWavelengthCalibration.h */
/*! \ingroup core */

#include <map>

class operaSpectralElements;
class operaSpectralOrder;

struct DetectionParameters {
	double RawLineWidth;
	double RawLineWidthErr;
	double LocalMaxFilterWidth;
	double MinPeakDepth;
	double DetectionThreshold;
};

class WavelengthSolutions {
private:
	std::map <int, Polynomial> polynomials;
	void readLineSet(string inputLineSetFilename, int order, operaVector& wavelength, operaVector& distance) const;
public:
	void SetFromSpectralOrders(const operaSpectralOrderVector& spectralOrders, int minorder, int maxorder);
	void CalculateFromLineSet(string inputLineSetFilename, int minorder, int maxorder, unsigned maxorderofpolynomial);
	bool HasSolutionForOrder(int order) const;
	const Polynomial& GetSolutionForOrder(int order) const;
};

class WavelengthCalibration {
public:
	WavelengthCalibration(operaSpectralOrder* spectralOrder, const operaSpectrum& atlasSpectrum, const operaSpectralLineList& atlasLines);

	void SetFromInitialSolution(const WavelengthSolutions& initialSolutions, int order);

	void CalculateWavelengthSolution(double ParRangeSizeInPerCent, DetectionParameters detectionParams, double nsigclip, bool parseSolution, unsigned NpointsPerPar, double initialAcceptableMismatch, double dampingFactor, unsigned minNumberOfLines, unsigned maxNIter);

	void InitializeDistanceLimitsFromGeometry();
	void FindAndSetComparisonAndAtlasLines(DetectionParameters detectionParams, double nsigclip);
	void FitSolutionPolynomialUsingMatchingLines(double acceptableMismatch);
	void RefineSolutionPolynomialFit(double acceptableMismatch, double dampingFactor, unsigned minNumberOfLines, double nsigclip, unsigned maxNIter);
	void FinishedShrinkingPolynomialFit();

	bool HasUncalLines() const;
	bool HasUncalSpectrum() const;
	operaSpectralLineList GetUncalLines();
	operaSpectralLineList GetLinesFromUncalSpectrum(DetectionParameters detection, double nsigclip);
	operaSpectralLineList GetAtlasLinesInRange(operaWavelengthRange wlrange, double rawLineWidth);
	operaSpectralLineList GetLinesFromAtlasSpectrum(operaWavelengthRange wlrange, DetectionParameters detection);
	void SetComparisonLines(const operaSpectralLineList& comparisonLines);
	void SetAtlasLines(const operaSpectralLineList& atlasLines);

	operaSpectralLineList DetectAtlasLines(operaWavelengthRange wlrange, DetectionParameters detection);

	void PrintSolution();
	void PrintPrecisionAndResolution();

	void WritePlotAtlasData(const operaSpectrum& atlasRegion);
	void WritePlotComparisonData(double rawlinewidth);
	void WritePlotLinesData();
	void WritePlotOrdersData();

	static ofstream fatlasdata;
	static ofstream fcompdata;
	static ofstream flinesdata;
	static ofstream fordersdata;
	static bool generate3DPlot;
	static bool skipPlots;

private:
	operaSpectralOrder* spectralOrder;
	const operaSpectrum& atlasSpectrumFull;
	const operaSpectralLineList& atlasLinesFull;
	bool hasComparisonLines;
	bool hasAtlasLines;
	double comparisonLineWidth;
	double comparisonLineWidthErr;
	unsigned order;
	operaWavelength* wavelength;
	double wl_central;
};

operaSpectralLineList readThoriumArgonAtlas(string atlas_lines);

operaSpectrum readAtlasSpectrum(string atlas_spectrum);

int DetermineOrderShift(operaSpectralOrderVector& spectralOrders, int referenceMinOrder, int referenceMaxOrder, const WavelengthSolutions& initialSolutions, int nOrdersToSearchAround, const operaSpectrum& atlasSpectrum, const operaSpectralLineList& atlasLines, double ParRangeSizeInPerCent, DetectionParameters detectionParams, double nsigclip, unsigned NpointsPerPar, double initialAcceptableMismatch, double dampingFactor, unsigned minNumberOfLines, unsigned maxNIter);

void GenerateWavelengthOrdersPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display);

void GenerateWavelengthSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, string linesdatafilename, bool subtractCentralWavelength, bool display);

void GenerateWavelength3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, bool subtractCentralWavelength, bool display);

void GenerateWavelengthSolutionPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned npolynomials, Polynomial *polynomials[], bool display);

#endif
