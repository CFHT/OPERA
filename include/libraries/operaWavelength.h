/*******************************************************************
 ****               		OPERA PIPELINE v1.0                     ****
 ********************************************************************
 Library name: operaWavelength
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 
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

#ifndef OPERAWAVELENGTH_H
#define OPERAWAVELENGTH_H

#include "libraries/operaLibCommon.h"			// for doubleValue_t
#include "libraries/Polynomial.h"	
#include "libraries/operaWavelength.h"

#define MAXORDEROFWAVELENGTHPOLYNOMIAL 10

/*! 
 * \sa class operaWavelength
 * \brief Encapsulation of Wavelength information.
 * \return none
 * \file operaWavelength.h
 * \ingroup libraries
 */

class operaSpectralLineList;
class operaSpectralOrder;

class operaWavelength {
	
private:
	Polynomial wavelengthPolynomial; // Polynomial lambda(d) representing the relation between distance and wavelength
    
	// Minimum and maximum distances from geometry information
	double dmin;
	double dmax;		
    
    // Distance and wavelength for each point to model the wavelength polynomial
    unsigned nDataPoints;
    operaVector distanceData;
    operaVector wavelengthData;
    operaVector wavelengthErrors;
    Vector<unsigned> matchAtlasindex; // The index of the atlas line that the wavelength came from
    Vector<unsigned> matchComparisonindex; // The index of the comparison line that the distance came from

    // Atlas lines of used for line matching, by wavelength
    unsigned nAtlasLines;
    operaVector atlasLinesflux;
    operaVector atlasLineswl;
    operaVector atlasLineswlError;

    // Comparison lines used for line matching, by distance
    unsigned nComparisonLines;
    operaVector comparisonLinesflux;
    operaVector comparisonLinespix;
    operaVector comparisonLinespixError;
    operaVector comparisonLineswl; // Comparison line wavelengths are calculated using the current polynomial

    // Precision/resolution values
    doubleValue_t resolutionElementInPixels;
    doubleValue_t spectralResolution;
    double radialVelocityPrecision;
    
    double xcorrelation; // Cross-correlation between simulated spectra made from the comparison and atlas lines
    
    operaVector createSimulatedSpectrum(const operaVector& wl, const operaVector& flux, unsigned nstepspersigma) const;
    void InsertLineMatch(unsigned atlasindex, unsigned compareindex); // Helper function to insert data points from matching atlas and comparison lines
	
public:
	operaWavelength();
	operaWavelength(unsigned Coeffs);
	
	Polynomial *getWavelengthPolynomial(void);
	const Polynomial *getWavelengthPolynomial(void) const;
	double evaluateWavelength(double distanceValue) const;
    operaVector evaluateWavelength(const operaVector& distancevalues) const;
    double convertPixelToWavelength(double DeltaDistanceInPixels) const;
    void applyRadialVelocityCorrection(double rvshift_InKPS);
    
	double getDmin(void) const;
	double getDmax(void) const;
	void setDmin(double Dmin);
	void setDmax(double Dmax);
	double getinitialWavelength(void) const;
	double getfinalWavelength(void) const;
	double getcentralWavelength(void) const;
	
	unsigned getnDataPoints(void) const;
	double getDistance(unsigned index) const;
	double getWavelength(unsigned index) const;
	double getWavelengthError(unsigned index) const;
	unsigned getMatchAtlasIndex(unsigned index) const;
    unsigned getMatchComparisonIndex(unsigned index) const;
    void createDataVectors(const operaVector& WavelengthData, const operaVector& WavelengthErrors, const operaVector& DistanceData);
    void resize(unsigned NDataPoints);
	void clear();
	
	unsigned getnAtlasLines(void) const;   
    double getatlasLinesflux(unsigned index) const;
    double getatlasLineswl(unsigned index) const;
    double getatlasLineswlError(unsigned index) const;    
    void setAtlasDataVectors(const operaSpectralLineList& atlasLines);
    
    unsigned getnComparisonLines(void) const;
    double getcomparisonLinesflux(unsigned index) const;
    double getcomparisonLinespix(unsigned index) const;
    double getcomparisonLinespixError(unsigned index) const;
    double getcomparisonLineswl(unsigned index) const;
    void setComparisonDataVectors(const operaSpectralLineList& comparisonLines);
    void recalculateComparisonLineswlVector(void);
    
    void matchAtlaswithComparisonLines(double acceptableMismatch);    
    double getPerCentageOfComparisonMatch(void) const;
    double getPerCentageOfAtlasMatch(void) const;
    void filterDataPointsBySigmaClip(double nsig);
    void filterDataPointsByErrorClip(double nsig);
    
	doubleValue_t getResolutionElementInPixels(void) const;
    void setResolutionElementInPixels(doubleValue_t ResolutionElementInPixels);
    doubleValue_t getSpectralResolution(void) const;
	void setSpectralResolution(doubleValue_t Resolution);
    double getRadialVelocityPrecision(void) const;
    void setRadialVelocityPrecision(double radialvelocityprecision);    
    void calculateSpectralResolution();
    void calculateRadialVelocityPrecision(void);
    double calculateWavelengthRMSPrecision(void);
    double calculateWavelengthMedianPrecision(void);   
    
    double getxcorrelation(void) const;
	void setxcorrelation(double Xcorrelation); 
    void calculateXCorrelation(void);
    
	operaVector createAtlasSimulatedSpectrumWithConstantFlux(unsigned nstepspersigma) const;
    operaVector createComparisonSimulatedSpectrumWithConstantFlux(unsigned nstepspersigma);
    operaVector createAtlasSimulatedSpectrum(unsigned nstepspersigma) const;
    operaVector createComparisonSimulatedSpectrum(unsigned nstepspersigma);
    
	void CalculateWavelengthSolution(unsigned maxcoeffs, bool witherrors);
    void RefineWavelengthSolution(unsigned ncoeffs, bool witherrors);    
	void refineWavelengthSolutionOfSecondOrderByXCorrelation(unsigned nPointsPerParameter, double parameterRangetoSearch);
    void refineWavelengthSolutionByFindingMaxMatching(unsigned NpointsPerPar, double ParRangeSizeInPerCent, double acceptableMismatch);
};

#endif
