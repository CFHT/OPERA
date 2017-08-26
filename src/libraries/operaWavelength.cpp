/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaWavelength.h"
#include "libraries/Gaussian.h"
#include "libraries/operaFit.h"	 // for operaMPFitPolynomial and operaLMFitPolynomial
#include "libraries/operaMath.h"	 // for DiffPolynomialFunction
#include "libraries/operaSpectralTools.h"		// for operaSpectralLineList

#define NPOINTPERSIGMA 4

/*!
 * operaWavelength
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the wavelength object.
 * \file operaWavelength.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaWavelength
 * \brief Encapsulation of Wavelength information.
 * \return none
 */

operaWavelength::operaWavelength() :
dmin(0.0),
dmax(0.0),
nDataPoints(0),
nAtlasLines(0),
nComparisonLines(0), 
radialVelocityPrecision(0.0),
xcorrelation(0.0)
{ }

operaWavelength::operaWavelength(unsigned Coeffs) :
dmin(0.0),
dmax(0.0),
nDataPoints(0),
nAtlasLines(0),
nComparisonLines(0), 
radialVelocityPrecision(0.0),
xcorrelation(0.0),
wavelengthPolynomial(Coeffs)
{ }

// Wavelength polynomial

Polynomial *operaWavelength::getWavelengthPolynomial(void) {
	return &wavelengthPolynomial;
}

const Polynomial *operaWavelength::getWavelengthPolynomial(void) const {
	return &wavelengthPolynomial;
}

double operaWavelength::evaluateWavelength(double distanceValue) const {
    return wavelengthPolynomial(distanceValue);
}

operaVector operaWavelength::evaluateWavelength(const operaVector& distancevalues) const {
	return Operation(wavelengthPolynomial, distancevalues);
}

double operaWavelength::convertPixelToWavelength(double DeltaDistanceInPixels) const {
	const operaVector& par = wavelengthPolynomial.getCoefficients();
    return DeltaDistanceInPixels * DiffPolynomialFunction(getcentralWavelength(), par.datapointer(), par.size());
}

void operaWavelength::applyRadialVelocityCorrection(double rvshift_InKPS) {
    operaVector coeffs = wavelengthPolynomial.getCoefficients();
    coeffs *= (1 + rvshift_InKPS/SPEED_OF_LIGHT_KMS);
    wavelengthPolynomial.setCoefficients(coeffs);
}

// Min and max distances

double operaWavelength::getDmin(void) const {
	return dmin;
}

double operaWavelength::getDmax(void) const {
	return dmax;
}

void operaWavelength::setDmin(double Dmin) {
	dmin = Dmin;
}

void operaWavelength::setDmax(double Dmax) {
	dmax = Dmax;
}

double operaWavelength::getinitialWavelength(void) const {
	return wavelengthPolynomial.Evaluate(dmin);
}

double operaWavelength::getfinalWavelength(void) const {
	return wavelengthPolynomial.Evaluate(dmax);
}

double operaWavelength::getcentralWavelength(void) const {
	return wavelengthPolynomial.Evaluate((dmin + dmax)/2);
}

// Distance-wavelength data points

unsigned operaWavelength::getnDataPoints(void) const {
    return nDataPoints;
}

double operaWavelength::getDistance(unsigned index) const {
	return distanceData[index];
}

double operaWavelength::getWavelength(unsigned index) const {
	return wavelengthData[index];
}

double operaWavelength::getWavelengthError(unsigned index) const {
	return wavelengthErrors[index];
}

unsigned operaWavelength::getMatchAtlasIndex(unsigned index) const {
	return matchAtlasindex[index];
}

unsigned operaWavelength::getMatchComparisonIndex(unsigned index) const {
	return matchComparisonindex[index];
}

void operaWavelength::createDataVectors(const operaVector& WavelengthData, const operaVector& WavelengthErrors, const operaVector& DistanceData) {
    distanceData = DistanceData;
    wavelengthData = WavelengthData;
    wavelengthErrors = WavelengthErrors;
    nDataPoints = DistanceData.size();
    matchAtlasindex.resize(nDataPoints);
    matchComparisonindex.resize(nDataPoints);
    for(unsigned i=0; i< nDataPoints; i++ ) {
        matchAtlasindex[i] = i;
        matchComparisonindex[i] = i;
    }
}

void operaWavelength::resize(unsigned NDataPoints) {
	distanceData.resize(NDataPoints);
	wavelengthData.resize(NDataPoints);
	wavelengthErrors.resize(NDataPoints);
	matchAtlasindex.resize(NDataPoints);
	matchComparisonindex.resize(NDataPoints);
    nDataPoints = NDataPoints;
}

void operaWavelength::clear() {
	distanceData.clear();
	wavelengthData.clear();
	wavelengthErrors.clear();
	matchAtlasindex.clear();
	matchComparisonindex.clear();
    nDataPoints = 0;
}

// Atlas wavelength data

unsigned operaWavelength::getnAtlasLines(void) const {
    return nAtlasLines;
}

double operaWavelength::getatlasLinesflux(unsigned index) const {
	return atlasLinesflux[index];
} 

double operaWavelength::getatlasLineswl(unsigned index) const {
	return atlasLineswl[index];
}

double operaWavelength::getatlasLineswlError(unsigned index) const {
	return atlasLineswlError[index];
}

void operaWavelength::setAtlasDataVectors(const operaSpectralLineList& atlasLines) {
    atlasLineswl = atlasLines.center;
    atlasLineswlError = atlasLines.centerError;
    atlasLinesflux = atlasLines.amplitude;
    nAtlasLines = atlasLines.size();
}

//  Comparison pixel data

unsigned operaWavelength::getnComparisonLines(void) const {
    return nComparisonLines;
}  

double operaWavelength::getcomparisonLinesflux(unsigned index) const {
	return comparisonLinesflux[index];
} 

double operaWavelength::getcomparisonLinespix(unsigned index) const {
	return comparisonLinespix[index];
} 

double operaWavelength::getcomparisonLinespixError(unsigned index) const {
 	return comparisonLinespixError[index];
} 

double operaWavelength::getcomparisonLineswl(unsigned index) const {
	return comparisonLineswl[index];
} 

void operaWavelength::setComparisonDataVectors(const operaSpectralLineList& comparisonLines){
    comparisonLinesflux = comparisonLines.amplitude;
    comparisonLinespix = comparisonLines.center;
    comparisonLinespixError = comparisonLines.centerError;
    comparisonLineswl = evaluateWavelength(comparisonLines.center);
    nComparisonLines = comparisonLines.size();
}

void operaWavelength::recalculateComparisonLineswlVector(void) {
    comparisonLineswl = evaluateWavelength(comparisonLinespix);
}

// Line matching and filtering

void operaWavelength::matchAtlaswithComparisonLines(double acceptableMismatch) {
    
    double wlc = getcentralWavelength();
    double lineSigma = wlc / spectralResolution.value;
    double acceptMismatchInwlUnits = acceptableMismatch*lineSigma;
    
    recalculateComparisonLineswlVector();
    
    clear();
    
    // Identify and select the set of lines that match both comparison and atlas.
    // The criteria for matching is the difference between centers must be < acceptableMismatch * sigma 
    unsigned nextfirstline = 0;
	
    for (unsigned i=0; i<getnComparisonLines(); i++) {
        
        unsigned bestAtlasMatchIndex = 0;
        double mindifference = acceptMismatchInwlUnits;
        
        for(unsigned l=nextfirstline;l<getnAtlasLines();l++) {
            
            double difference = fabs(comparisonLineswl[i] - atlasLineswl[l]);
            
            if(difference < mindifference) {
                if(comparisonLineswl[i] > atlasLineswl[l]) {
                    mindifference = difference;
                    bestAtlasMatchIndex = l;
                } else {
                    InsertLineMatch(l, i);
                    nextfirstline = l+1;
                    break;
                }
            } else if (comparisonLineswl[i] <= atlasLineswl[l] && difference > mindifference) {
                if(bestAtlasMatchIndex) {
                    InsertLineMatch(bestAtlasMatchIndex, i);
                    nextfirstline = bestAtlasMatchIndex+1;
                }
                break;
            }
            if(l==nAtlasLines-1) {
                if(bestAtlasMatchIndex) {
                    InsertLineMatch(bestAtlasMatchIndex, i);
                    nextfirstline = bestAtlasMatchIndex+1;
                }
            }
			if (nDataPoints >= getnAtlasLines()) {
				throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (nDataPoints >= getnComparisonLines()) {
				throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
			}
        }
    }
}

void operaWavelength::InsertLineMatch(unsigned atlasindex, unsigned compareindex) {
	distanceData.insert(comparisonLinespix[compareindex]);
	wavelengthData.insert(atlasLineswl[atlasindex]);
	wavelengthErrors.insert(atlasLineswlError[atlasindex]);
	matchAtlasindex.insert(atlasindex);
	matchComparisonindex.insert(compareindex);
	nDataPoints++;
}

double operaWavelength::getPerCentageOfComparisonMatch(void) const {
    return 100*(double)getnDataPoints()/(double)nComparisonLines;
}

double operaWavelength::getPerCentageOfAtlasMatch(void) const {
    return 100*(double)getnDataPoints()/(double)nAtlasLines;
}

void operaWavelength::filterDataPointsBySigmaClip(double nsig) {
    double sigmaclip = (nsig)*calculateWavelengthRMSPrecision();
    unsigned np = 0;
    for(unsigned i=0;i<getnDataPoints();i++) {
        if(fabs(wavelengthData[i] - evaluateWavelength(distanceData[i])) < sigmaclip) {
            distanceData[np] = distanceData[i];
            wavelengthData[np] = wavelengthData[i];
            wavelengthErrors[np] = wavelengthErrors[i];
            matchAtlasindex[np] = matchAtlasindex[i];
            matchComparisonindex[np] = matchComparisonindex[i];
            np++;
        }
    }
    resize(np);
}

void operaWavelength::filterDataPointsByErrorClip(double nsig) {
    double sigmaclip = (nsig)*calculateWavelengthRMSPrecision();
    unsigned np = 0;
    for(unsigned i=0;i<getnDataPoints();i++) {
        if(wavelengthErrors[i] < sigmaclip) {
            distanceData[np] = distanceData[i];
            wavelengthData[np] = wavelengthData[i];
            wavelengthErrors[np] = wavelengthErrors[i];
            matchAtlasindex[np] = matchAtlasindex[i];
            matchComparisonindex[np] = matchComparisonindex[i];
            np++;
        }
    }
    resize(np);
}

// Resolution and precision

doubleValue_t operaWavelength::getResolutionElementInPixels(void) const {
	return resolutionElementInPixels;
}
    
void operaWavelength::setResolutionElementInPixels(doubleValue_t ResolutionElementInPixels) {
	resolutionElementInPixels = ResolutionElementInPixels;
}

doubleValue_t operaWavelength::getSpectralResolution(void) const {
    return spectralResolution;
}

void operaWavelength::setSpectralResolution(doubleValue_t Resolution) {
    spectralResolution = Resolution;
}

double operaWavelength::getRadialVelocityPrecision(void) const {
    return radialVelocityPrecision;
}

void operaWavelength::setRadialVelocityPrecision(double radialvelocityprecision) {
    radialVelocityPrecision = radialvelocityprecision;
}

void operaWavelength::calculateSpectralResolution() {
    double deltawl = convertPixelToWavelength(resolutionElementInPixels.value);
    double errorwl = convertPixelToWavelength(resolutionElementInPixels.error);
    spectralResolution.value = getcentralWavelength() / deltawl;
    spectralResolution.error = errorwl * getcentralWavelength() / (deltawl*deltawl);
}

void operaWavelength::calculateRadialVelocityPrecision(void) {
    double speedoflight = 299792458; // in m/s
    const operaVector& tempwl = evaluateWavelength(distanceData);
    operaVector scaledResiduals = (wavelengthData - tempwl)/tempwl;
    radialVelocityPrecision = speedoflight * RMS(scaledResiduals);
}

double operaWavelength::calculateWavelengthRMSPrecision(void) {
    return RMS(wavelengthData - evaluateWavelength(distanceData));
}

double operaWavelength::calculateWavelengthMedianPrecision(void) {
    operaVector wlresiduals = Abs(wavelengthData - evaluateWavelength(distanceData));
    double medianResidual = MedianQuick(wlresiduals);
    return MedianStdDev(wlresiduals, medianResidual);
}

// Cross-correlation

double operaWavelength::getxcorrelation(void) const {
    return xcorrelation;
}

void operaWavelength::setxcorrelation(double Xcorrelation) {
    xcorrelation = Xcorrelation;
}

void operaWavelength::calculateXCorrelation(void) {
    operaVector atlasSimulSpectrum = createAtlasSimulatedSpectrum(NPOINTPERSIGMA);
    operaVector comparisonSimulSpectrum = createComparisonSimulatedSpectrum(NPOINTPERSIGMA);
    double crosscorrelation = operaCrossCorrelation(atlasSimulSpectrum, comparisonSimulSpectrum);
    setxcorrelation(crosscorrelation);
}

// Simulated spectra

operaVector operaWavelength::createAtlasSimulatedSpectrum(unsigned nstepspersigma) const {
    return createSimulatedSpectrum(atlasLineswl, atlasLinesflux, nstepspersigma);
}

operaVector operaWavelength::createComparisonSimulatedSpectrum(unsigned nstepspersigma) {
    recalculateComparisonLineswlVector();
    return createSimulatedSpectrum(comparisonLineswl, comparisonLinesflux, nstepspersigma);
}

operaVector operaWavelength::createAtlasSimulatedSpectrumWithConstantFlux(unsigned nstepspersigma) const {
    operaVector lineAmplitudes(nAtlasLines);
    lineAmplitudes = 1.0;
    return createSimulatedSpectrum(atlasLineswl, lineAmplitudes, nstepspersigma);
}

operaVector operaWavelength::createComparisonSimulatedSpectrumWithConstantFlux(unsigned nstepspersigma) {
    operaVector lineAmplitudes(nComparisonLines);
    lineAmplitudes = 1.0;
    recalculateComparisonLineswlVector();
    return createSimulatedSpectrum(comparisonLineswl, lineAmplitudes, nstepspersigma);
}

operaVector operaWavelength::createSimulatedSpectrum(const operaVector& wl, const operaVector& flux, unsigned nstepspersigma) const {
    unsigned nLines = wl.size();
    if (nLines >= MAXPOINTSINSIMULATEDSPECTRUM) {
        throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double wlc = getcentralWavelength();
    double minwl = getinitialWavelength();
    double maxwl = getfinalWavelength();
    
    double lineSigma = wlc / spectralResolution.value;
    
    operaVector sigmaVector(nLines);
    sigmaVector = lineSigma;
    
    Gaussian spectrumModel(nLines, flux.datapointer(), sigmaVector.datapointer(), wl.datapointer());
    
    double wlstep = lineSigma / nstepspersigma;
    unsigned npoints = (unsigned)ceil(fabs(maxwl - minwl)/wlstep);
    double wlcurrent = minwl;
    
    operaVector outputSpectrum;
    for(unsigned i=0;i<npoints;i++) {
        outputSpectrum.insert(spectrumModel.EvaluateGaussian(wlcurrent));
        wlcurrent += wlstep;
    }
    return outputSpectrum;
}

// Wavelength solution calculation

void operaWavelength::CalculateWavelengthSolution(unsigned maxcoeffs, bool witherrors) {
    Polynomial bestfitpoly;
    bestfitpoly.setChisqr(BIG);
    
    if (maxcoeffs > wavelengthPolynomial.getOrderOfPolynomial()) {
		maxcoeffs = wavelengthPolynomial.getOrderOfPolynomial();
	}
	
	for (unsigned currentfit=1; currentfit<=maxcoeffs; currentfit++) {
		operaVector coeffs(wavelengthPolynomial.getVector(), currentfit);
		operaVector errors(wavelengthPolynomial.getErrorVector(), currentfit);
		
		Polynomial polyfit;
		if (witherrors) {
			polyfit = PolynomialFit(distanceData, wavelengthData, wavelengthErrors, coeffs, errors);
		} else {
			polyfit = PolynomialFit(distanceData, wavelengthData, coeffs);
		}
		
		if (polyfit.getChisqr() < bestfitpoly.getChisqr()) {
			bestfitpoly = polyfit;
		}
	}
	
	wavelengthPolynomial = bestfitpoly;
}

void operaWavelength::RefineWavelengthSolution(unsigned ncoeffs, bool witherrors) {
    unsigned npar = wavelengthPolynomial.getOrderOfPolynomial();
    operaVector coeffs = wavelengthPolynomial.getCoefficients();
	operaVector errors = wavelengthPolynomial.getErrors();
	coeffs.resize(ncoeffs);
	errors.resize(ncoeffs);
    if(ncoeffs > npar) {
        for(unsigned i=npar; i<ncoeffs; i++) {
			coeffs[i] = 1.0;
			errors[i] = 0.0;
        }
    }
    
	Polynomial newpoly;
    if (witherrors) {
        newpoly = PolynomialFit(distanceData, wavelengthData, wavelengthErrors, coeffs, errors);
    } else {
        newpoly = PolynomialFit(distanceData, wavelengthData, coeffs);
    }
    if (newpoly.getChisqr() < wavelengthPolynomial.getChisqr()) {
		wavelengthPolynomial = newpoly;
    }
}

void operaWavelength::refineWavelengthSolutionOfSecondOrderByXCorrelation(unsigned nPointsPerParameter, double parameterRangetoSearch) {
    // EM May 25 2015 -- it only works to search solutions on a 2nd order polynomial (parabola).
    // I changed this function because the previous version was wrong and it wouldn't work. However it is not used by ESPaDOnS.
    
    int npar = 3;  // Force 2nd order polynomial
    wavelengthPolynomial.resize(npar);
    
    // Note that the input parameters determine the step and range for which the coefficients will be searched
    operaVector par = wavelengthPolynomial.getCoefficients();
    par[npar-1] = 0;
    operaVector range = Abs(par) * parameterRangetoSearch/100.0;
    range[npar-1] = 2.0*1e-5;
    operaVector delta = range / double(nPointsPerParameter);
    operaVector initial = par - (range/2.0);
    operaVector maxpar = par;
    
    // Attempt to find the coefficients that gives the highest correlation between the raw and atlas simulated spectra.
    double maxcorrelation = -1.0;
    
    par[2] = initial[2];
    for (unsigned k=0;k<nPointsPerParameter; k++) {
        par[1] = initial[1];
        for (unsigned j=0;j<nPointsPerParameter; j++) {
            par[0] = initial[0];
            for (unsigned i=0;i<nPointsPerParameter; i++) {
                recalculateComparisonLineswlVector();
                
                operaVector atlasSimulSpectrum = createAtlasSimulatedSpectrumWithConstantFlux(NPOINTPERSIGMA);
                operaVector comparisonSimulSpectrum = createComparisonSimulatedSpectrumWithConstantFlux(NPOINTPERSIGMA);
                
                double crosscorrelation = operaCrossCorrelation(atlasSimulSpectrum, comparisonSimulSpectrum);
                if(crosscorrelation > maxcorrelation) {
                    maxcorrelation = crosscorrelation;
                    maxpar = par;
                }
                par[0] += delta[0];
            }
            par[1] += delta[1];
        }
        par[2] += delta[2];
    }
    setxcorrelation(maxcorrelation);
    wavelengthPolynomial.setCoefficients(maxpar);
    
    recalculateComparisonLineswlVector();
}

void operaWavelength::refineWavelengthSolutionByFindingMaxMatching(unsigned NpointsPerPar, double ParRangeSizeInPerCent, double acceptableMismatch) {
    double coeff = wavelengthPolynomial.getCoefficient(0);
    double searchRange = fabs(coeff * ParRangeSizeInPerCent/100.0);
    double delta = searchRange / NpointsPerPar;
    
    double maxcoeff = coeff;
    double maxpercentage = 0;
    coeff -= searchRange/2.0;
    
    for(unsigned i = 0; i < NpointsPerPar; i++) {
		wavelengthPolynomial.setCoefficient(0, coeff);
        matchAtlaswithComparisonLines(acceptableMismatch);
        double MatchPercentage = (getPerCentageOfComparisonMatch() + getPerCentageOfAtlasMatch())/2;
        if(MatchPercentage > maxpercentage) {
            maxpercentage = MatchPercentage;
            maxcoeff = coeff;
        }
        coeff += delta;
    }
    wavelengthPolynomial.setCoefficient(0, maxcoeff);
}
