#ifndef OPERASPECTRALTOOLS_H
#define OPERASPECTRALTOOLS_H

/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
Library name: operaSpectralTools
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include <stdarg.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include "libraries/operaSpectralElements.h"
#include "libraries/operaSpectralLines.h"
#include "libraries/ladfit.h"						// for ladfit

#define MAXNUMBEROFWLRANGES 1000
#define MINNUMBEROFPOINTSINSIDEBIN 5


/*! 
 * \brief general library routines.
 * \file operaSpectralTools.h
 * \ingroup libraries
 */

using namespace std;

// Class to store multiple flux vectors of the same length and process them simultaneously.
class operaFluxVectors {
private:
	std::vector<operaFluxVector> vectors;
	void lengthcheck(const operaFluxVector& samelength) const { if (vectors.front().getlength() != samelength.getlength()) throw operaException("operaFluxVectors: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__); }
	void countcheck(const operaVector& samecount) const { if(vectors.size() != samecount.size()) throw operaException("operaFluxVectors: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__); }
	void countcheck(unsigned count) const { if(vectors.size() != count) throw operaException("operaFluxVectors: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__); }
	void rangecheck(unsigned index) const { if(index >= vectors.size()) throw operaException("operaFluxVectors: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__); }
public:
	operaFluxVectors(unsigned count) : vectors(count) { }
	operaFluxVectors(unsigned count, unsigned length) : vectors(count, operaFluxVector(length)) { }
	unsigned vectorcount() const { return vectors.size(); }
	unsigned length() const { return vectors.front().getlength(); }
    void clear() { for(unsigned v=0; v<vectors.size(); v++) vectors[v].clear(); }
	void resize(unsigned newlength) { for(unsigned v=0; v<vectors.size(); v++) vectors[v].resize(newlength); }
	void trim(operaIndexRange range) { for(unsigned v=0; v<vectors.size(); v++) vectors[v].trim(range); }
	void insert(const operaVector& newfluxes, const operaVector& newvariances) { countcheck(newfluxes); countcheck(newvariances); for(unsigned v=0; v<vectors.size(); v++) vectors[v].insert(newfluxes[v], newvariances[v]); }
	void insert(double newflux, double newvariance) { countcheck(1); vectors[0].insert(newflux, newvariance); }
	void reverse() { for(unsigned v=0; v<vectors.size(); v++) vectors[v].reverse(); }
	void reorder(const operaIndexMap& indexmap) { for(unsigned v=0; v<vectors.size(); v++) vectors[v].reorder(indexmap); }
    const operaFluxVector& getvector(unsigned v) const { rangecheck(v); return vectors[v]; }
	operaFluxVector& getvector(unsigned v) { rangecheck(v); return vectors[v]; } //this should hopefully go away eventually, but some things require a non-const pointer currently
	operaVector getfluxes(unsigned index) const { operaVector fluxes; for(unsigned v=0; v<vectors.size(); v++) fluxes.insert(vectors[v].getflux(index)); return fluxes; }
    operaVector getvariances(unsigned index) const { operaVector variances; for(unsigned v=0; v<vectors.size(); v++) variances.insert(vectors[v].getvariance(index)); return variances; }
    operaVector geterrors(unsigned index) const { operaVector errors; for(unsigned v=0; v<vectors.size(); v++) errors.insert(vectors[v].geterror(index)); return errors; }
    void setvector(const operaFluxVector& newvector, unsigned v) { rangecheck(v); lengthcheck(newvector); vectors[v] = newvector; }
	void setfluxes(const operaVector& newfluxes, unsigned index) { countcheck(newfluxes); for(unsigned v=0; v<vectors.size(); v++) vectors[v].setflux(newfluxes[v], index); }
    void setvariances(const operaVector& newvariances, unsigned index) { countcheck(newvariances); for(unsigned v=0; v<vectors.size(); v++) vectors[v].setvariance(newvariances[v], index); }
};

// Class to store and process a wavelength vector and one or more corresponding flux vectors of the same length, with or without variances.
class operaSpectrum {
protected:
	operaVector wavelength;
	operaFluxVectors fluxvectors;
public:
	operaSpectrum() : fluxvectors(1) { }
	operaSpectrum(unsigned presize) : wavelength(presize), fluxvectors(1, presize) { }
	operaSpectrum(unsigned vectors, unsigned presize) : wavelength(presize), fluxvectors(vectors, presize) { }
	operaSpectrum(const operaVector& wavelength) : wavelength(wavelength), fluxvectors(1, wavelength.size()) { }
	operaSpectrum(unsigned vectors, const operaVector& wavelength) : wavelength(wavelength), fluxvectors(vectors, wavelength.size()) { }
	unsigned size() const { return wavelength.size(); }
	unsigned fluxcount() const { return fluxvectors.vectorcount(); }
	bool empty() const { return wavelength.empty(); }
	void clear() { wavelength.clear(); fluxvectors.clear(); }
	void insert(double wl, double flux, double variance = 0) { wavelength.insert(wl); fluxvectors.insert(flux, variance); }
	void insert(double wl, const operaVector& fluxes, const operaVector& variances) { wavelength.insert(wl); fluxvectors.insert(fluxes, variances); }
	void insert(double wl, const operaVector& fluxes) { insert(wl, fluxes, operaVector(fluxes.size())); }
	void reverse() { wavelength.reverse(); fluxvectors.reverse(); }
	void resize(unsigned newsize) { wavelength.resize(newsize); fluxvectors.resize(newsize); }
	void sort() { operaIndexMap indexmap = wavelength.indexsort(); wavelength.reorder(indexmap); fluxvectors.reorder(indexmap); }
	void insertfrom(const operaSpectrum& input, unsigned i) { insert(input.getwavelength(i), input.getfluxes(i), input.getvariances(i)); }
	
	double getwavelength(unsigned i) const { return wavelength[i]; }
	double getflux(unsigned i) const { return fluxvectors.getvector(0).getflux(i); }
	double getvariance(unsigned i) const { return fluxvectors.getvector(0).getvariance(i); }
	operaVector getfluxes(unsigned i) const { return fluxvectors.getfluxes(i); }
	operaVector getvariances(unsigned i) const { return fluxvectors.getvariances(i); }
	
	double firstwl() const { return wavelength.first(); }
	double midwl() const { return wavelength[wavelength.size()/2]; }
	double lastwl() const { return wavelength.last(); }
	double medianwl () const { return Median(wavelength); }
	double meanwl () const { return Mean(wavelength); }
	double medianflux (unsigned v=0) const { return Median(fluxvector(v)); }
	double meanflux (unsigned v=0) const { return Mean(fluxvector(v)); }
	operaVector medianfluxes() const { operaVector temp; for(unsigned v = 0; v < fluxcount(); v++) temp.insert(medianflux(v)); return temp; }
	operaVector meanfluxes() const { operaVector temp; for(unsigned v = 0; v < fluxcount(); v++) temp.insert(meanflux(v)); return temp; }
	
	const operaVector& wavelengthvector() const { return wavelength; }
    const operaFluxVector& operafluxvector(unsigned v=0) const { return fluxvectors.getvector(v); }
    const operaVector& fluxvector(unsigned v=0) const { return fluxvectors.getvector(v).getflux(); }
    const operaVector& variancevector(unsigned v=0) const { return fluxvectors.getvector(v).getvariance(); }
    void setoperafluxvector(const operaFluxVector& newvector, unsigned v=0) { fluxvectors.setvector(newvector, v); }
    
    double* wavelength_ptr() { return wavelength.datapointer(); }
	double* flux_ptr(unsigned v=0) { return fluxvectors.getvector(v).getfluxpointer(); }
	double* variance_ptr(unsigned v=0) { return fluxvectors.getvector(v).getvariancepointer(); }
	const double* wavelength_ptr() const { return wavelength.datapointer(); }
	const double* flux_ptr(unsigned v=0) const { return fluxvectors.getvector(v).getfluxpointer(); }
	const double* variance_ptr(unsigned v=0) const { return fluxvectors.getvector(v).getvariancepointer(); }
};

// Class to hold the attributes of a list of spectral lines
class operaSpectralLineList {
private:
	void inserter(const operaVector& src, unsigned srcindex, operaVector& dest) { if(src.size() >= srcindex) dest.insert(src[srcindex]); }
public: // Everything public for simplicity's sake
	operaVector center; // wavelength
	operaVector centerError; // wavelength error
	operaVector sigma; // line width
	operaVector amplitude; // flux
	double medianWidth;
	double medianWidthError;
	unsigned size() const { return center.size(); }
	void sort() { operaIndexMap indexmap = center.indexsort(); center.reorder(indexmap); centerError.reorder(indexmap); sigma.reorder(indexmap); amplitude.reorder(indexmap); }
	double getwavelength(unsigned i) const { return center[i]; }
	double getflux(unsigned i) const { return amplitude[i]; }
	double getsigma(unsigned i) const { return sigma[i]; }
	void insertfrom(const operaSpectralLineList& src, unsigned i) { inserter(src.center, i, center); inserter(src.centerError, i, centerError); inserter(src.sigma, i, sigma); inserter(src.amplitude, i, amplitude); }
};

class operaWavelengthRange {
private:
	double wl0;
	double wlf;
public:
	operaWavelengthRange(double start, double end) : wl0(start), wlf(end) { }
    bool contains(double wavelength) { return wl0 <= wavelength && wavelength <= wlf; }
    double getwl0() const { return wl0; }
    double getwlf() const { return wlf; }
};

class operaWavelengthRanges {
private:
	vector <operaWavelengthRange> ranges;
public:
    void addWavelengthRange(double start, double end) { ranges.push_back(operaWavelengthRange(start, end)); }
    double wl0(unsigned i) const { return ranges[i].getwl0(); }
    double wlf(unsigned i) const { return ranges[i].getwlf(); }
    bool contains(double wavelength) { for(unsigned i = 0; i < ranges.size(); i++) if(ranges[i].contains(wavelength)) return true; return false; }
    bool contains(double wavelength, unsigned i) { return ranges[i].contains(wavelength); }
    operaWavelengthRange getrange(unsigned i) { return ranges[i];}
    unsigned int size() const { return ranges.size(); }
};

double operaCrossCorrelation(operaVector a, operaVector b);

/* 
 * calculateXCorrWithGaussian(const operaVector& wavelength, const operaFluxVector& flux, double sigma)
 * \brief This function calculates the cross-correlation between an input spectrum and gaussian function
 * \param const operaVector& wavelength
 * \param const operaFluxVector& flux
 * \param double sigma
 * \return operaVector
 */
operaVector calculateXCorrWithGaussian(const operaVector& wavelength, const operaVector& flux, double sigma);

/*
 * convolveSpectrumWithGaussian(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double sigma)
 * \brief This function calculates the convolution between an input spectrum and a gaussian function
 * \param unsigned np
 * \param double *wavelength
 * \param double *flux
 * \param return double *convolvedSpectrum
 * \param double sigma
 * \return void
 */
operaVector convolveSpectrumWithGaussian(const operaVector& wavelength, const operaVector& flux, double sigma);


/*
 * convolveSpectrumWithGaussianByResolution(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double spectralResolution)
 * \brief This function calculates the convolution between an input spectrum and a gaussian function using the spectral Resolution to calculate line width
 * \param unsigned np
 * \param double *wavelength
 * \param double *flux
 * \param return double *convolvedSpectrum
 * \param double spectralResolution
 * \return void
 */
operaVector convolveSpectrumWithGaussianByResolution(const operaVector& wavelength, const operaVector& flux, double spectralResolution);

/*
 * normalizeSpectrum(unsigned nLines, double *lineflux)
 * \brief This function normalize an input vector of fluxes by the maximum value in the array
 * \param unsigned nLines
 * \param double *lineflux
 * \return void
 */
void normalizeSpectrum(unsigned nLines, double *lineflux);

/*
 * normalizeSpectrum(unsigned nLines, double *lineflux)
 * \brief This function normalize an input vector of fluxes and respective variances by the maximum value in the array
 * \param unsigned nLines
 * \param double *lineflux
 * \param double *linevariance 
 * \return void
 */
void normalizeSpectrum(unsigned nLines, double *lineflux, double *linevariance);

/*
 * double convertVacuumToAirWavelength(double vac_wl)
 * \brief This function converts from vacuum to air wavelength using the IAU standard for conversion
 * \brief  from air to vacuum wavelengths as given in Morton (1991, ApJS, 77, 119)
 * \param double vac_wl (in Angstrom)
 * \return double air_wl (in Angstrom)
 */
double convertVacuumToAirWavelength(double vac_wl);

/*
 * double planck(double T, double wl)
 * \brief This function returns the spectral radiance of a black body in photons/(s m^2 dlambda).
 * \param Temperature is a double input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a double input that represents the wavelength at which the black body is observed in nanometers
 * \return double value for the spectral radiance of the black body in ph/(s m^2 dlambda)
 */
double PlanckFunction(double T, double wl);

/*
 * double calculateBlackBodyVFlux(double Temperature)
 * \brief This function calculates the flux (in ph/(m^2 s m)) for a Black Body
 * \param double Temperature (in Kelvin)
 * \return double flux (in ph/(m^2 s m))
 */
double calculateBlackBodyVFlux(double Temperature);

/*
 * double IntegrateSpectralElementOfBlackBody(double wl0, double wlf, double T)
 * \brief This function returns the integrated spectral flux of a black body in photons/(s m^2) for a given
 * \brief wavelength range using the simpson method.
 * \param Temperature is a double input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a double input that represents the wavelength at which the black body is observed in nanometers
 * \return double value for the spectral radiance of the black body in ph/(s m^2 dlambda)
 */
double IntegrateSpectralElementOfBlackBody(double wl0, double wlf, double T);

double getFactorToMatchFluxesBetweenElements(operaSpectralElements *refElements,operaSpectralElements *elementsToMatch, double delta_wl);

operaSpectrum calculateUniformSample(const operaSpectrum& input, unsigned npout);

float getFluxAtWavelength(unsigned np,float *wl,float *flux,float wavelengthForNormalization);
double getFluxAtWavelength(const operaVector& wavelength, const operaVector& flux, double target);
operaVector getFluxesAtWavelength(const operaSpectrum& spectrum, double target);

unsigned readContinuumWavelengthMask(string wavelength_mask, double *wl0, double *wlf);
operaWavelengthRanges readContinuumWavelengthMask(string wavelength_mask);
operaWavelengthRanges readIndexedWavelengthMask(string wavelength_mask);


unsigned getSpectrumWithinWLRange(operaSpectralElements *inputSpectrum, double wl0, double wlf, double *outputFlux, double *outputWavelength);

bool getOverlappingWLRange(operaSpectralElements *refElements, operaSpectralElements *elementsToMatch, double &wl0, double &wlf);

operaWavelengthRanges getWavelengthMaskAroundLines(const operaSpectrum sourceLines, double spectralResolution, double nsig);

double calculateDeltaRadialVelocityInKPS(double telluricWL, double observedWL);

/*!
 * \brief Resamples a spectrum to a new set of wavelengths.
 * \param inputWavelength The wavelengths of the input spectrum
 * \param inputFlux The flux of the input spectrum
 * \param outputWavelength The wavelengths to resample to
 * \return The flux of the fit spectrum at each point along outputWavelength
 */
operaVector fitSpectrum(const operaVector& inputWavelength, const operaVector& inputFlux, const operaVector& outputWavelength);

/*!
 * \brief Resamples a spectrum to a new set of wavelengths.
 * \param inputSpectrum Input operaSpectrum
 * \param outputWavelength The wavelengths to resample to
 * \return The fit spectrum at each point along outputWavelength
 */
operaSpectrum fitSpectrum(const operaSpectrum& inputSpectrum, const operaVector& outputWavelength);

/*!
 * \brief Returns the subset of an input spectrum, which lies within the mask ranges
 * \param string The wavelength ranges to return output spectrum
 * \param inputSpectrum Input operaSpectrum
 * \return Output operaSpectrum
 */
operaSpectrum getSpectrumWithinMask(string wavelength_mask, const operaSpectrum& inputSpectrum);

/*!
 * \brief Returns the subset of an input spectrum, which lies within a given opera wavelength range
 * \param operaWavelengthRange The wavelength range to return output spectrum
 * \param inputSpectrum Input operaSpectrum
 * \return Output operaSpectrum
 */
template <class T, class W>
void appendSpectrumWithinRange(W wlrange, const T& inputSpectrum, T& outputSpectrum) {
	for(unsigned i=0; i<inputSpectrum.size(); i++) {
        if(wlrange.contains(inputSpectrum.getwavelength(i))) {
            outputSpectrum.insertfrom(inputSpectrum, i);
        }
    }
}

template <class T, class W>
T getSpectrumWithinRange(W wlrange, const T& inputSpectrum) {
    T outputSpectrum;
    appendSpectrumWithinRange(wlrange, inputSpectrum, outputSpectrum);
    return outputSpectrum;
}

operaSpectralLineList getAllSpectralLines(const operaSpectralLines& spectralLines);

operaSpectralLineList getSpectralLinesInWavelengthRange(const operaSpectralLines& spectralLines, operaWavelengthRange wlrange);

operaSpectralLineList getSpectralLinesNearMedianWidth(const operaSpectralLines& spectralLines, double nsigclip);

operaSpectralLines DetectSpectralLines(operaSpectralElements* spectralElements, double linewidth, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, dispersionaxis_t dispersiontype);

/*!
 * \brief Returns a masked spectrum, where it removes points at +/- nsig x resolution element away from line centers, from list of lines provided
 * \param operaSpectrum The input spectrum
 * \param operaSpectrum telluricLines
 * \param double spectral resolution
 * \param double nsig -- define size of mask around each line in units of resolution element.
 * \return Output operaSpectrum
 */
operaSpectrum maskSpectrumAroundLines(const operaSpectrum inputSpectrum, const operaSpectrum telluricLines, double spectralResolution, double nsig);

operaVector convolveSpectrum(operaSpectrum inputSpectrum, double spectralResolution);

operaVector fitSpectrumToPolynomial(const operaVector& inputWavelength, const operaVector& inputFlux, const operaVector& outputWavelength, unsigned order);

void LinearFit(const operaVector& x, const operaVector& y, double& a, double& b, double& absdev);

void LinearFit(const operaVector& x, const operaVector& y, double& a, double& aError, double& b, double& bError, double& absdev);

Polynomial PolynomialFit(const operaVector& x, const operaVector& y, operaVector initcoeffs);

Polynomial PolynomialFit(const operaVector& x, const operaVector& y, const operaVector& yerr, operaVector initcoeffs, operaVector initcoefferrors);

unsigned findClosestInSortedRange(const operaVector& vector, double target, unsigned startindex, unsigned endindex);

void measureFluxContinuum(const operaFluxVector& uncalibratedFlux, unsigned binsize, unsigned nsigcut, operaVector& continuumElemSamples, operaVector& continuumFluxSamples);

#endif
