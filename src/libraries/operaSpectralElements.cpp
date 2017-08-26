/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralElements
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

#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralElements.h"

/*!
 * operaSpectralElements
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates a vector of spectral elements.
 * \file operaSpectralElements.cpp
 * \ingroup libraries
 */

using namespace std;

operaSpectralElements::operaSpectralElements() :
nSpectralElements(0),
elementHeight(0),
SpectrumType(None),
hasRawFlux(false),
hasStandardFlux(false),
hasOptimalFlux(false),
hasOperaOptimalFlux(false),
hasXCorrelation(false),
hasWavelength(false),
hasDistance(false),
hasFluxSNR(false),
hasExtendedBeamFlux(false)  
{
}

operaSpectralElements::operaSpectralElements(unsigned nElements):
nSpectralElements(0),
elementHeight(0),
SpectrumType(None),
hasRawFlux(false),
hasStandardFlux(false),
hasOptimalFlux(false),
hasOperaOptimalFlux(false),
hasXCorrelation(false),
hasWavelength(false),
hasDistance(false),
hasFluxSNR(false),
hasExtendedBeamFlux(false)
{
	resize(nElements);
}

operaSpectralElements::operaSpectralElements(unsigned nElements, operaSpectralOrder_t format, bool extended):
nSpectralElements(0),
elementHeight(0),
SpectrumType(None),
hasRawFlux(false),
hasStandardFlux(false),
hasOptimalFlux(false),
hasOperaOptimalFlux(false),
hasXCorrelation(false),
hasWavelength(false),
hasDistance(false),
hasFluxSNR(false),
hasExtendedBeamFlux(false)
{
	hasExtendedBeamFlux = extended;
	resize(nElements);
	SpectrumType = format;
	switch (format) {
		case RawSpectrum:
		case RawBeamSpectrum:
			setHasRawSpectrum(true);
			break;
		case StandardSpectrum:
		case StandardBeamSpectrum:
			setHasStandardSpectrum(true);
			break;
		case OptimalSpectrum:
		case OptimalBeamSpectrum:
			setHasOptimalSpectrum(true);
			break;
		case OperaOptimalSpectrum:
		case OperaOptimalBeamSpectrum:
			setHasOperaOptimalSpectrum(true);
			break;
		default:
			break;
	}
}

void operaSpectralElements::createExtendedVectors() {
	if(hasExtendedBeamFlux) throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    hasExtendedBeamFlux = true;
    tell.resize(nSpectralElements);
    rvel.resize(nSpectralElements);
	rawFlux.resize(nSpectralElements);
	normalizedFlux.resize(nSpectralElements);
	fcalFlux.resize(nSpectralElements);
}

void operaSpectralElements::resize(unsigned nElements) {
	fluxvector.resize(nElements);
	photoCenterX.resize(nElements);
	photoCenterY.resize(nElements);
	distd.resize(nElements);
	XCorrelation.resize(nElements);
	fluxSNR.resize(nElements);
	wavelength.resize(nElements);
	if(hasExtendedBeamFlux) {
		tell.resize(nElements);
		rvel.resize(nElements);
		rawFlux.resize(nElements);
		normalizedFlux.resize(nElements);
		fcalFlux.resize(nElements);
	}
	nSpectralElements = nElements;
}

void operaSpectralElements::trim(operaIndexRange range) {
	fluxvector.trim(range);
	photoCenterX.trim(range);
	photoCenterY.trim(range);
	distd.trim(range);
	XCorrelation.trim(range);
	fluxSNR.trim(range);
	wavelength.trim(range);
	if(hasExtendedBeamFlux) {
		tell.trim(range);
		rvel.trim(range);
		rawFlux.trim(range);
		normalizedFlux.trim(range);
		fcalFlux.trim(range);
	}
	nSpectralElements = range.size();
}

unsigned operaSpectralElements::getnSpectralElements(void) const { 
	return nSpectralElements;
}

const operaFluxVector& operaSpectralElements::getFluxVector(void) const {
	return fluxvector;
}

const operaVector& operaSpectralElements::getXCorrelation() const {
	return XCorrelation;
}

const operaVector& operaSpectralElements::getPhotoCenterX() const {
	return photoCenterX;
}

const operaVector& operaSpectralElements::getPhotoCenterY() const {
	return photoCenterY;
}

const operaVector& operaSpectralElements::getDistd() const {
	return distd;
}

const operaVector& operaSpectralElements::getWavelength(void) const {
	return wavelength;
}

const operaVector& operaSpectralElements::getFluxSNR() const {
	return fluxSNR;
}

const operaVector& operaSpectralElements::getRvel() const {
	return rvel;
}	

void operaSpectralElements::setFluxVector(const operaFluxVector &FluxVector) {
	if(fluxvector.getlength() != FluxVector.getlength()) throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	fluxvector = FluxVector;
}

void operaSpectralElements::setXCorrelation(const operaVector& XCorr) {
	if (XCorrelation.size() != XCorr.size()) throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	XCorrelation = XCorr;
}

void operaSpectralElements::setPhotoCenter(const operaVector& PhotoCenterX, const operaVector& PhotoCenterY) {
	if (photoCenterX.size() != PhotoCenterX.size() || photoCenterY.size() != PhotoCenterY.size()) throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	photoCenterX = PhotoCenterX; photoCenterY = PhotoCenterY;
}

void operaSpectralElements::setDistd(const operaVector& Distd) {
	if (distd.size() != Distd.size()) throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	distd = Distd;
}

void operaSpectralElements::setWavelength(const operaVector& Wavelength) {
	if (wavelength.size() != Wavelength.size()) throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	wavelength = Wavelength;
}

void operaSpectralElements::setFluxSNR(const operaVector& FluxSNR) {
	if (fluxSNR.size() != FluxSNR.size()) throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	fluxSNR = FluxSNR;
}

void operaSpectralElements::setRvel(const operaVector& Rvel) {
	if (rvel.size() != Rvel.size()) throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	rvel = Rvel;
}

void operaSpectralElements::calculateFluxSNR() {
	fluxSNR = Abs(fluxvector.getflux() / Sqrt(fluxvector.getvariance()));
}

void operaSpectralElements::setwavelengthsFromCalibration(const operaWavelength *Wavelength) {
    for(unsigned indexElem=0;indexElem<nSpectralElements;indexElem++) {
        wavelength[indexElem] = Wavelength->evaluateWavelength(distd[indexElem]);
    }
    setHasWavelength(true);
}

operaIndexRange operaSpectralElements::getContainedWavelengthSubrangeIndexes(double wl0, double wlf) const {
	return wavelength.subrange(wl0, wlf);
}

double operaSpectralElements::getFlux(unsigned indexElem) const {
	return fluxvector.getflux(indexElem);
}

void operaSpectralElements::setFlux(double Flux, unsigned indexElem) {
	fluxvector.setflux(Flux, indexElem);
}

double operaSpectralElements::getFluxVariance(unsigned indexElem) const {
	return fluxvector.getvariance(indexElem);
}

void operaSpectralElements::setFluxVariance(double FluxVariance, unsigned indexElem) {
	fluxvector.setvariance(FluxVariance, indexElem);
}

double operaSpectralElements::getFluxSNR(unsigned indexElem) const {
	return fluxSNR[indexElem];
}

void operaSpectralElements::setFluxSNR(double FluxSNR, unsigned indexElem) {
	fluxSNR[indexElem] = FluxSNR;
}

double operaSpectralElements::getphotoCenterX(unsigned indexElem) const {
	return photoCenterX[indexElem];
}
double operaSpectralElements::getphotoCenterY(unsigned indexElem) const {
	return photoCenterY[indexElem];
}
void operaSpectralElements::setphotoCenter(double x, double y, unsigned indexElem) {
	photoCenterX[indexElem] = x;
	photoCenterY[indexElem] = y;	
}

double operaSpectralElements::getdistd(unsigned indexElem) const {
	return distd[indexElem];
}
void operaSpectralElements::setdistd(double Distd, unsigned indexElem) {
	distd[indexElem] = Distd;
}

double operaSpectralElements::getwavelength(unsigned indexElem) const {
	return wavelength[indexElem];
}

void operaSpectralElements::setwavelength(double Wavelength, unsigned indexElem) {
	wavelength[indexElem] = Wavelength;
}

double operaSpectralElements::getXCorrelation(unsigned indexElem) const {
	return XCorrelation[indexElem];
}

void operaSpectralElements::setXCorrelation(double Xcorr, unsigned indexElem){
	XCorrelation[indexElem] = Xcorr;
}

double operaSpectralElements::gettell(unsigned indexElem) const { 
	return tell[indexElem];
}

void operaSpectralElements::settell(double value, unsigned indexElem) { 
	tell[indexElem] = value;
}

void operaSpectralElements::copyTOtell(void) {
	tell = wavelength;
}

void operaSpectralElements::copyFROMtell(void) {
	wavelength = tell;
}

double operaSpectralElements::getrvel(unsigned indexElem) const { 
	return rvel[indexElem];
}

void operaSpectralElements::setrvel(double value, unsigned indexElem) { 
	rvel[indexElem] = value;
}

double operaSpectralElements::getnormalizedFlux(unsigned indexElem) const { 
	return normalizedFlux.getflux(indexElem);
}

void operaSpectralElements::setnormalizedFlux(double value, unsigned indexElem) { 
	normalizedFlux.setflux(value, indexElem);
}

double operaSpectralElements::getnormalizedFluxVariance(unsigned indexElem) const { 
	return normalizedFlux.getvariance(indexElem);
}

void operaSpectralElements::setnormalizedFluxVariance(double value, unsigned indexElem) {
	normalizedFlux.setvariance(value, indexElem);
}

void operaSpectralElements::copyTOnormalizedFlux(void) { 
    normalizedFlux = fluxvector;
}

void operaSpectralElements::copyFROMnormalizedFlux(void) { 
    fluxvector = normalizedFlux;
}

double operaSpectralElements::getfcalFlux(unsigned indexElem) const { 
	return fcalFlux.getflux(indexElem);
}

void operaSpectralElements::setfcalFlux(double value, unsigned indexElem) {
	fcalFlux.setflux(value, indexElem);
}

double operaSpectralElements::getfcalFluxVariance(unsigned indexElem) const { 
	return fcalFlux.getvariance(indexElem);
}

void operaSpectralElements::setfcalFluxVariance(double value, unsigned indexElem) {
	fcalFlux.setvariance(value, indexElem);
}

void operaSpectralElements::copyTOfcalFlux(void) {
    fcalFlux = fluxvector;
}

void operaSpectralElements::copyFROMfcalFlux(void) { 
    fluxvector = fcalFlux;
}

double operaSpectralElements::getrawFlux(unsigned indexElem) const {
	return rawFlux.getflux(indexElem);
}

void operaSpectralElements::setrawFlux(double value, unsigned indexElem) {
	rawFlux.setflux(value, indexElem);
}

double operaSpectralElements::getrawFluxVariance(unsigned indexElem) const {
	return rawFlux.getvariance(indexElem);
}

void operaSpectralElements::setrawFluxVariance(double value, unsigned indexElem) {
	rawFlux.setvariance(value, indexElem);
}

void operaSpectralElements::copyTOrawFlux(void) {
    rawFlux = fluxvector;
}

void operaSpectralElements::copyFROMrawFlux(void) {
    fluxvector = rawFlux;
}
