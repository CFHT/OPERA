/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralEnergyDistribution
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
#include "libraries/operaSpectralTools.h"
#include "libraries/operaSpectralEnergyDistribution.h"  // for operaSpectralEnergyDistribution

/*!
 * operaSpectralEnergyDistribution
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the spectral energy distribution object.
 * \file operaSpectralEnergyDistribution.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaSpectralEnergyDistribution
 * \brief Encapsulation of Wavelength information.
 * \return none
 */

/*
 * Constructors
 */

operaSpectralEnergyDistribution::operaSpectralEnergyDistribution() :
nDataPoints(0),
wavelengthForNormalization(0),
hasFluxData(false),   
hasUncalibratedFlux(false),
hasCalibratedFlux(false),
hasFluxCalibration(false),
hasInstrumentThroughput(false) 
{
}

operaSpectralEnergyDistribution::operaSpectralEnergyDistribution(unsigned NDataPoints) :
nDataPoints(0),
wavelengthForNormalization(0),
hasFluxData(false),   
hasUncalibratedFlux(false),
hasCalibratedFlux(false),
hasFluxCalibration(false),
hasInstrumentThroughput(false) 
{
    if (NDataPoints == 0) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	resizeDataVectors(NDataPoints);
}

operaSpectralEnergyDistribution::operaSpectralEnergyDistribution(unsigned NDataPoints, unsigned nElements) :
nDataPoints(0),
wavelengthForNormalization(0),
hasFluxData(false),   
hasUncalibratedFlux(false),
hasCalibratedFlux(false),
hasFluxCalibration(false),
hasInstrumentThroughput(false) 
{
    if (NDataPoints == 0 || nElements == 0) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    resizeDataVectors(NDataPoints); 
    resizeCalibrationVectors(nElements);
}

/*
 * Methods for managing data
 */  

void operaSpectralEnergyDistribution::resizeDataVectors(unsigned NDataPoints) {
	distanceData.resize(NDataPoints);
	wavelengthData.resize(NDataPoints);
	fluxData.resize(NDataPoints);
	nDataPoints = NDataPoints;
}

void operaSpectralEnergyDistribution::resizeCalibrationVectors(unsigned nElements) {
	calibrationWavelength.resize(nElements);
	calibrationDist.resize(nElements);
	uncalibratedFlux.resize(nElements);
    calibratedFlux.resize(nElements);
    fluxCalibration.resize(nElements);
    instrumentThroughput.resize(nElements);
}

unsigned operaSpectralEnergyDistribution::getnDataPoints(void) const {
    return nDataPoints;
}

void operaSpectralEnergyDistribution::setwavelengthForNormalization(double WavelengthForNormalization) {
    wavelengthForNormalization = WavelengthForNormalization;
}

double operaSpectralEnergyDistribution::getwavelengthForNormalization(void) const {
    return wavelengthForNormalization;
}

void operaSpectralEnergyDistribution::setdistanceData(double Distance, unsigned index) {
    distanceData[index] = Distance;
}

void operaSpectralEnergyDistribution::setwavelengthData(double Wavelength, unsigned index) {
    wavelengthData[index] = Wavelength;
}

void operaSpectralEnergyDistribution::setfluxData(double Flux, unsigned index) {
    fluxData[index] =  Flux;
}

double operaSpectralEnergyDistribution::getdistanceData(unsigned index) const {
    return distanceData[index];
}

double operaSpectralEnergyDistribution::getwavelengthData(unsigned index) const {
    return wavelengthData[index];
}

double operaSpectralEnergyDistribution::getfluxData(unsigned index) const {
    return fluxData[index];
}

/*
 * Methods for managing flux calibration and throughput elements
 */

void operaSpectralEnergyDistribution::setCalibrationDist(const operaVector& CalibrationDist) {
	calibrationDist = CalibrationDist;
}
	
void operaSpectralEnergyDistribution::setCalibrationWavelength(const operaVector& CalibrationWavelength) {
	calibrationWavelength = CalibrationWavelength;
}

void operaSpectralEnergyDistribution::setUncalibratedFlux(const operaFluxVector& UncalibratedFlux) {
    uncalibratedFlux = UncalibratedFlux;
}

void operaSpectralEnergyDistribution::setCalibratedFlux(const operaFluxVector& CalibratedFlux) {
    calibratedFlux = CalibratedFlux;
}

void operaSpectralEnergyDistribution::setFluxCalibration(const operaFluxVector& FluxCalibration) {
    fluxCalibration = FluxCalibration;
}

void operaSpectralEnergyDistribution::setThroughput(const operaFluxVector& Throughput) {
    instrumentThroughput = Throughput;
}

/*
 * Other Methods
 */

void operaSpectralEnergyDistribution::measureUncalibratedContinuum(unsigned binsize, unsigned nsigcut) {
    if(!getHasUncalibratedFlux()) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    resizeDataVectors(0);
    operaVector continuumElemSample;
    measureFluxContinuum(uncalibratedFlux, binsize, nsigcut, continuumElemSample, fluxData);
    for(unsigned i=0; i<continuumElemSample.size(); i++){
        unsigned sampledIndex = (unsigned)continuumElemSample[i];
        distanceData.insert(calibrationDist[sampledIndex]);
        wavelengthData.insert(calibrationWavelength[sampledIndex]);
        setwavelengthForNormalization(wavelengthData.last()); // DT Sep 6 2013 is this right?
	}
	
    if(!fluxData.empty()) setHasFluxData(true);
}

void operaSpectralEnergyDistribution::populateUncalibratedFluxFromContinuumData() {
    if(!getHasFluxData()) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    uncalibratedFlux = fitSpectrum(wavelengthData, fluxData, calibrationWavelength);
    
    setHasUncalibratedFlux(true);
}

void operaSpectralEnergyDistribution::calculateUncalibratedFlux(unsigned binsize, unsigned nsigcut) {
    measureUncalibratedContinuum(binsize,nsigcut);
    populateUncalibratedFluxFromContinuumData();
}
