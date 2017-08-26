/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralOrder
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

#include <iostream>
#include <fstream>
#include <cmath>

#include "globaldefines.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaGeometry.h"
#include "libraries/operaWavelength.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaSpectralLines.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaInstrumentProfile.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/PixelSet.h"
#include "libraries/Gaussian.h"
#include "libraries/operaArgumentHandler.h"

#define DEBUG false
#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef SATURATIONLIMIT
#define SATURATIONLIMIT 65535  // this should be retrieved from the config/param file
#endif
#ifndef MAXNUMBEROFLINES
#define MAXNUMBEROFLINES 10000
#endif
#ifndef NSIGMACUT
#define NSIGMACUT 3
#endif
#ifndef DEFAULT_SNR_SMOOTHING
#define DEFAULT_SNR_SMOOTHING 5
#endif

#include "libraries/operaLib.h"     // for itos
#include "libraries/operaStats.h"
#include "libraries/operaFit.h"
#include "libraries/operaCCD.h"
#include "libraries/operaMath.h"
#include "libraries/ladfit.h"
#include "libraries/VLArray.h"

/*!
 * operaSpectralOrder
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the FITS image.
 * \file operaSpectralOrder.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \brief operaSpectralOrder
 * \details A spectral order (SO) consists of a data set containing 
 * \details the information concerned with a full spectral order.
 */

/*
 * Constructors
 */
operaSpectralOrder::operaSpectralOrder() : 
orderNumber(0),
SpectrumType(None),
SNR(0.0),
SpectralElements(NULL),
SkyElements(NULL),
Geometry(NULL),
Wavelength(NULL),
InstrumentProfile(NULL),
SpectralLines(NULL),
Polarimetry(NULL),
SpectralEnergyDistribution(NULL),
numberOfBeams(0),
hasSpectralElements(false), 
hasSkyElements(false),
hasGeometry(false), 
hasWavelength(false), 
hasInstrumentProfile(false), 
hasSpectralLines(false), 
hasExtractionApertures(false),
hasPolarimetry(false),
hasSNR(false),
hasCenterSNROnly(false),
hasSpectralEnergyDistribution(false),
hasWavelengthRange(false)
{
    tiltInDegrees.value = 0.0;
    tiltInDegrees.error = 0.0;
    
    for (unsigned backgroundIndex = 0 ; backgroundIndex < LEFTANDRIGHT ; backgroundIndex++) {
        BackgroundElements[backgroundIndex] = NULL;
        BackgroundApertures[backgroundIndex] = NULL;
    }
    
    for (unsigned beam = 0 ; beam < MAXNUMBEROFBEAMS ; beam++) {
        BeamElements[beam] = NULL;
        BeamProfiles[beam] = NULL;
        ExtractionApertures[beam] = NULL;
    }
}

operaSpectralOrder::operaSpectralOrder(unsigned order) : 
SpectrumType(None),
SNR(0.0),
SpectralElements(NULL),
SkyElements(NULL),
Geometry(NULL),
Wavelength(NULL),
InstrumentProfile(NULL),
SpectralLines(NULL),
Polarimetry(NULL),
SpectralEnergyDistribution(NULL),
numberOfBeams(0),
hasSpectralElements(false), 
hasSkyElements(false),
hasGeometry(false), 
hasWavelength(false), 
hasInstrumentProfile(false), 
hasSpectralLines(false), 
hasExtractionApertures(false),
hasPolarimetry(false),
hasSNR(false),
hasCenterSNROnly(false),
hasSpectralEnergyDistribution(false),
hasWavelengthRange(false)
{ 	
	orderNumber = order;
    tiltInDegrees.value = 0.0;
    tiltInDegrees.error = 0.0;
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
        BackgroundElements[backgroundIndex] = NULL;//new operaSpectralElements();
        BackgroundApertures[backgroundIndex] = NULL;//new operaExtractionAperture();  
    }
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
        BeamElements[beam] = NULL;//new operaSpectralElements();
        ExtractionApertures[beam] = NULL;//new operaExtractionAperture();
        BeamProfiles[beam] = NULL;//new operaInstrumentProfile();
		BeamSED[beam] = NULL;
    }
}

operaSpectralOrder::operaSpectralOrder(unsigned order, unsigned maxdatapoints, unsigned maxValues, unsigned nElements, operaSpectralOrder_t format)  : 
SpectrumType(None),
SNR(0.0),
SpectralElements(NULL),
SkyElements(NULL),
Geometry(NULL),
Wavelength(NULL),
InstrumentProfile(NULL),
SpectralLines(NULL),
Polarimetry(NULL),
SpectralEnergyDistribution(NULL),
numberOfBeams(0),
hasSpectralElements(false), 
hasSkyElements(false), 
hasGeometry(false), 
hasWavelength(false), 
hasInstrumentProfile(false), 
hasSpectralLines(false), 
hasExtractionApertures(false),
hasPolarimetry(false),
hasSNR(false),
hasCenterSNROnly(false),
hasSpectralEnergyDistribution(false),
hasWavelengthRange(false)
{
	if (maxdatapoints == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	// This is OK, it just means you don't need spectralElements...
	if (nElements == 0) {
		//throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (maxValues == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	orderNumber = order;
	SpectrumType = format;
    tiltInDegrees.value = 0.0;
    tiltInDegrees.error = 0.0;
    
	Geometry = new operaGeometry(maxdatapoints, maxValues);
	if (nElements) {
		SpectralElements = new operaSpectralElements(nElements, SpectrumType);
        SkyElements = new operaSpectralElements(nElements, SpectrumType);        
    }
	if (nElements)
		Polarimetry = new operaPolarimetry(nElements);
#if 0   
	if (nElements)
		for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
			BackgroundElements[backgroundIndex] = new operaSpectralElements(nElements, SpectrumType);
			BackgroundApertures[backgroundIndex] = new operaExtractionAperture<Line>();  
		}
    
	if (nElements)
		for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
			BeamElements[beam] = new operaSpectralElements(nElements, SpectrumType);
			ExtractionApertures[beam] = new operaExtractionAperture<Line>();  
			BeamProfiles[beam] = new operaInstrumentProfile();        
		}    
#else
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
        BackgroundElements[backgroundIndex] = NULL;
        BackgroundApertures[backgroundIndex] = NULL;  
    }
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
        BeamElements[beam] = NULL;
        ExtractionApertures[beam] = NULL;
        BeamProfiles[beam] = NULL;
		BeamSED[beam] = NULL;
    }
#endif
}

/*
 * Destructor
 */
operaSpectralOrder::~operaSpectralOrder() {
	deleteAll();    
}

/*
 * Common Methods
 */

void operaSpectralOrder::deleteAll() {
	
	if (Geometry && hasGeometry) 
		delete Geometry;
	Geometry = NULL;
	if (Wavelength && hasWavelength) 
		delete Wavelength;
	Wavelength = NULL;
	if (SpectralElements && hasSpectralElements) 
		delete SpectralElements;
	SpectralElements = NULL;
	if (SkyElements && hasSkyElements) 
		delete SkyElements;
	SkyElements = NULL;    
	if (InstrumentProfile && hasInstrumentProfile) 
		delete InstrumentProfile;
	InstrumentProfile = NULL;
	if (SpectralLines && hasSpectralLines) 
		delete SpectralLines;
	SpectralLines = NULL;
    if (Polarimetry && hasPolarimetry) 
		delete Polarimetry;
	Polarimetry = NULL;
    if (SpectralEnergyDistribution && hasSpectralEnergyDistribution) 
		delete SpectralEnergyDistribution;
	SpectralEnergyDistribution = NULL;
	
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
		if (BackgroundElements[backgroundIndex]) 
			delete BackgroundElements[backgroundIndex];
		BackgroundElements[backgroundIndex] = NULL;
		
		if (BackgroundApertures[backgroundIndex])
			delete BackgroundApertures[backgroundIndex];
		BackgroundApertures[backgroundIndex] = NULL;
    }
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
		if (BeamElements[beam]) 
			delete BeamElements[beam];
		BeamElements[beam] = NULL;
		if (ExtractionApertures[beam]) 
			delete ExtractionApertures[beam];
		ExtractionApertures[beam] = NULL;
		if (BeamProfiles[beam]) 
			delete BeamProfiles[beam];
		BeamProfiles[beam] = NULL;
        if (BeamSED[beam]) 
            delete BeamSED[beam];
        BeamSED[beam] = NULL;        
    }    
	hasSpectralElements = false; 
	hasSkyElements = false; 
	hasGeometry = false; 
	hasWavelength = false; 
	hasInstrumentProfile = false; 
	hasSpectralLines = false; 
	hasExtractionApertures = false;
	hasPolarimetry = false;
	hasSNR = false;
	hasSpectralEnergyDistribution = false;
	hasWavelengthRange = false;
}

void operaSpectralOrder::deleteBeamProfiles(void) {
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
		if (BeamProfiles[beam]) 
			delete BeamProfiles[beam];
		BeamProfiles[beam] = NULL;
    }       
}

void operaSpectralOrder::deleteApertures(void) {
    
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
		if (BackgroundApertures[backgroundIndex]) 
            delete BackgroundApertures[backgroundIndex];
        BackgroundApertures[backgroundIndex] = NULL;
    }    
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
		if (ExtractionApertures[beam]) 
			delete ExtractionApertures[beam];
		ExtractionApertures[beam] = NULL;
    }        
}


void operaSpectralOrder::deleteBeamsAndBackgroundElements(void) {
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
		if (BackgroundElements[backgroundIndex]) 
			delete BackgroundElements[backgroundIndex];
		BackgroundElements[backgroundIndex] = NULL;
    }
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
		if (BeamElements[beam]) 
			delete BeamElements[beam];
		BeamElements[beam] = NULL;
    }    
}

void operaSpectralOrder::createGeometry(unsigned maxdatapoints, unsigned maxValues) {
	if (Geometry) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	Geometry = new operaGeometry(maxdatapoints, maxValues);
}

void operaSpectralOrder::createWavelength(unsigned maxnumberofcoefficients) {
	if (Wavelength) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	Wavelength = new operaWavelength(maxnumberofcoefficients);    
}

void operaSpectralOrder::createSpectralElements(unsigned maxdatapoints, operaSpectralOrder_t SpectrumType, bool extended) {
	if (SpectralElements != NULL) {
		// CU Aug 20, 2015 -- Swapped which line is commented out. If this is causing a crash, that's sign that something else is wrong.
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		//SpectralElements->Resizevectors(maxdatapoints, SpectralElements->getSpectrumType());
	} else {
		SpectralElements = new operaSpectralElements(maxdatapoints, SpectrumType, extended);
		createSpectralEnergyDistribution();
	}
}

void operaSpectralOrder::createSkyElements(unsigned maxdatapoints, operaSpectralOrder_t SpectrumType) {
	if (SkyElements) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	SkyElements = new operaSpectralElements(maxdatapoints, SpectrumType);
}

void operaSpectralOrder::createBeamsAndBackgrounds(unsigned nElements, unsigned nBeams, operaSpectralOrder_t format, bool extendedbeams) {
	if (nElements == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	/*
     * E. Martioli Jan 26, 2016 - I commented this line because it was crashing when reading spc files
        When nBeams=0 should not be a problem (I think) since all loops containing beams will not initialize.
     if (nBeams == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	*/
	for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
		if (BackgroundElements[backgroundIndex] != NULL) {
			BackgroundElements[backgroundIndex]->resize(nElements);
			BackgroundElements[backgroundIndex]->setSpectrumType(format);
		} else {
			BackgroundElements[backgroundIndex] = new operaSpectralElements(nElements, format);
		}
    }
    
    for(unsigned beam=0;beam<nBeams; beam++) {
		if (BeamElements[beam] != NULL) {
			BeamElements[beam]->resize(nElements);
			BeamElements[beam]->setSpectrumType(format);
		} else {
			BeamElements[beam] = new operaSpectralElements(nElements, format, extendedbeams);
		}
    }    
}

void operaSpectralOrder::createSpectralEnergyDistribution() {
	if(SpectralEnergyDistribution) delete SpectralEnergyDistribution;
    SpectralEnergyDistribution = new operaSpectralEnergyDistribution;
    for(unsigned beam=0; beam < MAXNUMBEROFBEAMS; beam++) {
        if(BeamSED[beam]) delete BeamSED[beam];
        BeamSED[beam] = new operaSpectralEnergyDistribution;
    }
}


void operaSpectralOrder::deletePolarimetry(void) {
    if (Polarimetry) 
		delete Polarimetry;
	Polarimetry = NULL;
}

void operaSpectralOrder::createPolarimetry(unsigned nElements) {
	if (nElements == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    deletePolarimetry();
	Polarimetry = new operaPolarimetry(nElements);
}

unsigned operaSpectralOrder::getorder(void) const {
	return orderNumber;
}

void operaSpectralOrder::setorder(unsigned ordernumber) {
	orderNumber = ordernumber;
}

void operaSpectralOrder::sethasSpectralElements(bool HasSpectralElements) {
	hasSpectralElements = HasSpectralElements;
}

void operaSpectralOrder::sethasSkyElements(bool HasSkyElements) {
	hasSkyElements = HasSkyElements;
}

void operaSpectralOrder::sethasGeometry(bool HasGeometry) {
	hasGeometry = HasGeometry;
}

void operaSpectralOrder::sethasWavelength(bool HasWavelength) {
	hasWavelength = HasWavelength;
}

void operaSpectralOrder::sethasInstrumentProfile(bool HasInstrumentProfile) {
	hasInstrumentProfile = HasInstrumentProfile;
}

void operaSpectralOrder::sethasSpectralLines(bool HasSpectralLines) {
	hasSpectralLines = HasSpectralLines;
}

void operaSpectralOrder::sethasExtractionApertures(bool HasExtractionApertures) {
	hasExtractionApertures = HasExtractionApertures;
}

void operaSpectralOrder::sethasPolarimetry(bool HasPolarimetry) {
	hasPolarimetry = HasPolarimetry;
}

void operaSpectralOrder::sethasSNR(bool HasSNR) {
	hasSNR = HasSNR;
}

void operaSpectralOrder::sethasCenterSNROnly(bool HasCenterSNROnly) {
	hasCenterSNROnly = HasCenterSNROnly;
}

void operaSpectralOrder::sethasSpectralEnergyDistribution(bool HasSpectralEnergyDistribution) {
	hasSpectralEnergyDistribution = HasSpectralEnergyDistribution;
}

bool operaSpectralOrder::gethasSpectralElements(void) const {
	return hasSpectralElements;
}

bool operaSpectralOrder::gethasSkyElements(void) const {
	return hasSkyElements;
}

bool operaSpectralOrder::gethasGeometry(void) const {
	return hasGeometry;
}

bool operaSpectralOrder::gethasWavelength(void) const {
	return hasWavelength;
}

bool operaSpectralOrder::gethasInstrumentProfile(void) const {
	return hasInstrumentProfile;
}

bool operaSpectralOrder::gethasSpectralLines(void) const {
	return hasSpectralLines;
}

bool operaSpectralOrder::gethasExtractionApertures(void) const {
	return hasExtractionApertures;
}

bool operaSpectralOrder::gethasPolarimetry(void) const {
	return hasPolarimetry;
}

bool operaSpectralOrder::gethasSNR(void) const {
	return hasSNR;
}

bool operaSpectralOrder::gethasCenterSNROnly(void) const {
	return hasCenterSNROnly;
}

bool operaSpectralOrder::gethasSpectralEnergyDistribution(void) const {
	return hasSpectralEnergyDistribution;
}

operaSpectralOrder_t operaSpectralOrder::getSpectrumType(void) const { 
	return SpectrumType;
}

void operaSpectralOrder::setSpectrumType(operaSpectralOrder_t format) { 
	SpectrumType = format;
}

/*
 * Common Methods
 */
operaSpectralElements *operaSpectralOrder::getSpectralElements() {
	return SpectralElements;
}

const operaSpectralElements *operaSpectralOrder::getSpectralElements() const {
	return SpectralElements;
}

operaSpectralElements *operaSpectralOrder::getSkyElements() {
	return SkyElements;
}

const operaSpectralElements *operaSpectralOrder::getSkyElements() const {
	return SkyElements;
}

operaGeometry *operaSpectralOrder::getGeometry() {
	return Geometry;
}

const operaGeometry *operaSpectralOrder::getGeometry() const {
	return Geometry;
}

operaWavelength *operaSpectralOrder::getWavelength() {
	return Wavelength;
}

const operaWavelength *operaSpectralOrder::getWavelength() const {
	return Wavelength;
}

operaInstrumentProfile *operaSpectralOrder::getInstrumentProfile() {
	return InstrumentProfile;
}

const operaInstrumentProfile *operaSpectralOrder::getInstrumentProfile() const {
	return InstrumentProfile;
}

operaSpectralLines *operaSpectralOrder::getSpectralLines() {
	return SpectralLines;
}

const operaSpectralLines *operaSpectralOrder::getSpectralLines() const {
	return SpectralLines;
}

operaPolarimetry *operaSpectralOrder::getPolarimetry() {
	return Polarimetry;
}

const operaPolarimetry *operaSpectralOrder::getPolarimetry() const {
	return Polarimetry;
}

operaSpectralEnergyDistribution *operaSpectralOrder::getSpectralEnergyDistribution(void) {
	return SpectralEnergyDistribution;
}

const operaSpectralEnergyDistribution *operaSpectralOrder::getSpectralEnergyDistribution(void) const {
	return SpectralEnergyDistribution;
}

double operaSpectralOrder::getCenterSNR(void) const {
	return SNR;
}

void operaSpectralOrder::setCenterSNR(double Snr) {
	SNR = Snr;
}

unsigned operaSpectralOrder::getnumberOfBeams(void) const {
	return numberOfBeams;
}   

void operaSpectralOrder::setnumberOfBeams(unsigned NumberOfBeams) {
#ifdef RANGE_CHECK
	if (NumberOfBeams >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	numberOfBeams = NumberOfBeams;
}

doubleValue_t operaSpectralOrder::getTiltInDegrees(void) const {
	return tiltInDegrees;
}

void operaSpectralOrder::setTiltInDegrees(doubleValue_t TiltInDegrees) {
	tiltInDegrees = TiltInDegrees;
}

void operaSpectralOrder::setTiltInDegrees(double tilt, double error) {
	tiltInDegrees.value = tilt;
	tiltInDegrees.error = error;    
}

double operaSpectralOrder::getTiltInDegreesValue(void) const {
	return tiltInDegrees.value;    
}

double operaSpectralOrder::getTiltInDegreesError(void) const {
	return tiltInDegrees.error;    
}

operaSpectralElements *operaSpectralOrder::getBeamElements(unsigned beam) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamElements[beam];
}

const operaSpectralElements *operaSpectralOrder::getBeamElements(unsigned beam) const {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamElements[beam];
}

operaInstrumentProfile *operaSpectralOrder::getBeamProfiles(unsigned beam) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamProfiles[beam]; 
}

const operaInstrumentProfile *operaSpectralOrder::getBeamProfiles(unsigned beam) const {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamProfiles[beam]; 
}

operaSpectralElements *operaSpectralOrder::getBackgroundElements(unsigned LeftOrRight) {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BackgroundElements[LeftOrRight];
}

const operaSpectralElements *operaSpectralOrder::getBackgroundElements(unsigned LeftOrRight) const {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BackgroundElements[LeftOrRight];
}

operaExtractionAperture<Line> *operaSpectralOrder::getExtractionApertures(unsigned beam) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return ExtractionApertures[beam];
}

const operaExtractionAperture<Line> *operaSpectralOrder::getExtractionApertures(unsigned beam) const {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return ExtractionApertures[beam];
}

operaExtractionAperture<Line> *operaSpectralOrder::getBackgroundApertures(unsigned LeftOrRight) {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BackgroundApertures[LeftOrRight];
}

const operaExtractionAperture<Line> *operaSpectralOrder::getBackgroundApertures(unsigned LeftOrRight) const {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BackgroundApertures[LeftOrRight];
}

void operaSpectralOrder::setBeamElements(unsigned beam, operaSpectralElements *beamElements) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (BeamElements[beam]) {
		delete BeamElements[beam];
	}
	BeamElements[beam] = beamElements;
}

void operaSpectralOrder::setBeamProfiles(unsigned beam, operaInstrumentProfile *beamProfiles){
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
  	if (BeamProfiles[beam]) {
		delete BeamProfiles[beam];
	}
	BeamProfiles[beam] = beamProfiles;  
}

void operaSpectralOrder::setBackgroundElements(unsigned LeftOrRight, operaSpectralElements *backgroundElements) {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (BackgroundElements[LeftOrRight]) {
		delete BackgroundElements[LeftOrRight];
	}
	BackgroundElements[LeftOrRight] = backgroundElements;
}

void operaSpectralOrder::setExtractionApertures(unsigned beam, operaExtractionAperture<Line> *extractionApertures) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (ExtractionApertures[beam]) {
		delete ExtractionApertures[beam];
	}
	ExtractionApertures[beam] = extractionApertures;
}

void operaSpectralOrder::setBackgroundApertures(unsigned LeftOrRight, operaExtractionAperture<Line> *backgroundApertures) {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (BackgroundApertures[LeftOrRight]) {
		delete BackgroundApertures[LeftOrRight];
	}
	BackgroundApertures[LeftOrRight] = backgroundApertures;
}

double operaSpectralOrder::getminwavelength() const {
	return minwavelength;
}

double operaSpectralOrder::getmaxwavelength() const {
	return maxwavelength;
}

bool operaSpectralOrder::gethasWavelengthRange() const {
	return hasWavelengthRange;
}

void operaSpectralOrder::setminwavelength(double wl) {
	minwavelength = wl;
}

void operaSpectralOrder::setmaxwavelength(double wl) {
	maxwavelength = wl;
}

void operaSpectralOrder::sethasWavelengthRange(bool hasRange) {
	hasWavelengthRange = hasRange;
}

double operaSpectralOrder::getsnrSpectralBinSize() const {
	return snrSpectralBinSize;
}

void operaSpectralOrder::setsnrSpectralBinSize(double spectralbinsize) {
	snrSpectralBinSize = spectralbinsize;
}

void operaSpectralOrder::NormalizeFlat(operaFITSImage &flatMatrix, operaFITSImage &outputMatrix, unsigned nx, unsigned ny, unsigned binsize){
	
	// This function only works for spectralElementHeight=1; IPysize=1; IPyampling=1, so we enforce these below.
	if (InstrumentProfile->getysize()!=1 || InstrumentProfile->getYsampling()!=1 || SpectralElements->getelementHeight()!=1){
		throw operaException("operaSpectralOrder: instrument profile ysize, ysampling and spectralelement height must be 1, not "+itos(InstrumentProfile->getysize())+" and "+dtos(SpectralElements->getelementHeight())+" and "+dtos(SpectralElements->getelementHeight())+" (resp.) as given.",operaErrorInstrumentProfileImproperDimensions, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (SpectralElements->getnSpectralElements() == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned NXPoints = InstrumentProfile->getxsize()*InstrumentProfile->getXsampling();	
	
	unsigned NumberofElementsToBin = binsize;
	
	// DT Jan 213 moved out of loop
	float *SubPixElemSampleMedianCounts = (float *)malloc(NXPoints*sizeof(float));
	if (!SubPixElemSampleMedianCounts) {
		throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElementCounts = (float *)malloc(NumberofElementsToBin*sizeof(float));	
	if (!SubPixElementCounts) {
		throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixXcoords = (float *)malloc(NXPoints*sizeof(float));
	if (!SubPixXcoords) {
		throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *y2 = (float *)malloc(NXPoints*sizeof(float));
	if (!y2) {
		throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	// Below it loops over all spectral elements
	for (unsigned indexElem=0; indexElem < SpectralElements->getnSpectralElements(); indexElem++) {		
		
		// Figure out which is the first and last elements given the number of elements to bin
		unsigned firstElement = indexElem - (unsigned)ceil(NumberofElementsToBin/2);
		unsigned lastElement =  indexElem + (unsigned)ceil(NumberofElementsToBin/2);
		
		// Enforce posive elements
		if (indexElem < ceil(NumberofElementsToBin/2)) {
			firstElement = 0;
		}
		
		// Enforce last element not to be greater than nelements
		if (lastElement > SpectralElements->getnSpectralElements()) {
			lastElement = SpectralElements->getnSpectralElements();
		}
		
		for (unsigned i=0; i<NXPoints; i++) {	
			
			unsigned indexElemInsideBin=0;
			
			for(unsigned iElem=firstElement; iElem < lastElement; iElem++) {
				float xcenter = SpectralElements->getphotoCenterX(iElem);
				float XCenterOfSubPixel = xcenter + InstrumentProfile->getIPixXCoordinate(i);				
				unsigned xx = (unsigned)floor(XCenterOfSubPixel);				
				unsigned yy = (unsigned)round(SpectralElements->getphotoCenterY(iElem));
				
				SubPixElementCounts[indexElemInsideBin++] = (float)flatMatrix[yy][xx];
			} 
			SubPixElemSampleMedianCounts[i] = operaArrayMedianQuick(indexElemInsideBin,SubPixElementCounts);		
		}
		free(SubPixElementCounts);
		
		/***** Prep for interpolation **********/
		for (unsigned i=0; i<NXPoints; i++) {
			SubPixXcoords[i] = InstrumentProfile->getIPixXCoordinate(i);
		}
		
		float yp1 = (SubPixElemSampleMedianCounts[1] - SubPixElemSampleMedianCounts[0])/(SubPixXcoords[1] - SubPixXcoords[0]);
		float ypn = (SubPixElemSampleMedianCounts[NXPoints-1] - SubPixElemSampleMedianCounts[NXPoints-2])/(SubPixXcoords[NXPoints-1] - SubPixXcoords[NXPoints-2]);
		
		// Call cubicspline to get second derivatives 
		cubicspline(SubPixXcoords, SubPixElemSampleMedianCounts, NXPoints, yp1, ypn, y2);
		/**** End of prep for interpolation ***/
		
		unsigned xlocalmin = (unsigned)floor(SpectralElements->getphotoCenterX(indexElem) - Geometry->getapertureWidth()/2);
		unsigned xlocalmax = (unsigned)ceil(SpectralElements->getphotoCenterX(indexElem) + Geometry->getapertureWidth()/2);
		
		if (xlocalmin < 0) { xlocalmin = 0;}		
		if (xlocalmax > nx) { xlocalmax = nx;}			
		
		unsigned yy = (unsigned)round(SpectralElements->getphotoCenterY(indexElem));
		
		for(unsigned xx=xlocalmin; xx<xlocalmax; xx++) {
			float ExpectedFlux;
			float xcoordInIPPixunits = (float)xx + 0.5 - SpectralElements->getphotoCenterX(indexElem);
			// Call splineinterpolate for interpolations 
			splineinterpolate(SubPixXcoords, SubPixElemSampleMedianCounts, y2, NXPoints, xcoordInIPPixunits, &ExpectedFlux);			
			outputMatrix[yy][xx] = flatMatrix[yy][xx]/ExpectedFlux;
		}
	}
	// DT Jan 2013 moved outside of loop
	free(y2);
	free(SubPixXcoords);
	free(SubPixElemSampleMedianCounts);
	free(SubPixElementCounts);
}

void operaSpectralOrder::calculateSNR(void) {
    
	unsigned length = SpectralElements->getnSpectralElements();
	SpectralElements->calculateFluxSNR();
	SNR = getCentralSmoothedSNR(length/DEFAULT_SNR_SMOOTHING);
	//SNR = SpectralElements->getFluxSNR(length/2);
	hasSNR = true;
}

float operaSpectralOrder::getCentralSmoothedSNR(int upperlowerbound) const {
    
	float snrs[MAXREFWAVELENGTHSPERORDER];
	unsigned length = SpectralElements->getnSpectralElements();
    
	if (length == 0) {
		return 1.0;
	}
	if (upperlowerbound == 0) {
		return SpectralElements->getFluxSNR(length/2);
	}
	int start = max(length/2-upperlowerbound, 0);
	int stop = min(length/2+upperlowerbound, length-1);
	int range = stop-start+1;
	if ((range > MAXREFWAVELENGTHSPERORDER) || (range == 0)) {
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned i = 0;
	for (int indexElem=start; indexElem <= stop; indexElem++){
		snrs[i++]= SpectralElements->getFluxSNR(indexElem);
	}
	return operaArrayMedianQuick(range, snrs);
}

float operaSpectralOrder::getPeakSmoothedSNR(int upperlowerbound) const {
    
	float snrs[MAXREFWAVELENGTHSPERORDER];
	unsigned length = SpectralElements->getnSpectralElements();
	unsigned maxindex = length/2;
    float maxSNR = 0.0;
	
	for (unsigned indexElem=0; indexElem < length; indexElem++){
		float snr = SpectralElements->getFluxSNR(indexElem);
		if (snr > maxSNR) {
			maxindex = indexElem;
			maxSNR = snr;
		}
	}
	if (length == 0) {
		return 1.0;
	}
	if (upperlowerbound == 0) {
		return SpectralElements->getFluxSNR(length/2);
	}
	int start = max(maxindex-upperlowerbound, 0);
	int stop = min(maxindex+upperlowerbound, length-1);
	int range = stop-start+1;
	if ((range > MAXREFWAVELENGTHSPERORDER) || (range == 0)) {
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned i = 0;
	for (int indexElem=start; indexElem <= stop; indexElem++){
		snrs[i++]= SpectralElements->getFluxSNR(indexElem);
	}
	return operaArrayMedianQuick(range, snrs);
}

void operaSpectralOrder::extractRawSum(operaFITSImage &inputImage, ofstream &sout){
	
	float distanceInPixelUnits;
	
	for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
		unsigned xlocalmin = (unsigned)floor(SpectralElements->getphotoCenterX(indexElem) - getGeometry()->getapertureWidth()/2);
		unsigned xlocalmax = (unsigned)ceil(SpectralElements->getphotoCenterX(indexElem) + getGeometry()->getapertureWidth()/2);
		unsigned ylocalmin = (unsigned)round(SpectralElements->getphotoCenterY(indexElem) - SpectralElements->getelementHeight()/2);
		unsigned ylocalmax = (unsigned)ceil(SpectralElements->getphotoCenterY(indexElem) + SpectralElements->getelementHeight()/2);
		
		if (xlocalmin < 0) { xlocalmin = 0;}		
		if (xlocalmax > inputImage.getnaxis1()) { xlocalmax = inputImage.getnaxis1();}		
		if (ylocalmin < 0) { ylocalmin = 0;}
		if (ylocalmax > inputImage.getnaxis2()) { ylocalmax = inputImage.getnaxis2();}
		
#ifdef PRINT_DEBUG
		cout << "Elem # " << indexElem << 
		"; PhotoCenter: [" << SpectralElements->getphotoCenterX(indexElem) << ":" << SpectralElements->getphotoCenterY(indexElem) << "]" <<
		"\tXRANGE= [" << xlocalmin <<
		":" << xlocalmax - 1 << "]" <<
		"\tYRANGE= [" << ylocalmin <<
		":" << ylocalmax - 1 << "]" << endl;
#endif				
		
		float sumFlux = 0.0;
		for(unsigned y=ylocalmin;y<ylocalmax;y++) {
			for(unsigned x=xlocalmin;x<xlocalmax;x++) {
				sumFlux += inputImage[y][x];
			}
		}
		
		SpectralElements->setFlux(sumFlux,indexElem);
		distanceInPixelUnits = SpectralElements->getdistd(indexElem);	
		
#ifdef PRINT_DEBUG
		cout << orderNumber <<
		"\t" << indexElem << 
		"\t" << SpectralElements->getphotoCenterX(indexElem) << 
		"\t" << SpectralElements->getphotoCenterY(indexElem) <<
		"\t" << xlocalmin <<
		"\t" << xlocalmax + 1 <<
		"\t" << ylocalmin <<
		"\t" << ylocalmax + 1 <<				
		"\t" << distanceInPixelUnits <<		
		"\t" << SpectralElements->getFlux(indexElem)<< endl;
#endif				
		
		// DT made to be more like upena...
		sout << distanceInPixelUnits <<
		"\t" << SpectralElements->getFlux(indexElem) <<
		"\t" << getorder() << endl;				
	}
	setSpectrumType(RawSpectrum);
	SpectralElements->setSpectrumType(RawSpectrum);
	sethasSpectralElements(true);
}

void operaSpectralOrder::extractRawSum(operaFITSImage &inputImage, float noise, float gain){
	
	for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
		unsigned xlocalmin = (unsigned)floor(SpectralElements->getphotoCenterX(indexElem) - getGeometry()->getapertureWidth()/2);
		unsigned xlocalmax = (unsigned)ceil(SpectralElements->getphotoCenterX(indexElem) + getGeometry()->getapertureWidth()/2);
		unsigned ylocalmin = (unsigned)round(SpectralElements->getphotoCenterY(indexElem) - SpectralElements->getelementHeight()/2);
		unsigned ylocalmax = (unsigned)ceil(SpectralElements->getphotoCenterY(indexElem) + SpectralElements->getelementHeight()/2);
		
		if (xlocalmin < 0) { xlocalmin = 0;}		
		if (xlocalmax > inputImage.getnaxis1()) { xlocalmax = inputImage.getnaxis1();}		
		if (ylocalmin < 0) { ylocalmin = 0;}
		if (ylocalmax > inputImage.getnaxis2()) { ylocalmax = inputImage.getnaxis2();}
		
#ifdef PRINT_DEBUG
		cout << "Elem # " << indexElem << 
		"; PhotoCenter: [" << SpectralElements->getphotoCenterX(indexElem) << ":" << SpectralElements->getphotoCenterY(indexElem) << "]" <<
		"\tXRANGE= [" << xlocalmin <<
		":" << xlocalmax - 1 << "]" <<
		"\tYRANGE= [" << ylocalmin <<
		":" << ylocalmax - 1 << "]" << endl;
#endif				
		
		float sumFlux = 0.0;
		float sumVar = 0.0;
		for(unsigned y=ylocalmin;y<ylocalmax;y++) {
			for(unsigned x=xlocalmin;x<xlocalmax;x++) {
				float flux = inputImage[y][x];
				sumFlux += flux*gain;
				sumVar += ((noise*noise) + fabs(flux*gain)); // detector noise + photon noise ;
			}
		}
		
		SpectralElements->setFlux(sumFlux,indexElem);
		SpectralElements->setFluxVariance(sumVar,indexElem);
		
#ifdef PRINT_DEBUG
		float distanceInPixelUnits = SpectralElements->getdistd(indexElem);	
		cout << orderNumber <<
		"\t" << indexElem << 
		"\t" << SpectralElements->getphotoCenterX(indexElem) << 
		"\t" << SpectralElements->getphotoCenterY(indexElem) <<
		"\t" << xlocalmin <<
		"\t" << xlocalmax + 1 <<
		"\t" << ylocalmin <<
		"\t" << ylocalmax + 1 <<				
		"\t" << distanceInPixelUnits <<		
		"\t" << SpectralElements->getFlux(indexElem)<< endl;
#endif				
	}
	setSpectrumType(RawSpectrum);
	SpectralElements->setSpectrumType(RawSpectrum);
	sethasSpectralElements(true);
}

void operaSpectralOrder::CalculateWavelengthSolution(void) {
	if (gethasGeometry() && gethasWavelength()) {
		operaWavelength *wavelength = getWavelength();
		operaGeometry *geometry = getGeometry();
		operaSpectralElements *spectralElements = getSpectralElements();
        double dmin = 0.0;
		wavelength->setDmin(dmin);
        double dmax = (double)geometry->CalculateDistance(geometry->getYmin(),geometry->getYmax());
        wavelength->setDmax(dmax);
		for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
			spectralElements->setwavelength(wavelength->evaluateWavelength(spectralElements->getdistd(k)), k);
		}
		spectralElements->setHasWavelength(true);
		sethasWavelength(true);
	}
}

void operaSpectralOrder::measureInstrumentProfileAlongRows(operaFITSImage &masterFlatImage, unsigned binsize, unsigned sampleElementForPlot, ostream *pout){
	
	if (InstrumentProfile->getysize()!=1 || InstrumentProfile->getYsampling()!=1){
		throw operaException("operaSpectralOrder: instrument profile ysize and ysampling must be 1.",operaErrorInstrumentProfileImproperDimensions, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	unsigned NXPoints = InstrumentProfile->getNXPoints();		
	unsigned NYPoints = InstrumentProfile->getNYPoints();
	
	if (NXPoints == 0){
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElemSampleMedianCounts = new float[NXPoints];
	float *SubPixElemSampleMedSigCounts = new float[NXPoints];	
	unsigned NumberOfElementSamples;
	unsigned NumberofElementsToBin = binsize;
	
	NumberOfElementSamples = (unsigned)ceil((float)SpectralElements->getnSpectralElements()/(float)NumberofElementsToBin); 
	
	if (NumberOfElementSamples == 0){
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElementCounts = new float[NumberOfElementSamples];		
	float *SumSubPixElements = new float[NumberOfElementSamples];
	
	for(unsigned k=0;k<NumberOfElementSamples;k++){
		
		unsigned firstElement = NumberofElementsToBin*(k);
		unsigned lastElement =  NumberofElementsToBin*(k+1);
		if (lastElement > SpectralElements->getnSpectralElements()){
            lastElement = SpectralElements->getnSpectralElements();   
        }
		
		SumSubPixElements[k] = 0;
		
        float distdmidElem = (float)fabs(SpectralElements->getdistd(lastElement-1) + SpectralElements->getdistd(firstElement))/2;		
        
        for (unsigned i=0; i<NXPoints; i++) {	
            
            unsigned indexElemInsideBin=0;
            
            for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
                
                float xcenter =  SpectralElements->getphotoCenterX(indexElem);
                float XCenterOfSubPixel = xcenter + InstrumentProfile->getIPixXCoordinate(i);				
                unsigned xx = (unsigned)floor(XCenterOfSubPixel);				
                
                float ycenter =  SpectralElements->getphotoCenterY(indexElem);
                float YCenterOfSubPixel = ycenter + InstrumentProfile->getIPixYCoordinate(0);                
                unsigned yy = (unsigned)floor(YCenterOfSubPixel);
                
                SubPixElementCounts[indexElemInsideBin++] = (float)masterFlatImage[yy][xx]/((float)InstrumentProfile->getXsampling()*(float)InstrumentProfile->getYsampling());
            } 
            
            SubPixElemSampleMedianCounts[i] = operaArrayMedian(indexElemInsideBin,SubPixElementCounts);
            SumSubPixElements[k] += SubPixElemSampleMedianCounts[i];		
            SubPixElemSampleMedSigCounts[i] = operaArrayMedianSigma(indexElemInsideBin, SubPixElementCounts, SubPixElemSampleMedianCounts[i]);
        }
        
        for (unsigned j=0; j<NYPoints; j++) {			
            for (unsigned i=0; i<NXPoints; i++) {
                InstrumentProfile->setdataCubeValues(SubPixElemSampleMedianCounts[i]/SumSubPixElements[k],i,j,k);
            }
        }
        InstrumentProfile->setdistd(distdmidElem,k);
        
#ifdef PRINT_DEBUG        
        unsigned midElement = firstElement + fabs(lastElement - firstElement)/2;        
        unsigned MinimumBinSize = Geometry->CalculateMinimumYBinSize(SpectralElements->getphotoCenterY(midElement));
        cout << k << " " 
        << SpectralElements->getphotoCenterY(midElement) << " " 
        << binsize << " " 
        << MinimumBinSize
        << endl;        
#endif        
        
	}
	
	InstrumentProfile->normalizeCubeData();
	
	InstrumentProfile->FitPolyMatrixtoIPDataVector(3,false);
	
	if (pout) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
    
	delete[] SubPixElementCounts;	
	delete[] SubPixElemSampleMedianCounts;	
	delete[] SubPixElemSampleMedSigCounts;	
	delete[] SumSubPixElements;	
}

void operaSpectralOrder::measureInstrumentProfileAlongRowsInto2DWithGaussian(operaFITSImage &masterFlatImage, operaFITSImage &badpix, unsigned binsize, float gaussSig, float tiltInDegrees, bool witherrors, unsigned sampleElementForPlot, ostream *pout, unsigned minimumLines){
	
	int npar = 3;
	double par[3];
	par[0] = 1.0/((double)gaussSig*sqrt(2.0*M_PI));
	par[1] = 0;
	par[2] = (double)gaussSig;    
	
    float tangentOfTiltInDegrees = tan(tiltInDegrees*M_PI/180.0);
    
	unsigned NXPoints = InstrumentProfile->getNXPoints();		
	unsigned NYPoints = InstrumentProfile->getNYPoints();
	
	if (NXPoints == 0){
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElemSampleMedianCounts = new float[NXPoints];
    float *SubPixElemSampleMedSigCounts = new float[NXPoints];
	
	unsigned NumberOfElementSamples;
	unsigned NumberofElementsToBin = binsize;
	
	NumberOfElementSamples = (unsigned)ceil((float)SpectralElements->getnSpectralElements()/(float)NumberofElementsToBin); 
    
	if (NumberofElementsToBin == 0){
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElementCounts = new float[NumberofElementsToBin];
	
    unsigned ndataPoints = 0;
    
	for(unsigned k=0;k<NumberOfElementSamples;k++){
		
		unsigned firstElement = NumberofElementsToBin*(k);
		unsigned lastElement =  NumberofElementsToBin*(k+1);
		if (lastElement > SpectralElements->getnSpectralElements()){
            lastElement = SpectralElements->getnSpectralElements();   
        }
		
        float IPNormalizationFactor = 0;
        float fluxFractionLost = 0; 
        
		for (unsigned i=0; i<NXPoints; i++) {
            
			unsigned indexElemInsideBin=0;
			
			for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
				
                float xcenter =  SpectralElements->getphotoCenterX(indexElem);
				float XCenterOfSubPixel = xcenter + InstrumentProfile->getIPixXCoordinate(i);
                unsigned xx = (unsigned)floor(XCenterOfSubPixel);				
                
                unsigned nps = 0;
                float SubPixMeanCounts = 0;
                
                if (xx > 0 && xx < masterFlatImage.getnaxis1()) {
                    float ycenter =  SpectralElements->getphotoCenterY(indexElem);
                    unsigned ylocalmin = (unsigned)round(ycenter -  SpectralElements->getelementHeight()/2);
                    unsigned ylocalmax = (unsigned)ceil(ycenter +  SpectralElements->getelementHeight()/2);
                    
                    for(unsigned yy=ylocalmin; yy<ylocalmax; yy++){
                        if (yy > 0 && yy < masterFlatImage.getnaxis2() &&
                            masterFlatImage[yy][xx] < SATURATIONLIMIT && 
                            badpix[yy][xx] == 1 ) { 
                            
                            SubPixMeanCounts += (float)masterFlatImage[yy][xx]; // in units of ADU					
                            nps++;
                        }
                    }
                }
                if(nps) {
                    SubPixElementCounts[indexElemInsideBin++] = SubPixMeanCounts/(float)(nps*InstrumentProfile->getXsampling());
                }                
			} 
            
            if(indexElemInsideBin > 2) {
                SubPixElemSampleMedianCounts[i] = operaArrayMedian(indexElemInsideBin,SubPixElementCounts);
                IPNormalizationFactor += SubPixElemSampleMedianCounts[i];		                
                if(witherrors) {
                    SubPixElemSampleMedSigCounts[i] = operaArrayMedianSigma(indexElemInsideBin, SubPixElementCounts, SubPixElemSampleMedianCounts[i]);
                }
            } else if (indexElemInsideBin > 0 && indexElemInsideBin < 3) {
                SubPixElemSampleMedianCounts[i] = operaArrayMean(indexElemInsideBin,SubPixElementCounts);
                IPNormalizationFactor += SubPixElemSampleMedianCounts[i];                
                if(witherrors) {
                    SubPixElemSampleMedSigCounts[i] = operaArraySigma(indexElemInsideBin,SubPixElementCounts);
                }
            } else {
                SubPixElemSampleMedianCounts[i] = NAN; 
                fluxFractionLost += 1.0/(float)(InstrumentProfile->getXsampling()*NXPoints);                
                if(witherrors) {
                    SubPixElemSampleMedSigCounts[i] = NAN;
                }
            }
		}
        
        if(IPNormalizationFactor) {
            float distdmidElem = (float)fabs(SpectralElements->getdistd(lastElement-1) + SpectralElements->getdistd(firstElement))/2;		
			
            for (unsigned j=0; j<NYPoints; j++) {	
                double IPPixYCoord = (double)(InstrumentProfile->getIPixYCoordinate(j));				
                for (unsigned i=0; i<NXPoints; i++) {
                    float ipvalue = 0.0;
                    if(!isnan(SubPixElemSampleMedianCounts[i])) { 
                        ipvalue = (SubPixElemSampleMedianCounts[i]/IPNormalizationFactor)*(1+fluxFractionLost/IPNormalizationFactor);
                    } else {
                        ipvalue = 1.0/(float)(InstrumentProfile->getXsampling()*NXPoints);
                    }
                    par[1]  = tangentOfTiltInDegrees*(double)(InstrumentProfile->getIPixXCoordinate(i));
                    InstrumentProfile->setdataCubeValues(ipvalue*(float)GaussianFunction(IPPixYCoord,par,npar),i,j,ndataPoints);
                }
            }        
            
            InstrumentProfile->setdistd(distdmidElem,ndataPoints);
            ndataPoints++;
        }
	}
    
    InstrumentProfile->setnDataPoints(ndataPoints);
    
    if(ndataPoints > minimumLines) {
        InstrumentProfile->normalizeCubeData();        
        InstrumentProfile->FitPolyMatrixtoIPDataVector(3,witherrors);
        
        if (pout != NULL) {
            InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
        }
        sethasInstrumentProfile(true);
    } else {
        sethasInstrumentProfile(false);
    }
    
    delete[] SubPixElementCounts;
	delete[] SubPixElemSampleMedianCounts;
    delete[] SubPixElemSampleMedSigCounts;
}

void operaSpectralOrder::setSpectralElementsByHeight(double Height) {
	
    if (!gethasGeometry()){
        throw operaException("operaSpectralOrder::setSpectralElementsByHeight: ",operaErrorHasNoGeometry, __FILE__, __FUNCTION__, __LINE__);
	}
    if(Height <= 0.0) {
        throw operaException("operaSpectralOrder:",operaErrorDivideByZeroError, __FILE__, __FUNCTION__, __LINE__);
    } 
	unsigned nElements = (unsigned)(ceil(fabs(Geometry->getYmax() - Geometry->getYmin())/Height));
	
    if(nElements == 0) {
        throw operaException("operaSpectralOrder:",operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
    } 
	if (SpectralElements == NULL) {
		SpectralElements = new operaSpectralElements(nElements, getSpectrumType());	
	} else {
		SpectralElements->resize(nElements);
		SpectralElements->setSpectrumType(getSpectrumType());
	}
	
	SpectralElements->setelementHeight(Height);
	
	double ymin = Geometry->getYmin();
	
	double ycenter = ymin + Height/2.0;
	
	for(unsigned indexElem=0; indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
        double d = Geometry->CalculateDistance(ymin, ycenter);
		SpectralElements->setphotoCenter(Geometry->getCenterPolynomial()->Evaluate(ycenter), ycenter, indexElem);
		SpectralElements->setdistd(d, indexElem);
		ycenter += Height;		
	}
	Geometry->CalculateAndSetOrderLength();
	SpectralElements->setHasDistance(true);
	sethasSpectralElements(true);
}

void operaSpectralOrder::setSpectralElementsByStitchingApertures(double effectiveApertureFraction) {
    
    if (!gethasGeometry()){
        throw operaException("operaSpectralOrder::setSpectralElementsByHeight: ",operaErrorHasNoGeometry, __FILE__, __FUNCTION__, __LINE__);
	}
    
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::setSpectralElementsByHeight: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double elementHeight = ExtractionApertures[0]->getApertureShape()->getLineYWidth()*effectiveApertureFraction;
    
    if (SpectralElements == NULL) {
		SpectralElements = new operaSpectralElements(0, getSpectrumType());
	} else {
		SpectralElements->resize(0);
		SpectralElements->setSpectrumType(getSpectrumType());
	}
    
    SpectralElements->setelementHeight(elementHeight);
	
    // First remember that the order polynomial is defined as:
    // x(y) = aa + bb*y + cc*y*y
    double aa = (double)Geometry->getCenterPolynomial()->getCoefficient(0);
    double bb = (double)Geometry->getCenterPolynomial()->getCoefficient(1);
    double cc = (double)Geometry->getCenterPolynomial()->getCoefficient(2);
    
    double ymin = (double)Geometry->getYmin();
    double ymax = (double)Geometry->getYmax();
	
    // Pick y coordinate of first element:
    double firstYCenter = ymin + elementHeight/2.0;
    double firstXCenter = (double)Geometry->getCenterPolynomial()->Evaluate((double)firstYCenter);
    // Then calculate the inital aperture line in the image reference system:
    
    // And the aperture central line is defined as: y(x) = dd*x + ff
    double dd = (double)ExtractionApertures[0]->getApertureShape()->getSlope();
    //double ff = (double)ExtractionApertures[0]->getApertureShape()->getIntercept();
    
	//    Geometry->getCenterPolynomial()->printEquation(&cout);
    
    double xcenter = firstXCenter;
    double ycenter = firstYCenter;
    
    SpectralElements->resize(1);
    SpectralElements->setphotoCenter(xcenter, ycenter, 0);
    SpectralElements->setdistd(Geometry->CalculateDistance(ymin, ycenter), 0);
	
    for(unsigned indexElem=1; ycenter + elementHeight/2 < ymax; indexElem++) {
		SpectralElements->resize(indexElem+1);
        
        // Now let's bring the aperture line to the image reference system:
        // y(x) - (firstYCenter + indexElem*elementHeight) = dd*(x - firstXcenter) + ff
        // And let's invert this equation and rearrange the terms:
        // x(y) = ([y - (firstYCenter + indexElem*elementHeight + ff)] / dd) + (firstXcenter)
        // Then we have the new line given by: x(y) = (1/dd)*y + firstXcenter - (firstYCenter + indexElem*elementHeight + ff)/dd
		// And intersect = firstXcenter - (firstYCenter + (double)indexElem*elementHeight + ff)/dd;
        
        double slope = (1.0/dd);
        double intersect = xcenter - (ycenter + elementHeight)/dd; // Don't acutally include ff, or we'd shift by that amount each iteration... 
		
#ifdef PRINT_DEBUG
		//      Below we print the line equation
		//        cout << "g" << indexElem << "(x)="<< slope << "*x+" << intersect << endl << endl;
        
        // Now we need to find a point {xcenter,ycenter} for the intersection
        // between the order polynomial and the aperture central line
		
        // and let's solve for the y-coordinate of the crossing point:
        // y*slope + intersect = aa + bb*y + cc*y*y
        // rearranging the terms:
        // alpha + beta*y + cc*y*y = 0, where alpha = (aa-intersect), beta = (bb-slope)
#endif
        double alpha = aa - intersect;
        double beta = bb - slope;
		
        // Applying Bhaskara's equation we have:
        // y = (-beta +/- sqrt(beta^2 - 4*alpha*cc)/(2*alpha)
		
        double ycenter_plus, ycenter_minus;
        if((-beta + sqrt(beta*beta - 4*alpha*cc)) != 0 && (beta*beta - 4*alpha*cc) >= 0) {
            ycenter_plus = (2*alpha)/(-beta + sqrt(beta*beta - 4*alpha*cc));
            ycenter_minus = (2*alpha)/(-beta - sqrt(beta*beta - 4*alpha*cc));
        } else {
            //throw operaException("operaSpectralOrder:",operaErrorDivideByZeroError, __FILE__, __FUNCTION__, __LINE__);
			// DT May 20 2014 -- We are likely to run in to problems when indexing through MAXSPECTRALELEMENTSPERORDER
			// so, just break rather than throwing an exception and thus aborting...
			break;
            //throw operaException("operaSpectralOrder:",operaErrorDivideByZeroError, __FILE__, __FUNCTION__, __LINE__);
		}
		//        cout << "ycenter_plus=" << ycenter_plus << " ycenter_minus=" << ycenter_minus << endl;
        
        // pick the branch with closest y to the previous point:
        double dist_yplus = fabs(ycenter_plus - ycenter);
        double dist_yminus = fabs(ycenter_minus - ycenter);
		//        cout << "dist_yplus=" << dist_yplus << " dist_yminus=" << dist_yminus << endl;
		
        if(dist_yminus < dist_yplus) {
            ycenter = ycenter_minus;
        } else {
            ycenter = ycenter_plus;
        }
        
        //ycenter = (2*alpha)/(-beta - sqrt(beta*beta - 4*alpha*cc));
        xcenter = (double)Geometry->getCenterPolynomial()->Evaluate((double)ycenter);
        
        SpectralElements->setphotoCenter(xcenter, ycenter, indexElem);
        SpectralElements->setdistd((double)Geometry->CalculateDistance(ymin, ycenter), indexElem);
    }
    Geometry->CalculateAndSetOrderLength();
	SpectralElements->setHasDistance(true);
	sethasSpectralElements(true);
}

void operaSpectralOrder::setApertureElements(operaSpectralOrder_t format) {
    if(!hasSpectralElements) {
        throw operaException("operaSpectralOrder:",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);	
    } 
    
	// only reallocate if we need to...
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
		if (BeamElements[beam] == NULL) {
			BeamElements[beam] = new operaSpectralElements(SpectralElements->getnSpectralElements(), format);    			
		} else {
			BeamElements[beam]->resize(SpectralElements->getnSpectralElements());
		}
    }
    for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
		if (BackgroundElements[backgroundIndex] == NULL) {
			BackgroundElements[backgroundIndex] = new operaSpectralElements(SpectralElements->getnSpectralElements(), format);    
		} else {
			BackgroundElements[backgroundIndex]->resize(SpectralElements->getnSpectralElements());
		}
    }
	
	for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {     
		
        double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
        double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            double beamXcenter = elemXcenter + (double)ExtractionApertures[beam]->getApertureShape()->getMidPoint().getXcoord();
            double beamYcenter = elemYcenter + (double)ExtractionApertures[beam]->getApertureShape()->getMidPoint().getYcoord();
			
            BeamElements[beam]->setphotoCenter(beamXcenter,beamYcenter,indexElem);    
        }
        for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
            double backgroundXcenter = elemXcenter + (double)BackgroundApertures[backgroundIndex]->getApertureShape()->getMidPoint().getXcoord();
            double backgroundYcenter = elemYcenter + (double)BackgroundApertures[backgroundIndex]->getApertureShape()->getMidPoint().getYcoord();
            
            BackgroundElements[backgroundIndex]->setphotoCenter(backgroundXcenter,backgroundYcenter,indexElem);                
        }           
    }    
}

operaExtractionAperture<Line> *operaSpectralOrder::calculateMainApertureFromExtractionBeams(bool useIP) {
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::calculateMainApertureFromExtractionBeams: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaPoint point(0,0);
    float mainSlope = ExtractionApertures[0]->getApertureShape()->getSlope();
    float mainWidth = ExtractionApertures[0]->getApertureShape()->getWidth();
    float mainLength = 0;
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        mainLength += ExtractionApertures[beam]->getApertureShape()->getLength();
    }
	
    Line MainExtractionLineAperture(mainSlope, mainWidth, mainLength, point);
    
    operaExtractionAperture<Line> *MainExtractionAperture;
    
    if(useIP) {
        if(!gethasInstrumentProfile()) {
            throw operaException("operaSpectralOrder::calculateMainApertureFromExtractionBeams: ",operaErrorHasNoInstrumentProfile, __FILE__, __FUNCTION__, __LINE__);
        }
		MainExtractionAperture = new operaExtractionAperture<Line>(&MainExtractionLineAperture,InstrumentProfile);
    } else {
        unsigned xsampling = ExtractionApertures[0]->getXsampling();
        unsigned ysampling = ExtractionApertures[0]->getYsampling();
        MainExtractionAperture = new operaExtractionAperture<Line>(&MainExtractionLineAperture, xsampling, ysampling);
    }
    
    return MainExtractionAperture;
}


void operaSpectralOrder::deleteInstrumentProfile(void) {
	
	if (InstrumentProfile) {
		delete InstrumentProfile;
	}
	InstrumentProfile = NULL;
	sethasInstrumentProfile(false);
}

//
// Caution -- replaces the instrumentProfile pointer!!!!
//
void operaSpectralOrder::setInstrumentProfileVector(unsigned IPxsize, unsigned IPxsampling, unsigned IPysize, unsigned IPysampling, unsigned NDataPoints) {
	
	if (InstrumentProfile) {
		delete InstrumentProfile;
	}
	InstrumentProfile = new operaInstrumentProfile(IPxsize,IPxsampling,IPysize,IPysampling, NDataPoints);
	sethasInstrumentProfile(true);
}

/*
 * interface used by read/write orders
 */
void operaSpectralOrder::setSpectralLines(operaSpectralLines *spectralLines) {
	SpectralLines = spectralLines;
	sethasSpectralLines(true);
}

void operaSpectralOrder::setSpectralLines(operaFITSImage &masterCompImage, operaFITSImage &badpix, operaFITSImage &bias, float noise, float gain, float ReferenceLineWidth,float DetectionThreshold, float LocalMaxFilterWidth, float MinPeakDepth) {
	
    if(!gethasSpectralElements()) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    if(!gethasInstrumentProfile()) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoInstrumentProfile, __FILE__, __FUNCTION__, __LINE__);
    }    
    if(!SpectralElements->getHasXCorrelation()) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoXCorrelation, __FILE__, __FUNCTION__, __LINE__);
    }
	
    if(!SpectralElements->getHasRawSpectrum()) { 
        
        unsigned NXPoints = InstrumentProfile->getNXPoints();		
        unsigned NYPoints = InstrumentProfile->getNYPoints();
        unsigned nMaxDataPoints = SpectralElements->getnSpectralElements();
        
        for(unsigned indexElem=0;indexElem < nMaxDataPoints; indexElem++) {    
            
            float xcenter =  SpectralElements->getphotoCenterX(indexElem);
            float ycenter =  SpectralElements->getphotoCenterY(indexElem);
            float distdElem = SpectralElements->getdistd(indexElem);
            
            float my = 0;
            float myvar = 0; 
            float fluxFractionLost = 0;
            float IPNormalizationFactor = 0;  
            
            for (unsigned j=0; j<NYPoints; j++) {	
                float YCenterOfSubPixel = ycenter + InstrumentProfile->getIPixYCoordinate(j);				
                unsigned yy = (unsigned)floor(YCenterOfSubPixel);             
                for (unsigned i=0; i<NXPoints; i++) {
                    float XCenterOfSubPixel = xcenter + InstrumentProfile->getIPixXCoordinate(i);				
                    unsigned xx = (unsigned)floor(XCenterOfSubPixel); 
                    if (xx > 0 && xx < masterCompImage.getnaxis1() &&
                        yy > 0 && yy < masterCompImage.getnaxis2() &&
                        masterCompImage[yy][xx] < SATURATIONLIMIT && 
                        badpix[yy][xx] == 1 && 
                        (float)masterCompImage[yy][xx] > 0 ) {                    
                        my +=  (float)(masterCompImage[yy][xx]  - bias[yy][xx]) * InstrumentProfile->getipDataFromPolyModel(distdElem,i,j);
                        myvar += (noise/gain)*(noise/gain)/float(NXPoints*NYPoints) + fabs((float)masterCompImage[yy][xx] * InstrumentProfile->getipDataFromPolyModel(distdElem,i,j));
                    } else {
                        fluxFractionLost += InstrumentProfile->getipDataFromPolyModel(distdElem,i,j);
                    }
                    IPNormalizationFactor += InstrumentProfile->getipDataFromPolyModel(distdElem,i,j); //In case IP is not properly normalized. This could occurs after the polynomial fit
                }
            }      
            
            if(my && IPNormalizationFactor) {
                my = (my/IPNormalizationFactor)*(1+fluxFractionLost/IPNormalizationFactor);
                myvar = (myvar/IPNormalizationFactor)*(1+fluxFractionLost/IPNormalizationFactor);
            } else {
                my = NAN;
                myvar = NAN;
            }
            
            SpectralElements->setFlux(my,indexElem);
            SpectralElements->setFluxVariance(myvar,indexElem);
            
#ifdef PRINT_DEBUG	
            //print out extracted cross-correlation spectrum 
            cout << SpectralElements->getphotoCenterY(indexElem) <<
            "\t" << SpectralElements->getFlux(indexElem) <<
            "\t" << sqrt(SpectralElements->getFluxVariance(indexElem)) <<
            "\t" << getorder() << endl;			
#endif	        
        }     
        setSpectrumType(RawSpectrum);   
        SpectralElements->setSpectrumType(RawSpectrum);   
        SpectralElements->setHasRawSpectrum(true);   
    }
    if (SpectralLines) {
		delete SpectralLines;
	}
    SpectralLines = new operaSpectralLines(SpectralElements,(double)ReferenceLineWidth,y_distance_disp);
    
    SpectralLines->detectSpectralFeatures((double)DetectionThreshold,(double)LocalMaxFilterWidth,(double)MinPeakDepth); 
    
    if(SpectralLines->getNFeatures() == 0) { 
        sethasSpectralLines(false);
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    } else {
        sethasSpectralLines(true);	// Note that if the exception was thrown then hasspectrallines is false
    }
}

void operaSpectralOrder::calculateXCorrBetweenIPandImage(operaFITSImage &Image, operaFITSImage &badpix, ostream *pout) {
    if(!gethasInstrumentProfile()) {
        throw operaException("operaSpectralOrder::calculateXCorrBetweenIPandImage: ",operaErrorHasNoInstrumentProfile, __FILE__, __FUNCTION__, __LINE__);
    }
    if(!gethasSpectralElements()) {
        throw operaException("operaSpectralOrder::calculateXCorrBetweenIPandImage: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }    
    
    const unsigned NXPoints = InstrumentProfile->getNXPoints();		
    const unsigned NYPoints = InstrumentProfile->getNYPoints();
    const unsigned nMaxDataPoints = SpectralElements->getnSpectralElements();
    
    for(unsigned indexElem=0;indexElem < nMaxDataPoints; indexElem++) {    
        float xcenter =  SpectralElements->getphotoCenterX(indexElem);
        float ycenter =  SpectralElements->getphotoCenterY(indexElem);
        float distdElem = SpectralElements->getdistd(indexElem);
		
		VLArray<int> xCoords(NXPoints);
        VLArray<int> yCoords(NYPoints);
		for (unsigned j=0; j<NYPoints; j++) yCoords(j) = (int)(ycenter + InstrumentProfile->getIPixYCoordinate(j));
		for (unsigned i=0; i<NXPoints; i++) xCoords(i) = (int)(xcenter + InstrumentProfile->getIPixXCoordinate(i));
        
        unsigned iMin, jMin, iMax, jMax;
        for (iMin=0; iMin < NXPoints; iMin++) if(xCoords(iMin) > 0) break; //find lowest index where xCoords[i] > 0
        for (jMin=0; jMin < NYPoints; jMin++) if(yCoords(jMin) > 0) break; //find lowest index where yCoords[j] > 0
		for (iMax=NXPoints; iMax > 0; iMax--) if(unsigned(xCoords(iMax-1)) < Image.getnaxis1() || xCoords(iMax-1) < 0) break; //find highest index where xCoords[i] < naxis1
		for (jMax=NYPoints; jMax > 0; jMax--) if(unsigned(yCoords(jMax-1)) < Image.getnaxis2() || yCoords(jMax-1) < 0) break; //find highest index where yCoords[j] < naxis2
		
		VLArray<bool> validCoords(NYPoints, NXPoints);
		for (unsigned j=jMin; j<jMax; j++) {
            for (unsigned i=iMin; i<iMax; i++) {
				validCoords(j, i) = Image[yCoords(j)][xCoords(i)] > 0 && Image[yCoords(j)][xCoords(i)] < SATURATIONLIMIT && badpix[yCoords(j)][xCoords(i)];
			}
		}
        
        VLArray<float> ipvals(NYPoints, NXPoints);
		float avgImg = 0, meanIP = 0;
        unsigned npImg = 0;
		for (unsigned j=jMin; j<jMax; j++) {	
            for (unsigned i=iMin; i<iMax; i++) {
                if (validCoords(j, i)) {                    
                    avgImg += (float)Image[yCoords(j)][xCoords(i)];
                    ipvals(j, i) = InstrumentProfile->getipDataFromPolyModel(distdElem,i,j);
                    meanIP += ipvals(j, i);
                    npImg++;
                }
            }
        }
        avgImg /= (float)npImg;
        meanIP /= (float)npImg;
        
        float Xcorr = 0, imgsqr = 0, ipsqr = 0;
		for (unsigned j=jMin; j<jMax; j++) {	
            for (unsigned i=iMin; i<iMax; i++) {
                if (validCoords(j, i)) {
					const float ip = ipvals(j, i) - meanIP;
					const float imgval = (float)Image[yCoords(j)][xCoords(i)] - avgImg;
                    Xcorr +=  imgval * ip;
                    imgsqr += imgval * imgval;
                    ipsqr += ip * ip;
                }
            }
        }      
        
        if(Xcorr) Xcorr /= sqrt(imgsqr*ipsqr);
        else Xcorr = NAN;
        SpectralElements->setXCorrelation((double)Xcorr,indexElem);
        
        if (pout != NULL) {
            //print out cross-correlation 
            *pout << getorder() << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "            
            << SpectralElements->getXCorrelation(indexElem) << " "
            << endl;			
        }
    }
    SpectralElements->setHasXCorrelation(true);
}

void operaSpectralOrder::measureInstrumentProfile(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, unsigned minimumLines) {
	
    if(!gethasSpectralLines()) {
        throw operaException("operaSpectralOrder::measureInstrumentProfile: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double *LinePositionVector = new double[SpectralLines->getnLines()];
    double *LineSigmaVector = new double[SpectralLines->getnLines()];
    double *LineAmplitudeVector = new double[SpectralLines->getnLines()];
    
    unsigned nSelectedLines = SpectralLines->selectLines(MaxContamination,nSigCut,amplitudeCutOff,LinePositionVector,LineSigmaVector,LineAmplitudeVector);
    
    if(operaArgumentHandler::debug) {
		for(unsigned k=0; k<nSelectedLines; k++) {   
			cout << k << " " << LinePositionVector[k] << " " << LineAmplitudeVector[k] << " " << LineSigmaVector[k] << endl;
		}
	}
	
    unsigned NXPoints = InstrumentProfile->getNXPoints();		
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    unsigned cubeDataIndex = 0;
    
    for(unsigned k=0; k<nSelectedLines; k++) {
        
        float ycenter =  LinePositionVector[k];
        
		float xcenter = (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        
		float distd = (float)Geometry->CalculateDistance(Geometry->getYmin(), ycenter);        
        
        InstrumentProfile->setdataCubeValues(masterCompImage,badpix,xcenter,ycenter,cubeDataIndex);        
		
        InstrumentProfile->subtractOuterFrame(cubeDataIndex);
		
        float FractionOfFluxLost = 0;
        float SumOfUsefulFluxWithinLine = 0;        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                if(!isnan(InstrumentProfile->getdataCubeValues(i,j,cubeDataIndex))){
                    SumOfUsefulFluxWithinLine += InstrumentProfile->getdataCubeValues(i,j,cubeDataIndex);                    
                } else {
                    FractionOfFluxLost += InstrumentProfile->getipDataFromPolyModel(distd,i,j);
                }  
            }
        }   
        
        float NormalizationFactor = 0;
        
        if(FractionOfFluxLost<1 && SumOfUsefulFluxWithinLine)
            NormalizationFactor = SumOfUsefulFluxWithinLine/(1.0 - FractionOfFluxLost);
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                float ipmatrixvalue = 0;
                if(NormalizationFactor && !isnan(InstrumentProfile->getdataCubeValues(i,j,cubeDataIndex))){
                    ipmatrixvalue = InstrumentProfile->getdataCubeValues(i,j,cubeDataIndex)/NormalizationFactor;
                } else {
                    ipmatrixvalue = InstrumentProfile->getipDataFromPolyModel(distd,i,j);  
                }
                InstrumentProfile->setdataCubeValues(ipmatrixvalue,i,j,cubeDataIndex);
            }
        }  
        
        InstrumentProfile->normalizeCubeData(cubeDataIndex);        
		
        InstrumentProfile->setdistd(distd,cubeDataIndex);            
        
        cubeDataIndex++;  
		
    }
    
    InstrumentProfile->setnDataPoints(cubeDataIndex);  
    
    if(cubeDataIndex < minimumLines) {
#ifdef PRINT_DEBUG
		cerr << "operaSpectralOrder::measureInstrumentProfile: Order " << getorder() << " IP measurements not possible." << endl;
#endif
		sethasInstrumentProfile(false);
    } else if(cubeDataIndex >= minimumLines && cubeDataIndex <= 3*minimumLines) {
        InstrumentProfile->FitMediantoIPDataVector();
    } else if (cubeDataIndex > 3*minimumLines) {
        InstrumentProfile->FitPolyMatrixtoIPDataVector(3,false);        
    }    
	
	if (pout != NULL) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
    
    delete[] LinePositionVector;
    delete[] LineSigmaVector;
    delete[] LineAmplitudeVector;
}

void operaSpectralOrder::measureInstrumentProfileWithBinning(operaFITSImage &masterCompImage, operaFITSImage &badpix, double binsize, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, unsigned minimumLines) {
    
    if(!gethasSpectralLines()) {
        throw operaException("operaSpectralOrder::measureInstrumentProfileWithBinning: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double *LinePositionVector = new double[SpectralLines->getnLines()];
    double *LineSigmaVector = new double[SpectralLines->getnLines()];
    double *LineAmplitudeVector = new double[SpectralLines->getnLines()];
    
    unsigned nSelectedLines = SpectralLines->selectLines(MaxContamination,nSigCut,amplitudeCutOff,LinePositionVector,LineSigmaVector,LineAmplitudeVector);
    
    if(nSelectedLines == 0) {
		delete[] LinePositionVector;
		delete[] LineSigmaVector;
		delete[] LineAmplitudeVector;
        throw operaException("operaSpectralOrder::measureInstrumentProfileWithBinning: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
	if(operaArgumentHandler::debug) {
		for(unsigned k=0; k<nSelectedLines; k++) {   
			cout << k << " " << LinePositionVector[k] << " " << LineAmplitudeVector[k] << " " << LineSigmaVector[k] << endl;
		}
	}
    
    unsigned NXPoints = InstrumentProfile->getNXPoints();		
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    unsigned cubeDataIndex = 0;   
    
    unsigned ipxsize = InstrumentProfile->getxsize();	
    unsigned ipysize = InstrumentProfile->getysize();	
    unsigned ipxsampling = InstrumentProfile->getXsampling();  
    unsigned ipysampling = InstrumentProfile->getYsampling();  
    
    operaInstrumentProfile tempIP(ipxsize,ipxsampling,ipysize,ipysampling,nSelectedLines);
    float dlimit = binsize;
    
    unsigned npInBin = 0;
    float averageDistd = 0;
    
    for(unsigned k=0; k<nSelectedLines; k++) {
        
        float ycenter =  LinePositionVector[k];
		float xcenter = (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        
		float distd = Geometry->CalculateDistance(Geometry->getYmin(), ycenter);        
        averageDistd += distd;
        
        tempIP.setdistd(distd,npInBin);
        
        tempIP.setdataCubeValues(masterCompImage,badpix,xcenter,ycenter,npInBin);
        tempIP.subtractOuterFrame(npInBin);
        
        float FractionOfFluxLost = 0;
        float SumOfUsefulFluxWithinLine = 0;        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                if(!isnan(tempIP.getdataCubeValues(i,j,npInBin))){
                    SumOfUsefulFluxWithinLine += tempIP.getdataCubeValues(i,j,npInBin);                    
                } else {
                    FractionOfFluxLost += InstrumentProfile->getipDataFromPolyModel(distd,i,j);
                }  
            }
        }   
        
        float NormalizationFactor = 0;
        
        if(FractionOfFluxLost<1 && SumOfUsefulFluxWithinLine)
            NormalizationFactor = SumOfUsefulFluxWithinLine/(1.0 - FractionOfFluxLost);
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                float ipmatrixvalue;
                if(NormalizationFactor && !isnan(tempIP.getdataCubeValues(i,j,npInBin))){
                    ipmatrixvalue = tempIP.getdataCubeValues(i,j,npInBin)/NormalizationFactor;
                } else {
                    ipmatrixvalue = InstrumentProfile->getipDataFromPolyModel(distd,i,j);  
                }
                tempIP.setdataCubeValues(ipmatrixvalue,i,j,npInBin);
            }
        }  
        
        tempIP.normalizeCubeData(npInBin);        
        
        npInBin++;
        
        if((npInBin >= 3 && distd > dlimit) || k == nSelectedLines-1)  {
            
			tempIP.setnDataPoints(npInBin);            
            tempIP.FitMediantoIPDataVector();            
            averageDistd /= (float)npInBin;
            
            InstrumentProfile->setdistd(averageDistd,cubeDataIndex);
            
            for (unsigned j=0; j<NYPoints; j++) {
                for (unsigned i=0; i<NXPoints; i++) {
                    float ipmatrixvalue = tempIP.getipDataFromPolyModel(averageDistd,i,j);
                    InstrumentProfile->setdataCubeValues(ipmatrixvalue,i,j,cubeDataIndex);
                }
            }
            InstrumentProfile->normalizeCubeData(cubeDataIndex); 
			
            dlimit += binsize;
            averageDistd = 0;
            cubeDataIndex++;
            npInBin = 0;
        } else if (npInBin < 3 && distd > dlimit) {
            dlimit += binsize;  
        }      
    }
    
    InstrumentProfile->setnDataPoints(cubeDataIndex);        
    
    if(cubeDataIndex < minimumLines) {
#ifdef PRINT_DEBUG	    
		cerr << "operaSpectralOrder::measureInstrumentProfileWithBinning: Order " << getorder() << " IP measurements not possible." << endl;
#endif
        sethasInstrumentProfile(false);
    } else if(cubeDataIndex >= minimumLines) {
        InstrumentProfile->FitPolyMatrixtoIPDataVector(3,false);        
    }    
    
	if (pout != NULL) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
	
    delete[] LinePositionVector;
    delete[] LineSigmaVector;
    delete[] LineAmplitudeVector;
}

void operaSpectralOrder::measureInstrumentProfileUsingMedian(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, unsigned minimumLines) {
    
    if(!gethasSpectralLines()) {
        throw operaException("operaSpectralOrder::measureInstrumentProfileUsingMedian: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double *LinePositionVector = new double[SpectralLines->getnLines()];
    double *LineSigmaVector = new double[SpectralLines->getnLines()];
    double *LineAmplitudeVector = new double[SpectralLines->getnLines()];
    
    unsigned nSelectedLines = SpectralLines->selectLines(MaxContamination,nSigCut,amplitudeCutOff,LinePositionVector,LineSigmaVector,LineAmplitudeVector);
    
	if(operaArgumentHandler::debug) {
		if(nSelectedLines == 0) { 
			cout << "operaSpectralOrder::measureInstrumentProfileUsingMedian: Order " << getorder() << " no spectral lines." << endl;
		}
		for(unsigned k=0; k<nSelectedLines; k++) {   
			cout << k << " " << LinePositionVector[k] << " " << LineAmplitudeVector[k] << " " << LineSigmaVector[k] << endl;
		}
	}
    
    unsigned NXPoints = InstrumentProfile->getNXPoints();		
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    for(unsigned k=0; k<nSelectedLines; k++) {
        
        float ycenter =  LinePositionVector[k];
		float xcenter = (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        
		float distd = Geometry->CalculateDistance(Geometry->getYmin(), ycenter);        
        
        InstrumentProfile->setdistd(distd,k);
        InstrumentProfile->setdataCubeValues(masterCompImage,badpix,xcenter,ycenter,k);
        InstrumentProfile->subtractOuterFrame(k);
        
        float FractionOfFluxLost = 0;
        float SumOfUsefulFluxWithinLine = 0;  
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                if(!isnan(InstrumentProfile->getdataCubeValues(i,j,k))){
                    SumOfUsefulFluxWithinLine += InstrumentProfile->getdataCubeValues(i,j,k);                    
                } else {
                    FractionOfFluxLost += InstrumentProfile->getipDataFromPolyModel(distd,i,j);
                }  
            }
        }   
        
        float NormalizationFactor = 0;
        
        if(FractionOfFluxLost<1 && SumOfUsefulFluxWithinLine)
            NormalizationFactor = SumOfUsefulFluxWithinLine/(1.0 - FractionOfFluxLost);
		
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                float ipmatrixvalue = 0;
                if(NormalizationFactor && !isnan(InstrumentProfile->getdataCubeValues(i,j,k))){
                    ipmatrixvalue = InstrumentProfile->getdataCubeValues(i,j,k)/NormalizationFactor;
                } else {
                    ipmatrixvalue = InstrumentProfile->getipDataFromPolyModel(distd,i,j);  
                }
                InstrumentProfile->setdataCubeValues(ipmatrixvalue,i,j,k);
            }
        }  
        
        InstrumentProfile->normalizeCubeData(k);        
    }
	
    if(nSelectedLines > minimumLines) {
        InstrumentProfile->FitMediantoIPDataVector();
        sethasInstrumentProfile(true);
    } else {
        sethasInstrumentProfile(false);
    }
	if (pout != NULL) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
	
    delete[] LinePositionVector;
    delete[] LineSigmaVector;
    delete[] LineAmplitudeVector;
}

void operaSpectralOrder::measureInstrumentProfileUsingWeightedMean(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout) {
    
    if(!gethasSpectralLines()) {
        throw operaException("operaSpectralOrder::measureInstrumentProfileUsingWeightedMean: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
    if(SpectralLines->getnLines() == 0) {
        throw operaException("operaSpectralOrder::measureInstrumentProfileUsingWeightedMean: ",operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double *LinePositionVector = new double[SpectralLines->getnLines()];
    double *LineSigmaVector = new double[SpectralLines->getnLines()];
    double *LineAmplitudeVector = new double[SpectralLines->getnLines()];
    
    unsigned nSelectedLines = SpectralLines->selectLines(MaxContamination,nSigCut,amplitudeCutOff,LinePositionVector,LineSigmaVector,LineAmplitudeVector);
    
    double TotalWeight = 0;  
    
    for(unsigned k=0; k<nSelectedLines; k++) {  
        TotalWeight += LineAmplitudeVector[k];
    }  
    
    if(operaArgumentHandler::debug) {
		for(unsigned k=0; k<nSelectedLines; k++) {   
			cout << k << " " << LinePositionVector[k] << " " << LineAmplitudeVector[k] << " " << LineSigmaVector[k] << endl;
		}
	}
    
    unsigned NXPoints = InstrumentProfile->getNXPoints();
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    unsigned ipxsize = InstrumentProfile->getxsize();
    unsigned ipysize = InstrumentProfile->getysize();
    unsigned ipxsampling = InstrumentProfile->getXsampling();
    unsigned ipysampling = InstrumentProfile->getYsampling();
    
    operaInstrumentProfile tempIP(ipxsize,ipxsampling,ipysize,ipysampling,1);
    
    for (unsigned j=0; j<NYPoints; j++) {	 
        for (unsigned i=0; i<NXPoints; i++) {
            tempIP.setdataCubeValues(0.0,i,j,0);
        }
    }
    
    for(unsigned k=0; k<nSelectedLines; k++) {
        
        float ycenter =  LinePositionVector[k];
		float xcenter = (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        
		float distd = Geometry->CalculateDistance(Geometry->getYmin(), ycenter);        
        
        InstrumentProfile->setdistd(distd,k);
        InstrumentProfile->setdataCubeValues(masterCompImage,badpix,xcenter,ycenter,k);
        InstrumentProfile->subtractOuterFrame(k);
        
        float FractionOfFluxLost = 0;
        float SumOfUsefulFluxWithinLine = 0;  
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                if(!isnan(InstrumentProfile->getdataCubeValues(i,j,k))){
                    SumOfUsefulFluxWithinLine += InstrumentProfile->getdataCubeValues(i,j,k);                    
                } else {
                    FractionOfFluxLost += InstrumentProfile->getipDataFromPolyModel(distd,i,j);
                }  
            }
        }   
		
        float NormalizationFactor = 0;
        
        if(FractionOfFluxLost<1 && SumOfUsefulFluxWithinLine)
            NormalizationFactor = SumOfUsefulFluxWithinLine/(1.0 - FractionOfFluxLost);
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                float tempipmatrixvalue = tempIP.getdataCubeValues(i,j,0);
                if(NormalizationFactor && !isnan(InstrumentProfile->getdataCubeValues(i,j,k))){
                    tempipmatrixvalue += (InstrumentProfile->getdataCubeValues(i,j,k)/NormalizationFactor)*(float)(LineAmplitudeVector[k]/TotalWeight);
                } else {
                    tempipmatrixvalue += InstrumentProfile->getipDataFromPolyModel(distd,i,j)*(float)(LineAmplitudeVector[k]/TotalWeight);
                }
                tempIP.setdataCubeValues(tempipmatrixvalue,i,j,0);
            }
        }  
    }
    
    InstrumentProfile->setnDataPoints(1);
    
    for (unsigned j=0; j<NYPoints; j++) {
        for (unsigned i=0; i<NXPoints; i++) {
            InstrumentProfile->setdataCubeValues(tempIP.getdataCubeValues(i,j,0),i,j,0);
        }
    }
    
    unsigned zero = 0;
    InstrumentProfile->normalizeCubeData(zero);
    
    InstrumentProfile->FitPolyMatrixtoIPDataVector(1,false); 
    
	if (pout != NULL) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
    
    delete[] LinePositionVector;
    delete[] LineSigmaVector;
    delete[] LineAmplitudeVector;
}

void operaSpectralOrder::recenterOrderPosition(void) {
	
    for(unsigned i=0; i<Geometry->getNdatapoints(); i++) {    
        
        float ycenter =  Geometry->getCenterY(i);
        float xcenter =  (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        float d = Geometry->CalculateDistance(Geometry->getYmin(), ycenter);
        
        // Set error as the size of an IP sub-pixel. In fact the error can be calculated,
        // however geometry is not using errors to trace order position, so we leave like this for now.
        
        double xerror = 1.0/(double)InstrumentProfile->getXsampling(); 
		
        double x = (double)(xcenter + InstrumentProfile->getIPphotoCenterX(d));
        
        double y = (double)(ycenter + InstrumentProfile->getIPphotoCenterY(d));
        
        Geometry->resetCenter(x,y,xerror,i);
    }
    
    double chisqr;
    bool witherrors = false;
    unsigned maxOrderofTracePolynomial = 3;	// was 7, too high... DT Apr 25 2013
    
    Geometry->traceOrder(maxOrderofTracePolynomial, chisqr, witherrors);
	
}

/*
 * Print out beam spectra. Useful for plotting.
 */

void operaSpectralOrder::printBeamSpectrum(ostream *pout) {
    if (pout != NULL && gethasSpectralElements()) {
        unsigned NumberofElements = SpectralElements->getnSpectralElements();
        
        for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
            *pout << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getXCorrelation(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *pout << BeamElements[beam]->getphotoCenterX(indexElem) << " "
                << BeamElements[beam]->getphotoCenterY(indexElem) << " "
                << BeamElements[beam]->getFlux(indexElem) << " "
                << BeamElements[beam]->getFluxVariance(indexElem) << " ";
            }
            *pout << endl;
        }
        *pout << endl;
    }
}

void operaSpectralOrder::printBeamSpectrum(string addFirstColumnEntry, ostream *pout) {
    if (pout != NULL && gethasSpectralElements()) {
        unsigned NumberofElements = SpectralElements->getnSpectralElements();
        
        for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
            *pout << addFirstColumnEntry << " "
            << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getXCorrelation(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *pout << BeamElements[beam]->getphotoCenterX(indexElem) << " "
                << BeamElements[beam]->getphotoCenterY(indexElem) << " "
                << BeamElements[beam]->getFlux(indexElem) << " "
                << BeamElements[beam]->getFluxVariance(indexElem) << " ";
            }
            *pout << endl;
        }
        *pout << endl;
    }
}

// Extracts the flux per subpixel at a specified spectral element using a given extraction aperture pixelset, skipping bad pixels
operaFluxVector operaSpectralOrder::extractFluxElement(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, unsigned indexElem, const PixelSet *aperturePixels) {
	operaFluxVector fluxVector;
	operaFluxVector temp = extractSubpixelFlux(objectImage, nflatImage, biasImage, badpix, 0, gainBiasNoise, indexElem, aperturePixels);
	for(unsigned i=0; i<temp.getlength(); i++) {
		if(!isnan(temp.getflux(i))) fluxVector.insert(temp.getflux(i), temp.getvariance(i));
	}
	return fluxVector;
}

// Extracts the flux per subpixel at a specified spectral element using a given extraction aperture pixelset, including bad pixels as NaN
operaFluxVector operaSpectralOrder::extractSubpixelFlux(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, unsigned badpixValue, GainBiasNoise &gainBiasNoise, unsigned indexElem, const PixelSet *aperturePixels, Vector<unsigned>* pixcol, Vector<unsigned> *pixrow) {
	operaFluxVector fluxVector;
	
	double subpixelArea = aperturePixels->getSubpixelArea(); // Area of each subpixel in pixels
	double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
	double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
	
	for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {
		// Select image col and row of subpixel
		unsigned col = (unsigned)floor(elemXcenter + aperturePixels->getXcenter(pix));
		unsigned row = (unsigned)floor(elemYcenter + aperturePixels->getYcenter(pix));
        if(pixcol) pixcol->insert(col);
		if(pixrow) pixrow->insert(row);
        if(col < objectImage.getnaxis1() && row < objectImage.getnaxis2() && objectImage[row][col] < SATURATIONLIMIT && badpix[row][col] > badpixValue && nflatImage[row][col]) {
			double gain = gainBiasNoise.getGain(col,row);
			double noise = gainBiasNoise.getNoise(col,row);
			
			double pixelFlux = gain*(objectImage[row][col] - biasImage[row][col])/nflatImage[row][col]; // Measured flux in e-/pixel
			double pixelFluxVar = 2.0*noise*noise; // Just count detector noise for now, we can add the rest depending on the extraction method (factor of 2 in detector noise due to bias subtraction)
			
			fluxVector.insert(pixelFlux * subpixelArea, pixelFluxVar * subpixelArea); // Convert from e-/pixel to e-/subpixel, since we have flux values per subpixel
		} else {
			fluxVector.insert(NAN, NAN);
		}
	}
	return fluxVector;
}

// Extracts the background flux per spectral using median binning on extracted fluxes followed by a spline fit.
operaFluxVector operaSpectralOrder::extractBackground(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, unsigned NumberofElementsToBin) {
	// Vectors to hold the background flux and dist for each bin
	operaVector BackgroundDistd;
    operaFluxVector BackgroundFlux;
    
    for(unsigned startindex=0; startindex<SpectralElements->getnSpectralElements(); startindex+=NumberofElementsToBin) {
		// Make sure bin doesn't run past the end of the elements
		unsigned endindex = startindex+NumberofElementsToBin;
		if (endindex > SpectralElements->getnSpectralElements()) endindex = SpectralElements->getnSpectralElements();
        
        // Loop through all left and right background elements in the current bin and put the extracted flux of each subpixel into a vector (extracted variance is ignored)
        operaVector fBackgroundFlux;
		for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
            const PixelSet *backgroundPixels = BackgroundApertures[backgroundIndex]->getSubpixels();
            for(unsigned indexElem=startindex; indexElem < endindex; indexElem++) {
                fBackgroundFlux.append(extractFluxElement(objectImage, nflatImage, biasImage, badpix, gainBiasNoise, indexElem, backgroundPixels).getflux());
            }
        }
        
		// Take the median of the bin -- could use a destrutive median since the vector is never used again...
        if (fBackgroundFlux.size() > 0) {
            BackgroundDistd.insert(fabs(SpectralElements->getdistd(endindex-1) + SpectralElements->getdistd(startindex))/2); //Center dist of the bin
            double BackgroundFluxVal = Median(fBackgroundFlux);
            double BackgroundFluxError = MedianStdDev(fBackgroundFlux, BackgroundFluxVal);
            BackgroundFlux.insert(BackgroundFluxVal, BackgroundFluxError*BackgroundFluxError);
        }
#ifdef PRINT_OUTOFBOUNDS
		else {
			cerr << "operaSpectralOrder::extractStandardSpectrumWarning: NPointsInBkg (" << fBackgroundFlux.size() << ") == 0" << endl;
		}
#endif
	}
	
	// Create a model for the background flux by resampling to the binned flux to the SpectralElements, with a constant variance (mean of all variance bins)
    operaFluxVector BackgroundModelFlux(SpectralElements->getnSpectralElements());
    if(BackgroundFlux.getlength() > 0) {
        BackgroundModelFlux.setflux(fitSpectrum(BackgroundDistd, BackgroundFlux.getflux(), SpectralElements->getDistd()));
        BackgroundModelFlux.setvariance(Mean(BackgroundFlux.getvariance()));
    }
    return BackgroundModelFlux;
}

// Extracts the flux per spectral element for each beam as well as the combined flux aperture. Sets the flux vectors of the spectral elements and beam elements.
void operaSpectralOrder::extractSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, const operaFluxVector& backgroundModelFlux) {
	for (unsigned indexElem=0; indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
        // Total extracted flux and number of points for combined aperture
        double objFlux = 0;
        double objFluxVar = 0;
        unsigned NTotalUsefulPoints = 0;
        unsigned NTotalPoints = 0;
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            const PixelSet *aperturePixels = ExtractionApertures[beam]->getSubpixels();
            
            // Extract the flux in the beam aperture, subtract the background flux element, add up all subpixels in aperture
            operaFluxVector pixelFlux = extractFluxElement(objectImage, nflatImage, biasImage, badpix, gainBiasNoise, indexElem, aperturePixels);
            double objBeamFlux = Sum(pixelFlux.getflux() - backgroundModelFlux.getflux(indexElem));
            double objBeamFluxVar = Sum(pixelFlux.getvariance() + Abs(pixelFlux.getflux()) + backgroundModelFlux.getvariance(indexElem)); // Using total noise = detector noise + photon noise
            
            // Keep running totals for the combined aperture
            objFlux += objBeamFlux;
            objFluxVar += objBeamFluxVar;
            NTotalUsefulPoints += pixelFlux.getlength();
            NTotalPoints += aperturePixels->getNPixels();
            
            // Weight the flux according to the number of extracted subpixels so that missing subpixels won't lower the flux
            objBeamFlux *= (double)aperturePixels->getNPixels()/(double)pixelFlux.getlength();
            objBeamFluxVar *= (double)aperturePixels->getNPixels()/(double)pixelFlux.getlength();            
            
            BeamElements[beam]->setFlux(objBeamFlux,indexElem);
            BeamElements[beam]->setFluxVariance(objBeamFluxVar,indexElem);
        }
        
        // Weight the flux according to the number of extracted subpixels so that missing subpixels won't lower the flux
        objFlux *= (double)NTotalPoints/(double)NTotalUsefulPoints;
        objFluxVar *= (double)NTotalPoints/(double)NTotalUsefulPoints;
        
        SpectralElements->setFlux(objFlux,indexElem);
        SpectralElements->setFluxVariance(objFluxVar,indexElem);
	}
}

void operaSpectralOrder::normalizeFluxToAperture() {
	double npixels = 0;
	for(unsigned beam = 0; beam < numberOfBeams; beam++) npixels += ExtractionApertures[beam]->getSubpixels()->getNPixels();
	npixels /= numberOfBeams;
	for(unsigned beam = 0; beam < numberOfBeams; beam++) {
		const PixelSet *aperturePixels = ExtractionApertures[beam]->getSubpixels();
		BeamElements[beam]->setFluxVector(BeamElements[beam]->getFluxVector()*(npixels/aperturePixels->getNPixels()));
	}
}

// Extract Raw Spectrum - this function uses extraction apertures.
void operaSpectralOrder::extractRawSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction) {
    
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::extractRawSpectrum: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    
    setSpectralElementsByStitchingApertures(effectiveApertureFraction); // Set master spectral elements.
    
    setApertureElements(SpectralElements->getSpectrumType()); // Set spectral elements for each beam aperture and for the background apertures
    
    operaFluxVector backgroundModel(SpectralElements->getnSpectralElements()); // Use empty background for extraction
    
    extractSpectrum(objectImage, nflatImage, biasImage, badpix, gainBiasNoise, backgroundModel);
    
    setSpectrumType(RawBeamSpectrum);
    SpectralElements->setSpectrumType(RawBeamSpectrum);
    SpectralElements->setHasRawSpectrum(true);
    sethasSpectralElements(true);
}


// Extract Raw Spectrum using empirical IP to reject bad pixels
void operaSpectralOrder::extractRawSpectrumRejectingBadpixels(operaFITSImage &objectImage, operaFITSImage &flatImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize, double minSigmaClip, unsigned iterations, bool onTargetProfile, bool usePolynomialFit, bool removeBackground, bool verbose, bool calculateXCorrelation, ostream *pout) {
    
    if (!hasSpectralElements || !SpectralElements->getHasStandardSpectrum()) {
        if(removeBackground) {
            if(verbose) cerr << "operaSpectralOrder::extractRawSpectrumRejectingBadpixels: calculating standard spectrum ..." << endl;
            extractStandardSpectrum(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,BackgroundBinsize);
        } else {
            if(verbose) cerr << "operaSpectralOrder::extractRawSpectrumRejectingBadpixels: calculating standard spectrum without background ..." << endl;
            extractStandardSpectrumNoBackground(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction);
        }
    }
    
    if(calculateXCorrelation) {
        calculateXCorrBetweenIPandImage(objectImage,badpix,NULL);
    }
    
    // Construct spatial profile
    if(onTargetProfile) {
        if(verbose) cerr << "operaSpectralOrder::extractRawSpectrumRejectingBadpixels: measuring spatial profile on object image ..." << endl;
        measureBeamSpatialProfiles(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,usePolynomialFit);
    } else {
        if(verbose) cerr << "operaSpectralOrder::extractRawSpectrumRejectingBadpixels: measuring spatial profile on flat-field image ..." << endl;
        measureBeamSpatialProfiles(flatImage,nflatImage,biasImage,badpix,gainBiasNoise,usePolynomialFit);
    }
    
    ostream *current_pout = NULL;
    
    for(unsigned iter = 0; iter < iterations; iter++){
        if(iter == iterations-1) current_pout = pout;
        if(verbose) cerr << "operaSpectralOrder::extractRawSpectrumRejectingBadpixels: iter="<<iter<<" measuring raw spectrum..." << endl;
        updateBadPixelsToRejectCosmicRays(objectImage,nflatImage,biasImage,badpix,gainBiasNoise, minSigmaClip);
        measureRawSpectrumRejectingBadpixels(objectImage,nflatImage,biasImage,badpix,gainBiasNoise);
		printBeamSpectrum(current_pout);
    }
    
    // get rid of stuff we don't need anymore
    for (unsigned beam=0; beam < numberOfBeams; beam++) {
        if (BeamProfiles[beam]) {
            delete BeamProfiles[beam];
            BeamProfiles[beam] = NULL;
        }
    }
    setSpectrumType(RawBeamSpectrum);
    SpectralElements->setSpectrumType(RawBeamSpectrum);
    SpectralElements->setHasRawSpectrum(true);
    sethasSpectralElements(true);
}

// Extract Standard Spectrum - raw extraction with background subtraction
void operaSpectralOrder::extractStandardSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize) {
    
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::extractStandardSpectrum: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    
    setSpectralElementsByStitchingApertures(effectiveApertureFraction); // Set master spectral elements.
    
    setApertureElements(SpectralElements->getSpectrumType()); // Set spectral elements for each beam aperture and for the background apertures
	
    operaFluxVector backgroundModel = extractBackground(objectImage, nflatImage, biasImage, badpix, gainBiasNoise, BackgroundBinsize); // Extract binned background flux then fit to spectral elements
    
    // Set background elements from background model
    for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
        BackgroundElements[backgroundIndex]->setFluxVector(backgroundModel);
    }
    
    extractSpectrum(objectImage, nflatImage, biasImage, badpix, gainBiasNoise, backgroundModel);
    
    setSpectrumType(StandardBeamSpectrum);
    SpectralElements->setSpectrumType(StandardBeamSpectrum);
    SpectralElements->setHasOptimalSpectrum(true);
    sethasSpectralElements(true);
}

// Extract Standard Spectrum without background subtraction - basically the same as raw extraction, but it also sets dummy background elements so optimal extraction won't get confused.
void operaSpectralOrder::extractStandardSpectrumNoBackground(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction) {
    
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::extractStandardSpectrumNoBackground: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    
    setSpectralElementsByStitchingApertures(effectiveApertureFraction); // Set master spectral elements.
    
    setApertureElements(SpectralElements->getSpectrumType()); // Set spectral elements for each beam aperture and for the background apertures
    
    operaFluxVector backgroundModel(SpectralElements->getnSpectralElements()); // Use empty background for extraction
    
    // Set background elements from background model
    for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
        BackgroundElements[backgroundIndex]->setFluxVector(backgroundModel);
    }
    
    extractSpectrum(objectImage, nflatImage, biasImage, badpix, gainBiasNoise, backgroundModel);
    
    setSpectrumType(StandardBeamSpectrum);
    SpectralElements->setSpectrumType(StandardBeamSpectrum);
    SpectralElements->setHasOptimalSpectrum(true);
    sethasSpectralElements(true);
}

// Extract Optimal Spectrum - this function uses extraction apertures.
void operaSpectralOrder::extractOptimalSpectrum(operaFITSImage &objectImage, operaFITSImage &flatImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize, double minSigmaClip, double sigmaClipRange, unsigned iterations, bool onTargetProfile, bool usePolynomialFit, bool removeBackground, bool verbose, bool calculateXCorrelation, ostream *pout) {
    
    // Extract the standard spectrum if it doesn't exist. STEPS #1 thru #4 K. Horne, 1986
	if(verbose) cerr << "operaSpectralOrder::extractOptimalSpectrum: calculating optimal spectrum using apertures ..." << endl;
	
    if (!hasSpectralElements || !SpectralElements->getHasStandardSpectrum()) {
        if(removeBackground) {
            if(verbose) cerr << "operaSpectralOrder::extractOptimalSpectrum: calculating standard spectrum ..." << endl;
            extractStandardSpectrum(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,BackgroundBinsize);
        } else {
            if(verbose) cerr << "operaSpectralOrder::extractOptimalSpectrum: calculating standard spectrum without background ..." << endl;
            extractStandardSpectrumNoBackground(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction);
        }
    }
    
    if(calculateXCorrelation) {
        calculateXCorrBetweenIPandImage(objectImage,badpix,NULL);
    }
    
    // Construct spatial profile. STEP #5 K. Horne, 1986
    if(onTargetProfile) {
        if(verbose) cerr << "operaSpectralOrder::extractOptimalSpectrum: measuring spatial profile on object image ..." << endl;
        measureBeamSpatialProfiles(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,usePolynomialFit);
    } else {
        if(verbose) cerr << "operaSpectralOrder::extractOptimalSpectrum: measuring spatial profile on flat-field image ..." << endl;
        measureBeamSpatialProfiles(flatImage,nflatImage,biasImage,badpix,gainBiasNoise,usePolynomialFit);
    }
    
    /*
     * The method measureOptimalSpectrum executes the following operations:
     * 1. Revise variance estimates. STEP #6 K. Horne, 1986
     * 2. Mask cosmic ray hits.      STEP #7 K. Horne, 1986
     * 3. Extract optimal spectrum.  STEP #8 K. Horne, 1986     
     */       
    ostream *current_pout = NULL;
    
    for(unsigned iter = 0; iter < iterations; iter++){
        if(iter == iterations-1) current_pout = pout;
        if(verbose) cerr << "operaSpectralOrder::extractOptimalSpectrum: iter="<<iter<<" measuring optimal spectrum..." << endl;
        measureOptimalSpectrum(objectImage,nflatImage,biasImage,badpix,gainBiasNoise, minSigmaClip, sigmaClipRange);
		printBeamSpectrum(current_pout);
    }
    
    // get rid of stuff we don't need anymore
    for (unsigned beam=0; beam < numberOfBeams; beam++) {
        if (BeamProfiles[beam]) {
            delete BeamProfiles[beam];
            BeamProfiles[beam] = NULL;
        }
    }
    if(verbose) cerr << "operaSpectralOrder::extractOptimalSpectrum: optimal extraction run successfully! exiting..." << endl;            
    
    setSpectrumType(OptimalBeamSpectrum);
    SpectralElements->setSpectrumType(OptimalBeamSpectrum);
    SpectralElements->setHasOptimalSpectrum(true);
    sethasSpectralElements(true);
}

void operaSpectralOrder::measureBeamSpatialProfiles(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, bool usePolynomialFit) {
    
    unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
#ifdef PRINT_OUTOFBOUNDS
	if (NumberofElements == 0)
		cerr << "operaSpectralOrder::measureSpatialProfileWithinAperture Warning: NumberofElements (" << NumberofElements << ") == 0" << endl;
#endif
	
    // The IP dimension is given by 1 times the number of subpixels in pixelset
    unsigned NXPoints = 0;
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        unsigned nbeampixels = ExtractionApertures[beam]->getSubpixels()->getNPixels();
        NXPoints += nbeampixels;
        if(BeamProfiles[beam]) delete BeamProfiles[beam];
        BeamProfiles[beam] = new operaInstrumentProfile(nbeampixels,1,1,1, NumberofElements);
    }
    if(InstrumentProfile) delete InstrumentProfile;
    InstrumentProfile = new operaInstrumentProfile(NXPoints,1,1,1, NumberofElements);
        
    for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
        
        double BackgroundFlux = 0;
        for(unsigned background=0;background<LEFTANDRIGHT;background++) {
            BackgroundFlux += BackgroundElements[background]->getFlux(indexElem)/(double)LEFTANDRIGHT;
        }
        
        float distd = (float)SpectralElements->getdistd(indexElem);
        
        unsigned beamstart = 0;
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            operaVector pixelFlux = extractSubpixelFlux(inputImage, nflatImage, biasImage, badpix, 0, gainBiasNoise, indexElem, ExtractionApertures[beam]->getSubpixels()).getflux() - BackgroundFlux;
            for(unsigned pix=0; pix<pixelFlux.size(); pix++) {
                double pixip = pixelFlux[pix];
                BeamProfiles[beam]->setdataCubeValues(pixip,pix,0,indexElem);
                InstrumentProfile->setdataCubeValues(pixip,beamstart+pix,0,indexElem);
            }
            beamstart += pixelFlux.size();
            
            BeamProfiles[beam]->setdistd(distd,indexElem);
            BeamProfiles[beam]->normalizeCubeData(indexElem);
        }
        InstrumentProfile->setdistd(distd,indexElem);
        InstrumentProfile->normalizeCubeData(indexElem);
    }
    
    if(usePolynomialFit) {
        InstrumentProfile->FitPolyMatrixtoIPDataVector(3,false);
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            BeamProfiles[beam]->FitPolyMatrixtoIPDataVector(3,false);
        }
    } else {
        InstrumentProfile->FitMediantoIPDataVector();
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            BeamProfiles[beam]->FitMediantoIPDataVector();
        }
    }
	
    InstrumentProfile->setdataCubeFromPolyModel();
    InstrumentProfile->normalizeCubeData();
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        BeamProfiles[beam]->setdataCubeFromPolyModel();
        BeamProfiles[beam]->normalizeCubeData();
    }
}

void operaSpectralOrder::updateBadPixelsToRejectCosmicRays(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double minSigmaClip) {
    unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
    for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
        double BackgroundFlux=0;
        for(unsigned background=0;background<LEFTANDRIGHT;background++) {
            BackgroundFlux += BackgroundElements[background]->getFlux(indexElem)/(double)LEFTANDRIGHT;
        }
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            double OldBeamFlux = BeamElements[beam]->getFlux(indexElem);
            
            Vector<unsigned> pixelcols, pixelrows; // Store the pixel positions of each extracted subpixel
            operaFluxVector pixelFlux = extractSubpixelFlux(inputImage, nflatImage, biasImage, badpix, 0, gainBiasNoise, indexElem, ExtractionApertures[beam]->getSubpixels(), &pixelcols, &pixelrows) - BackgroundFlux;
            
            for(unsigned pix=0; pix<pixelFlux.getlength(); pix++) {
                const double beamip = BeamProfiles[beam]->getdataCubeValues(pix, 0, indexElem);
                
                if(!isnan(pixelFlux.getflux(pix)) && !isnan(beamip)) {
                    double BeamResidual = pixelFlux.getflux(pix) - OldBeamFlux*beamip;
                    double RevisedBeamVariance = pixelFlux.getvariance(pix) + fabs(BackgroundFlux) + fabs(OldBeamFlux*beamip);
                    if(BeamResidual*BeamResidual < minSigmaClip*RevisedBeamVariance) {
                        badpix[pixelrows[pix]][pixelcols[pix]] += 1.0; // Mark this pixel as good on our mask
                    }
                }
            }
        }
    }
}

void operaSpectralOrder::measureRawSpectrumRejectingBadpixels(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise) {
    
    unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
    for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
        double BackgroundFlux=0;
        for(unsigned background=0;background<LEFTANDRIGHT;background++) {
            BackgroundFlux += BackgroundElements[background]->getFlux(indexElem)/(double)LEFTANDRIGHT;
        }
        
        double OldFlux = SpectralElements->getFlux(indexElem);
        double rawFluxWithoutBadpixelsAllBeams = 0;
        double SumOfUsefulFluxWithinApertureAllBeams = 0;
        double rawFluxVarianceAllBeams = 0;
        
        unsigned beamstart = 0;
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            double OldBeamFlux = BeamElements[beam]->getFlux(indexElem);
            double rawFluxWithoutBadpixels = 0;
            double SumOfUsefulBeamFluxWithinAperture = 0;
            double rawFluxVariance = 0;
            
            operaFluxVector pixelFlux = extractSubpixelFlux(inputImage, nflatImage, biasImage, badpix, 1, gainBiasNoise, indexElem, ExtractionApertures[beam]->getSubpixels()) - BackgroundFlux;
            for(unsigned pix=0; pix<pixelFlux.getlength(); pix++) {
                const double beamip = BeamProfiles[beam]->getdataCubeValues(pix, 0, indexElem);
                const double fullip = InstrumentProfile->getdataCubeValues(beamstart+pix, 0, indexElem);
                if(!isnan(pixelFlux.getflux(pix))) {
                    if(!isnan(beamip)) {
                        double RevisedBeamVariance = pixelFlux.getvariance(pix) + fabs(BackgroundFlux) + fabs(OldBeamFlux*beamip);
                        rawFluxWithoutBadpixels += pixelFlux.getflux(pix);
                        SumOfUsefulBeamFluxWithinAperture += beamip;
                        rawFluxVariance += RevisedBeamVariance;
                    }
                    if(!isnan(fullip)) {
                        double RevisedVariance = pixelFlux.getvariance(pix) + fabs(BackgroundFlux) + fabs(OldFlux*fullip);
                        rawFluxWithoutBadpixelsAllBeams += pixelFlux.getflux(pix);
                        SumOfUsefulFluxWithinApertureAllBeams += fullip;
                        rawFluxVarianceAllBeams += RevisedVariance;
#ifdef PRINT_DEBUG
                        // Uncomment this part to plot the IP-related values for an specific indexElem.
                        if (indexElem==NumberofElements/2) {
                            cout << pix << " " << pixelFlux.getflux(pix) << " " << pixelFlux.getvariance(pix) << " " << Residual << " " << RevisedVariance << " "<< OldFlux << " " << fullip << endl;
                        }
#endif
                    }
                }
            }
            beamstart += pixelFlux.getlength();
            
            if(SumOfUsefulBeamFluxWithinAperture) {
                BeamElements[beam]->setFlux(rawFluxWithoutBadpixels/SumOfUsefulBeamFluxWithinAperture,indexElem);
                BeamElements[beam]->setFluxVariance(rawFluxVariance/SumOfUsefulBeamFluxWithinAperture,indexElem);
            }
        }
        
        if(SumOfUsefulFluxWithinApertureAllBeams) {
            SpectralElements->setFlux(rawFluxWithoutBadpixelsAllBeams/SumOfUsefulFluxWithinApertureAllBeams,indexElem);
            SpectralElements->setFluxVariance(rawFluxVarianceAllBeams/SumOfUsefulFluxWithinApertureAllBeams,indexElem);
        }
    }
}

void operaSpectralOrder::measureOptimalSpectrum(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double minSigmaClip ,double sigmaClipRange) {
    
    unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
    for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
        double BackgroundFlux=0;
        for(unsigned background=0;background<LEFTANDRIGHT;background++) {
            BackgroundFlux += BackgroundElements[background]->getFlux(indexElem)/(double)LEFTANDRIGHT;
        }
        
        double OldFlux = SpectralElements->getFlux(indexElem);
        double optimalFluxDenominatorAllBeams = 0;
        double optimalFluxNumeratorAllBeams = 0;
        double SumOfUsefulFluxWithinApertureAllBeams = 0;
        double maxSigSq = 0;
        
        unsigned beamstart = 0;
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            double OldBeamFlux = BeamElements[beam]->getFlux(indexElem);
            double optimalBeamFluxDenominator = 0;
            double optimalBeamFluxNumerator = 0;
            double SumOfUsefulBeamFluxWithinAperture = 0;
            double maxBeamSigSq = 0;
            
            operaFluxVector pixelFlux = extractSubpixelFlux(inputImage, nflatImage, biasImage, badpix, 0, gainBiasNoise, indexElem, ExtractionApertures[beam]->getSubpixels()) - BackgroundFlux;
            for(unsigned pix=0; pix<pixelFlux.getlength(); pix++) {
                const double beamip = BeamProfiles[beam]->getdataCubeValues(pix, 0, indexElem);
                const double fullip = InstrumentProfile->getdataCubeValues(beamstart+pix, 0, indexElem);
                if(!isnan(pixelFlux.getflux(pix))) {
                    if(!isnan(beamip)) {
                        double BeamResidual = pixelFlux.getflux(pix) - OldBeamFlux*beamip;
                        double RevisedBeamVariance = pixelFlux.getvariance(pix) + fabs(BackgroundFlux) + fabs(OldBeamFlux*beamip);
                        double BeamSigmaSq = BeamResidual*BeamResidual / RevisedBeamVariance;
                        if(BeamSigmaSq > maxBeamSigSq) {
                            maxBeamSigSq = BeamSigmaSq;
                        }
                    }
                    if(!isnan(fullip)) {
                        double Residual = pixelFlux.getflux(pix) - OldFlux*fullip;
                        double RevisedVariance = pixelFlux.getvariance(pix) + fabs(BackgroundFlux) + fabs(OldFlux*fullip);
                        double SigmaSq = Residual*Residual/RevisedVariance;
                        if(SigmaSq > maxSigSq) {
                            maxSigSq = SigmaSq;
                        }
                    }
                }
            }
            for(unsigned pix=0; pix<pixelFlux.getlength(); pix++) {
                const double beamip = BeamProfiles[beam]->getdataCubeValues(pix, 0, indexElem);
                const double fullip = InstrumentProfile->getdataCubeValues(beamstart+pix, 0, indexElem);
                if(!isnan(pixelFlux.getflux(pix))) {
                    if(!isnan(beamip)) {
                        double BeamResidual = pixelFlux.getflux(pix) - OldBeamFlux*beamip;
                        if(BeamResidual < 0) BeamResidual = 0;
                        double RevisedBeamVariance = pixelFlux.getvariance(pix) + fabs(BackgroundFlux) + fabs(OldBeamFlux*beamip);
                        double BeamSigmaSq = BeamResidual*BeamResidual / RevisedBeamVariance;
                        if(BeamSigmaSq < minSigmaClip || BeamSigmaSq < maxBeamSigSq/sigmaClipRange) {
                            optimalBeamFluxNumerator += beamip*pixelFlux.getflux(pix)/RevisedBeamVariance;
                            optimalBeamFluxDenominator += beamip*beamip/RevisedBeamVariance;
                            SumOfUsefulBeamFluxWithinAperture += beamip;
                        } else {
                            BeamProfiles[beam]->setdataCubeValues(NAN, pix, 0, indexElem);
                        }
                    }
                    if(!isnan(fullip)) {
                        double Residual = pixelFlux.getflux(pix) - OldFlux*fullip;
                        if(Residual < 0) Residual = 0;
                        double RevisedVariance = pixelFlux.getvariance(pix) + fabs(BackgroundFlux) + fabs(OldFlux*fullip);
                        double SigmaSq = Residual*Residual/RevisedVariance;
                        if(SigmaSq < minSigmaClip || SigmaSq < maxSigSq/sigmaClipRange) {
                            optimalFluxNumeratorAllBeams += fullip*pixelFlux.getflux(pix)/RevisedVariance;
                            optimalFluxDenominatorAllBeams += fullip*fullip/RevisedVariance;
                            SumOfUsefulFluxWithinApertureAllBeams += fullip;
                        } else {
                            InstrumentProfile->setdataCubeValues(NAN, beamstart+pix, 0, indexElem);
                        }
                    }
                }
            }
            beamstart += pixelFlux.getlength();
            
            if(optimalBeamFluxDenominator) {
                BeamElements[beam]->setFlux(optimalBeamFluxNumerator/optimalBeamFluxDenominator,indexElem);
                BeamElements[beam]->setFluxVariance(SumOfUsefulBeamFluxWithinAperture/optimalBeamFluxDenominator,indexElem);
            }
        }
        
        if(optimalFluxDenominatorAllBeams) {
            SpectralElements->setFlux(optimalFluxNumeratorAllBeams/optimalFluxDenominatorAllBeams,indexElem);
            SpectralElements->setFluxVariance(SumOfUsefulFluxWithinApertureAllBeams/optimalFluxDenominatorAllBeams,indexElem);
        }
    }
}

void operaSpectralOrder::setWavelengthsFromCalibration() {
	if(gethasWavelength() && gethasSpectralElements()) {
		SpectralElements->setwavelengthsFromCalibration(Wavelength);
		for (unsigned beam = 0; beam < numberOfBeams; beam++) {
			BeamElements[beam]->setWavelength(SpectralElements->getWavelength());
		}
	}
}

void operaSpectralOrder::deleteSpectralEnergyDistribution(void) {
    if(SpectralEnergyDistribution) delete SpectralEnergyDistribution;
	SpectralEnergyDistribution = NULL;
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        if(BeamSED[beam]) delete BeamSED[beam];            
        BeamSED[beam] = NULL;
    } 
}

operaSpectralEnergyDistribution *operaSpectralOrder::getBeamSED(unsigned beam) {
#ifdef RANGE_CHECK
	if (beam >= numberOfBeams) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamSED[beam];
}

const operaSpectralEnergyDistribution *operaSpectralOrder::getBeamSED(unsigned beam) const {
#ifdef RANGE_CHECK
	if (beam >= numberOfBeams) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamSED[beam];
}

operaSpectralElements& operaSpectralOrder::MainAndBeamElements(unsigned index) {
	if (index) return *BeamElements[index-1];
	return *SpectralElements;
}

operaSpectralEnergyDistribution& operaSpectralOrder::MainAndBeamSED(unsigned index) {
	if (index) return *BeamSED[index-1];
	return *SpectralEnergyDistribution;
}

unsigned operaSpectralOrder::MainAndBeamCount() {
	return numberOfBeams+1;
}

operaVector operaSpectralOrder::getMainAndBeamFluxes(unsigned index) {
	operaVector fluxes;
	for (unsigned b = 0; b<numberOfBeams+1; b++) fluxes.insert(MainAndBeamElements(b).getFlux(index));
	return fluxes;
}

void operaSpectralOrder::fitSEDUncalibratedFluxToSample(const operaSpectrum& uniformSample) {
	// Fit a spline to the uniform sample along the wavelength of our SpectralElements and store the result in uncalibratedModel.
	operaSpectrum uncalibratedModel = fitSpectrum(uniformSample, SpectralElements->getWavelength());
	
	// Set the flux of UncalibratedFlux of the SEDs to the values in uncalibratedModel.
	for (unsigned i = 0; i < MainAndBeamCount() && i < uncalibratedModel.fluxcount(); i++) {
		MainAndBeamSED(i).setCalibrationWavelength(SpectralElements->getWavelength());
		MainAndBeamSED(i).setUncalibratedFlux(uncalibratedModel.fluxvector(i));
		MainAndBeamSED(i).setHasUncalibratedFlux(true);
	}
}

void operaSpectralOrder::fitSEDFcalAndThroughputToFlatResp(const operaSpectrum& flatresp) {
	// Fit a spline to the flat response along the wavelength of our SpectralElements and store the result in flatResponseModel.
	operaSpectrum flatResponseModel = fitSpectrum(flatresp, SpectralElements->getWavelength());
	
	operaFluxVector fluxcal = 1.0/flatResponseModel.operafluxvector();
	
	for (unsigned i = 0; i < MainAndBeamCount(); i++) {
		MainAndBeamSED(i).setCalibrationWavelength(SpectralElements->getWavelength());
		MainAndBeamSED(i).setFluxCalibration(fluxcal);
		MainAndBeamSED(i).setThroughput(fluxcal);
		//MainAndBeamSED(i).sethasCalibrationWavelength(true);
		MainAndBeamSED(i).setHasFluxCalibration(true);
		MainAndBeamSED(i).setHasInstrumentThroughput(true);
	}
	sethasSpectralEnergyDistribution(true);
}

void operaSpectralOrder::applyNormalization(unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial, bool overwriteUncalFlux, bool normalizeBeams) {
	if (binsize == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    
	for (unsigned i = 0; i < MainAndBeamCount(); i++) {
		unsigned numberOfElements = MainAndBeamElements(i).getnSpectralElements();
		if (numberOfElements == 0) {
			throw operaException("operaSpectralOrder: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
		}
		operaFluxVector normalizedFlux(numberOfElements);
		operaFluxVector continuum(numberOfElements);
		normalizeSpectrum(MainAndBeamElements(i).getFluxVector(), normalizedFlux, continuum, binsize, orderOfPolynomial, usePolynomial);
		if (overwriteUncalFlux) MainAndBeamElements(i).setFluxVector(normalizedFlux);
		
		if (!normalizeBeams) break;
    }
	
    if (gethasWavelength() && SpectralElements->getHasDistance()) {
        SpectralElements->setwavelengthsFromCalibration(getWavelength());
    }
}

void operaSpectralOrder::applyNormalizationForEmissionSpectrum(unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial, bool overwriteUncalFlux, bool normalizeBeams) {
	if (binsize == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    
    for (unsigned i = 0; i < MainAndBeamCount(); i++) {
		double maxflux = Max(MainAndBeamElements(i).getFluxVector().getflux());
		operaFluxVector flux = 1.0 - (MainAndBeamElements(i).getFluxVector() / maxflux);
		
		unsigned numberOfElements = MainAndBeamElements(i).getnSpectralElements();
		operaFluxVector normalizedFlux(numberOfElements);
		operaFluxVector continuum(numberOfElements);
		normalizeSpectrum(flux, normalizedFlux, continuum, binsize, orderOfPolynomial, usePolynomial);
        
		if (overwriteUncalFlux) MainAndBeamElements(i).setFluxVector(1.0 - normalizedFlux);
		else MainAndBeamElements(i).setFluxVector(1.0 - flux);
		
		if (!normalizeBeams) break;
    }
	
    if (gethasWavelength() && SpectralElements->getHasDistance()) {
        SpectralElements->setwavelengthsFromCalibration(getWavelength());
    }
}

void operaSpectralOrder::applyNormalizationFromExistingContinuum(bool normalizeBeams) {
	for (unsigned i = 0; i < MainAndBeamCount(); i++) {
		if (MainAndBeamElements(i).getnSpectralElements() == 0) {
			throw operaException("operaSpectralOrder: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
		}
		MainAndBeamElements(i).setFluxVector(MainAndBeamElements(i).getFluxVector() / MainAndBeamSED(i).getUncalibratedFlux());
		if (!normalizeBeams) break;
	}
}

void operaSpectralOrder::normalizeSpectrum(const operaFluxVector &uncalibratedFlux, operaFluxVector &normalizedFlux, operaFluxVector &outputContinuum, unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial) {
	measureContinuum(uncalibratedFlux,outputContinuum, binsize,NSIGMACUT,orderOfPolynomial,usePolynomial);    
    normalizedFlux = uncalibratedFlux/outputContinuum;
}

void operaSpectralOrder::measureContinuum(const operaFluxVector &uncalibratedFlux, operaFluxVector &outputContinuum, unsigned binsize, unsigned nsigcut, unsigned orderOfPolynomial, bool usePolynomial) {
	
	operaVector continuumElemSample, continuumFluxSample;
    measureFluxContinuum(uncalibratedFlux, binsize, nsigcut, continuumElemSample, continuumFluxSample);
    
	operaVector continuumModelElem;
    for(unsigned i = 0; i < uncalibratedFlux.getlength(); i++) continuumModelElem.insert(i);
    
    if(usePolynomial) {
        const operaVector& continuumModelFlux = fitSpectrumToPolynomial(continuumElemSample, continuumFluxSample, continuumModelElem, orderOfPolynomial);
        outputContinuum.setflux(continuumModelFlux);
        outputContinuum.setvariance(uncalibratedFlux.getvariance());
    } else {
        const operaVector& continuumModelFlux = fitSpectrum(continuumElemSample, continuumFluxSample, continuumModelElem);
		outputContinuum.setflux(continuumModelFlux);
		outputContinuum.setvariance(uncalibratedFlux.getvariance());
    }
#ifdef PRINT_DEBUG           
    for(unsigned i=0;i<NumberofPoints;i++) {           
        cout << i << " " 
        << uncalibratedFlux.getflux(i) << " "
        << outputContinuum.getflux(i) << " " 
        << outputContinuum.getvariance(i) << " " << endl;   
    }
#endif
}

void operaSpectralOrder::calculateContinuum(unsigned binsize, unsigned nsigcut) {
	if (binsize == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    // use info in SpectralElements and binsize to create spectralEnergyDistribution
    for (unsigned i = 0; i < MainAndBeamCount(); i++) {
		MainAndBeamSED(i).setCalibrationDist(MainAndBeamElements(i).getDistd());
		MainAndBeamSED(i).setCalibrationWavelength(MainAndBeamElements(i).getWavelength());
		MainAndBeamSED(i).setUncalibratedFlux(MainAndBeamElements(i).getFluxVector());
		MainAndBeamSED(i).calculateUncalibratedFlux(binsize, nsigcut);
    }

    sethasSpectralEnergyDistribution(true);
}

void operaSpectralOrder::normalizeSpectralElementsByConstant(const operaVector& maxFluxForNormalization) {
	if (maxFluxForNormalization.size() != MainAndBeamCount()) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    
	for (unsigned i = 0; i < MainAndBeamCount(); i++) {
        if (maxFluxForNormalization[i] == 0) {
            throw operaException("operaSpectralOrder: ", operaErrorDivideByZeroError, __FILE__, __FUNCTION__, __LINE__);
        }
		MainAndBeamElements(i).setFluxVector(MainAndBeamElements(i).getFluxVector() / maxFluxForNormalization[i]);
    }
}

void operaSpectralOrder::divideSpectralElementsBySEDElements(bool useThroughput) {
    if(!gethasSpectralEnergyDistribution()) {
        throw operaException("operaSpectralOrder::divideSpectralElementsBySEDElements: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    
	for (unsigned i = 0; i < MainAndBeamCount(); i++) {
		MainAndBeamElements(i).setFluxVector(MainAndBeamElements(i).getFluxVector() / MainAndBeamSED(i).getCalibration(useThroughput));
    }
}

void operaSpectralOrder::multiplySpectralElementsBySEDElements(bool useThroughput, const operaVector& spectralBinConstants) {
    if(!gethasSpectralEnergyDistribution()) {
        throw operaException("operaSpectralOrder::multiplySpectralElementsBySEDElements: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    
	for (unsigned i = 0; i < spectralBinConstants.size(); i++) {
		MainAndBeamElements(i).setFluxVector(MainAndBeamElements(i).getFluxVector() * MainAndBeamSED(i).getCalibration(useThroughput) / spectralBinConstants[i]);
    }
}

/*
 * Star+Sky Mode, subtract the sky beam
 */
void operaSpectralOrder::calculateStarAndSkyElements(bool starplusskyInvertSkyFiber, double skyOverStarFiberAreaRatio) {  
    // Test if number of beams is 2.
    if(numberOfBeams != 2) { 
        throw operaException("operaSpectralOrder:calculateStarAndSkyElements: numberOfBeams must be 2;", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
	
	// Below it assumes there is only two beams, by default beam=0 is star+sky and beam=1 is sky
	operaSpectralElements* starPlusSkyElements = BeamElements[0];
	operaSpectralElements* skyElements = BeamElements[1];
    if(starplusskyInvertSkyFiber) {
		skyElements = BeamElements[0];
		starPlusSkyElements = BeamElements[1];
    }
    
	SpectralElements->setFluxVector(starPlusSkyElements->getFluxVector() - skyElements->getFluxVector() / skyOverStarFiberAreaRatio);
    
	createSkyElements(SpectralElements->getnSpectralElements(), SpectrumType);
	SkyElements->setPhotoCenter(SpectralElements->getPhotoCenterX(),SpectralElements->getPhotoCenterY());
	SkyElements->setDistd(SpectralElements->getDistd());
	SkyElements->setWavelength(SpectralElements->getWavelength());
    SkyElements->setFluxVector(skyElements->getFluxVector());
    sethasSkyElements(true);
}

/*
 * Radial Velocity Wavelength Correction
 */
void operaSpectralOrder::applyRvelWavelengthCorrection(double RVcorrectionInKmPerSecond) {
    SpectralElements->setWavelength(SpectralElements->getWavelength() * (1.0 + (RVcorrectionInKmPerSecond / SPEED_OF_LIGHT_KMS)));
}

void operaSpectralOrder::setExtendedRvelWavelengthCorrection(double RVcorrectionInKmPerSecond) {
    if(!SpectralElements->getHasExtendedBeamFlux()) throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    SpectralElements->setRvel(SpectralElements->getWavelength()*(RVcorrectionInKmPerSecond / SPEED_OF_LIGHT_KMS));
}

void operaSpectralOrder::applyWavelengthCorrectionFromExtendedRvel(void) {
    if(!SpectralElements->getHasExtendedBeamFlux()) throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	SpectralElements->setWavelength(SpectralElements->getWavelength() + SpectralElements->getRvel());
}

void operaSpectralOrder::CreateExtendedVectors() {
	for (unsigned i = 0; i < MainAndBeamCount(); i++) MainAndBeamElements(i).createExtendedVectors();
}

void operaSpectralOrder::CopyFluxVectorIntoRawFlux() {
	for (unsigned i = 0; i < MainAndBeamCount(); i++) MainAndBeamElements(i).copyTOrawFlux();
}

void operaSpectralOrder::CopyRawFluxIntoFluxVector() {
	for (unsigned i = 0; i < MainAndBeamCount(); i++) MainAndBeamElements(i).copyFROMrawFlux();
}

void operaSpectralOrder::CopyFluxVectorIntoFcalFlux() {
	for (unsigned i = 0; i < MainAndBeamCount(); i++) MainAndBeamElements(i).copyTOfcalFlux();
}

void operaSpectralOrder::CopyFcalFluxIntoFluxVector() {
	for (unsigned i = 0; i < MainAndBeamCount(); i++) MainAndBeamElements(i).copyFROMfcalFlux();
}

void operaSpectralOrder::CopyFluxVectorIntoNormalizedFlux() {
	for (unsigned i = 0; i < MainAndBeamCount(); i++) MainAndBeamElements(i).copyTOnormalizedFlux();
}

void operaSpectralOrder::CopyNormalizedFluxIntoFluxVector() {
	for (unsigned i = 0; i < MainAndBeamCount(); i++) MainAndBeamElements(i).copyFROMnormalizedFlux();
}

void operaSpectralOrder::TrimOrderToWavelengthRange() {
	operaIndexRange indexrange = SpectralElements->getContainedWavelengthSubrangeIndexes(minwavelength, maxwavelength);
	for (unsigned i = 0; i < MainAndBeamCount(); i++) MainAndBeamElements(i).trim(indexrange);
	if (gethasPolarimetry()) Polarimetry->trim(indexrange);
}

operaSpectralLineList operaSpectralOrder::getRawLinesFromUncalSpectrum(double linewidth, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, double nsigclip, dispersionaxis_t dispersiontype) {
	
	if (!SpectralElements->getHasXCorrelation()) {
		const operaVector& xvals = (dispersiontype == wavelength_disp ? SpectralElements->getWavelength() : SpectralElements->getDistd());
		const operaVector& compXcorr = calculateXCorrWithGaussian(xvals, SpectralElements->getFluxVector().getflux(), linewidth);
		SpectralElements->setXCorrelation(compXcorr);
		SpectralElements->setHasXCorrelation(true);
	}
	
	operaSpectralLines compLines = DetectSpectralLines(SpectralElements, linewidth, LocalMaxFilterWidth, MinPeakDepth, DetectionThreshold, dispersiontype);
	
	operaSpectralLineList filteredLines = getSpectralLinesNearMedianWidth(compLines, nsigclip);
	
	if (filteredLines.size() > 0) {
		SpectralElements->setHasWavelength(true);
	}
	return filteredLines;
}

operaSpectralLineList operaSpectralOrder::detectSpectralLines(double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, double nsigclip, double spectralResolution, bool emissionSpectrum) {
    if (!hasSpectralElements || !hasWavelength) {
		throw operaException("detectSpectralLinesInSpectralOrder: order has no elements/wavelength. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
	}
	
	SpectralElements->setwavelengthsFromCalibration(Wavelength);
	
	operaFluxVector saveFlux;
	if (!emissionSpectrum) {
		saveFlux = SpectralElements->getFluxVector();
		operaFluxVector invertedFlux = (1.0 + DetectionThreshold) - saveFlux;
		for(unsigned i = 0; i < invertedFlux.getlength(); i++) if(isnan(invertedFlux.getflux(i)) || invertedFlux.getflux(i) < 0) invertedFlux.setflux(0.0, i);
		SpectralElements->setFluxVector(invertedFlux);
	}

	double linewidth = Wavelength->getcentralWavelength() / spectralResolution;
	
	operaSpectralLineList filteredLines = getRawLinesFromUncalSpectrum(linewidth, LocalMaxFilterWidth, MinPeakDepth, DetectionThreshold, nsigclip, wavelength_disp);
	
	// Transform spectrum back to absorption
	if (!emissionSpectrum) {
		SpectralElements->setFluxVector(saveFlux);
		filteredLines.amplitude = 1.0 - filteredLines.amplitude;
	}

	doubleValue_t Resolution;
	Resolution.value = Wavelength->getcentralWavelength() / filteredLines.medianWidth;
	Resolution.error = filteredLines.medianWidthError * Wavelength->getcentralWavelength() / (filteredLines.medianWidth*filteredLines.medianWidth);
	Wavelength->setSpectralResolution(Resolution);
	
	return filteredLines;
}
