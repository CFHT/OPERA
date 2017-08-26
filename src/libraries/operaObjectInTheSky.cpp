/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaObjectInTheSky
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Feb/2013
 
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

#include <math.h>       // for pow

#include "operaError.h"
#include "globaldefines.h"
#include "libraries/operaException.h"
#include "libraries/operaObjectInTheSky.h"      // for operaObjectInTheSky
#include "libraries/operaSpectralTools.h"       // for PlanckFunction and calculateBlackBodyVFlux

#include "libraries/operaLibCommon.h"       // for VEGA_PHOTONFLUX_AT_548NM

/*!
 * operaObjectInTheSky
 * \author Eder Martioli
 * \brief This class encapsulates the information of an object in the sky.
 * \file operaObjectInTheSky.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaObjectInTheSky
 * \brief Encapsulation of information about an astronomical source.
 * \return none
 */

/*
 * Constructors
 */
operaObjectInTheSky::operaObjectInTheSky() :
RA(0.0),
Dec(0.0),
ProperMotionRA(0.0),
ProperMotionDec(0.0),
Parallax(0.0),
V_magnitude(0.0),
EffectiveTemperature(0.0),
RadialVelocity(0.0),
SpectralType(A_type)
{
    
}

/*
 * Destructor
 */

operaObjectInTheSky::~operaObjectInTheSky() {
  
}

/*
 * Methods for managing data
 */

void operaObjectInTheSky::setsourceID(string SourceID) {
    sourceID = SourceID;
}

string operaObjectInTheSky::getsourceID(void) {
    return sourceID;
}


void operaObjectInTheSky::setRA(double alpha) {
    RA = alpha;
}

double operaObjectInTheSky::getRA(void) {
    return RA;
}

void operaObjectInTheSky::setDec(double delta) {
    Dec = delta;
}

double operaObjectInTheSky::getDec(void) {
    return Dec;
}

void operaObjectInTheSky::setProperMotionRA(double muRA) {
    ProperMotionRA = muRA;
}

double operaObjectInTheSky::getProperMotionRA(void) {
    return ProperMotionRA;
}

void operaObjectInTheSky::setProperMotionDec(double muDec) {
    ProperMotionDec = muDec;
}

double operaObjectInTheSky::getProperMotionDec(void) {
    return ProperMotionDec;
}

void operaObjectInTheSky::setParallax(double pi) {
    Parallax = pi;
}

double operaObjectInTheSky::getParallax(void) {
    return Parallax;
}

void operaObjectInTheSky::setV_magnitude(double Vmag) {
    V_magnitude = Vmag;
}

double operaObjectInTheSky::getV_magnitude(void) {
    return V_magnitude;
}

void operaObjectInTheSky::setEffectiveTemperature(double Teff) {
    EffectiveTemperature = Teff;
}

double operaObjectInTheSky::getEffectiveTemperature(void) {
    return EffectiveTemperature;
}

void operaObjectInTheSky::setRadialVelocity(double RV) {
    RadialVelocity = RV;
}

double operaObjectInTheSky::getRadialVelocity(void) {
    return RadialVelocity;
}

void operaObjectInTheSky::setSpectralType(operaSpectralType_t SpecTy) {
    SpectralType = SpecTy;
}

operaSpectralType_t operaObjectInTheSky::getSpectralType(void) {
    return SpectralType;
}

// Calculate the flux in the V-Jonhnson band in ph/(m^2 s nm)
double operaObjectInTheSky::getBlackBodyVFlux(void) {
    return calculateBlackBodyVFlux(EffectiveTemperature);
}


// The phi factor is given by (Radius/distance)^2
double operaObjectInTheSky::calculatePhiFactor(void) {
    double phiVega = VEGA_PHOTONFLUX_AT_548NM/PlanckFunction(VEGA_EFF_TEMPERATURE_K,548e-9); // calculate factor phi ~ R^2/d^2 for Vega
    double phi = pow(10,-(V_magnitude - VEGA_VBAND_MAGNITUDE)/2.5)*(phiVega*calculateBlackBodyVFlux(VEGA_EFF_TEMPERATURE_K))/getBlackBodyVFlux();
    return phi;
}

double operaObjectInTheSky::getSpectralBinFlux(double wl0, double wlf) {
    return (calculatePhiFactor()*IntegrateSpectralElementOfBlackBody(wl0, wlf,EffectiveTemperature));
}

