/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaTelescope
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

#include "operaError.h"
#include "globaldefines.h"
#include "libraries/operaException.h"
#include "libraries/operaTelescope.h"  // for operaTelescope

/*!
 * operaTelescope
 * \author Eder Martioli
 * \brief Encapsulation of telescope information.
 * \file operaTelescope.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaTelescope
 * \brief Encapsulation of Wavelength information.
 * \return none
 */

/*
 * Constructors
 */

operaTelescope::operaTelescope() :
latitude(0.0),
longitude(0.0),
elevation(0.0),
CollectingArea(0.0),
Aperture(0.0),
FocalRatio(1.0),
TelescopeMount(undefinedMount),
OpticalCoating(undefinedCoating)
{
}

/*
 * Destructor
 */

operaTelescope::~operaTelescope() {
  
}

/*
 * Methods for managing data
 */

void operaTelescope::setlatitude(double Latitude) {
    latitude = Latitude;
}

double operaTelescope::getlatitude(void) {
    return latitude;
}

void operaTelescope::setlongitude(double Longitude) {
    longitude = Longitude;
}

double operaTelescope::getlongitude(void) {
    return longitude;
}

void operaTelescope::setelevation(double Elevation) {
    elevation = Elevation;
}

double operaTelescope::getelevation(void) {
    return elevation;
}

void operaTelescope::setCollectingArea(double area) {
    CollectingArea = area;
}

double operaTelescope::getCollectingArea(void) {
    return CollectingArea;
}

void operaTelescope::setAperture(double diameter) {
    Aperture = diameter;
}

double operaTelescope::getAperture(void) {
    return Aperture;
}

void operaTelescope::setFocalRatio(double fnumber) {
    FocalRatio = fnumber;
}

double operaTelescope::getFocalRatio(void){
    return FocalRatio;
}

void operaTelescope::setTelescopeMount(telescopeMount_t telescopeMountType) {
    TelescopeMount = telescopeMountType;
}

telescopeMount_t operaTelescope::getTelescopeMount(void) {
    return TelescopeMount;
}

void operaTelescope::setOpticalCoating(opticalCoating_t coatingType) {
    OpticalCoating = coatingType;
}

opticalCoating_t operaTelescope::getOpticalCoating(void) {
    return OpticalCoating;
}

