/*******************************************************************
 ****               		OPERA PIPELINE v1.0                     ****
 ********************************************************************
 Library name: GainBiasNoise
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
#include "libraries/GainBiasNoise.h"

/*!
 * GainBiasNoise
 * \author Doug Teeple / Eder Martioli
 * \brief GainBiasNoise  calculation
 * \details Calculate the gain / bias / noise in an image.
 * \file GainBiasNoise.cpp
 * \ingroup libraries
 */

/*
 * Constructors / Destructors
 */

GainBiasNoise::GainBiasNoise() :
namps(MAXAMPS)
{
	for (unsigned i=0; i<MAXAMPS; i++) {
		setGain(i, 0.0);
		setNoise(i, 0.0);
		ampsDataSec[i].x1 = 0;
		ampsDataSec[i].y1 = 0;
		ampsDataSec[i].x2 = 0;
		ampsDataSec[i].y2 = 0;
	}
}

GainBiasNoise::GainBiasNoise(unsigned Namps)
{
	if (Namps > MAXAMPS) {
		throw operaException("GainBiasNoise ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	namps = Namps;
	for (unsigned i=0; i<MAXAMPS; i++) {
		setGain(i, 0.0);
		setNoise(i, 0.0);
		ampsDataSec[i].x1 = 0;
		ampsDataSec[i].y1 = 0;
		ampsDataSec[i].x2 = 0;
		ampsDataSec[i].y2 = 0;
	}
}

GainBiasNoise::~GainBiasNoise() {

}

/* 
 * \sa method double getDatasec();
 * \brief returns the datasec ot amp
 */
void GainBiasNoise::getDatasec(unsigned amp, DATASEC_t &datasec) const {
	if (amp >= namps) {
		throw operaException("GainBiasNoise ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	memcpy((void *)&datasec, (void *)&ampsDataSec[amp], sizeof(DATASEC_t));
}
/* 
 * \sa method void setDatasec(unsigned amp, DATASEC_t &datasec);
 * \brief sets the datsec of amp
 */
void GainBiasNoise::setDatasec(unsigned amp, DATASEC_t &datasec) {
	if (amp >= namps) {
		throw operaException("GainBiasNoise ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	memcpy((void *)&ampsDataSec[amp], (void *)&datasec, sizeof(DATASEC_t));
}

/*!
 * \sa method double getNoise();
 * \brief returns the Noise of amp
 */
double GainBiasNoise::getNoise(unsigned x, unsigned y) const {
    unsigned current_amp = 0;
    for(unsigned amp=0;amp<namps; amp++) {
        if(x >= ampsDataSec[amp].x1 && x <= ampsDataSec[amp].x2 &&
           y >= ampsDataSec[amp].y1 && y <= ampsDataSec[amp].y2) {
            current_amp=amp;
            break;
        }
    }
    return noises[current_amp];
}

/*!
 * \sa method double getGain();
 * \brief returns the Gain of amp
 */
double GainBiasNoise::getGain(unsigned x, unsigned y) const {
    unsigned current_amp = 0;
    for(unsigned amp=0;amp<namps; amp++) {
        if(x >= ampsDataSec[amp].x1 && x <= ampsDataSec[amp].x2 &&
           y >= ampsDataSec[amp].y1 && y <= ampsDataSec[amp].y2) {
            current_amp=amp;
            break;
        }
    }
    return gains[current_amp];
}

