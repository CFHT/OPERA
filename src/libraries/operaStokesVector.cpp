/*******************************************************************
 ****               	    OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaStokesVector
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaStokesVector.h"

/*!
 * \file operaStokesVector.cpp
 * \brief This file holds the implementation of the class operaStokesVector.
 */

string StokesName(stokes_parameter_t StokesParameter)
{
    string StokesParameterName;
    
    switch (StokesParameter) {
        case StokesI:
            StokesParameterName = "Stokes I";
            break;
        case StokesQ:
            StokesParameterName = "Stokes Q";
            break;
        case StokesU:
            StokesParameterName = "Stokes U";
            break;
        case StokesV:
            StokesParameterName = "Stokes V";
            break;
        default:
            break;
    }
    
    return StokesParameterName;
}

operaStokesVector::operaStokesVector() : length(0), stokesVector(4) { }

operaStokesVector::operaStokesVector(unsigned Length) : length(Length)
{
	for(unsigned i = 0; i < 4; i++) stokesVector.push_back(operaFluxVector(length));
}

operaStokesVector::operaStokesVector(const operaFluxVector &StokesVectorI, const operaFluxVector &StokesVectorQ, const operaFluxVector &StokesVectorU, const operaFluxVector &StokesVectorV) : length(0)
{
	stokesVector.push_back(StokesVectorI);
	stokesVector.push_back(StokesVectorQ);
	stokesVector.push_back(StokesVectorU);
	stokesVector.push_back(StokesVectorV);
	length = StokesVectorI.getlength();
}

void operaStokesVector::resize(unsigned Length)
{
	for(unsigned i = 0; i < 4; i++) stokesVector[i].resize(Length);
    length = Length;
}

void operaStokesVector::trim(operaIndexRange range)
{
	for(unsigned i = 0; i < 4; i++) stokesVector[i].trim(range);
    length = range.size();
}

unsigned operaStokesVector::getLength(void) const
{
    return length;
}

void operaStokesVector::setStokesParameters(const operaFluxVector &StokesIFluxVector, const operaFluxVector &StokesQFluxVector, const operaFluxVector &StokesUFluxVector, const operaFluxVector &StokesVFluxVector)
{
    stokesVector[StokesI] = StokesIFluxVector;
    stokesVector[StokesQ] = StokesQFluxVector;
    stokesVector[StokesU] = StokesUFluxVector;
    stokesVector[StokesV] = StokesVFluxVector;
    length = StokesIFluxVector.getlength();
}

void operaStokesVector::setStokesParameter(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector)
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    stokesVector[StokesIndex] = FluxVector;
}

void operaStokesVector::setStokesParameter(stokes_parameter_t StokesIndex, double StokesValue, double Variance, unsigned index)
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    stokesVector[StokesIndex].setflux(StokesValue, index);
    stokesVector[StokesIndex].setvariance(Variance, index);
}

const operaFluxVector &operaStokesVector::getStokesParameter(stokes_parameter_t StokesIndex) const
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector[StokesIndex];
}

double operaStokesVector::getStokesParameterFlux(stokes_parameter_t StokesIndex, unsigned index) const
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector[StokesIndex].getflux(index);
}

double operaStokesVector::getStokesParameterVariance(stokes_parameter_t StokesIndex, unsigned index) const
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector[StokesIndex].getvariance(index);
}
