/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaPolarimetry
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
#include "libraries/operaPolarimetry.h"

/*!
 * \file operaPolarimetry.cpp
 */

operaPolarimetry::operaPolarimetry() :
length(0),
hasStokesI(false),
hasStokesQ(false),
hasStokesU(false),
hasStokesV(false),
hasDegreeOfStokesI(false),
hasDegreeOfStokesQ(false),
hasDegreeOfStokesU(false),
hasDegreeOfStokesV(false),
hasContinuumRemoved(false),
hasFirstNullPolarization(false),
hasSecondNullPolarization(false),
hasWavelength(false)
{ }

operaPolarimetry::operaPolarimetry(unsigned Length) :
length(Length),
stokesParameter(Length),
degreeOfPolarization(Length),
continuumRemoved(Length),
firstNullPolarization(Length),
secondNullPolarization(Length),
wavelength(Length),
hasStokesI(false),
hasStokesQ(false),
hasStokesU(false),
hasStokesV(false),
hasDegreeOfStokesI(false),
hasDegreeOfStokesQ(false),
hasDegreeOfStokesU(false),
hasDegreeOfStokesV(false),
hasContinuumRemoved(false),
hasFirstNullPolarization(false),
hasSecondNullPolarization(false),
hasWavelength(false)
{
	if (length == 0) {
		throw operaException("operaPolarimetry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
}

void operaPolarimetry::resize(unsigned Length)
{
	stokesParameter.resize(Length);
    degreeOfPolarization.resize(Length);
    continuumRemoved.resize(Length);
    firstNullPolarization.resize(Length);
    secondNullPolarization.resize(Length);
    wavelength.resize(Length);
	length = Length;
}

void operaPolarimetry::trim(operaIndexRange range)
{
	stokesParameter.trim(range);
    degreeOfPolarization.trim(range);
    continuumRemoved.trim(range);
    firstNullPolarization.trim(range);
    secondNullPolarization.trim(range);
    wavelength.trim(range);
	length = range.size();
}

unsigned operaPolarimetry::getLength(void) const
{
    return length;
}

double operaPolarimetry::getwavelength(unsigned indexElem) const
{
	return wavelength[indexElem];
}

void operaPolarimetry::setwavelength(double Wavelength, unsigned indexElem)
{
	wavelength[indexElem] = Wavelength;
}

const operaFluxVector &operaPolarimetry::getStokesParameter(stokes_parameter_t StokesIndex) const
{
    return stokesParameter.getStokesParameter(StokesIndex);
}

void operaPolarimetry::setStokesParameter(stokes_parameter_t StokesIndex, const operaFluxVector& FluxVector)
{
#ifdef RANGE_CHECK
    if (FluxVector.size() != length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    stokesParameter.setStokesParameter(StokesIndex, FluxVector);
    setHasStokes(StokesIndex,true);
}

double operaPolarimetry::getStokesParameterFlux(stokes_parameter_t StokesIndex, unsigned index) const {
	return stokesParameter.getStokesParameterFlux(StokesIndex, index);
}

double operaPolarimetry::getStokesParameterVariance(stokes_parameter_t StokesIndex, unsigned index) const {
	return stokesParameter.getStokesParameterVariance(StokesIndex, index);
}

void operaPolarimetry::setStokesParameter(stokes_parameter_t StokesIndex, double StokesValue, double Variance, unsigned index)
{
    stokesParameter.setStokesParameter(StokesIndex, StokesValue, Variance, index);
    setHasStokes(StokesIndex,true);
}

const operaFluxVector &operaPolarimetry::getDegreeOfPolarization(stokes_parameter_t StokesIndex) const
{
    return degreeOfPolarization.getStokesParameter(StokesIndex);
}

void operaPolarimetry::setDegreeOfPolarization(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector)
{
#ifdef RANGE_CHECK
    if (FluxVector.size() != length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    degreeOfPolarization.setStokesParameter(StokesIndex, FluxVector);
    setHasDegreeOfStokes(StokesIndex,true);
}

double operaPolarimetry::getDegreeOfPolarizationFlux(stokes_parameter_t StokesIndex, unsigned index) const
{
	return degreeOfPolarization.getStokesParameterFlux(StokesIndex, index);
}

double operaPolarimetry::getDegreeOfPolarizationVariance(stokes_parameter_t StokesIndex, unsigned index) const
{
	return degreeOfPolarization.getStokesParameterVariance(StokesIndex, index);
}

void operaPolarimetry::setDegreeOfPolarization(stokes_parameter_t StokesIndex, double DegreeOfPolarizationValue, double Variance, unsigned index)
{
    degreeOfPolarization.setStokesParameter(StokesIndex, DegreeOfPolarizationValue, Variance, index);
    setHasDegreeOfStokes(StokesIndex,true);
}

const operaFluxVector &operaPolarimetry::getContinuumRemoved(stokes_parameter_t StokesIndex) const
{
    return continuumRemoved.getStokesParameter(StokesIndex);
}

void operaPolarimetry::setContinuumRemoved(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector)
{
#ifdef RANGE_CHECK
    if (FluxVector.size() != length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    continuumRemoved.setStokesParameter(StokesIndex, FluxVector);
    hasContinuumRemoved = true;
}

double operaPolarimetry::getContinuumRemovedFlux(stokes_parameter_t StokesIndex, unsigned index) const
{
	return continuumRemoved.getStokesParameterFlux(StokesIndex, index);
}

double operaPolarimetry::getContinuumRemovedVariance(stokes_parameter_t StokesIndex, unsigned index) const
{
	return continuumRemoved.getStokesParameterVariance(StokesIndex, index);
}

void operaPolarimetry::setContinuumRemoved(stokes_parameter_t StokesIndex, double ContinuumRemovedValue, double Variance, unsigned index)
{
    continuumRemoved.setStokesParameter(StokesIndex, ContinuumRemovedValue, Variance, index);
    hasContinuumRemoved = true;
}

const operaFluxVector &operaPolarimetry::getFirstNullPolarization(stokes_parameter_t StokesIndex) const
{
    return firstNullPolarization.getStokesParameter(StokesIndex);
}

void operaPolarimetry::setFirstNullPolarization(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector)
{
#ifdef RANGE_CHECK
    if (FluxVector.size() != length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    firstNullPolarization.setStokesParameter(StokesIndex, FluxVector);
    hasFirstNullPolarization = true;
}

double operaPolarimetry::getFirstNullPolarizationFlux(stokes_parameter_t StokesIndex, unsigned index) const
{
	return firstNullPolarization.getStokesParameterFlux(StokesIndex, index);
}

double operaPolarimetry::getFirstNullPolarizationVariance(stokes_parameter_t StokesIndex, unsigned index) const
{
	return firstNullPolarization.getStokesParameterVariance(StokesIndex, index);
}

void operaPolarimetry::setFirstNullPolarization(stokes_parameter_t StokesIndex, double FirstNullPolarizationValue, double Variance, unsigned index)
{
    firstNullPolarization.setStokesParameter(StokesIndex, FirstNullPolarizationValue, Variance, index);
    hasFirstNullPolarization = true;
}

const operaFluxVector &operaPolarimetry::getSecondNullPolarization(stokes_parameter_t StokesIndex) const
{
    return secondNullPolarization.getStokesParameter(StokesIndex);
}

void operaPolarimetry::setSecondNullPolarization(stokes_parameter_t StokesIndex, const operaFluxVector &FluxVector)
{
#ifdef RANGE_CHECK
    if (FluxVector.size() != length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    secondNullPolarization.setStokesParameter(StokesIndex, FluxVector);
    hasSecondNullPolarization = true;
}

double operaPolarimetry::getSecondNullPolarizationFlux(stokes_parameter_t StokesIndex, unsigned index) const
{
	return secondNullPolarization.getStokesParameterFlux(StokesIndex, index);
}

double operaPolarimetry::getSecondNullPolarizationVariance(stokes_parameter_t StokesIndex, unsigned index) const
{
	return secondNullPolarization.getStokesParameterVariance(StokesIndex, index);
}

void operaPolarimetry::setSecondNullPolarization(stokes_parameter_t StokesIndex, double SecondNullPolarizationValue, double Variance, unsigned index)
{
    secondNullPolarization.setStokesParameter(StokesIndex, SecondNullPolarizationValue, Variance, index);    
    hasSecondNullPolarization = true;
}


bool operaPolarimetry::getHasStokes(stokes_parameter_t StokesIndex) const
{
    switch (StokesIndex) {
        case StokesI:
            return hasStokesI;
        case StokesQ:
            return hasStokesQ;
        case StokesU:
            return hasStokesU;
        case StokesV:
            return hasStokesV;
        default:
#ifdef RANGE_CHECK
            throw operaException("operaPolarimetry: unrecognized Stokes parameter. ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
#endif
            return false;
    }
}

void operaPolarimetry::setHasStokes(stokes_parameter_t StokesIndex, bool HasStokes)
{
    switch (StokesIndex) {
        case StokesI:
            hasStokesI = HasStokes;
            break;
        case StokesQ:
            hasStokesQ = HasStokes;
            break;
        case StokesU:
            hasStokesU = HasStokes;
            break;
        case StokesV:
            hasStokesV = HasStokes;
            break;
        default:
#ifdef RANGE_CHECK
            throw operaException("operaPolarimetry: unrecognized Stokes parameter. ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
#endif
			break;
    }
}

bool operaPolarimetry::getHasDegreeOfStokes(stokes_parameter_t StokesIndex) const
{
    switch (StokesIndex) {
        case StokesI:
            return hasDegreeOfStokesI;
        case StokesQ:
            return hasDegreeOfStokesQ;
        case StokesU:
            return hasDegreeOfStokesU;
        case StokesV:
            return hasDegreeOfStokesV;
        default:
#ifdef RANGE_CHECK
            throw operaException("operaPolarimetry: unrecognized Stokes parameter. ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
#endif
            return false;
    }
}

void operaPolarimetry::setHasDegreeOfStokes(stokes_parameter_t StokesIndex, bool HasDegreeOfStokes)
{
    switch (StokesIndex) {
        case StokesI:
            hasDegreeOfStokesI = HasDegreeOfStokes;
            break;
        case StokesQ:
            hasDegreeOfStokesQ = HasDegreeOfStokes;
            break;
        case StokesU:
            hasDegreeOfStokesU = HasDegreeOfStokes;
            break;
        case StokesV:
            hasDegreeOfStokesV = HasDegreeOfStokes;
            break;
        default:
#ifdef RANGE_CHECK
            throw operaException("operaPolarimetry: unrecognized Stokes parameter. ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);
#endif
            break;
    }
}

void operaPolarimetry::calculatePolarization(void)
{
	//Don't bother with StokesI since the degree of polarization should always be 1.0
	calculatePolarization(StokesQ);
	calculatePolarization(StokesU);
	calculatePolarization(StokesV);
}

void operaPolarimetry::calculatePolarization(stokes_parameter_t StokesIndex)
{
    if (hasStokesI) {
        if (getHasDegreeOfStokes(StokesIndex)) {
            operaFluxVector fluxprod = degreeOfPolarization.getStokesParameter(StokesIndex) * stokesParameter.getStokesParameter(StokesI);
            stokesParameter.setStokesParameter(StokesIndex, fluxprod);
            setHasStokes(StokesIndex,true);
        }
    } else {
        throw operaException("operaPolarimetry: missing Stokes I.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
}

void operaPolarimetry::calculateDegreeOfPolarization(void)
{
    calculateDegreeOfPolarization(StokesI); //This should always result in a degree of polarization of 1.0
    calculateDegreeOfPolarization(StokesQ);
    calculateDegreeOfPolarization(StokesU);
    calculateDegreeOfPolarization(StokesV);
}

void operaPolarimetry::calculateDegreeOfPolarization(stokes_parameter_t StokesIndex)
{
    if (hasStokesI) {
        if (getHasStokes(StokesIndex)) {
			operaFluxVector fluxdiv = stokesParameter.getStokesParameter(StokesIndex) / stokesParameter.getStokesParameter(StokesI);
            degreeOfPolarization.setStokesParameter(StokesIndex, fluxdiv);
            setHasStokes(StokesIndex,true);
        }
    } else {
        throw operaException("operaPolarimetry: missing Stokes I.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
}

void operaPolarimetry::calculateDegreeOfPolarization(stokes_parameter_t StokesIndex, operaFluxVector *iE[4], operaFluxVector *iA[4], unsigned NumberOfExposures) {
    if(NumberOfExposures != 2 && NumberOfExposures != 4) {
        throw operaException("operaPolarimetry: NumberOfExposures not valid.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    if (StokesIndex != StokesQ && StokesIndex != StokesU && StokesIndex != StokesV) {
        throw operaException("operaPolarimetry: unrecognized Stokes parameter.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
	
	operaFluxVector PoverI(length);
    operaFluxVector N1(length),N2(length);
    
    if(method == Difference || method == DifferenceWithBeamSwapped) {
        operaFluxVector G[4];
        
        for(unsigned i=0;i<NumberOfExposures;i++) {
            if(method == DifferenceWithBeamSwapped) {
				if(i==0 || i==2) {
                    G[i] = (*iE[i] - *iE[i+1])/(*iE[i] + *iE[i+1]);
                } else if (i==1 || i==3) {
                    G[i] = (*iA[i-1] - *iA[i])/(*iA[i-1] + *iA[i]);
                }
            } else {
				G[i] = (*iE[i] - *iA[i])/(*iE[i] + *iA[i]);
            }
        }
		
        if(NumberOfExposures==2) {
            operaFluxVector D1 = G[0] - G[1];
            PoverI = D1 / NumberOfExposures;
        } else if (NumberOfExposures==4) {
            operaFluxVector D1 = G[0] - G[1];
            operaFluxVector D2 = G[2] - G[3];
            operaFluxVector D1s = G[0] - G[3];
            operaFluxVector D2s = G[2] - G[1];
            
            PoverI = (D1 + D2) / NumberOfExposures;
            N1 = (D1 - D2) / NumberOfExposures;
            N2 = (D1s - D2s) / NumberOfExposures;
        }
    } else if (method == Ratio) {
        operaFluxVector r[4];
        
        for(unsigned i=0;i<NumberOfExposures;i++) {
            r[i] = (*iE[i]) / (*iA[i]);
        }
		
		if(NumberOfExposures==2) {
            operaFluxVector R1 = r[0] / r[1];
            operaFluxVector R = Pow(R1 , 1.0/NumberOfExposures);
            PoverI = (R - 1.0) / (R + 1.0);
        } else if (NumberOfExposures==4) {
            operaFluxVector R1 = r[0] / r[1];
            operaFluxVector R2 = r[2] / r[3]; // changed from r[0] / r[1] Apr 15 2013 DT
            operaFluxVector R1s = r[0] / r[3];
            operaFluxVector R2s = r[2] / r[1];
            operaFluxVector R = Pow(R1 * R2, 1.0/NumberOfExposures);
            operaFluxVector RN1 = Pow(R1 / R2 , 1.0/NumberOfExposures);
            operaFluxVector RN2 = Pow(R1s / R2s , 1.0/NumberOfExposures);
            PoverI = (R - 1.0) / (R + 1.0);
            N1 = (RN1 - 1.0) / (RN1 + 1.0);
            N2 = (RN2 - 1.0) / (RN2 + 1.0);
        }
    } else if (method == NewMethod) {
        PoverI = 0.0;
        N1 = 0.0;
        N2 = 0.0;
    } else {
        throw operaException("operaPolarimetry: unrecognized method to calculate polarization.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);        
    }
    
    // The fix below is necessary since the angles of the analyzer with respect to the reference system for ESPaDOnS don't follow the same order as described in the literature. Added Dec 04 2013 EM
    if(StokesIndex == StokesQ) {
        PoverI *= -1.0;
    }
    
    setDegreeOfPolarization(StokesIndex, PoverI);
    if (NumberOfExposures==4) {
        setFirstNullPolarization(StokesIndex, N1);
        setSecondNullPolarization(StokesIndex, N2);
    }
}

void operaPolarimetry::calculateStokesParameter(stokes_parameter_t StokesIndex, operaFluxVector *iE[4], operaFluxVector *iA[4], unsigned NumberOfExposures) {
    
    if(NumberOfExposures != 2 && NumberOfExposures != 4) {
        throw operaException("operaPolarimetry: NumberOfExposures not valid.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaFluxVector Intensity(length);
    
    Intensity = 0.0;
    for(unsigned i=0;i<NumberOfExposures;i++) {
        Intensity += (*(iE[i]) + *(iA[i])) / ((double)NumberOfExposures); // DT Apr 26 2013 added *2.0 Deleted by LMa Feb 12 2016
    }
    setStokesParameter(StokesI, Intensity);
    
    if (StokesIndex == StokesQ || StokesIndex == StokesU || StokesIndex == StokesV) {
        if(!getHasDegreeOfStokes(StokesIndex)) {
            calculateDegreeOfPolarization(StokesIndex,iE,iA,NumberOfExposures);
        }
        calculatePolarization(StokesIndex);
    }
}
