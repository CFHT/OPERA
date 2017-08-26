/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: PixelSet
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
#include "libraries/operaLib.h"		// for itos
#include "libraries/operaException.h"
#include "libraries/PixelSet.h"
#include "libraries/operaVectorOperations.h"

/*! 
 * PixelSet
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates a set of pixels or subpixels
 * \file PixelSet.cpp
 * \ingroup libraries
 */

using namespace std;


/* 
 * PixelSet
 * \brief It defines a set of pixels or subpixels where  
 * \brief the flux information is extracted from
 * \return none
 */

/*
 * PixelSet Constructor
 */

PixelSet::PixelSet(void) : nPixels(0), subpixelArea(1) { }

PixelSet::PixelSet(float SubPixelArea) : nPixels(0), subpixelArea(SubPixelArea) { }

PixelSet::PixelSet(unsigned NPixels, float SubPixelArea) :
nPixels(NPixels),
xcenter(NPixels),
ycenter(NPixels),
iIndex(NPixels),
jIndex(NPixels),
redundancy(NPixels),
pixelValue(NPixels),
subpixelArea(SubPixelArea)
{
}

/*
 * PixelSet Setters/Getters
 */

unsigned PixelSet::getNPixels(void) const {
    return nPixels;
}

float PixelSet::getXcenter(unsigned index) const {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return xcenter[index];
}
float PixelSet::getYcenter(unsigned index) const {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return ycenter[index];
}
int PixelSet::getiIndex(unsigned index) const {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return iIndex[index];    
}
int PixelSet::getjIndex(unsigned index) const {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return jIndex[index];    
}
unsigned PixelSet::getredundancy(unsigned index) const {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return redundancy[index];    
}     
float PixelSet::getPixelValue(unsigned index) const {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return pixelValue[index];
}
float PixelSet::getSubpixelArea(void) const {
    return subpixelArea;
}

void PixelSet::setXcenter(float Xcenter, unsigned index) {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    xcenter[index] = Xcenter;
}

void PixelSet::setYcenter(float Ycenter, unsigned index) {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    ycenter[index] = Ycenter;
}

void PixelSet::setiIndex(int iindex, unsigned k) {
#ifdef RANGE_CHECK
    if (k >= nPixels) {
		throw operaException("PixelSet: k="+itos(k)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    iIndex[k] = iindex;
}

void PixelSet::setjIndex(int jindex, unsigned k) {
#ifdef RANGE_CHECK
    if (k >= nPixels) {
		throw operaException("PixelSet: k="+itos(k)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    jIndex[k] = jindex;
}

void PixelSet::setredundancy(unsigned Redundancy, unsigned k) {
#ifdef RANGE_CHECK
    if (k >= nPixels) {
		throw operaException("PixelSet: k="+itos(k)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    redundancy[k] = Redundancy;
}

void PixelSet::setPixelValue(float PixelValue, unsigned index) {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	pixelValue[index] = PixelValue;
}

void PixelSet::setSubpixelArea(float SubpixelArea) {
    subpixelArea = SubpixelArea;
}

/*
 * PixelSet Methods
 */

void PixelSet::resize(unsigned newsize) {
	xcenter.resize(newsize);
    ycenter.resize(newsize);
    iIndex.resize(newsize);
    jIndex.resize(newsize);
    redundancy.resize(newsize);
    pixelValue.resize(newsize);
    nPixels = newsize;
}

float PixelSet::getMinXcoord(void) const {    
    return Min(xcenter);
}

float PixelSet::getMaxXcoord(void) const {
    return Max(xcenter);               
}

float PixelSet::getMinYcoord(void) const {
    return Min(ycenter);               
}

float PixelSet::getMaxYcoord(void) const {
    return Max(ycenter);               
}
