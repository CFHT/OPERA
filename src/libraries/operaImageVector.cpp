/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaImageVector
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
 
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
#include "libraries/operaGeometricShapes.h"		// for Box
#include "libraries/operaImageVector.h"
#include "libraries/operaFITSImage.h"

/*!
 * operaImageVector
 * \author Doug Teeple
 * \brief The operaImageVector stores a vector of image indices and the associated image
 * \file operaImageVector.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * where
 * \brief where function
 * \brief return an index vector of count size of all image values at index that satisfy condition.
 * \param b - operaFITSImage*
 * \param count - count of 1's
 * \return operaFITSImage&
 */	
operaImageVector &where(operaFITSImage* b, unsigned *count) {
	float *bp = (float *)b->getpixels(); 
	unsigned c = 0, d = 0;
	unsigned n = b->getnpixels(); 
	unsigned nn = n; 
	while (n--) if (*bp++ != 0.0) c++;
	*count = c;
	operaImageVector *ivector = new operaImageVector(c);
	bp = (float *)b->getpixels();
	while (nn--) if (*bp++ != 0.0) ivector->setIndex(d++, nn);
	//if (b->getIstemp()) delete b;
	return *ivector;
}

operaImageVector &where(operaFITSImage& b, unsigned *count) {
	float *bp = (float *)b.getpixels(); 
	unsigned c = 0, d = 0;
	unsigned n = b.getnpixels(); 
	unsigned nn = n; 
	while (n--) if (*bp++ != 0.0) c++;
	*count = c;
	operaImageVector *ivector = new operaImageVector(c);
	bp = (float *)b.getpixels();
	while (nn--) if (*bp++ != 0.0) ivector->setIndex(d++, nn);
	//if (b.getIstemp()) delete &b;
	return *ivector;
}

operaImageVector &where(operaFITSImage* b) {
	unsigned count;
	return where(b, &count);
}

operaImageVector &where(operaFITSImage& b) {
	unsigned count;
	return where(b, &count);
}

void SetImage(operaImageVector *vector, operaFITSImage *image) {
	vector->setFITSImage(image);
}
 /*
  * the image vector
 */

operaImageVector::operaImageVector(void) :
length(0),
index(NULL),
associatedImage(NULL)
{
}

operaImageVector::operaImageVector(unsigned Length) :
length(0),
index(NULL),
associatedImage(NULL)
{
	length = Length;
	index = new unsigned[Length];
}


operaImageVector::operaImageVector(Box *boxes, unsigned count) :
length(0),
index(NULL),
associatedImage(NULL)
{
	for (unsigned i=0; i<count; i++) {
		length += boxes[i].getSize();
	}
	index = new unsigned[length];
	for (unsigned k=0; k<count; k++) {
		unsigned i = 0;
		Box box = boxes[k];
		for (unsigned x=box.getX1(); x<=box.getX2(); x++) {
			for (unsigned y=box.getY1(); y<=box.getY2(); y++) {
				index[i++] = x + y * 2048;
			}
		}		
	}
}

operaImageVector::operaImageVector(Box &box) :
length(0),
index(NULL),
associatedImage(NULL)
{
	length = box.getSize();
	index = new unsigned[length];
	unsigned i = 0;
	for (unsigned x=box.getX1(); x<=box.getX2(); x++) {
		for (unsigned y=box.getY1(); y<=box.getY2(); y++) {
			index[i++] = x + y * 2048;
		}
	}
}

void operaImageVector::setIndex(unsigned Index, unsigned IndexValue) {
#ifdef RANGE_CHECK
    if (Index > length) {
		throw operaException("operaImageVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	index[Index] = IndexValue;
}

operaFITSImage *operaImageVector::getFITSImage(void) {
	return associatedImage;
}

void operaImageVector::setFITSImage(operaFITSImage *image) {
	associatedImage = image;
}

unsigned operaImageVector::getlength(void) {
	return length;
}

operaImageVector::~operaImageVector(void)
{
	if (index) {
		delete[] index;
	}
	index = NULL;
	associatedImage = NULL;	// do not delete...
}

