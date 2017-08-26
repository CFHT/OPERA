#ifndef OPERAIMAGEVECTOR_H
#define OPERAIMAGEVECTOR_H
/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name:  operaImageVector
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include "libraries/operaGeometricShapes.h"		// for box
#include "libraries/operaFITSImage.h"			// for operaFITSImage

/*! 
 * \sa class operaImageVector
 * \brief The operaImageVector stores a vector of image indices and the associated image.
 * \file operaImageVector.h
 * \ingroup libraries
 */
class operaFITSImage;
class operaFITSSubImage;

class operaImageVector {
	
private:
	unsigned length;
	unsigned *index;
	operaFITSImage *associatedImage;
	
public:
    operaImageVector(void);
	
    operaImageVector(unsigned length);
		
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return a pixel value.
	 * \param i - index into row
	 * \note usage: myfitsimage[y][x] = 0.0; sets the image pixel value at coords x, y to zero
	 * \return float pointer to row 
	 */
	unsigned operator[](unsigned i) {return index[i];};
	
	/*! 
	 * \brief operator =
	 * \brief copy value.
	 * \return operaImageVector^ 
	 */
	operaImageVector& operator=(float f) {
		float *pixels = (float *)associatedImage->getpixels();
		unsigned i = 0;
		while (i < length) {
			pixels[index[i++]] = f;
		}
		if (associatedImage->getIstemp()) {
			delete associatedImage;
			associatedImage = NULL;
		}
		return *this;
	};
	
	operaImageVector(Box *boxes, unsigned count);	
	
	operaImageVector(Box &box);
	
	operaFITSImage *getFITSImage(void);
	
	void setFITSImage(operaFITSImage *image);
	
	void setIndex(unsigned Index, unsigned indexvalue);
	
	unsigned getIndex(unsigned Index) {return index[Index];};
	
	unsigned getlength(void);
	
	~operaImageVector(void);
 };

void SetImage(operaImageVector *vector, operaFITSImage *image);

operaImageVector &where(operaFITSImage* b, unsigned *count);

operaImageVector &where(operaFITSImage& b, unsigned *count);

operaImageVector &where(operaFITSImage* b);

operaImageVector &where(operaFITSImage& b);

operaImageVector &where(operaFITSSubImage* b, unsigned *count);

operaImageVector &where(operaFITSSubImage& b, unsigned *count);

operaImageVector &where(operaFITSSubImage* b);

operaImageVector &where(operaFITSSubImage& b);

#endif
