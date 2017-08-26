#ifndef OPERAFITSCUBE_H
#define OPERAFITSCUBE_H
/*******************************************************************
 ****               	OPERA PIPELINE v1.0                     ****
 *******************************************************************
 Class: operaFITSCube 
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

#include <string>
#include <fitsio.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"			// for operaFITSImage
#include "libraries/operaStats.h"

/*
 * operaFITSCube class
 * \author Doug Teeple
 * \brief This class encapsulates the FITS Cube.
 * \file operaFITSCube.h
 * \ingroup libraries
 */
class operaFITSCube : public operaFITSImage {
	
private:
		
protected:	
	
public:
	/*!
	 * \sa class operaFITSCube()
	 * \brief Basic operaFITSCube constructor.
	 */
	operaFITSCube();														// simply construct a default image
	/*!
	 * \sa class operaFITSCube
	 * \brief construct an in-memory FITSImage object
	 * \brief operaFITSCube(unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, edatatype Datatype=tushort, unsigned Compression=0)
	 * \brief Create an in-memory image of given dimensions.
	 */
	operaFITSCube(unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, edatatype Datatype = tushort);	// simply construct an image of a given size
	/*! 
	 * \sa class operaFITSCube
	 * \brief create a writeable file image of an in memory FITSImage object
	 * \brief operaFITSCube(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, edatatype Datatype, unsigned Compression, bool IsLazy = false)
	 * \brief Constructor to create a new FITS file.
	 */
	operaFITSCube(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, edatatype Datatype, unsigned Compression = 0, bool IsLazy = false);	// create a new empty FITSImage
	/*! 
	 * \sa class operaFITSCube
	 * \brief create a FITSIMage object from a FITS file
	 * \brief operaFITSCube(string Filename, int mode=READWRITE, unsigned Compression = 0, bool IsLazy = false)
	 * \brief Constructor to create a FITSImage from a FITS file.
	 */
	operaFITSCube(string Filename, edatatype Datatype, int mode=READWRITE/*READONLY*/, unsigned Compression = 0, bool IsLazy = false);		// read an existing FITSImage from file
    /*! 
     * operaFITSCube(operaFITSCube &imageIn, bool ViewOnly = false, bool AddHeader = false)
     * \brief Clone a operaFITSCube object.
     */
	operaFITSCube(operaFITSCube &image, bool ViewOnly = false, bool AddHeader = false);
	
	~operaFITSCube();   // destructor	
    
    /* OPERATORS */
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return a operaFITSImage*.
	 */
	operaFITSImage& operator[](unsigned i) {
#ifdef FAST_SLICES
		setSlice(i);
		return (operaFITSImage& )*this;
#else
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (i > this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		operaFITSCube *t = new operaFITSCube(*this, true);
		t->setSlice(i);
		t->npixels = t->npixels_per_slice;
		t->naxis3 = 1;
		t->naxis = 2;
		return *t;
#endif
    };
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return an ImageVector.
	 * \note usage: myfitsimage[where(input>saturation)] = 0.0; sets the image pixel values where true to zero
	 */
	operaImageVector *operator[](operaImageVector *vector) {
        SetImage(vector, this);
        return vector;
    };
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return an ImageVector.
	 * \note usage: myfitsimage[where(input>saturation)] = 0.0; sets the image pixel values where true to zero
	 */
	operaImageVector &operator[](operaImageVector &vector) {
        SetImage(&vector, this);
        return vector;
    };
    /*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \note usage: operaFITSImage a = operaFITSImage b; copies the pixel values from b to a
	 */	
	operaFITSCube& operator=(operaFITSCube* b) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		float *bp = (float *)b->pixptr+(b->current_slice-1)*b->naxis1*b->naxis2; 
		unsigned n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b->naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		while (n--) *p++ = *bp++;
		if (b->istemp && !b->viewOnly) delete b;
		return *this;
	};
	/*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \note usage: operaFITSImage a = operaFITSImage b; copies the pixel values from b to a
	 */	
	operaFITSCube& operator=(operaFITSCube& b) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		while (n--) *p++ = *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief Median operator %
	 * \note usage: operaFITSCube& a % 8; median combines 8 slices of a
	 */	
	operaFITSCube& operator%(unsigned i) {
		this->medianCollapse(i);
		return *this;
	};
	/*! 
	 * \brief operator ==
	 * \brief equality operator.
	 * \note usage: operaFITSCube a = operaFITSCube b * operaFITSCube c; multiplies the pixel values  b * c and assigns to a
	 */	
	operaFITSCube& operator==(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) {*tp++ = (*p++ == *bp++?1.0:0.0);}
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
    /*! 
     * \brief operator ==
     * \brief equality operator.
     * \note usage:
     */
    operaFITSCube& operator==(float f) {
		operaFITSCube *t = new operaFITSCube(*this);
 		float *tp = (float *)t->pixptr; 
        float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		unsigned n = npixels; 
		while (n--) {*tp++ = (*p++ == f?1.0:0.0);}
        return *this;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \note usage: operaFITSCube a += operaFITSCube b; adds the pixel values from b to a
	 */	
	operaFITSCube& operator+=(operaFITSCube& b) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		while (n--) *p++ += *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \note usage: operaFITSCube a += 100.0; adds the float value to a
	 */	
	operaFITSCube& operator+=(float f) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		unsigned n = npixels; 
		while (n--) *p++ += f;
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \note usage: operaFITSCube a -= operaFITSCube b; subtracts the pixel values from b to a
	 */	
	operaFITSCube& operator-=(operaFITSCube& b) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		while (n--) *p++ -= *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \note usage: operaFITSCube a -= 100.0; subtracts the float value from a
	 */	
	operaFITSCube& operator-=(float f) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		unsigned n = npixels; 
		while (n--) *p++ -= f;
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \note usage: operaFITSCube a *= operaFITSCube b; multiplies the pixel values from b to a
	 */	
	operaFITSCube& operator*=(operaFITSCube& b) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		while (n--) *p++ *= *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \note usage: operaFITSCube a *= 100.0; multiplies the float value times a
	 */	
	operaFITSCube& operator*=(float f) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		unsigned n = npixels; 
		while (n--) *p++ *= f;
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \note usage: operaFITSCube a /= operaFITSCube b; divides the pixel values from b into a
	 */	
	operaFITSCube& operator/=(operaFITSCube& b) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		while (n--) *p++ /= (*bp!=0.0?*bp++:(bp++,1.0));
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \note usage: operaFITSCube a /= 100.0; divides the float value from a
	 */	
	operaFITSCube& operator/=(float f) {
		float *p = (float *)pixptr+(current_slice-1)*naxis1*naxis2; 
		unsigned n = npixels; 
		while (n--) *p++ /= (f!=0.0?f:1.0);
		return *this;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \note usage: operaFITSCube a = operaFITSCube b * operaFITSCube c; multiplies the pixel values  b * c and assigns to a
	 */	
	operaFITSCube& operator*(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ * *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \note usage: operaFITSCube a = operaFITSCube b * 10.0; multiplies the pixel values  b * 10.0 and assigns to a
	 */	
	operaFITSCube& operator*(float f) {
		operaFITSCube *t = new operaFITSCube(*this, this->viewOnly);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ * f;
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \note usage: operaFITSCube a = operaFITSCube b / operaFITSCube c; divides the pixel values  b / c and assigns to a
	 */	
	operaFITSCube& operator/(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, this->viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ / (*bp!=0.0?*bp++:(bp++,1.0));
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \note usage: operaFITSCube a = operaFITSCube b / 100.0; divides the pixel values  b / 100.0 and assigns to a
	 */	
	operaFITSCube& operator/(float f) {
		operaFITSCube *t = new operaFITSCube(*this, this->viewOnly);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ / (f!=0.0?f:1.0);
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \note usage: operaFITSCube a = operaFITSCube b + operaFITSCube c; adds the pixel values  b + c and assigns to a
	 */	
	operaFITSCube& operator+(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ + *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \note usage: operaFITSCube a = operaFITSCube b + 100.0; adds the pixel values  b + 100.0 and assigns to a
	 */	
	operaFITSCube& operator+(float f) {
		operaFITSCube *t = new operaFITSCube(*this, this->viewOnly);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ + f;
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \note usage: operaFITSCube a = operaFITSCube b - operaFITSCube c; subtracts the pixel values  b - c and assigns to a
	 */	
	operaFITSCube& operator-(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ - *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \note usage: operaFITSCube a = operaFITSCube b - 100.0; subtracts the pixel values  b - 100.0 and assigns to a
	 */	
	operaFITSCube& operator-(float f) {
		operaFITSCube *t = new operaFITSCube(*this, this->viewOnly);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ - f;
		return *t;
	};
	/*! 
	 * \brief operator !
	 * \brief invert operator.
	 * \note usage: operaFITSCube a = !operaFITSCube b; inverts the pixel values of b and assigns to a (0.0 becomes 1.0, non-zero becomes 0.0)
	 */	
	operaFITSCube& operator!() {
		float *tp;
		operaFITSCube *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++==0.0?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator >
	 * \brief greater than operator.
	 * \note usage: operaFITSCube a = operaFITSCube b > operaFITSCube c; creates a mask of 1.0 if b > c or 0.0 if b <= c and assigns to a
	 */	
	operaFITSCube& operator>(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++>*bp++?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator >
	 * \brief greater than operator.
	 * \note usage: operaFITSCube a = operaFITSCube b > 100.0; creates a mask of 1.0 if b > 100.0 or 0.0 if b <= 100.0 and assigns to a
	 */	
	operaFITSCube& operator>(float f) {
		operaFITSCube *t = new operaFITSCube(*this, this->viewOnly);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++>f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator >=
	 * \brief greater than or equal operator.
	 * \note usage: operaFITSCube a = operaFITSCube b >= operaFITSCube c; creates a mask of 1.0 if b >= c or 0.0 if b < c and assigns to a
	 */	
	operaFITSCube& operator>=(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++>=*bp++?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator >=
	 * \brief greater than or equal operator.
	 * \note usage: operaFITSCube a = operaFITSCube b >= 100.0; creates a mask of 1.0 if b >= 100.0 or 0.0 if b < 100.0 and assigns to a
	 */	
	operaFITSCube& operator>=(float f) {
		operaFITSCube *t = new operaFITSCube(*this, this->viewOnly);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++>=f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator <
	 * \brief less than operator.
	 * \note usage: operaFITSCube a = operaFITSCube b < operaFITSCube c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 */	
	operaFITSCube& operator<(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++<*bp++?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator <
	 * \brief less than operator.
	 * \note usage: operaFITSCube a = operaFITSCube b < operaFITSCube c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 */	
	operaFITSCube& operator<(float f) {
		operaFITSCube *t = new operaFITSCube(*this, this->viewOnly);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++<f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator <=
	 * \brief less than or equal operator.
	 * \note usage: operaFITSCube a = operaFITSCube b <= operaFITSCube c; creates a mask of 1.0 if b <= c or 0.0 if b > c and assigns to a
	 */	
	operaFITSCube& operator<=(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++<=*bp++?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator <=
	 * \brief less than or equal operator.
	 * \note usage: operaFITSCube a = operaFITSCube b <= 100.0; creates a mask of 1.0 if b <= 100.0 or 0.0 if b > 100.0 and assigns to a
	 */	
	operaFITSCube& operator<=(float f) {
		operaFITSCube *t = new operaFITSCube(*this, this->viewOnly);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++<f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator &&
	 * \brief less than operator.
	 * \note usage: operaFITSCube a = operaFITSCube b && operaFITSCube c; creates a mask of 1.0 where b && c != 0.0
	 */	
	operaFITSCube& operator&&(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++!=0.0&&*bp++!=0.0?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator ||
	 * \brief less than operator.
	 * \note usage: operaFITSCube a = operaFITSCube b && operaFITSCube c; creates a mask of 1.0 where b && c != 0.0
	 */	
	operaFITSCube& operator||(operaFITSCube& b) {
		float *tp;
		operaFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSCube(*this, viewOnly);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr+(b.current_slice-1)*b.naxis1*b.naxis2; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++!=0.0||*bp++!=0.0?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * operaFITSCubeSave() 
	 * \brief Saves the current image to disk.
	 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
	 */
	void operaFITSCubeSave();
	/*! 
	 * operaFITSCubeeSaveAs(string newFilename) 
	 * \brief Saves the current image to disk, with the given filename.
	 */
	void operaFITSCubeSaveAs(string newFilename);
	/*! 
	 * void operaFITSCubeConvertImage(edatatype todatattype)
	 * \brief Convert image type and create new image values of that type.
	 */
	void operaFITSCubeConvertImage(edatatype todatatype);
    /*!
     * void medianCollapse()
     * \brief collapse the cube into a single image
     * WARNING: that this method utilizes operaArrayMedianQuick, so the input image will be scrambled.
     * WARNING: Takes a very long time to complete
     */
    void medianCollapse();
	/*!
	 * void medianCollapse()
	 * \brief collapse n slices of the cube into a single slice
	 */
	void medianCollapse(unsigned n);
    /*!
     * bool isCube()
     * \brief Is this image a cube?
     */    
    bool isCube(void);
    /*! 
	 * void setSlice(unsigned slice)
	 * \brief sets the current slice
	 */
	void setSlice(unsigned slice);
    /*! 
     * void getSlice()
     * \brief get the current slice
     */
    unsigned getSlice();
	/*! 
	 * unsigned getslices(void)
	 * \brief gets the number of slices
	 */
	unsigned getslices(void);
	/*! 
	 * void readSlice(unsigned slice)
	 * \brief read a slice into memory
	 */
	void readSlice(unsigned slice);
	/*! 
	 * void saveSlice(unsigned slice)
	 * \brief save image data
	 */
	void saveSlice(unsigned slice);
};

/*! 
 * void map(operaFITSCube &, float from, float to)
 * \brief map all pixels from into to.
 * usage: map(image, F_NAN, 0.0);
 */
static inline void map(operaFITSCube &image, float from, float to) {
	float *p = (float *)image.getpixels(); 
	unsigned n = image.getnpixels(); 
	while (n--) {*p=(*p==from?to:*p); p++;}
}

#endif

