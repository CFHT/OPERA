/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
Library name: liboperaFITSSubImage
Class: operaFITSSubImage
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

#ifndef OPERAFITSSUBIMAGE_H
#define OPERAFITSSUBIMAGE_H

#include "libraries/operaEspadonsImage.h"

using namespace std;

class operaFITSImage;
class operaFITSSubImage;
class operaImageVector;
class Box;

/*!
 * operaFITSSubImage class
 * \author Doug Teeple
 * \brief This class encapsulates a FITS subimage.
 * \file operaFITSSubImage.h
 * \ingroup libraries
 */
class operaFITSSubImage {

	friend class operaEspadonsSubImage;
	friend class operaEspadonsImage;
	friend class operaFITSImage;
	friend class operaMultiExtensionFITSImage;

protected:
	unsigned nx;		// x-dimension from arguments (ncols) 
	unsigned ny;		// y-dimension arguments (nrows)
	unsigned npixels;	// number of subimage pixels	
	unsigned x;			// start x, not used after the initial copy
	unsigned y;			// start y, not used after the initial copy	
	void *pixptr;		// pixel base
	bool istemp;		// set if this instance is a temp created in an expression
	
public:
	/*!
	 * \class operaFITSSubImage()
	 * \brief Basic operaFITSSubImage constructor.
	 */
	operaFITSSubImage();														// simply construct a default image
	/*!
	 * \class operaFITSSubImage(unsigned nx, unsigned ny)
	 * \brief Basic operaFITSSubImage constructor that allocates pixels.
	 */
	operaFITSSubImage(unsigned nx, unsigned ny);								// simply construct a default image
	/*! 
	 * operaFITSSubImage(operaFITSImage &from, ImageIndexVector v) 
	 * \brief Creates a sub image from an operaFITSImage.
	 * \param from - Image from which to clone the sub image
	 * \param v - ImageIndexVector 
	 * \note - the subimage is one dimensional, in the x direction
	 * \return void
	 */
	operaFITSSubImage(operaFITSSubImage &from, operaImageVector &v);
	/*! 
	 * operaFITSSubImage(operaFITSSubImage &from, ImageIndexVector v) 
	 * \brief Creates a sub image from an operaFITSSubImage.
	 * \param from - Image from which to clone the sub image
	 * \param v - ImageIndexVector 
	 * \note - the subimage is one dimensional, in the x direction
	 * \return void
	 */
	operaFITSSubImage(operaFITSImage &from, operaImageVector &v);
	
	/*! 
	 * operaFITSSubImage(string filename, unsigned X, unsigned Y, unsigned NX, unsigned NY) 
	 * \brief Creates a sub image by reading a tile from filename.
	 * \param filename - the FITS file to read
	 * \param X - beginning x location
	 * \param Y - beginning y location
	 * \param NX - x width
	 * \param NY - y width
	 */
	operaFITSSubImage(string filename, unsigned X, unsigned Y, unsigned NX, unsigned NY);

	/*! 
	 * operaFITSSubImage(operaFITSImage &from, unsigned X, unsigned Y, unsigned NX, unsigned NY) 
	 * \brief Creates a sub image from an operaFITSImage.
	 * \param from - Image from which to clone the sub image
	 * \param X - beginning x location
	 * \param Y - beginning y location
	 * \param NX - x width
	 * \param NY - y width
	 */
	operaFITSSubImage(operaFITSImage &from, unsigned X, unsigned Y, unsigned NX, unsigned NY);
	
	/*! 
	 * operaFITSSubImage(operaFITSImage &from, , Box box) 
	 * \brief Creates a sub image from an operaFITSImage.
	 * \param from - Image from which to clone the sub image
	 * \param Box box
	 */
	operaFITSSubImage(operaFITSImage &from, Box &box);
	
	/*! 
	 * \brief operator []
	 * \brief assignment operator.
	 * \param i - index into row
	 * \note usage: myfitssubimage[y][x] = 0.0; sets the image pixel value at coords x, y to zero
	 * \return float pointer to row 
	 */	
	float *operator[](unsigned i) { return (float *)pixptr+(nx*i);};
	/*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b; copies the pixel values from b to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator=(operaFITSSubImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ = *bp++;
		if (b.istemp) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = 0.0; copies the float value to every pixel in = a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator=(float f) {
		float *p = (float *)pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ = f;
		return *this;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a += operaFITSSubImage b; adds the pixel values from b to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator+=(operaFITSSubImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ += *bp++;
		if (b.istemp) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a += 100.0; adds the float value to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator+=(float f) {
		float *p = (float *)pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ += f;
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a -= operaFITSSubImage b; subtracts the pixel values from b to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator-=(operaFITSSubImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ -= *bp++;
		if (b.istemp) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a -= 100.0; subtracts the float value from a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator-=(float f) {
		float *p = (float *)pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ -= f;
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a *= operaFITSSubImage b; multiplies the pixel values from b to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator*=(operaFITSSubImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ *= *bp++;
		if (b.istemp) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a *= 100.0; multiplies the float value times a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator*=(float f) {
		float *p = (float *)pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ *= f;
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a /= operaFITSSubImage b; divides the pixel values from b into a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator/=(operaFITSSubImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ /= (*bp!=0.0?*bp++:(bp++,1.0));
		if (b.istemp) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a /= 100.0; divides the float value from a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator/=(float f) {
		float *p = (float *)pixptr; 
		unsigned n = npixels; 
		while (n--) *p++ /= (f!=0.0?f:1.0);
		return *this;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b * operaFITSSubImage c; multiplies the pixel values  b * c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator*(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ * *bp++;
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b * 10.0; multiplies the pixel values  b * 10.0 and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator*(float f) {
		operaFITSSubImage *t = new operaFITSSubImage(this->nx, this->ny);
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
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b / operaFITSSubImage c; divides the pixel values  b / c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator/(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ / (*bp!=0.0?*bp++:(bp++,1.0));
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b / 100.0; divides the pixel values  b / 100.0 and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator/(float f) {
		operaFITSSubImage *t = new operaFITSSubImage(this->nx, this->ny);
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
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b + operaFITSSubImage c; adds the pixel values  b + c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator+(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ + *bp++;
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b + 100.0; adds the pixel values  b + 100.0 and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator+(float f) {
		operaFITSSubImage *t = new operaFITSSubImage(this->nx, this->ny);
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
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b - operaFITSSubImage c; subtracts the pixel values  b - c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator-(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = *p++ - *bp++;
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b - 100.0; subtracts the pixel values  b - 100.0 and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator-(float f) {
		operaFITSSubImage *t = new operaFITSSubImage(this->nx, this->ny);
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
	 * \note usage: operaFITSSubImage a = !operaFITSSubImage b; inverts the pixel values of b and assigns to a (0.0 becomes 1.0, non-zero becomes 0.0)
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator!() {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
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
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b > operaFITSSubImage c; creates a mask of 1.0 if b > c or 0.0 if b <= c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator>(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++>*bp++?1.0:0.0);
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator >
	 * \brief greater than operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b > 100.0; creates a mask of 1.0 if b > 100.0 or 0.0 if b <= 100.0 and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator>(float f) {
		operaFITSSubImage *t = new operaFITSSubImage(this->nx, this->ny);
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
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b >= operaFITSSubImage c; creates a mask of 1.0 if b >= c or 0.0 if b < c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator>=(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++>=*bp++?1.0:0.0);
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator >=
	 * \brief greater than or equal operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b >= 100.0; creates a mask of 1.0 if b >= 100.0 or 0.0 if b < 100.0 and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator>=(float f) {
		operaFITSSubImage *t = new operaFITSSubImage(this->nx, this->ny);
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
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b < operaFITSSubImage c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator<(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++<*bp++?1.0:0.0);
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator <
	 * \brief less than operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b < operaFITSSubImage c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator<(float f) {
		operaFITSSubImage *t = new operaFITSSubImage(this->nx, this->ny);
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
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b <= operaFITSSubImage c; creates a mask of 1.0 if b <= c or 0.0 if b > c and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator<=(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p<=*bp?(bp++,*p++):(p++,*bp++));
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator <=
	 * \brief less than or equal operator.
	 * \param f - float
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b <= 100.0; creates a mask of 1.0 if b <= 100.0 or 0.0 if b > 100.0 and assigns to a
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator<=(float f) {
		operaFITSSubImage *t = new operaFITSSubImage(this->nx, this->ny);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++<=f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator &&
	 * \brief less than operator.
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b && operaFITSSubImage c; creates a mask of 1.0 where b && c != 0.0
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator&&(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++!=0.0&&*bp++!=0.0?1.0:0.0);
		if (b.istemp) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator ||
	 * \brief less than operator.
	 * \param b - operaFITSSubImage&
	 * \note usage: operaFITSSubImage a = operaFITSSubImage b && operaFITSSubImage c; creates a mask of 1.0 where b && c != 0.0
	 * \return operaFITSSubImage&
	 */	
	operaFITSSubImage& operator||(operaFITSSubImage& b) {
		float *tp;
		operaFITSSubImage *t;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSSubImage(this->nx, this->ny);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned n = npixels; 
		while (n--) *tp++ = (*p++!=0.0||*bp++!=0.0?1.0:0.0);
		if (b.istemp) delete &b;
		return *t;
	};
	
	/*! 
	 * float* operaFITSSubImage::operaFITSSubImageClonePixels()
	 * \brief Clone float pixel data.
	 * \return pixels*
	 */
	float* operaFITSSubImageClonePixels();
	/*! 
	 * void operaFITSSubImage::operaFITSSubImageSetData(data)
	 * \brief set the iamge data pointer to a buffer of data.
	 * \param data pointer to the data
	 * \return void
	 */
	void operaFITSSubImageSetData(float* data);									// set pixptr to a data buffer
	
	/*! 
	 * operaFITSSubImage::getpixel(unsigned x, unsigned y)
	 * \brief get a pixel value at coordinates x,y.
	 * \param x - the x coordinate
	 * \param y - the y coordinate
	 * \return unsigned short
	 */
	inline float getpixel(unsigned x, unsigned y) {return ((float *)pixptr)[(nx*y)+x];};
	
	/*! 
	 * operaFITSSubImage::setpixel(float value, unsigned x, unsigned y)
	 * \brief get a pixel value at coordinates x,y.
	 * \param value - the float pixel value
	 * \param x - the x coordinate
	 * \param y - the y coordinate
	 * \return unsigned short
	 */
	inline void setpixel(float value, unsigned x, unsigned y) {((float *)pixptr)[(nx*y)+x] = value;};
	
	/*! 
	 * operaFITSSubImage::transpose()
	 * \brief transpose y for x.
	 * \return void
	 */
	void transpose();
	
	/* getters and setters */

	/*! 
	 * void *getpixels()
	 * \brief get the image array base.
	 * \return void *pixels pointer
	 */
	void *getpixels();
	/*! 
	 * unsigned getnx()
	 * \brief get the image array length of y dimension.
	 * \return unsigned length of axis x
	 */
	unsigned getnx();
	/*! 
	 * unsigned getnx()
	 * \brief get the image array length of y dimension.
	 * \return unsigned length of axis y
	 */
	unsigned getny();
	/*! 
	 * unsigned getnpixels()
	 * \brief get the image array number of pixels.
	 * \return unsigned number of pixels
	 */
	unsigned getnpixels();
	/*! 
	 * bool getIstemp()
	 * \brief return true if this is a temp.
	 * \return unsigned number of pixels
	 */
	bool getIstemp() {return istemp;};
};

#endif
