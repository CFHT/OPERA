/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaFITSSubImage
 Version: 1.0
 Description: class encapsulates a FITS subimage.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 Contact: opera@cfht.hawaii.edu
 
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
 
#include "fitsio.h"
#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaImage.h"
#include "libraries/operaGeometricShapes.h"	// for Box
#include "libraries/operaImageVector.h"

/*!
 * \brief This class encapsulates the FITS subimage.
 * \file operaFITSSubIMage.cpp
 */

using namespace std;

/*!
 * operaFITSSubIMage
 * \author Doug Teeple
 * \brief This class encapsulates the FITS sub image.
 * \ingroup libraries
 */

/* 
 * \class operaFITSSubImage()
 * \brief Basic operaFITSSubImage constructor.
 * \note SUBIMAGES ARE FLOATS ONLY!!!!!!
 * \return void
 */
operaFITSSubImage::operaFITSSubImage() : 
nx(0),				// x-dimension to be figured out from arguments (ncols) 
ny(0),				// y-dimension to be figured out from arguments (nrows)
npixels(0),			// number of pixels	
pixptr(NULL),		// pixel data values
istemp(false)		// set if this instance is a temp created in an expression
{
}

/*
 * \class operaFITSSubImage(unsigned nx, unsigned ny)
 * \brief Basic operaFITSSubImage constructor that allocates pixels.
 * \param NX - x size
 * \param NY - y size
 */
operaFITSSubImage::operaFITSSubImage(unsigned NX, unsigned NY) : 
nx(0),				// x-dimension to be figured out from arguments (ncols) 
ny(0),				// y-dimension to be figured out from arguments (nrows)
npixels(0),			// number of pixels	
pixptr(NULL),		// pixel data values
istemp(false)		// set if this instance is a temp created in an expression
{
	nx = NX;
	ny = NY;
	npixels = nx * ny;
	pixptr = malloc(npixels*sizeof(float));
}
/* 
 * operaFITSSubImage(string filename, unsigned X, unsigned Y, unsigned NX, unsigned NY) 
 * \brief Creates a sub image by reading a tile from filename.
 * \param filename - the FITS file to read
 * \param X - beginning x location
 * \param Y - beginning y location
 * \param NX - x width
 * \param NY - y width
 */
operaFITSSubImage::operaFITSSubImage(string filename, unsigned X, unsigned Y, unsigned NX, unsigned NY)  : 
nx(0),				// x-dimension to be figured out from arguments (ncols) 
ny(0),				// y-dimension to be figured out from arguments (nrows)
npixels(0),			// number of pixels	
pixptr(NULL),		// pixel data values
istemp(false)		// set if this instance is a temp created in an expression
{
	int status = 0;
	int datatype;
	void *tmpptr;
	fitsfile *fptr;		// FITS file pointer
	long fpixel[2] = {1+X, 1+Y};
	long lpixel[2] = {1+X+NX,1+Y+NY};
	long inc = 1;
	
	nx = NX;
	ny = NY;
	npixels = nx * ny;
	pixptr = malloc(npixels*sizeof(float));
	if (fits_open_diskfile(&fptr, filename.c_str(), READONLY, &status)) {
		throw operaException("operaFITSSubImage: cfitsio error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_get_img_equivtype(fptr, &datatype, &status)) {
		throw operaException("operaFITSSubImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned i = npixels;
	switch ((edatatype)datatype) {
		case tshort: {
			tmpptr = malloc(sizeof(short)*npixels);
			if (!tmpptr) {
				throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_read_subset(fptr, tshort, fpixel, lpixel, &inc, NULL, (short *)tmpptr, NULL, &status)) { // workaround!!!! FIXMERIGHT
				throw operaException("operaFITSSubImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}				
			float *topp = (float *)pixptr;
			short *from = (short *)tmpptr;
			while (i--) {
				*topp++ = *from++;
			}
		}
			break;
		case tushort: {	
			tmpptr = malloc(sizeof(unsigned short)*npixels);
			if (!tmpptr) {
				throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_read_subset(fptr, tushort, fpixel, lpixel, &inc, NULL, (unsigned short *)pixptr, NULL, &status)) {
				throw operaException("operaFITSSubImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}				
			float *topp = (float *)pixptr;
			unsigned short *from = (unsigned short *)tmpptr;
			while (i--) {
				*topp++ = *from++;
			}
		}
			break;
		case tfloat: {
			tmpptr = malloc(sizeof(float)*npixels);
			if (!tmpptr) {
				throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_read_subset(fptr, tfloat, fpixel, lpixel, &inc, NULL, (float *)pixptr, NULL, &status)) {
				throw operaException("operaFITSSubImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}				
			float *topp = (float *)pixptr;
			float *from = (float *)tmpptr;
			while (i--) {
				*topp++ = *from++;
			}
		}
			break;
		case tdouble: {	
			tmpptr = malloc(sizeof(double)*npixels);
			if (!tmpptr) {
				throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_read_subset(fptr, tdouble, fpixel, lpixel, &inc, NULL, (double *)pixptr, NULL, &status)) {
				throw operaException("operaFITSSubImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}				
			float *topp = (float *)pixptr;
			double *from = (double *)tmpptr;
			while (i--) {
				*topp++ = *from++;
			}
		}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;		
			break;
	}
	free(tmpptr);
	fits_close_file(fptr, &status);	// for some reason most cfitsio programs don't check the status...
}
/* 
 * operaFITSSubImage(operaFITSImage &from) 
 * \brief Creates a sub image from an operaFITSImage.
 * \param from - Image from which to clone the sub image
 * \param X - beginning x location
 * \param Y - beginning y location
 * \param NX - x width
 * \param NY - y width
 * \return void
 */
operaFITSSubImage::operaFITSSubImage(operaFITSImage &from, unsigned X, unsigned Y, unsigned NX, unsigned NY) :
pixptr(NULL),		// pixel data values
istemp(false)		// set if this instance is a temp created in an expression
{	
	nx = NX;
	ny = NY;
	x = X;
	y = Y;
	npixels = nx * ny;
	
	pixptr = malloc(from.getelementsize()*npixels); 
	if (!pixptr) {
		throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	switch (from.getdatatype()) {
		case tushort:
			for (unsigned j = 0; j < ny; j++) {
				for (unsigned i = 0; i < nx; i++) {
					setpixel((float)from.getpixelUSHORT(x+i,y+j), i, j); 
				}
			}
			break;
		case tfloat:
			for (unsigned j = 0; j < ny; j++) {
				for (unsigned i = 0; i < nx; i++) {
					setpixel((float)from.getpixel(x+i,y+j), i, j); 
				}
			}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;		
			break;
	}
}

/* 
 * operaFITSSubImage(operaFITSImage &from, ImageIndexVector v) 
 * \brief Creates a sub image from an operaFITSImage.
 * \param from - Image from which to clone the sub image
 * \param v - ImageIndexVector 
 * \note - the subimage is one dimensional, in the x direction
 * \return void
 */
operaFITSSubImage::operaFITSSubImage(operaFITSImage &from, operaImageVector &v) :
pixptr(NULL),		// pixel data values
istemp(false)		// set if this instance is a temp created in an expression
{	
	nx = v.getlength();
	ny = 1;
	x = 0;
	y = 0;
	npixels = nx * ny;
	
	pixptr = malloc(nx*sizeof(float));
	if (!pixptr) {
		throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *fromptr = (float *)from.getpixels();
	float *toptr = (float *)pixptr;
	unsigned i = nx;
	while (--i) {
		*toptr++ = fromptr[v[i]-1];
	}
}

/* 
 * operaFITSSubImage(operaFITSSubImage &from, ImageIndexVector v) 
 * \brief Creates a sub image from an operaFITSImage.
 * \param from - SubImage from which to clone the sub image
 * \param v - ImageIndexVector 
 * \note - the subimage is one dimensional, in the x direction
 * \return void
 */
operaFITSSubImage::operaFITSSubImage(operaFITSSubImage &from, operaImageVector &v) :
pixptr(NULL),		// pixel data values
istemp(false)		// set if this instance is a temp created in an expression
{	
	
	nx = v.getlength();
	ny = 1;
	x = 0;
	y = 0;
	npixels = nx * ny;
	
	pixptr = malloc(nx*sizeof(float));
	if (!pixptr) {
		throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *fromptr = (float *)from.getpixels();
	float *toptr = (float *)pixptr;
	unsigned i = nx;
	while (--i) {
		*toptr++ = fromptr[v[i]-1];
	}
}


/* 
 * operaFITSSubImage(operaFITSImage &from, , Box box) 
 * \brief Creates a sub image from an operaFITSImage.
 * \param from - Image from which to clone the sub image
 * \param Box box
 */
operaFITSSubImage::operaFITSSubImage(operaFITSImage &from, Box &box) :
pixptr(NULL),		// pixel data values
istemp(false)		// set if this instance is a temp created in an expression
{	
	nx = box.getDX();
	ny = box.getDY();
	x = box.getX1();
	y = box.getY1();
	npixels = nx * ny;
	
	pixptr = malloc(from.getelementsize()*npixels); 
	if (!pixptr) {
		throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	switch (from.getdatatype()) {
		case tushort:
			for (unsigned j = 0; j < ny; j++) {
				for (unsigned i = 0; i < nx; i++) {
					setpixel((float)from.getpixelUSHORT(x+i,y+j), i, j); 
				}
			}
			break;
		case tfloat:
			for (unsigned j = 0; j < ny; j++) {
				for (unsigned i = 0; i < nx; i++) {
					setpixel((float)from.getpixel(x+i,y+j), i, j); 
				}
			}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;		
			break;
	}
}

/* 
 * float* operaFITSSubImage::operaFITSSubImageClonePixels()
 * \brief Clone float pixel data.
 * \return pixels*
 */
float* operaFITSSubImage::operaFITSSubImageClonePixels() {
	long npixels = nx * ny;
	long size = sizeof(float)*npixels;
	float *p = (float *)malloc(size); 
	if (!p) {
		throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (pixptr)
		memcpy(p, pixptr, size);
	else
		memset(p, 0, size);
	return p;
}

/* 
 * void operaFITSSubImage::operaFITSSubImageSetData(float* data)
 * \brief set the iamge data pointer to a buffer of data.
 * \param data pointer to the data
 * \return void
 */
void operaFITSSubImage::operaFITSSubImageSetData(float* data) {
	pixptr = (void *)data;
}

/* 
 * operaFITSSubImage::transpose()
 * \brief transpose y for x.
 * \return void
 */
void operaFITSSubImage::transpose() {
	
	float *newpixptr = (float *)malloc(sizeof(float)*npixels);
	if (!newpixptr)
		throw operaException("operaFITSSubImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			((float *)newpixptr)[(ny*i)+j] = (float)getpixel(i,j); 
		}
	}
	unsigned tmp = ny;
	ny = nx;
	nx = tmp;
	if (pixptr)
		free(pixptr);
	pixptr = newpixptr;
}

/* 
 * void *getpixels()
 * \brief get the image array base.
 * \return void *pixels pointer
 */
void *operaFITSSubImage::getpixels() {
	return pixptr;
}
/* 
 * unsigned getnx()
 * \brief get the image array length of x dimension.
 * \return unsigned length of axis x
 */
unsigned operaFITSSubImage::getnx() {
	return nx;
}
/* 
 * unsigned getny()
 * \brief get the image array length of y dimension.
 * \return unsigned length of axis y
 */
unsigned operaFITSSubImage::getny() {
	return ny;
}
/* 
 * unsigned getnpixels()
 * \brief get the image array number of pixels.
 * \return unsigned number of pixels
 */
unsigned operaFITSSubImage::getnpixels() {
	return npixels;
}
