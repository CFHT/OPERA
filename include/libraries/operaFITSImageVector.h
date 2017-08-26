#ifndef OPERAFITSIMAGEVECTOR_H
#define OPERAFITSIMAGEVECTOR_H
/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaFITSImageVector   
 Class: operaFITSImageVector 
 Version: 1.0  
 Author(s): CFHT OPERA team 
 Affiliation: Canada France Hawaii Telescope  
 Location: Hawaii USA  
 Date: Aug/2011 
 
 Copyright (C) 2012  Opera Pipeline team, Canada France Hawaii Telescope
 
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

 /*! 
 * operaFITSImageVector
 * \author Doug Teeple
 * \brief A template class for FITS Image Vectors.
 * \file operaFITSImageVector.h
 * \ingroup libraries
 */

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"

/*!
 * operaFITSImageVector class
 * \author Doug Teeple
 * \brief A template class for FITS Image Vectors.
 * \note ONLY works with operaFITSImage at the moment.
 * \note All images must have the same dimensionality.
 * \file operaFITSImageVector.h
 * \package operaFITSImageVector
 * \ingroup libraries
 */
#define MAXFITSVECTORIMAGES 100

template <class FITSImage>
class operaFITSImageVector {
    
private:
    unsigned imagecount;
	FITSImage *images[MAXFITSVECTORIMAGES];
	float etimes[MAXFITSVECTORIMAGES];
	edatatype datatype;
	int mode;
	unsigned compression;
	bool isLazy;
	unsigned long naxis1;
	unsigned long naxis2;
	
public:
	/*!
	 * \sa operaFITSImageVector(string Filenames[], bool IsLazy = false)
	 * \brief Basic operaFITSImageVector constructor from FITS files.
	 * \return void
	 */
	operaFITSImageVector(string Filenames[], bool IsLazy = false);

	/*!
	 * \sa addImage(string Filename)
	 * \brief add an image to the image list.
	 * \return void
	 */
	void addImage(string Filename);
	
	/*!
	 * \sa FITSImage *getImage(unsigned index);
	 * \brief get an image * at index or NULL.
	 * \return FITSImage *
	 */
	FITSImage *getImage(unsigned index);
	
	/*!
	 * \sa unsigned getCount(void);
	 * \brief return the count of input images.
	 * \return unsigned
	 */
	unsigned getCount(void);
	
	/*!
	 * \sa ufloat getMeanEtime(void);
	 * \brief return the mean etime of all input images.
	 * \return float
	 */
	float getMeanEtime(void);
	
	/*!
	 * \sa void removebias(float bias)
	 * \brief remove the bias from the input images.
	 * \return void
	 */
	void removebias(float bias);
	
	/*!
	 * \sa FITSImage *median(void)
	 * \brief return a pointer to a FITSImage containing the median of the input images.
	 * \return FITSImage pointer
	 */
	FITSImage *median(void);
	
	/*!
	 * \sa FITSImage *median(void)
	 * \brief return a pointer to a FITSImage containing the stack of the input images.
	 * \return FITSImage pointer
	 */
	FITSImage *stack(void);
	/*!
	 * \sa FITSImage *mean(void)
	 * \brief return a pointer to a FITSImage containing the mean of the input images.
	 * \return FITSImage pointer
	 */
	FITSImage *mean(void);
	/*!
	 * \sa FITSImage *meanunsaturated(unsigned long saturation)
	 * \brief return a pointer to a FITSImage containing the mean of the input images
	 * \brief neglecting unsaturated pixels.
	 * \return FITSImage pointer
	 */
	FITSImage *meanunsaturated(unsigned long saturation);
	
	/*!
	 * \sa FITSImage *meanunsaturatedbyetime(unsigned long saturation)
	 * \brief return a pointer to a FITSImage containing the mean of the input images
	 * \brief neglecting unsaturated pixels and accounting for etime.
	 * \return FITSImage pointer
	 */
	FITSImage *meanunsaturatedbyetime(unsigned long saturation);
																	   
	/*!
	 * \sa void close(void)
	 * \brief close the input images.
	 * \return void
	 */
	void close(void);
};

template <class FITSImage>
operaFITSImageVector<FITSImage>::operaFITSImageVector(string Filenames[], bool IsLazy) :
imagecount(0),
datatype(tfloat),
mode(READONLY),
compression(cNone),
isLazy(IsLazy),
naxis1(0),
naxis2(0)
{
	while (!Filenames[imagecount].empty()) {
		images[imagecount] = new FITSImage(Filenames[imagecount], datatype, mode, compression, isLazy);
		if ((naxis1 != 0 && naxis1 != images[imagecount]->getnaxis1()) ||
			(naxis2 != 0 && naxis2 != images[imagecount]->getnaxis2())) {
			throw operaException("operaFITSImageVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		} else {
			naxis1 = images[imagecount]->getnaxis1();
			naxis2 = images[imagecount]->getnaxis2();
			etimes[imagecount] = images[imagecount]->operaFITSGetFloatHeaderValue("EXPTIME");
		}
		imagecount++;
	}
}

template <class FITSImage>
void operaFITSImageVector<FITSImage>::addImage(string Filename)
{
	images[imagecount] = new FITSImage(Filename, datatype, mode, compression, isLazy);
	etimes[imagecount] = images[imagecount]->operaFITSGetFloatHeaderValue("EXPTIME");
	imagecount++;
}

template <class FITSImage>
FITSImage *operaFITSImageVector<FITSImage>::getImage(unsigned index)
{
	if (index < imagecount) {
		return images[index];
	} else {
		return NULL;
	}
}

template <class FITSImage>
unsigned operaFITSImageVector<FITSImage>::getCount(void)
{
	return imagecount;
}

template <class FITSImage>
float operaFITSImageVector<FITSImage>::getMeanEtime(void)
{
	float meanetime = 0.0;
	if (imagecount > 0) {
		unsigned i = imagecount;
		while (i--) {
			meanetime += etimes[i];
		}
		meanetime /= imagecount;
	}
	return meanetime;
}

template <class FITSImage>
FITSImage *operaFITSImageVector<FITSImage>::median(void)
{
	FITSImage *imagemedian = NULL;
	if (imagecount > 0) {
		imagemedian = new FITSImage(naxis1, naxis2, datatype);
		float *pixelStack = new float[imagecount];
		for (unsigned y=0; y<naxis2; y++) {
			for (unsigned x=0; x<naxis1; x++) {
				for (unsigned i=0; i<imagecount; i++) {
					FITSImage *image = (FITSImage *)images[i];
					pixelStack[i] = (*image)[y][x];
				}
				(*imagemedian)[y][x] = operaArrayMedianQuick(imagecount, (float *)pixelStack);
			}
		}
		delete[] pixelStack;
	}
	return imagemedian;
}

template <class FITSImage>
void operaFITSImageVector<FITSImage>::removebias(float bias)
{
	unsigned i = imagecount;
	while (i--) {
		*images[i] -= bias;
	}
}

template <class FITSImage>
FITSImage *operaFITSImageVector<FITSImage>::stack(void)
{
	FITSImage *imagestack = NULL;
	if (imagecount > 0) {
		imagestack = new FITSImage(naxis1, naxis2, datatype);
		*imagestack = (float)0.0;
		unsigned i = imagecount;
		while (i--) {
			*imagestack += *images[i];
		}		
	}
	return imagestack;
}

template <class FITSImage>
FITSImage *operaFITSImageVector<FITSImage>::mean(void)
{
	FITSImage *imagestack = NULL;
	if (imagecount > 0) {
		imagestack = new FITSImage(naxis1, naxis2, datatype);
		*imagestack = (float)0.0;
		unsigned i = imagecount;
		while (i--) {
			*imagestack += *images[i];
		}
		*imagestack /= imagecount;
	}
	return imagestack;
}

template <class FITSImage>
FITSImage *operaFITSImageVector<FITSImage>::meanunsaturated(unsigned long saturation)
{
	FITSImage *imagemeanunsaturated = NULL;
	if (imagecount > 0) {
		FITSImage divisors(naxis1, naxis2, datatype);
		divisors = (float)1.0;
		imagemeanunsaturated = new FITSImage(naxis1, naxis2, datatype);
		*imagemeanunsaturated = (float)0.0;
		for (unsigned y=0; y<naxis2; y++) {
			for (unsigned x=0; x<naxis1; x++) {
				for (unsigned i=0; i<imagecount; i++) {
					FITSImage *image = (FITSImage *)images[i];
					float pixelvalue = (*image)[y][x];
					if (pixelvalue < (float)saturation) {
						(*imagemeanunsaturated)[y][x] += pixelvalue;
						divisors[y][x]++;
					}
				}
			}
		}
		*imagemeanunsaturated /= divisors;
	}
	return imagemeanunsaturated;
}

template <class FITSImage>
FITSImage *operaFITSImageVector<FITSImage>::meanunsaturatedbyetime(unsigned long saturation)
{
	FITSImage *imagemeanunsaturated = NULL;
	if (imagecount > 0) {
		FITSImage divisors(naxis1, naxis2, datatype);
		divisors = (float)1.0;
		imagemeanunsaturated = new FITSImage(naxis1, naxis2, datatype);
		*imagemeanunsaturated = (float)0.0;
		for (unsigned y=0; y<naxis2; y++) {
			for (unsigned x=0; x<naxis1; x++) {
				for (unsigned i=0; i<imagecount; i++) {
					FITSImage *image = (FITSImage *)images[i];
					float pixelvalue = (*image)[y][x];
					if (pixelvalue < (float)saturation) {
						(*imagemeanunsaturated)[y][x] += pixelvalue / etimes[i];
						divisors[y][x]++;
					}
				}
			}
		}
		*imagemeanunsaturated /= divisors;
		*imagemeanunsaturated *= getMeanEtime();
	}
	return imagemeanunsaturated;
}

template <class FITSImage>
void operaFITSImageVector<FITSImage>::close(void)
{
	while (imagecount--) {
		images[imagecount]->operaFITSImageClose();
	}
}

#endif

