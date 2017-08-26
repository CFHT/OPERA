/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaFITSCube
 Version: 1.0
 Description: class encapsulates a Cube FITS image.
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaImage.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSCube.h"
#include "libraries/operaException.h"

/*! \file operaFITSCube.cpp */

using namespace std;

/*! 
 * operaFITSCube
 * \author Megan Tannock
 * \brief This class encapsulates the Cube FITS image.
 * \ingroup libraries
 */

/*
 * operaFITSCube()
 * \brief Basic operaFITSCube constructor.
 * \return void
 */
operaFITSCube::operaFITSCube() : operaFITSImage()
{
	imageType = FITSCube;
    naxis = 3;
}

/* 
 * \class operaFITSCube
 * \brief construct an in-memory FITSImage object
 * \brief operaFITSCube(unsigned Naxis1, unsigned Naxis2, edatatype Datatype=tushort, unsigned Compression=0)
 * \brief Create an in-memory image of given dimensions.
 * \param Naxis1 - x ccd dimension
 * \param Naxis2 - y ccd dimension
 * \param Datatype optional datatype defaults to tshort
 * \param Compression optional compression, defaults to none
 * \return void
 */
operaFITSCube::operaFITSCube(unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, edatatype Datatype)
{
	mode = READWRITE;
    naxis = 3;
    naxis1 = Naxis1;
    naxis2 = Naxis2;
    naxis3 = Naxis3;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
    npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
	npixels_per_slice = naxis1*naxis2;
	npixels_per_extension = naxis1*naxis2*naxis3;
	compression = cNone;
    datatype = Datatype;
	imageType = FITSCube;
    hdu = 0;	// signals the a header may need to be created...
	isLazy = false;
	
	bitpix = tobitpix(Datatype);
	size_t size = toSize(bitpix, npixels);
	pixptr = malloc(size); 
	memset(pixptr, 0, size);		
	AllExtensions = true;
	AllSlices = true;    
}

/* 
 * \class operaFITSCube
 * \brief create a writeable file image of an in memory FITSImage object
 * \brief operaFITSCube(string Filename, unsigned Naxis1, unsigned Naxis2, edatatype Datatype, unsigned Compression)
 * \brief Constructor to create a new FITS file.
 * \param Filename to create (file is deleted if it exists)
 * \param Naxis1 - dimensions
 * \param Naxis2 - dimensions
 * \param Datatype defaults to tshort
 * \param Compression, defaults to none
 * \throws operaException cfitsio error code
 * \return void
 */
operaFITSCube::operaFITSCube(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, edatatype Datatype, unsigned Compression, bool IsLazy)
{ 
    int status = 0;
    
    filename = Filename;
    compression = Compression;
    mode = READWRITE;
    naxis = 3;
    naxis1 = Naxis1;
    naxis2 = Naxis2;
    naxis3 = Naxis3;
    naxes[0] = naxis1;
    naxes[1] = naxis2;
    naxes[2] = naxis3;
    npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
    npixels_per_slice = naxis1*naxis2;
    npixels_per_extension = naxis1*naxis2*naxis3;
    compression = Compression;
    datatype = Datatype;
	imageType = FITSCube;
    hdu = 0;	// signals the a header may need to be created...
	isLazy = IsLazy;
    
    // remove existing file - cfitsio returns an error if it exists...
    remove(filename.c_str());
    
    bitpix = tobitpix(Datatype);
    
    // Open FITS file for output
    if (fits_create_file(&fptr, filename.c_str(), &status)) {
        throw operaException("operaFITSCube: create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
    }
    if (fits_set_compression_type(fptr, compression, &status)) {
        throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
    }
    size_t size = toSize(bitpix, npixels);
    pixptr = malloc(size); 
    if (!pixptr) {
        throw operaException("operaFITSCube: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
    }
    memset(pixptr, 0, size);
    if (fptr == NULL || hdu == 0) {
        // Create headers (including primary) and images
        if (fits_create_img(fptr, bitpix, 0, naxes, &status)) {
            throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
    }
    if (!isLazy) {
        if (extensions > 0) {
            AllExtensions = true;
            AllSlices = true;
        }
    }    
}

/* 
 * \class operaFITSCube
 * \brief create a FITSIMage object from a FITS file
 * \brief operaFITSCube(string Filename, int mode=READWRITE)
 * \brief Constructor to create a FITSImage from a FITS file.
 * \param Filename
 * \param mode
 * \return void
 */
operaFITSCube::operaFITSCube(string Filename, edatatype Datatype, int mode/*=READWRITE|READONLY*/, unsigned Compression, bool IsLazy) : operaFITSImage(Filename, Datatype, mode, Compression, IsLazy)
{
    // should this be filled in more? eg: Naxis3 = something, naxis = 3, etc
}
/* 
 * operaFITSImage* operaFITSCube::operaFITSCubeClone(operaFITSImage &imageIn, bool ViewOnly, bool AddHeader)
 * \brief Clone a  FITSImage object.
 * \param imageIn - pointer to image to clone
 * \return operaFITSImage* 
 */
operaFITSCube::operaFITSCube(operaFITSCube &imageIn, bool ViewOnly, bool AddHeader) : operaFITSImage() {
    // for cfitsio routines:
	int status = 0;
	int hdutype = ANY_HDU;
    
	naxis = imageIn.naxis;
	naxis1 = imageIn.naxis1;
	naxis2 = imageIn.naxis2; 
	naxis3 = imageIn.naxis3; 
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	extensions = imageIn.extensions;
	npixels = imageIn.npixels;
	npixels_per_slice = imageIn.npixels_per_slice;
	npixels_per_extension = imageIn.npixels_per_extension;
	current_extension = imageIn.current_extension;
	current_slice = imageIn.current_slice;
	current_slice = imageIn.current_slice;
    mode = imageIn.mode;
    	
	fptr = imageIn.fptr;
	isClone = true;				// i.e. do not close fptr
	datatype = imageIn.datatype;
	imageType = FITSCube;
	bitpix = imageIn.bitpix;
	isLazy = imageIn.isLazy;
	viewOnly = ViewOnly;
	AllExtensions = imageIn.AllExtensions;
	AllSlices = imageIn.AllSlices;
	
	// The memcpy may not be needed
	if (ViewOnly) {
		pixptr = imageIn.pixptr;
	} else {
		size_t size = imageIn.getsize();
		pixptr = malloc(size); 
		if (!isLazy) 
			memcpy(pixptr, imageIn.getpixels(), size);
	}
    // add the header values
    if (AddHeader) {
        // Move to the current extension HDU (the current_extension should have been set before calling 
        // this constructor. The primary header for a FITS Cube is extension 1 (extension 0 means no header)
        if ( fits_movabs_hdu(imageIn.fptr, cfitsioextension(0), &hdutype, &status) ) {
            throw operaException("operaFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if ( fits_movrel_hdu(imageIn.fptr, imageIn.current_extension, &hdutype, &status) ) {
            throw operaException("operaFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
        }
	} else {
        hdu = 0; // no header to be copied. This will give a default header in save or saveAs
    }
}

/* 
 * operaFITSCubeSave() 
 * \brief Saves the current image to disk.
 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
 * \throws operaException operaErrorCodeNoFilename
 * \throws operaException cfitsio error code
 * \return void 
 */
void operaFITSCube::operaFITSCubeSave() {
	int status = 0;
	long fnaxes[MAXFITSDIMENSIONS];			// FITS image dimension in file
	fitsfile *newfptr;		// FITS file pointer for the new image
	const long fpixel = 1;
	const int morekeys = 0;
	
	if (mode == READONLY) {
        throw operaException("operaFITSCube:: "+filename+" ", operaErrorCodeChangeREADONLYError, __FILE__, __FUNCTION__, __LINE__);	
    }	

	if (filename.empty()) {
		throw operaException("operaFITSCube: ", operaErrorCodeNoFilename, __FILE__, __FUNCTION__, __LINE__);
    }
	
	// remove existing file - cfitsio returns an error if it exists...
	remove(filename.c_str());
	
	if (fits_create_file(&newfptr, filename.c_str(), &status)) {
		throw operaException("operaFITSCube: cfitsio create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fptr == NULL || hdu == 0) {
		if (fits_create_img(newfptr, bitpix, naxis, naxes, &status)) {
			throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	//
	// check for an image resize, if hdu exists
	//
	if (hdu != 0 && fptr != 0) {
		//	get img size from FITS file
		if (fits_get_img_size(fptr, 2, fnaxes, &status)) {
			throw operaException("operaFITSCube: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}	
		if ((long)naxis1 != fnaxes[0] || (long)naxis2 != fnaxes[1]) {
			if (fits_resize_img(newfptr, bitpix, 2, naxes, &status)) {
				throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	//
	// be sure we have the basic headers for floats if no header existed
	//
	if (hdu == 0) {
		if (datatype == tfloat || datatype == tdouble) {
			float bzero = 0.0, bscale = 1.0;
			if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
				throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (datatype == tfloat) {
				if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
					throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {
				if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
					throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} 
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	if (fits_write_img(newfptr, datatype, fpixel, npixels, pixptr, &status)) {
		throw operaException("operaFITSCube: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	fits_close_file(newfptr, &status);
	newfptr = NULL;
	// note that fptr is left open in case the user wants to save again.
	// It will be closed by a close call if made, or the destructor otherwise
}
/* 
 * operaFITSCubeeSaveAs(string newFilename) 
 * \brief Saves the current image to disk, with the given filename.
 * \param newFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSCube::operaFITSCubeSaveAs(string newFilename) {
	int status = 0;
	fitsfile *newfptr;		// FITS file pointer for the new image
	const long fpixel = 1;
	const int morekeys = 0;
	
	// remove existing file - cfitsio returns an error if it exists...
	remove(newFilename.c_str());
	
	if (fits_create_file(&newfptr, newFilename.c_str(), &status)) {
		throw operaException("operaFITSCube: cfitsio create error: "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fptr == NULL || hdu == 0) {
		if (fits_create_img(newfptr, bitpix, naxis, naxes, &status)) {
			throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	//
	// be sure we have the basic headers for floats if no header existed
	//
	if (hdu == 0) {
		if (datatype == tfloat || datatype == tdouble) {
			float bzero = 0.0, bscale = 1.0;
			if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
				throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (datatype == tfloat) {
				if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
					throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {
				if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
					throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	if (fits_write_img(newfptr, datatype, fpixel, npixels, pixptr, &status)) {
		throw operaException("operaFITSCube: cfitsio error "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	fits_close_file(newfptr, &status);
	newfptr = NULL;
	// note that fptr is left open in case the user wants to save again.
	// It will be closed by a close call if made, or the destructor otherwise
}
/*
 * void medianCollapse()
 * \brief collapse the cube into a single slice
 */
void operaFITSCube::medianCollapse() {
	operaFITSImage &output = (operaFITSImage &)*this;
    const unsigned maxx = output.getXDimension();
    const unsigned maxy = output.getYDimension();
    const unsigned maxz = output.getZDimension();
    
    if (maxz == 1) {
        return;         // no collapse to be done because there is only one slice
    } else {            // if more than one slice, need to collapse slices
		float *pixelStack = new float[maxz];
		for (unsigned y=0; y<maxy; y++) {
            for (unsigned x=0; x<maxx; x++) {
                for (unsigned z=1; z<=maxz; z++) {
					float *base = (float *)output.getpixels() + (z-1) * npixels_per_slice;
                    pixelStack[z] = *(base + (y*naxis1) + x);
                }
                output[y][x] = operaArrayMedianQuick(maxz, (float *)pixelStack);
            }
        }
		naxes[2] = naxis3 = 1;
		current_slice = 1;
		delete[] pixelStack;
    }
}
/*
 * void medianCollapse()
 * \brief collapse n slices of the cube into a single slice
 */
void operaFITSCube::medianCollapse(unsigned n) {
	operaFITSImage &output = (operaFITSImage &)*this;
    const unsigned maxx = output.getXDimension();
    const unsigned maxy = output.getYDimension();
    const unsigned maxz = output.getZDimension();
    
    if (maxz == 1) {
        return;         // no collapse to be done because there is only one slice
    } else {            // if more than one slice, need to collapse slices
		float *pixelStack = new float[maxz];
		for (unsigned y=0; y<maxy; y++) {
            for (unsigned x=0; x<maxx; x++) {
                for (unsigned z=1; z<=n; z++) {
					float *base = (float *)output.getpixels() + (z-1) * npixels_per_slice;
                    pixelStack[z] = *(base + (y*naxis1) + x);
                }
                output[y][x] = operaArrayMedianQuick(maxz, (float *)pixelStack);
            }
        }
		naxes[2] = naxis3 = 1;
		current_slice = 1;
		delete[] pixelStack;
    }
}
/*
 * bool isCube()
 * \brief Is this image a cube?
 * \return bool
 */    
bool operaFITSCube::isCube(void) {
	return getZDimension() > 1;
}
/* 
 * void setSlice(int slice)
 * \brief sets the current slice
 * \param int slice
 */
void operaFITSCube::setSlice(unsigned slice) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (slice > getslices()) {
        throw operaException("operaFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	current_slice = slice;
}
/* 
 * void getSlice()
 * \brief get the current slice
 * \param unsigned slice
 */
unsigned operaFITSCube::getSlice() {
    return current_slice;
}
/* 
 * unsigned getslices(void)
 * \brief gets the number of slices
 */
unsigned operaFITSCube::getslices(void) {
	return naxis3;
}
operaFITSCube::~operaFITSCube() {
}
