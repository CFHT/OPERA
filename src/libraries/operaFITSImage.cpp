/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaFITSIMage 
 Version: 1.0
 Description: class encapsulates a FITS image..
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

#include <fstream>

#include "fitsio.h"
#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaImage.h"
#include "libraries/operaLib.h"					// trimFITSKeyword
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaImageVector.h"
#include "libraries/operaGeometricShapes.h"		// Box
#include "libraries/operaException.h"

using namespace std;

/*!
 * operaFITSImage
 * \author Doug Teeple
 * \brief This class encapsulates the FITS image.
 * \file operaFITSIMage.cpp
 * \ingroup libraries
 */

/* 
 * \class operaFITSImage()
 * \brief Basic operaFITSImage constructor.
 * \return void
 */
operaFITSImage::operaFITSImage(bool IsLazy) : 
fptr(NULL),			// FITS file pointer
bitpix(ushort_img),	// BITPIX keyword value (BYTE_IMG, SHORT_IMG, USHORT_IMG, FLOAT_IMG, DOUBLE_IMG)
hdu(0),				// current active hdu
nhdus(1),			// nummber of hdus
naxis(2),			// FITS image dimension
naxis1(0),			// x-dimension to be figured out from NAXIS1 (ncols) 
naxis2(0),			// y-dimension to be figured out from NAXIS2 (nrows)
naxis3(1),			// z-dimension to be figured out from NAXIS3 (nslices)
npixels(0),			// number of pixels	
npixels_per_slice(0),	// total number of pixels in all slices and extensions
npixels_per_extension(0),	// number of ccd pixels per extension
compression(0),		// no compression
datatype(tushort),	// (TSHORT, TUSHORT, TFLOAT, TDOUBLE)		
istemp(false),		// set if this instance is a temp created in an expression
mode(0),			// READWRITE / READONLY
pixptr(NULL),		// pixel data values
current_extension(1),// current active extension
current_slice(1),	// current active slice
extensions(0),		// number of extensions
isLazy(false),		// Lazy read
viewOnly(false),	// Is this a view of somebody else's pixels? BEWARE of deletion!
isClone(false),		// is this a clone of somebody else's fptr? If so do not close!
AllExtensions(false), // all extensions in memory
AllSlices(true),	// all slices in memory
imageType(FITS),	// Kind of image
super(NULL)			// the super class that created this instance
{
	isLazy = IsLazy;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	setHasBeenRead(0, false);
	setHasBeenWritten(0, false);
}

/*
 * \class operaFITSImage(string Filename)
 * \brief Basic operaFITSImage constructor from a FITS file which retains the original datatype.
 * \return void
 */
operaFITSImage::operaFITSImage(string Filename, bool IsLazy) :
fptr(NULL),			// FITS file pointer
bitpix(ushort_img),	// BITPIX keyword value (BYTE_IMG, SHORT_IMG, USHORT_IMG, FLOAT_IMG, DOUBLE_IMG)
hdu(0),				// current active extension
nhdus(1),			// nummber of hdus
naxis(2),			// FITS image dimension
naxis1(0),			// x-dimension to be figured out from NAXIS1 (ncols) 
naxis2(0),			// y-dimension to be figured out from NAXIS2 (nrows)
naxis3(1),			// z-dimension to be figured out from NAXIS3 (nslices)
npixels(0),			// number of pixels	
npixels_per_slice(0),	// total number of pixels in all slices and extensions
npixels_per_extension(0),	// number of ccd pixels per extension
compression(0),		// no compression
datatype(tushort),	// (TSHORT, TUSHORT, TFLOAT, TDOUBLE)		
istemp(false),		// set if this instance is a temp created in an expression
mode(READONLY),		// READWRITE / READONLY
pixptr(NULL),		// pixel data values
varptr(NULL),		// variances
current_extension(1),// current active extension
current_slice(1),	// current active slice
extensions(0),		// number of extensions
isLazy(false),		// Lazy read
viewOnly(false),	// Is this a view of somebody else's pixels? BEWARE of deletion!
isClone(false),		// is this a clone of somebody else's fptr? If so do not close!
AllExtensions(false), // all extensions in memory
AllSlices(true),	// all slices in memory
imageType(FITS),	// Kind of image
super(NULL)			// the super class that created this instance
{
	isLazy = IsLazy;
	filename = Filename;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	setHasBeenRead(0, false);
	setHasBeenWritten(0, false);
	ifstream ifile(filename.c_str());
	if (ifile.good()) {
		openFITSfile(Filename, READONLY);
		readFITSHeaderInfo();
		if (!isLazy) {
			AllExtensions = true;
			AllSlices = true;
			readFITSArray();
		}
	} else if (mode == READONLY) {
		throw operaException("operaFITSImage: ", operaErrorCodeFileDoesNotExistError, __FILE__, __FUNCTION__, __LINE__);	
	}
}

/* 
 * \class operaFITSImage
 * \brief construct an in-memory FITSImage object
 * \brief operaFITSImage(unsigned Naxis1, unsigned Naxis2, edatatype Datatype=tushort, unsigned Compression=0)
 * \brief Create an in-memory image of given dimensions.
 * \param Naxis1 - x ccd dimension
 * \param Naxis2 - y ccd dimension
 * \param Datatype optional datatype defaults to tshort
 * \param Compression optional compression, defaults to none
 * \return void
 */
operaFITSImage::operaFITSImage(unsigned Naxis1, unsigned Naxis2, edatatype Datatype)  : 
fptr(NULL),			// FITS file pointer
bitpix(ushort_img),	// BITPIX keyword value (BYTE_IMG, SHORT_IMG, USHORT_IMG, FLOAT_IMG, DOUBLE_IMG)
bzero(0.0),			// bzero
bscale(1.0),		// bscale
hdu(0),				// current active extension
nhdus(1),			// nummber of hdus
naxis(2),			// FITS image dimension
naxis1(0),			// x-dimension to be figured out from NAXIS1 (ncols) 
naxis2(0),			// y-dimension to be figured out from NAXIS2 (nrows)
naxis3(1),			// z-dimension to be figured out from NAXIS3 (nslices)
npixels(0),			// number of pixels	
npixels_per_slice(0),	// total number of pixels in all slices and extensions
npixels_per_extension(0),	// number of ccd pixels per extension
compression(0),		// no compression
datatype(tushort),	// (TSHORT, TUSHORT, TFLOAT, TDOUBLE)		
istemp(false),		// set if this instance is a temp created in an expression
mode(READWRITE),	// READWRITE / READONLY
pixptr(NULL),		// pixel data values
varptr(NULL),		// variances
current_extension(1),// current active extension
current_slice(1),	// current active slice
extensions(0),		// number of extensions
isLazy(false),		// Lazy read
viewOnly(false),	// Is this a view of somebody else's pixels? BEWARE of deletion!
isClone(false),		// is this a clone of somebody else's fptr? If so do not close!
AllExtensions(true), // all extensions in memory
AllSlices(true),	// all slices in memory
imageType(FITS),	// Kind of image
super(NULL)			// the super class that created this instance
{
	naxis1 = Naxis1;
	naxis2 = Naxis2;
	npixels = naxis1*naxis2;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	setHasBeenRead(0, false);
	setHasBeenWritten(0, false);
	datatype = Datatype;
	if (datatype == tushort) {
		bzero = 32768.0;
	}
	bitpix = tobitpix(Datatype);
	size_t size = toSize(bitpix, npixels);
	pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
	memset(pixptr, 0, size);		
}

operaFITSImage::operaFITSImage(string Filename, edatatype Datatype, int Mode/*=READWRITE|READONLY*/, unsigned Compression, bool IsLazy) :
fptr(NULL),			// FITS file pointer
bitpix(ushort_img),	// BITPIX keyword value (BYTE_IMG, SHORT_IMG, USHORT_IMG, FLOAT_IMG, DOUBLE_IMG)
bzero(0.0),			// bzero
bscale(1.0),		// bscale
hdu(0),				// current active extension
nhdus(1),			// nummber of hdus
naxis(2),			// FITS image dimension
naxis1(0),			// x-dimension to be figured out from NAXIS1 (ncols) 
naxis2(0),			// y-dimension to be figured out from NAXIS2 (nrows)
naxis3(0),			// z-dimension to be figured out from NAXIS3 (nslices)
npixels(0),			// number of pixels	
npixels_per_slice(0),	// total number of pixels in all slices and extensions
npixels_per_extension(0),	// number of ccd pixels per extension
compression(0),		// no compression
datatype(tushort),	// (TSHORT, TUSHORT, TFLOAT, TDOUBLE)		
istemp(false),		// set if this instance is a temp created in an expression
mode(0),			// READWRITE / READONLY
pixptr(NULL),		// pixel data values
varptr(NULL),		// variances
current_extension(1),// current active extension
current_slice(1),	// current active slice
extensions(0),		// number of extensions
isLazy(false),		// Lazy read
viewOnly(false),	// Is this a view of somebody else's pixels? BEWARE of deletion!
isClone(false),		// is this a clone of somebody else's fptr? If so do not close!
AllExtensions(false), // all extensions in memory
AllSlices(true),	// all slices in memory
imageType(FITS),	// Kind of image
super(NULL)			// the super class that created this instance
{
    isLazy = IsLazy;
	setHasBeenRead(0, false);
	setHasBeenWritten(0, false);
	filename = Filename;
	compression = Compression;
	mode = Mode;
	ifstream ifile(filename.c_str());
	if (ifile.good()||true) {
		openFITSfile(Filename, mode);
		readFITSHeaderInfo();
		if (isLazy) {
			npixels = npixels_per_extension = naxis1*naxis2*naxis3;
			bitpix = tobitpix(Datatype);
			size_t size = toSize(bitpix, npixels);
			pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
			memset(pixptr, 0, size);		
		} else {
			AllExtensions = true;
			AllSlices = true;
			readFITSArray();
			if (Datatype != datatype) {
				operaFITSImageConvertImageInPlace(datatype, Datatype);
				if (Datatype == tfloat) {
					bzero = 0.0;
					bscale = 1.0;
				}
			}
		}
		datatype = Datatype;
		bitpix = tobitpix(Datatype);
	} else if (mode == READONLY) {
		throw operaException("operaFITSImage: ", operaErrorCodeFileDoesNotExistError, __FILE__, __FUNCTION__, __LINE__);	
	}
}
/* 
 * \class operaFITSImage
 * \brief create a writeable file image of an in memory FITSImage object
 * \brief operaFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, edatatype Datatype, unsigned Compression)
 * \brief Constructor to create a new FITS file.
 * \param Filename to create (file is deleted if it exists)
 * \param Naxis1 - dimensions
 * \param Naxis2 - dimensions
 * \param Datatype defaults to tshort
 * \param Compression, defaults to none
 * \throws operaException cfitsio error code
 * \return void
 */
operaFITSImage::operaFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, edatatype Datatype, unsigned Compression, bool IsLazy) :
fptr(NULL),			// FITS file pointer
bitpix(ushort_img),	// BITPIX keyword value (BYTE_IMG, SHORT_IMG, USHORT_IMG, FLOAT_IMG, DOUBLE_IMG)
bzero(0.0),			// bzero
bscale(1.0),		// bscale
hdu(0),				// current active extension
nhdus(1),			// nummber of hdus
naxis(2),			// FITS image dimension
naxis1(0),			// x-dimension to be figured out from NAXIS1 (ncols) 
naxis2(0),			// y-dimension to be figured out from NAXIS2 (nrows)
naxis3(0),			// z-dimension to be figured out from NAXIS3 (nslices)
npixels(0),			// number of pixels	
npixels_per_slice(0),	// total number of pixels in all slices and extensions
npixels_per_extension(0),	// number of ccd pixels per extension
compression(0),		// no compression
datatype(tushort),	// (TSHORT, TUSHORT, TFLOAT, TDOUBLE)		
istemp(false),		// set if this instance is a temp created in an expression
mode(READWRITE),	// READWRITE / READONLY
pixptr(NULL),		// pixel data values
varptr(NULL),		// variances
current_extension(1),// current active extension
current_slice(1),	// current active slice
extensions(0),		// number of extensions
isLazy(false),		// Lazy read
viewOnly(false),	// Is this a view of somebody else's pixels? BEWARE of deletion!
isClone(false),		// is this a clone of somebody else's fptr? If so do not close!
AllExtensions(true), // all extensions in memory
AllSlices(true),	// all slices in memory
imageType(FITS),	// Kind of image
super(NULL)			// the super class that created this instance
{
	isLazy = IsLazy;
	int status = 0;
	setHasBeenRead(0, false);
	setHasBeenWritten(0, false);
	filename = Filename;
	compression = Compression;
	naxis = 2;
    naxis1 = Naxis1;
    naxis2 = Naxis2;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
    npixels = Naxis1*Naxis2;
	compression = Compression;
    datatype = Datatype;
	if (datatype == tushort) {
		bzero = 32768.0;
	}
 	bitpix = tobitpix(Datatype);
    hdu = 0;	// signals the a header may need to be created...
    
	// remove existing file - cfitsio returns an error if it exists...
	remove(filename.c_str());
	
	// Open FITS file for output
    if (fits_create_file(&fptr, filename.c_str(), &status)) {
		throw operaException("operaFITSImage: create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
	}
	if (fits_set_compression_type(fptr, compression, &status)) {
		throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	size_t size = toSize(bitpix, npixels);
	pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
	if (!pixptr) {
		throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	memset(pixptr, 0, size);
}

/* 
 * operaFITSImage(operaFITSImage &imageIn)
 * \brief Clone a FITSImage object.
 * \param imageIn - pointer to image to clone
 * \Note - typical usage does not require the memcpy, as the clone is used as a temp.
 * So if this is too slow, the memcpy can be removed.
 * \return operaFITSImage*
 */
operaFITSImage::operaFITSImage(operaFITSImage &imageIn, bool ViewOnly, bool AddHeader) :
fptr(NULL),			// FITS file pointer
bitpix(ushort_img),	// BITPIX keyword value (BYTE_IMG, SHORT_IMG, USHORT_IMG, FLOAT_IMG, DOUBLE_IMG)
bzero(0.0),			// bzero
bscale(1.0),		// bscale
hdu(0),				// current active extension
nhdus(1),			// nummber of hdus
naxis(2),			// FITS image dimension
naxis1(0),			// x-dimension to be figured out from NAXIS1 (ncols) 
naxis2(0),			// y-dimension to be figured out from NAXIS2 (nrows)
naxis3(0),			// z-dimension to be figured out from NAXIS3 (nslices)
npixels(0),			// number of pixels	
npixels_per_slice(0),	// total number of pixels in all slices and extensions
npixels_per_extension(0),	// number of ccd pixels per extension
compression(0),		// no compression
datatype(tushort),	// (TSHORT, TUSHORT, TFLOAT, TDOUBLE)		
istemp(false),		// set if this instance is a temp created in an expression
mode(READWRITE),	// READWRITE / READONLY
pixptr(NULL),		// pixel data values
varptr(NULL),		// variances
current_extension(1),// current active extension
current_slice(1),	// current active slice
extensions(0),		// number of extensions
isLazy(false),		// Lazy read
viewOnly(false),	// Is this a view of somebody else's pixels? BEWARE of deletion!
isClone(false),		// is this a clone of somebody else's fptr? If so do not close!
AllExtensions(true), // all extensions in memory
AllSlices(true),	// all slices in memory
imageType(FITS),	// Kind of image
super(NULL)			// the super class that created this instance
{
	filename = imageIn.filename;
	hdu = imageIn.hdu;	// this is wrong?
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
    mode = imageIn.mode;
    imageType = imageIn.imageType;
	super = &imageIn;
	
	fptr = imageIn.fptr;
	bzero = imageIn.bzero;
	bscale = imageIn.bscale;
	datatype = imageIn.datatype;
	bitpix = imageIn.bitpix;
	isLazy = imageIn.isLazy;
    isClone = true;				// i.e.do not close fptr
	viewOnly = ViewOnly;
	AllExtensions = imageIn.AllExtensions;
	if (AllExtensions) {
		//current_extension = 1; -- DT Aug 2013 - this would undo setting the extension...
	}
	AllSlices = imageIn.AllSlices;
	if (AllSlices) {
		//current_slice = 1; -- DT Aug 2013 - this would undo setting the slice...
	}
	
	// The memcpy may not be needed
	if (ViewOnly && imageIn.pixptr != NULL) {
		pixptr = imageIn.pixptr;
	} else {
		size_t size = getsize();
        pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
		if (imageIn.getpixels()) {
			memcpy(pixptr, (void *)((float *)imageIn.pixptr+(imageIn.current_extension-1)*naxis1*naxis2*naxis3+(imageIn.current_slice-1)*naxis1*naxis2), size);
		}
		if (imageIn.pixptr == NULL) {
			imageIn.pixptr = pixptr;
		}
	}
	if (isLazy) {
		// in the case of a lazy read, we must check in the FITSImage clone operator
		// if the extension has been read and proactively read it in...
		int status = 0;
		setHasBeenRead(imageIn);
		setHasBeenWritten(imageIn);
		if (!extensionHasBeenRead[current_extension]) {
			long fpixel[MAXFITSDIMENSIONS] = {1,1,1};
			int filedatatype;
			if (fits_get_img_equivtype(fptr, &filedatatype, &status)) {
				throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			float mybzero = bzero;
			if (filedatatype == USHORT_IMG) {
				mybzero = 32768.0;
			}
			if (fits_read_pix(fptr, todatatype(ebitpix(filedatatype), mybzero, bscale), fpixel, npixels, NULL, pixptr, NULL, &status)) {
				//throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
				status = 0;
			}
			if (todatatype(ebitpix(filedatatype), mybzero, bscale) != datatype) {
				operaFITSImageConvertImageInPlace(todatatype(ebitpix(filedatatype), bzero, bscale), datatype);
			}
			setHasBeenRead(current_extension);
		}
	}
	if (AddHeader) {
		int status = 0;
		int hdutype = 0;
        // Move to the current extension HDU (before calling this constructor, the current_extension must
        // have been set. A specific extension for MEF Images, or the primary header (extension 1) for a FITS Cube
        if ( fits_movabs_hdu(imageIn.fptr, cfitsioextension(0), &hdutype, &status) ) {
            throw operaException("operaFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if ( fits_movrel_hdu(imageIn.fptr, imageIn.current_extension, &hdutype, &status) ) {
            throw operaException("operaFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
        }
		hdu = fits_get_hdu_num(imageIn.fptr, &status);
	} else {
        hdu = 0; // no header to be copied. This will give a default header in save or saveAs
    }
}

/* 
 * \class operaFITSImage()
 * \brief Destructor releases the image array memory and closes the opened fitsfile.
 * \return void
 */
operaFITSImage::~operaFITSImage(){
    if (!super && !viewOnly && !isClone) {
		if (pixptr) free(pixptr);
		pixptr = NULL;
		if (varptr) free(varptr);
		varptr = NULL;
	}
}

/*
 * Utility routines
 */

/* 
 * operaFITSImageSave() 
 * \brief Saves the current image to disk.
 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
 * \throws operaException operaErrorCodeNoFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSImageSave() {
	int status = 0;
	long fnaxes[MAXFITSDIMENSIONS];			// FITS image dimension in file
	fitsfile *newfptr;		// FITS file pointer for the new image
	const long fpixel = 1;
	const int morekeys = 0;

	if (mode == READONLY)
		throw operaException("operaFITSImage:: "+filename+" ", operaErrorCodeChangeREADONLYError, __FILE__, __FUNCTION__, __LINE__);	
	
	if (filename.empty())
		throw operaException("operaFITSImage: ", operaErrorCodeNoFilename, __FILE__, __FUNCTION__, __LINE__);	
	
	// remove existing file - cfitsio returns an error if it exists...
	remove(filename.c_str());
	
	if (fits_create_file(&newfptr, filename.c_str(), &status)) {
		throw operaException("operaFITSImage: cfitsio create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fptr == NULL || hdu == 0) {
		if (fits_create_img(newfptr, bitpix, naxis, naxes, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	//
	// check for an image resize, if hdu exists
	//
	//cerr << "resizeing? " << naxes[0] << ' ' << naxes[1] << ' ' << naxis1 << ' ' << naxis2 << endl;
	if (hdu != 0 && fptr != 0) {
		//	get img size from FITS file
		if (fits_get_img_size(fptr, 2, fnaxes, &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}	
		//cerr << "resizeing,, " << fnaxes[0] << ' ' << fnaxes[1] << endl;
		if ((long)naxis1 != fnaxes[0] || (long)naxis2 != fnaxes[1]) {
			//cout << "resized " << naxes[0] << ' ' << naxes[1] << endl;
			if (fits_resize_img(newfptr, bitpix, 2, naxes, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
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
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (datatype == tfloat) {
				if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
					throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {
				if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
					throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	if (fits_update_key_lng(newfptr, "NAXIS1", naxis1, (char *)"Number of pixel columns", &status)) {
		throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key_lng(newfptr, "NAXIS2", naxis2, (char *)"Number of pixel rows", &status)) {
		throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_write_img(newfptr, datatype, fpixel, npixels, pixptr, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	fits_close_file(newfptr, &status);
	newfptr = NULL;
	// note that fptr is left open in case the user wants to save again.
	// It will be closed by a close call if made, or the destructor otherwise
}

/* 
 * operaFITSImageSaveAs(string newFilename) 
 * \brief Saves the current image to disk, with the given filename.
 * \param newFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSImageSaveAs(string newFilename) {
	int status = 0;
	long fnaxes[2];			// FITS image dimension in file
	fitsfile *newfptr;		// FITS file pointer for the new image
	long  fpixel = 1;
	
	// remove existing file - cfitsio returns an error if it exists...
	remove(newFilename.c_str());
	
    if (fits_create_file(&newfptr, newFilename.c_str(), &status)) {
		throw operaException("operaFITSImage: cfitsio error "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}		
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}	
	if (fptr == NULL || hdu == 0) {
		if (fits_create_img(newfptr, bitpix, naxis, naxes, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
		if (fits_copy_hdu(fptr, newfptr, 0, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	//
	// check for an image resize
	//
	if (hdu != 0) {
		//	get img size from FITS file if it exists in a header   
		if (fits_get_img_size(fptr, 2, fnaxes, &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}	
		if ((long)naxis1 != fnaxes[0] || (long)naxis2 != fnaxes[1]) {
			if (fits_resize_img(newfptr, bitpix, 2, naxes, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
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
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			if (datatype == tfloat) {
				if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
					throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {
				if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
					throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	if (fits_write_img(newfptr, datatype, fpixel, npixels, pixptr, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
	}
	fits_close_file(newfptr, &status);
	newfptr = NULL;
}

/* 
 * operaFITSImageCopyHeader(operaFITSImage *from) 
 * \brief Copies all of the header information from image.
 * \param from
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSImageCopyHeader(operaFITSImage *from) {
	int status = 0;
	int junk = 0;
	
	if (fits_copy_hdu(from->getfitsfileptr(), fptr, 0, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	//
	// update the current HDU
	//
	hdu = fits_get_hdu_num(fptr, &junk);
	//
	// now update the headers in case they differ
	//
	if (datatype == tfloat || datatype == tdouble) {
		float bzero = 0.0, bscale = 1.0;
		if (fits_resize_img(fptr, bitpix, naxis, naxes, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_set_bscale(fptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (datatype == tfloat) {
			if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		} else {
			if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
		if (fits_update_key_flt(fptr, "BZERO", bzero, -1, NULL, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_flt(fptr, "BSCALE", bscale, -1, NULL, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}


/*
 * operaFITSImageCopyHeader(operaFITSImage *from, unsigned inputhdu)
 * \brief Copies all of the header information from image.
 * \param from
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSImageCopyHeaderFromDifferentHDU(operaFITSImage *from, unsigned inputabshdu) {
	int status = 0;
	int junk = 0;
	int hdutype = 0;
    
    if (fits_movabs_hdu(from->getfitsfileptr(), inputabshdu, &hdutype, &status)) {
        throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
    }
	if (fits_copy_hdu(from->getfitsfileptr(), fptr, 0, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
	}
	//
	// update the current HDU
	//
	hdu = fits_get_hdu_num(fptr, &junk);
	//
	// now update the headers in case they differ
	//
	if (datatype == tfloat || datatype == tdouble) {
		float bzero = 0.0, bscale = 1.0;
		if (fits_resize_img(fptr, bitpix, naxis, naxes, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (fits_set_bscale(fptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (datatype == tfloat) {
			if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
		} else {
			if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
		}
		if (fits_update_key_flt(fptr, "BZERO", bzero, -1, NULL, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (fits_update_key_flt(fptr, "BSCALE", bscale, -1, NULL, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
	}
}



/*! 
 * operaFITSImageCopyNonStructuralHeader(operaFITSImage *from) 
 * \brief Copies all of the non-structural header information from image.
 * \param from
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSImageCopyNonStructuralHeader(operaFITSImage *from) {
	int status = 0;
	int nkeys = 0;
	char card[81];

	if (!fptr) {
		throw operaException("operaFITSImage: ", operaErrorCodeNULLfptr, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (!from->getfitsfileptr()) {
		throw operaException("operaFITSImage: ", operaErrorCodeNULLfptr, __FILE__, __FUNCTION__, __LINE__);	
	}
	/* copy all the user keywords (not the structural keywords) */
	if (fits_get_hdrspace(from->getfitsfileptr(), &nkeys, NULL, &status)) {
		throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	} 
	for (int ii = 1; ii <= nkeys; ii++) {
		if (fits_read_record(from->getfitsfileptr(), ii, card, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
			if (fits_write_record(fptr, card, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
}
/* 
 * string operaFITSGetFilename(void) 
 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
 * \param keyword
 * \return string keyword value or empty string
 */
string operaFITSImage::operaFITSGetFilename(void) {
	return filename;
}
/* 
 * void *row()
 * \brief get base of ith row.
 * \return void *
 */
float* operaFITSImage::row(unsigned i) {
	return (*this)[i];
}

/* 
 * string operaFITSFindComment(string contains) 
 * \brief returns the value of a FITS comment containing string.
 * \param keyword
 * \return string keyword value or empty string
 */
string operaFITSImage::operaFITSFindComment(string contains) {
	string result;
	int status = 0;
	int nkeys = 0;
	char card[81];
	
	if (fits_get_hdrspace(fptr, &nkeys, NULL, &status)) {
		throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	} 
	for (int ii = 1; ii <= nkeys; ii++) {
		if (fits_read_record(fptr, ii, card, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		string comment(card);
		if (comment.find(contains) != string::npos) {
			result = comment;
			break;
		}
	}
	return result;
}

/* 
 * string operaFITSGetHeaderValue(string keyword) 
 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
 * \param keyword
 * \return string keyword value or empty string
 */
string operaFITSImage::operaFITSGetHeaderValue(string keyword) {
	string raw = operaFITSGetRawHeaderValue(keyword);
	return trimFITSKeyword(raw.c_str());
}
/* 
 * float operaFITSGetFloatHeaderValue(string keyword) 
 * \brief returns the float value of a FITS keyword.
 * \param keyword
 * \return float keyword value 
 */
float operaFITSImage::operaFITSGetFloatHeaderValue(string keyword) {
	string raw = operaFITSGetRawHeaderValue(keyword);
	string clean = trimFITSKeyword(raw.c_str());
	return atof(clean.c_str());
}
/* 
 * int operaFITSGetIntHeaderValue(string keyword) 
 * \brief returns the int value of a FITS keyword.
 * \param keyword
 * \return int keyword value
 */
int operaFITSImage::operaFITSGetIntHeaderValue(string keyword) {
	string raw = operaFITSGetRawHeaderValue(keyword);
	string clean = trimFITSKeyword(raw.c_str());
	return atoi(clean.c_str());
}
/* 
 * string operaFITSGetRawHeaderValue(string keyword) 
 * \brief returns the value of a FITS keyword verbatim.
 * \param keyword
 * \throws operaException cfitsio error code
 * \return string keyword value or empty string
 */
string operaFITSImage::operaFITSGetRawHeaderValue(string keyword) {
	char zvalue[FLEN_VALUE],comment[FLEN_COMMENT];
	int status = 0;
	if (fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, comment, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" "+keyword+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	return string(zvalue);
}

/* 
 * void getFITSImageInformation(string Filename, &XDimension, &YDimension, &ZDimension, &Extensions, &Datatype, &npixels)
 * \brief get usefule information about a FITS file.
 * \param Filename
 * \param XDimension
 * \param YDimension
 * \param ZDimension
 * \param Extensions
 * \param Datatype
 * \param Npixels
 * \return void
 */
void getFITSImageInformation(string Filename, unsigned *XDimension, unsigned *YDimension, unsigned *ZDimension, unsigned *Extensions, edatatype *Datatype, long *Npixels){
	
	operaFITSImage info(Filename, tfloat, READONLY, cNone, true);
    
	*Datatype = info.getdatatype();
	*XDimension = info.getXDimension();
	*YDimension = info.getYDimension();
	*ZDimension = info.getZDimension();
	*Extensions = info.getNExtensions();
	*Npixels = (info.getNExtensions()>0?info.getnpixels()*info.getNExtensions():info.getnpixels());
	info.operaFITSImageClose();
}
/* 
 * operaFITSSetHeaderValue(string keyword, string value, string comment)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSSetHeaderValue(string keyword, string value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tstring, (char*)keyword.c_str(), (char*)value.c_str(), (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tstring, (char*)keyword.c_str(), (char*)value.c_str(), (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}

/* 
 * operaFITSSetHeaderValue(string keyword, unsigned short value, string comment)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \return void
 */
void operaFITSImage::operaFITSSetHeaderValue(string keyword, unsigned short value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tushort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tushort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}
/* 
 * operaFITSSetHeaderValue(string keyword, unsigned short value, string comment)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \return void
 */
void operaFITSImage::operaFITSSetHeaderValue(string keyword, short value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tshort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tshort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}
/* 
 * operaFITSSetHeaderValue(string keyword, float value, string comment)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSSetHeaderValue(string keyword, float value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tfloat, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tfloat, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}
/* 
 * operaFITSSetHeaderValue(string keyword, double value, string comment)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSSetHeaderValue(string keyword, double value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tdouble, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tdouble, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}

/* 
 * operaFITSDeleteHeaderKey(string keyword)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSDeleteHeaderKey(string keyword) {
	int status = 0;
	
	fits_delete_key(fptr, (char*)keyword.c_str(), &status);
}
/* 
 * operaFITSAddComment(string comment)
 * \brief adds a comment to header.
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSAddComment(string comment) {
	int status = 0;
	
	if (fits_write_comment(fptr, (char*)comment.c_str(), &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
}

/* 
 * unsigned short* operaFITSImage::operaFITSImageClonePixelsUSHORT()
 * \brief Clone pixel image data.
 * \return pixels*
 */
unsigned short* operaFITSImage::operaFITSImageClonePixelsUSHORT() {
	long npixels = naxis1 * naxis2;
	long size = sizeof(unsigned short)*npixels;
	unsigned short *p = (unsigned short *)malloc(size); 
	if (!p) {
		throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (pixptr)
		memcpy(p, pixptr, size);
	else
		memset(p, 0, size);
	return p;
}

/* 
 * unsigned short* operaFITSImage::operaFITSImageClonePixelsUSHORT(unsigned x, unsigned y, unisgned nx, unsigned ny)
 * \brief Clone a pixel image subwindow.
 * \param x
 * \param y
 * \param nx
 * \param ny
 * \return pixels*
 */
unsigned short* operaFITSImage::operaFITSImageClonePixelsUSHORT(unsigned x, unsigned y, unsigned nx, unsigned ny) {
	long npixels = nx * ny;
	long size = sizeof(float)*npixels;
	unsigned short *p = (unsigned short *)malloc(size);
	if (!p) {
		throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned short *pp = p;
	unsigned short *pptr = (unsigned short *)pixptr;
	for (unsigned j = y; j < ny; j++) {
		for (unsigned i = x; i < nx; i++) {
			if (pixptr)
				*p++ = (unsigned short)pptr[(naxis1*j)+i];
			else
				*p++ = 0;
		}
	}
	return pp;
}

/* 
 * float* operaFITSImage::operaFITSImageClonePixels()
 * \brief Clone float pixel data.
 * \return pixels*
 */
float* operaFITSImage::operaFITSImageClonePixels() {
	long npixels = naxis1 * naxis2;
	long size = sizeof(float)*npixels;
	float *p = (float *)malloc(size); 
	if (!p) {
		throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (pixptr)
		memcpy(p, pixptr, size);
	else
		memset(p, 0, size);
	return p;
}

/* 
 * float* operaFITSImage::operaFITSImageClonePixels(unsigned x, unsigned y, unsigned nx, unsigned ny)
 * \brief Clone a float pixel image subwindow.
 * \param x
 * \param y
 * \param nx
 * \param ny
 * \return pixels*
 */
float* operaFITSImage::operaFITSImageClonePixels(unsigned x, unsigned y, unsigned nx, unsigned ny) {
	long npixels = nx * ny;
	long size = sizeof(float)*npixels;
	float *p = (float *)malloc(size);
	if (!p) {
		throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *pp = p;
	switch (datatype) {
		case tshort: {
			short *pptr = (short *)pixptr;
			for (unsigned j = y; j < ny; j++) {
				for (unsigned i = x; i < nx; i++) {
					if (pixptr)
						*p++ = (float)pptr[(naxis1*j)+i];
					else
						*p++ = 0.0;
				}
			}
		}
			break;
		case tushort: {	
			unsigned short *pptr = (unsigned short *)pixptr;
			for (unsigned j = y; j < ny; j++) {
				for (unsigned i = x; i < nx; i++) {
					if (pixptr)
						*p++ = (float)pptr[(naxis1*j)+i];
					else
						*p++ = 0.0;
				}
			}
		}
			break;
		case tfloat: {
			float *pptr = (float *)pixptr;
			for (unsigned j = y; j < ny; j++) {
				for (unsigned i = x; i < nx; i++) {
					if (pixptr)
						*p++ = pptr[(naxis1*j)+i];
					else
						*p++ = 0.0;
				}
			}
		}
			break;
		case tdouble: {	
			double *pptr = (double *)pixptr;
			for (unsigned j = y; j < ny; j++) {
				for (unsigned i = x; i < nx; i++) {
					if (pixptr)
						*p++ = (float)pptr[(naxis1*j)+i];
					else
						*p++ = 0.0;
				}
			}
		}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;		
			break;
	}
	return pp;
}

/* 
 * void operaFITSImage::operaFITSImageSetData(operaFITSSubImage &subImage, unsigned X, unsigned Y)
 * \brief copy a subimage in to an operaFITSImage.
 * \param subImage the subImage
 * \param X - beginning x location in image
 * \param Y - beginning y location in image
 * \return void
 */
void operaFITSImage::operaFITSImageSetData(operaFITSSubImage &subImage, unsigned X, unsigned Y) {
	unsigned nx = subImage.nx;
	unsigned ny = subImage.ny;
	for (unsigned j=0; j<ny; j++) {
		for (unsigned i=0; i<nx; i++) {
			setpixel(subImage.getpixel(i, j), i+X, j+Y);
		}
	}
}

/* 
 * void operaFITSImage::operaFITSImageWriteVirtual(operaFITSSubImage &subImage, unsigned X, unsigned Y)
 * \brief Write a subimage on to avirtual operaFITSImage.
 * \brief This is used for deep stacking. 
 * \brief The virtual image is Lazy, and can be very large.
 * \brief The actual image exists only on disk.
 * \param subImage the subImage
 * \param X - beginning x location in image
 * \param Y - beginning y location in image
 * \return void
 */
void operaFITSImage::operaFITSImageWriteVirtual(operaFITSSubImage &subImage, unsigned long X, unsigned long Y) {
	unsigned nx = subImage.nx;
	unsigned ny = subImage.ny;
	int status = 0;
	long fpixel[MAXFITSDIMENSIONS], lpixel[MAXFITSDIMENSIONS];
	
	fpixel[0] = X; 
	fpixel[1] = Y;
	fpixel[2] = 1;
	lpixel[0] = nx; 
	lpixel[1] = ny;
	lpixel[2] = 1;
	
	if (isLazy) {
		if (fits_write_subset(fptr, datatype, fpixel, lpixel, subImage.getpixels(), &status))
			throw operaException("operaFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
	} else {
		throw operaException("operaFITSImage: "+filename+" ", operaErrorNotLazy, __FILE__, __FUNCTION__, __LINE__);	
	}
}

/* 
 * void operaFITSImage::operaFITSImageReadVirtual(operaFITSSubImage &subImage, unsigned X, unsigned Y)
 * \brief Read a subimage from a virtual operaFITSImage.
 * \brief This is used for deep stacking. 
 * \brief The virtual image is Lazy, and can be very large.
 * \brief The actual image exists only on disk.
 * \param subImage the subImage
 * \param X - beginning x location in image
 * \param Y - beginning y location in image
 * \return void
 */
void operaFITSImage::operaFITSImageReadVirtual(operaFITSSubImage &subImage, unsigned long X, unsigned long Y) {
	unsigned nx = subImage.nx;
	unsigned ny = subImage.ny;
	int status = 0;
	long fpixel[MAXFITSDIMENSIONS], lpixel[MAXFITSDIMENSIONS], inc[MAXFITSDIMENSIONS];
	
	inc[0] = 1; 
	inc[1] = 1;
	inc[2] = 1;
	fpixel[0] = X; 
	fpixel[1] = Y;
	fpixel[2] = 1;
	lpixel[0] = nx; 
	lpixel[1] = ny;
	lpixel[2] = 1;
	
	if (isLazy) {
		if (fits_read_subset(fptr, datatype, fpixel, lpixel, inc, NULL, subImage.getpixels(), NULL, &status))
			throw operaException("operaFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	} else {
		throw operaException("operaFITSImage: "+filename+" ", operaErrorNotLazy, __FILE__, __FUNCTION__, __LINE__);	
	}
}

/* 
 * void operaFITSImage::operaFITSImageSetData(unsigned short* data)
 * \brief set the image data pointer to a buffer of data.
 * \param data pointer to the data
 * \return void
 */
void operaFITSImage::operaFITSImageSetData(unsigned short* data) {
	pixptr = (void *)data;
}

/* 
 * void operaFITSImage::operaFITSImageSetData(float* data)
 * \brief set the image data pointer to a buffer of data.
 * \param data pointer to the data
 * \return void
 */
void operaFITSImage::operaFITSImageSetData(float* data) {
	pixptr = (void *)data;
}

/* 
 * void operaFITSImage::operaFITSImageConvertImage(edatatype todatattype)
 * \brief Convert image type and create new image values of that type.
 * \param todatatype edatatype to convert to
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSImageConvertImage(edatatype todatatype) {
	int status = 0;
	int i = npixels;
	edatatype olddatatype = datatype;
	datatype = todatatype;
	bitpix = tobitpix(todatatype); 
	size_t newsize = toSize(bitpix, npixels);
	void *newpixptr = malloc(newsize);
	if (!newpixptr) {
		throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	switch (todatatype) {
		case tshort:
			switch (olddatatype) {
				case tshort: {
					short *topp = (short *)newpixptr;
					short *from = (short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tushort: {
					short *topp = (short *)newpixptr;
					unsigned short *from = (unsigned short *)pixptr;
					while (i--) {
						*topp++ = (unsigned short)(*from++);
					}
				}
					break;
				case tfloat: {
					short *topp = (short *)newpixptr;
					float *from = (float *)pixptr;
					while (i--) {
						*topp++ = (short)(*from++);
					}
				}
					break;
				case tdouble: {
					short *topp = (short *)newpixptr;
					double *from = (double *)pixptr;
					while (i--) {
						*topp++ = (short)(*from++);
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			break;
		case tushort:
			switch (olddatatype) {
				case tshort: {
					unsigned short *topp = (unsigned short *)newpixptr;
					short *from = (short *)pixptr;
					while (i--) {
						*topp++ = (unsigned short)(*from++);
					}
				}
					break;
				case tushort: {
					unsigned short *topp = (unsigned short *)newpixptr;
					unsigned short *from = (unsigned short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tfloat: {
					unsigned short *topp = (unsigned short *)newpixptr;
					float *from = (float *)pixptr;
					while (i--) {
						*topp++ = (unsigned short)(*from++);
					}
				}
					break;
				case tdouble: {
					unsigned short *topp = (unsigned short *)newpixptr;
					double *from = (double *)pixptr;
					while (i--) {
						*topp++ = (unsigned short)(*from++);
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			break;
		case tfloat:
			switch (olddatatype) {
				case tshort: {
					float *topp = (float *)newpixptr;
					short *from = (short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tushort: {
					float *topp = (float *)newpixptr;
					unsigned short *from = (unsigned short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tfloat: {
					float *topp = (float *)newpixptr;
					float *from = (float *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tdouble: {
					float *topp = (float *)newpixptr;
					double *from = (double *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}
			if (fits_set_bscale(fptr, 1.0, 0.0, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			break;
		case tdouble:
			switch (olddatatype) {
				case tshort: {
					double *topp = (double *)newpixptr;
					short *from = (short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tushort: {
					double *topp = (double *)newpixptr;
					unsigned short *from = (unsigned short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tfloat: {
					double *topp = (double *)newpixptr;
					float *from = (float *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tdouble: {
					double *topp = (double *)newpixptr;
					double *from = (double *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			if (fits_set_bscale(fptr, 1.0, 0.0, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;
			break;
	}
	if (pixptr) free(pixptr);
	pixptr = newpixptr;
}
/* 
 * void operaFITSImage::operaFITSImageConvertImage(edatatype fromdatatype, edatatype todatatype)
 * \brief Convert image type and create new image values of that type.
 * \param todatatype edatatype to convert to
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSImageConvertImage(edatatype fromdatatype, edatatype todatatype) {
	int status = 0;
	int i = npixels;
	edatatype olddatatype = fromdatatype;
	datatype = todatatype;
	bitpix = tobitpix(todatatype); 
	size_t newsize = toSize(bitpix, npixels);
	void *newpixptr = malloc(newsize);
	if (!newpixptr) {
		throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	switch (todatatype) {
		case tshort:
			switch (olddatatype) {
				case tshort: {
					short *topp = (short *)newpixptr;
					short *from = (short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tushort: {
					short *topp = (short *)newpixptr;
					unsigned short *from = (unsigned short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tfloat: {
					short *topp = (short *)newpixptr;
					float *from = (float *)pixptr;
					while (i--) {
						*topp++ = (short)(*from++);
					}
				}
					break;
				case tdouble: {
					short *topp = (short *)newpixptr;
					double *from = (double *)pixptr;
					while (i--) {
						*topp++ = (short)(*from++);
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			break;
		case tushort:
			switch (olddatatype) {
				case tshort: {
					unsigned short *topp = (unsigned short *)newpixptr;
					short *from = (short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tushort: {
					unsigned short *topp = (unsigned short *)newpixptr;
					unsigned short *from = (unsigned short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tfloat: {
					unsigned short *topp = (unsigned short *)newpixptr;
					float *from = (float *)pixptr;
					while (i--) {
						*topp++ = (unsigned short)(*from++);
					}
				}
					break;
				case tdouble: {
					unsigned short *topp = (unsigned short *)newpixptr;
					double *from = (double *)pixptr;
					while (i--) {
						*topp++ = (unsigned short)(*from++);
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			break;
		case tfloat:
			switch (olddatatype) {
				case tshort: {
					float *topp = (float *)newpixptr;
					short *from = (short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tushort: {
					float *topp = (float *)newpixptr;
					unsigned short *from = (unsigned short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tfloat: {
					float *topp = (float *)newpixptr;
					float *from = (float *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tdouble: {
					float *topp = (float *)newpixptr;
					double *from = (double *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}
			if (fits_set_bscale(fptr, 1.0, 0.0, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			break;
		case tdouble:
			switch (olddatatype) {
				case tshort: {
					double *topp = (double *)newpixptr;
					short *from = (short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tushort: {
					double *topp = (double *)newpixptr;
					unsigned short *from = (unsigned short *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tfloat: {
					double *topp = (double *)newpixptr;
					float *from = (float *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				case tdouble: {
					double *topp = (double *)newpixptr;
					double *from = (double *)pixptr;
					while (i--) {
						*topp++ = *from++;
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			if (fits_set_bscale(fptr, 1.0, 0.0, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;
			break;
	}
	if (pixptr) free(pixptr);
	pixptr = newpixptr;
}
/* 
 * void operaFITSImage::operaFITSImageConvertImageInPlace(edatatype fromdatatype, edatatype todatatype)
 * \brief Convert image type and create new image values of that type.
 * \note - this assumes enough storage is available for the target type.
 * \param todatatype edatatype to convert to
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::operaFITSImageConvertImageInPlace(edatatype fromdatatype, edatatype todatatype) {
	int status = 0;
	int i = npixels;
	edatatype olddatatype = fromdatatype;
	datatype = todatatype;
	bitpix = tobitpix(todatatype); 
	
	switch (todatatype) {
		case tshort:
			switch (olddatatype) {
				case tshort: {
					short *topp = (short *)pixptr+npixels-1;
					short *from = (short *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tushort: {
					short *topp = (short *)pixptr+npixels-1;
					unsigned short *from = (unsigned short *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tfloat: {
					short *topp = (short *)pixptr+npixels-1;
					float *from = (float *)pixptr+npixels-1;
					while (i--) {
						*topp-- = (short)(*from--);
					}
				}
					break;
				case tdouble: {
					short *topp = (short *)pixptr+npixels-1;
					double *from = (double *)pixptr+npixels-1;
					while (i--) {
						*topp-- = (short)(*from--);
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			break;
		case tushort:
			switch (olddatatype) {
				case tshort: {
					unsigned short *topp = (unsigned short *)pixptr+npixels-1;
					short *from = (short *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tushort: {
					unsigned short *topp = (unsigned short *)pixptr+npixels-1;
					unsigned short *from = (unsigned short *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tfloat: {
					unsigned short *topp = (unsigned short *)pixptr+npixels-1;
					float *from = (float *)pixptr+npixels-1;
					while (i--) {
						*topp-- = (unsigned short)(*from--);
					}
				}
					break;
				case tdouble: {
					unsigned short *topp = (unsigned short *)pixptr+npixels-1;
					double *from = (double *)pixptr+npixels-1;
					while (i--) {
						*topp-- = (unsigned short)(*from--);
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			break;
		case tfloat:
			switch (olddatatype) {
				case tshort: {
					float *topp = (float *)pixptr+npixels-1;
					short *from = (short *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tushort: {
					float *topp = (float *)pixptr+npixels-1;
					unsigned short *from = (unsigned short *)pixptr+npixels-1;;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tfloat: {
					float *topp = (float *)pixptr+npixels-1;
					float *from = (float *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tdouble: {
					float *topp = (float *)pixptr+npixels-1;
					double *from = (double *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}
			if (fits_set_bscale(fptr, 1.0, 0.0, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			break;
		case tdouble:
			switch (olddatatype) {
				case tshort: {
					double *topp = (double *)pixptr+npixels-1;
					short *from = (short *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tushort: {
					double *topp = (double *)pixptr+npixels-1;
					unsigned short *from = (unsigned short *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tfloat: {
					double *topp = (double *)pixptr+npixels-1;
					float *from = (float *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				case tdouble: {
					double *topp = (double *)pixptr+npixels-1;
					double *from = (double *)pixptr+npixels-1;
					while (i--) {
						*topp-- = *from--;
					}
				}
					break;
				default:
					throw operaErrorCodeDatatypeNotSupported;
					break;
			}	    
			if (fits_set_bscale(fptr, 1.0, 0.0, &status)) {
				throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;
			break;
	}
}
/* 
 * operaFITSImage::rotate90()
 * \brief rotate 90 degrees.
 * \return void
 */
void operaFITSImage::rotate90() {
	int status = 0;
	switch (datatype) {
		case tushort: {
			unsigned short *clone = operaFITSImageClonePixelsUSHORT();
			for (unsigned j = 0; j < naxis2; j++) {
				for (unsigned i = 0; i < naxis1; i++) {
					setpixel(getpixelUSHORT(clone, i, j, naxis1), j, naxis1-i-1, naxis2); 
				}
			}
			free(clone);
		}
			break;
		case tfloat: {
			float *clone = operaFITSImageClonePixels();
			for (unsigned j = 0; j < naxis2; j++) {
				for (unsigned i = 0; i < naxis1; i++) {
					setpixel(getpixel(clone, i, j, naxis1), j, naxis1-i-1, naxis2); 
				}
			}
			unsigned tmp = naxis2;
			naxis2 = naxis1;
			naxis1 = tmp;
			free(clone);
		}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;
			break;
	}
	//
	// be sure to update the header keyword if there is an fptr
	unsigned tmp = naxis2;
	naxis2 = naxis1;
	naxis1 = tmp;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	//cerr << "rotate90 " << naxes[0] << ' ' << naxes[1] << ' ' << naxis1 << ' ' << naxis2 << endl;
	if (fptr) {
		if (fits_update_key_lng(fptr, "NAXIS1", naxis1, (char *)"Number of pixel columns", &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(fptr, "NAXIS2", naxis2, (char *)"Number of pixel rows", &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}

/*
 * operaFITSImage::mirrorColumns()
 * \brief mirror x.
 * \return void
 */
void operaFITSImage::mirrorColumns(void) {
	switch (datatype) {
		case tushort: {
			unsigned short *clone = operaFITSImageClonePixelsUSHORT();
			for (unsigned j = 0; j < naxis2; j++) {
				for (unsigned i = 0; i < naxis1; i++) {
					setpixel(getpixelUSHORT(clone, i, j, naxis1),naxis1-1-i,j);
				}
			}
			free(clone);
		}
			break;
		case tfloat: {
			float *clone = operaFITSImageClonePixels();
			for (unsigned j = 0; j < naxis2; j++) {
				for (unsigned i = 0; i < naxis1; i++) {
					setpixel(getpixel(clone, i, j, naxis1),naxis1-1-i,j);
				}
			}
			free(clone);
		}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;
			break;
	}
}

/*
 * operaFITSImage::mirrorRows()
 * \brief mirror y.
 * \return void
 */
void operaFITSImage::mirrorRows(void) {
	switch (datatype) {
		case tushort: {
			unsigned short *clone = operaFITSImageClonePixelsUSHORT();
			for (unsigned j = 0; j < naxis2; j++) {
				for (unsigned i = 0; i < naxis1; i++) {
					setpixel(getpixelUSHORT(clone, i, j, naxis1),i,naxis2-1-j);
				}
			}
			free(clone);
		}
			break;
		case tfloat: {
			float *clone = operaFITSImageClonePixels();
			for (unsigned j = 0; j < naxis2; j++) {
				for (unsigned i = 0; i < naxis1; i++) {
					setpixel(getpixel(clone, i, j, naxis1),i,naxis2-1-j);
				}
			}
			free(clone);
		}
			break;
		default:
			throw operaErrorCodeDatatypeNotSupported;
			break;
	}
}

/* 
 * operaFITSImage::assignVariances(float gain)
 * \brief Assign variances to each pixel in an image in units of ADU.
 */
void operaFITSImage::assignVariances(float gain) {
	if (varptr == NULL) {
		varptr = malloc(sizeof(float)*npixels);
	}
	float *pixels = (float *)pixptr;
	float *vars = (float *)varptr;
	unsigned long n = npixels;
	while (n--) {
		float noise = sqrt(*pixels++ * gain) / gain;
		*vars++ = pow(noise,2);
	}
}

/* 
 * operaFITSImage::transpose()
 * \brief transpose y for x.
 * \return void
 */
void operaFITSImage::transpose(operaFITSImage &inImage) {
	
	if (inImage.naxis1 != naxis2 || inImage.naxis2 != naxis1) {
		throw operaException("operaFITSImage: ",(operaErrorCode)operaErrorInputValuesDisagree, __FILE__, __FUNCTION__, __LINE__);		
	}
	unsigned ii = inImage.naxis1-1;
	for (unsigned j = 0; j < naxis2; j++) {
		unsigned jj = inImage.naxis2-1;
		for (unsigned i = 0; i < naxis1; i++) {
			(*this)[j][i] = inImage[jj][ii];
			jj--;
		}
		ii--;
	}
}

/* 
 * operaFITSImage::operaFITSImageClose()
 * \brief Close a FITSImage object file.
 * \brief Be careful if this is a view or a clone to not close the fptr.
 * \return void
 */
void operaFITSImage::operaFITSImageClose() {
    if (fptr && !viewOnly && !isClone) {
        int status = 0;
        fits_close_file(fptr, &status) ;
        fptr = NULL;
    }
}

/*
 * getters and setters
 */

/* 
 * setfilename(string Filename)
 * \brief set the image filename.
 * \param Filename
 * \return void
 */
void operaFITSImage::setfilename(string Filename) {
	filename = Filename;
}
/* 
 * string getfilename()
 * \brief set the image filename.
 * \return string filename
 */
string operaFITSImage::getfilename() {
	return filename;
}
/* 
 * fitsfile *getfitsfileptr()
 * \brief get the image FITS file point.
 * \return fitsfile *fileptr
 */
fitsfile *operaFITSImage::getfitsfileptr() {
	return fptr;
}
/* 
 * void setfitsfileptr(fitsfile *newfptr)
 * \brief set the image FITS file pointer.
 * \return void
 */
void operaFITSImage::setfitsfileptr(fitsfile *newfptr) {
	fptr = newfptr;
}
/* 
 * void *getpixels()
 * \brief get the image array base.
 * \return void *pixels pointer
 */
void *operaFITSImage::getpixels() {
	return pixptr;
}
/* 
 * void *getvars()
 * \brief get the variance array base.
 */
void *operaFITSImage::getvars() {
	return varptr;
}
/* 
 * unsigned getnaxis()
 * \brief get the image array number of axes (2 in usually).
 * \return unsigned number of axes
 */
unsigned operaFITSImage::getnaxis() {
	return naxis;
}
/* 
 * unsigned getnaxis1()
 * \brief get the image array length of x dimension.
 * \return unsigned length of axis 1
 */
unsigned operaFITSImage::getnaxis1() {
	return naxis1;
}
/* 
 * unsigned getnaxis2()
 * \brief get the image array length of y dimension.
 * \return unsigned length of axis 2
 */
unsigned operaFITSImage::getnaxis2() {
	return naxis2;
}
/* 
 * unsigned getnaxis3()
 * \brief get the image array length of y dimension.
 * \return unsigned length of axis 3
 */
unsigned operaFITSImage::getnaxis3() {
	return naxis3;
}
/* 
 * unsigned getXDimension()
 * \brief get the image array length of x dimension.
 * \return unsigned length of axis 1
 */
unsigned operaFITSImage::getXDimension() {
	return naxis1;
}
/* 
 * unsigned getYDimension()
 * \brief get the image array length of y dimension.
 * \return unsigned length of axis 2
 */
unsigned operaFITSImage::getYDimension() {
	return naxis2;
}
/* 
 * unsigned getZDimension()
 * \brief get the image array length of z dimension.
 * \return unsigned length of axis 3
 */
unsigned operaFITSImage::getZDimension() {
	return naxis3;
}
/* 
 * unsigned getNExtensions()
 * \brief get the number of extensions in the file.
 * \return unsigned number of extensions
 */
unsigned operaFITSImage::getNExtensions() {
	return extensions;
}
/* 
 * unsigned getExtension()
 * \brief get the current extension.
 * \return unsigned extension
 */
unsigned operaFITSImage::getExtension() {
	return current_extension;
}
/*
 * bool memoryAvailable()
 * \brief Is there enough memory available to read the file in?
 * \return bool
 */
bool operaFITSImage::memoryAvailable(long bytes){
	return true;
}
/* 
 * unsigned getnpixels()
 * \brief get the image array number of pixels.
 * \return unsigned number of pixels
 */
unsigned operaFITSImage::getnpixels() {
	return npixels;
}
/* 
 * unsigned getNPixelsPerExtension()
 * \brief get the image array number of pixels in an extension.
 * \return unsigned number of pixels
 */
unsigned operaFITSImage::getNPixelsPerExtension() {
	return npixels_per_extension;
}
/* 
 * unsigned long getsize()
 * \brief get the size of an image.
 * \return unsigned long sie of the image
 */
unsigned long operaFITSImage::getsize() {
	switch (datatype) {
		case tshort: {
			return sizeof(short)*npixels;				
		}
			break;
		case tushort: {	
			return sizeof(unsigned short)*npixels;				
		}
			break;
		case tfloat: {
			return sizeof(float)*npixels;				
		}
			break;
		case tdouble: {	
			return sizeof(double)*npixels;				
		}
			break;
		default:
            throw operaException("operaFITSImage: "+filename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);		
			break;
	}
}
/* 
 * unsigned getelementsize()
 * \brief get the size of a pixel.
 * \return unsigned long sie of the image
 */
unsigned operaFITSImage::getelementsize() {
	switch (datatype) {
		case tshort: {
			return sizeof(short);				
		}
			break;
		case tushort: {	
			return sizeof(unsigned short);				
		}
			break;
		case tfloat: {
			return sizeof(float);				
		}
			break;
		case tdouble: {	
			return sizeof(double);				
		}
			break;
		default:
            throw operaException("operaFITSImage: "+filename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);		
			break;
	}
}
/* 
 * unsigned getCompression()
 * \brief get the current compression method.
 * \return 0 = none or the FITS compression constant (usually RICE_1)
 */
unsigned operaFITSImage::getCompression() {
	return compression;
}
/* 
 * void setCompression(nsigned Compression)
 * \brief set the current compression method.
 * \param Compression 0 = none or the FITS compression constant (usually RICE_1)
 * \return void
 */
void operaFITSImage::setCompression(unsigned Compression) {
	compression = Compression;
}
/* 
 * edatatype getdatatype()
 * \brief get the current datatype
 * \return datatype
 */
edatatype operaFITSImage::getdatatype() {
	return datatype;
}
/* 
 * ebitpix getbitpix()
 * \brief get the current bitpix
 * \return bitpix
 */
ebitpix operaFITSImage::getbitpix() {
	return bitpix;
}

/*! 
 * setHasBeenRead(operaFITSImage &image, bool Read)
 * \brief Copy hasBeenRead property.
 */
void operaFITSImage::setHasBeenRead(operaFITSImage &image) {
	for (unsigned ext=0; ext<=MAXFITSEXTENSIONS; ext++) {
		extensionHasBeenRead[ext] = image.extensionHasBeenRead[ext];
	}		
}

/*! 
 * setHasBeenRead(unsigned index, bool Read)
 * \brief Set hasBeenRead property index = 0 means all.
 */
void operaFITSImage::setHasBeenRead(unsigned index, bool Read) {
	if (index == 0) {
		for (unsigned ext=0; ext<=MAXFITSEXTENSIONS; ext++) {
			extensionHasBeenRead[ext] = Read;
		}		
	} else {
		extensionHasBeenRead[index] = Read;
	}
}

/*! 
 * setHasBeenRead(unsigned index, bool Read)
 * \brief Set hasBeenRead property for this extension, resetting all others.
 */
void operaFITSImage::setHasBeenRead(unsigned index) {
	setHasBeenRead(0, false);
	extensionHasBeenRead[index] = true;
}

/*! 
 * setHasBeenWritten(operaFITSImage &image, bool Read)
 * \brief Copy hasBeenWritten property.
 */
void operaFITSImage::setHasBeenWritten(operaFITSImage &image) {
	for (unsigned ext=0; ext<=MAXFITSEXTENSIONS; ext++) {
		extensionHasBeenWritten[ext] = image.extensionHasBeenWritten[ext];
	}		
}

/*! 
 * setHasBeenWritten(unsigned index, bool Read)
 * \brief Set hasBeenWritten property index = 0 means all.
 */
void operaFITSImage::setHasBeenWritten(unsigned index, bool Read) {
	if (index == 0) {
		for (unsigned ext=0; ext<=MAXFITSEXTENSIONS; ext++) {
			extensionHasBeenWritten[ext] = Read;
		}		
	} else {
		extensionHasBeenWritten[index] = Read;
	}
}

/*! 
 * setHasBeenWritten(unsigned index, bool Read)
 * \brief Set hasBeenWritten property for this extension, resetting all others.
 */
void operaFITSImage::setHasBeenWritten(unsigned index) {
	setHasBeenWritten(0, false);
	extensionHasBeenWritten[index] = true;
}


operaImageVector *operaFITSImage::where(Box &box) {
	unsigned npixels = (box.getX2() - box.getX1()) * (box.getY2() - box.getY1());
	operaImageVector *ivector = new operaImageVector(npixels);
	unsigned d = 0;
	for (unsigned j=box.getY1(); j<=box.getY2(); j++) {
		for (unsigned i=box.getX1(); i<=box.getX2(); i++) {
			ivector->setIndex(d++, i+naxis1*j);
		}
	}
	return ivector;
}

operaImageVector *operaFITSImage::where(unsigned x, unsigned y, unsigned dx, unsigned dy) {
	unsigned npixels = dx * dy;
	operaImageVector *ivector = new operaImageVector(npixels);
	unsigned d = 0;
	for (unsigned j=y; j<=y+dy; j++) {
		for (unsigned i=x; i<=x+dx; i++) {
			ivector->setIndex(d++, i+naxis1*j);
		}
	}
	return ivector;
}

/*
 * \method resize(unsigned x0, unsigned xf, unsigned y0, unsigned yf)
 * \brief resize a product image to the given rows, cols, retaining any existing data
 * \note Resizing rows does not require copying data, where resizing columns does
 * \return void
 */
void operaFITSImage::resize(unsigned x0, unsigned xf, unsigned y0, unsigned yf) {
    
    if((x0 == 0 && xf == 0 && y0 == 0 && yf == 0) ||
       x0 >= xf  || y0 >= yf || xf > getnaxis1() || yf > getnaxis2()) {
        throw operaException("operaFITSImage: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    unsigned cols = xf-x0;
    unsigned rows = yf-y0;
    unsigned oldcols = naxis1;
    unsigned oldrows = naxis2;
    
	int status = 0;
	
	if (cols && rows && (oldrows != rows || oldcols != cols)) {
        switch (datatype) {
            case tushort: {
                unsigned short *clone = operaFITSImageClonePixelsUSHORT();
                for (unsigned j = 0; j < rows; j++) {
                    unsigned y = y0 + j;
                    for (unsigned i = 0; i < cols; i++) {
                        unsigned x = x0 + i;
                        setpixel(getpixelUSHORT(clone, x, y, oldcols),i,j,cols);
                    }
                }
                free(clone);
            }
                break;
            case tfloat: {
                float *clone = operaFITSImageClonePixels();
                for (unsigned j = 0; j < rows; j++) {
                    unsigned y = y0 + j;
                    for (unsigned i = 0; i < cols; i++) {
                        unsigned x = x0 + i;
                        setpixel(getpixel(clone, x, y, oldcols),i,j,cols);
                    }
                }
                free(clone);
            }
                break;
            default:
                throw operaErrorCodeDatatypeNotSupported;
                break;
        }

        naxis1 = cols;
        naxis2 = rows;
        npixels = naxis1*naxis2;
        naxes[0] = naxis1;
        naxes[1] = naxis2;
        
        if (fptr) {
            if (fits_update_key_lng(fptr, "NAXIS1", naxis1, (char *)"Number of pixel columns", &status)) {
                throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
            }
            if (fits_update_key_lng(fptr, "NAXIS2", naxis2, (char *)"Number of pixel rows", &status)) {
                throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
            }
        }
    }
}



/* 
 * \method resize(unsigned cols, unsigned rows)
 * \brief resize a product image to the given rows, cols, retaining any existing data
 * \note Resizing rows does not require copying data, where resizing columns does
 * \return void
 */
void operaFITSImage::resize(unsigned cols, unsigned rows) {
	int status = 0;
	unsigned oldcols = naxis1;
	unsigned oldrows = naxis2;
	naxis1 = cols;
	naxis2 = rows;
	npixels = naxis1*naxis2;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	
	if (npixels && (oldrows != rows || oldcols != cols)) {
		bitpix = tobitpix(datatype);
		size_t size = toSize(bitpix, npixels);	
		if (oldcols == cols) {	// let realloc do the work
			pixptr = realloc(pixptr, size);
		} else {
			void *newpixptr = malloc(size); 
			if (!newpixptr) {
				throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
			}
			memset(newpixptr, 0, size);
			if (pixptr)
				free(pixptr);
			pixptr = newpixptr;
		}
		// update the headers
		if (fits_resize_img(fptr, bitpix, 2, naxes, &status)) {
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}
/*
 * PRIVATE
 */
/* 
 * ebitpix tobitpix(edatatype Datatype)
 * \brief returns the bitpix value from given datatype.
 * \note PRIVATE
 * \param Datatype
 * \return ebitpix
 */
ebitpix operaFITSImage::tobitpix(edatatype Datatype) {
	
	switch (Datatype) {
		case tshort:
			return short_img;
			break;
		case tushort:
			return ushort_img;
			break;
		case tfloat:
			return float_img;
			break;
		case tdouble:
			return double_img;
			break;
		case tbyte:
			return byte_img;
			break;
		default:
			throw operaException("operaFITSImage: cfitsio error ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
			break;
	}	    
}	

/* 
 * edatatype todatatype(ebitpix Bitpix, float bzero, float bscale)
 * \brief returns the datatype value from given bitpix, bzero and bscale.
 * \note if bzero == 32768. and bscale == 1. then datatype = tushort else tshort
 * \note PRIVATE
 * \param Bitpix
 * \param bzero
 * \param bscale
 * \return edatatype
 */
//	datatype = todatatype(bitpix, bzero, bscale); <- call from readFITSHeaderInfo
edatatype operaFITSImage::todatatype(ebitpix Bitpix, float mybzero, float mybscale) {

	switch (Bitpix) {
		case short_img:
		case ushort_img:
			if ((mybzero == 32768. || mybzero == 65535.) && (mybscale == 1.))
				return tushort;
			else
				return tshort;
			break;
		case float_img:
			return tfloat;
			break;
		case double_img:
			return tdouble;
			break;
		case byte_img:
			return tbyte;
			break;
		default:
			throw operaException("operaFITSImage: datatype conversion error ("+itos(Bitpix)+") ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
			break;
	}
}

/* 
 * long toSize(ebitpix bitpix, long npixels)
 * \brief returns the size of an image array given the bitpix and npixels.
 * \note PRIVATE
 * \param bitpix
 * \param npixels
 * \return size
 */
size_t operaFITSImage::toSize(ebitpix bitpix, long npixels) {
	switch (bitpix) {
		case short_img:
			return npixels * sizeof(short);
			break;
		case ushort_img:
			return npixels * sizeof(unsigned short);
			break;
		case float_img:
			return npixels * sizeof(float);
			break;
		case double_img:
			return npixels * sizeof(double);
			break;
		case byte_img:
			return npixels * sizeof(char);
			break;
		default:
			throw operaException("operaFITSImage: datatype conversion error ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
			break;
	}
}

/* 
 * void openFITSfile(string Filename, int Mode=READWRITE)
 * \brief Opens a fitsfile of a given name and sets the filename into the object.
 * \note PRIVATE
 * \param Filename
 * \param mode (READONLY or READWRITE, defaults to latter
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::openFITSfile(string Filename, int Mode) {
	int status = 0;
	mode = Mode;
	
	filename = Filename;
	//	open FITS file	
	if (fits_open_diskfile(&fptr, filename.c_str(), mode, &status)) {
		throw operaException("operaFITSImage: cfitsio error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	status = 0;
	if (fits_is_compressed_image(fptr, &status)) {
		if (fits_get_compression_type(fptr, (int*)&compression, &status))
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
}

/* 
 * void readFITSHeaderInfo()
 * \brief reads the FITS header into the class object.
 * \note PRIVATE
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::readFITSHeaderInfo() {
	int status = 0;
	int hdutype = 0;
	int filedatatype;
	char zbitpix[FLEN_VALUE],zbzero[FLEN_VALUE],zbscale[FLEN_VALUE],zextend[FLEN_VALUE],comment[FLEN_COMMENT];
	
	naxes[0] = 1;
	naxes[1] = 1;
	naxes[2] = 1;
	//	get current hdunum
	fits_get_hdu_num(fptr, (int *)&hdu);
	
	//	get img dimension from FITS file 
	if (fits_get_img_dim(fptr, (int *)&naxis, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	//	get bitpix (datatype) from FITS file    
	if (fits_read_keyword(fptr, "BITPIX", zbitpix, comment, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	bitpix = (ebitpix)atoi(zbitpix); 
	extensions = 0;
	nhdus = 1;
	// if naxis == 0 then this is a PHU
	if (naxis == 0) {	// wircam
						//	get number of extensions    
		if (fits_read_keyword(fptr, "NEXTEND", zextend, comment, &status)) {
			status = 0;
			int hdus = 0;
			if (fits_get_num_hdus(fptr, &hdus, &status)) {
				throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			extensions = nhdus = (unsigned)hdus - 1;	// skip the primary
		} else {
			nhdus = atoi(zextend);
			extensions = nhdus;
		}
		hdu = 1;
		if (fits_movrel_hdu(fptr, hdu, &hdutype, &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		//	get img dimension from FITS file 
		if (fits_get_img_dim(fptr, (int *)&naxis, &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		//	get bitpix (datatype) from FITS file  -- Note that the extension bitpix may differ from PHU  
		if (fits_read_keyword(fptr, "BITPIX", zbitpix, comment, &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		bitpix = (ebitpix)atoi(zbitpix); 

		
	} 
	//	get img size from FITS file    
	if (fits_get_img_size(fptr, MAXFITSDIMENSIONS, naxes, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}	
	
	//	get img type from FITS file   - note you can't just use bitpix, because the
	// distinction between tushort and tshort also depends on bscale and bzero
	if (fits_get_img_equivtype(fptr, &filedatatype, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	naxis1 = naxes[0];
	naxis2 = naxes[1];
	naxis3 = ((naxis==2)?1:naxes[2]);
	current_extension = 1;
	current_slice = 1;
	npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
	npixels_per_slice = naxis1*naxis2;
	npixels_per_extension = naxis1*naxis2*naxis3;
	
	bzero = 0.0;
	//	get bitpix (datatype) from FITS file    
	fits_read_keyword(fptr, "BZERO", zbzero, comment, &status);
	if (status == 0) {
		bzero = atof(zbzero);    
	}
	status = 0;
	bscale = 1.0;
	//	get bitpix (datatype) from FITS file    
	fits_read_keyword(fptr, "BSCALE", zbscale, comment, &status);
	if (status == 0) {
		bzero = atof(zbzero);    
	}
	datatype = todatatype(bitpix, bzero, bscale);		
}

/* 
 * void readFITSArray()
 * \brief reads the image array into the class object (allocates memory).
 * \note PRIVATE
 * \throws operaException cfitsio error code
 * \return void
 */
void operaFITSImage::readFITSArray() {
	int status = 0;
	int hdutype = ANY_HDU;
	int filedatatype;
	long fpixel[MAXFITSDIMENSIONS] = {1,1,1};
	size_t size = toSize(bitpix, npixels);
	
	pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
	/*
	 * If we are reading all extensions one adjacent to the other then loop through each
	 * extension, incrementing the base pointer.
	 */
	if (AllExtensions && extensions > 0) {
		void *readpointer = malloc(npixels_per_extension*sizeof(double));	// biggest size
		size_t extensionsize = toSize(bitpix, npixels_per_extension);
		void *basepointer = pixptr;
		for (unsigned extension=first_extension; extension<=extensions; extension++) {
			/* select extension */
			if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
				throw operaException("operaFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if ( fits_movrel_hdu(fptr, extension, &hdutype, &status) ) {
				throw operaException("operaFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_get_img_equivtype(fptr, &filedatatype, &status)) {
				throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_read_pix(fptr,todatatype(ebitpix(filedatatype), bzero, bscale), fpixel, npixels_per_extension, NULL, readpointer, NULL, &status)) {
				throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			memcpy(basepointer, readpointer, extensionsize); 
			switch (datatype) {
				case tshort: 
				case short_img: {
					basepointer = (void *)((short *)basepointer + npixels_per_extension);
				}
					break;
				case tushort: {	
					basepointer = (void *)((unsigned short *)basepointer + npixels_per_extension);
				}
					break;
				case tfloat: 
				{
					basepointer = (void *)((float *)basepointer + npixels_per_extension);
				}
					break;
				case tdouble: 
				{	
					basepointer = (void *)((double *)basepointer + npixels_per_extension);
				}
					break;
				default:
					throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
					break;
			}
			/* select extension */
			if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
				throw operaException("operaFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
			}
			if ( fits_movrel_hdu(fptr, 1, &hdutype, &status) ) {
				throw operaException("operaFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
			}
		}
		free(readpointer);
	} else {
		if (fits_get_img_equivtype(fptr, &filedatatype, &status)) {
			throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		switch (todatatype(ebitpix(filedatatype), bzero, bscale)) {
			case tshort: {
				short *s_pixptr;
				s_pixptr = (short *)malloc(npixels * sizeof(short));		
				if (!s_pixptr) {
					throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
				}
				if (fits_read_pix(fptr, datatype, fpixel, npixels, NULL, s_pixptr, NULL, &status)) {
					throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				memcpy(pixptr, s_pixptr, size); 
				free(s_pixptr);
			}
				break;
			case tushort: {	
				unsigned short *us_pixptr = (unsigned short *)malloc(size);
				if (!us_pixptr) {
					throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
				}
				if (fits_read_pix(fptr, datatype, fpixel, npixels, NULL, us_pixptr, NULL, &status)) {
					throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				memcpy(pixptr, us_pixptr, size); 
				free(us_pixptr);
			}
				break;
			case tfloat: 
			{
				float *f_pixptr = (float *)malloc(size);
				if (!f_pixptr) {
					throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
				}
				if (fits_read_pix(fptr, datatype, fpixel, npixels, NULL, f_pixptr, NULL, &status)) {
					throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}				
				memcpy(pixptr, f_pixptr, size); 
				free(f_pixptr);
			}
				break;
			case tdouble: 
			{	
				double *d_pixptr = (double *)malloc(size);
				if (!d_pixptr) {
					throw operaException("operaFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
				}
				if (fits_read_pix(fptr, datatype, fpixel, npixels, NULL, d_pixptr, NULL, &status)) {
					throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}				
				memcpy(pixptr, d_pixptr, size); 
				free(d_pixptr);
			}
				break;
			default:
				throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
				break;
		}
	}
}
