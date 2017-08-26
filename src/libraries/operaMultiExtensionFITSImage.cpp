/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaMultiExtensionFITSImage
 Version: 1.0
 Description: class encapsulates a MEF FITS image.
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
#include "libraries/operaMultiExtensionFITSImage.h"
#include "libraries/operaMultiExtensionFITSCube.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSSubImage.h"

/*! \file operaMultiExtensionFITSImage.cpp */

using namespace std;

/*
 * \author Megan Tannock
 * \brief This class encapsulates multi extension FITS images.
 * \ingroup libraries
 */

/* 
 * operaMultiExtensionFITSImage()
 * \brief Basic operaMultiExtensionFITSImage constructor.
 * \return void
 */
operaMultiExtensionFITSImage::operaMultiExtensionFITSImage() : 
operaFITSImage()
{
}

/* 
 * \class operaMultiExtensionFITSImage
 * \brief construct an in-memory FITSImage object
 * \brief operaMultiExtensionFITSImage(unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype=tushort)
 * \brief Create an in-memory image of given dimensions.
 * \param Naxis1 - x ccd dimension
 * \param Naxis2 - y ccd dimension
 * \param Extensions - number of extensions
 * \param Datatype optional datatype defaults to tshort
 * \throws operaException operaErrorNoMemory
 * \return void
 */
operaMultiExtensionFITSImage::operaMultiExtensionFITSImage(unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype) : 
operaFITSImage(Naxis1, Naxis2, Datatype)
{
	mode = READWRITE;
    naxis = 2;
    naxis1 = Naxis1;
    naxis2 = Naxis2;
    naxis3 = 1;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	extensions = Extensions;
    npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
	npixels_per_slice = naxis1*naxis2;
	npixels_per_extension = naxis1*naxis2*naxis3;
	compression = cNone;
    datatype = Datatype;
	imageType = MEF;
    hdu = 0;	// signals the a header may need to be created...
	isLazy = false;
	
	bitpix = tobitpix(Datatype);
	if (pixptr == NULL) {
		size_t size = toSize(bitpix, npixels);
		pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
		memset(pixptr, 0, size);		
	}
	AllExtensions = true;
	AllSlices = true;
}

/* 
 * \class operaMultiExtensionFITSImage
 * \brief create a writeable file image of an in memory FITSImage object
 * \brief operaMultiExtensionFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype, unsigned Compression, bool isLazy)
 * \brief Constructor to create a new FITS file.
 * \param Filename to create (file is deleted if it exists)
 * \param Naxis1 - dimensions
 * \param Naxis2 - dimensions
 * \param Extensions - number of extensions
 * \param Datatype defaults to tshort
 * \param Compression, defaults to none
 * \param isLazy, defaults to true
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoMemory
 * \return void
 */
operaMultiExtensionFITSImage::operaMultiExtensionFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype, unsigned Compression, bool IsLazy) {
	int status = 0;
	unsigned x1 = 1;
	unsigned x2 = x1+naxis1;
	unsigned y1 = 1;
	unsigned y2 = y1+naxis2;
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int hdutype = ANY_HDU;
	
	filename = Filename;
	compression = Compression;
	mode = READWRITE;
	naxis = 2;
	naxis1 = Naxis1;
    naxis2 = Naxis2;
    naxis3 = 1;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	extensions = Extensions;
	npixels_per_slice = naxis1*naxis2;
	npixels_per_extension = naxis1*naxis2*naxis3;
	isLazy = IsLazy;
	if (isLazy) {
		npixels = npixels_per_extension;
	} else {
		npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
	}
	compression = Compression;
    datatype = Datatype;
	imageType = MEF;
    hdu = 0;	// signals the a header may need to be created...
    
	// remove existing file - cfitsio returns an error if it exists...
	remove(filename.c_str());
	
	bitpix = tobitpix(Datatype);
	
	// Open FITS file for output
    if (fits_create_file(&fptr, filename.c_str(), &status)) {
		throw operaException("operaMultiExtensionFITSImage: create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
	}
	if (fits_set_compression_type(fptr, compression, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	size_t size = toSize(bitpix, npixels);
	pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
	if (!pixptr) {
		throw operaException("operaMultiExtensionFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	memset(pixptr, 0, size);
    if (hdu == 0) {
        // Create headers (including primary) and images
		if (fits_create_img(fptr, bitpix, /*naxis=0*/0, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		//	get current hdunum
		fits_get_hdu_num(fptr, (int *)&hdu);
		//if (fits_set_hdrsize(fptr, 600, &status)) {
		//	throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		//}
		if (fits_write_key_lng(fptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		for (unsigned extension=1; extension<=extensions; extension++) {
            if (fits_create_img(fptr, bitpix, naxis, naxes, &status)) {
                throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
			fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
			status= 0;
#if 1
 			float bzero = 0.0, bscale = 1.0;
			if (datatype == tfloat) {
				if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {
				if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
			if (fits_update_key_flt(fptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(fptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
#endif
			fits_read_keyword(fptr, "DETSIZE", zvalue, zcomment, &status);
			if (status) {	// keyword does not exist, add the MEF base keywords
				status = 0;
				char buff[1024];
				sprintf(buff, "[%u:%u,%u:%ld]", x1, x2, y1, naxis2*extensions+1);
				// DETSIZE and DETSEC keywords both REQUIRED for file to open as a MEF.
				if (fits_update_key(fptr, tstring, "DETSIZE", buff, (char *)"Total data pixels in full mosaic", &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				sprintf(buff, "[%u:%u,%u:%u]", x1, x2, y1, y2);
				// DETSEC changes for each extension.
				if (fits_update_key(fptr, tstring, "DETSEC", buff, (char *)"Mosaic area of the detector", &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				y1 = y2+1;
				y2 += naxis2;
				
			}
		}
    }
	/* select extension 1 */
	if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if ( fits_movrel_hdu(fptr, 1, &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (!isLazy) {
		if (extensions > 0) {
			AllExtensions = true;
			AllSlices = true;
		}
	}    
}

/* 
 * \class operaMultiExtensionFITSImage
 * \brief create a multi extension FITS image object from a multi extension FITS file
 * \brief operaMultiExtensionFITSImage(string Filename, edatatype Datatype, int mode, unsigned Compression, bool isLazy)
 * \brief Constructor to create a MEF FITS Image from a MEF FITS file.
 * \param Filename
 * \param Datatype, defaults to tshort
 * \param mode, defaults to READWRITE (READWRITE | READONLY)
 * \param Compression, defaults to none
 * \param isLazy, defaults to true
 * \return void
 */
operaMultiExtensionFITSImage::operaMultiExtensionFITSImage(string Filename, edatatype Datatype, int mode/*=READWRITE|READONLY*/, unsigned Compression, bool IsLazy) : 
operaFITSImage(Filename, Datatype, mode, Compression, IsLazy)
{
	npixels_per_slice = naxis1*naxis2;
	npixels_per_extension = naxis1*naxis2*naxis3;
	imageType = MEF;
    if (IsLazy) {   // if lazy read...
        npixels = npixels_per_extension;
    } else { 
        npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
    }
}
/* 
 * operaMultiExtensionFITSImage(operaMultiExtensionFITSImage &imageIn, bool ViewOnly, bool AddHeader)
 * \brief Clone a Multi Extension FITSImage object.
 * \param imageIn - address of image to clone
 * \param ViewOnly, defaults to false
 * \param AddHeader, defaults to false
 * \return operaMultiExtensionFITSImage*
 */
operaMultiExtensionFITSImage::operaMultiExtensionFITSImage(operaMultiExtensionFITSImage &imageIn, bool ViewOnly, bool AddHeader) : operaFITSImage(imageIn, ViewOnly, AddHeader)
{
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
	
	filename = imageIn.filename;
	fptr = imageIn.fptr;
	isClone = true;				// i.e.do not close fptr
	datatype = imageIn.datatype;
	bzero = imageIn.bzero;
	bscale = imageIn.bscale;
	imageType = imageIn.imageType;
	bitpix = imageIn.bitpix;
	isLazy = imageIn.isLazy;
	AllExtensions = imageIn.AllExtensions;
	AllSlices = imageIn.AllSlices;
	
    if (isLazy) {
        for (unsigned ext=0; ext<=MAXFITSEXTENSIONS; ext++) {
            // Arrays have slots 0 -4, but there are only 4 extensions (numbered 1,2,3,4)
            // First slot (0) will be left as false and not used.
            extensionHasBeenRead[ext] = imageIn.extensionHasBeenRead[ext];
            extensionHasBeenWritten[ext] = imageIn.extensionHasBeenWritten[ext];
        }
    }
	// The memcpy may not be needed
	if (ViewOnly && imageIn.pixptr != NULL) {
		pixptr = imageIn.pixptr;
	} else {
		size_t size = imageIn.getsize();
        pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
		if (imageIn.getpixels()) {
			memcpy(pixptr, imageIn.getpixels(), size);
		}
		if (imageIn.pixptr == NULL) {
			imageIn.pixptr = pixptr;
		}
	}
    if (AddHeader){
        //
    } else {
        hdu = 0;
    }
}
/* 
 * operaMultiExtensionFITSImage(operaMultiExtensionFITSCube &imageIn, unsigned slice, bool AddHeader)
 * \brief Clone a Multi Extension FITSImage object.
 * \param imageIn - address of image to clone
 * \param slice, which slice
 * \return operaMultiExtensionFITSImage*
 */
operaMultiExtensionFITSImage::operaMultiExtensionFITSImage(operaMultiExtensionFITSCube &imageIn, unsigned slice, bool AddHeader)
{
	naxis = 2;
	naxis1 = imageIn.naxis1;
	naxis2 = imageIn.naxis2;
	naxis3 = 0;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	extensions = imageIn.extensions;
	npixels = naxis1*naxis2*extensions;
	npixels_per_slice = imageIn.npixels_per_slice;
	npixels_per_extension = naxis1*naxis2;
	current_extension = 1;
	current_slice = 1; 
    mode = MEF;
	
	filename = "";
	fptr = imageIn.fptr;
	isClone = true;				// i.e.do not close fptr
	datatype = imageIn.datatype;
	bzero = imageIn.bzero;
	bscale = imageIn.bscale;
	imageType = imageIn.imageType;
	bitpix = imageIn.bitpix;
	isLazy = false;
	AllExtensions = true;
	AllSlices = true;
	
    if (isLazy) {
        for (unsigned ext=0; ext<=MAXFITSEXTENSIONS; ext++) {
            // Arrays have slots 0 -4, but there are only 4 extensions (numbered 1,2,3,4)
            // First slot (0) will be left as false and not used.
            extensionHasBeenRead[ext] = imageIn.extensionHasBeenRead[ext];
            extensionHasBeenWritten[ext] = imageIn.extensionHasBeenWritten[ext];
        }
    }
	size_t size = getsize();
	pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
	if (imageIn.getpixels()) {
		float *basepointer = (float *)imageIn.getpixels() + (slice-1)*npixels_per_slice;
		float *newbasepointer = (float *)pixptr;
		for (unsigned extension = 1; extension <= extensions; extension++) {
			memcpy(newbasepointer, basepointer, size/extensions);
			basepointer += imageIn.npixels_per_extension;
			newbasepointer += npixels_per_extension;
		}
	}
    if (AddHeader){
        //
    } else {
        hdu = 0;
    }
}

operaMultiExtensionFITSImage::~operaMultiExtensionFITSImage() {
}

/* 
 * bool isInstantiated(operaMultiExtensionFITSImage &Image)
 * \brief Is this image in memory yet?
 * \param operaMultiExtensionFITSImage &Image
 * \return bool
 */
bool operaMultiExtensionFITSImage::isInstantiated(operaMultiExtensionFITSImage &Image) {
	return Image.getpixels() != NULL;
}
/* 
 * operaMultiExtensionFITSImageSave() 
 * \brief Saves the current image to disk.
 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
 * \throws operaException operaErrorCodeNoFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageSave() {
	int status = 0;
	unsigned x1 = 1;
	unsigned x2 = x1+naxis1;
	unsigned y1 = 1;
	unsigned y2 = y1+naxis2;
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int hdutype = ANY_HDU;
	fitsfile *newfptr;         // FITS file pointer for the new image
	const long fpixel = 1;
	const int morekeys = 0;
	
	// lazy writes have already happened...
	if (isLazy) {
		return;
	}
    // check if there is a filename 
	if (filename.empty() || filename == "")
		throw operaException("operaMultiExtensionFITSImage: ", operaErrorCodeNoFilename, __FILE__, __FUNCTION__, __LINE__);	
    // check if there is a filename 
	if (mode == READONLY)
		throw operaException("operaMultiExtensionFITSImage: "+filename+" ", operaErrorCodeChangeREADONLYError, __FILE__, __FUNCTION__, __LINE__);	
	// remove existing file - cfitsio returns an error if it exists...
	remove(filename.c_str());
    // creates a fits file with filename filename.c_str()
	if (fits_create_file(&newfptr, filename.c_str(), &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
    // if fptr == NULL, image was only created in memory
    if (fptr == NULL) {
        // Create headers (including primary) and images
		if (fits_create_img(newfptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
        // copy the current HDU from file associated with fptr and appends to the
        // end of the file associated with newfptr (in this case the end is the primary)
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}  
    }
	if (fits_update_key_lng(newfptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key_lng(newfptr, "NAXIS", (unsigned short)0, (char *)"Number of axes", &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	fits_delete_key(newfptr, "NAXIS1", &status);
	status = 0;
	fits_delete_key(newfptr, "NAXIS2", &status); 
	status = 0;
    // set the compression type for the file associated with newfptr
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	void *basepointer = pixptr;
	// save each extension with header
	for (unsigned extension=1; extension<=extensions; extension++) {
		// move to the specified HDU number  (select extension)
		if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if ( fits_movrel_hdu(fptr, extension, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
		// copy the current HDU from the file associated with fptr to the file associated with newfptr
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		fits_read_keyword(fptr, "DETSIZE", zvalue, zcomment, &status);
		if (status) {	// keyword does not exist, add the MEF base keywords
			status = 0;
			char buff[1024];
			sprintf(buff, "[%u:%u,%u:%ld]", x1, x2, y1, naxis2*extensions+1);
			// DETSIZE and DETSEC keywords both REQUIRED for file to open as a MEF.
			if (fits_update_key(newfptr, tstring, "DETSIZE", buff, (char *)"Total data pixels in full mosaic", &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			sprintf(buff, "[%u:%u,%u:%u]", x1, x2, y1, y2);
			// DETSEC changes for each extension.
			if (fits_update_key(newfptr, tstring, "DETSEC", buff, (char *)"Mosaic area of the detector", &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			y1 = y2+1;
			y2 += naxis2;
			
		}
		status = 0;
#if 0
		// be sure we have the basic headers for floats if no header existed
		if (datatype == tfloat || datatype == tdouble) {
			float bzero = 0.0, bscale = 1.0;
			if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
#endif
		if (AllExtensions) {
			// write elements into the fits file
			if (fits_write_img(newfptr, datatype, fpixel, npixels_per_extension, basepointer, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			// move the basepointer to next extension
			switch (datatype) {
				case tshort: 
				case short_img:
				{
					basepointer = (void *)((short *)basepointer + npixels_per_extension);
				}
					break;
				case tushort:
				{	
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
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
					break;
			}
		}
	}
	fits_close_file(newfptr, &status);
	newfptr = NULL;
	// note that fptr is left open in case the user wants to save again, it will be closed
    // by a close call if made, or the destructor otherwise
}
/* 
 * operaMultiExtensionFITSImageSaveAs(string newFilename) 
 * \brief Saves the current image to disk, with the given filename.
 * \param newFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageSaveAs(string newFilename) {
	int status = 0;
	unsigned x1 = 1;
	unsigned x2 = x1+naxis1;
	unsigned y1 = 1;
	unsigned y2 = y1+naxis2;
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int hdutype = ANY_HDU;
	fitsfile *newfptr;         // FITS file pointer for the new image
	const long fpixel = 1;
	
	// remove existing file - cfitsio returns an error if it exists...
	remove(newFilename.c_str());
    // creates a fits file with filename newFilename.c_str()
	if (fits_create_file(&newfptr, newFilename.c_str(), &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio create error: "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
    // if fptr == NULL, image was only created in memory
    if (fptr == NULL) {
        // Create headers (including primary) and images
		if (fits_create_img(newfptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
        // copy the current HDU from file associated with fptr and appends to the
        // end of the file associated with newfptr (in this case the end is the primary)
		if (fits_create_img(newfptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(newfptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
#if 1
		/* copy all the user keywords (not the structural keywords) */
		int nkeys = 0;
		char card[81];
		if (fits_get_hdrspace(fptr, &nkeys, NULL, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
		for (int ii = 1; ii <= nkeys; ii++) {
			if (fits_read_record(fptr, ii, card, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
				if (fits_write_record(newfptr, card, &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
		}
#endif	
	}
	if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
#if 0
	// be sure we have the basic headers for floats if no header existed
	if (datatype == tfloat || datatype == tdouble) {
		//if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
		//	throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		//}
		if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
#endif
	// it is OK if the keyword isn't there i.e. fptr was NULL
	if (fptr != NULL) {
		//if (fits_delete_key(newfptr, "NAXIS3", &status)) {
		//	throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		//}
	}
    // set the compression type for the file associated with newfptr
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}    
	void *basepointer = pixptr;
	// save each extension with header
	for (unsigned extension=1; extension<=extensions; extension++) {
		// move to the specified HDU number  (select extension)
		if (fptr != NULL) {
            if (fits_create_img(newfptr, bitpix, naxis, naxes, &status)) {
                throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if ( fits_movrel_hdu(fptr, extension, &hdutype, &status) ) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
			}
			// copy the current HDU from the file associated with fptr to the file associated with newfptr
#if 1
			/* copy all the user keywords (not the structural keywords) */
			int nkeys = 0;
			char card[81];
			if (fits_get_hdrspace(fptr, &nkeys, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			} 
			for (int ii = 1; ii <= nkeys; ii++) {
				if (fits_read_record(fptr, ii, card, &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
					if (fits_write_record(newfptr, card, &status)) {
						throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
					}
				}
			}
#endif	
			if (fits_update_key_lng(newfptr, "NAXIS", (unsigned short)naxis, (char *)"Number of axes", &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		} else {
            if (fits_create_img(newfptr, bitpix, naxis, naxes, &status)) {
                throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
			fits_movrel_hdu(newfptr, 1, NULL, &status);  /* try to move to next HDU */
			status= 0;
			if (fits_update_key_lng(newfptr, "NAXIS", (unsigned short)naxis, (char *)"Number of axes", &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_lng(newfptr, "NAXIS1", (unsigned short)naxis1, (char *)"Length of data axis 1", &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_lng(newfptr, "NAXIS2", (unsigned short)naxis2, (char *)"Length of data axis 2", &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		// be sure we have the basic headers for floats if no header existed
		if (datatype == tfloat || datatype == tdouble) {
			//if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
			//	throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			//}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
		if (fptr != NULL) {
			fits_read_keyword(fptr, "DETSIZE", zvalue, zcomment, &status);
		}
		if (status != 0 || fptr == NULL) {	// keyword does not exist, add the MEF base keywords
			status = 0;
			char buff[1024];
			sprintf(buff, "[%u:%u,%u:%ld]", x1, x2, y1, naxis2*extensions+1);
			// DETSIZE and DETSEC keywords both REQUIRED for file to open as a MEF.
			if (fits_update_key(newfptr, tstring, "DETSIZE", buff, (char *)"Total data pixels in full mosaic", &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			sprintf(buff, "[%u:%u,%u:%u]", x1, x2, y1, y2);
			// DETSEC changes for each extension.
			if (fits_update_key(newfptr, tstring, "DETSEC", buff, (char *)"Mosaic area of the detector", &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			y1 = y2+1;
			y2 += naxis2;
			
		}
		status = 0;
		if (AllExtensions) {
			// write elements into the fits file
			if (fits_write_img(newfptr, datatype, fpixel, npixels_per_extension, basepointer, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
			}
			// move the basepointer to next extension
			switch (datatype) {
				case tshort: 
				case short_img:
				{
					basepointer = (void *)((short *)basepointer + npixels_per_extension);
				}
					break;
				case tushort:
				{	
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
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+newFilename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
					break;
			}
		}
	}
	fits_close_file(newfptr, &status);
	newfptr = NULL;
	// note that fptr is left open in case the user wants to save again, it will be closed
    // by a close call if made, or the destructor otherwise
}
/*
 * void addExtension(string name, unsigned Columns, unsigned Rows)
 * \brief adds an extension 
 * \param unsigned Columns
 * \param unsigned Rows
 */
void operaMultiExtensionFITSImage::addExtension(string name, unsigned Columns, unsigned Rows, DATASEC_t &detsec, bool reuse) {
	int status = 0;
	int hdutype = IMAGE_HDU;

	if (!reuse) {
		extensions++;
		naxis = 2;
		naxis1 = Columns;
		naxis2 = Rows;
		naxis3 = 0;
		naxes[0] = naxis1;
		naxes[1] = naxis2;
		naxes[2] = naxis3;
		npixels_per_extension = npixels_per_slice = npixels = naxis1 * naxis2;
		/*
		 * we always need a new pixptr since the size of each extension may differ
		 */
		if (pixptr) {
			free(pixptr);
		}
		size_t size = toSize(bitpix, npixels);
		pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
		if (!pixptr) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
		}
		memset(pixptr, 0, size);
		if (fits_create_img(fptr, bitpix, naxis, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	if (fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key_lng(fptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_movrel_hdu(fptr, extensions, NULL, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (datatype == tfloat) {
		bzero = 0.0;
		bscale = 1.0;
	}
	if (datatype == tfloat) {
		if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
		if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	if (fits_update_key_flt(fptr, "BZERO", bzero, -1, NULL, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key_flt(fptr, "BSCALE", bscale, -1, NULL, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key(fptr, tstring, "EXTNAME", (void *)name.c_str(), (char *)"Extension name", &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned x1 = detsec.x1;
	unsigned x2 = detsec.x2;
	unsigned y1 = detsec.y1;
	unsigned y2 = detsec.y2;
	char buff[80];
	sprintf(buff, "[%u:%u,%u:%u]", x1, x2, y1, y2);
	// DETSEC changes for each extension.
	if (fits_update_key(fptr, tstring, "DETSEC", buff, (char *)"Mosaic area of the detector", &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	x1 = 1;
	x2 = detsec.x2;
	y1 = 1;
	y2 = detsec.y2;
	sprintf(buff, "[%u:%u,%u:%u]", x1, x2, y1, y2);
	// DETSIZE and DETSEC keywords both REQUIRED for file to open as a MEF.
	if (fits_update_key(fptr, tstring, "DETSIZE", buff, (char *)"Total data pixels in full mosaic", &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
}
/* 
 * void setExtension(unsigned extension)
 * \brief sets the current extension
 * \param unsigned new_extension
 */
void operaMultiExtensionFITSImage::setExtension(unsigned new_extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (new_extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	int status = 0;
	int hdutype = ANY_HDU;
	
	if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
	}
	if ( fits_movrel_hdu(fptr, new_extension, &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
	}
	current_extension = new_extension;
}
/* 
 * void getExtension()
 * \brief get the current extension
 * \param unsigned extension
 */
int operaMultiExtensionFITSImage::getExtension(){
    return current_extension;
}

/* 
 * operaMultiExtensionFITSImageCopyHeader(operaMultiExtensionFITSImage *from, unsigned extension) 
 * \brief Copies all of the header information from image. 
 * \param from
 * \param extension
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageCopyHeader(operaMultiExtensionFITSImage *from) {
	int status = 0;
	const int morekeys = 0;
	int hdutype = ANY_HDU;
	
	if (!fptr) {
		throw operaException("operaMultiExtensionFITSImage: ", operaErrorCodeNULLfptr, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (!from->getfitsfileptr()) {
		throw operaException("operaMultiExtensionFITSImage: ", operaErrorCodeNULLfptr, __FILE__, __FUNCTION__, __LINE__);	
	}
	// if hdu == NULL, image was only created in memory
	if (hdu == 0) {
		// Create header (primary)
		if (fits_create_img(fptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		//	get current hdunum
		fits_get_hdu_num(fptr, (int *)&hdu);
	}
#if 0
	/* copy all the user keywords (not the structural keywords) */
	int nkeys = 0;
	char card[81];
	if (fits_get_hdrspace(from->getfitsfileptr(), &nkeys, NULL, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	} 
	for (int ii = 1; ii <= nkeys; ii++) {
		if (fits_read_record(from->getfitsfileptr(), ii, card, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
			if (fits_write_record(fptr, card, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
#endif	
	for (unsigned ext=0; ext<=extensions; ext++) {
		if ( fits_delete_hdu(fptr, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	for (unsigned ext=0; ext<=extensions; ext++) {
		if ( fits_movabs_hdu(from->getfitsfileptr(), cfitsioextension(0), &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		if ( fits_movrel_hdu(from->getfitsfileptr(), ext, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		// copy the current HDU from file associated with from and appends to the
		// end of the file associated with fptr (in this case the end is the primary)
		if (fits_copy_hdu(from->getfitsfileptr(), fptr, morekeys, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (datatype == tfloat || datatype == tdouble) {
			float bzero = 0.0, bscale = 1.0;
			//if (fits_set_bscale(fptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
			//	throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			//}
			if (fits_update_key_flt(fptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(fptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if ( fits_movrel_hdu(fptr, 1, &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	// special case of a single extension only, with a primary header...
	if (from->getZDimension() == 1 && ( from->getNExtensions() == 2 || from->getNExtensions() == 1 ) ) {
		// keep header info....
	} else {
		if ( fits_delete_hdu(fptr, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}
/* 
 * operaMultiExtensionFITSImageCopyNonStructuralHeader(operaMultiExtensionFITSImage *from, unsigned extension) 
 * \brief Copies all of the non-structural header information from image. 
 * \param from
 * \param extension
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageCopyNonStructuralHeader(operaMultiExtensionFITSImage *from) {
	int status = 0;
	int hdutype = ANY_HDU;
	int nkeys = 0;
	char card[81];
	
	if (!fptr) {
		throw operaException("operaMultiExtensionFITSImage: ", operaErrorCodeNULLfptr, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (!from->getfitsfileptr()) {
		throw operaException("operaMultiExtensionFITSImage: ", operaErrorCodeNULLfptr, __FILE__, __FUNCTION__, __LINE__);	
	}
	/* copy all the user keywords (not the structural keywords) */
	if (fits_get_hdrspace(from->getfitsfileptr(), &nkeys, NULL, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	} 
	for (int ii = 1; ii <= nkeys; ii++) {
		if (fits_read_record(from->getfitsfileptr(), ii, card, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
			if (fits_write_record(fptr, card, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	for (unsigned ext=0; ext<=extensions; ext++) {
		if ( fits_delete_hdu(fptr, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	for (unsigned ext=0; ext<=extensions; ext++) {
		if ( fits_movabs_hdu(from->getfitsfileptr(), cfitsioextension(0), &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		if ( fits_movrel_hdu(from->getfitsfileptr(), ext, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (fits_get_hdrspace(from->getfitsfileptr(), &nkeys, NULL, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
		for (int ii = 1; ii <= nkeys; ii++) {
			if (fits_read_record(from->getfitsfileptr(), ii, card, &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
				if (fits_write_record(fptr, card, &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
		}
	}
	if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if ( fits_movrel_hdu(fptr, 1, &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if ( fits_delete_hdu(fptr, &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
}
/* 
 * void readExtension(unsigned extension)
 * \brief read an extension into memory
 * \param operaMultiExtension FITS Image &Image 
 * \param unsigned extension
 */
void operaMultiExtensionFITSImage::readExtension(unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif    
	int status = 0;
	int hdutype = ANY_HDU;
	long fpixel[MAXFITSDIMENSIONS] = {1,1,1};
    int filedatatype;
    
    // Must be a lazy read to use this method
    if (!isLazy) {
		throw operaException("operaMultiExtensionFITSImage: ", operaErrorNotLazy, __FILE__, __FUNCTION__, __LINE__);	
	}
	// Bail if we have already read the extension -> it is already in memory...
	if (extensionHasBeenRead[extension] && extension == current_extension) {
		return;
	}
	setHasBeenRead(extension);
    if (super) {
        super->setHasBeenRead(extension);
	}
    current_extension = extension;
	AllExtensions = false;
    // read from original file
	// Note that different extensions may have different dimensions!
    // select extension 
    if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
        throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
    }
    if ( fits_movrel_hdu(fptr, current_extension, &hdutype, &status) ) {
        throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
    }
	if (fits_get_img_size(fptr, MAXFITSDIMENSIONS, naxes, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}	
    if (pixptr == NULL || naxis1 != (unsigned long)naxes[0] || naxis2 != (unsigned long)naxes[1]) {
		naxis1 = naxes[0];
		naxis2 = naxes[1];
		naxis3 = ((naxis==2)?1:naxes[2]);
		npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
		npixels_per_slice = naxis1*naxis2;
		npixels_per_extension = naxis1*naxis2*naxis3;
        size_t size = toSize(bitpix, npixels);
		if (pixptr) {
			free(pixptr);
		}
        pixptr = malloc(MAX(size, toSize(float_img, npixels_per_extension))); 
        memset(pixptr, 0, size);		
		if (super) {
			super->pixptr = pixptr;
		}
    }
    if (fits_get_img_equivtype(fptr, &filedatatype, &status)) {
        throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
    }
	float mybzero = bzero;
	if (filedatatype == USHORT_IMG) {
		mybzero = 32768.0;
	}
    if (fits_read_pix(fptr, todatatype(ebitpix(filedatatype), mybzero, bscale), fpixel, npixels_per_extension, NULL, pixptr, NULL, &status)) {
        throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
    }
    if (todatatype(ebitpix(filedatatype), mybzero, bscale) != datatype) {
        operaFITSImageConvertImageInPlace(todatatype(ebitpix(filedatatype), bzero, bscale), datatype);
		if (datatype == tfloat) {
			bzero = 0.0;
			bscale = 1.0;
			if (super) {
				super->bzero = bzero;
				super->bscale = bscale;
			}
		}
    }
}
/* 
 * void saveExtension(unsigned extension)
 * \brief save image data
 * \param unsigned extension
 */ 
void operaMultiExtensionFITSImage::saveExtension(unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    int status = 0;
	long fpixel = 1;
    
    // Must be READWRITE mode to use this method
    if (mode == READONLY) {
		throw operaException("operaMultiExtensionFITSImage: ", operaErrorCodeChangeREADONLYError, __FILE__, __FUNCTION__, __LINE__);	
	}
    // Must be a lazy read to use this method
    if (isLazy) {
		if (current_extension != extension) {
			setExtension(extension);
		}
		setHasBeenWritten(extension);
		if (super) {
			super->setHasBeenWritten(extension);
		}
		if (fits_write_img(fptr, datatype, fpixel, npixels_per_extension, pixptr, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}			
	} else {	// just write it....
		if (fits_write_img(fptr, datatype, fpixel, npixels_per_extension, pixptr, &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}			
	}
}
/* 
 * string operaFITSGetHeaderValue(string keyword) 
 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
 * \param keyword
 * \return string keyword value or empty string
 */
string operaMultiExtensionFITSImage::operaFITSGetHeaderValue(string keyword) {
	string raw = operaFITSGetRawHeaderValue(keyword);
	return trimFITSKeyword(raw.c_str());
}
/* 
 * float operaFITSGetFloatHeaderValue(string keyword) 
 * \brief returns the float value of a FITS keyword.
 * \param keyword
 * \return float keyword value 
 */
float operaMultiExtensionFITSImage::operaFITSGetFloatHeaderValue(string keyword) {
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
int operaMultiExtensionFITSImage::operaFITSGetIntHeaderValue(string keyword) {
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
string operaMultiExtensionFITSImage::operaFITSGetRawHeaderValue(string keyword) {
	char zvalue[FLEN_VALUE],comment[FLEN_COMMENT];
	int status = 0;
	if (fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, comment, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" "+keyword+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	return string(zvalue);
}

/* 
 * string operaFITSGetHeaderValue(string keyword, unsigned extension) 
 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
 * \param keyword
 * \param extension
 * \return string keyword value or empty string
 */
string operaMultiExtensionFITSImage::operaFITSGetHeaderValue(string keyword, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    setExtension(extension);
	string raw = operaFITSGetRawHeaderValue(keyword, extension);
	return trimFITSKeyword(raw.c_str());
}
/* 
 * float operaFITSGetFloatHeaderValue(string keyword, unsigned extension) 
 * \brief returns the float value of a FITS keyword.
 * \param keyword
 * \param extension
 * \return float keyword value 
 */
float operaMultiExtensionFITSImage::operaFITSGetFloatHeaderValue(string keyword, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    setExtension(extension);
	string raw = operaFITSGetRawHeaderValue(keyword, extension);
	string clean = trimFITSKeyword(raw.c_str());
	return atof(clean.c_str());
}
/* 
 * int operaFITSGetIntHeaderValue(string keyword, unsigned extension) 
 * \brief returns the int value of a FITS keyword.
 * \param keyword
 * \param extension
 * \return int keyword value
 */
int operaMultiExtensionFITSImage::operaFITSGetIntHeaderValue(string keyword, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    setExtension(extension);
	string raw = operaFITSGetRawHeaderValue(keyword, extension);
	string clean = trimFITSKeyword(raw.c_str());
	return atoi(clean.c_str());
}
/* 
 * string operaFITSGetRawHeaderValue(string keyword, unsigned extension) 
 * \brief returns the value of a FITS keyword verbatim.
 * \param keyword
 * \param extension
 * \throws operaException cfitsio error code
 * \return string keyword value or empty string
 */
string operaMultiExtensionFITSImage::operaFITSGetRawHeaderValue(string keyword, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    setExtension(extension);
	char zvalue[FLEN_VALUE],comment[FLEN_COMMENT];
	int status = 0;
	if (fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, comment, &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" "+keyword+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	return string(zvalue);
}
/* 
 * operaFITSSetHeaderValue(string keyword, string value, string comment, unsigned comment)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, string value, string comment, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	setExtension(extension);
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tstring, (char*)keyword.c_str(), (char*)value.c_str(), (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tstring, (char*)keyword.c_str(), (char*)value.c_str(), (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}

/* 
 * operaFITSSetHeaderValue(string keyword, unsigned short value, string comment, unsigned extension)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \return void
 */
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, unsigned short value, string comment, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	setExtension(extension);
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tushort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tushort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}
/* 
 * operaFITSSetHeaderValue(string keyword, short value, string comment, unsigned extension)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \return void
 */
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, short value, string comment, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	setExtension(extension);
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tshort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tshort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}
/* 
 * operaFITSSetHeaderValue(string keyword, float value, string comment, unsigned extension)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, float value, string comment, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	setExtension(extension);
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tfloat, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tfloat, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}
/* 
 * operaFITSSetHeaderValue(string keyword, double value, string comment, unsigned extension)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \param value
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, double value, string comment, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	setExtension(extension);
	fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
	if (status) {	// keyword does not exist
		status = 0;
		if (fits_write_key(fptr, tdouble, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {		// keyword exists
		if (fits_update_key(fptr, tdouble, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
			throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}

/* 
 * operaFITSDeleteHeaderKey(string keyword, unsigned extension)
 * \brief sets the given keyword to value with comment.
 * \param keyword
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaFITSDeleteHeaderKey(string keyword, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	int status = 0;
	
	setExtension(extension);
	fits_delete_key(fptr, (char*)keyword.c_str(), &status);
}
/* 
 * operaFITSAddComment(string comment, unsigned extension)
 * \brief adds a comment to header.
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaFITSAddComment(string comment, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	int status = 0;
	
	setExtension(extension);
	if (fits_write_comment(fptr, (char*)comment.c_str(), &status)) {
		throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
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
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, string value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	if (AllExtensions) {
		for (unsigned extension=1; extension<=extensions; extension++) {
			setExtension(extension);
			fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
			if (status) {	// keyword does not exist
				status = 0;
				if (fits_write_key(fptr, tstring, (char*)keyword.c_str(), (char*)value.c_str(), (char*)comment.c_str(), &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {		// keyword exists
				if (fits_update_key(fptr, tstring, (char*)keyword.c_str(), (char*)value.c_str(), (char*)comment.c_str(), &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
		}
	} else {
		fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
		if (status) {	// keyword does not exist
			status = 0;
			if (fits_write_key(fptr, tstring, (char*)keyword.c_str(), (char*)value.c_str(), (char*)comment.c_str(), &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		} else {		// keyword exists
			if (fits_update_key(fptr, tstring, (char*)keyword.c_str(), (char*)value.c_str(), (char*)comment.c_str(), &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
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
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, unsigned short value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	if (AllExtensions) {
		for (unsigned extension=1; extension<=extensions; extension++) {
			setExtension(extension);
			fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
			if (status) {	// keyword does not exist
				status = 0;
                if (fits_write_key(fptr, tushort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {		// keyword exists
                if (fits_update_key(fptr, tushort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
		}
	} else {
        fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
		if (status) {	// keyword does not exist
			status = 0;
            if (fits_write_key(fptr, tushort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		} else {		// keyword exists
            if (fits_update_key(fptr, tushort, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
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
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, float value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	if (AllExtensions) {
		for (unsigned extension=1; extension<=extensions; extension++) {
			setExtension(extension);
            fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
			if (status) {	// keyword does not exist
				status = 0;
                if (fits_write_key(fptr, tfloat, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {		// keyword exists
                if (fits_update_key(fptr, tfloat, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
		}
	} else {
        fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
		if (status) {	// keyword does not exist
			status = 0;
            if (fits_write_key(fptr, tfloat, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		} else {		// keyword exists
            if (fits_update_key(fptr, tfloat, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
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
void operaMultiExtensionFITSImage::operaFITSSetHeaderValue(string keyword, double value, string comment) {
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int status = 0;
	
	if (AllExtensions) {
		for (unsigned extension=1; extension<=extensions; extension++) {
			setExtension(extension);
            fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
			if (status) {	// keyword does not exist
				status = 0;
                if (fits_write_key(fptr, tdouble, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {		// keyword exists
                if (fits_update_key(fptr, tdouble, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
					throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
		}
	} else {
        fits_read_keyword(fptr, (char*)keyword.c_str(), zvalue, zcomment, &status);
		if (status) {	// keyword does not exist
			status = 0;
            if (fits_write_key(fptr, tdouble, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		} else {		// keyword exists
            if (fits_update_key(fptr, tdouble, (char*)keyword.c_str(), &value, (char*)comment.c_str(), &status)) {
				throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
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
void operaMultiExtensionFITSImage::operaFITSDeleteHeaderKey(string keyword) {
	int status = 0;
	
	if (AllExtensions) {
		for (unsigned extension=1; extension<=extensions; extension++) {
			setExtension(extension);
            fits_delete_key(fptr, (char*)keyword.c_str(), &status);
		}
	} else {
        fits_delete_key(fptr, (char*)keyword.c_str(), &status);
	}
}
/* 
 * operaFITSAddComment(string comment)
 * \brief adds a comment to header.
 * \param comment
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSImage::operaFITSAddComment(string comment) {
	int status = 0;
	
	if (AllExtensions) {
		for (unsigned extension=1; extension<=extensions; extension++) {
			setExtension(extension);
            if (fits_write_comment(fptr, (char*)comment.c_str(), &status)) {
                throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
		}
	} else {
        if (fits_write_comment(fptr, (char*)comment.c_str(), &status)) {
            throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
	}
}
/*
 * bool isMEF()
 * \brief Is this image a MEF?
 * \return bool
 */    
bool operaMultiExtensionFITSImage::isMEF(void){
    return extensions > 1;
}

/* 
 * void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageWriteVirtual(operaFITSSubImage &subImage, unsigned X, unsigned Y, unsigned extension)
 * \brief Write a subimage on to avirtual operaMultiExtensionFITSImage.
 * \brief This is used for deep stacking. 
 * \brief The virtual image is Lazy, and can be very large.
 * \brief The actual image exists only on disk.
 * \param subImage the subImage
 * \param X - beginning x location in image
 * \param Y - beginning y location in image
 * \return void
 */
void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageWriteVirtual(operaFITSSubImage &subImage, unsigned long X, unsigned long Y, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    current_extension = extension;
    setExtension(current_extension); // calls fits_movabs_hdu
    
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
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
	} else {
		throw operaException("operaMultiExtensionFITSImage: "+filename+" ", operaErrorNotLazy, __FILE__, __FUNCTION__, __LINE__);	
	}
}

/* 
 * void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageReadVirtual(operaFITSSubImage &subImage, unsigned X, unsigned Y, unsigned extension)
 * \brief Read a subimage from a virtual operaMultiExtensionFITSImage.
 * \brief This is used for deep stacking. 
 * \brief The virtual image is Lazy, and can be very large.
 * \brief The actual image exists only on disk.
 * \param subImage the subImage
 * \param X - beginning x location in image
 * \param Y - beginning y location in image
 * \return void
 */
void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageReadVirtual(operaFITSSubImage &subImage, unsigned long X, unsigned long Y, unsigned extension) {
#ifdef RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSImage: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    current_extension = extension;
    setExtension(current_extension); // calls fits_movabs_hdu
    
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
			throw operaException("operaMultiExtensionFITSImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	} else {
		throw operaException("operaMultiExtensionFITSImage: "+filename+" ", operaErrorNotLazy, __FILE__, __FUNCTION__, __LINE__);	
	}
}
#if 0
/* 
 * Thread Support to process all extensions in parallel
 */
#include <pthread.h>

typedef struct thread_args {
	operaMultiExtensionFITSImage &a;		// input value LHS
	operaMultiExtensionFITSImage &b;		// input value RHS
	operaMultiExtensionFITSImage &c;		// output value
	unsigned extension;						// which extension is this
	ops_t op;								// the op to perform
} thread_args_t;

void *paropthunkthread(void *argument) {
	thread_args_t *thread_args_s = (thread_args_t *)argument;
	operaMultiExtensionFITSImage &a = thread_args_s->a;
	operaMultiExtensionFITSImage &b = thread_args_s->b;
	operaMultiExtensionFITSImage &c = thread_args_s->c;
	unsigned extension = thread_args_s->extension;
	a.setExtension(extension);
	b.setExtension(extension);
	c.setExtension(extension);
	switch (thread_args_s->op) {
		case op_mul:
			c = a * b;
			break;
		case op_div:
			c = a / b;
			break;
		case op_add:
			c = a + b;
			break;
		case op_sub:
			c = a - b;
			break;
		default:
			break;
	}
	return NULL;
}
operaMultiExtensionFITSImage::operaMultiExtensionFITSImage &par(operaMultiExtensionFITSImage &a, operaMultiExtensionFITSImage &b, operaMultiExtensionFITSImage &c, ops_t Op) {
	unsigned long nthreads = b.getNExtensions();
	// set up the thread vector
	pthread_t *threads = (pthread_t *)calloc(nthreads, sizeof(pthread_t*));
	// set up the argument vector
	thread_args_t *thread_args = (thread_args_t *)calloc(nthreads, sizeof(thread_args_t));
	// create all threads
	for (unsigned extension=1; extension<=nthreads; extension++) {
		thread_args[extension].op = Op;
		thread_args[extension].extension = extension;
		thread_args[extension].a = a;
		thread_args[extension].b = b;
		thread_args[extension].c = c;
		if (pthread_create(&threads[extension], NULL, paropthunkthread, (void *) &thread_args[extension]) != 0)
			return c;
	}
	// wait for all threads to complete
	for (unsigned extension=1; extension<=nthreads; extension++) {
		if (pthread_join(threads[extension], NULL) != 0) {
			return c;
		}
	}
	return c;
}
#endif

