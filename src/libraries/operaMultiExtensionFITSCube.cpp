/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaMultiExtensionFITSCube 
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

#include <pthread.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaImage.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMultiExtensionFITSImage.h"
#include "libraries/operaMultiExtensionFITSCube.h"
#include "libraries/operaException.h"
#include "libraries/operaStats.h"

/*! \file operaMultiExtensionFITSCube.cpp */

using namespace std;

/*! 
 * operaMultiExtensionFITSCube
 * \author Megan Tannock
 * \brief This class encapsulates the MEF FITS image.
 * \ingroup libraries
 */

/* 
 * operaMultiExtensionFITSCube()
 * \brief Basic operaMultiExtensionFITSCube constructor.
 * \return void
 */
operaMultiExtensionFITSCube::operaMultiExtensionFITSCube() : operaMultiExtensionFITSImage()
{
    naxis = 3;
}
/* 
 * \class operaMultiExtensionFITSCube
 * \brief construct an in-memory FITSImage object
 * \brief operaMultiExtensionFITSCube(unsigned Naxis1, unsigned Naxis2, edatatype Datatype=tushort, unsigned Compression=0)
 * \brief Create an in-memory image of given dimensions.
 * \param Naxis1 - x ccd dimension
 * \param Naxis2 - y ccd dimension
 * \param Datatype optional datatype defaults to tshort
 * \param Compression optional compression, defaults to none
 * \return void
 */
operaMultiExtensionFITSCube::operaMultiExtensionFITSCube(unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, unsigned Extensions, edatatype Datatype)
{
	mode = READWRITE;
    naxis = 3;
    naxis1 = Naxis1;
    naxis2 = Naxis2;
    naxis3 = Naxis3;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	extensions = Extensions;
	if (naxis3 == 1) {			// it is not really a cube...
		naxis = 2;
	}
    npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
	npixels_per_slice = naxis1*naxis2;
	npixels_per_extension = naxis1*naxis2*naxis3;
	compression = cNone;
    datatype = Datatype;
	imageType = MEFCube;
    hdu = 0;	// signals the a header may need to be created...
	isLazy = false;
	
	bitpix = tobitpix(Datatype);
	size_t size = toSize(bitpix, npixels);
	pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
	memset(pixptr, 0, size);		
	AllExtensions = true;
	AllSlices = true;    
}
/* 
 * \class operaMultiExtensionFITSCube
 * \brief create a writeable file image of an in memory FITSImage object
 * \brief operaMultiExtensionFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype, unsigned Compression)
 * \brief Constructor to create a new FITS file.
 * \param Filename to create (file is deleted if it exists)
 * \param Naxis1 - dimensions
 * \param Naxis2 - dimensions
 * \param Datatype defaults to tshort
 * \param Compression, defaults to none
 * \throws operaException cfitsio error code
 * \return void
 */
operaMultiExtensionFITSCube::operaMultiExtensionFITSCube(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, unsigned Extensions, edatatype Datatype, unsigned Compression, bool IsLazy) 
{
    int status = 0;
	unsigned x1 = 1;
	unsigned x2 = x1+naxis1;
	unsigned y1 = 1;
	unsigned y2 = y1+naxis2;
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
    int hdutype = ANY_HDU;
 	const long fpixel = 1;
   
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
	extensions = Extensions;
	if (naxis3 == 1) {			// it is not really a cube...
		naxis = 2;
	}
    if (isLazy) {   // if lazy read...
        npixels = npixels_per_extension;
    } else { 
        npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
    }
    npixels_per_slice = naxis1*naxis2;
    npixels_per_extension = naxis1*naxis2*naxis3;
    compression = Compression;
    datatype = Datatype;
	imageType = MEFCube;
    hdu = 0;	// signals the a header may need to be created...
    isLazy = IsLazy;
    
    // remove existing file - cfitsio returns an error if it exists...
    remove(filename.c_str());
    
    bitpix = tobitpix(Datatype);
    
    // Open FITS file for output
    if (fits_create_file(&fptr, filename.c_str(), &status)) {
        throw operaException("operaMultiExtensionFITSCube: create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
    }
    if (fits_set_compression_type(fptr, compression, &status)) {
        throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
    }
    size_t size = toSize(bitpix, npixels);
	pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
    if (!pixptr) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
    }
    memset(pixptr, 0, size);
    if (hdu == 0) {
        // Create headers (including primary) and images
        if (fits_create_img(fptr, bitpix, 0, naxes, &status)) {
            throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        /* select extension */
        if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
            throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_update_key_lng(fptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
            throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
		for (unsigned extension=1; extension<=extensions; extension++) {
 			float bzero = 0.0, bscale = 1.0;
            if (fits_create_img(fptr, bitpix, naxis, naxes, &status)) {
                throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
			fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
			status= 0;
			if (datatype == tfloat) {
				if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*4 (complex, stored as float)", &status)) {
					throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			} else {
				if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"Real*8 (complex, stored as float)", &status)) {
					throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
			}
			if (fits_update_key_flt(fptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(fptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			fits_read_keyword(fptr, "DETSIZE", zvalue, zcomment, &status);
			if (status) {	// keyword does not exist, add the MEF base keywords
				status = 0;
				char buff[1024];
				sprintf(buff, "[%u:%u,%u:%ld]", x1, x2, y1, naxis2*extensions+1);
				// DETSIZE and DETSEC keywords both REQUIRED for file to open as a MEF.
				if (fits_update_key(fptr, tstring, "DETSIZE", buff, (char *)"Total data pixels in full mosaic", &status)) {
					throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				sprintf(buff, "[%u:%u,%u:%u]", x1, x2, y1, y2);
				// DETSEC changes for each extension.
				if (fits_update_key(fptr, tstring, "DETSEC", buff, (char *)"Mosaic area of the detector", &status)) {
					throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				y1 = y2+1;
				y2 += naxis2;
				
			}
			//for (unsigned i=0; i<500; i++) {
			//	if (fits_write_comment(fptr, "COMMENT  Reserved space.  This line can be used to add a new FITS card.", &status)) {
			//		throw operaException("operaMultiExtensionFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			//	}			
			//}
		}
    }
	if (fits_write_img(fptr, datatype, fpixel, npixels, pixptr, &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	/* select extension */
	if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if ( fits_movrel_hdu(fptr, 1, &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (!isLazy) {
		if (extensions > 0) {
			AllExtensions = true;
			AllSlices = true;
		}
	}    
}

operaMultiExtensionFITSCube::operaMultiExtensionFITSCube(string Filename, edatatype Datatype, int mode/*=READWRITE|READONLY*/, unsigned Compression, bool IsLazy) : 
operaMultiExtensionFITSImage(Filename, Datatype, mode, Compression, IsLazy)
{
    isLazy = IsLazy;
	npixels_per_slice = naxis1*naxis2;
	npixels_per_extension = naxis1*naxis2*naxis3;
	imageType = MEF;
    if (isLazy) {   // if lazy read...
        npixels = npixels_per_extension;
    } else { 
        npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
    }
}
/* 
 * ooperaMultiExtensionFITSCube(operaMultiExtensionFITSCube &imageIn, bool ViewOnly, bool AddHeader)
 * \brief Clone a Multi Extension FITSImage object.
 * \param imageIn - pointer to image to clone
 * \return operaMultiExtensionFITSImage*
 */
operaMultiExtensionFITSCube::operaMultiExtensionFITSCube(operaMultiExtensionFITSCube &imageIn, bool ViewOnly, bool AddHeader) {
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
	imageType = MEFCube;
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
    // add the header values
    if (AddHeader){
        //
    } else {
        hdu = 0;
    }
}

operaMultiExtensionFITSCube::~operaMultiExtensionFITSCube() {
}

/* 
 * operaMultiExtensionFITSCubeCopyHeader(operaMultiExtensionFITSCube *from, unsigned extension) 
 * \brief Copies all of the header information from image. 
 * If AllExtensions copies all headers, including primary.
 * \param from
 * \param extension
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSCube::operaMultiExtensionFITSCubeCopyHeader(operaMultiExtensionFITSCube *from) {
	int status = 0;
	const int morekeys = 100;
	int hdutype = ANY_HDU;
	
	if (!fptr) {
		throw operaException("operaMultiExtensionFITSCube: ", operaErrorCodeNULLfptr, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (!from->getfitsfileptr()) {
		throw operaException("operaMultiExtensionFITSCube: ", operaErrorCodeNULLfptr, __FILE__, __FUNCTION__, __LINE__);	
	}
	// if hdu == NULL, image was only created in memory
	if (hdu == 0) {
		// Create header (primary)
		if (fits_create_img(fptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		//	get current hdunum
		fits_get_hdu_num(fptr, (int *)&hdu);
	}
	if ( fits_movabs_hdu(from->getfitsfileptr(), cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
#if 0
	/* copy all the user keywords (not the structural keywords) */
	int nkeys = 0;
	char card[81];
	if (fits_get_hdrspace(from->getfitsfileptr(), &nkeys, NULL, &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	} 
	for (int ii = 1; ii <= nkeys; ii++) {
		if (fits_read_record(from->getfitsfileptr(), ii, card, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
			if (fits_write_record(fptr, card, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
#endif
	fits_delete_key(fptr, "NAXIS1", &status);
	status = 0;
	fits_delete_key(fptr, "NAXIS2", &status);
	status = 0;
	
	for (unsigned ext=0; ext<=extensions; ext++) {
		if ( fits_delete_hdu(fptr, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
	for (unsigned ext=0; ext<=extensions; ext++) {
		if ( fits_movabs_hdu(from->getfitsfileptr(), cfitsioextension(0), &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		if ( fits_movrel_hdu(from->getfitsfileptr(), ext, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		// copy the current HDU from file associated with from and appends to the
		// end of the file associated with fptr (in this case the end is the primary)
		if (fits_copy_hdu(from->getfitsfileptr(), fptr, morekeys, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(fptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (datatype == tfloat || datatype == tdouble) {
			float bzero = 0.0, bscale = 1.0;
			//if (fits_set_bscale(fptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
			//	throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			//}
			if (fits_update_key_flt(fptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(fptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
	}
	if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if ( fits_movrel_hdu(fptr, 1, &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if ( fits_delete_hdu(fptr, &hdutype, &status) ) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
	}
}

/* 
 * operaMultiExtensionFITSCubeSave() 
 * \brief Saves the current image to disk.
 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
 * \throws operaException operaErrorCodeNoFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSCube::operaMultiExtensionFITSCubeSave() {
	int status = 0;
	unsigned x1 = 1;
	unsigned x2 = x1+naxis1;
	unsigned y1 = 1;
	unsigned y2 = y1+naxis2;
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int hdutype = ANY_HDU;
	fitsfile *newfptr;         // FITS file pointer for the new image
	const long fpixel = 1;
	const int morekeys = 100;
	
    // check if there is a filename 
	if (filename.empty())
		throw operaException("operaMultiExtensionFITSCube: ", operaErrorCodeNoFilename, __FILE__, __FUNCTION__, __LINE__);	
    // check if there is a filename 
	if (mode == READONLY)
		throw operaException("operaMultiExtensionFITSCube: "+filename+" ", operaErrorCodeChangeREADONLYError, __FILE__, __FUNCTION__, __LINE__);	
	// remove existing file - cfitsio returns an error if it exists...
	remove(filename.c_str());
    // creates a fits file with filename filename.c_str()
	if (fits_create_file(&newfptr, filename.c_str(), &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
    // if fptr == NULL, image was only created in memory
    if (fptr == NULL) {
        // Create headers (including primary) and images
		if (fits_create_img(newfptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
        // copy the current HDU from file associated with fptr and appends to the
        // end of the file associated with newfptr (in this case the end is the primary)
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}  
    }
	if (fits_update_key_lng(newfptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key_lng(newfptr, "NAXIS", (unsigned short)0, (char *)"Number of axes", &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	fits_delete_key(newfptr, "NAXIS1", &status);
	status = 0;
	fits_delete_key(newfptr, "NAXIS2", &status); 
	status = 0;
    // set the compression type for the file associated with newfptr
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	void *basepointer = pixptr;
	// save each extension with header
	for (unsigned extension=1; extension<=extensions; extension++) {
		// move to the specified HDU number  (select extension)
		if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if ( fits_movrel_hdu(fptr, extension, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
		// copy the current HDU from the file associated with fptr to the file associated with newfptr
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		fits_read_keyword(fptr, "DETSIZE", zvalue, zcomment, &status);
		if (status) {	// keyword does not exist, add the MEF base keywords
			status = 0;
			char buff[1024];
			sprintf(buff, "[%u:%u,%u:%ld]", x1, x2, y1, naxis2*extensions+1);
			// DETSIZE and DETSEC keywords both REQUIRED for file to open as a MEF.
			if (fits_update_key(newfptr, tstring, "DETSIZE", buff, (char *)"Total data pixels in full mosaic", &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			sprintf(buff, "[%u:%u,%u:%u]", x1, x2, y1, y2);
			// DETSEC changes for each extension.
			if (fits_update_key(newfptr, tstring, "DETSEC", buff, (char *)"Mosaic area of the detector", &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
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
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
#endif
		if (AllExtensions) {
			// write elements into the fits file
			if (fits_write_img(newfptr, datatype, fpixel, npixels_per_extension, basepointer, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
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
					throw operaException("operaMultiExtensionFITSCube: cfitsio error "+filename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
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
 * operaMultiExtensionFITSCubeSaveAs(string newFilename) 
 * \brief Saves the current image to disk, with the given filename.
 * \param newFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaMultiExtensionFITSCube::operaMultiExtensionFITSCubeSaveAs(string newFilename) {
	int status = 0;
	unsigned x1 = 1;
	unsigned x2 = x1+naxis1;
	unsigned y1 = 1;
	unsigned y2 = y1+naxis2;
	char zvalue[FLEN_VALUE],zcomment[FLEN_COMMENT];
	int hdutype = ANY_HDU;
	fitsfile *newfptr;         // FITS file pointer for the new image
	const long fpixel = 1;
	const int morekeys = 100;
	
	// remove existing file - cfitsio returns an error if it exists...
	remove(newFilename.c_str());
	// creates a fits file with filename newFilename.c_str()
	if (fits_create_file(&newfptr, newFilename.c_str(), &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio create error: "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	// if fptr == NULL, image was only created in memory
	if (fptr == NULL) {
		// Create headers (including primary) and images
		if (fits_create_img(newfptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
		// copy the current HDU from file associated with fptr and appends to the
		// end of the file associated with newfptr (in this case the end is the primary)
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
#if 0
	// be sure we have the basic headers for floats if no header existed
	if (datatype == tfloat || datatype == tdouble) {
		float bzero = 0.0, bscale = 1.0;
		if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
#endif
	if (fits_update_key_lng(newfptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key_lng(newfptr, "NAXIS", (unsigned short)0, (char *)"Number of axes", &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_delete_key(newfptr, "NAXIS1", &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_delete_key(newfptr, "NAXIS2", &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	// set the compression type for the file associated with newfptr
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}    
	void *basepointer = pixptr;
	// save each extension with header
	for (unsigned extension=1; extension<=extensions; extension++) {
		// move to the specified HDU number  (select extension)
		if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if ( fits_movrel_hdu(fptr, extension, &hdutype, &status) ) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);	
		}
		// copy the current HDU from the file associated with fptr to the file associated with newfptr
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
#if 0
		// be sure we have the basic headers for floats if no header existed
		if (datatype == tfloat || datatype == tdouble) {
			float bzero = 0.0, bscale = 1.0;
			if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
#endif
		fits_read_keyword(fptr, "DETSIZE", zvalue, zcomment, &status);
		if (status) {	// keyword does not exist, add the MEF base keywords
			status = 0;
			char buff[1024];
			sprintf(buff, "[%u:%u,%u:%ld]", x1, x2, y1, naxis2*extensions+1);
			// DETSIZE and DETSEC keywords both REQUIRED for file to open as a MEF.
			if (fits_update_key(newfptr, tstring, "DETSIZE", buff, (char *)"Total data pixels in full mosaic", &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			sprintf(buff, "[%u:%u,%u:%u]", x1, x2, y1, y2);
			// DETSEC changes for each extension.
			if (fits_update_key(newfptr, tstring, "DETSEC", buff, (char *)"Mosaic area of the detector", &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			y1 = y2+1;
			y2 += naxis2;
			
		}
		status = 0;
		if (AllExtensions) {
			// write elements into the fits file
			if (fits_write_img(newfptr, datatype, fpixel, npixels_per_extension, basepointer, &status)) {
				throw operaException("operaMultiExtensionFITSCube: cfitsio error "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
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
					throw operaException("operaMultiExtensionFITSCube: cfitsio error "+newFilename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
					break;
			}
		}
	}
	fits_close_file(newfptr, &status);
	newfptr = NULL;
	// note that fptr is left open in case the user wants to save again, it will be closed
	// by a close call if made, or the destructor otherwise
}

typedef struct thread_args {
	operaMultiExtensionFITSImage *output;
	unsigned extension;
	unsigned maxx;
	unsigned maxy;
	unsigned maxz;
} thread_args_t;

static pthread_t *threads = NULL;

static void collapse(operaMultiExtensionFITSImage *output, const unsigned extension, const unsigned maxx, const unsigned maxy, const unsigned maxz) {
	unsigned long npixels_per_extension = output->getnaxis1() * output->getnaxis2() * output->getnaxis3();
	unsigned long npixels_per_slice = output->getnaxis1() * output->getnaxis2();
	unsigned int npixels_per_row = output->getnaxis1();
	float *base = (float *)output->getpixels() + (extension-1)*npixels_per_extension;
	float *baseout = (float *)output->getpixels() + (extension-1)*npixels_per_extension;
	float *pixelStack = new float[maxz];
	for (unsigned int y=0; y<maxy; y++) {
		for (unsigned int x=0; x<maxx; x++) {
			for (unsigned int z=0; z<maxz; z++) {
				pixelStack[z] = *(base + z*npixels_per_slice + (y*npixels_per_row) + x);
			}
			*(baseout + (y*npixels_per_row) + x) = operaArrayMedianQuick(maxz, (float *)pixelStack);
		}
	}
	delete[] pixelStack;
}

static void *collapsethread(void *argument) {
	thread_args_t *thread_args_s = (thread_args_t *)argument;
	collapse(thread_args_s->output, thread_args_s->extension, thread_args_s->maxx, thread_args_s->maxy, thread_args_s->maxz);
	return NULL;
}

static void processSingleExtension(const unsigned extension, thread_args_t *thread_args, operaMultiExtensionFITSImage *output, const unsigned maxx, const unsigned maxy, const unsigned maxz) {
	thread_args[0].extension = extension;
	thread_args[0].output = output;
	thread_args[0].maxx = maxx;
	thread_args[0].maxy = maxy;
	thread_args[0].maxz = maxz;
    collapsethread((void *) &thread_args[0]);
}

static bool spawnthreads(unsigned int extension, unsigned int maxextension, unsigned int count, thread_args_t *thread_args, operaMultiExtensionFITSImage *output, const unsigned maxx, const unsigned maxy, const unsigned maxz) {
	unsigned int j = 0;
    for (unsigned int i=extension; i<=maxextension; i++) {
		thread_args[i].extension = i;
		thread_args[i].output = output;
		thread_args[i].maxx = maxx;
		thread_args[i].maxy = maxy;
		thread_args[i].maxz = maxz;
		if (pthread_create(&threads[i], NULL, collapsethread, (void *) &thread_args[i]) != 0)
			return false;
        if (++j >= count)
            break;
	}
    return true;
}

static bool waitthreads(unsigned int extension, unsigned int maxextension, unsigned int count) {
	unsigned int j = 0;
	for (unsigned int i=extension; i<=maxextension; i++) {
		if (pthread_join(threads[i], NULL) != 0)
			return false;
        if (++j >= count)
            break;
	}
    return true;
}

static bool processExtensions(unsigned int minextension, unsigned int maxextension, unsigned int maxthreads, thread_args_t *thread_args, operaMultiExtensionFITSImage *output, const unsigned maxx, const unsigned maxy, const unsigned maxz) {
	bool worked = true;
	for (unsigned int extension=minextension; extension<=maxextension; extension+=maxthreads) {
        worked &= spawnthreads(extension, maxextension, maxthreads, thread_args, output, maxx, maxy, maxz);
        worked &= waitthreads(extension, maxextension, maxthreads);
	}
	return worked;
}
/*
 * void medianCollapse()
 * \brief collapse the cube into a single image
 * Collapses entire image if AllExtensions = true
 * Collapses current_extension if AllExtensions = false
 * To collapse a single extension, use medianCollapseSingleExtension(extension)
 * Note: indexing for a MEF Cube is [extension][slice][y][x]
 */
void operaMultiExtensionFITSCube::medianCollapse() {
    operaMultiExtensionFITSImage &output = *this;
    const unsigned maxx = output.getXDimension();
    const unsigned maxy = output.getYDimension();
    const unsigned maxz = output.getZDimension();
    const unsigned extensions = output.getNExtensions();
    const int current_extension = output.getExtension();
    
    if (maxz == 1) {
        return;                   // no collapse to be done because there is only one slice
    } else {                      // if more than one slice, need to collapse slices
        if (AllExtensions) {      // all extensions read in
            for (unsigned extension=1; extension<=extensions; extension++) {
				collapse(this, extension, maxx, maxy, maxz);
				// Now, we have to move the extensions down to make them contiguous
				if (extension > 1) {
					float *newextensionptr = (float *)pixptr+naxis1*naxis2*(extension-1);
					float *oldextensionptr = (float *)pixptr+naxis1*naxis2*naxis3*(extension-1);
					memcpy(newextensionptr, oldextensionptr, sizeof(float *)*naxis1*naxis2);
				}
            }
        } else {                  // collapses only current_extension
			collapse(this, current_extension, maxx, maxy, maxz);
            // if only reading in one extension at a time, save after changes
			saveExtension(current_extension, first_slice);
        }
		naxes[2] = naxis3 = 1;
		current_slice = 1;
		npixels_per_extension = naxis1*naxis2*naxis3;
		npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
    }
}
/*
 * void medianCollapse(unsigned int n)
 * \brief collapse n slices of the cube into a single image
 * \details Collapses entire image if AllExtensions = true
 * Collapses current_extension if AllExtensions = false
 * To collapse a single extension, use medianCollapseSingleExtension(extension)
 * WARNING: that this method uses operaArrayMedianQuick, so the input image will be scrambled.
 */
void operaMultiExtensionFITSCube::medianCollapse(unsigned int n) {
	operaMultiExtensionFITSCube &output = *this;
	const unsigned maxx = output.getXDimension();
	const unsigned maxy = output.getYDimension();
	const unsigned maxz = n;
	const unsigned extensions = output.getNExtensions();
	const int current_extension = output.getExtension();

	if (maxz == 1) {
		return;                   // no collapse to be done because there is only one slice
	} else {                      // if more than one slice, need to collapse slices
		if (AllExtensions) {      // all extensions read in
			for (unsigned extension=1; extension<=extensions; extension++) {
				collapse(this, extension, maxx, maxy, n);
				// Now, we have to move the extensions down to make them contiguous
				if (extension > 1) {
					float *newextensionptr = (float *)pixptr+naxis1*naxis2*(extension-1);
					float *oldextensionptr = (float *)pixptr+naxis1*naxis2*naxis3*(extension-1);
					memcpy(newextensionptr, oldextensionptr, sizeof(float *)*naxis1*naxis2);
				}
			}
		} else {                  // collapses only current_extension
			collapse(this, current_extension, maxx, maxy, n);
			// if only reading in one extension at a time, save after changes
			saveExtension(current_extension, first_slice);
		}
		naxes[2] = naxis3 = 1;
		current_slice = 1;
		npixels_per_extension = naxis1*naxis2*naxis3;
		npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
	}
}
/*
 * void medianCollapseParallel(unsigned int maxthreads)
 * \brief extension-based parallel collapse the cube into a single image
 * Collapses entire image if AllExtensions = true
 * Collapses current_extension if AllExtensions = false
 * To collapse a single extension, use medianCollapseSingleExtension(extension)
 * Note: indexing for a MEF Cube is [extension][slice][y][x]
 */
void operaMultiExtensionFITSCube::medianCollapseParallel(unsigned int maxthreads) {
    operaMultiExtensionFITSCube &output = *this;
    const unsigned maxx = output.getXDimension();
    const unsigned maxy = output.getYDimension();
    const unsigned maxz = output.getZDimension();
    const unsigned extensions = output.getNExtensions();
    const int current_extension = output.getExtension();
    
	unsigned long nthreads = maxthreads+1;
	threads = (pthread_t *)calloc(nthreads, sizeof(pthread_t*));
	thread_args_t *thread_args = (thread_args_t *)calloc(nthreads, sizeof(thread_args_t));
	
    if (maxz == 1) {
        return;                   // no collapse to be done because there is only one slice
    } else {                      // if more than one slice, need to collapse slices
        if (AllExtensions) {      // all extensions read in
			if (maxthreads > 1) {
				processExtensions(1, extensions, maxthreads, thread_args, this, maxx, maxy, maxz);
			} else {
				for (unsigned int extension=1; extension<=extensions; extension++) {
					processSingleExtension(extension, thread_args, this, maxx, maxy, maxz);
				}
			}
			// Now, we have to move the extensions down to make them contiguous
			for (unsigned int extension=1; extension<=extensions; extension++) {
				if (extension > 1) {
					float *newextensionptr = (float *)pixptr+naxis1*naxis2*(extension-1);
					float *oldextensionptr = (float *)pixptr+naxis1*naxis2*naxis3*(extension-1);
					memcpy(newextensionptr, oldextensionptr, sizeof(float *)*naxis1*naxis2);
				}
			}
        } else {                  // collapses only current_extension
			collapse(this, current_extension, maxx, maxy, maxz);
            // if only reading in one extension at a time, save after changes
			saveExtension(current_extension, first_slice);
        }
		naxes[2] = naxis3 = 1;
		current_slice = 1;
		npixels_per_extension = naxis1*naxis2*naxis3;
		npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
    }
}
/*
 * void medianCollapseParallel(unsigned int maxthreads, unsigned int n)
 * \brief extension-based parallel collapse n slices of the cube into a single image
 * Collapses entire image if AllExtensions = true
 * Collapses current_extension if AllExtensions = false
 * To collapse a single extension, use medianCollapseSingleExtension(extension)
 * Note: indexing for a MEF Cube is [extension][slice][y][x]
 */
void operaMultiExtensionFITSCube::medianCollapseParallel(unsigned int maxthreads, unsigned int n) {
    operaMultiExtensionFITSImage &output = *this;
    const unsigned maxx = output.getXDimension();
    const unsigned maxy = output.getYDimension();
    const unsigned maxz = n;
    const unsigned extensions = output.getNExtensions();
    const int current_extension = output.getExtension();
    
	unsigned long nthreads = maxthreads+1;
	threads = (pthread_t *)calloc(nthreads, sizeof(pthread_t*));
	thread_args_t *thread_args = (thread_args_t *)calloc(nthreads, sizeof(thread_args_t));
    if (maxz == 1) {
        return;                   // no collapse to be done because there is only one slice
    } else {                      // if more than one slice, need to collapse slices
        if (AllExtensions) {      // all extensions read in
			if (maxthreads > 1) {
				processExtensions(1, extensions, maxthreads, thread_args, this, maxx, maxy, maxz);
			} else {
				for (unsigned int extension=1; extension<=extensions; extension++) {
					processSingleExtension(extension, thread_args, this, maxx, maxy, maxz);
				}
			}
			// Now, we have to move the extensions down to make them contiguous
			for (unsigned int extension=1; extension<=extensions; extension++) {
				if (extension > 1) {
					float *newextensionptr = (float *)pixptr+naxis1*naxis2*(extension-1);
					float *oldextensionptr = (float *)pixptr+naxis1*naxis2*naxis3*(extension-1);
					memcpy(newextensionptr, oldextensionptr, sizeof(float *)*naxis1*naxis2);
				}
			}
        } else {                  // collapses only current_extension
			collapse(this, current_extension, maxx, maxy, maxz);
            // if only reading in one extension at a time, save after changes
			saveExtension(current_extension, first_slice);
        }
		naxes[2] = naxis3 = 1;
		current_slice = 1;
		npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);
		npixels_per_extension = naxis1*naxis2*naxis3;
    }
}
/*
 * void medianCollapseSingleExtension(int extension)
 * \brief collapse the cube into a single image
 * Collapses a single extension
 * To collapse all extensions, use medianCollapse in non-lazy read (AllExtensions = true)
 * \param extension 
 * WARNING: that this method utilizes operaArrayMedianQuick, so the input image will be scrambled.
 */
void operaMultiExtensionFITSCube::medianCollapseSingleExtension(int extension) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > (int)getNExtensions()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    setExtension(extension);
    
    operaMultiExtensionFITSCube &output = *this;
    const unsigned maxx = output.getXDimension();
    const unsigned maxy = output.getYDimension();
    const unsigned maxz = output.getZDimension();
    
    if (maxz == 1) {
        return;                         // no collapse to be done because there is only one slice
    } else {                            // collapses only current_extension
		collapse(this, extension, maxx, maxy, maxz);
        if (!AllExtensions) {
            saveExtension(current_extension, first_slice);
        }
    }       
}
/* 
 * void getSlice()
 * \brief get the current slice
 * \param unsigned slice
 */
unsigned operaMultiExtensionFITSCube::getSlice() {
    return current_slice;
}
/* 
 * unsigned getslices(void)
 * \brief gets the number of slices
 */
unsigned operaMultiExtensionFITSCube::getslices(void) {
	return naxis3;
}
/*  
 * void setSlice(int slice)
 * \brief sets the current slice
 * \param int slice
 */
void operaMultiExtensionFITSCube::setSlice(unsigned slice) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
    if (slice > getslices()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	current_slice = slice;
}
/* 
 * void readSlice(int slice)
 * \brief read a slice into memory
 * \param int slice
 * \Note The slice is read in as the first slice
 */
void operaMultiExtensionFITSCube::readSlice(unsigned slice) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
    if (slice > getslices()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
	int status = 0;
	long fpixel[MAXFITSDIMENSIONS], lpixel[MAXFITSDIMENSIONS], inc[MAXFITSDIMENSIONS];
	inc[0] = 1; 
	inc[1] = 1;
	inc[2] = 1;
	fpixel[0] = 1; 
	fpixel[1] = 1;
	fpixel[2] = slice;
	lpixel[0] = naxes[0]; 
	lpixel[1] = naxes[1];
	lpixel[2] = slice;
	
	/* load image slice into the current extension */
	if (isLazy && pixptr == NULL) {
		size_t size = toSize(bitpix, npixels);
        pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
		memset(pixptr, 0, size);		
		 
	}
	if (fits_read_subset(fptr, datatype, fpixel, lpixel, inc, NULL, pixptr, NULL, &status))
		throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);			
}
/* 
 * void saveSlice(int slice)
 * \brief save image data
 * \param int slice
 */
void operaMultiExtensionFITSCube::saveSlice(unsigned slice) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
    if (slice > getslices()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    int status = 0;
    long fpixel[MAXFITSDIMENSIONS], lpixel[MAXFITSDIMENSIONS];
    fpixel[0] = 1; 
    fpixel[1] = 1;
    fpixel[2] = slice;
    lpixel[0] = naxes[0]; 
    lpixel[1] = naxes[1];
    lpixel[2] = slice;
    
    if (fits_write_subset(fptr, datatype, fpixel, lpixel, pixptr, &status))
        throw operaException("operaMultiExtensionFITSCube: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
}
/* 
 * void setExtension(int extension)
 * \brief sets the new_extension to the current extension
 * \param unsigned new_extension
 */
void operaMultiExtensionFITSCube::setExtension(unsigned new_extension) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
    if (new_extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    operaMultiExtensionFITSImage &output = *this;
    output.AllSlices = true;
	output.setExtension(new_extension);
}
/* 
 * void setExtension(int extension, int slice)
 * \brief sets the newe_xtension to the current extension
 * \param unsigned new_extension
 * \param int slice
 */
void operaMultiExtensionFITSCube::setExtension(unsigned new_extension, unsigned slice) {
#ifdef FITS_RANGE_CHECK  // Check that the given values are within the actual range 
    if (new_extension > getNExtensions() || slice > getslices()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    operaMultiExtensionFITSImage &output = *this;
	output.setExtension(new_extension);
	setSlice(slice);
}
/* 
 * void readExtension(int extension)
 * \brief read an extension into memory
 * \param operaMultiExtensionFITSCube &Image
 * \param int extension
 */
void operaMultiExtensionFITSCube::readExtension(unsigned extension) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    operaMultiExtensionFITSImage &output = *this;
    output.AllSlices = true;
    output.readExtension(extension);
}
/* 
 * void readExtension(int extension, int slice)
 * \brief read an extension into memory
 * \param operaMultiExtensionFITSCube &Image
 * \param int extension
 * \param int slice
 */
void operaMultiExtensionFITSCube::readExtension(unsigned extension, unsigned slice) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions() || slice > getslices()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    operaMultiExtensionFITSImage &output = *this;
    output.readExtension(extension);
	readSlice(slice);
}
/* 
 * void saveExtension(int extension) 
 * \brief save image data
 * \param int extension
 * \param int slice
 */
void operaMultiExtensionFITSCube::saveExtension(unsigned extension) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
    if (extension > getNExtensions()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    operaMultiExtensionFITSImage &output = *this;
    output.AllSlices = true;
    output.setExtension(extension);
}
/* 
 * void saveExtension(int extension, int slice)
 * \brief save image data
 * \param int extension
 * \param int slice
 */
void operaMultiExtensionFITSCube::saveExtension(unsigned extension, unsigned slice) {
#ifdef FITS_RANGE_CHECK  // Check that the given values are within the actual range 
    if (extension > getNExtensions() || slice > getslices()) {
        throw operaException("operaMultiExtensionFITSCube: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
    }
#endif
    operaMultiExtensionFITSImage &output = *this;
    output.setExtension(extension);
	saveSlice(slice);
}
/*
 * bool isMEFCube()
 * \brief Is this image a MEFcube?
 * \return bool
 */    
bool operaMultiExtensionFITSCube::isMEFCube(void){
    return extensions > 1 && getZDimension() > 1;
}
