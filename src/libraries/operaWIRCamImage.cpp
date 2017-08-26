/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaWIRCamImage
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
#include "libraries/operaException.h"
#include "libraries/operaImageVector.h"				// for operaImageVector
#include "libraries/operaImage.h"
#include "libraries/operaLib.h"						// for itos
#include "libraries/operaWIRCamImage.h"

#include "libraries/ladfit.h"

/*! \file operaWIRCamImage.cpp */

using namespace std;

/*! 
 * operaWIRCamImage
 * \author Megan Tannock
 * \brief This class encapsulates the WIRCam FITS image.
 * \ingroup libraries
 */

operaWIRCamImage::operaWIRCamImage() : operaMultiExtensionFITSCube(),
ExposureTime(WIRCAM_MINIMUM_EXPOSURE_TIME),
ChipBias(WIRCAM_CHIPBIAS),
haveHeaderSkyLevels(false)
{
}
/*
 * \class operaWIRCamImage
 * \brief construct an in-memory FITSImage object
 * \brief operaWIRCamImage(unsigned Naxis1, unsigned Naxis2, edatatype Datatype=tushort)
 * \brief Create an in-memory image of given dimensions.
 * \param Naxis1 - x ccd dimension
 * \param Naxis2 - y ccd dimension
 * \param Datatype optional datatype defaults to tshort
 * \param Compression optional compression, defaults to none
 * \return void
 */
operaWIRCamImage::operaWIRCamImage(unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, unsigned Extensions, edatatype Datatype) : 
operaMultiExtensionFITSCube(Naxis1, Naxis2, Naxis3, Extensions, Datatype),
ExposureTime(WIRCAM_MINIMUM_EXPOSURE_TIME),
ChipBias(WIRCAM_CHIPBIAS),
haveHeaderSkyLevels(false)
{
}
/* 
 * \class operaWIRCamImage
 * \brief create a writeable file image of an in memory FITSImage object
 * \brief operaWIRCamImage(string Filename, unsigned Naxis1, unsigned Naxis2, edatatype Datatype, unsigned Compression)
 * \brief Constructor to create a new FITS file.
 * \param Filename to create (file is deleted if it exists)
 * \param Naxis1 - dimensions
 * \param Naxis2 - dimensions
 * \param Datatype defaults to tshort
 * \param Compression, defaults to none
 * \throws operaException cfitsio error code
 * \return void
 */
operaWIRCamImage::operaWIRCamImage(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, unsigned Extensions, edatatype Datatype, unsigned Compression, bool IsLazy) : 
operaMultiExtensionFITSCube(Filename, Naxis1, Naxis2, Naxis3, Extensions, Datatype, Compression, IsLazy),
ExposureTime(WIRCAM_MINIMUM_EXPOSURE_TIME),
ChipBias(WIRCAM_CHIPBIAS),
haveHeaderSkyLevels(false)
{
    // Note: Use constructors below to fill in guide cubes.
    try {
        for (unsigned extension=1; extension<=getNExtensions(); extension++) {
            for (unsigned slice=1; slice<=getZDimension(); slice++) {
                char skyheader[12];
                sprintf(skyheader, "SKYLVL%02d", slice);
                skyBackground[extension-1][slice-1].skyLevel = operaFITSGetFloatHeaderValue(skyheader, extension);
                sprintf(skyheader, "SKYDEV%02d", slice);
                skyBackground[extension-1][slice-1].skyDeviation = operaFITSGetFloatHeaderValue(skyheader, extension);
                skyBackground[extension-1][slice-1].skyRate = skyBackground[extension-1][slice-1].skyLevel / ExposureTime;
            }
        }
        haveHeaderSkyLevels = true;
    }
    catch (...) {
        // it is OK, just note that the levels could not be found
        haveHeaderSkyLevels = false;
    }
    
}
/* 
 * \class operaWIRCamImage
 * \brief create a FITSIMage object from a FITS file
 * \brief operaWIRCamImage(string Filename, int mode=READWRITE)
 * \brief Constructor to create a FITSImage from a FITS file.
 * \param Filename
 * \param mode
 * \return void
 */
operaWIRCamImage::operaWIRCamImage(string Filename, edatatype Datatype, int mode, unsigned Compression, bool IsLazy) : 
operaMultiExtensionFITSCube(Filename, Datatype, mode, Compression, IsLazy), 
ExposureTime(WIRCAM_MINIMUM_EXPOSURE_TIME),
ChipBias(WIRCAM_CHIPBIAS),
haveHeaderSkyLevels(false)
{
	// fill in the guide cubes....
    operaWIRCamImage &output = *this;
	output.setExtension(1);
	int dx = operaFITSGetIntHeaderValue("WCSIZEX");
	int dy = operaFITSGetIntHeaderValue("WCSIZEY");
	GuideWindows[0].setX1(operaFITSGetIntHeaderValue("WCPOSX1"));
	GuideWindows[0].setY1(operaFITSGetIntHeaderValue("WCPOSY1"));
	GuideWindows[0].setX2(GuideWindows[0].getX1()+dx);
	GuideWindows[0].setY2(GuideWindows[0].getY1()+dy);
	
	GuideWindows[1].setX1(operaFITSGetIntHeaderValue("WCPOSX2"));
	GuideWindows[1].setX1(operaFITSGetIntHeaderValue("WCPOSY2"));
	GuideWindows[1].setX2(GuideWindows[1].getX1()+dx);
	GuideWindows[1].setY2(GuideWindows[1].getY1()+dy);
	
	GuideWindows[2].setX1(operaFITSGetIntHeaderValue("WCPOSX3"));
	GuideWindows[2].setX1(operaFITSGetIntHeaderValue("WCPOSY3"));
	GuideWindows[2].setX2(GuideWindows[2].getX1()+dx);
	GuideWindows[2].setY2(GuideWindows[2].getY1()+dy);
	
	GuideWindows[3].setX1(operaFITSGetIntHeaderValue("WCPOSX4"));
	GuideWindows[3].setX1(operaFITSGetIntHeaderValue("WCPOSY4"));
	GuideWindows[3].setX2(GuideWindows[3].getX1()+dx);
	GuideWindows[3].setY2(GuideWindows[3].getY1()+dy);
	
	ChipBias = operaFITSGetFloatHeaderValue("CHIPBIAS");
	ExposureTime = operaFITSGetFloatHeaderValue("EXPTIME");
	
	try {
		for (unsigned extension=1; extension<=getNExtensions(); extension++) {
			for (unsigned slice=1; slice<=getZDimension(); slice++) {
				char skyheader[12];
				sprintf(skyheader, "SKYLVL%02d", slice);
				skyBackground[extension-1][slice-1].skyLevel = operaFITSGetFloatHeaderValue(skyheader, extension);
				sprintf(skyheader, "SKYDEV%02d", slice);
				skyBackground[extension-1][slice-1].skyDeviation = operaFITSGetFloatHeaderValue(skyheader, extension);
				skyBackground[extension-1][slice-1].skyRate = skyBackground[extension-1][slice-1].skyLevel / ExposureTime;
			}
		}
		haveHeaderSkyLevels = true;
	}
	catch (...) {
		// it is OK, just note that the levels could not be found
		haveHeaderSkyLevels = false;
	}
}
/* 
 * \class operaWIRCamImage
 * \brief create a operaWIRCamImage object from another operaWIRCamImage
 * \brief operaWIRCamImage(operaWIRCamImage &image)
 * \param operaWIRCamImage
 */
operaWIRCamImage::operaWIRCamImage(operaWIRCamImage &imageIn, bool ViewOnly, bool AddHeader)  : 
operaMultiExtensionFITSCube(),
ExposureTime(WIRCAM_MINIMUM_EXPOSURE_TIME),
ChipBias(WIRCAM_CHIPBIAS),
haveHeaderSkyLevels(false)
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
	imageType = MEFCube;
	bitpix = imageIn.bitpix;
	isLazy = imageIn.isLazy;
	AllExtensions = imageIn.AllExtensions;
	AllSlices = imageIn.AllSlices;
	
	ChipBias = imageIn.ChipBias;
	ExposureTime = imageIn.ExposureTime;
	haveHeaderSkyLevels = imageIn.haveHeaderSkyLevels;
	
	for (unsigned extension=1; extension<=getNExtensions(); extension++) {
		for (unsigned slice=1; slice<=getZDimension(); slice++) {
			skyBackground[extension-1][slice-1].skyLevel = imageIn.skyBackground[extension-1][slice-1].skyLevel;
			skyBackground[extension-1][slice-1].skyDeviation = imageIn.skyBackground[extension-1][slice-1].skyDeviation;
			skyBackground[extension-1][slice-1].skyRate = imageIn.skyBackground[extension-1][slice-1].skyRate;
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
}

/* 
 * \class operaWIRCamImage
 * \brief destructor
 */
operaWIRCamImage::~operaWIRCamImage() {
}

/* 
 * ~operaWIRCamImageSave() 
 * \brief Saves the current image to disk.
 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
 * \throws operaException operaErrorCodeNoFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaWIRCamImage::operaWIRCamImageSave() {
	int status = 0;
	int hdutype = ANY_HDU;
	fitsfile *newfptr;         // FITS file pointer for the new image
	const long fpixel = 1;
	const int morekeys = 100;
	
    // check if there is a filename 
	if (filename.empty())
		throw operaException("operaWIRCamImage: ", operaErrorCodeNoFilename, __FILE__, __FUNCTION__, __LINE__);	
    // check if there is a filename 
	if (mode == READONLY)
		throw operaException("operaWIRCamImage: "+filename+" ", operaErrorCodeChangeREADONLYError, __FILE__, __FUNCTION__, __LINE__);	
	// remove existing file - cfitsio returns an error if it exists...
	remove(filename.c_str());
    // creates a fits file with filename filename.c_str()
	if (fits_create_file(&newfptr, filename.c_str(), &status)) {
		throw operaException("operaWIRCamImage: cfitsio create error: "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
    // if fptr == NULL, image was only created in memory
    if (fptr == NULL) {
        // Create headers (including primary) and images
		if (fits_create_img(newfptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
        // copy the current HDU from file associated with fptr and appends to the
        // end of the file associated with newfptr (in this case the end is the primary)
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        // DEFAULT HEADER KEYWORDS FOR WIRCAM
        if (fits_write_key(newfptr, tstring, "DETECTOR", (char *)WIRCAM_DETECTOR, NULL, &status)) {
            throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
		fits_delete_key(newfptr, "DETSIZE", &status);
        if (fits_write_key(newfptr, tstring, "DETSIZE", (char *)WIRCAM_DETSIZE, (char *)"Total data pixels in full mosaic", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CCDBIN1", WIRCAM_CCDBIN1, -1, (char *)"Binning factor along first axis", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CCDBIN2", WIRCAM_CCDBIN2, -1, (char *)"Binning factor along second axis", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSIZE1", WIRCAM_PIXSIZE1, -1, (char *)"Pixel size for axis 1 (microns)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSIZE2", WIRCAM_PIXSIZE2, -1, (char *)"Pixel size for axis 2 (microns)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSCAL1", WIRCAM_PIXSCAL1, -1, (char *)"Pixel scale for axis 1 (arcsec/pixel)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSCAL2", WIRCAM_PIXSCAL2, -1, (char *)"Pixel scale for axis 2 (arcsec/pixel)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "IRFPA1", (char *)WIRCAM_EXTNAME1, (char *)"Name of the detector in ext 1", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA2", (char *)WIRCAM_EXTNAME2, (char *)"Name of the detector in ext 2", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA3", (char *)WIRCAM_EXTNAME3, (char *)"Name of the detector in ext 3", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA4", (char *)WIRCAM_EXTNAME4, (char *)"Name of the detector in ext 4", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_update_key_flt(newfptr, "CHIPBIAS", ChipBias, -1, (char *)"Science frame chip bias added to CDS value", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "DETSECA", (char *)WIRCAM_DETSECA, (char *)"Guider Mosaic Area of IRFPA1", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_write_key(newfptr, tstring, "DETSECB", (char *)WIRCAM_DETSECB, (char *)"Guider Mosaic Area of IRFPA2", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
        if (fits_write_key(newfptr, tstring, "DETSECC", (char *)WIRCAM_DETSECC, (char *)"Guider Mosaic Area of IRFPA3", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
        if (fits_write_key(newfptr, tstring, "DETSECD", (char*)WIRCAM_DETSECD, (char *)"Guider Mosaic Area of IRFPA4", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        // END OF DEFAULT HEADER KEYWORDS FOR WIRCAM
    }
	if (fits_update_key_lng(newfptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
		throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key_lng(newfptr, "NAXIS", (unsigned short)0, (char *)"Number of axes", &status)) {
		throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	fits_delete_key(newfptr, "NAXIS1", &status);
	status = 0;
	fits_delete_key(newfptr, "NAXIS2", &status); 
	status = 0;
    // set the compression type for the file associated with newfptr
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	void *basepointer = pixptr;
	// save each extension with header
	for (unsigned extension=1; extension<=extensions; extension++) {
		// move to the specified HDU number  (select extension)
		if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
			throw operaException("operaWIRCamImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		if ( fits_movrel_hdu(fptr, extension, &hdutype, &status) ) {
			throw operaException("operaWIRCamImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		// copy the current HDU from the file associated with fptr to the file associated with newfptr
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
#if 0
		// be sure we have the basic headers for floats if no header existed
		if (datatype == tfloat || datatype == tdouble) {
			float bzero = 0.0, bscale = 1.0;
			if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
				throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
#endif
        // DEFAULT HEADER KEYWORDS FOR WIRCAM
        if (extension == 1){
            // EXTNAME and DETSEC change for each extension
            if (fits_write_key(newfptr, tstring, "EXTNAME", (char *)WIRCAM_EXTNAME1, (char *)"Extension name", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
            if (fits_update_key(newfptr, tstring, "DETSEC", (char *)WIRCAM_DETSEC1, (char *)"Mosaic area of the detector", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
        if (extension == 2){
            // EXTNAME and DETSEC change for each extension
            if (fits_write_key(newfptr, tstring, "EXTNAME", (char *)WIRCAM_EXTNAME2, (char *)"Extension name", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
            if (fits_update_key(newfptr, tstring, "DETSEC", (char *)WIRCAM_DETSEC2, (char *)"Mosaic area of the detector", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
        if (extension == 3){
            // EXTNAME and DETSEC change for each extension
            if (fits_write_key(newfptr, tstring, "EXTNAME", (char *)WIRCAM_EXTNAME3, (char *)"Extension name", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
            if (fits_update_key(newfptr, tstring, "DETSEC", (char *)WIRCAM_DETSEC3, (char *)"Mosaic area of the detector", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
        if (extension == 4){
            // EXTNAME and DETSEC change for each extension
            if (fits_write_key(newfptr, tstring, "EXTNAME", (char *)WIRCAM_EXTNAME4, (char *)"Extension name", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
            if (fits_update_key(newfptr, tstring, "DETSEC", (char *)WIRCAM_DETSEC4, (char *)"Mosaic area of the detector", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
        if (fits_write_key(newfptr, tstring, "DETECTOR", (char *)WIRCAM_DETECTOR, NULL, &status)) {
            throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
		fits_delete_key(newfptr, "DETSIZE", &status);
        if (fits_write_key(newfptr, tstring, "DETSIZE", (char *)WIRCAM_DETSIZE, (char *)"Total data pixels in full mosaic", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CCDBIN1", WIRCAM_CCDBIN1, -1, (char *)"Binning factor along first axis", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CCDBIN2", WIRCAM_CCDBIN2, -1, (char *)"Binning factor along second axis", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSIZE1", WIRCAM_PIXSIZE1, -1, (char *)"Pixel size for axis 1 (microns)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSIZE2", WIRCAM_PIXSIZE2, -1, (char *)"Pixel size for axis 2 (microns)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSCAL1", WIRCAM_PIXSCAL1, -1, (char *)"Pixel scale for axis 1 (arcsec/pixel)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSCAL2", WIRCAM_PIXSCAL2, -1, (char *)"Pixel scale for axis 2 (arcsec/pixel)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "IRFPA1", (char *)WIRCAM_EXTNAME1, (char *)"Name of the detector in ext 1", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA2", (char *)WIRCAM_EXTNAME2, (char *)"Name of the detector in ext 2", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA3", (char *)WIRCAM_EXTNAME3, (char *)"Name of the detector in ext 3", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA4", (char *)WIRCAM_EXTNAME4, (char *)"Name of the detector in ext 4", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }    
        if (fits_write_key(newfptr, tstring, "CHIPSIZE", (char *)WIRCAM_CHIPSIZE, (char *)"Detector imaging area sizes", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "DATASEC", (char *)WIRCAM_DATASEC, (char *)"Imaging area of the detector", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CHIPBIAS", ChipBias, -1, (char *)"Science frame chip bias added to CDS value", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "DETSECA", (char *)WIRCAM_DETSECA, (char *)"Guider Mosaic Area of IRFPA1", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_write_key(newfptr, tstring, "DETSECB", (char *)WIRCAM_DETSECB, (char *)"Guider Mosaic Area of IRFPA2", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
        if (fits_write_key(newfptr, tstring, "DETSECC", (char *)WIRCAM_DETSECC, (char *)"Guider Mosaic Area of IRFPA3", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
        if (fits_write_key(newfptr, tstring, "DETSECD", (char*)WIRCAM_DETSECD, (char *)"Guider Mosaic Area of IRFPA4", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}       
        // END DEFAULT HEADER KEYWORDS FOR WIRCAM         
		if (AllExtensions) {
			// write elements into the fits file
			if (fits_write_img(newfptr, datatype, fpixel, npixels_per_extension, basepointer, &status)) {
				throw operaException("operaWIRCamImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
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
					throw operaException("operaWIRCamImage: cfitsio error "+filename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
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
 * operaWIRCamImageSaveAs(string newFilename) 
 * \brief Saves the current image to disk, with the given filename.
 * \param newFilename
 * \throws operaException cfitsio error code
 * \return void
 */
void operaWIRCamImage::operaWIRCamImageSaveAs(string newFilename) {
	int status = 0;
	int hdutype = ANY_HDU;
	fitsfile *newfptr;         // FITS file pointer for the new image
	const long fpixel = 1;
	const int morekeys = 100;
	
	// remove existing file - cfitsio returns an error if it exists...
	remove(newFilename.c_str());
    // creates a fits file with filename newFilename.c_str()
	if (fits_create_file(&newfptr, newFilename.c_str(), &status)) {
		throw operaException("operaWIRCamImage: cfitsio create error: "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
    // if fptr == NULL, image was only created in memory
    if (fptr == NULL) {
        // Create headers (including primary) and images
		if (fits_create_img(newfptr, bitpix, 0, naxes, &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
        // copy the current HDU from file associated with fptr and appends to the
        // end of the file associated with newfptr (in this case the end is the primary)
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {	// use the HDU
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        // DEFAULT HEADER KEYWORDS FOR WIRCAM
        if (fits_write_key(newfptr, tstring, "DETECTOR", (char *)WIRCAM_DETECTOR, NULL, &status)) {
            throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
		fits_delete_key(newfptr, "DETSIZE", &status);
        if (fits_write_key(newfptr, tstring, "DETSIZE", (char *)WIRCAM_DETSIZE, (char *)"Total data pixels in full mosaic", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CCDBIN1", WIRCAM_CCDBIN1, -1, (char *)"Binning factor along first axis", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CCDBIN2", WIRCAM_CCDBIN2, -1, (char *)"Binning factor along second axis", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSIZE1", WIRCAM_PIXSIZE1, -1, (char *)"Pixel size for axis 1 (microns)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSIZE2", WIRCAM_PIXSIZE2, -1, (char *)"Pixel size for axis 2 (microns)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSCAL1", WIRCAM_PIXSCAL1, -1, (char *)"Pixel scale for axis 1 (arcsec/pixel)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSCAL2", WIRCAM_PIXSCAL2, -1, (char *)"Pixel scale for axis 2 (arcsec/pixel)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "IRFPA1", (char *)WIRCAM_EXTNAME1, (char *)"Name of the detector in ext 1", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA2", (char *)WIRCAM_EXTNAME2, (char *)"Name of the detector in ext 2", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA3", (char *)WIRCAM_EXTNAME3, (char *)"Name of the detector in ext 3", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA4", (char *)WIRCAM_EXTNAME4, (char *)"Name of the detector in ext 4", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_update_key_flt(newfptr, "CHIPBIAS", ChipBias, -1, (char *)"Science frame chip bias added to CDS value", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "DETSECA", (char *)WIRCAM_DETSECA, (char *)"Guider Mosaic Area of IRFPA1", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_write_key(newfptr, tstring, "DETSECB", (char *)WIRCAM_DETSECB, (char *)"Guider Mosaic Area of IRFPA2", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
        if (fits_write_key(newfptr, tstring, "DETSECC", (char *)WIRCAM_DETSECC, (char *)"Guider Mosaic Area of IRFPA3", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
        if (fits_write_key(newfptr, tstring, "DETSECD", (char*)WIRCAM_DETSECD, (char *)"Guider Mosaic Area of IRFPA4", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        // END OF DEFAULT HEADER KEYWORDS FOR WIRCAM
	}
#if 0
	// be sure we have the basic headers for floats if no header existed
	if (datatype == tfloat || datatype == tdouble) {
		float bzero = 0.0, bscale = 1.0;
		if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
#endif
	if (fits_update_key_lng(newfptr, "NEXTEND", (unsigned short)extensions, (char *)"Number of Extensions", &status)) {
		throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_update_key_lng(newfptr, "NAXIS", (unsigned short)0, (char *)"Number of axes", &status)) {
		throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_delete_key(newfptr, "NAXIS1", &status)) {
		throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (fits_delete_key(newfptr, "NAXIS2", &status)) {
		throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}
    // set the compression type for the file associated with newfptr
	if (fits_set_compression_type(newfptr, compression, &status)) {
		throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}    
	void *basepointer = pixptr;
	// save each extension with header
	for (unsigned extension=1; extension<=extensions; extension++) {
		// move to the specified HDU number  (select extension)
		if ( fits_movabs_hdu(fptr, cfitsioextension(0), &hdutype, &status) ) {
			throw operaException("operaWIRCamImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		if ( fits_movrel_hdu(fptr, extension, &hdutype, &status) ) {
			throw operaException("operaWIRCamImage: cfitsio error ", status, __FILE__, __FUNCTION__, __LINE__);
		}
		// copy the current HDU from the file associated with fptr to the file associated with newfptr
		if (fits_copy_hdu(fptr, newfptr, morekeys, &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
		}
		if (fits_update_key_lng(newfptr, "BITPIX", bitpix, (char *)"number of bits per data pixel", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
#if 0
		// be sure we have the basic headers for floats if no header existed
		if (datatype == tfloat || datatype == tdouble) {
			float bzero = 0.0, bscale = 1.0;
			if (fits_set_bscale(newfptr, bscale, bzero, &status)) {	// ??? doesn't seem to work...
				throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BZERO", bzero, -1, NULL, &status)) {
				throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (fits_update_key_flt(newfptr, "BSCALE", bscale, -1, NULL, &status)) {
				throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
		}
#endif
        // DEFAULT HEADER KEYWORDS FOR WIRCAM
        if (extension == 1){
            // EXTNAME and DETSEC change for each extension
            if (fits_write_key(newfptr, tstring, "EXTNAME", (char *)WIRCAM_EXTNAME1, (char *)"Extension name", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
            if (fits_update_key(newfptr, tstring, "DETSEC", (char *)WIRCAM_DETSEC1, (char *)"Mosaic area of the detector", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
        if (extension == 2){
            // EXTNAME and DETSEC change for each extension
            if (fits_write_key(newfptr, tstring, "EXTNAME", (char *)WIRCAM_EXTNAME2, (char *)"Extension name", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
            if (fits_update_key(newfptr, tstring, "DETSEC", (char *)WIRCAM_DETSEC2, (char *)"Mosaic area of the detector", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
        if (extension == 3){
            // EXTNAME and DETSEC change for each extension
            if (fits_write_key(newfptr, tstring, "EXTNAME", (char *)WIRCAM_EXTNAME3, (char *)"Extension name", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
            if (fits_update_key(newfptr, tstring, "DETSEC", (char *)WIRCAM_DETSEC3, (char *)"Mosaic area of the detector", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
        if (extension == 4){
            // EXTNAME and DETSEC change for each extension
            if (fits_write_key(newfptr, tstring, "EXTNAME", (char *)WIRCAM_EXTNAME4, (char *)"Extension name", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
            if (fits_update_key(newfptr, tstring, "DETSEC", (char *)WIRCAM_DETSEC4, (char *)"Mosaic area of the detector", &status)) {
                throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
            }
        }
        if (fits_write_key(newfptr, tstring, "DETECTOR", (char *)WIRCAM_DETECTOR, NULL, &status)) {
            throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
		fits_delete_key(newfptr, "DETSIZE", &status);
        if (fits_write_key(newfptr, tstring, "DETSIZE", (char *)WIRCAM_DETSIZE, (char *)"Total data pixels in full mosaic", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CCDBIN1", WIRCAM_CCDBIN1, -1, (char *)"Binning factor along first axis", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CCDBIN2", WIRCAM_CCDBIN2, -1, (char *)"Binning factor along second axis", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSIZE1", WIRCAM_PIXSIZE1, -1, (char *)"Pixel size for axis 1 (microns)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSIZE2", WIRCAM_PIXSIZE2, -1, (char *)"Pixel size for axis 2 (microns)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSCAL1", WIRCAM_PIXSCAL1, -1, (char *)"Pixel scale for axis 1 (arcsec/pixel)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "PIXSCAL2", WIRCAM_PIXSCAL2, -1, (char *)"Pixel scale for axis 2 (arcsec/pixel)", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "IRFPA1", (char *)WIRCAM_EXTNAME1, (char *)"Name of the detector in ext 1", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA2", (char *)WIRCAM_EXTNAME2, (char *)"Name of the detector in ext 2", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA3", (char *)WIRCAM_EXTNAME3, (char *)"Name of the detector in ext 3", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }
        if (fits_write_key(newfptr, tstring, "IRFPA4", (char *)WIRCAM_EXTNAME4, (char *)"Name of the detector in ext 4", &status)) {
            throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
        }    
        if (fits_write_key(newfptr, tstring, "CHIPSIZE", (char *)WIRCAM_CHIPSIZE, (char *)"Detector imaging area sizes", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "DATASEC", (char *)WIRCAM_DATASEC, (char *)"Imaging area of the detector", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_update_key_flt(newfptr, "CHIPBIAS", ChipBias, -1, (char *)"Science frame chip bias added to CDS value", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (fits_write_key(newfptr, tstring, "DETSECA", (char *)WIRCAM_DETSECA, (char *)"Guider Mosaic Area of IRFPA1", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (fits_write_key(newfptr, tstring, "DETSECB", (char *)WIRCAM_DETSECB, (char *)"Guider Mosaic Area of IRFPA2", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
        if (fits_write_key(newfptr, tstring, "DETSECC", (char *)WIRCAM_DETSECC, (char *)"Guider Mosaic Area of IRFPA3", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		} 
        if (fits_write_key(newfptr, tstring, "DETSECD", (char*)WIRCAM_DETSECD, (char *)"Guider Mosaic Area of IRFPA4", &status)) {
			throw operaException("operaWIRCamImage: cfitsio error ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
		}       
        // END DEFAULT HEADER KEYWORDS FOR WIRCAM         
		if (AllExtensions) {
			// write elements into the fits file
			if (fits_write_img(newfptr, datatype, fpixel, npixels_per_extension, basepointer, &status)) {
				throw operaException("operaWIRCamImage: cfitsio error "+newFilename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);
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
					throw operaException("operaWIRCamImage: cfitsio error "+newFilename+" ", (operaErrorCode)operaErrorCodeDatatypeNotSupported, __FILE__, __FUNCTION__, __LINE__);	
					break;
			}
		}
	}
	fits_close_file(newfptr, &status);
	newfptr = NULL;
	// note that fptr is left open in case the user wants to save again, it will be closed
    // by a close call if made, or the destructor otherwise
}

/*! 
 * applyNonLinearCorrection(unsigned extension, Polynomial &polynomial) 
 * \brief Applies the non-linear correction polynomial to chip.
 * \return void
 */
void operaWIRCamImage::applyNonLinearCorrection(unsigned extension, Polynomial &polynomial) {
	unsigned long n = this->npixels_per_extension;
	float *pixels = (float *)this->getpixels()+(extension-1)*n;
	while (n--) {
		*pixels = *pixels*polynomial.Evaluate(*pixels);
		pixels++;
	}
}

/* 
 * void detrend(operaWIRCamImage &image, operaWIRCamImage &dark, operaWIRCamImage &badpixelmask, operaWIRCamImage *flat)
 * \brief detrend a WIRCam image. all extensions in parallel
 */
void operaWIRCamImage::detrend(operaWIRCamImage &image, operaWIRCamImage &dark, operaWIRCamImage &badpixelmask, operaWIRCamImage *flat) {
	operaWIRCamImage &output = *this;
	if (flat == NULL) {
		output = image - dark;
	} else {
		output = (image - dark) / (*flat);
	}
	output.operaFITSSetHeaderValue("BPIXNAME", badpixelmask.operaFITSGetFilename(),	"Badpixelmask Name");
	output.operaFITSSetHeaderValue("BDPIXVAL", 0.0,					"Bad pixel value");
	output.operaFITSSetHeaderValue("DARKNAME", dark.operaFITSGetFilename(),	"Dark Name");
	output.operaFITSSetHeaderValue("DARKSUB",  "yes",				"Dark Subtraction done?");
	output.operaFITSSetHeaderValue("FLATCORR", "no",				"Flat field corrected (for horiz. lines)?");
	if (flat == NULL) {
		output.operaFITSSetHeaderValue("FLATDIV", "no", "Flat field division done?");
	} else {
		output.operaFITSSetHeaderValue("FLATDIV", "yes", "Flat field division done?");
	}
}

/* 
 * void weightmap(operaWIRCamImage &image, operaWIRCamImage &badpixelmask)
 * \brief create a weightmap of bad and saturated pixxels of a WIRCam image. all extensions in parallel
 */
void operaWIRCamImage::weightmap(operaWIRCamImage &image, operaWIRCamImage &badpixelmask) {
	operaWIRCamImage &output = *this;
	output = badpixelmask || (image > WIRCAM_SATURATION) || (image < -WIRCAM_SATURATION);
}

/* 
 * void skySubtraction(operaWIRCamImage &image, operaWIRCamImage &sky, float exptime)
 * \brief sky subtract a WIRCam image, set exptime to 1.0 if not wanted.
 * The sky is normalized to 1.0
 */
void operaWIRCamImage::skySubtraction(operaWIRCamImage &image, operaWIRCamImage &sky, float exptime) {
	operaWIRCamImage &output = *this;
	if (image.isMEFCube()) {
		image.medianCollapse();
	}
	float imageMedian = operaArrayMedian(image.getnpixels(), (float *)image.getpixels());
	output -= sky * imageMedian * exptime;
	output.operaFITSSetHeaderValue("SKYSUBTR", "yes", "Sky subtraction done?");
	output.operaFITSSetHeaderValue("SKYNAME", sky.operaFITSGetFilename(), "Sky name");
}

/* 
 * void maskGuideWindows(operaWIRCamImage &image)
 * \brief mask the entire rows and columns at the guide window positions in a  WIRCam image
 */
void operaWIRCamImage::maskGuideWindows(operaWIRCamImage &image) {
	operaFITSImage &output = *this;
	output = image;
	output[new operaImageVector(GuideWindows,4)] = 0.0;
}

/* 
 * void masterFlat(operaWIRCamImage &image[], operaFITSImage *weight, unsigned count)
 * \brief create a master flat from a medianCombined cube
 * \brief optionally creating a weight map
 * \param operaWIRCamImage images is an array of pointers to the flats images
 * \param operaWIRCamImage images is a pointer to the weight map image
 * \param unsigned count is the number of twilight flats to be combined
 */
void operaWIRCamImage::masterFlat(operaWIRCamImage *images[], operaWIRCamImage *weight, unsigned count) {
	operaWIRCamImage &output = *this;
	/*
	 * First, median combine. Most likely this would have already
	 * been done, so this is a no-op, but it needs to be done...
	 * Next subtract the chipbias
	 * Next calculate the median and median deviation
	 */
	const unsigned maxx = images[0]->getXDimension();
	const unsigned maxy = images[0]->getYDimension();
	float *median = new float[count];
	float *mediandeviation = new float[count];
	for (unsigned i=0; i<count; i++) {
		operaWIRCamImage *image = images[i];
		image->medianCollapse();
		*image -= ChipBias;
		median[i] = operaArrayMedian(maxx*maxy, (float *)(image->getpixels()));
		mediandeviation[i] = operaArrayMedianSigma(maxx*maxy, (float *)(image->getpixels()), median[i]);
	}
	
	/* 
	 * For each pixel, fit a straight line in the plot median(flux) vs. flux[y][x]
	 */
	float *pixelStack = new float[count];
	float *xx = new float[count];
	float *yy = new float[count];
	float a = 0.0;	// offset
	float b = 0.0;	//slope 
	float absoluteDeviation = 0.0;
	float finalDeviation = 0.0;
	
	for (unsigned y=0;  y<maxy; y++) {
		for (unsigned x=0; x<maxx; x++) {
			for (unsigned i=0; i<count; i++) {
				operaFITSImage *image = (operaFITSImage *)images[i];
				pixelStack[i] = *image[y][x];
			}
			ladfit(median, pixelStack, count, &a, &b, &absoluteDeviation);
			// Now refine the fit
			float yfit = 0.0;
			float total = 0.0;
			unsigned j = 0;
			for (unsigned i=0; i<maxx; i++) {
				yfit = fabs((median[i] * b) + a - pixelStack[i]);
				if (yfit < FiveSigma*absoluteDeviation){
					xx[j] = median[i];
					yy[j] = pixelStack[i];
					j++;
					total += yfit;
				}
			}
			if (j < maxx && total != 0.0) {
				ladfit(xx, yy, j, &finalDeviation, &a, &b);
			} else {
				finalDeviation = absoluteDeviation;
			}
			
			// now assign the slope to the output
			((operaFITSImage &)output)[y][x] = b;
			
			// and calculate the weight if asked
			if (weight) 
				((operaFITSImage &)*weight)[y][x] = 1.0 / pow(finalDeviation/0.8, 2);
		}
	}
	delete[] xx;
	delete[] yy;
	delete[] pixelStack;
	delete[] median;
	delete[] mediandeviation;
}

/* 
 * void masterFlatParallel(operaWIRCamImage &image[], operaFITSImage *weight, unsigned count)
 * \brief create a master flat from a medianCombined cube, uses threading to speed up the process
 * \brief optionally creating a weight map
 * \param operaWIRCamImage images is an array of pointers to the flats images
 * \param operaWIRCamImage images is a pointer to the weight map image
 * \param unsigned count is the number of twilight flats to be combined
 */
void operaWIRCamImage::masterFlatParallel(operaWIRCamImage *images[], operaWIRCamImage *weight, unsigned count) {
	operaWIRCamImage &output = *this;
	/*
	 * First, median combine. Most likely this would have already
	 * been done, so this is a no-op, but it needs to be done...
	 * Next subtract the chipbias
	 * Next calculate the median and median deviation
	 */
	unsigned int extensions = images[0]->getNExtensions();
	const unsigned maxx = images[0]->getXDimension();
	const unsigned maxy = images[0]->getYDimension();
	float *median = new float[count];
	float *mediandeviation = new float[count];
	for (unsigned i=0; i<count; i++) {
		operaWIRCamImage *image = images[i];
		image->medianCollapseParallel(extensions);
		*image -= ChipBias;
		median[i] = operaArrayMedian(maxx*maxy, (float *)(image->getpixels()));
		mediandeviation[i] = operaArrayMedianSigma(maxx*maxy, (float *)(image->getpixels()), median[i]);
	}
	
	/* 
	 * For each pixel, fit a straight line in the plot median(flux) vs. flux[y][x]
	 */
	float *pixelStack = new float[count];
	float *xx = new float[count];
	float *yy = new float[count];
	float a = 0.0;	// offset
	float b = 0.0;	//slope
	float absoluteDeviation = 0.0;
	float finalDeviation = 0.0;
	
	for (unsigned y=0;  y<maxy; y++) {
		for (unsigned x=0; x<maxx; x++) {
			for (unsigned i=0; i<count; i++) {
				operaFITSImage *image = (operaFITSImage *)images[i];
				pixelStack[i] = *image[y][x];
			}
			ladfit(median, pixelStack, count, &a, &b, &absoluteDeviation);
			// Now refine the fit
			float yfit = 0.0;
			float total = 0.0;
			unsigned j = 0;
			for (unsigned i=0; i<maxx; i++) {
				yfit = fabs((median[i] * b) + a - pixelStack[i]);
				if (yfit < FiveSigma*absoluteDeviation){
					xx[j] = median[i];
					yy[j] = pixelStack[i];
					j++;
					total += yfit;
				}
			}
			if (j < maxx && total != 0.0) {
				ladfit(xx, yy, j, &finalDeviation, &a, &b);
			} else {
				finalDeviation = absoluteDeviation;
			}
			
			// now assign the slope to the output
			((operaFITSImage &)output)[y][x] = b;
			
			// and calculate the weight if asked
			if (weight) 
				((operaFITSImage &)*weight)[y][x] = 1.0 / pow(finalDeviation/0.8, 2);
		}
	}
	delete[] xx;
	delete[] yy;
	delete[] pixelStack;
	delete[] median;
	delete[] mediandeviation;
}

/* 
 * void medianStack(operaWIRCamImage *images[], unsigned count)
 * \brief create a master dark by median combining the stack
 * \param operaWIRCamImage images is an array of pointers to the darks images
 * \param unsigned count is the number of darks to be combined
 */
void operaWIRCamImage::medianStack(operaWIRCamImage *images[], unsigned count) {
	operaFITSImage &output = *this;
	/*
	 * First, median combine. Most likely this would have already
	 * been done, so this is a no-op, but it needs to be done...
	 * This means slices == 1
	 */
	for (unsigned i=0; i<count; i++) {
		operaWIRCamImage *image = images[i];
		image->medianCollapse();
	}
	float *pixelStack = new float[count];
	const unsigned maxx = images[0]->getXDimension();
	const unsigned maxy = images[0]->getYDimension();
	for (unsigned y=0; y<maxy; y++) {
		for (unsigned x=0; x<maxx; x++) {
			for (unsigned i=0; i<count; i++) {
				operaFITSImage *image = (operaFITSImage *)images[i];
				pixelStack[i] = *image[y][x];
			}
			output[y][x] = operaArrayMedianQuick(count, (float *)pixelStack);
		}
	}
	delete[] pixelStack;
}

/* 
 * void medianStackParallel(operaWIRCamImage *images[], unsigned count)
 * \brief create a master dark by median combining the stack, uses threading to collapse individual darks
 * \param operaWIRCamImage images is an array of pointers to the darks images
 * \param unsigned count is the number of darks to be combined
 */
void operaWIRCamImage::medianStackParallel(operaWIRCamImage *images[], unsigned count) {
	operaFITSImage &output = *this;
	/*
	 * First, median combine. Most likely this would have already
	 * been done, so this is a no-op, but it needs to be done...
	 * This means slices == 1
	 */
	unsigned int extensions = images[0]->getNExtensions();
	for (unsigned i=0; i<count; i++) {
		operaWIRCamImage *image = images[i];
		image->medianCollapseParallel(extensions); // uses threading
	}
	float *pixelStack = new float[count];
	const unsigned maxx = images[0]->getXDimension();
	const unsigned maxy = images[0]->getYDimension();
	for (unsigned y=0; y<maxy; y++) {
		for (unsigned x=0; x<maxx; x++) {
			for (unsigned i=0; i<count; i++) {
				operaFITSImage *image = (operaFITSImage *)images[i];
				pixelStack[i] = *image[y][x];
			}
			output[y][x] = operaArrayMedianQuick(count, (float *)pixelStack);
		}
	}
	delete[] pixelStack;
}

/* 
 * void createReferencePixelImage(operaWIRCamImage &image)
 * \brief Create a bias image from a WIRCam image.
 * \details A band of pixels from the top and bottom edge of
 * and image is used as a bias estimate. The median value of
 * the top and bottom band is calculated and then propagated
 * through an entire column. Each column through the entire
 * image (or cube) is propagated independently. Each extension
 * is calculated independently.
 */
#define XBAND 8
void operaWIRCamImage::createReferencePixelImage(operaWIRCamImage &image) {
	operaWIRCamImage &output = *this;
	const unsigned edgeRows[XBAND] = WIRCAM_EDGE_ROWS;
	const unsigned maxx = image.getXDimension();
	const unsigned maxy = image.getYDimension();
	const unsigned maxz = image.getZDimension();
	float medianvector[XBAND];
	// do an extension at a time, and process slices at the same time
	for (unsigned extension=1; extension<=image.getNExtensions(); extension++) {
		operaFITSImage inimage(image[extension], true);		// true means a view not a copy
		operaFITSImage outfitsimage(output[extension], true);	// true means a view not a copy
		// for each row in the EDGEROWS, replace the individual row values
		// with the median of the whole row...
		for (unsigned i=0; i<XBAND; i++) {
			unsigned y = edgeRows[i];
			float *row = inimage[y];
			float rowmedian = operaArrayMedianQuick(maxx, row);
			for (unsigned x=0; x<maxx; x++) {
				*row++ -= rowmedian;
			}
		}
		// Now, take the median of the edge rows
		// and set the entire row to the median
		for (unsigned x=0; x<maxx; x++) {
			for (unsigned i=0; i<XBAND; i++) {
				unsigned y = edgeRows[i];
				medianvector[i] = inimage[y][x];
			}
			float median = operaArrayMedianQuick(XBAND, medianvector);// destructive version
			for (unsigned y=0; y<maxy*maxz; y++) {				// Note maxy*maxz does all slices...
				outfitsimage[y][x] = median;
			}
		}
		// Now, take a flaoting average across the columns
		// to smooth the result
		const unsigned floatwidth = REFERENCE_PIXEL_MOVING_AVERAGE_WIDTH;
		unsigned int y = 0;
		float *row = outfitsimage[y];
		for (unsigned x=0; x<maxx; x++) {				// Note maxy*maxz does all slices...
			float rowmovingaverage = 0.0;
			unsigned int count = 0;
			for (unsigned k=0; k<floatwidth; k++) {		// Note maxy*maxz does all slices...
				if ((x+k) < maxx) {
					rowmovingaverage += row[x+k];
					count++;
				}
			}
			if (count > 0) {
				rowmovingaverage /= count;
				for (unsigned yy=0; yy<maxy*maxz; yy++) {				// Note maxy*maxz does all slices...
					outfitsimage[yy][x] = rowmovingaverage;
				}
			}
		}
		inimage.operaFITSImageClose();
		outfitsimage.operaFITSImageClose();
	}
	output.setExtension(1);		// indexing above set the extension to 4
}

/* 
 * void calculateSkyBackground(operaWIRCamImage &image)
 * \brief Calculate the Sky Background - the median / exposure time
 * Could also be gotten from the headers SKYLVL, SKYDEV -- if they exist
 */
void operaWIRCamImage::calculateSkyBackground(operaWIRCamImage &image) {
	const unsigned maxx = image.getXDimension();
	const unsigned maxy = image.getYDimension();
	for (unsigned extension=1; extension<=image.getNExtensions(); extension++) {
		for (unsigned slice=1; slice<=getZDimension(); slice++) {
			skyBackground[extension-1][slice-1].skyLevel = operaArrayMedian(maxx*maxy, (float *)(image.getpixels()));
			skyBackground[extension-1][slice-1].skyDeviation = operaArrayMedianSigma(maxx*maxy, (float *)(image.getpixels()), skyBackground[extension-1][slice-1].skyLevel);
			skyBackground[extension-1][slice-1].skyRate = skyBackground[extension-1][slice-1].skyLevel / ExposureTime;
		}
	}
}

/* 
 * void calculateQERatio(operaWIRCamImage *images[], unsigned count)
 * \brief Calculate the QE Ratio - the median sky rate over the images
 * Could also be gotten from the headers SKYLVL, SKYDEV -- if they exist
 */
void operaWIRCamImage::calculateQERatio(operaWIRCamImage *images[], unsigned count) {
	float *qeratio = new float[count];
	float sum = 0.0;
	const unsigned extensions = images[0]->getNExtensions();
	for (unsigned extension=1; extension<=extensions; extension++) {
		for (unsigned slice=1; slice<=getZDimension(); slice++) {
			for (unsigned i=0; i<count; i++) {
				operaWIRCamImage *image = (operaWIRCamImage *)images[i];
				qeratio[i] = image->skyBackground[extension-1][slice-1].skyRate;
			}
			skyBackground[extension-1][slice-1].QERatio = operaArrayMedianQuick(count, (float *)qeratio);
			sum += skyBackground[extension-1][slice-1].QERatio;
		}
	}
	sum /= extensions;
	for (unsigned extension=1; extension<=extensions; extension++) {
		for (unsigned slice=1; slice<=getZDimension(); slice++) {
			skyBackground[extension-1][slice-1].QERatio /= sum;
		}
	}
	delete[] qeratio;
}

/* 
 * void createSky(operaWIRCamImage *images[], unsigned count)
 * \brief Create a master sky image
 */
void operaWIRCamImage::createSky(operaWIRCamImage *images[], unsigned count) {
	operaWIRCamImage &output = *this;
	
	// First, do a variance correction
	for (unsigned i=0; i<count; i++) {
		operaWIRCamImage *image = (operaWIRCamImage *)images[i];
		for (unsigned extension=1; extension<=extensions; extension++) {
			for (unsigned slice=1; slice<=getZDimension(); slice++) {
				(operaFITSImage &)image[extension][slice] /= image->ExposureTime * image->skyBackground[extension-1][slice-1].skyDeviation;
			}
		}
	}
	
	// next take a median stack
	float *pixelStack = new float[count];
	const unsigned maxx = images[0]->getXDimension();
	const unsigned maxy = images[0]->getYDimension();
	for (unsigned y=0; y<maxy; y++) {
		for (unsigned x=0; x<maxx; x++) {
			for (unsigned i=0; i<count; i++) {
				operaFITSImage *image = (operaFITSImage *)images[i];
				pixelStack[i] = (*image)[y][x];
			}
			output[y][x] = operaArrayMedianQuick(count, (float *)pixelStack);
		}
	}
	delete[] pixelStack;
	
	// next normalize by the median calue
	for (unsigned extension=1; extension<=extensions; extension++) {
		for (unsigned slice=1; slice<=getZDimension(); slice++) {
			float skyMedian = operaArrayMedian(maxx*maxy, (float *)(output[extension][slice].getpixels()));
			output[extension][slice] /= skyMedian;
		}
	}
	
}
/* 
 * void correctSkyLevel(operaWIRCamImage *images[], unsigned count)
 * \brief Correct sky level
 */
void operaWIRCamImage::correctSkyLevel(operaWIRCamImage *images[], unsigned count) {
	float *correctedMean = new float[count];
	float *correctedVariance = new float[count];
	const unsigned extensions = images[0]->getNExtensions();
	float total = 0.0;
	for (unsigned i=0; i<count; i++) {
		operaWIRCamImage *image = (operaWIRCamImage *)images[i];
		float sum = 0.0;
		for (unsigned extension=1; extension<=extensions; extension++) {
			for (unsigned slice=1; slice<=getZDimension(); slice++) {
				sum += skyBackground[extension-1][slice-1].skyRate / image->skyBackground[extension-1][slice-1].QERatio;
			}
		}
		correctedMean[i] = sum / extensions;
		total += correctedMean[i];
	}
	float mean = total / count;
	for (unsigned i=0; i<count; i++) {
		correctedVariance[i] = correctedMean[i] / mean;
	}
}
/* 
 * Box &getGuideWindow(unsigned index)
 * \brief Retrieve the nth Guide Window Box, may be indexed by chip number.
 */
Box *operaWIRCamImage::getGuideWindow(unsigned index) {
	switch (index) {
		case ChipOne:
		case 1:
			index = 0;
			break;
		case ChipTwo:
		case 2:
			index = 1;
			break;
		case ChipThree:
		case 3:
			index = 2;
			break;
		case ChipFour:
		case 4:
			index = 4;
			break;
		default:
			throw operaException("operaWircamImage: "+itos(index)+" ", operaErrorInvalidGuideWindowIndex, __FILE__, __FUNCTION__, __LINE__);	
			break;
	}
	return &GuideWindows[index];
}
/*
 * edgeExtend(Box &b, int Direction)
 * \brief extend the edge of a box in the 'X' (HORIZONTAL), 'Y' (VERTICAL), or 'N' (N = NONE, no extension) direction
 *  to the edge of the image
 * \param Box &b
 * \param int Direction
 */
Box operaWIRCamImage::edgeExtend(Box &b, char Direction/*=X|Y|N*/) {
	if (Direction == 'X') {
		// extend edges in 'X' (HORIZONTAL) direction 
		b.setX1(0);
		b.setX2(naxis1);
	} else if (Direction == 'Y') {
		// extend edges in 'Y' (VERTICAL) direction
		b.setY1(0);
		b.setY2(naxis2);
	} else if (Direction == 'N') {
		// do not extend edges
		// do nothing, no extension required
	} else {
		// invalid input, throws excception
		throw operaException("operaWircamImage error: invalid direction, must be 'X' (HORIZONTAL), 'Y' (VERTICAL), or 'N' (NONE)",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
	}
	return b;
}
/* 
 * float getChipBias(void)
 * \brief Retrieve the WIRCam ChipBias.
 */
float operaWIRCamImage::getChipBias(void) {
	return ChipBias;
}

/* 
 * float getAverageSkyLevel(unsigned extension)
 * \brief Retrieve the average sky level for all extensions.
 */
float operaWIRCamImage::getAverageSkyLevel(void) {
	return (skyBackground[0][0].skyLevel+skyBackground[1][0].skyLevel+skyBackground[2][0].skyLevel+skyBackground[3][0].skyLevel)/4.0;
}

/* 
 * float getSkyLevel(unsigned extension)
 * \brief Retrieve the sky level for an extension.
 */
float operaWIRCamImage::getSkyLevel(unsigned extension) {
	return skyBackground[extension-1][0].skyLevel;
}

/* 
 * float getSkyRate(unsigned extension)
 * \brief Retrieve the sky rate for an extension.
 */
float operaWIRCamImage::getSkyRate(unsigned extension) {
	return skyBackground[extension-1][0].skyRate;
}

