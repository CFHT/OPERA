#ifndef OPERAWIRCAMIMAGE_H
#define OPERAWIRCAMIMAGE_H
/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
Class: operaWIRCamImage
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

/*
 * ASL11/ASL22 causes the FITSImage [] operator to use the << 11 operation when calculating
 * addresses - since the x dimension is known to be a fixed 2048 == 2^11
 * and since the y dimension is also 2048 ASL22 = naxis1*naxis2
 */
#define ASL11 11
#define ASL22 22

#include "globaldefines.h"							// for FITS_RANGE_CHECK
#include "libraries/Polynomial.h"					// Polynomial
#include "libraries/operaGeometricShapes.h"			// for Box
#include "libraries/operaImageVector.h"				// for operaImageVector
#include "libraries/operaFITSImage.h"				// for operaFITSImage
#include "libraries/operaMultiExtensionFITSImage.h"	// for operaMultiExtensionFITSImage
#include "libraries/operaMultiExtensionFITSCube.h"	// for operaMultiExtensionFITSCube


#define WIRCAM_CHIPBIAS 7000.0
#define WIRCAM_SATURATION 35000.0

#define WIRCAM_CCDBIN1 1
#define WIRCAM_CCDBIN2 1
#define WIRCAM_PIXSIZE1 18.0
#define WIRCAM_PIXSIZE2 18.0
#define WIRCAM_PIXSCAL1 0.3
#define WIRCAM_PIXSCAL2 0.3
#define WIRCAM_MAXPIXEL 65535

#define WIRCAM_CHIPSIZE "[1:2048,1:2048]"
#define WIRCAM_DATASEC "[1:2048,1:2048]"
#define WIRCAM_DETSIZE "[1:4096,1:4096]"

#define WIRCAM_DETSEC1 "[2194:4241,2194:4241]"
#define WIRCAM_DETSEC2 "[2194:4241,1:2048]"
#define WIRCAM_DETSEC3 "[1:2048,1:2048]"
#define WIRCAM_DETSEC4 "[1:2048,2194:4241]"

#define WIRCAM_DETSECA "[15:28,15:28]"
#define WIRCAM_DETSECB "[15:28,1:14]"
#define WIRCAM_DETSECC "[1:14,1:14]"
#define WIRCAM_DETSECD "[1:14,15:28]"

#define WIRCAM_EXTNAME1 "HAWAII-2RG-#77"
#define WIRCAM_EXTNAME2 "HAWAII-2RG-#52"
#define WIRCAM_EXTNAME3 "HAWAII-2RG-#54"
#define WIRCAM_EXTNAME4 "HAWAII-2RG-#60"
#define WIRCAM_DETECTOR "WIRCam  "

#define WIRCAM_NAXIS1 2048
#define WIRCAM_NAXIS2 2048
#define WIRCAM_EXTENSIONS 4
#define WIRCAM_MAX_SLICES 18
#define WIRCAM_MINIMUM_EXPOSURE_TIME 2.375			// seconds

#define WIRCAM_EDGE_ROWS {0, 1, 2, 3, 2044, 2045, 2046, 2047} /* for referene pixel sibtraction */
#define REFERENCE_PIXEL_MOVING_AVERAGE_WIDTH 4

const int last_extension = primary_extension+WIRCAM_EXTENSIONS;
const unsigned ChipOne = 77;
const unsigned ChipTwo = 52;
const unsigned ChipThree = 54;
const unsigned ChipFour = 60;

/*
 * operaWIRCamImage class
 * \author Doug Teeple, Megan Tannock
 * \brief This class encapsulates the WIRCam image.
 * \file operaWIRCamImage.h
 * \ingroup libraries
 */
class SkyBackGround {
public:
	float skyLevel;
	float skyDeviation;
	float skyRate;
	float QERatio;
};

class operaWIRCamImage : public operaMultiExtensionFITSCube  {

private:
	Box GuideWindows[WIRCAM_EXTENSIONS];
	SkyBackGround skyBackground[WIRCAM_EXTENSIONS][WIRCAM_MAX_SLICES];
	
	float ExposureTime;
	float ChipBias;
	bool haveHeaderSkyLevels;
	
public:
	/*!
	 * \class operaWIRCamImage()
	 * \brief Basic operaWIRCamImage constructor.
	 * \return void
	 */
	operaWIRCamImage();     // simply construct a default image

	~operaWIRCamImage();    // destructor	

	/*!
	 * \class operaWIRCamImage
	 * \brief construct an in-memory FITSImage object
	 * \brief operaWIRCamImage(unsigned Naxis1, unsigned Naxis2, edatatype Datatype=tushort, unsigned Compression=0)
	 * \brief Create an in-memory image of given dimensions.
	 * \param Naxis1 - x ccd dimension
	 * \param Naxis2 - y ccd dimension
	 * \param Datatype optional datatype defaults to tshort
	 * \param Compression optional compression, defaults to none
	 * \return void
	 */
	operaWIRCamImage(unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, unsigned Extensions, edatatype Datatype);	// simply construct an image of a given size
	/*! 
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
	operaWIRCamImage(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, unsigned Extensions, 
				   edatatype Datatype, unsigned Compression = 0, bool isLazy = true);	// create a new empty FITSImage
	/*! 
	 * \class operaWIRCamImage
	 * \brief create a FITSIMage object from a FITS file
	 * \brief operaWIRCamImage(string Filename, int mode=READWRITE)
	 * \brief Constructor to create a FITSImage from a FITS file.
	 * \param Filename
	 * \param mode
	 * \return void
	 */
	operaWIRCamImage(string Filename, edatatype Datatype, int mode=READWRITE/*READONLY*/, unsigned Compression = 0, bool isLazy = true);		// read an existing FITSImage from file
	/*! 
	 * \class operaWIRCamImage
	 * \brief create a operaWIRCamImage object from another operaWIRCamImage
	 * \brief operaWIRCamImage(operaWIRCamImage &image)
	 * \param operaWIRCamImage
	 * \return void
	 */
	operaWIRCamImage(operaWIRCamImage &image, bool ViewOnly = false, bool AddHeader = false);
	
    /**********************************************************************************/
    /**********************************************************************************/
    /**** NOTE WELL:                                                               ****/
    /**** The operators and Matrix View and column selector are defined for the    ****/
    /**** float datatype only.                                                     ****/
    /**********************************************************************************/
    /**********************************************************************************/
	
    /* OPERATORS */
    
    /*! 
	 * \brief operator []
	 * \brief indexing operator to return an operaFITSImage&.
	 * \param i - index into extension
	 * \return operaFITSImage& 
	 */
	operaFITSCube& operator[](unsigned i) {
#ifdef FAST_EXTENSIONS
		setExtension((unsigned long)i);
		npixels = npixels_per_extension;
		return (operaFITSCube&)*this;
#else
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (i > this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, true);
		t->npixels = npixels_per_extension;
		t->naxis = 3;
		t->naxes[0] = t->naxis1;
		t->naxes[1] = t->naxis2;
		t->naxes[2] = t->naxis3;
		t->current_extension = (unsigned long)i;
		t->current_slice = 1;
		t->extensions = extensions;
        setExtension(current_extension);
		t->super = this;							// record where we came from
		return (operaFITSCube&)*t;			
#endif
    };
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return an ImageVector.
	 * \param i - index into row
	 * \note usage: myfitsimage[where(input>saturation)] = 0.0; sets the image pixel values where true to zero
	 * \return float pointer to row 
	 */
	operaImageVector *operator[](operaImageVector *vector) {
        SetImage(vector, this); 
        return vector;
    };
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return an ImageVector.
	 * \param i - index into row
	 * \note usage: myfitsimage[where(input>saturation)] = 0.0; sets the image pixel values where true to zero
	 * \return float pointer to row 
	 */
	operaImageVector &operator[](operaImageVector &vector) {
        SetImage(&vector, this);
        return vector;
    };
    /*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b; copies the pixel values from b to a
	 * \return operaFITSImage&
	 */	
	operaWIRCamImage& operator=(operaWIRCamImage* b) {
        float *p = (float *)pixptr; 
        float *bp = (float *)b->pixptr;
		float *vp = (float *)varptr; 
		float *bvp = (float *)b->varptr; 
		unsigned long n = npixels, nv = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b->extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b->naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
			// we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
				setExtension(ext);
				b->setExtension(ext);
				b->extensionHasBeenRead[ext] = false;
				b->readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b->pixptr;
                while (n--) *p++ = *bp++;
				if (vp && bvp) {
					while (nv--) *vp++ = *bvp++;
				}
				saveExtension(ext);
                if (b->super)
					b->super->setHasBeenRead(ext);
                if (super)
					super->setHasBeenWritten(ext);
            }
		} else {
			// we need to copy all extensions that are stacked in memory
            while (n--) *p++ = *bp++;
			if (vp && bvp) {
				while (nv--) *vp++ = *bvp++;
			}
			npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);	// in case we were indexed...
        }
		if (b->istemp && !b->viewOnly) delete b;
		return *this;
	};
	/*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b; copies the pixel values from b to a
	 * \return operaFITSImage&
	 */	
	operaWIRCamImage& operator=(operaWIRCamImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		float *vp = (float *)varptr; 
		float *bvp = (float *)b.varptr; 
		unsigned long n = npixels, nv = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
				setExtension(ext);
				b.setExtension(ext);
				b.extensionHasBeenRead[ext] = false;
				b.readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
                while (n--) *p++ = *bp++;
				if (vp && bvp) {
					while (nv--) *vp++ = *bvp++;
				}
                saveExtension(ext);
                if (b.super)
					b.super->setHasBeenRead(ext);
                if (super)
					super->setHasBeenWritten(ext);
            }
		} else {
            while (n--) *p++ = *bp++;
			if (vp && bvp) {
				while (nv--) *vp++ = *bvp++;
			}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = 0.0; copies the float value to every pixel in = a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator=(float f) {
		float *p = (float *)pixptr;
		float *vp = (float *)varptr; 
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                //readExtension(ext);
				setExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				vp = (float *)varptr;
                while (n--) *p++ = f;
				if (vp) {
					while (nv--) *vp++ = 0.0;
				}
                saveExtension(ext);
                if (super)
					super->setHasBeenWritten(ext);
            }
		} else {
            while (n--) *p++ = f;
			if (vp) {
				while (nv--) *vp++ = 0.0;
			}
        }
		return *this;
	};
	/*! 
	 * \brief Median operator %
	 * \note usage: operaMultiExtensionFITSImage a % 8; median combines 8 slices of a
	 */	
	operaMultiExtensionFITSImage& operator%(unsigned i) {
		this->medianCollapseParallel(extensions, i);
		if (varptr) {
			for (unsigned k=0; k<i; k++) {
				float *vp = (float *)varptr; 
				unsigned long n = npixels; 
				while (n--) { *vp = sqrt(*vp); vp++; }
			}
		}
		return *this;
	};
    /*! 
	 * \brief operator >>
	 * \brief write image to named file (saveas).
	 * \param const string filename
	 * \note usage: operaWIRcamimage a(); a >> "foo.fits"; 
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator>>(const string Filename) {
		operaWIRCamImageSaveAs(Filename);
		return *this;
	};
	/*! 
	 * \brief operator >>
	 * \brief write image to another image (headers and pixels).
	 * \param operaWIRCamImage& image
	 * \note usage: operaWIRcamimage a(); operaWIRcamimage b(); a >> b; 
	 * \return operaWIRCamImage&
	 */	
	const operaWIRCamImage& operator>>(operaWIRCamImage& image) {
		image.operaMultiExtensionFITSImageCopyHeader(this);
		image = this;
		return image;
	};
	/*! 
	 * \brief operator <<
	 * \brief read an image from a named file.
	 * \param const string filename
	 * \note usage: operaWIRcamimage a(); a << "foo.fits"; 
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator<<(const string Filename) {
		filename = Filename;
		naxes[0] = naxis1;
		naxes[1] = naxis2;
		naxes[2] = naxis3;
		if (isLazy) {
			setHasBeenRead(0, false);
			openFITSfile(Filename, READONLY);
			readFITSHeaderInfo();
			if (!isLazy) {
				AllExtensions = true;
				AllSlices = true;
				readFITSArray();					
			}
		}
		return *this;
	};
	/*! 
	 * \brief operator <<
	 * \brief write image to another image (headers and pixels).
	 * \param operaWIRCamImage& image
	 * \note usage: operaWIRcamimage a(); operaWIRcamimage b(); a << b; 
	 * \return operaWIRCamImage&
	 */	
	const operaWIRCamImage& operator<<(operaWIRCamImage& image) {
		operaMultiExtensionFITSImageCopyHeader(&image);
		*this = image;
		return *this;
	};
    /*! 
	 * \brief operator ==
	 * \brief equality operator.
	 * \param b - operaWIRCamImage&
	 * \note usage:
	 * \return operaWIRCamImage&
     */
    operaWIRCamImage& operator==(operaWIRCamImage& b) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		float *tp;
 		operaWIRCamImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr;
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels;
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
  				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
                while (n--) *tp++ = *p++ * *bp++;
                saveExtension(ext);
            }
        } else {
            while (n--) {*tp++ = (*p++ == *bp++? 1.0: 0.0);}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
    /*! 
     * \brief operator ==
     * \brief equality operator.
     * \param b - operaWIRCamImage&
     * \note usage:
     * \return operaWIRCamImage&
     */
    operaWIRCamImage& operator==(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
 		unsigned long n = npixels;
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
  				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
				while (n--) {*tp++ = (*p++ == f?1.0:0.0);}
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ * f;
        }
		return *t;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a += operaWIRCamImage b; adds the pixel values from b to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator+=(operaWIRCamImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		float *vp = (float *)varptr; 
		float *bvp = (float *)b.varptr; 
		unsigned long n = npixels, nv = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
                while (n--) *p++ += *bp++;
				if (vp && bvp) {
					while (nv--) *vp++ += *bvp++;
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ += *bp++;
			if (vp && bvp) {
				while (nv--) *vp++ += *bvp++;
			}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a += 100.0; adds the float value to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator+=(float f) {
		float *p = (float *)pixptr; 
		float *vp = (float *)varptr; 
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr;
				vp = (float *)varptr;
                while (n--) *p++ += f;
				if (vp) {
					while (nv--) *vp++ = 0.0;
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ += f;
			if (vp) {
				while (nv--) *vp++ = 0.0;
			}
        }
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a -= operaWIRCamImage b; subtracts the pixel values from b to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator-=(operaWIRCamImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr;
		float *vp = (float *)varptr; 
		float *bvp = (float *)b.varptr; 
		unsigned long n = npixels, nv = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
				n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
				while (n--) *p++ -= *bp++;
				if (vp && bvp) {
					while (nv--) *vp++ += *bvp++;
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ -= *bp++;
			if (vp && bvp) {
				while (nv--) *vp++ += *bvp++;
			}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a -= 100.0; subtracts the float value from a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator-=(float f) {
		float *p = (float *)pixptr;
		float *vp = (float *)varptr; 
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				vp = (float *)varptr;
                while (n--) *p++ -= f;
				if (vp) {
					while (nv--) *vp++ = 0.0;
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ -= f;
			if (vp) {
				while (nv--) *vp++ = 0.0;
			}
        }
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a *= operaWIRCamImage b; multiplies the pixel values from b to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator*=(operaWIRCamImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr;
		float *vp = (float *)varptr; 
		float *bvp = (float *)b.varptr; 
		unsigned long n = npixels, nv = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
                while (n--) *p++ *= *bp++;
				if (vp && bvp) {
					p = (float *)pixptr; 
					bp = (float *)b.pixptr;
					while (nv--) { double td = *vp; *vp++ = ( pow(*p,2) * *bvp++ ) + ( pow(*bp,2) * td ); p++; bp++; }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ *= *bp++;
			if (vp && bvp) {
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				while (nv--) { double td = *vp; *vp++ = ( pow(*p,2) * *bvp++ ) + ( pow(*bp,2) * td ); p++; bp++; }
			}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a *= 100.0; multiplies the float value times a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator*=(float f) {
		float *p = (float *)pixptr;
		float *vp = (float *)varptr; 
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr;
				vp = (float *)varptr;
                while (n--) *p++ *= f;
				if (vp) {
					p = (float *)pixptr; 
					while (nv--) { *vp++ *= pow(f,2); }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ *= f;
			if (vp) {
				p = (float *)pixptr; 
				while (nv--) { *vp++ *= pow(f,2); }
			}
        }
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a /= operaWIRCamImage b; divides the pixel values from b into a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator/=(operaWIRCamImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr;
		float *vp = (float *)varptr; 
		float *bvp = (float *)b.varptr; 
		unsigned long n = npixels, nv = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
                while (n--) *p++ /= *bp++;
				if (vp && bvp) {
					p = (float *)pixptr; 
					bp = (float *)b.pixptr;
					while (nv--) { double td = *vp; *vp++ = ( pow(*bp,-2) * td) + ( pow(*p / pow(*bp,2),2) * *bvp++ ); p++; bp++; }
				}
                saveExtension(ext);
            }
        } else {
            while (n--) *p++ /= *bp++;
			if (vp && bvp) {
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
                while (nv--) { double td = *vp; *vp++ = ( pow(*bp,-2) * td) + ( pow(*p / pow(*bp,2),2) * *bvp++ ); p++; bp++; }
			}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a /= 100.0; divides the float value from a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator/=(float f) {
		float *p = (float *)pixptr;
		float *vp = (float *)varptr; 
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr;
				vp = (float *)varptr;
                while (n--) *p++ /= f;
				if (vp) {
					p = (float *)pixptr; 
					while (nv--) { *vp++ /= pow(f,2); }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ /= f;
			if (vp) {
				p = (float *)pixptr; 
				while (nv--) { *vp++ /= pow(f,2); }
			}
        }
		return *this;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b * operaWIRCamImage c; multiplies the pixel values  b * c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator*(operaWIRCamImage& b) {
		float *tp, *tvp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
			tvp = (float *)this->varptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
			tvp = (float *)t->varptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr;
		float *vp = (float *)varptr; 
		float *bvp = (float *)b.varptr; 
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
				tvp = (float *)t->varptr; 
                while (n--) *tp++ = *p++ * *bp++;
				if (vp && bvp && tvp) {
					p = (float *)pixptr; 
					bp = (float *)b.pixptr;
					while (nv--) { double td = *vp; *tvp++ = ( pow(*p,2) * *bvp++ ) + ( pow(*bp,2) * td ); p++; bp++; }
				}
                saveExtension(ext);
            }
        } else {
            while (n--) *tp++ = *p++ * *bp++;
 			if (vp && bvp && tvp) {
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				while (nv--) { double td = *vp; *tvp++ = ( pow(*p,2) * *bvp++ ) + ( pow(*bp,2) * td ); p++; bp++; }
			}
		}
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b * 10.0; multiplies the pixel values  b * 10.0 and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator*(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *tvp = (float *)t->varptr; 
		float *p = (float *)this->pixptr;
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
				tvp = (float *)t->varptr;
                while (n--) *tp++ = *p++ * f;
				if (tvp) {
					while (nv--) { *tvp++ *= pow(f,2); }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ * f;
 			if (tvp) {
				while (nv--) { *tvp++ *= pow(f,2); }
			}
        }
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b / operaWIRCamImage c; divides the pixel values  b / c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator/(operaWIRCamImage& b) {
		float *tp, *tvp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tvp = (float *)this->varptr; 
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
			tvp = (float *)t->varptr; 
		}
		float *p = (float *)this->pixptr; 
		float *vp = (float *)this->varptr; 
		float *bp = (float *)b.pixptr;
		float *bvp = (float *)b.varptr;
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
				tvp = (float *)t->varptr; 
                while (n--) *tp++ = *p++ / *bp++;
				if (vp && bvp && tvp) {
					p = (float *)pixptr; 
					bp = (float *)b.pixptr;
					while (nv--) { *tvp++ = ( pow(*bp,-2) * *vp++ ) + ( pow(*p / pow(*bp,2),2) * *bvp++ ); p++; bp++;  }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ / *bp++;
 			if (vp && bvp && tvp) {
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
                while (nv--) { *tvp++ = ( pow(*bp,-2) * *vp++ ) + ( pow(*p / pow(*bp,2),2) * *bvp++ ); p++; bp++;  }
			}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b / 100.0; divides the pixel values  b / 100.0 and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator/(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *tvp = (float *)t->varptr; 
		float *p = (float *)this->pixptr;
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
				tvp = (float *)t->varptr;
                while (n--) *tp++ = *p++ / f;
				if (tvp) {
					while (nv--) { *tvp++ /= pow(f,2); }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ / f;
 			if (tvp) {
				while (nv--) { *tvp++ /= pow(f,2); }
			}
        }
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b + operaWIRCamImage c; adds the pixel values  b + c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator+(operaWIRCamImage& b) {
		float *tp, *tvp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
			tvp = (float *)this->varptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
			tvp = (float *)t->varptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr;
		float *vp = (float *)varptr; 
		float *bvp = (float *)b.varptr; 
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
				tvp = (float *)t->varptr; 
                while (n--) *tp++ = *p++ + *bp++;
				if (vp && bvp && tvp) {
					while (nv--) { *tvp++ = *vp++ + *bvp++; }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ + *bp++;
 			if (vp && bvp && tvp) {
				while (nv--) { *tvp++ = *vp++ + *bvp++; }
			}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b + 100.0; adds the pixel values  b + 100.0 and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator+(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *tvp = (float *)t->varptr; 
		float *p = (float *)this->pixptr;
		float *vp = (float *)this->varptr;
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
				tvp = (float *)t->varptr;
                while (n--) *tp++ = *p++ + f;
				if (tvp && vp) {
					while (nv--) { *tvp++ = *vp++ + 0.0; }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ + f;
 			if (tvp && vp) {
				while (nv--) { *tvp++ = *vp++ + 0.0; }
			}
        }
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b - operaWIRCamImage c; subtracts the pixel values  b - c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator-(operaWIRCamImage& b) {
		float *tp, *tvp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
			tvp = (float *)this->varptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
			tvp = (float *)t->varptr; 
		}
		float *p = (float *)this->pixptr; 
		float *vp = (float *)this->varptr; 
		float *bp = (float *)b.pixptr;
		float *bvp = (float *)b.varptr;
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
				vp = (float *)varptr; 
				bvp = (float *)b.varptr; 
				tvp = (float *)t->varptr; 
                while (n--) *tp++ = *p++ - *bp++;
				if (vp && bvp && tvp) {
					while (nv--) { *tvp++ = *vp++ + *bvp++; }
				}
                saveExtension(ext);
            }
        } else {
            while (n--) *tp++ = *p++ - *bp++;
 			if (vp && bvp && tvp) {
				while (nv--) { *tvp++ = *vp++ + *bvp++; }
			}
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b - 100.0; subtracts the pixel values  b - 100.0 and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator-(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *tvp = (float *)t->varptr; 
		float *p = (float *)this->pixptr;
		unsigned long n = npixels, nv = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = nv = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
				tvp = (float *)t->varptr;
                while (n--) *tp++ = *p++ - f;
				if (tvp) {
					while (nv--) { *tvp++ /= pow(f,2); }
				}
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ - f;
 			if (tvp) {
				while (nv--) { *tvp++ /= pow(f,2); }
			}
        }
		return *t;
	};
	/*! 
	 * \brief operator !
	 * \brief invert operator.
	 * \note usage: operaWIRCamImage a = !operaWIRCamImage b; inverts the pixel values of b and assigns to a (0.0 becomes 1.0, non-zero becomes 0.0)
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator!() {
		float *tp;
		operaWIRCamImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
        unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++==0.0?1.0:0.0);
                saveExtension(ext);
            }
        } else {
            while (n--) *tp++ = (*p++==0.0?1.0:0.0);
        }
		return *t;
	};
	/*! 
	 * \brief operator >
	 * \brief greater than operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b > operaWIRCamImage c; creates a mask of 1.0 if b > c or 0.0 if b <= c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator>(operaWIRCamImage& b) {
		float *tp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++>*bp++?1.0:0.0);
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = (*p++>*bp++?1.0:0.0);
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator >
	 * \brief greater than operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b > 100.0; creates a mask of 1.0 if b > 100.0 or 0.0 if b <= 100.0 and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator>(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++>f?1.0:0.0);
                saveExtension(ext);
            }
        } else {
            while (n--) *tp++ = (*p++>f?1.0:0.0);
        }
		return *t;
	};
	/*! 
	 * \brief operator >=
	 * \brief greater than or equal operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b >= operaWIRCamImage c; creates a mask of 1.0 if b >= c or 0.0 if b < c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator>=(operaWIRCamImage& b) {
		float *tp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
				n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
				while (n--) *tp++ = (*p++>=*bp++?1.0:0.0);
                saveExtension(ext);
            }
        } else {
            while (n--) *tp++ = (*p++>=*bp++?1.0:0.0);
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator >=
	 * \brief greater than or equal operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b >= 100.0; creates a mask of 1.0 if b >= 100.0 or 0.0 if b < 100.0 and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator>=(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++>=f?1.0:0.0);
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = (*p++>=f?1.0:0.0);
        }
		return *t;
	};
	/*! 
	 * \brief operator <
	 * \brief less than operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b < operaWIRCamImage c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator<(operaWIRCamImage& b) {
		float *tp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++<*bp++?1.0:0.0);
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = (*p++<*bp++?1.0:0.0);
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator <
	 * \brief less than operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b < operaWIRCamImage c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator<(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++<f?1.0:0.0);
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = (*p++<f?1.0:0.0);
        }
		return *t;
	};
	/*! 
	 * \brief operator <=
	 * \brief less than or equal operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b <= operaWIRCamImage c; creates a mask of 1.0 if b <= c or 0.0 if b > c and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator<=(operaWIRCamImage& b) {
		float *tp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++<=*bp++?1.0:0.0);
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = (*p++<=*bp++?1.0:0.0);
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator <=
	 * \brief less than or equal operator.
	 * \param f - float
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b <= 100.0; creates a mask of 1.0 if b <= 100.0 or 0.0 if b > 100.0 and assigns to a
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator<=(float f) {
		operaWIRCamImage *t = new operaWIRCamImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++<f?1.0:0.0);
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = (*p++<f?1.0:0.0);
        }
		return *t;
	};
	/*! 
	 * \brief operator &&
	 * \brief less than operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b && operaWIRCamImage c; creates a mask of 1.0 where b && c != 0.0
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator&&(operaWIRCamImage& b) {
		float *tp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++!=0.0&&*bp++!=0.0?1.0:0.0);
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = (*p++!=0.0&&*bp++!=0.0?1.0:0.0);
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	
	/*! 
	 * \brief operator ||
	 * \brief less than operator.
	 * \param b - operaWIRCamImage&
	 * \note usage: operaWIRCamImage a = operaWIRCamImage b && operaWIRCamImage c; creates a mask of 1.0 where b && c != 0.0
	 * \return operaWIRCamImage&
	 */	
	operaWIRCamImage& operator||(operaWIRCamImage& b) {
		float *tp;
		operaWIRCamImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaWIRCamImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaWIRCamImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaWIRCamImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				tp = (float *)t->pixptr;
                while (n--) *tp++ = (*p++!=0.0||*bp++!=0.0?1.0:0.0);
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = (*p++!=0.0||*bp++!=0.0?1.0:0.0);
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	
    /*! 
	 * operaWIRCamImageSave() 
	 * \brief Saves the current image to disk.
	 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
	 * \throws operaException operaErrorCodeNoFilename
	 * \throws operaException cfitsio error code
	 * \return void
	 */
	void operaWIRCamImageSave();
	
	/*! 
	 * operaWIRCamImageSaveAs(string newFilename) 
	 * \brief Saves the current image to disk, with the given filename.
	 * \param newFilename
	 * \throws operaException cfitsio error code
	 * \return void
	 */
	void operaWIRCamImageSaveAs(string newFilename);
	
	/*! 
	 * applyNonLinearCorrection(unsigned extension, Polynomial &polynomial) 
	 * \brief Applies the non-linerar correction polynomial to chip.
	 * \throws operaException cfitsio error code
	 * \return void
	 */
	void applyNonLinearCorrection(unsigned extension, Polynomial &polynomial);

	/*! 
	 * void detrend(operaWIRCamImage &output, operaWIRCamImage &image, operaWIRCamImage &dark, operaWIRCamImage &mask, operaWIRCamImage *flat)
	 * \brief detrend a WIRCam image
	 */
	void detrend(operaWIRCamImage &image, operaWIRCamImage &dark, operaWIRCamImage &mask, operaWIRCamImage *flat);
	
	/*! 
	 * void weightmap(operaWIRCamImage &image, operaWIRCamImage &badpixelmask)
	 * \brief create a weightmap of bad and saturated pixxels of a WIRCam image. all extensions in parallel
	 */
	void weightmap(operaWIRCamImage &image, operaWIRCamImage &badpixelmask);

	/* 
	 * void skySubtraction(operaWIRCamImage &image, operaWIRCamImage &sky, float exptime)
	 * \brief sky subtract a WIRCam image, set exptime to 1.0 if not wanted
	 */
	void skySubtraction(operaWIRCamImage &image, operaWIRCamImage &sky, float exptime);
		
	/*! 
	 * void maskGuideWindows(operaWIRCamImage &image)
	 * \brief mask the entire rows and columns at the guide window positions in a  WIRCam image
	 */
	void maskGuideWindows(operaWIRCamImage &image);
	
	/* 
	 * void masterFlat(operaWIRCamImage &image[], operaFITSImage *weight, unsigned count) 
	 * \brief create a master flat from a medianCombined cube
	 * \brief optionally creating a weight map
	 * \param operaWIRCamImage images is an array of pointers to the flats images
	 * \param operaWIRCamImage images is a pointer to the weight map image
	 * \param unsigned count is the number of twilight flats to be combined
	 */
	void masterFlat(operaWIRCamImage *images[], operaWIRCamImage *weight, unsigned count);
	
	/* 
	 * void masterFlatParallel(operaWIRCamImage &image[], operaFITSImage *weight, unsigned count)
	 * \brief create a master flat from a medianCombined cube, uses threading to speed up the process
	 * \brief optionally creating a weight map
	 * \param operaWIRCamImage images is an array of pointers to the flats images
	 * \param operaWIRCamImage images is a pointer to the weight map image
	 * \param unsigned count is the number of twilight flats to be combined
	 */
	void masterFlatParallel(operaWIRCamImage *images[], operaWIRCamImage *weight, unsigned count);
	
	/* 
	 * void medianStack(operaWIRCamImage *images[], unsigned count)
	 * \brief create a master dark by median combining the stack
	 * \param operaWIRCamImage images is an array of pointers to the darks images
	 * \param unsigned count is the number of darks to be combined
	 */
	void medianStack(operaWIRCamImage *images[], unsigned count);
	
	/* 
	 * void medianStackParallel(operaWIRCamImage *images[], unsigned count)
	 * \brief create a master dark by median combining the stack, uses threading  to collapse individual darks
	 * \param operaWIRCamImage images is an array of pointers to the darks images
	 * \param unsigned count is the number of darks to be combined
	 */
	void medianStackParallel(operaWIRCamImage *images[], unsigned count);
	
	/* 
	 * void createReferencePixelImage(operaWIRCamImage &image)
	 * \brief Create a bias image from a WIRCam image
	 * \details A band of pixels from the top and bottom edge of
	 * and image is used as a bias estimate. The median value of
	 * the top and bottom band is calculated and then propagated
	 * through an entire column. Each column through the entire
	 * image (or cube) is propagated independently. Each extension
	 * is calculated independently.
	 */
	void createReferencePixelImage(operaWIRCamImage &image);
	
	/*! 
	 * void calculateSkyBackground(operaWIRCamImage &image)
	 * \brief Calculate the Sky Background - the median / exposure time
	 * Could also be gotten from the headers SKYLVL, SKYDEV -- if they exist
	 */
	void calculateSkyBackground(operaWIRCamImage &image);
	
	/*! 
	 * void calculateQERatio(operaWIRCamImage *images[], unsigned count)
	 * \brief Calculate the QE Ratio - the median sky rate over the images
	 * Could also be gotten from the headers SKYLVL, SKYDEV -- if they exist
	 */
	void calculateQERatio(operaWIRCamImage *images[], unsigned count);
	
	/*! 
	 * void createSky(operaWIRCamImage *images[], unsigned count)
	 * \brief Create a master sky image
	 */
	void createSky(operaWIRCamImage *images[], unsigned count);
	
	/*! 
	 * void correctSkyLevel(operaWIRCamImage *images[], unsigned count)
	 * \brief Correct sky level
	 */ 
	void correctSkyLevel(operaWIRCamImage *images[], unsigned count);
	
	/*! 
	 * Box &getGuideWindow(unsigned index)
	 * \brief Retrieve the nth Guide Window Box, may be indexed by chip number.
	 */
	Box *getGuideWindow(unsigned index);
    /*
     * edgeExtend(Box &b, int Direction)
     * \brief extend the edge of a box in the 'X' (HORIZONTAL), 'Y' (VERTICAL), or 'N' (N = NONE, no extension) direction
     *  to the edge of the image
     * \param Box &b
     * \param int Direction
     */
    Box edgeExtend(Box &b, char Direction/*=X|Y|N*/);
	/*! 
	 * float getChipBias(void)
	 * \brief Retrieve the WIRCam ChipBias.
	 */
    float getChipBias(void);
	/*!
	* float getAverageSkyLevel(unsigned extension);
	* \brief Retrieve the average sky level for all extensions.
	*/
		float getAverageSkyLevel(void);		
	/*!
	 * float getSkyLevel(unsigned extension)
	 * \brief Retrieve the sky level for an extension.
	 */
		float getSkyLevel(unsigned extension);		
	/*!
	 * float getSkyRate(unsigned extension)
	 * \brief Retrieve the sky rate for an extension.
	 */
		float getSkyRate(unsigned extension);
};
#endif

