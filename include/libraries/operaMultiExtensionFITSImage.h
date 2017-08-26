#ifndef OPERAMULTIEXTENSIONFITSIMAGE_H
#define OPERAMULTIEXTENSIONFITSIMAGE_H
/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
Class: operaMultiExtensionFITSImage 
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

#include "libraries/operaFITSImage.h"			// for operaFITSImage
#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"

/*!
 * \file operaMultiExtensionFITSImage.h
 */


// These definitions are for the eventual paralllelization of extension operations
typedef enum {
	op_mul,op_div,op_add,op_sub
} ops_t;

class operaMultiExtensionFITSCube;

/*!
 * \author Megan Tannock
 * \brief This class encapsulates multi extension FITS images.
 * \ingroup libraries
 * \sa class operaFITSImage, class operaMultiExtensionFITSCube, class operaFITSCube
 * \details This class is designed to process multi extension FITS images. The operators are defined
 * to work with an entire image, and goes pixel by pixel. This simplifies calling the 
 * operators for the user, and they do not need to use any loops. The methods included are
 * specific to processing a multi extension FITS Image.
 */
class operaMultiExtensionFITSImage : public operaFITSImage  {

	friend class operaMultiExtensionFITSCube;
	
private:	
	
protected:	
	
public:
    
    /*
	 * Constructors / Destructors
	 */
    
	/*!
	 * \sa class operaMultiExtensionFITSImage()
	 * \brief Basic operaMultiExtensionFITSImage constructor.
	 */
	operaMultiExtensionFITSImage();     // simply construct a default image    
	/*!
	 * \sa class operaMultiExtensionFITSImage
	 * \brief construct an in-memory FITSImage object
	 * \details operaMultiExtensionFITSImage(unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype=tushort)
	 * \details Create an in-memory image of given dimensions.
     * \throws operaException operaErrorNoMemory
	 */
	operaMultiExtensionFITSImage(unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype = tushort);	// simply construct an image of a given size
	/*! 
	 * \sa class operaMultiExtensionFITSImage
	 * \brief create a writeable file image of an in memory FITSImage object
	 * \details operaMultiExtensionFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype, unsigned Compression, bool isLazy)
	 * \details Constructor to create a new FITS file.
	 * \throws operaException cfitsio error code
     * \throws operaException operaErrorNoMemory
	 */
	operaMultiExtensionFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Extensions, edatatype Datatype, unsigned Compression = 0, bool isLazy = true);	// create a new empty FITSImage
	/*! 
     * \sa class operaMultiExtensionFITSImage
     * \brief create a multi extension FITS image object from a multi extension FITS file
     * \details operaMultiExtensionFITSImage(string Filename, edatatype Datatype, int mode, unsigned Compression, bool isLazy)
     * \details Constructor to create a MEF FITS Image from a MEF FITS file.
     * \return void
	 */ 
	operaMultiExtensionFITSImage(string Filename, edatatype Datatype, int mode=READWRITE/*READONLY*/, unsigned Compression = 0, bool isLazy = true);		// read an existing FITSImage from file
    /*! 
     * \sa class operaMultiExtensionFITSImage
     * \brief Clone a Multi Extension FITSImage object.
     * \details operaMultiExtensionFITSImage(operaFITSImage &imageIn, bool ViewOnly, bool AddHeader)
     */
	operaMultiExtensionFITSImage(operaMultiExtensionFITSImage &image, bool ViewOnly = false, bool AddHeader = false);
	
	/* 
	 * operaMultiExtensionFITSImage(operaMultiExtensionFITSCube &imageIn, unsigned slice)
	 * \brief Clone a Multi Extension FITSImage object.
	 */
	operaMultiExtensionFITSImage(operaMultiExtensionFITSCube &imageIn, unsigned slice, bool AddHeader = false);
	
    ~operaMultiExtensionFITSImage();		// destructor	

    
    /**********************************************************************************/
    /**********************************************************************************/
    /**** NOTE WELL:                                                               ****/
    /**** The operators are defined for the float datatype only.                   ****/
    /**********************************************************************************/
    /**********************************************************************************/
	
	/* OPERATORS */
    
	/*! 
	 * \brief Indexing operator, returns an operaFITSImage&
     * \details Indexes into a single extension in a multi extension FITS image
     * \throws operaException operaErrorLengthMismatch
	 */
	operaFITSImage& operator[](unsigned i) {
#ifdef FAST_EXTENSIONS
		if (current_extension != i)
			setExtension((unsigned long)i);
		npixels = npixels_per_extension;
		return (operaFITSImage&)*this;
#else
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (i > this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this, true);
		setExtension((unsigned long)i);
		t->npixels = npixels_per_extension;
		t->naxis = 2;
		t->naxis3 = 1;
		t->naxes[0] = t->naxis1;
		t->naxes[1] = t->naxis2;
		t->naxes[2] = t->naxis3;
		t->current_extension = (unsigned long)i;
		t->current_slice = 1;
		t->extensions = extensions;
		t->super = this;							// record where we came from
		return *t;			
#endif
	};
    /*! 
	 * \brief Indexing operator [], returns operaImageVector pointer
	 * \brief indexing operator to return a pixel value.
	 * \note usage:  
	 */
	operaImageVector *operator[](operaImageVector *vector) {
        SetImage(vector, this);
        return vector;
    };
    /*! 
	 * \brief Assignment operator =
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b; copies the pixel values from b to a
	 */	
	operaMultiExtensionFITSImage& operator=(operaMultiExtensionFITSImage* b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b->pixptr;
		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b->extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b->naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
			// we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
				setExtension(ext);
				b->extensionHasBeenRead[ext] = false;
				b->setExtension(ext);
				b->readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b->pixptr;
				// don't read the this, because we will write over it anyway...
                while (n--) *p++ = *bp++;
                saveExtension(ext);
                if (b->super)
					b->super->setHasBeenRead(ext);
                if (super)
					super->setHasBeenWritten(ext);
            }
		} else {
			// we need to copy all extensions that are stacked in memory
            while (n--) *p++ = *bp++;
			npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);	// in case we were indexed...
        }
		if (b->istemp && !b->viewOnly) delete b;
        return *this;
	};
	/*! 
	 * \brief Assignment operator =
	 * \brief assignment operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b; copies the pixel values from b to a
	 */	
	operaMultiExtensionFITSImage& operator=(operaMultiExtensionFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr;
		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
				setExtension(ext);
				b.extensionHasBeenRead[ext] = false;
				b.setExtension(ext);
                b.readExtension(ext);
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
                n = npixels_per_extension; // has to be reset each time we go through the loop
                // don't read the this, because we will write over it anyway...
                while (n--) *p++ = *bp++;
                saveExtension(ext);
                if (b.super)
					b.super->setHasBeenRead(ext);
                if (super)
					super->setHasBeenWritten(ext);
            }
		} else {
            while (n--) *p++ = *bp++;
			npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);	// in case we were indexed...
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief Assignment operator =
	 * \note usage: operaMultiExtensionFITSImage a = 0.0; copies the float value to every pixel in = a
	 */	
	operaMultiExtensionFITSImage& operator=(float f) {
		float *p = (float *)pixptr;
		unsigned long n = npixels;
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                // don' do this, it will be written over anyway... readExtension(ext);
				setExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
                while (n--) *p++ = f;
                saveExtension(ext);
                if (super)
					super->setHasBeenWritten(ext);
            }
		} else {
            while (n--) *p++ = f;
        }
		return *this;
	};
    /*! 
	 * \brief operator >>
	 * \brief write image to named file (saveas).
	 * \note usage: operamultiextensionFITSimage a(); a >> "foo.fits"; 
	 */	
	operaMultiExtensionFITSImage& operator>>(const string Filename) {
		operaMultiExtensionFITSImageSaveAs(Filename);
		return *this;
	};
	/*! 
	 * \brief operator >>
	 * \brief write image to another image (headers and pixels).
	 * \note usage: operamultiextensionFITSimage a(); operamultiextensionFITSimage b(); a >> b; 
	 */	
	const operaMultiExtensionFITSImage& operator>>(operaMultiExtensionFITSImage& image) {
		image.operaMultiExtensionFITSImageCopyHeader(this);
		image = this;
		return image;
	};
	/*! 
	 * \brief operator <<
	 * \brief read an image from a named file.
	 * \note usage: operamultiextensionFITSimage a(); a << "foo.fits"; 
	 */	
	operaMultiExtensionFITSImage& operator<<(const string Filename) {
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
	 * \note usage: operamultiextensionFITSimage a(); operamultiextensionFITSimage b(); a << b; 
	 */	
	const operaMultiExtensionFITSImage& operator<<(operaMultiExtensionFITSImage& image) {
		operaMultiExtensionFITSImageCopyHeader(&image);
		*this = image;
		return *this;
	};
    /*! 
	 * \brief operator ==
     * \brief equality operator, returns 1.0 for each pixel equal in the other image, else 0.0.
	 * \note usage: operamultiextensionFITSimage a(); operamultiextensionFITSimage b(); a = a == b; 
     */
    operaMultiExtensionFITSImage& operator==(operaMultiExtensionFITSImage& b) {
		float *tp;
 		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
     * \brief equality operator, returns 1.0 for each pixel equal in the other image, else 0.0.
	 * \note usage: operamultiextensionFITSimage a(); operamultiextensionFITSimage b(); a = a == b; 
     */
    operaMultiExtensionFITSImage& operator==(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a += operaMultiExtensionFITSImage b; adds the pixel values from b to a
	 */	
	operaMultiExtensionFITSImage& operator+=(operaMultiExtensionFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
 				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
                while (n--) *p++ += *bp++;
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ += *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \note usage: operaMultiExtensionFITSImage a += 100.0; adds the float value to a
	 */	
	operaMultiExtensionFITSImage& operator+=(float f) {
		float *p = (float *)pixptr; 
		unsigned long n = npixels;
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
 				p = (float *)pixptr; 
                while (n--) *p++ += f;
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ += f;
        }
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \note usage: operaMultiExtensionFITSImage a -= operaMultiExtensionFITSImage b; subtracts the pixel values from b to a
	 */	
	operaMultiExtensionFITSImage& operator-=(operaMultiExtensionFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
 				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
                while (n--) *p++ -= *bp++;
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ -= *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \note usage: operaMultiExtensionFITSImage a -= 100.0; subtracts the float value from a
	 */	
	operaMultiExtensionFITSImage& operator-=(float f) {
		float *p = (float *)pixptr; 
		unsigned long n = npixels;
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
 				p = (float *)pixptr; 
                while (n--) *p++ -= f;
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ -= f;
        }
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \note usage: operaMultiExtensionFITSImage a *= operaMultiExtensionFITSImage b; multiplies the pixel values from b to a
	 */	
	operaMultiExtensionFITSImage& operator*=(operaMultiExtensionFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
 				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
                while (n--) *p++ *= *bp++;
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ *= *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \note usage: operaMultiExtensionFITSImage a *= 100.0; multiplies the float value times a
	 */	
	operaMultiExtensionFITSImage& operator*=(float f) {
		float *p = (float *)pixptr; 
 		unsigned long n = npixels;
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
 				p = (float *)pixptr; 
                while (n--) *p++ *= f;
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ *= f;
        }
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \note usage: operaMultiExtensionFITSImage a /= operaMultiExtensionFITSImage b; divides the pixel values from b into a
	 */	
	operaMultiExtensionFITSImage& operator/=(operaMultiExtensionFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
 		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
 				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
                while (n--) *p++ /= *bp++;
                saveExtension(ext);
            }
        } else {
            while (n--) *p++ /= *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \note usage: operaMultiExtensionFITSImage a /= 100.0; divides the float value from a
	 */	
	operaMultiExtensionFITSImage& operator/=(float f) {
		float *p = (float *)pixptr; 
 		unsigned long n = npixels; 
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
 				p = (float *)pixptr; 
                while (n--) *p++ /= f;
                saveExtension(ext);
            }
		} else {
            while (n--) *p++ /= f;
        }
		return *this;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b * operaMultiExtensionFITSImage c; multiplies the pixel values  b * c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator*(operaMultiExtensionFITSImage& b) {
		float *tp;
 		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
            while (n--) *tp++ = *p++ * *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b * 10.0; multiplies the pixel values  b * 10.0 and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator*(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
               while (n--) *tp++ = *p++ * f;
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ * f;
        }
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b / operaMultiExtensionFITSImage c; divides the pixel values  b / c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator/(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
                while (n--) *tp++ = *p++ / *bp++;
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ / *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b / 100.0; divides the pixel values  b / 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator/(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
                while (n--) *tp++ = *p++ / f;
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ / f;
        }
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b + operaMultiExtensionFITSImage c; adds the pixel values  b + c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator+(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
                while (n--) *tp++ = *p++ + *bp++;
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ + *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b + 100.0; adds the pixel values  b + 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator+(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
                while (n--) *tp++ = *p++ + f;
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ + f;
        }
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b - operaMultiExtensionFITSImage c; subtracts the pixel values  b - c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator-(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
                while (n--) *tp++ = *p++ - *bp++;
                saveExtension(ext);
            }
        } else {
            while (n--) *tp++ = *p++ - *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b - 100.0; subtracts the pixel values  b - 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator-(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
                while (n--) *tp++ = *p++ - f;
                saveExtension(ext);
            }
		} else {
            while (n--) *tp++ = *p++ - f;
        }
		return *t;
	};
	/*! 
	 * \brief operator !
	 * \brief invert operator.
	 * \note usage: operaMultiExtensionFITSImage a = !operaMultiExtensionFITSImage b; inverts the pixel values of b and assigns to a (0.0 becomes 1.0, non-zero becomes 0.0)
	 */	
	operaMultiExtensionFITSImage& operator!() {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b > operaMultiExtensionFITSImage c; creates a mask of 1.0 if b > c or 0.0 if b <= c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator>(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b > 100.0; creates a mask of 1.0 if b > 100.0 or 0.0 if b <= 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator>(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b >= operaMultiExtensionFITSImage c; creates a mask of 1.0 if b >= c or 0.0 if b < c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator>=(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b >= 100.0; creates a mask of 1.0 if b >= 100.0 or 0.0 if b < 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator>=(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b < operaMultiExtensionFITSImage c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator<(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b < operaMultiExtensionFITSImage c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator<(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b <= operaMultiExtensionFITSImage c; creates a mask of 1.0 if b <= c or 0.0 if b > c and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator<=(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b <= 100.0; creates a mask of 1.0 if b <= 100.0 or 0.0 if b > 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSImage& operator<=(float f) {
		operaMultiExtensionFITSImage *t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b && operaMultiExtensionFITSImage c; creates a mask of 1.0 where b && c != 0.0
	 */	
	operaMultiExtensionFITSImage& operator&&(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
	 * \note usage: operaMultiExtensionFITSImage a = operaMultiExtensionFITSImage b && operaMultiExtensionFITSImage c; creates a mask of 1.0 where b && c != 0.0
	 */	
	operaMultiExtensionFITSImage& operator||(operaMultiExtensionFITSImage& b) {
		float *tp;
		operaMultiExtensionFITSImage *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSImage: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSImage(*this);
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
	 * bool isInstantiated(operaMultiExtensionFITSImage &Image)
	 * \brief Is this image in memory yet?
	 * \param operaMultiExtensionFITSImage &Image
	 * \return bool
	 */
	bool isInstantiated(operaMultiExtensionFITSImage &Image);
	/*! 
	 * operaMultiExtensionFITSImageSave() 
	 * \brief Saves the current image to disk.
	 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
	 * \throws operaException operaErrorCodeNoFilename
	 * \throws operaException cfitsio error code
	 */
	void operaMultiExtensionFITSImageSave();
	/*! 
	 * operaMultiExtensionFITSImageSaveAs(string newFilename) 
	 * \brief Saves the current image to disk, with the given filename.
	 * \throws operaException cfitsio error code
	 */
	void operaMultiExtensionFITSImageSaveAs(string newFilename);
	
	/*! 
	 * void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageConvertImage(edatatype todatattype)
	 * \brief Convert image type and create new image values of that type.
	 * \throws operaException cfitsio error code
	 */
	void operaMultiExtensionFITSImageConvertImage(edatatype todatatype);
    /*! 
     * void getExtension()
     * \brief get the current extension number
     */
    int getExtension();
    /*! 
	 * void addExtension(unsigned Columns, unsigned Rows)
	 * \briefadds an extension 
	 */
	void addExtension(string name, unsigned Columns, unsigned Rows, DATASEC_t &datasec, bool reuse = false);
    /*! 
     * void saveExtension()
     * \brief save image data
     */
    void saveExtension(unsigned extension);
    /*! 
     * void readExtension(unsigned extension)
     * \brief read an extension into memory
     */
    void readExtension(unsigned extension);
	/*! 
	 * void setExtension(unsigned extension)
	 * \brief sets the current extension
	 * \param unsigned new_extension
	 */
	void setExtension(unsigned new_extension);
    /*! 
     * operaFITSImageCopyHeader(operaFITSImage *from) 
     * \brief Copies all of the header information from image.
     * \throws operaException cfitsio error code
     * \return void
     */
    void operaMultiExtensionFITSImageCopyHeader(operaMultiExtensionFITSImage *from);    //copy header unit
	/* 
	 * operaMultiExtensionFITSImageCopyNonStructuralHeader(operaMultiExtensionFITSImage *from, unsigned extension) 
	 * \brief Copies all of the non-structural header information from image. 
	 * \throws operaException cfitsio error code
	 */
	void operaMultiExtensionFITSImageCopyNonStructuralHeader(operaMultiExtensionFITSImage *from);
	/*! 
	 * string operaFITSGetHeaderValue(string keyword) 
	 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
	 * \throws operaException cfitsio error code
	 */
	string operaFITSGetHeaderValue(string keyword);							// return header kweyword value as string
	/*! 
	 * float operaFITSGetFloatHeaderValue(string keyword) 
	 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
	 * \throws operaException cfitsio error code
	 */
	float operaFITSGetFloatHeaderValue(string keyword);
	/*! 
	 * int operaFITSGetFloatHeaderValue(string keyword) 
	 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
	 * \throws operaException cfitsio error code
	 */
	int operaFITSGetIntHeaderValue(string keyword);
	/*! 
	 * string operaFITSGetRawHeaderValue(string keyword) 
	 * \brief returns the value of a FITS keyword verbatim.
	 * \throws operaException cfitsio error code
	 */
	string operaFITSGetRawHeaderValue(string keyword);						// return raw header keyword value as string
	/*! 
     * string operaFITSGetHeaderValue(string keyword, unsigned extension) 
     * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
	 * \throws operaException cfitsio error code
     */
    string operaFITSGetHeaderValue(string keyword, unsigned extension);
    /*! 
     * float operaFITSGetFloatHeaderValue(string keyword, unsigned extension) 
     * \brief returns the float value of a FITS keyword.
	 * \throws operaException cfitsio error code
     */
    float operaFITSGetFloatHeaderValue(string keyword, unsigned extension);
    /*! 
     * int operaFITSGetIntHeaderValue(string keyword, unsigned extension) 
     * \brief returns the int value of a FITS keyword.
	 * \throws operaException cfitsio error code
     */
    int operaFITSGetIntHeaderValue(string keyword, unsigned extension);
    /*! 
     * string operaFITSGetRawHeaderValue(string keyword, unsigned extension) 
     * \brief returns the value of a FITS keyword verbatim.
     */
    string operaFITSGetRawHeaderValue(string keyword, unsigned extension);

	/*! 
	 * operaFITSSetHeaderValue(string keyword, string value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, string value, string comment, unsigned extension); // set header kweyword value as string
	/*! 
	 * operaFITSSetHeaderValue(string keyword, short value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, short value, string comment, unsigned extension); // set header keyword value as ushort
	/*! 
	 * operaFITSSetHeaderValue(string keyword, unsigned short value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, unsigned short value, string comment, unsigned extension); // set header keyword value as ushort
	/*! 
	 * operaFITSSetHeaderValue(string keyword, float value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, float value, string comment, unsigned extension); // set header keyword value as float
	/*! 
	 * operaFITSSetHeaderValue(string keyword, double value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, double value, string comment, unsigned extension); // set header keyword value as double
	/*! 
	 * operaFITSAddComment(string comment)
	 * \brief adds a comment to header.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSAddComment(string comment, unsigned extension);
	/*! 
	 * operaFITSDeleteHeaderKey(string keyword)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSDeleteHeaderKey(string keyword, unsigned extension);
	/*! 
	 * operaFITSSetHeaderValue(string keyword, string value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, string value, string comment); // set header kweyword value as string
	/*! 
	 * operaFITSSetHeaderValue(string keyword, unsigned short value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, unsigned short value, string comment); // set header keyword value as ushort
	/*! 
	 * operaFITSSetHeaderValue(string keyword, float value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, float value, string comment); // set header keyword value as float
	/*! 
	 * operaFITSSetHeaderValue(string keyword, double value, string comment)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSSetHeaderValue(string keyword, double value, string comment); // set header keyword value as double
	/*! 
	 * operaFITSAddComment(string comment)
	 * \brief adds a comment to header.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSAddComment(string comment);
	/*! 
	 * operaFITSDeleteHeaderKey(string keyword)
	 * \brief sets the given keyword to value with comment.
	 * \throws operaException cfitsio error code
	 */
	void operaFITSDeleteHeaderKey(string keyword);
    /*!
     * bool isMEF()
     * \brief Is this image a MEF?
     */    
    bool isMEF(void); 
    /* 
     * void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageWriteVirtual(operaFITSSubImage &subImage, unsigned X, unsigned Y, unsigned extension)
     * \brief Write a subimage on to avirtual operaMultiExtensionFITSImage.
     * \details This is used for deep stacking. 
     * \details The virtual image is Lazy, and can be very large.
     * \details The actual image exists only on disk.
     */
    void operaMultiExtensionFITSImageWriteVirtual(operaFITSSubImage &subImage, unsigned long X, unsigned long Y, unsigned extension);    
    /* 
     * void operaMultiExtensionFITSImage::operaMultiExtensionFITSImageReadVirtual(operaFITSSubImage &subImage, unsigned X, unsigned Y)
     * \brief Read a subimage from a virtual operaFITSImage.
     * \details This is used for deep stacking. 
     * \details The virtual image is Lazy, and can be very large.
     * \details The actual image exists only on disk.
     */
    void operaMultiExtensionFITSImageReadVirtual(operaFITSSubImage &subImage, unsigned long X, unsigned long Y, unsigned extension);    
    
    
	operaMultiExtensionFITSImage &parmul(operaMultiExtensionFITSImage &b);

};

/*! 
 * void map(operaMultiExtensionFITSImage &, float from, float to)
 * \brief map all pixels from into to.
 * usage: map(image, F_NAN, 0.0);
 */
static inline void map(operaMultiExtensionFITSImage &image, float from, float to) {
	float *p = (float *)image.getpixels(); 
	unsigned n = image.getnpixels(); 
	while (n--) {*p=(*p==from?to:*p); p++;}
}
#endif

