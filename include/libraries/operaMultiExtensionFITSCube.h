#ifndef OPERAMULTIEXTENSIONFITSCUBE_H
#define OPERAMULTIEXTENSIONFITSCUBE_H
/*******************************************************************
 ****               	OPERA PIPELINE v1.0                     ****
 *******************************************************************
 Class: operaMultiExtensionFITSCube  
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

#include "libraries/operaFITSImage.h"				// for operaFITSImage
#include "libraries/operaFITSCube.h"				// for operaFITSCube
#include "libraries/operaMultiExtensionFITSImage.h"	// for operaMultiExtensionFITSImage
#include "libraries/operaImageVector.h"				// for operaMultiExtensionFITSImage
#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"

/*
 * operaMultiExtensionFITSCube class
 * \author Doug Teeple
 * \brief This class encapsulates the MEF FITS image.
 * \file operaMultiExtensionFITSCube.h
 * \ingroup libraries
 */
class operaMultiExtensionFITSCube : public operaMultiExtensionFITSImage {
	
private:
    
protected:
	
public:
	/*!
	 * \class operaMultiExtensionFITSCube()
	 * \brief Basic operaMultiExtensionFITSCube constructor.
	 */
	operaMultiExtensionFITSCube();														// simply construct a default image
	 
	/*!
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
	operaMultiExtensionFITSCube(unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, unsigned Extensions, edatatype Datatype);	// simply construct an image of a given size
	/*! 
	 * \class operaMultiExtensionFITSCube
	 * \brief create a writeable file image of an in memory FITSImage object
	 * \brief operaMultiExtensionFITSCube(string Filename, unsigned Naxis1, unsigned Naxis2, edatatype Datatype, unsigned Compression)
	 * \brief Constructor to create a new FITS file.
	 * \param Filename to create (file is deleted if it exists)
	 * \param Naxis1 - dimensions
	 * \param Naxis2 - dimensions
	 * \param Datatype defaults to tshort
	 * \param Compression, defaults to none
	 */
	operaMultiExtensionFITSCube(string Filename, unsigned Naxis1, unsigned Naxis2, unsigned Naxis3, unsigned Extensions, edatatype Datatype, unsigned Compression = 0, bool isLazy = true);	// create a new empty FITSImage
	/*! 
	 * \class operaMultiExtensionFITSCube
	 * \brief create a FITSIMage object from a FITS file
	 * \brief operaMultiExtensionFITSCube(string Filename, int mode=READWRITE)
	 * \brief Constructor to create a FITSImage from a FITS file.
	 * \param Filename
	 * \param mode
	 * \return void 
	 */
	operaMultiExtensionFITSCube(string Filename, edatatype Datatype, int mode=READWRITE/*READONLY*/, unsigned Compression = 0, bool isLazy = true);		// read an existing FITSImage from file

    /*! 
     * operaMultiExtensionFITSCube(operaFITSImage &imageIn, bool ViewOnly = false, bool AddHeader = false)
     * \brief Clone a Multi Extension FITSImage object.
     * \param imageIn - pointer to image to clone
     */
	operaMultiExtensionFITSCube(operaMultiExtensionFITSCube &image, bool ViewOnly = false, bool AddHeader = false);
	
	~operaMultiExtensionFITSCube();			// destructor	

    /**********************************************************************************/
    /**********************************************************************************/
    /**** NOTE WELL:                                                               ****/
    /**** The operators are defined for thefloat datatype only.                    ****/
    /**********************************************************************************/
    /**********************************************************************************/
	
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return a operaFITSCube&.
	 */
	operaFITSCube& operator[](unsigned i) {
#ifdef FAST_EXTENSIONS
		setExtension((unsigned long)i);
		npixels = npixels_per_extension;
		return (operaFITSCube&)*this;
#else
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (i > this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
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
        setExtension(t->current_extension);
		t->super = this;							// record where we came from
		return (operaFITSCube&)*t;
#endif
    };
 	/* 
	 * \brief operator []
	 * \brief indexing operator to return a pixel value.
	 */
	operaImageVector *operator[](operaImageVector *vector) {
        SetImage(vector, this);
        return vector;
    };
    /*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b; copies the pixel values from b to a
	 */	
	operaMultiExtensionFITSCube& operator=(operaMultiExtensionFITSCube* b) {
        float *p = (float *)pixptr; 
		float *bp = (float *)b->pixptr;
		unsigned long n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b->extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b->naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
				// don't read the this, because we will write over it anyway...
				setExtension(ext);
				b->setExtension(ext);
				b->extensionHasBeenRead[ext] = false;
                b->readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b->pixptr;
                while (n--) *p++ = *bp++;
                saveExtension(ext);
                if (b->super)
					b->super->setHasBeenRead(ext);
                if (super)
					super->setHasBeenWritten(ext);
            }
		} else {
            while (n--) *p++ = *bp++;
			npixels = naxis1*naxis2*naxis3*(extensions==0?1:extensions);	// in case we were indexed...
        }
        if (b->istemp && !b->viewOnly) delete b;
		return *this;
	};
	/*!  
	 * \brief operator =
	 * \brief assignment operator.
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b; copies the pixel values from b to a
	 */	
	operaMultiExtensionFITSCube& operator=(operaMultiExtensionFITSCube& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
                // don't read the this, because we will write over it anyway...
				setExtension(ext);
				b.setExtension(ext);
				b.extensionHasBeenRead[ext] = false;
                b.readExtension(ext);
                n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
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
	 * \brief assignment operator.
	 * \note usage: operaMultiExtensionFITSCube a = 0.0; copies the float value to every pixel in = a
	 */	
	operaMultiExtensionFITSCube& operator=(float f) {
		float *p = (float *)pixptr;
		unsigned long n = npixels;
        if (isLazy) {
            // we need to read/write all extensions...
            for (unsigned ext=1; ext<=extensions; ext++) {
				setExtension(ext);
                // will be written over anyway... readExtension(ext);
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
	 * \brief Median operator %
	 * \note usage: operaMultiExtensionFITSCube& a % 8; median combines 8 slices of a
	 */	
	operaMultiExtensionFITSCube& operator%(unsigned int i) {
		this->medianCollapseParallel(extensions, i);
		return *this;
	};
    /*! 
	 * \brief operator >>
	 * \brief write image to named file (saveas).
	 * \note usage: operamultiextensionFITScube a(); a >> "foo.fits"; 
	 */	
	operaMultiExtensionFITSCube& operator>>(const string Filename) {
		operaMultiExtensionFITSCubeSaveAs(Filename);
		return *this;
	};
	/*! 
	 * \brief operator >>
	 * \brief write image to another image (headers and pixels).
	 * \note usage: operamultiextensionFITScube a(); operamultiextensionFITScube b(); a >> b; 
	 */	
	const operaMultiExtensionFITSCube& operator>>(operaMultiExtensionFITSCube& image) {
		image.operaMultiExtensionFITSCubeCopyHeader(this);
		image = this;
		return image;
	};
	/*! 
	 * \brief operator <<
	 * \brief read an image from a named file.
	 * \note usage: operamultiextensionFITScube a(); a << "foo.fits"; 
	 */	
	operaMultiExtensionFITSCube& operator<<(const string Filename) {
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
	 * \note usage: operamultiextensionFITScube a(); operamultiextensionFITScube b(); a << b; 
	 */	
	const operaMultiExtensionFITSCube& operator<<(operaMultiExtensionFITSCube& image) {
		operaMultiExtensionFITSCubeCopyHeader(&image);
		*this = image;
		return *this;
	};
    /*! 
	 * \brief operator ==
	 * \brief equality operator.
	 * \note usage: operamultiextensionFITScube a(); operamultiextensionFITScube b(); a = a == b;
     */
    operaMultiExtensionFITSCube& operator==(operaMultiExtensionFITSCube& b) {
		float *tp;
 		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this);
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
     * \note usage: operamultiextensionFITScube a(); operamultiextensionFITScube b(); a = a == b;
     */
    operaMultiExtensionFITSCube& operator==(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this);
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
	 * \note usage: operaMultiExtensionFITSCube a += operaMultiExtensionFITSCube b; adds the pixel values from b to a
	 */	
	operaMultiExtensionFITSCube& operator+=(operaMultiExtensionFITSCube& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
        if (isLazy) {
            for (unsigned ext=1; ext<=extensions; ext++) {
                readExtension(ext);
				n = npixels_per_extension; // has to be reset each time we go through the loop
				p = (float *)pixptr; 
				bp = (float *)b.pixptr;
				while (n--) *p++ += *bp++;
                saveExtension(ext);
            }
		} else {
            // we need to read/write all extensions...
            while (n--) *p++ += *bp++;
        }
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \note usage: operaMultiExtensionFITSCube a += 100.0; adds the float value to a
	 */	
	operaMultiExtensionFITSCube& operator+=(float f) {
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
	 * \note usage: operaMultiExtensionFITSCube a -= operaMultiExtensionFITSCube b; subtracts the pixel values from b to a
	 */	
	operaMultiExtensionFITSCube& operator-=(operaMultiExtensionFITSCube& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
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
	 * \note usage: operaMultiExtensionFITSCube a -= 100.0; subtracts the float value from a
	 */	
	operaMultiExtensionFITSCube& operator-=(float f) {
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
	 * \note usage: operaMultiExtensionFITSCube a *= operaMultiExtensionFITSCube b; multiplies the pixel values from b to a
	 */	
	operaMultiExtensionFITSCube& operator*=(operaMultiExtensionFITSCube& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels; 
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
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
	 * \note usage: operaMultiExtensionFITSCube a *= 100.0; multiplies the float value times a
	 */	
	operaMultiExtensionFITSCube& operator*=(float f) {
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
	 * \note usage: operaMultiExtensionFITSCube a /= operaMultiExtensionFITSCube b; divides the pixel values from b into a
	 */	
	operaMultiExtensionFITSCube& operator/=(operaMultiExtensionFITSCube& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
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
	 * \note usage: operaMultiExtensionFITSCube a /= 100.0; divides the float value from a
	 */	
	operaMultiExtensionFITSCube& operator/=(float f) {
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b * operaMultiExtensionFITSCube c; multiplies the pixel values  b * c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator*(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b * 10.0; multiplies the pixel values  b * 10.0 and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator*(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, this->viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b / operaMultiExtensionFITSCube c; divides the pixel values  b / c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator/(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b / 100.0; divides the pixel values  b / 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator/(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, this->viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b + operaMultiExtensionFITSCube c; adds the pixel values  b + c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator+(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b + 100.0; adds the pixel values  b + 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator+(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, this->viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b - operaMultiExtensionFITSCube c; subtracts the pixel values  b - c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator-(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b - 100.0; subtracts the pixel values  b - 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator-(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, this->viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = !operaMultiExtensionFITSCube b; inverts the pixel values of b and assigns to a (0.0 becomes 1.0, non-zero becomes 0.0)
	 */	
	operaMultiExtensionFITSCube& operator!() {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b > operaMultiExtensionFITSCube c; creates a mask of 1.0 if b > c or 0.0 if b <= c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator>(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b > 100.0; creates a mask of 1.0 if b > 100.0 or 0.0 if b <= 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator>(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, this->viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b >= operaMultiExtensionFITSCube c; creates a mask of 1.0 if b >= c or 0.0 if b < c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator>=(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b >= 100.0; creates a mask of 1.0 if b >= 100.0 or 0.0 if b < 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator>=(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, this->viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b < operaMultiExtensionFITSCube c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator<(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b < operaMultiExtensionFITSCube c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator<(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, this->viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b <= operaMultiExtensionFITSCube c; creates a mask of 1.0 if b <= c or 0.0 if b > c and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator<=(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b <= 100.0; creates a mask of 1.0 if b <= 100.0 or 0.0 if b > 100.0 and assigns to a
	 */	
	operaMultiExtensionFITSCube& operator<=(float f) {
		operaMultiExtensionFITSCube *t = new operaMultiExtensionFITSCube(*this, this->viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b && operaMultiExtensionFITSCube c; creates a mask of 1.0 where b && c != 0.0
	 */	
	operaMultiExtensionFITSCube& operator&&(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * \note usage: operaMultiExtensionFITSCube a = operaMultiExtensionFITSCube b && operaMultiExtensionFITSCube c; creates a mask of 1.0 where b && c != 0.0
	 */	
	operaMultiExtensionFITSCube& operator||(operaMultiExtensionFITSCube& b) {
		float *tp;
		operaMultiExtensionFITSCube *t = NULL;
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if (b.extensions != this->extensions) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorExtensionOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (b.naxis3 != this->naxis3) {
			throw operaException("operaMultiExtensionFITSCube: ", operaErrorSliceOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaMultiExtensionFITSCube(*this, viewOnly);
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
	 * operaMultiExtensionFITSCubeCopyHeader(operaMultiExtensionFITSCube *from) 
	 * \brief Copies all of the header information from image. 
	 */
	void operaMultiExtensionFITSCubeCopyHeader(operaMultiExtensionFITSCube *from);
	/*! 
	 * operaMultiExtensionFITSCubeSave() 
	 * \brief Saves the current image to disk.
	 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
	 */
	void operaMultiExtensionFITSCubeSave();
	
	/*! 
	 * operaMultiExtensionFITSCubeSaveAs(string newFilename) 
	 * \brief Saves the current image to disk, with the given filename.
	 */
	void operaMultiExtensionFITSCubeSaveAs(string newFilename);
    /*! 
     * void getSlice()
     * \brief get the current slice
     */
    unsigned getSlice();
	/*! 
	 * unsigned getslices(void)
	 * \brief gets the number of slices
	 */
	unsigned getslices(void);
    /*! 
	 * void setSlice(int slice)
	 * \brief sets the current slice
	 */
	void setSlice(unsigned slice);
	/*! 
	 * void readSlice(int slice)
	 * \brief read a slice into memory
	 */
	void readSlice(unsigned slice);
	/*! 
	 * void saveSlice(int slice)
	 * \brief save slice data
	 */
	void saveSlice(unsigned slice);
    /*!
     * void medianCollapse()
     * \brief collapse the cube into a single image
     * \details Collapses entire image if AllExtensions = true
     * Collapses current_extension if AllExtensions = false
     * To collapse a single extension, use medianCollapseSingleExtension(extension)
     * WARNING: that this method uses operaArrayMedianQuick, so the input image will be scrambled.
     */
    void medianCollapse();
    /*!
     * void medianCollapse(unsigned int n)
     * \brief collapse n slices of the cube into a single image
     * \details Collapses entire image if AllExtensions = true
     * Collapses current_extension if AllExtensions = false
     * To collapse a single extension, use medianCollapseSingleExtension(extension)
     * WARNING: that this method uses operaArrayMedianQuick, so the input image will be scrambled.
     */
    void medianCollapse(unsigned int n);
	/*!
	 * void medianCollapseParallel(unsigned int maxthreads)
	 * \brief extension-based parallel collapse the cube into a single image
	 * Collapses entire image if AllExtensions = true
	 * Collapses current_extension if AllExtensions = false
	 * To collapse a single extension, use medianCollapseSingleExtension(extension)
	 * Note: indexing for a MEF Cube is [extension][slice][y][x]
	 */
	void medianCollapseParallel(unsigned int maxthreads);
	/*!
	 * void medianCollapseParallel(unsigned int maxthreads, unsigned n slices)
	 * \brief extension-based parallel collapse n slices of the cube into a single image
	 * Collapses entire image if AllExtensions = true
	 * Collapses current_extension if AllExtensions = false
	 * To collapse a single extension, use medianCollapseSingleExtension(extension)
	 * Note: indexing for a MEF Cube is [extension][slice][y][x]
	 */
	void medianCollapseParallel(unsigned int maxthreads, unsigned int n);
    /*!
     * void medianCollapseSingleExtension(int extension)
     * \brief collapse the cube into a single image
     * \details Collapses a single extension
     * To collapse all extensions, use medianCollapse in non-lazy read (AllExtensions = true)
     * WARNING: that this method uses operaArrayMedianQuick, so the input image will be scrambled.
     */
    void medianCollapseSingleExtension(int extension);
    /*! 
     * void setExtension(int extension)
     * \brief sets the new_extension to the current extension
     */
    void setExtension(unsigned new_extension);	
    /*! 
     * void setExtension(int extension, int slice)
     * \brief sets the new_extension to the current extension
     */
    void setExtension(unsigned new_extension, unsigned slice);
    /*! 
	 * void readExtension(int extension)
	 * \brief read an extension into memory
	 */
	void readExtension(unsigned extension);
	/*! 
	 * void readExtension(int extension, int slice)
	 * \brief read an extension into memory
	 */
	void readExtension(unsigned extension, unsigned slice);
    /*! 
	 * void saveExtension(int extension)
	 * \brief save image data
	 */
	void saveExtension(unsigned extension);
    /*! 
	 * void saveExtension(int extension, int slice)
	 * \brief save image data
	 */
	void saveExtension(unsigned extension, unsigned slice);
    /*!
     * bool isMEFCube()
     * \brief Is this image a MEFcube?
     */    
    bool isMEFCube(void);

};

/*! 
 * void map(operaMultiExtensionFITSCube &, float from, float to)
 * \brief map all pixels from into to.
 * usage: map(image, F_NAN, 0.0);
 */
static inline void map(operaMultiExtensionFITSCube &image, float from, float to) {
	float *p = (float *)image.getpixels(); 
	unsigned n = image.getnpixels(); 
	while (n--) {*p=(*p==from?to:*p); p++;}
}
#endif

