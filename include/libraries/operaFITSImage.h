#ifndef OPERAFITSIMAGE_H
#define OPERAFITSIMAGE_H
/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaGITSImage   
 Class: operaFITSImage 
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
#include <fitsio.h>								// usually in /usr/local/include/

/*! 
 * operaFITSImage
 * \author Doug Teeple
 * \brief This class encapsulates the FITS image.
 * \file operaFITSImage.h
 * \ingroup libraries
 */

class operaFITSSubImage;
class Box;
class operaImageVector;

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaLib.h"					// for trimFITSKeyword

using namespace std;

#define MAXFITSDIMENSIONS 3
#define MAXFITSEXTENSIONS 20
#define cfitsioextension(e) (e+1)

typedef struct DATASEC {
	unsigned x1;
	unsigned x2;
	unsigned y1;
	unsigned y2;
} DATASEC_t;

const int primary_extension = 0;
const int first_extension = primary_extension+1;
const int first_slice = 1;
const bool lazyRead = true;

//					8				16						20						-32					-64
enum ebitpix {byte_img=BYTE_IMG,short_img=SHORT_IMG, ushort_img=USHORT_IMG, float_img=FLOAT_IMG, double_img=DOUBLE_IMG};

//					8				21				20									42			82
enum edatatype {tbyte=TBYTE, tshort=TSHORT, tushort=TUSHORT, tstring=TSTRING, tfloat=TFLOAT, tdouble=TDOUBLE};

//								21				11				41					31
enum eCompression {cNone=0, cGZIP=GZIP_1, cRICE=RICE_1, cHCOMPRESS=HCOMPRESS_1, cPLIO=PLIO_1};

enum eImageType {UNK, MEF, MEFCube, FITSCube, FITS};
class operaImageVector;
class operaFITSImage;

void SetImage(operaImageVector *vector, operaFITSImage *image);

/*! 
 * void getFITSImageInformation(string Filename, int *XDimension, int *YDimension, int *ZDimension, int *Extensions, edatatype *Datatype, long *Npixels)
 * \brief get useful information about a FITS file.
 */
void getFITSImageInformation(string Filename, unsigned *XDimension, unsigned *YDimension, unsigned *ZDimension, unsigned *Extensions, edatatype *Datatype, long *Npixels);

/*!
 * operaFITSIMage class
 * \author Doug Teeple
 * \brief This class encapsulates the FITS image.
 * \file operaFITSIMage.h
 * \package operaFITSIMage
 * \ingroup libraries
 */
class operaFITSImage {
    
	friend class operaEspadonsImage;
	friend class operaFITSSubImage;
	friend class operaFITSCube;
	friend class operaMultiExtensionFITSCube;
	friend class operaMultiExtensionFITSImage;
	friend class operaWIRCamImage;
    
private:
	void openFITSfile(string Filename, int Mode=READWRITE/*READONLY*/);
	void readFITSHeaderInfo();
	void readFITSArray();
    
protected:	
	string filename;					// filename
	fitsfile *fptr;						// FITS file pointer
	ebitpix bitpix;						// BITPIX keyword value (SHORT_IMG, USHORT_IMG, FLOAT_IMG, DOUBLE_IMG)
	float bzero;						// if bzero == 32768. then this is tushort else tshort
	float bscale;						// if bscale == 1. then this is tushort else tshort
	unsigned hdu;						// current active extension hdu
	unsigned nhdus;						// number of hdus
	long naxes[MAXFITSDIMENSIONS];		// FITS image dimensions
	long naxis;							// FITS image dimensionality
	unsigned long naxis1;				// x-dimension to be figured out from NAXIS1 (ncols) 
	unsigned long naxis2;				// y-dimension to be figured out from NAXIS2 (nrows)
	unsigned long naxis3;				// z-dimension to be figured out from NAXIS3 (nslices)
	unsigned long npixels;				// total number of pixels in all slices and extensions
	unsigned long npixels_per_slice;	// number of ccd pixels per slice	
	unsigned long npixels_per_extension;// number of ccd pixels per extension	
	unsigned compression;				// NULL, GZIP_1, RICE_1, HCOMPRESS_1 or PLIO_1
	edatatype datatype;					// (TSHORT, TUSHORT, TFLOAT, TDOUBLE)		
	bool istemp;						// set if this instance is a temp created in an expression
	int mode;							// READWRITE or READONLY - retained for efficiency iin setExetnsion
	void *pixptr;						// pixel base
	void *varptr;						// variance base for propagating errors
    
	ebitpix tobitpix(edatatype Datatype);		// get BITPIX from datatype
	edatatype todatatype(ebitpix Bitpix, float bzero, float bscale);	// get datatype from BITPIX, BZERO, SCALE
	size_t toSize(ebitpix bitpix, long npixels);	// get size in bytes from BITPIX and npixels
	unsigned current_extension;			// current active extension
	unsigned current_slice;				// current active slice in cube
	unsigned extensions;				// number of extensions
    bool isLazy;						// is this a lazy read?
    bool viewOnly;						// is this a view of somebody else's pixels? If so be careful about deletion!
    bool isClone;						// is this a clone of somebody else's fptr? If so do not close!
    bool AllExtensions;					// are all extensions in memory?
    bool AllSlices;						// are all slices in memory?
	eImageType imageType;				// is this MEF, MEFCube, Cube, FITS?
	bool extensionHasBeenRead[MAXFITSEXTENSIONS+1];		// in the case of a lazy extension by extension read, have we read in that extension yet?
	bool extensionHasBeenWritten[MAXFITSEXTENSIONS+1];	// in the case of a lazy extension by extension read, have we saved that extension yet?
    operaFITSImage *super;				// the super class that created this instance
	
public:
	/*!
	 * \sa class operaFITSImage()
	 * \brief Basic operaFITSImage constructor.
	 */
	operaFITSImage(bool isLazy = false);	// simply construct a default image
	/*!
	 * \sa class operaFITSImage(string Filename)
	 * \brief Basic operaFITSImage constructor from a FITS file.
	 */
	operaFITSImage(string Filename, bool isLazy = false);
	
	/*!
	 * \sa class operaFITSImage
	 * \brief construct an in-memory FITSImage object
	 * \brief operaFITSImage(unsigned Naxis1, unsigned Naxis2, edatatype Datatype=tushort, unsigned Compression=0)
	 * \brief Create an in-memory image of given dimensions.
	 */
	operaFITSImage(unsigned Naxis1, unsigned Naxis2, edatatype Datatype = tushort);	// simply construct an image of a given size
	/*! 
	 * \sa class operaFITSImage
	 * \brief create a writeable file image of an in memory FITSImage object
	 * \brief operaFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, edatatype Datatype, unsigned Compression)
	 * \brief Constructor to create a new FITS file.
	 */
	operaFITSImage(string Filename, unsigned Naxis1, unsigned Naxis2, 
				   edatatype Datatype, unsigned Compression = 0, bool isLazy = false);	// create a new empty FITSImage
	/*! 
	 * \sa class operaFITSImage
	 * \brief create a FITSIMage object from a FITS file
	 * \brief operaFITSImage(string Filename, int Mode=READWRITE)
	 * \brief Constructor to create a FITSImage from a FITS file.
	 */
	operaFITSImage(string Filename, edatatype Datatype, int Mode/*READWRITE/READONLY*/, unsigned Compression = 0, bool isLazy = false);		// read an existing FITSImage from file
	/*! 
	 * operaFITSImage* operaFITSImage(operaFITSImage &imageIn, bool ViewOnly)
	 * \brief Clone a FITSImage object.
	 */
	operaFITSImage(operaFITSImage &imageIn, bool ViewOnly = false, bool AddHeader = false);			// clone image object
	/*! 
	 * operaFITSImage* operaFITSImage(operaFITSImage &imageIn)
	 * \brief Clone a FITSImage object.
	 */
	~operaFITSImage();		// destructor
	
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return a pixel value.
	 */
	float *operator[](unsigned i) {
#ifdef FITS_RANGE_CHECK  // Check that the given value is within the actual range 
		if ((unsigned long)i > this->naxis2) {
			throw operaException("operaFITSImage: ", MatrixInvalidDimensions, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		// in the case of a lazy read, we must check in the FITSImage [] operator
		// if the extension has been read... This is the case of a reference without an assignment
		if (isLazy) {
			int status = 0;
			if (!(extensionHasBeenRead[current_extension]||extensionHasBeenWritten[current_extension])) {
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
					throw operaException("operaFITSImage: cfitsio error "+filename+" ("+itos(filedatatype)+") ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				if (todatatype(ebitpix(filedatatype), mybzero, bscale) != datatype) {
					operaFITSImageConvertImageInPlace(todatatype(ebitpix(filedatatype), mybzero, bscale), datatype);
				}
				setHasBeenRead(current_extension);
                if (super)
					super->setHasBeenRead(current_extension);;
			}
		}
#ifdef ASL11
		if ((current_extension == 1 && current_slice == 1) || isLazy)
			return (float *)pixptr+(i<<ASL11);
		else
			return ((float *)pixptr+((current_extension-1)<<ASL22)*naxis3)+((current_slice-1)<<ASL22)+(i<<ASL11);
#else
		if ((current_extension == 1 && current_slice == 1) || isLazy)
			return (float *)pixptr+(naxis1*i);
		else
			return ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2+(naxis1*i);
#endif
	};
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return an ImageVector.
	 * \note usage: myfitsimage[where(input>saturation)] = 0.0; sets the image pixel values where true to zero
	 */
	operaImageVector &operator[](operaImageVector *vector) {
        SetImage(vector, this); 
        return *vector;
    };
	/*! 
	 * \brief operator []
	 * \brief indexing operator to return an ImageVector.
	 * \param i - index into row
	 * \note usage: myfitsimage[where(input>saturation)] = 0.0; sets the image pixel values where true to zero
	 */
	operaImageVector &operator[](operaImageVector &vector) {
        SetImage(&vector, this);
        return vector;
    };
	
	/*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \note usage: operaFITSImage a = operaFITSImage b; copies the pixel values from b to a
	 */	
	operaFITSImage& operator=(operaFITSImage* b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b->pixptr;
		unsigned long n = npixels; 
        if (isLazy) {
			int status = 0;
			long fpixel = 1;
			if (!b->extensionHasBeenRead[current_extension]) {
				long fpixel[MAXFITSDIMENSIONS] = {1,1,1};
				int filedatatype;
				if (fits_get_img_equivtype(b->fptr, &filedatatype, &status)) {
					throw operaException("operaFITSImage: cfitsio error "+b->filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				float mybzero = b->bzero;
				if (filedatatype == USHORT_IMG) {
					mybzero = 32768.0;
				}
				if (fits_read_pix(b->fptr, b->todatatype(ebitpix(filedatatype), mybzero, b->bscale), fpixel, b->npixels, NULL, b->pixptr, NULL, &status)) {
					throw operaException("operaFITSImage: cfitsio error "+b->filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				if (b->todatatype(ebitpix(filedatatype), mybzero, b->bscale) != b->datatype) {
					b->operaFITSImageConvertImageInPlace(b->todatatype(ebitpix(filedatatype), mybzero, b->bscale), b->datatype);
				}
				b->setHasBeenRead(current_extension);
				if (b->super)
					b->super->setHasBeenRead(current_extension);
			}
			while (n--) *p++ = *bp++;
			if (fits_write_img(this->fptr, datatype, fpixel, npixels, this->pixptr, &status)) {
				throw operaException("operaFITSImage: cfitsio error "+this->filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			setHasBeenWritten(current_extension);
			if (super)
				super->setHasBeenWritten(current_extension);
		} else {
            if (current_extension > 1 || current_slice > 1)
                p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
            if (b->current_extension > 1 || b->current_slice > 1)
                bp = ((float *)b->pixptr+(b->current_extension-1)*b->naxis1*b->naxis2*b->naxis3)+(b->current_slice-1)*b->naxis1*b->naxis2;			
			// we need to copy all extensions that are stacked in memory
            while (n--) *p++ = *bp++;
        }
        if (b->istemp && !b->viewOnly) delete b;
		return *this;
	};
	/*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \note usage: operaFITSImage a = operaFITSImage b; copies the pixel values from b to a
	 */	
	operaFITSImage& operator=(operaFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		unsigned long n = npixels; 
        if (isLazy) {
			int status = 0;
			long fpixel = 1;
			if (!b.extensionHasBeenRead[b.current_extension]) {
				long fpixel[MAXFITSDIMENSIONS] = {1,1,1};
				int filedatatype;
				if (fits_get_img_equivtype(b.fptr, &filedatatype, &status)) {
					throw operaException("operaFITSImage: cfitsio error "+b.filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
				}
				float mybzero = b.bzero;
				if (filedatatype == USHORT_IMG) {
					mybzero = 32768.0;
				}
				if (fits_read_pix(b.fptr, b.todatatype(ebitpix(filedatatype), mybzero, b.bscale), fpixel, b.npixels, NULL, b.pixptr, NULL, &status)) {
					//throw operaException("operaFITSImage: cfitsio error "+b.filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
					status = 0;
				}
				if (b.todatatype(ebitpix(filedatatype), mybzero, b.bscale) != b.datatype) {
					b.operaFITSImageConvertImageInPlace(b.todatatype(ebitpix(filedatatype), mybzero, b.bscale), b.datatype);
				}
				b.setHasBeenRead(current_extension);
				if (b.super)
					b.super->setHasBeenRead(current_extension);
			}
			while (n--) *p++ = *bp++;
			if (fits_write_img(fptr, datatype, fpixel, npixels, pixptr, &status)) {
				throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			setHasBeenWritten(current_extension);
			if (super)
				super->setHasBeenWritten(current_extension);
		} else {
            if (current_extension > 1 || current_slice > 1)
                p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
            if (b.current_extension > 1 || b.current_slice > 1)
                bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;			
			// we need to copy all extensions that are stacked in memory
            while (n--) *p++ = *bp++;
        }
        if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator =
	 * \brief assignment operator.
	 * \note usage: operaFITSImage a = 0.0; copies the float value to every pixel in = a
	 */	
	operaFITSImage& operator=(float f) {
		float *p = (float *)pixptr;
		unsigned long n = npixels;
        if (isLazy) {
            int status = 0;
			long fpixel = 1;
            while (n--) *p++ = f;
			if (fits_write_img(fptr, datatype, fpixel, npixels, this->pixptr, &status)) {
				throw operaException("operaFITSImage: cfitsio error "+this->filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
			}
			setHasBeenWritten(current_extension);
        } else {
			if (current_extension > 1 || current_slice > 1)
                p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		while (n--) *p++ = f;
		return *this;
	};
	/*! 
	 * \brief operator >>
	 * \brief write image to named file (saveas).
	 * \note usage: operafitsmage a(); a >> "foo.fits"; 
	 */	
	operaFITSImage& operator>>(const string Filename) {
		operaFITSImageSaveAs(Filename);
		return *this;
	};
	/*! 
	 * \brief operator >>
	 * \brief write image to another image (headers and pixels).
	 * \note usage: operafitsmage a(); operafitsmage b(); a >> b; 
	 */	
	const operaFITSImage& operator>>(operaFITSImage& image) {
		image.operaFITSImageCopyHeader(this);
		image = this;
		return image;
	};
	/*! 
	 * \brief operator <<
	 * \brief read an image from a named file.
	 * \note usage: operafitsmage a(); a << "foo.fits"; 
	 */	
	operaFITSImage& operator<<(const string Filename) {
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
	 * \note usage: operafitsmage a(); operafitsmage b(); a << b; 
	 */	
	const operaFITSImage& operator<<(operaFITSImage& image) {
		operaFITSImageCopyHeader(&image);
		*this = image;
		return *this;
	};
	/*! 
	 * \brief operator ==
	 * \brief equality operator.
	 * \note usage:
	 */
	operaFITSImage& operator==(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) {*tp++ = (*p++ == *bp++?1.0:0.0);}
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator ==
	 * \brief equality operator.
	 * \note usage:
	 */
	operaFITSImage& operator==(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) {*tp++ = (*p++ == f?1.0:0.0);}
		return *t;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \note usage: operaFITSImage a += operaFITSImage b; adds the pixel values from b to a
	 */	
	operaFITSImage& operator+=(operaFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *p++ += *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator +=
	 * \brief add/assignment operator.
	 * \note usage: operaFITSImage a += 100.0; adds the float value to a
	 */	
	operaFITSImage& operator+=(float f) {
		float *p = (float *)pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *p++ += f;
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \note usage: operaFITSImage a -= operaFITSImage b; subtracts the pixel values from b to a
	 */	
	operaFITSImage& operator-=(operaFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *p++ -= *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator -=
	 * \brief subtract/assignment operator.
	 * \note usage: operaFITSImage a -= 100.0; subtracts the float value from a
	 */	
	operaFITSImage& operator-=(float f) {
		float *p = (float *)pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *p++ -= f;
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \note usage: operaFITSImage a *= operaFITSImage b; multiplies the pixel values from b to a
	 */	
	operaFITSImage& operator*=(operaFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *p++ *= *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator *=
	 * \brief multiply/assignment operator.
	 * \note usage: operaFITSImage a *= 100.0; multiplies the float value times a
	 */	
	operaFITSImage& operator*=(float f) {
		float *p = (float *)pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *p++ *= f;
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \note usage: operaFITSImage a /= operaFITSImage b; divides the pixel values from b into a
	 */	
	operaFITSImage& operator/=(operaFITSImage& b) {
		float *p = (float *)pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *p++ /= *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *this;
	};
	/*! 
	 * \brief operator /=
	 * \brief divide/assignment operator.
	 * \note usage: operaFITSImage a /= 100.0; divides the float value from a
	 */	
	operaFITSImage& operator/=(float f) {
		float *p = (float *)pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *p++ /= f;
		return *this;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \note usage: operaFITSImage a = operaFITSImage b * operaFITSImage c; multiplies the pixel values  b * c and assigns to a
	 */	
	operaFITSImage& operator*(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = *p++ * *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator *
	 * \brief multiply operator.
	 * \note usage: operaFITSImage a = operaFITSImage b * 10.0; multiplies the pixel values  b * 10.0 and assigns to a
	 */	
	operaFITSImage& operator*(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = *p++ * f;
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \note usage: operaFITSImage a = operaFITSImage b / operaFITSImage c; divides the pixel values  b / c and assigns to a
	 */	
	operaFITSImage& operator/(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = *p++ / *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator /
	 * \brief divide operator.
	 * \note usage: operaFITSImage a = operaFITSImage b / 100.0; divides the pixel values  b / 100.0 and assigns to a
	 */	
	operaFITSImage& operator/(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}unsigned long n = npixels; 
		while (n--) *tp++ = *p++ / f;
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \note usage: operaFITSImage a = operaFITSImage b + operaFITSImage c; adds the pixel values  b + c and assigns to a
	 */	
	operaFITSImage& operator+(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = *p++ + *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator +
	 * \brief add operator.
	 * \note usage: operaFITSImage a = operaFITSImage b + 100.0; adds the pixel values  b + 100.0 and assigns to a
	 */	
	operaFITSImage& operator+(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = *p++ + f;
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \note usage: operaFITSImage a = operaFITSImage b - operaFITSImage c; subtracts the pixel values  b - c and assigns to a
	 */	
	operaFITSImage& operator-(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = *p++ - *bp++;
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator -
	 * \brief subtract operator.
	 * \note usage: operaFITSImage a = operaFITSImage b - 100.0; subtracts the pixel values  b - 100.0 and assigns to a
	 */	
	operaFITSImage& operator-(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = *p++ - f;
		return *t;
	};
	/*! 
	 * \brief operator !
	 * \brief invert operator
	 * \note usage: operaFITSImage a = !operaFITSImage b; inverts the pixel values of b and assigns to a (0.0 becomes 1.0, non-zero becomes 0.0)
	 */	
	operaFITSImage& operator!() {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++==0.0?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator >
	 * \brief greater than operator.
	 * \note usage: operaFITSImage a = operaFITSImage b > operaFITSImage c; creates a mask of 1.0 if b > c or 0.0 if b <= c and assigns to a
	 */	
	operaFITSImage& operator>(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++>*bp++?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator >
	 * \brief greater than operator.
	 * \note usage: operaFITSImage a = operaFITSImage b > 100.0; creates a mask of 1.0 if b > 100.0 or 0.0 if b <= 100.0 and assigns to a
	 */	
	operaFITSImage& operator>(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++>f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator >=
	 * \brief greater than or equal operator.
	 * \note usage: operaFITSImage a = operaFITSImage b >= operaFITSImage c; creates a mask of 1.0 if b >= c or 0.0 if b < c and assigns to a
	 */	
	operaFITSImage& operator>=(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++>=*bp++?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator >=
	 * \brief greater than or equal operator.
	 * \note usage: operaFITSImage a = operaFITSImage b >= 100.0; creates a mask of 1.0 if b >= 100.0 or 0.0 if b < 100.0 and assigns to a
	 */	
	operaFITSImage& operator>=(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++>=f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator <
	 * \brief less than operator.
	 * \note usage: operaFITSImage a = operaFITSImage b < operaFITSImage c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 */	
	operaFITSImage& operator<(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++<*bp++?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator <
	 * \brief less than operator.
	 * \note usage: operaFITSImage a = operaFITSImage b < operaFITSImage c; creates a mask of 1.0 if b < c or 0.0 if b >= c and assigns to a
	 */	
	operaFITSImage& operator<(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++<f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator <=
	 * \brief less than or equal operator.
	 * \note usage: operaFITSImage a = operaFITSImage b <= operaFITSImage c; creates a mask of 1.0 if b <= c or 0.0 if b > c and assigns to a
	 */	
	operaFITSImage& operator<=(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++<=*bp++?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator <=
	 * \brief less than or equal operator.
	 * \note usage: operaFITSImage a = operaFITSImage b <= 100.0; creates a mask of 1.0 if b <= 100.0 or 0.0 if b > 100.0 and assigns to a
	 */	
	operaFITSImage& operator<=(float f) {
		operaFITSImage *t = new operaFITSImage(*this);
		t->istemp = true;
		float *tp = (float *)t->pixptr; 
		float *p = (float *)this->pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++<f?1.0:0.0);
		return *t;
	};
	/*! 
	 * \brief operator &&
	 * \brief less than operator.
	 * \note usage: operaFITSImage a = operaFITSImage b && operaFITSImage c; creates a mask of 1.0 where b && c != 0.0
	 */	
	operaFITSImage& operator&&(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++!=0.0&&*bp++!=0.0?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * \brief operator ||
	 * \brief less than operator.
	 * \note usage: operaFITSImage a = operaFITSImage b && operaFITSImage c; creates a mask of 1.0 where b && c != 0.0
	 */	
	operaFITSImage& operator||(operaFITSImage& b) {
		float *tp;
		operaFITSImage *t = NULL;
		if (this->istemp) {
			t = this;
			tp = (float *)this->pixptr; 
		} else {
			t = new operaFITSImage(*this);
			t->istemp = true;
			tp = (float *)t->pixptr; 
		}
		float *p = (float *)this->pixptr; 
		float *bp = (float *)b.pixptr; 
		if (!isLazy) {
			if (current_extension > 1 || current_slice > 1)
				p = ((float *)pixptr+(current_extension-1)*naxis1*naxis2*naxis3)+(current_slice-1)*naxis1*naxis2;
			if (b.current_extension > 1 || b.current_slice > 1)
				bp = ((float *)b.pixptr+(b.current_extension-1)*b.naxis1*b.naxis2*b.naxis3)+(b.current_slice-1)*b.naxis1*b.naxis2;
		}
		unsigned long n = npixels; 
		while (n--) *tp++ = (*p++!=0.0||*bp++!=0.0?1.0:0.0);
		if (b.istemp && !b.viewOnly) delete &b;
		return *t;
	};
	/*! 
	 * string operaFITSGetFilename(void) 
	 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
	 */
	string operaFITSGetFilename(void);
	/*! 
	 * operaFITSImageSave() 
	 * \brief Saves the current image to disk.
	 * \note This is tricky and does some fancy fptr work so that the HDU is retained.
	 */
	void operaFITSImageSave();												// save a FITSImage
	/*! 
	 * operaFITSImageSaveAs(string newFilename) 
	 * \brief Saves the current image to disk, with the given filename.
	 */
	void operaFITSImageSaveAs(string newFilename);							// save as
	/*! 
	 * unsigned short* operaFITSImageClonePixelsUSHORT()
	 * \brief Clone pixel image data.
	 */
	unsigned short* operaFITSImageClonePixelsUSHORT();						// clone ccd pixels
	/*! 
	 * unsigned short* operaFITSImageClonePixelsUSHORT(unsigned x, unsigned y, unisgned nx, unsigned ny)
	 * \brief Clone a pixel image subwindow.
	 */
	unsigned short* operaFITSImageClonePixelsUSHORT(unsigned x, unsigned y, unsigned nx, unsigned ny); // clone subimage
	/*! 
	 * float* operaFITSImageClonePixels()
	 * \brief Clone float pixel data.
	 * \return pixels*
	 */
	float* operaFITSImageClonePixels();										// clone image pixels
	/*! 
	 * float* operaFITSImageClonePixels(unsigned x, unsigned y, unsigned nx, unsigned ny
	 * \brief Clone pixels.
	 */
	float* operaFITSImageClonePixels(unsigned x, unsigned y, unsigned nx, unsigned ny); // clone subwindow
	/*! 
	 * void operaFITSImage::operaFITSImageSetData(operaFITSSubImage &subImage, unsigned X, unsigned Y)
	 * \brief copy a subimage in to an operaFITSImage.
	 */
	void operaFITSImageSetData(operaFITSSubImage &subImage, unsigned X, unsigned Y);
	/* 
	 * void operaFITSImage::operaFITSImageWriteVirtual(operaFITSSubImage &subImage, unsigned X, unsigned Y)
	 * \brief Write a subimage on to avirtual operaFITSImage.
	 * \brief This is used for deep stacking. 
	 * \brief The virtual image is Lazy, and can be very large.
	 * \brief The actual image exists only on disk.
	 */
	void operaFITSImageWriteVirtual(operaFITSSubImage &subImage, unsigned long X, unsigned long Y);
	/* 
	 * void operaFITSImage::operaFITSImageReadVirtual(operaFITSSubImage &subImage, unsigned X, unsigned Y)
	 * \brief Read a subimage from a virtual operaFITSImage.
	 * \brief This is used for deep stacking. 
	 * \brief The virtual image is Lazy, and can be very large.
	 * \brief The actual image exists only on disk.
	 */
	void operaFITSImageReadVirtual(operaFITSSubImage &subImage, unsigned long X, unsigned long Y);
	/*! 
	 * void operaFITSImage::operaFITSImageSetData(unsigned short* data)
	 * \brief set the iamge data pointer to a buffer of data.
	 * \param data pointer to the data
	 */
	void operaFITSImageSetData(unsigned short* data);						// set pixptr to a data buffer
	/*! 
	 * void operaFITSImage::operaFITSImageSetData(float* data)
	 * \brief set the iamge data pointer to a buffer of data.
	 */
	void operaFITSImageSetData(float* data);								// set pixptr to a data buffer
	/* 
	 * string operaFITSFindComment(string contains) 
	 * \brief returns the value of a FITS comment containing string.
	 * \param keyword
	 * \return string keyword value or empty string
	 */
	string operaFITSFindComment(string contains);
	/*! 
	 * string operaFITSGetHeaderValue(string keyword) 
	 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
	 */
	string operaFITSGetHeaderValue(string keyword);							// return header kweyword value as string
	/*! 
	 * float operaFITSGetFloatHeaderValue(string keyword) 
	 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
	 */
	float operaFITSGetFloatHeaderValue(string keyword);
	
	/*! 
	 * int operaFITSGetFloatHeaderValue(string keyword) 
	 * \brief returns the value of a FITS keyword with leading trailing quotes and trailing spaces removed.
	 */
	int operaFITSGetIntHeaderValue(string keyword);
	
	/*! 
	 * string operaFITSGetRawHeaderValue(string keyword) 
	 * \brief returns the value of a FITS keyword verbatim.
	 */
	string operaFITSGetRawHeaderValue(string keyword);						// return raw header kweyword value as string
	/*! 
	 * operaFITSSetHeaderValue(string keyword, string value, string comment)
	 * \brief sets the given keyword to value with comment.
	 */
	void operaFITSSetHeaderValue(string keyword, string value, string comment); // set header kweyword value as string
	/*! 
	 * operaFITSSetHeaderValue(string keyword, unsigned short value, string comment)
	 * \brief sets the given keyword to value with comment.
	 */
	void operaFITSSetHeaderValue(string keyword, unsigned short value, string comment); // set header keyword value as ushort
	/*! 
	 * operaFITSSetHeaderValue(string keyword, unsigned short value, string comment)
	 * \brief sets the given keyword to value with comment.
	 */
	void operaFITSSetHeaderValue(string keyword, short value, string comment); // set header keyword value as ushort
	/*! 
	 * operaFITSSetHeaderValue(string keyword, float value, string comment)
	 * \brief sets the given keyword to value with comment.
	 */
	void operaFITSSetHeaderValue(string keyword, float value, string comment); // set header keyword value as float
	/*! 
	 * operaFITSSetHeaderValue(string keyword, double value, string comment)
	 * \brief sets the given keyword to value with comment.
	 */
	void operaFITSSetHeaderValue(string keyword, double value, string comment); // set header keyword value as double
	/*! 
	 * operaFITSAddComment(string comment)
	 * \brief adds a comment to header.
	 */
	void operaFITSAddComment(string comment);
	/*! 
	 * operaFITSDeleteHeaderKey(string keyword)
	 * \brief sets the given keyword to value with comment.
	 */
	void operaFITSDeleteHeaderKey(string keyword);
	/*! 
	 * operaFITSImageCopyHeader(operaFITSImage *from) 
	 * \brief Copies all of the header information from image.
	 */
	void operaFITSImageCopyHeader(operaFITSImage *from);					// copy header unit
    /*
     * operaFITSImageCopyHeader(operaFITSImage *from, unsigned inputhdu)
     * \brief Copies all of the header information from image.
     */
    void operaFITSImageCopyHeaderFromDifferentHDU(operaFITSImage *from, unsigned inputabshdu);
	
    /*!
	 * operaFITSImageCopyNonStructuralHeader(operaFITSImage *from) 
	 * \brief Copies all of the non-structural header information from image.
	 */
	void operaFITSImageCopyNonStructuralHeader(operaFITSImage *from);					// copy header unit
	/*! 
	 * void operaFITSImage::operaFITSImageConvertImage(edatatype todatattype)
	 * \brief Convert image type and create new image values of that type.
	 */
	void operaFITSImageConvertImage(edatatype todatatype);					// convert image type and values
	/*! 
	 * void operaFITSImageConvertImage(edatatype fromdatatype, edatatype todatatype)
	 * \brief Convert image type and create new image values of that type.
	 */
	void operaFITSImageConvertImage(edatatype fromdatatype, edatatype todatatype);	// convert image type and values
	/*!
	 * void operaFITSImageConvertImageInPlace(edatatype fromdatatype, edatatype todatatype)
	 * \brief Convert image type and create new image values of that type.
	 * \note - this assumes enough storage is available for the target type.
	 */
	void operaFITSImageConvertImageInPlace(edatatype fromdatatype, edatatype todatatype);
	/*! 
	 * operaFITSImage::rotate90()
	 * \brief rotate 90 degrees.
	 */
	void rotate90();
    /*!
     * operaFITSImage::mirrorColumns()
     * \brief mirror x.
     */
    void mirrorColumns(void);
    
    /*!
     * operaFITSImage::mirrorRows()
     * \brief mirror y.
     */
    void mirrorRows(void);
    
	/*!
	 * operaFITSImage::assignVariances(float gain)
	 * \brief Assign variances to each pixel in an image.
	 */
	void assignVariances(float gain);
	/*! 
	 * operaFITSImage::transpose(operaFITSImage &inImage)
	 * \brief transpose y for x.
	 * \return void
	 */
	void transpose(operaFITSImage &inImage);
	/*! 
	 * operaFITSImageClose()
	 * \brief Close a FITSImage object file.
	 */
	void operaFITSImageClose();												// close
	
	/*! 
	 * getpixelUSHORT(unsigned x, unsigned y, unsigned long naxis1)
	 * \brief get an unsigned short pixel value at coordinates x,y.
	 */
	inline unsigned short getpixelUSHORT(unsigned x, unsigned y, unsigned long naxis1) {return ((unsigned short *)pixptr)[(naxis1*y)+x];};
	
	/*! 
	 * getpixelUSHORT(unsigned short*p, unsigned x, unsigned y, unsigned long naxis1)
	 * \brief get an unsigned short pixel value at coordinates x,y.
	 */
	inline unsigned short getpixelUSHORT(unsigned short*p, unsigned x, unsigned y, unsigned long naxis1) {return ((unsigned short *)p)[(naxis1*y)+x];};
	
	/*! 
	 * getpixelUSHORT(unsigned x, unsigned y)
	 * \brief get an unsigned short pixel value at coordinates x,y.
	 */
	inline unsigned short getpixelUSHORT(unsigned x, unsigned y) {return ((unsigned short *)pixptr)[(naxis1*y)+x];};
	
	/*! 
	 * getpixel(unsigned x, unsigned y, unsigned long naxis1)
	 * \brief get a pixel value at coordinates x,y.
	 */
	inline float getpixel(unsigned x, unsigned y, unsigned long naxis1) {return ((float *)pixptr)[(naxis1*y)+x];};
	
	/*! 
	 * getpixel(float *p, unsigned x, unsigned y, unsigned long naxis1)
	 * \brief get a pixel value at coordinates x,y.
	 */
	inline float getpixel(float *p, unsigned x, unsigned y, unsigned long naxis1) {return ((float *)p)[(naxis1*y)+x];};
	
	/*! 
	 * getpixel(unsigned x, unsigned y)
	 * \brief get a pixel value at coordinates x,y.
	 */
	inline float getpixel(unsigned x, unsigned y) {return ((float *)pixptr)[(naxis1*y)+x];};
	
	/*! 
	 * setpixel(float value, unsigned x, unsigned y)
	 * \brief get a pixel value at coordinates x,y.
	 */
	inline void setpixel(unsigned short value, unsigned x, unsigned y) {((unsigned short *)pixptr)[(naxis1*y)+x] = value;};
	inline void setpixel(unsigned short value, unsigned x, unsigned y, unsigned long naxis1) {((unsigned short *)pixptr)[(naxis1*y)+x] = value;};
	inline void setpixel(unsigned short *p, unsigned short value, unsigned x, unsigned y, unsigned long naxis1) {((unsigned short *)p)[(naxis1*y)+x] = value;};
	inline void setpixel(float value, unsigned x, unsigned y) {((float *)pixptr)[(naxis1*y)+x] = value;};
	inline void setpixel(float value, unsigned x, unsigned y, unsigned long naxis1) {((float *)pixptr)[(naxis1*y)+x] = value;};
	inline void setpixel(float  *p, float value, unsigned x, unsigned y, unsigned long naxis1) {((float *)p)[(naxis1*y)+x] = value;};
	
	/*! 
	 * setHasBeenRead(operaFITSImage &image)
	 * \brief Copy hasBeenRead property.
	 */
	void setHasBeenRead(operaFITSImage &image);
	
	/*! 
	 * setHasBeenRead(unsigned index, bool Read)
	 * \brief Set hasBeenRead property index = 0 means all.
	 */
	void setHasBeenRead(unsigned index, bool Read);
	
	/*! 
	 * setHasBeenRead(unsigned index)
	 * \brief Set hasBeenRead property for this extension, resetting all others.
	 */
	void setHasBeenRead(unsigned index);
	
	/*! 
	 * setHasBeenWritten(operaFITSImage &image
	 * \brief Copy hasBeenWritten property.
	 */
	void setHasBeenWritten(operaFITSImage &image);
	
	/*! 
	 * setHasBeenWritten(unsigned index, bool Read)
	 * \brief Set hasBeenWritten property index = 0 means all.
	 */
	void setHasBeenWritten(unsigned index, bool Read);
	
	/*! 
	 * setHasBeenWritten(unsigned index)
	 * \brief Set hasBeenWritten property for this extension, resetting all others.
	 */
	void setHasBeenWritten(unsigned index);
	
	operaImageVector *where(Box &box);
	
	operaImageVector *where(unsigned x, unsigned y, unsigned dx, unsigned dy);
	
    /*
	 * \sa method resize(unsigned x0, unsigned xf, unsigned y0, unsigned yf)
	 * \brief resize a product image to the given rows, cols, retaining any existing data
     */
    void resize(unsigned x0, unsigned xf, unsigned y0, unsigned yf);

	/*! 
	 * \sa method resize(unsigned cols, unsigned rows)
	 * \brief resize a product image to the given rows, cols, retaining any existing data
	 */
	void resize(unsigned cols, unsigned rows);
	
	/* getters and setters */
	
	/*! 
	 * setfilename(string Filename)
	 * \brief set the image filename.
	 */
	void setfilename(string name);
	/*! 
	 * string getfilename()
	 * \brief set the image filename.
	 */
	string getfilename();
	/*! 
	 * row(unsigned Row)
	 * \brief get base of nth row.
	 */
	float* row(unsigned Row);
	/*! 
	 * fitsfile *getfitsfileptr()
	 * \brief get the image FITS file pointer.
	 */
	fitsfile *getfitsfileptr();
	/*! 
	 * void setfitsfileptr(fitsfile *newfptr)
	 * \brief set the image FITS file pointer.
	 */
	void setfitsfileptr(fitsfile *newfptr);
	/*! 
	 * void *getpixels()
	 * \brief get the image array base.
	 */
	void *getpixels();	
	/*! 
	 * void *getvars()
	 * \brief get the variance array base.
	 */
	void *getvars();	
	/*! 
	 * unsigned getnaxis()
	 * \brief get the image array number of axes (2 usually).
	 */
	unsigned getnaxis();
	/*! 
	 * unsigned getnaxis1()
	 * \brief get the image array length of x dimension.
	 */
	unsigned getnaxis1();
	/*! 
	 * unsigned getnaxis2()
	 * \brief get the image array length of y dimension.
	 */
	unsigned getnaxis2();
	/*! 
	 * unsigned getnaxis3()
	 * \brief get the image array length of y dimension.
	 * \return unsigned length of axis 3
	 */
	unsigned getnaxis3();
	/*! 
	 * unsigned getXDimension()
	 * \brief get the image array length of x dimension.
	 */
	unsigned getXDimension();
	/*! 
	 * unsigned getYDimension()
	 * \brief get the image array length of y dimension.
	 */
	unsigned getYDimension();
	/*! 
	 * unsigned getZDimension()
	 * \brief get the image array length of z dimension.
	 */
	unsigned getZDimension();
	/*! 
	 * unsigned getNExtensions()
	 * \brief get the number of extensions in the file.
	 */
	unsigned getNExtensions();
	/* 
	 * unsigned getNPixelsPerExtension()
	 * \brief get the image array number of pixels in an extension.
	 * \return unsigned number of pixels
	 */
	unsigned getNPixelsPerExtension();
	/* 
	 * unsigned getExtension()
	 * \brief get the current extension.
	 * \return unsigned extension
	 */
	unsigned getExtension();
	/*! 
	 * unsigned getnpixels()
	 * \brief get the image array number of pixels.
	 */
	unsigned getnpixels();
	/*! 
	 * unsigned long getsize()
	 * \brief get the size of an image.
	 */
	unsigned long getsize();
	/*! 
	 * unsigned long getelementsize()
	 * \brief get the size of a pixel.
	 */
	unsigned getelementsize();
	/*! 
	 * unsigned getCompression()
	 * \brief get the current compression method.
	 */
	unsigned getCompression();
	/*! 
	 * void getCompression(nsigned Compression)
	 * \brief set the current compression method.
	 */
	void setCompression(unsigned Compression);
	/*! 
	 * edatatype getdatatype()
	 * \brief get the current datatype
	 */
	edatatype getdatatype();
	/*! 
	 * ebitpix getbitpix()
	 * \brief get the current bitpix
	 */
	ebitpix getbitpix();
	/*! 
	 * bool getIstemp()
	 * \brief return true if this is a temp.
	 */
	bool getIstemp() {return istemp;};
	/*!
	 * bool memoryAvailable()
	 * \brief Is there enough memory available to read the file in?
	 */
	bool memoryAvailable(long bytes);
	
};

/*! 
 * void map(operaFITSImage &, float from, float to)
 * \brief map all pixels from into to.
 * usage: map(image, F_NAN, 0.0);
 */
static inline void map(operaFITSImage &image, float from, float to) {
	float *p = (float *)image.getpixels(); 
	unsigned n = image.getnpixels(); 
	while (n--) {*p=(*p==from?to:*p);p++;}
}

#endif

