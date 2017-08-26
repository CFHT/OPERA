/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaEspadonsIMage
 Version: 1.0
 Description: class encapsulates a FITS image.
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


/*! operaEspadonsImage
 * \file operaEspadonsImage.cpp
 * \author Doug Teeple
 * \brief This class extands the FITS image to espadons image.
 * \note  Keywards relevant to espadons images are made first-class objects.
 * \ingroup libraries
 */

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaEspadonsImage.h"

#include "libraries/operaImage.h"

using namespace std;

static string imagetypestrings[]		= {"GENERIC", "BIAS", "FLAT", "COMPARISON", "ALIGN", "OBJECT", "DARK"};
static string detectorstrings[]			= {"UNKNOWN", "EEV1", "OLAPA"};
static string amplifierstrings[]		= {"UNKNOWN", "a", "ab"};
static string instrumentmodestrings[]	= {"UNKNOWN", "STAR_ONLY", "STAR_PLUS_SKY", "POLAR"};
static string speedstrings[]			= {"UNKNOWN", "FAST", "NORMAL", "SLOW"};
static string stokesstrings[]			= {"STOKES_UNKNOWN", "STOKES_I", "STOKES_U", "STOKES_Q", "STOKES_V"};
static string polarquadstrings[]		= {"UNKNOWN", "1", "2", "3", "4"};

/* 
 * \class operaEspadonsImage()
 * \brief Basic operaEspadonsImage constructor.
 * \extends operaFITSImage
 * \return none
 */
operaEspadonsImage::operaEspadonsImage() : operaFITSImage(),
	x2(0),							// last col to be figured out from DATASEC 
	y2(0),							// last row to be figured out from DATASEC
	x1(0),							// start col to be figured out from DATASEC 
	y1(0),							// start row to be figured out from DATASEC
	ndpixels(0),					// number of DATASEC pixels	
	datasecSubImage(NULL),			// datasec sub Image	
	imagetype(IMTYPE_GENERIC),		// (IMTYPE_BIAS, IMTYPE_FLAT, etc)
	detector(DETECTOR_UNKNOWN),		// Detector EEV1/OLAPA
	amplifier(AMPLIFIER_UNKNOWN),	// a ab
	mode(MODE_UNKNOWN),				// sp1 sp2 pol
	speed(SPEED_UNKNOWN),			// sequence # of polar image
	stokes(STOKES_UNKNOWN),			// (Stokes_(I, U, Q, V))
	sequence(POLAR_QUAD_UNKNOWN)	// sequence # of polar image
{
}

/* 
 * \class operaEspadonsImage
 * \brief operaEspadonsImage(string Filename, int mode)
 * \brief Constructor for readng a FITS file and creating the corresponding object.
 * \param Filename
 * \param mode - READWRITE or READONLY, defaults to READWRITE
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaEspadonsImage::operaEspadonsImage(string Filename, int mode) : operaFITSImage(Filename, mode),
	x2(0),							// last col to be figured out from DATASEC 
	y2(0),							// last row to be figured out from DATASEC
	x1(0),							// start col to be figured out from DATASEC 
	y1(0),							// start row to be figured out from DATASEC
	ndpixels(0),					// number of DATASEC pixels	
	datasecSubImage(NULL),			// datasec sub Image	
	imagetype(IMTYPE_GENERIC),		// (IMTYPE_BIAS, IMTYPE_FLAT, etc)
	detector(DETECTOR_UNKNOWN),		// Detector EEV1/OLAPA
	amplifier(AMPLIFIER_UNKNOWN),	// a ab
	mode(MODE_UNKNOWN),				// sp1 sp2 pol
	speed(SPEED_UNKNOWN),			// sequence # of polar image
	stokes(STOKES_UNKNOWN),			// (Stokes_(I, U, Q, V))
	sequence(POLAR_QUAD_UNKNOWN)	// sequence # of polar image
{
	int status = 0;
	string value;
	char zvalue[FLEN_VALUE],zdatasec[FLEN_VALUE],comment[FLEN_COMMENT];
	
	//	get DATASEC from FITS file    
	if (fits_read_keyword(fptr, "DATASEC", zdatasec, comment, &status)) {
		x1 = 1;
		y1 = 1;
		x2 = naxes[0];
		y2 = naxes[1];
	} else {
		sscanf(zdatasec, "'[%u:%u,%u:%u]'", &x1, &x2, &y1, &y2);
	}
	ndpixels = (x2-x1+1) * (y2-y1+1);
	
	datasecSubImage = new operaFITSSubImage((operaFITSImage&)*this, x1, y1, x2-x1, y2-y1);
	
	// read header stuff here
	value = operaFITSGetHeaderValue("DETECTOR");
	if (value == "OLAPA") {
		detector = DETECTOR_OLAPA;
		value = operaFITSGetHeaderValue("AMPLIST");
		if (value == string("a"))
			amplifier = AMPLIFIER_a;
		else if (value == string("a,b"))
			amplifier = AMPLIFIER_ab;
		else
			throw operaException("operaEspadonsImage: missing AMPLIST header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	} else if (value == "EEV1") {
		detector = DETECTOR_EEV1;
		amplifier = AMPLIFIER_a;
	} else
		throw operaException("operaEspadonsImage: missing DETECTOR header keyword "+value+" ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	
	value = operaFITSGetHeaderValue("OBSTYPE");
	if (value == string("OBJECT")) 
		imagetype = IMTYPE_OBJECT;
	else if (value == string("FLAT")) 
		imagetype = IMTYPE_FLAT;
	else if (value == string("BIAS")) 
		imagetype = IMTYPE_BIAS;
	else if (value == string("COMPARISON")) 
		imagetype = IMTYPE_COMPARISON;
	else if (value == string("ALIGN")) 
		imagetype = IMTYPE_ALIGN;
	else if (value == string("DARK")) 
		imagetype = IMTYPE_DARK;
	else
		throw operaException("operaEspadonsImage: missing OBSTYPE header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	
	value = operaFITSGetHeaderValue("INSTMODE");
	if (value == string("Polarimetry, R=65,000"))
		mode = MODE_POLAR;
	else if (value == string("Spectroscopy, star only, R=80,000"))
		mode = MODE_STAR_ONLY;
	else if (value == string("Spectroscopy, star+sky, R=65,000"))
		mode = MODE_STAR_PLUS_SKY;
	else
		throw operaException("operaEspadonsImage: missing INSTMODE header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	
	value = operaFITSGetHeaderValue("EREADSPD");
	if (string("Fast:") <= value)
		speed = SPEED_FAST;
	else if (string("Normal:") <= value)
		speed = SPEED_NORMAL;
	else if (string("Slow:") <= value)
		speed = SPEED_SLOW;
	else
		throw operaException("operaEspadonsImage: missing EREADSPD header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	

	if ( mode == MODE_POLAR ) {
		// Stokes CMMTSEQ = V exposure 1, sequence 1 of 1
		if (fits_read_keyword(fptr, "CMMTSEQ", zvalue, comment, &status)) {
			throw operaException("operaEspadonsImage: missing CMMTSEQ header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
		}	
		// Polar Sequence CMMTSEQ = V exposure 1, sequence 1 of 1
		unsigned exp, seq, of;
		char stokes_c;
		if (sscanf(zvalue, "'%c exposure %u, sequence %u of %u'", &stokes_c, &exp, &seq, &of)) {
			switch (stokes_c) {
				case 'I':
					stokes = STOKES_I;
					break;
				case 'U':
					stokes = STOKES_U;
					break;
				case 'V':
					stokes = STOKES_V;
					break;
				case 'Q':
					stokes = STOKES_Q;
					break;
				default:
					break;
			}
			if (of == 4) {
				sequence = (polarquad_t)seq;
			}
		}
	}
}

/* 
 * \class operaEspadonsImage
 * \brief operaEspadonsImage(string Filename, int mode=READWRITE)
 * \brief create a espadons Image object from a FITS file
 * \brief Constructor to create a espadons Image from a FITS file.
 * \param Filename
 * \param Datatype - convert image to this datatype
 * \param mode
 * \return void
 */
operaEspadonsImage::operaEspadonsImage(string Filename, edatatype Datatype, int mode/*=READWRITE | READONLY*/) : operaFITSImage(Filename, Datatype, mode),		// read an existing FITSImage from file
	x2(0),							// last col to be figured out from DATASEC 
	y2(0),							// last row to be figured out from DATASEC
	x1(0),							// start col to be figured out from DATASEC 
	y1(0),							// start row to be figured out from DATASEC
	ndpixels(0),					// number of DATASEC pixels	
	datasecSubImage(NULL),			// datasec sub Image	
	imagetype(IMTYPE_GENERIC),		// (IMTYPE_BIAS, IMTYPE_FLAT, etc)
	detector(DETECTOR_UNKNOWN),		// Detector EEV1/OLAPA
	amplifier(AMPLIFIER_UNKNOWN),	// a ab
	mode(MODE_UNKNOWN),				// sp1 sp2 pol
	speed(SPEED_UNKNOWN),			// sequence # of polar image
	stokes(STOKES_UNKNOWN),			// (Stokes_(I, U, Q, V))
	sequence(POLAR_QUAD_UNKNOWN)	// sequence # of polar image
{
	int status = 0;
	string value;
	char zvalue[FLEN_VALUE],zdatasec[FLEN_VALUE],comment[FLEN_COMMENT];
	
	//	get DATASEC from FITS file    
	if (fits_read_keyword(fptr, "DATASEC", zdatasec, comment, &status)) {
		x1 = 1;
		y1 = 1;
		x2 = naxes[0];
		y2 = naxes[1];
	} else {
		sscanf(zdatasec, "'[%u:%u,%u:%u]'", &x1, &x2, &y1, &y2);
	}
	ndpixels = (x2-x1+1) * (y2-y1+1);
	
	datasecSubImage = new operaFITSSubImage((operaFITSImage&)*this, x1, y1, x2-x1, y2-y1);
	
	// read header stuff here
	if (value == "OLAPA") {
		detector = DETECTOR_OLAPA;
		value = operaFITSGetHeaderValue("AMPLIST");
		if (value == string("a"))
			amplifier = AMPLIFIER_a;
		else if (value == string("a,b"))
			amplifier = AMPLIFIER_ab;
		else
			throw operaException("operaEspadonsImage: missing AMPLIST header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	} else if (value == "EEV1") {
		detector = DETECTOR_EEV1;
		amplifier = AMPLIFIER_a;
	} else
		throw operaException("operaEspadonsImage: missing DETECTOR header keyword "+value+" ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	
	value = operaFITSGetHeaderValue("OBSTYPE");
	if (value == string("OBJECT")) 
		imagetype = IMTYPE_OBJECT;
	else if (value == string("FLAT")) 
		imagetype = IMTYPE_FLAT;
	else if (value == string("BIAS")) 
		imagetype = IMTYPE_BIAS;
	else if (value == string("COMPARISON")) 
		imagetype = IMTYPE_COMPARISON;
	else if (value == string("ALIGN")) 
		imagetype = IMTYPE_ALIGN;
	else if (value == string("DARK")) 
		imagetype = IMTYPE_DARK;
	else
		throw operaException("operaEspadonsImage: missing OBSTYPE header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	
	value = operaFITSGetHeaderValue("INSTMODE");
	if (value == string("Polarimetry, R=65,000"))
		mode = MODE_POLAR;
	else if (value == string("Spectroscopy, star only, R=80,000"))
		mode = MODE_STAR_ONLY;
	else if (value == string("Spectroscopy, star+sky, R=65,000"))
		mode = MODE_STAR_PLUS_SKY;
	else
		throw operaException("operaEspadonsImage: missing INSTMODE header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	
	value = operaFITSGetHeaderValue("EREADSPD");
	if (string("Fast:") <= value)
		speed = SPEED_FAST;
	else if (string("Normal:") <= value)
		speed = SPEED_NORMAL;
	else if (string("Slow:") <= value)
		speed = SPEED_SLOW;
	else
		throw operaException("operaEspadonsImage: missing EREADSPD header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	
	// Stokes CMMTSEQ = V exposure 1, sequence 1 of 1
	if (fits_read_keyword(fptr, "CMMTSEQ", zvalue, comment, &status)) {
		throw operaException("operaEspadonsImage: missing CMMTSEQ header keyword ", operaErrorHeaderProblem, __FILE__, __FUNCTION__, __LINE__);	
	}	
	// Polar Sequence CMMTSEQ = V exposure 1, sequence 1 of 1
	unsigned exp, seq, of;
	char stokes_c;
	if (sscanf(zvalue, "'%c exposure %u, sequence %u of %u'", &stokes_c, &exp, &seq, &of)) {
		switch (stokes_c) {
			case 'I':
				stokes = STOKES_I;
				break;
			case 'U':
				stokes = STOKES_U;
				break;
			case 'V':
				stokes = STOKES_V;
				break;
			case 'Q':
				stokes = STOKES_Q;
				break;
			default:
				break;
		}
		if (of == 4) {
			sequence = (polarquad_t)seq;
		}
	}
}

/* 
 * \class operaEspadonsImage
 * \brief operaEspadonsImage(string Filename, unsigned Naxis1, unsigned Naxis2, DATASEC_t &datasec, edatatype Datatype, unsigned Compression)
 * \brief Constructor to create a new FITS file.
 * \param Filename to create (file is deleted if it exists)
 * \param Naxis1 - x dimension
 * \param Naxis2 - y dimension
 * \param datasec - the datasec where the valid data pixels are
 * \param Datatype defaults to tshort
 * \param Compression, defaults to none
 * \throws operaException cfitsio error code
 * \return void
 */
operaEspadonsImage::operaEspadonsImage(string Filename, unsigned Naxis1, unsigned Naxis2, DATASEC_t &datasec,
									   edatatype Datatype, unsigned Compression) : 
	operaFITSImage(Filename, Naxis1, Naxis2, Datatype, Compression),
	x2(0),							// last col to be figured out from DATASEC 
	y2(0),							// last row to be figured out from DATASEC
	x1(0),							// start col to be figured out from DATASEC 
	y1(0),							// start row to be figured out from DATASEC
	ndpixels(0),					// number of DATASEC pixels	
	datasecSubImage(NULL),			// datasec sub Image	
	imagetype(IMTYPE_GENERIC),		// (IMTYPE_BIAS, IMTYPE_FLAT, etc)
	detector(DETECTOR_UNKNOWN),		// Detector EEV1/OLAPA
	amplifier(AMPLIFIER_UNKNOWN),	// a ab
	mode(MODE_UNKNOWN),				// sp1 sp2 pol
	speed(SPEED_UNKNOWN),			// sequence # of polar image
	stokes(STOKES_UNKNOWN),			// (Stokes_(I, U, Q, V))
	sequence(POLAR_QUAD_UNKNOWN)	// sequence # of polar image
{
	
	x1 = datasec.x1;
	y1 = datasec.y1;
	x2 = datasec.x2;
	y2 = datasec.y2;
	ndpixels = (x2-x1+1) * (y2-y1+1);
	
	datasecSubImage = new operaFITSSubImage((operaFITSImage&)*this, x1, y1, x2-x1, y2-y1);
}

/* 
 * \class operaEspadonsImage
 * \brief construct an in-memory espadons Image object
 * \brief operaEspadonsImage(unsigned Naxis1, unsigned Naxis2, edatatype Datatype=tushort, unsigned Compression=0)
 * \brief Create an in-memory image of given dimensions.
 * \param Naxis1 - x ccd dimension
 * \param Naxis2 - y ccd dimension
 * \param datasec DATASEC_t - espadons image datasec
 * \param Datatype optional datatype defaults to tshort
 * \param Compression optional compression, defaults to none
 * \return void
 */
operaEspadonsImage::operaEspadonsImage(unsigned Naxis1, unsigned Naxis2, DATASEC_t &datasec,
									   edatatype Datatype) : 
	operaFITSImage(Naxis1, Naxis2, Datatype),
	x2(0),							// last col to be figured out from DATASEC 
	y2(0),							// last row to be figured out from DATASEC
	x1(0),							// start col to be figured out from DATASEC 
	y1(0),							// start row to be figured out from DATASEC
	ndpixels(0),					// number of DATASEC pixels	
	datasecSubImage(NULL),			// datasec sub Image	
	imagetype(IMTYPE_GENERIC),		// (IMTYPE_BIAS, IMTYPE_FLAT, etc)
	detector(DETECTOR_UNKNOWN),		// Detector EEV1/OLAPA
	amplifier(AMPLIFIER_UNKNOWN),	// a ab
	mode(MODE_UNKNOWN),				// sp1 sp2 pol
	speed(SPEED_UNKNOWN),			// sequence # of polar image
	stokes(STOKES_UNKNOWN),			// (Stokes_(I, U, Q, V))
	sequence(POLAR_QUAD_UNKNOWN)	// sequence # of polar image
{
	
	x1 = datasec.x1;
	y1 = datasec.y1;
	x2 = datasec.x2;
	y2 = datasec.y2;
	ndpixels = (x2-x1+1) * (y2-y1+1);
	
	datasecSubImage = new operaFITSSubImage((operaFITSImage&)*this, x1, y1, x2-x1, y2-y1);
}

/* 
 * operaEspadonsImage()
 * \brief Destructor, frees image memory and closes FITS file.
 * \return none
 */
operaEspadonsImage::~operaEspadonsImage() {	
}

/* 
 * operaEspadonsImage* operaEspadonsImage::operaEspadonsImageClone(operaEspadonsImage &imageIn)
 * \brief Clone an espadons Image object.
 * \param imageIn - pointer to image to clone
 * \return operaEspadonsImage*
 */
operaEspadonsImage* operaEspadonsImage::operaEspadonsImageClone(operaEspadonsImage &imageIn) {
	DATASEC_t datasec = {imageIn.x1, imageIn.x2, imageIn.y1, imageIn.y2};
	operaEspadonsImage *image = new operaEspadonsImage(imageIn.naxis1, imageIn.naxis2, datasec, imageIn.datatype);
	
	image->imagetype = imageIn.imagetype;
	image->detector = imageIn.detector;
	image->amplifier = imageIn.amplifier;
	image->mode = imageIn.mode;
	image->speed = imageIn.speed;
	image->stokes = imageIn.stokes;
	image->sequence = imageIn.sequence;
	return image;
}

/* 
 * unsigned short* operaEspadonsImage::operaEspadonsImageCloneDatasecUSHORT()
 * \brief Clone a pixel datasec.
 * \return pixels*
 */
unsigned short* operaEspadonsImage::operaEspadonsImageCloneDatasecUSHORT() {
	long ndpixels = (x2-x1+1) * (y2-y1+1);
	long size = sizeof(unsigned short)*ndpixels;
	unsigned short *p = (unsigned short *)malloc(size); 
	if (!p) {
		throw operaException("operaEspadonsImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (pixptr)
		memcpy(p, pixptr, size);
	else
		memset(p, 0, size);
	return p;
}

/* 
 * float* operaEspadonsImage::operaEspadonsImageCloneDatasec()
 * \brief Clone a pixel datasec.
 * \return pixels*
 */
float* operaEspadonsImage::operaEspadonsImageCloneDatasec() {
	long ndpixels = (x2-x1+1) * (y2-y1+1);
	long size = sizeof(float)*ndpixels;
	float *p = (float *)malloc(size); 
	if (!p) {
		throw operaException("operaEspadonsImage: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (pixptr)
		memcpy(p, pixptr, size);
	else
		memset(p, 0, size);
	return p;
}

/* 
 * operaEspadonsImageCopyHeader(operaEspadonsImage *from) 
 * \brief Copies all of the header information from image.
 * \param from
 * \throws operaException cfitsio error code
 * \return void
 */
void operaEspadonsImage::operaEspadonsImageCopyHeader(operaEspadonsImage *from) {
	int status = 0;
	
	if (fits_copy_hdu(from->getfitsfileptr(), fptr, 0, &status)) {
		throw operaException("operaFITSImage: cfitsio error "+filename+" ", (operaErrorCode)status, __FILE__, __FUNCTION__, __LINE__);	
	}	
}

/* 
 * void operaEspadonsImage::operaFITSImageSetData(operaFITSSubImage &subImage)
 * \brief copy a subimage back in to an operaEspadonsImage datasec.
 * \param datasecSubImage the datasec sub image
 * \return void
 */
void operaEspadonsImage::operaFITSImageSetData(operaFITSSubImage &datasecSubImage) {
	unsigned nx = datasecSubImage.nx;
	unsigned ny = datasecSubImage.ny;
	for (unsigned j=0; j<ny; j++) {
		for (unsigned i=0; i<nx; i++) {
			setpixel(datasecSubImage.getpixel(i, j), i+x1, j+y1);
		}
	}
}

/*
 * getters / setters
 */
/* 
 * operaFITSSubImage *getDatasecSubImage()
 * \brief get the datasec subImage.
 * \return subImage pointers
 */
operaFITSSubImage *operaEspadonsImage::getDatasecSubImage() {
	return datasecSubImage;
}
/* 
 * unsigned getnx()
 * \brief get the DATASEC array length of dimension 1.
 * \return unsigned length of axis 1
 */
unsigned operaEspadonsImage::getnx() {
	return x2-x1;
}
/* 
 * unsigned getny()
 * \brief get the DATASEC array length of dimension 2.
 * \return unsigned length of axis 2
 */
unsigned operaEspadonsImage::getny() {
	return y2-y1;
}
/* 
 * unsigned getndpixels()
 * \brief get the image array number of pixels.
 * \return unsigned number of pixels
 */
unsigned operaEspadonsImage::getndpixels() {
	return ndpixels;
}
/* 
 * imtype_t getimtype()
 * \brief Returns the image type (OBJECT, FLAT, ...).
 * \return image type
 */
imtype_t operaEspadonsImage::getimtype() {
	return imagetype;
}

/* 
 * imtype_t getimtypestring()
 * \brief Returns the image type (OBJECT, FLAT, ...) as a string.
 * \return image type
 */
string operaEspadonsImage::getimtypestring() {
	return imagetypestrings[imagetype];
}

/* 
 * detector_t getdetector()
 * \brief Returns te detector type (EEV1, OLAPA, ...).
 * \return detector type
 */
detector_t operaEspadonsImage::getdetector() {
	return detector;
}
/* 
 * detector_t getdetector()
 * \brief Returns te detector type (EEV1, OLAPA, ...) as a string.
 * \return detector type
 */
string operaEspadonsImage::getdetectorstring() {
	return detectorstrings[detector];
}

/* 
 * amplifier_t getamplifier()
 * \brief Returns te amplifier type (a, ab ...).
 * \return amplifier type
 */
amplifier_t operaEspadonsImage::getamplifier() {
	return amplifier;
}
/* 
 * amplifier_t getamplifier()
 * \brief Returns te amplifier type (a, ab ...)  as a string.
 * \return amplifier type
 */
string operaEspadonsImage::getamplifierstring() {
	return amplifierstrings[amplifier];
}

/* 
 * instrumentmode_t getmode()
 * \brief Returns the mode (sp1, sp2, pol ...).
 * \return mode
 */
instrumentmode_t operaEspadonsImage::getmode() {
	return mode;
}
/* 
 * string getmodestring()
 * \brief Returns the mode (sp1, sp2, pol ...) as a string.
 * \return mode
 */
string operaEspadonsImage::getmodestring() {
	return instrumentmodestrings[mode];
}

/* 
 * speed_t getspeed()
 * \brief Returns the speed (Fast, Slow, Normal ...).
 * \return speed
 */
speed_t operaEspadonsImage::getspeed() {
	return speed;
}
/* 
 * string getspeedstring()
 * \brief Returns the speed (Fast, Slow, Normal ...) as a string.
 * \return speed
 */
string operaEspadonsImage::getspeedstring() {
	return speedstrings[speed];
}

/* 
 * speed_t getstokes()
 * \brief Returns the stokes parameter (U, I, Vl ...).
 * \return stokes
 */
stokes_t operaEspadonsImage::getstokes() {
	return stokes;
}
/* 
 * string getstokesstring()
 * \brief Returns the stokes parameter (U, I, Vl ...) as a string.
 * \return stokes
 */
string operaEspadonsImage::getstokesstring() {
	return stokesstrings[stokes];
}

/* 
 * polarquad_t getpolarquad()
 * \brief Returns the polar quad sequence number (1,2,3,4).
 * \return polar sequence number
 */
polarquad_t operaEspadonsImage::getpolarquad() {
	return sequence;
}
/* 
 * string getpolarquadstring()
 * \brief Returns the polar quad sequence number (1,2,3,4) as a string.
 * \return polar sequence number
 */
string operaEspadonsImage::getpolarquadstring() {
	return polarquadstrings[sequence];
}

