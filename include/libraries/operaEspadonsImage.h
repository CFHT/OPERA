#ifndef OPERAESPADONSIMAGE_H
#define OPERAESPADONSIMAGE_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaEspadonsImage
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jul/2011
 
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

#include <fitsio.h>
#include "operaFITSImage.h"
#include "operaFITSSubImage.h"

#define ESPADONS_DEFAULT_NAXIS1 2080
#define ESPADONS_DEFAULT_NAXIS2 4640

using namespace std;

/*! 
 * operaEspadonsImage
 * \author Doug Teeple
 * \brief This class extands the FITS image to espadons image.
 * \note  Keywards relevant to espadons images are made first-class objects.
 * \file operaEspadonsImage.h
 * \ingroup libraries
 */
typedef enum { IMTYPE_GENERIC, IMTYPE_BIAS, IMTYPE_FLAT, IMTYPE_COMPARISON, IMTYPE_ALIGN, IMTYPE_OBJECT, IMTYPE_DARK } imtype_t ;
typedef enum { DETECTOR_UNKNOWN, DETECTOR_EEV1, DETECTOR_OLAPA} detector_t;
typedef enum { AMPLIFIER_UNKNOWN, AMPLIFIER_a, AMPLIFIER_ab} amplifier_t;
typedef enum { MODE_UNKNOWN, MODE_STAR_ONLY, MODE_STAR_PLUS_SKY, MODE_POLAR} instrumentmode_t;
typedef enum { SPEED_UNKNOWN, SPEED_FAST, SPEED_NORMAL, SPEED_SLOW} speed_t;
typedef enum { STOKES_UNKNOWN, STOKES_I, STOKES_Q, STOKES_U, STOKES_V} stokes_t;
typedef enum { POLAR_QUAD_UNKNOWN, POLAR_QUAD_1, POLAR_QUAD_2, POLAR_QUAD_3, POLAR_QUAD_4} polarquad_t;

#define ESPADONS_DEFAULT_DATASEC {21,2068,1,4608}

class operaEspadonsImage : public operaFITSImage {
	
private:
	typedef operaFITSImage& super;	// a way of referring to the super class
	unsigned x2;			// last DATASEC pixel in x direction 
	unsigned y2;			// last DATASEC pixel in y direction
	unsigned x1;			// first DATASEC pixel in x direction 
	unsigned y1;			// first DATASEC pixel in y direction
	unsigned ndpixels;		// number of DATASEC pixels	
	operaFITSSubImage *datasecSubImage;	// the datasec subImage
	imtype_t imagetype;		// (IMTYPE_BIAS, IMTYPE_FLAT, etc)  
	detector_t detector;	// Detector EEV1/OLAPA
	amplifier_t amplifier;	// a ab
	instrumentmode_t mode;	// sp1 sp2 pol
	speed_t speed;			// sequence # of polar image
	stokes_t stokes;		// (Stokes_(I, U, Q, V))
	polarquad_t sequence;	// sequence # of polar image
	
public:
	/*! 
	 * \class operaEspadonsImage()
	 * \brief Basic operaEspadonsImage constructor.
	 * \extends operaFITSImage
	 * \return none
	 */
	operaEspadonsImage();
	
	/*! 
	 * \class operaEspadonsImage
	 * \brief Constructor for readng a FITS file and creating the corresponding object.
	 * \param Filename
	 * \throws operaException operaErrorHeaderProblem
	 * \return none
	 */
	operaEspadonsImage(string Filename, int mode=READWRITE/*READONLY*/);			// constructor to read an existing FITSImage
	/*! 
	 * \class operaEspadonsImage
	 * \brief create a espadons Image object from a FITS file
	 * \brief operaEspadonsImage(string Filename, int mode=READWRITE)
	 * \brief Constructor to create a espadons Image from a FITS file.
	 * \param Filename
	 * \param Datatype - convert image to this datatype
	 * \param mode
	 * \return void
	 */
	operaEspadonsImage(string Filename, edatatype Datatype, int mode=READWRITE/*READONLY*/);		// read an existing FITSImage from file
	/*! 
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
	operaEspadonsImage(string Filename, unsigned Naxis1, unsigned Naxis2, DATASEC_t &datasec,
					   edatatype Datatype = tushort, unsigned Compression = 0);			// constructor to create an in-memory espadons image with a datasec
	/*! 
	 * \class operaEspadonsImage
	 * \brief construct an in-memory espadons Image object
	 * \brief operaEspadonsImage(unsigned Naxis1, unsigned Naxis2, edatatype Datatype=tushort, unsigned Compression=0)
	 * \brief Create an in-memory image of given dimensions.
	 * \param Naxis1 - x ccd dimension
	 * \param Naxis2 - y ccd dimension
	 * \param Datatype optional datatype defaults to tshort
	 * \param Compression optional compression, defaults to none
	 * \return void
	 */
	operaEspadonsImage(unsigned Naxis1, unsigned Naxis2, DATASEC_t &datasec,
				   edatatype Datatype = tushort);	// simply construct an image of a given size
	/*! 
	 * ~operaEspadonsImage()
	 * \brief Destructor, frees image memory and closes FITS file.
	 * \return none
	 */
	~operaEspadonsImage();													// destructor
	
	/*! 
	 * unsigned short* operaEspadonsImage::operaEspadonsImageCloneDatasecUSHORT()
	 * \brief Clone a pixel datasec.
	 * \return pixels*
	 */
	unsigned short* operaEspadonsImageCloneDatasecUSHORT();					// clone the datasec as USHURT
	/*! 
	 * float* operaEspadonsImage::operaEspadonsImageCloneDatasec()
	 * \brief Clone a pixel datasec.
	 * \return pixels*
	 */
	float* operaEspadonsImageCloneDatasec();								// clone the datasec as float
	
	/*! 
	 *  getpixel(unsigned x, unsigned y)
	 * \brief get a pixel value at coordinates x,y.
	 * \param x - the x coordinate
	 * \param y - the y coordinate
	 * \note This function does a shift left by 11 as the image width is known to be 2048 pixels
	 * \note and a shift is faster than a multiply. If the espadons image width ever changes
	 * \note this will have to be revisited.
	 * \return unsigned short
	 */
	inline unsigned short getpixelUSHORT(unsigned x, unsigned y) {return ((unsigned short *)pixptr)[((x2-x1)<<11)+x];};
	inline float getpixel(unsigned x, unsigned y) {return ((float *)pixptr)[((x2-x1)<<11)+x];};
	
	/*! 
	 * setpixel(value, unsigned x, unsigned y)
	 * \brief get a pixel value at coordinates x,y.
	 * \param value - the pixel value
	 * \param x - the x coordinate
	 * \param y - the y coordinate
	 * \note This method does a shift left by 11 as the image width is known to be 2048 pixels
	 * \note and a shift is faster than a multiply. If the espadons image width ever changes
	 * \note this will have to be revisited.
	 * \return unsigned short
	 */
	inline void setpixel(unsigned short value, unsigned x, unsigned y) {((unsigned short *)pixptr)[((x2-x1)<<11)+x] = value;};
	inline void setpixel(float value, unsigned x, unsigned y) {((float *)pixptr)[((x2-x1)<<11)+x] = value;};
	
	/*! 
	 * operaEspadonsImageCopyHeader(operaEspadonsImage *from) 
	 * \brief Copies all of the header information from image.
	 * \param from
	 * \throws operaException cfitsio error code
	 * \return void
	 */
	void operaEspadonsImageCopyHeader(operaEspadonsImage *from);					// copy header unit
	
	/*! 
	 * ~operaEspadonsImage* operaEspadonsImage::operaEspadonsImageClone(operaEspadonsImage &imageIn)
	 * \brief Clone an espadons Image object.
	 * \param imageIn - pointer to image to clone
	 * \return operaEspadonsImage*
	 */
	operaEspadonsImage *operaEspadonsImageClone(operaEspadonsImage &imageIn);			// clone image object
	/*
	 * getters / setters
	 */
	/*! 
	 * operaFITSSubImage *getDatasecSubImage()
	 * \brief get the datasec subImage.
	 * \return subImage pointers
	 */
	operaFITSSubImage *getDatasecSubImage();
	/*! 
	 * void operaEspadonsImage::operaFITSImageSetData(operaFITSSubImage &subImage)
	 * \brief copy a subimage back in to an operaEspadonsImage datasec.
	 * \param datasecSubImage the datasec sub image
	 * \return void
	 */
	void operaFITSImageSetData(operaFITSSubImage &datasecSubImage);
	/*! 
	 * unsigned getnx()
	 * \brief get the DATASEC array length of dimension 1.
	 * \return unsigned length of axis 1
	 */
	unsigned getnx();
	/*! 
	 * unsigned getny()
	 * \brief get the DATASEC array length of dimension 2.
	 * \return unsigned length of axis 2
	 */
	unsigned getny();
	/*! 
	 * unsigned getndpixels()
	 * \brief get the image array number of pixels.
	 * \return unsigned number of pixels
	 */
	unsigned getndpixels();
	
	/*! 
	 * imtype_t getimtype()
	 * \brief Returns the image type (OBJECT, FLAT, ...).
	 * \return image type
	 */
	
	/*! 
	 * unsigned long getsize()
	 * \brief get the datasec size of an image.
	 * \return unsigned long sie of the image
	 */
	unsigned long getsize();
	
	imtype_t getimtype();
	/*! 
	 * string getimtypestring()
	 * \brief Returns the image type (OBJECT, FLAT, ...) as a string.
	 * \return image type
	 */
	string getimtypestring();
	
	/*! 
	 * detector_t getdetector()
	 * \brief Returns te detector type (EEV1, OLAPA, ...).
	 * \return detector type
	 */
	detector_t getdetector();
	/*! 
	 * string getdetectorstring()
	 * \brief Returns te detector type (EEV1, OLAPA, ...) as a string.
	 * \return detector type
	 */
	string getdetectorstring();
	
	/*! 
	 * amplifier_t getamplifier()
	 * \brief Returns te amplifier type (a, ab ...).
	 * \return amplifier type
	 */
	amplifier_t getamplifier();
	/*! 
	 * amplifier_t getamplifier()
	 * \brief Returns te amplifier type (a, ab ...)  as a string.
	 * \return amplifier type
	 */
	string getamplifierstring();
	
	/*! 
	 * instrumentmode_t getmode()
	 * \brief Returns the mode (sp1, sp2, pol ...).
	 * \return mode
	 */
	instrumentmode_t getmode();
	/*! 
	 * string getmodestring()
	 * \brief Returns the mode (sp1, sp2, pol ...) as a string.
	 * \return mode
	 */
	string getmodestring();
	
	/*! 
	 * speed_t getspeed()
	 * \brief Returns the speed (Fast, Slow, Normal ...).
	 * \return speed
	 */
	speed_t getspeed();
	/*! 
	 * string getspeedstring()
	 * \brief Returns the speed (Fast, Slow, Normal ...) as a string.
	 * \return speed
	 */
	string getspeedstring();
	
	/*! 
	 * speed_t getstokes()
	 * \brief Returns the stokes parameter (U, I, Vl ...).
	 * \return stokes
	 */
	stokes_t getstokes();
	/*! 
	 * string getstokesstring()
	 * \brief Returns the stokes parameter (U, I, Vl ...) as a string.
	 * \return stokes
	 */
	string getstokesstring();
	
	/*! 
	 * polarquad_t getpolarquad()
	 * \brief Returns the polar quad sequence number (1,2,3,4).
	 * \return polar sequence number
	 */
	polarquad_t getpolarquad();
	/*! 
	 * string getpolarquadstring()
	 * \brief Returns the polar quad sequence number (1,2,3,4) as a string.
	 * \return polar sequence number
	 */
	string getpolarquadstring();
};

#endif
