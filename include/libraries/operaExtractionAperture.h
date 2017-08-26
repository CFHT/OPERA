#ifndef OPERAEXTRACTIONAPERTURE_H
#define OPERAEXTRACTIONAPERTURE_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaExtractionAperture
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

#include <ostream>
#include <math.h>

#include "libraries/operaLibCommon.h"			// for MAXNPOLYGONSIDES
#include "libraries/operaStats.h"				// for operaArrayIndexSort

#include "libraries/operaGeometricShapes.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaInstrumentProfile.h"
#include "libraries/PixelSet.h"

/*! 
 * \sa class operaExtractionAperture
 * \brief The Extraction Aperture.
 * \details The Extraction Aperture defines the set of pixels or subpixels where  
 * \details the flux information is extracted from. The aperture shape can be defined 
 * \details as one of the available options.
 * \return none
 * \file operaExtractionAperture.h
 * \ingroup libraries
 */

class PixelBox;

template <class Shape>
class operaExtractionAperture {
	
private:    
    Shape apertureShape;    // operaGeometricShape
    BoundingBox boundingBox;  // Box size that englobes aperture shape
    unsigned xsampling;		// sampling factor for the extraction aperture in x direction
	unsigned ysampling;		// sampling factor for the extraction aperture in y direction
    PixelSet subpixels;     // set of subpixels with coordinates, values, and area of each subpixel
    float fluxFraction;     // fraction of flux (with respect to the entire IP) contained in aperture  
    
    // Helper functions
    void setSubpixelPositions(PixelBox pixelrange);
    void setSubpixelPositionsWithRedundancy(PixelBox pixelrange);
    void setSubpixelPositions(PixelBox pixelrange, operaInstrumentProfile *instrumentProfile);
    void setSubpixelValues(operaFITSImage &Image);
    void setSubpixelValues(operaInstrumentProfile *instrumentProfile, float d);
    void setSubpixels(bool withRedundancy = false);
    void setSubpixels(operaFITSImage &Image);
    void setSubpixels(int naxis1, int naxis2, bool withRedundancy = false);
    void setSubpixels(operaInstrumentProfile *instrumentProfile, float d); 
    void setSubpixels(operaInstrumentProfile *instrumentProfile);
    
    void shift(float xshift, float yshift);
    void recenter(const operaPoint &NewCenter);
    
public:
	
	/*
	 * Constructor
	 */
	
    operaExtractionAperture(void);
    
    operaExtractionAperture(Shape *ApertureShape, unsigned XSampling, unsigned YSampling);
    
    operaExtractionAperture(Shape *ApertureShape, unsigned XSampling, unsigned YSampling, operaFITSImage &Image);
    
    operaExtractionAperture(Shape *ApertureShape, operaInstrumentProfile *instrumentProfile, float distd);
    
    operaExtractionAperture(Shape *ApertureShape, operaInstrumentProfile *instrumentProfile);
	
    /*
	 * Setters/Getters
	 */      
    
	void setSubpixels(const PixelSet& Subpixels); 
    
    const PixelSet* getSubpixels(void) const;
    
    void setApertureShape(const Shape& ApertureShape);
    
    const Shape* getApertureShape(void) const;
    
    void setSampling(unsigned Xsampling, unsigned Ysampling);  
    
    unsigned getXsampling(void) const;
    
    unsigned getYsampling(void) const;

    float getFluxFraction(void) const;
    
    void setFluxFraction(float FluxFraction);
    
	/*
	 * Methods
	 */    

    void shiftAperture(float xshift, float yshift); 
    
    void shiftAperture(float xshift, float yshift, operaFITSImage &Image);
    
    void recenterAperture(const operaPoint &NewCenter, bool withRedundancy = false);
    
    void recenterAperture(const operaPoint &NewCenter, int naxis1, int naxis2, bool withRedundancy = false);
    
    void recenterAperture(const operaPoint &NewCenter, operaFITSImage &Image); 
};

#endif
