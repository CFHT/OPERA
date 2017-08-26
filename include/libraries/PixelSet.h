#ifndef PIXELSET_H
#define PIXELSET_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: PixelSet
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

/*! 
 * \sa class PixelSet
 * \brief Encapsulation of the pixel set.
 * \details Defines the set of pixels or subpixels where  
 * \details the flux information is extracted from. 
 * \file PixelSet.h
 * \ingroup libraries
 */

#include <vector>

class PixelSet {
private:
    unsigned nPixels;                   // number of subpixels
    std::vector<float> xcenter;         // subpixel x-center positions
    std::vector<float> ycenter;         // subpixel y-center positions
    std::vector<int> iIndex;            // original image pixel i-index (col) from which the subpixel was obtained
    std::vector<int> jIndex;            // original image pixel j-index (row) from which the subpixel was obtained
    std::vector<unsigned> redundancy;   // redundancy  
    std::vector<float> pixelValue;      // pixel value
    float subpixelArea;                 // area of subpixel in pixel^2, value depends on pixelation
public:
	
	/*
	 * Constructor
	 */
    PixelSet(void);
    PixelSet(float SubPixelArea);  
    PixelSet(unsigned NPixels, float SubPixelArea);
    
    /*
	 * Setters/Getters
	 */    
	unsigned getNPixels(void) const;
	float getXcenter(unsigned index) const;
    float getYcenter(unsigned index) const;
    int getiIndex(unsigned index) const;
    int getjIndex(unsigned index) const;
    unsigned getredundancy(unsigned index) const;
    float getPixelValue(unsigned index) const;
    float getSubpixelArea(void) const;
    void setXcenter(float Xcenter, unsigned index);
	void setYcenter(float Ycenter, unsigned index);
	void setiIndex(int iindex, unsigned k);
	void setjIndex(int jindex, unsigned k);
	void setredundancy(unsigned Redundancy, unsigned k);
	void setPixelValue(float PixelValue, unsigned index);
	void setSubpixelArea(float SubpixelArea);
    
	/*
	 * Methods
	 */      
    void resize(unsigned newsize);
	float getMinXcoord(void) const;
    float getMaxXcoord(void) const;
    float getMinYcoord(void) const;
    float getMaxYcoord(void) const;
};

#endif
