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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaExtractionAperture.h"
#include "libraries/PixelSet.h"
#include "libraries/operaStats.h"

/*!
 * operaExtractionAperture
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the Extraction Aperture.
 * \file operaExtractionAperture.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaExtractionAperture
 * \brief The Extraction Aperture defines the set of pixels or subpixels where  
 * \brief the flux information is extracted from
 * \return none
 */

/*
 * operaExtractionAperture Constructors/Destructor
 */

template <class T>
T BoundBetween(T var, T min, T max) {
	if(var < min) return min;
	else if(var > max) return max;
	return var;
}

struct PixelBox {
public:
	int x0;
    int xn;
    int y0;
    int yn;
	PixelBox(BoundingBox box) {
		x0 = floor(box.getMinX());
		xn = ceil(box.getMaxX());
		y0 = floor(box.getMinY());
		yn = ceil(box.getMaxY());
	}
	PixelBox(BoundingBox box, unsigned naxis1, unsigned naxis2) {
		x0 = BoundBetween<int>(floor(box.getMinX()), 0, naxis1);
		xn = BoundBetween<int>(ceil(box.getMaxX()), 0, naxis1);
		y0 = BoundBetween<int>(floor(box.getMinY()), 0, naxis2);
		yn = BoundBetween<int>(ceil(box.getMaxY()), 0, naxis2);
	}
	PixelBox(BoundingBox box, operaInstrumentProfile* instrumentProfile) {
		x0 = (int)instrumentProfile->getIPixiIndex(box.getMinX());
		xn = (int)instrumentProfile->getIPixiIndex(box.getMaxX());
		if(xn < (int)instrumentProfile->getNXPoints()) xn++;
		y0 = (int)instrumentProfile->getIPixjIndex(box.getMinY());
		yn = (int)instrumentProfile->getIPixjIndex(box.getMaxY());
		if(yn < (int)instrumentProfile->getNYPoints()) yn++;
	}
	unsigned getWidth() const {
		return (xn-x0);
	}
	unsigned getHeight() const {
		return (yn-y0);
	}
	unsigned getArea() const {
		return (xn-x0)*(yn-y0);
	}
};	

template <class Shape>
operaExtractionAperture<Shape>::operaExtractionAperture(void) : xsampling(1), ysampling(1) { }

template <class Shape>
operaExtractionAperture<Shape>::operaExtractionAperture(Shape *ApertureShape, unsigned XSampling, unsigned YSampling) : xsampling(XSampling), ysampling(YSampling)
{
    setApertureShape(*ApertureShape);
    subpixels.setSubpixelArea(1.0/(float)(xsampling*ysampling));
    setSubpixels();
}

template <class Shape>
operaExtractionAperture<Shape>::operaExtractionAperture(Shape *ApertureShape, unsigned XSampling, unsigned YSampling, operaFITSImage &Image) : xsampling(XSampling), ysampling(YSampling)
{
    setApertureShape(*ApertureShape);
    subpixels.setSubpixelArea(1.0/(float)(xsampling*ysampling));
    setSubpixels(Image);
}

template <class Shape>
operaExtractionAperture<Shape>::operaExtractionAperture(Shape *ApertureShape, operaInstrumentProfile *instrumentProfile, float distd)
{
    setApertureShape(*ApertureShape);
    xsampling = instrumentProfile->getXsampling();
    ysampling = instrumentProfile->getYsampling();
    subpixels.setSubpixelArea(1.0/(float)(xsampling*ysampling));
    setSubpixels(instrumentProfile,distd);
}

template <class Shape>
operaExtractionAperture<Shape>::operaExtractionAperture(Shape *ApertureShape, operaInstrumentProfile *instrumentProfile)
{
    setApertureShape(*ApertureShape);
    xsampling = instrumentProfile->getXsampling();
    ysampling = instrumentProfile->getYsampling();
    subpixels.setSubpixelArea(1.0/(float)(xsampling*ysampling));
    setSubpixels(instrumentProfile);
}

/*
 * operaExtractionAperture Setters/Getters
 */

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixels(const PixelSet& Subpixels) {
    subpixels = Subpixels;
}

template <class Shape>
const PixelSet* operaExtractionAperture<Shape>::getSubpixels(void) const {
    return &subpixels;
}

template <class Shape>
void operaExtractionAperture<Shape>::setApertureShape(const Shape& ApertureShape) {
	apertureShape = ApertureShape;
	boundingBox = apertureShape.getBoundingBox();
}

template <class Shape>
const Shape* operaExtractionAperture<Shape>::getApertureShape(void) const {
    return &apertureShape;
}

template <class Shape>
void operaExtractionAperture<Shape>::setSampling(unsigned Xsampling, unsigned Ysampling) {
    xsampling = Xsampling;
    ysampling = Ysampling;    
}

template <class Shape>
unsigned operaExtractionAperture<Shape>::getXsampling(void) const {
    return xsampling;
}

template <class Shape>
unsigned operaExtractionAperture<Shape>::getYsampling(void) const {
    return ysampling;
}

template <class Shape>
float operaExtractionAperture<Shape>::getFluxFraction(void) const {
    return fluxFraction;
}

template <class Shape>
void operaExtractionAperture<Shape>::setFluxFraction(float FluxFraction) {
    fluxFraction = FluxFraction;
}

/*
 * operaExtractionAperture Methods
 */

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixelPositions(PixelBox pixelrange) {
    unsigned nPixels = 0;
	for(int i=pixelrange.x0; i<pixelrange.xn; i++){
		for(unsigned ii=0;ii<xsampling;ii++) {
			float xsubpix = (float)i + (float)(0.5+ii)/float(xsampling);
			for(int j=pixelrange.y0; j<pixelrange.yn; j++){
				for(unsigned jj=0;jj<ysampling;jj++) {
					float ysubpix = (float)j + (float)(0.5+jj)/float(ysampling);
					operaPoint TestPoint(xsubpix,ysubpix);
					if(apertureShape.pointInShape(TestPoint)) {
						subpixels.setXcenter(xsubpix, nPixels);
						subpixels.setYcenter(ysubpix, nPixels);
						subpixels.setiIndex(i, nPixels);
						subpixels.setjIndex(j, nPixels);
						subpixels.setredundancy(1,nPixels);
						nPixels++;
					}
				}
			}
		}
	}
	subpixels.resize(nPixels);
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixelPositionsWithRedundancy(PixelBox pixelrange) {
    unsigned nPixels = 0;
	for(int i=pixelrange.x0; i< pixelrange.xn; i++){
		for(int j=pixelrange.y0; j<pixelrange.yn; j++){
			unsigned redundancy = 0;
			for(unsigned ii=0;ii<xsampling;ii++) {
				float xsubpix = (float)i + (float)(0.5+ii)/float(xsampling);
				for(unsigned jj=0;jj<ysampling;jj++) {
					float ysubpix = (float)j + (float)(0.5+jj)/float(ysampling);
					operaPoint TestPoint(xsubpix,ysubpix);
					if(apertureShape.pointInShape(TestPoint)) redundancy++;
				}
			}
			if(redundancy) {
				subpixels.setXcenter((float)i+0.5, nPixels);
				subpixels.setYcenter((float)j+0.5, nPixels);
				subpixels.setiIndex(i, nPixels);
				subpixels.setjIndex(j, nPixels);
				subpixels.setredundancy(redundancy,nPixels);
				nPixels++;
			}
		}
	}
	subpixels.resize(nPixels);
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixelPositions(PixelBox pixelrange, operaInstrumentProfile *instrumentProfile) {
    unsigned nPixels = 0;
    for(int i=pixelrange.x0;i<pixelrange.xn;i++){
		float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
		for(int j=pixelrange.y0;j<pixelrange.yn;j++){
			float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                           
			operaPoint TestPoint(xsubpix,ysubpix);
			if(apertureShape.pointInShape(TestPoint)) {
				subpixels.setXcenter(xsubpix, nPixels);
				subpixels.setYcenter(ysubpix, nPixels);
				subpixels.setiIndex(i, nPixels);
				subpixels.setjIndex(j, nPixels);
				subpixels.setredundancy(1,nPixels);
				nPixels++;
			}
		}
	}
	subpixels.resize(nPixels);
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixelValues(operaFITSImage &Image) {
    unsigned nPixels = subpixels.getNPixels();
	for(unsigned pix=0; pix<nPixels; pix++) {
		int i = subpixels.getiIndex(pix);
		int j = subpixels.getjIndex(pix);
		subpixels.setPixelValue((float)Image[j][i], pix);
	}
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixelValues(operaInstrumentProfile *instrumentProfile, float d) {
	unsigned nPixels = subpixels.getNPixels();
	for(unsigned pix=0; pix<nPixels; pix++) {
		int i = subpixels.getiIndex(pix);
		int j = subpixels.getjIndex(pix);
		subpixels.setPixelValue(instrumentProfile->getipDataFromPolyModel(d,(unsigned)i,(unsigned)j), pix);
	}
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixels(bool withRedundancy) {
    PixelBox pixelrange(boundingBox);
    subpixels.resize(pixelrange.getArea()*xsampling*ysampling);
	if(withRedundancy) setSubpixelPositionsWithRedundancy(pixelrange);
	else setSubpixelPositions(pixelrange);
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixels(int naxis1, int naxis2, bool withRedundancy) {
    PixelBox pixelrange(boundingBox, naxis1, naxis2);
    subpixels.resize(pixelrange.getArea()*xsampling*ysampling);
	if(withRedundancy) setSubpixelPositionsWithRedundancy(pixelrange);
	else setSubpixelPositions(pixelrange);
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixels(operaInstrumentProfile *instrumentProfile) {
    PixelBox pixelrange(boundingBox, instrumentProfile);
    subpixels.resize(pixelrange.getArea()*xsampling*ysampling);
    setSubpixelPositions(pixelrange, instrumentProfile);
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixels(operaFITSImage &Image) {
    setSubpixels(Image.getnaxis1(), Image.getnaxis2());
    setSubpixelValues(Image);
}

template <class Shape>
void operaExtractionAperture<Shape>::setSubpixels(operaInstrumentProfile *instrumentProfile, float d) {
    setSubpixels(instrumentProfile);
	setSubpixelValues(instrumentProfile, d);
}

template <class Shape>
void operaExtractionAperture<Shape>::shift(float xshift, float yshift) {
    apertureShape.shift(xshift, yshift);
    boundingBox.shift(xshift, yshift);
}

template <class Shape>
void operaExtractionAperture<Shape>::recenter(const operaPoint &NewCenter) {
    float xshift = boundingBox.getCenter().getXcoord() - NewCenter.getXcoord();
    float yshift = boundingBox.getCenter().getYcoord() - NewCenter.getYcoord();
    shift(xshift, yshift);
}

template <class Shape>
void operaExtractionAperture<Shape>::shiftAperture(float xshift, float yshift) {
    shift(xshift,yshift);
    setSubpixels();
}

template <class Shape>
void operaExtractionAperture<Shape>::shiftAperture(float xshift, float yshift, operaFITSImage &Image) {
    shift(xshift,yshift);
    setSubpixels(Image);
}

template <class Shape>
void operaExtractionAperture<Shape>::recenterAperture(const operaPoint &NewCenter, bool withRedundancy) {
    recenter(NewCenter);
    setSubpixels(withRedundancy);
}

template <class Shape>
void operaExtractionAperture<Shape>::recenterAperture(const operaPoint &NewCenter, operaFITSImage &Image) {
    recenter(NewCenter);
    setSubpixels(Image);
}

template <class Shape>
void operaExtractionAperture<Shape>::recenterAperture(const operaPoint &NewCenter, int naxis1, int naxis2, bool withRedundancy) {
    recenter(NewCenter);
    setSubpixels(naxis1,naxis2,withRedundancy);
}

template class operaExtractionAperture<Line>;
