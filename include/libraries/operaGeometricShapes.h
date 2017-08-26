#ifndef OPERAGEOMETRICSHAPES_H
#define OPERAGEOMETRICSHAPES_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaGeometricShapes
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
 
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

class operaPoint;
class Rectangle;
class Circle;
class Polygon;
class Line;

/*! 
 * \sa class operaGeometricShapes
 * \brief geometric shapes containers.
 * \details The operaGeometricShapes defines geometric shapes that can be used  
 * \details by operaExtractionAperture in order to obtain information from 
 * \details specific regions of FITS images.
 * \return none
 * \file operaGeometricShapes.h
 * \ingroup libraries
 */

/*
 * class operaPoint
 */
class operaPoint {
private:
    float Xcoord, Ycoord;
    
public:
    operaPoint(void);
    operaPoint(float xp, float yp);    
    void setPoint(float xp, float yp);
    float getXcoord(void) const;
    float getYcoord(void) const;
    void shift(float xshift, float yshift);
    void rotate(float angle);
};

/*
 * Simple boxes for things like guide windows in WIRCam
 */
class Box {
private:
	unsigned x1, x2, y1, y2;
	unsigned xDim, yDim;
	
public:
	Box(void);
	Box(unsigned X1, unsigned X2, unsigned Y1, unsigned Y2);
	Box(unsigned X1, unsigned X2, unsigned Y1, unsigned Y2, unsigned XDim, unsigned YDim);
	unsigned getX1(void) const;
	unsigned getX2(void) const;
	unsigned getY1(void) const;
	unsigned getY2(void) const;
	unsigned getDX(void) const;
	unsigned getDY(void) const;
	unsigned getXDim(void) const;
	unsigned getYDim(void) const;
	void setX1(unsigned X1);
	void setX2(unsigned X2);
	void setY1(unsigned Y1);
	void setY2(unsigned Y2);
	void setXDim(unsigned XDim);
	void setYDim(unsigned YDim);
	unsigned getSize(void) const;
};

/*
 * Axis aligned bounding box
 */
class BoundingBox {
private:
    float minx;
    float maxx;
    float miny;
    float maxy;
    
public:
	BoundingBox();
	BoundingBox(float MinX, float MaxX, float MinY, float MaxY);
    BoundingBox(float Width, float Height, const operaPoint &Center);
    float getMinX() const;
    float getMaxX() const;
    float getMinY() const;
    float getMaxY() const;
    float getWidth(void) const;
    float getHeight(void) const;
    const operaPoint getCenter(void) const;
    void shift(float xshift, float yshift);
    bool pointInShape(const operaPoint &testPoint) const;
};

/*
 * class Rectangle
 */
class Rectangle {
private:
    float width;
    float height;
    float angle; // 0 <= angle < 180 degrees
    operaPoint center;
    operaPoint corners[FOURSIDES];
    
public:
    Rectangle(void);
    Rectangle(float Width, float Height, float Angle);
	Rectangle(float Width, float Height, float Angle, const operaPoint &Center);   
    const operaPoint& getCorner(unsigned index)  const;
    float getWidth(void) const;
    float getHeight(void) const;
    float getAngle(void) const;   
	const operaPoint& getCenter(void) const;
    void rotate(float Angle);
    void shift(float xshift, float yshift);
    void printCorners(void); 
	bool pointInShape(const operaPoint &testPoint) const;
	BoundingBox getBoundingBox() const;
};

/*
 * class Circle
 */
class Circle {
private:
    float radius;
    operaPoint center;   
    
public:
    Circle(void);     
	Circle(float Radius);
    Circle(float Radius, const operaPoint &Center);
    float getRadius(void) const;
    const operaPoint &getCenter(void) const;
    void shift(float xshift, float yshift);
    bool pointInShape(const operaPoint &testPoint) const;
    BoundingBox getBoundingBox() const;
};

/*
 * class Polygon
 */
class Polygon {
private:
    unsigned nSides;
    operaPoint vertices[MAXNPOLYGONSIDES];
    
public:
    Polygon(void);
    Polygon(unsigned NSides, const operaPoint Vertices[]);
 	unsigned getNSides(void) const;
    const operaPoint& getVertex(unsigned index) const;
    void simplePolygonization(void);
    void shift(float xshift, float yshift);
    void printVertexCoordinates(void);
    bool pointInShape(const operaPoint &testPoint) const;
    BoundingBox getBoundingBox() const;
};

/*
 * class Line
 */
typedef enum LinePosition {
	line_duplicate, line_top, line_left, line_right, line_bottom, line_perpendicular
} LinePosition_t;

class Line {
private:
    float slope;
    float intercept;    
    float width;
    float length;
    operaPoint midPoint;
    
public:
    Line(void);
    Line(float Slope);
    Line(float Slope, const operaPoint &SamplePoint);
    Line(float Slope, float Intercept);
    Line(float Slope, float Width, float Length);
    Line(float Slope, float Intercept, float Width, float Length);
	Line(float Slope, float Width, float Length, const operaPoint &MidPoint);
    float getSlope(void) const;
    float getWidth(void) const;
    float getLength(void) const;
    const operaPoint& getMidPoint(void) const;
    void setMidPoint(const operaPoint &MidPoint);    
	float getIntercept(void) const;
    void printLineEquation(void); 
    float getYcoord(float x) const;
    float getXcoord(float y) const;
    void shift(float xshift, float yshift);
    float getLineYWidth(void) const;
    float getLineXLength(void) const;
    operaPoint getIntersectionPoint(const Line &inputLine) const;
    Line getPerpendicularLine() const;
    Line getTopLine() const;
    Line getBottomLine() const;
    Line getLeftLine() const;
    Line getRightLine() const;
    bool pointOnLine(const operaPoint &TestPoint) const;
    bool pointInShape(const operaPoint &testPoint) const;
    BoundingBox getBoundingBox() const;
};

#endif
