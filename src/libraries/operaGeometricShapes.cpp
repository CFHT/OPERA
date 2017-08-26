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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaGeometricShapes.h"
#include "libraries/operaException.h"

#include "libraries/operaLibCommon.h"
#include "libraries/operaStats.h"   // for operaArrayIndexSort

/*!
 * operaGeometricShapes
 * \author Doug Teeple / Eder Martioli
 * \brief GeometricShapes  
 * \details The operaGeometricShapes defines geometric shapes that can be used  
 * \details by operaExtractionAperture in order to obtain information from 
 * \details specific regions of FITS images.
 * \file operaGeometricShapes.cpp
 * \ingroup libraries
 */

using namespace std;

/*
 * operaPoint class 
 */

operaPoint::operaPoint(void) : Xcoord(0), Ycoord(0) { }

operaPoint::operaPoint(float xp, float yp) : Xcoord(xp), Ycoord(yp) { }

void operaPoint::setPoint(float xp, float yp) {
    Xcoord = xp;
    Ycoord = yp;
}

float operaPoint::getXcoord(void) const {
    return Xcoord;
}

float operaPoint::getYcoord(void) const {
    return Ycoord;
}

void operaPoint::shift(float xshift, float yshift) {
	Xcoord += xshift;
	Ycoord += yshift;
}

void operaPoint::rotate(float angle) {
	float angleInRadians = angle*M_PI/180.0;
    float s = sin(angleInRadians);
    float c = cos(angleInRadians);
    setPoint(Xcoord * c - Ycoord * s, Xcoord * s + Ycoord * c);
}

/*
 * Simple boxes for things like guide windows in WIRCam
 */

Box::Box(void) : x1(0), x2(0), y1(0), y2(0), xDim(0), yDim(0) { }

Box::Box(unsigned X1, unsigned X2, unsigned Y1, unsigned Y2) : x1(X1), x2(X2), y1(Y1), y2(Y2), xDim(0), yDim(0) { }

Box::Box(unsigned X1, unsigned X2, unsigned Y1, unsigned Y2, unsigned XDim, unsigned YDim) : x1(X1), x2(X2), y1(Y1), y2(Y2), xDim(XDim), yDim(YDim) { }

unsigned Box::getX1(void) const {
	return x1;
}

unsigned Box::getX2(void) const {
	return x2;
}

unsigned Box::getY1(void) const {
	return y1;
}

unsigned Box::getY2(void) const {
	return y2;
}

unsigned Box::getDX(void) const {
	return x2-x1;
}

unsigned Box::getDY(void) const {
	return y2-y1;
}

void Box::setX1(unsigned X1) {
	x1 = X1;
}

void Box::setX2(unsigned X2) {
	x2 = X2;
}

void Box::setY1(unsigned Y1) {
	y1 = Y1;
}

void Box::setY2(unsigned Y2) {
	y2 = Y2;
}

unsigned Box::getSize(void) const {
	return (y2-y1)*(x2-x1);
}

/*
 * BoundingBox class
 */

BoundingBox::BoundingBox() : minx(0), maxx(0), miny(0), maxy(0) { }

BoundingBox::BoundingBox(float MinX, float MaxX, float MinY, float MaxY) : minx(MinX), maxx(MaxX), miny(MinY), maxy(MaxY) { }

BoundingBox::BoundingBox(float Width, float Height, const operaPoint &Center) {
	minx = Center.getXcoord() - Width/2;
	maxx = Center.getXcoord() + Width/2;
	miny = Center.getYcoord() - Height/2;
	maxy = Center.getYcoord() + Height/2;
}

float BoundingBox::getMinX() const {
	return minx;
}

float BoundingBox::getMaxX() const {
	return maxx;
}

float BoundingBox::getMinY() const {
	return miny;
}

float BoundingBox::getMaxY() const {
	return maxy;
}

float BoundingBox::getWidth(void) const {
	return maxx-minx;
}

float BoundingBox::getHeight(void) const {
	return maxy-miny;
}

const operaPoint BoundingBox::getCenter(void) const {
	return operaPoint((minx+maxx)/2, (miny+maxy)/2);
}

void BoundingBox::shift(float xshift, float yshift) {
	minx += xshift;
	maxx += xshift;
	miny += yshift;
	maxy += yshift;
}

bool BoundingBox::pointInShape(const operaPoint &testPoint) const {
	return testPoint.getXcoord() >= minx && testPoint.getXcoord() <= maxx && testPoint.getYcoord() >= miny && testPoint.getYcoord() <= maxy;
}

/*
 * Rectangle class 
 */

Rectangle::Rectangle(void) : width(0), height(0), angle(0), center(0,0) { }

Rectangle::Rectangle(float Width, float Height, float Angle) : width(Width), height(Height), angle(Angle), center(0,0) {
    if(Width <= 0 || Height <= 0){
        throw operaException("Rectangle: error: sides must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }   
    
    corners[0].setPoint(-width/2,-height/2);
    corners[1].setPoint(+width/2,-height/2);    
    corners[2].setPoint(+width/2,+height/2);    
    corners[3].setPoint(-width/2,+height/2);    
    
    rotate(Angle);
}

Rectangle::Rectangle(float Width, float Height, float Angle, const operaPoint &Center) : width(Width), height(Height), angle(Angle), center(0,0) {   
    if(Width <= 0 || Height <= 0){
        throw operaException("Rectangle: error: sides must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    corners[0].setPoint(-width/2,-height/2);
    corners[1].setPoint(+width/2,-height/2);
    corners[2].setPoint(+width/2,+height/2);
    corners[3].setPoint(-width/2,+height/2);
    
    rotate(Angle);
    shift(Center.getXcoord(), Center.getYcoord());
}

const operaPoint& Rectangle::getCorner(unsigned index) const {
    return corners[index];
}

float Rectangle::getWidth(void) const {
	return width;   
}

float Rectangle::getHeight(void) const {
    return height;    
}

float Rectangle::getAngle(void) const {
    return angle;    
}

const operaPoint& Rectangle::getCenter(void) const {
    return center;
}

void Rectangle::rotate(float Angle) {
	for(unsigned i=0;i<FOURSIDES;i++) {
        corners[i].shift(-center.getXcoord(), -center.getYcoord());
        corners[i].rotate(angle);
        corners[i].shift(center.getXcoord(), center.getYcoord());
    }
}

void Rectangle::shift(float xshift, float yshift) {
	center.shift(xshift, yshift);
	for(unsigned i=0;i<FOURSIDES;i++) {
        corners[i].shift(xshift, yshift);
    }
}

void Rectangle::printCorners(void) {
    for (unsigned i=0; i<FOURSIDES; i++) {
        cout << "(" << corners[i].getXcoord() << "," << corners[i].getYcoord() << ")\t";
    }
    cout << endl;
}

bool Rectangle::pointInShape(const operaPoint &TestPoint) const {
	float lineSlope = 0;
    if(angle != 90 && angle != -90 && angle != 270 && angle != -270) {
        lineSlope = tan(angle*M_PI/180.0);
    }
	Line rectangleCenterLine(lineSlope, height, width, center);
    return rectangleCenterLine.pointInShape(TestPoint);
}

BoundingBox Rectangle::getBoundingBox(void) const {
	float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    for (unsigned i=0; i<FOURSIDES; i++) {
        if(corners[i].getXcoord() < minx) minx = corners[i].getXcoord();
        if(corners[i].getXcoord() > maxx) maxx = corners[i].getXcoord();
        if(corners[i].getYcoord() < miny) miny = corners[i].getYcoord();
        if(corners[i].getYcoord() > maxy) maxy = corners[i].getYcoord();
	}
    return BoundingBox(minx, maxx, miny, maxy);
}

/*
 * Circle class 
 */

Circle::Circle(void) : radius(0), center(0,0) { }

Circle::Circle(float Radius) : radius(Radius), center(0,0) {
    if(Radius <= 0){
        throw operaException("Circle: error: radius must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }
}

Circle::Circle(float Radius, const operaPoint &Center) : radius(Radius), center(Center)  {
    if(Radius <= 0){
        throw operaException("Circle: error: radius must be greater than zero",operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }
}

float Circle::getRadius(void) const {
    return radius;
}

const operaPoint& Circle::getCenter(void) const {
    return center;
}

void Circle::shift(float xshift, float yshift) {
	center.shift(xshift, yshift);
}

bool Circle::pointInShape(const operaPoint &TestPoint) const {
	float testX = TestPoint.getXcoord() - center.getXcoord();
    float testY = TestPoint.getYcoord() - center.getYcoord();
	return testX*testX + testY*testY < radius;
}

BoundingBox Circle::getBoundingBox(void) const {
	return BoundingBox(2*radius, 2*radius, center);
}

/*
 * Polygon class 
 */

Polygon::Polygon(void) : nSides(0) { }

Polygon::Polygon(unsigned NSides, const operaPoint Vertices[]) : nSides(NSides) {
    if(NSides < 3){
        throw operaException("Polygon: error: number of sides must be 3 or greater", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    for(unsigned i=0; i<nSides; i++) {
        vertices[i] = Vertices[i];
    }
    
    simplePolygonization(); // organize vertices to make polygon simple, i.e. with no intersecting lines     
}

unsigned Polygon::getNSides(void) const {
    return nSides;     
}

const operaPoint& Polygon::getVertex(unsigned index) const {
    return vertices[index];
}

void Polygon::simplePolygonization(void) {
    unsigned anchorPoint=0;
    float minYCoord = 1e30;
    
    // First determine the anchor point as the point with minimum Y coordinate value. 
    for(unsigned i=0;i<nSides;i++){
        if(vertices[i].getYcoord() < minYCoord) {
            minYCoord = vertices[i].getYcoord();
            anchorPoint = i;
        }
    }
	
    float polarAngle[MAXNPOLYGONSIDES];
    
    unsigned np=0;
    for(unsigned i=0;i<nSides;i++){
        if(i != anchorPoint) {
            float x = vertices[i].getXcoord() - vertices[anchorPoint].getXcoord();
            float y = vertices[i].getYcoord() - vertices[anchorPoint].getYcoord();
            polarAngle[np] = atan2(y,x);
			
#ifdef PRINT_DEBUG    
            cerr << "Polygon::simplePolygonization: anchorPoint = " << anchorPoint << endl;
            cerr << "Polygon::simplePolygonization: i = "<< i << " vertexXcoords,vertexYcoords = " << vertices[i]->getXcoord() << "," << vertices[i]->getYcoord() << endl;
            cerr << "Polygon::simplePolygonization: x,y = " << x << "," << y << endl;
            cerr << "Polygon::simplePolygonization: np = " << np << " polarAngle = " <<  polarAngle[np] << endl;
#endif            
            np++;
        }
    }
    int sindex[MAXNPOLYGONSIDES];
    operaArrayIndexSort((int)np,polarAngle,sindex);
    
    operaPoint NewVertices[MAXNPOLYGONSIDES]; 
    
    NewVertices[0].setPoint(vertices[anchorPoint].getXcoord(), vertices[anchorPoint].getYcoord());
    
    np=0;
    for(unsigned i=0;i<nSides;i++) {
        if(i != anchorPoint) {
            NewVertices[np+1].setPoint(vertices[sindex[np]].getXcoord(), vertices[sindex[np]].getYcoord());
            np++;
        }
    }
    
    for(unsigned i=0;i<nSides;i++) {
        vertices[i] = NewVertices[i];
    }
}

void Polygon::shift(float xshift, float yshift) {
	for(unsigned i=0;i<nSides;i++) {
        vertices[i].shift(xshift, yshift);
    }
}

void Polygon::printVertexCoordinates(void) {
    for (unsigned i=0; i<nSides; i++) {
        cout << "(" << vertices[i].getXcoord() << "," << vertices[i].getYcoord() << ")\t";
    }
    cout << endl;
}

/* 
 * Polygon::pointInPolygon(operaPoint testPoint)
 * \author Doug Teeple / Eder Martioli
 * \brief This function will return TRUE if the test point is inside the polygon, or
 * \brief FALSE if it is not.  If the point is exactly on the edge of the polygon,
 * \brief then the function may return TRUE or FALSE.
 * \Notes Source: http://alienryderflex.com/polygon/
 * \ingroup libraries
 */

bool Polygon::pointInShape(const operaPoint &testPoint) const {
    
    // Comment: Eder, April 10 2012. 
    // This function still doesn't work properly
	// Eder addedparens for ambiguity - please check that this is what you want...
	// DT changed ^= (bitwise XOR) to |= (or equals)
    
    unsigned j=nSides-1 ;
    bool oddNodes = false;
    
    float testXcoord = testPoint.getXcoord();
    float testYcoord = testPoint.getYcoord();
    
    for (unsigned i=0; i<nSides; i++) {
        if (((vertices[i].getYcoord() < testYcoord && vertices[j].getYcoord()>=testYcoord)
             ||   (vertices[j].getYcoord() < testYcoord && vertices[i].getYcoord()>=testYcoord))
            &&  (vertices[i].getXcoord()<=testXcoord || vertices[j].getXcoord()<=testXcoord)) {
            
            oddNodes |= (vertices[i].getXcoord()+(testYcoord-vertices[i].getYcoord())/(vertices[j].getYcoord()-vertices[i].getYcoord())*(vertices[j].getXcoord()-vertices[i].getXcoord())) < testXcoord;
        }
        j=i;
    }
    
    return oddNodes; 
}

BoundingBox Polygon::getBoundingBox(void) const {
	float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    for (unsigned i=0; i<nSides; i++) {
        if(vertices[i].getXcoord() < minx) minx = vertices[i].getXcoord();
        if(vertices[i].getXcoord() > maxx) maxx = vertices[i].getXcoord();
        if(vertices[i].getYcoord() < miny) miny = vertices[i].getYcoord();
        if(vertices[i].getYcoord() > maxy) maxy = vertices[i].getYcoord();
    }
    operaPoint center((maxx+minx)/2 , (maxy+miny)/2);
    return BoundingBox(minx, maxx, miny, maxy);
}

/*
 * Line class 
 */

Line::Line(void) : slope(0), intercept(0), width(0), length(BIG), midPoint(0,0) { }

Line::Line(float Slope) : slope(Slope), intercept(0), width(0), length(BIG), midPoint(0,0) { }

Line::Line(float Slope, float Intercept) : slope(Slope), intercept(Intercept), width(0), length(BIG), midPoint(0,Intercept) { }

Line::Line(float Slope, const operaPoint &SamplePoint) : slope(Slope), intercept(0), width(0), length(BIG), midPoint(SamplePoint) {
    intercept = midPoint.getYcoord() - slope*midPoint.getXcoord();
}

Line::Line(float Slope, float Width, float Length) : slope(Slope), intercept(0), width(Width), length(Length), midPoint(0,0) {
    if(Width < 0 || Length <= 0){
        throw operaException("Line: error: line dimensions must be greater than zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }
}

Line::Line(float Slope, float Intercept, float Width, float Length) : slope(Slope), intercept(Intercept), width(Width), length(Length), midPoint(0,intercept) {
    if(Width < 0 || Length <= 0){
        throw operaException("Line: error: line dimensions must be greater than zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }
}

Line::Line(float Slope, float Width, float Length, const operaPoint &MidPoint) : slope(Slope), intercept(0), width(Width), length(Length), midPoint(MidPoint) {
    if(Width < 0 || Length <= 0){
        throw operaException("Line: error: line dimensions must be greater than zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
    }          
    intercept = midPoint.getYcoord() - slope*midPoint.getXcoord();
}

float Line::getSlope(void) const {
	return slope;     
}

float Line::getWidth(void) const {
	return width;     
}

float Line::getLength(void) const {
	return length;     
}

const operaPoint& Line::getMidPoint(void) const {
	return midPoint;     
}

void Line::setMidPoint(const operaPoint &MidPoint) {
    midPoint = MidPoint;
}

float Line::getIntercept(void) const {
    return intercept;
}

void Line::printLineEquation(void) {
    cout << "y = " << slope << "*x + " << intercept << endl;
}

float Line::getYcoord(float x) const {
    return slope*x + intercept;
}

float Line::getXcoord(float y) const {
    if(slope == 0) {
        throw operaException("Line: error: x is undetermined because slope is zero", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);        
    }
    return (y - intercept)/slope;
}

void Line::shift(float xshift, float yshift) {
	midPoint.shift(xshift, yshift);
	intercept = midPoint.getYcoord() - slope*midPoint.getXcoord();
}

float Line::getLineYWidth(void) const {
    return width/cos(atan(slope));
}

float Line::getLineXLength(void) const {
    return length*cos(atan(slope));
}

operaPoint Line::getIntersectionPoint(const Line &inputLine) const {
    if(inputLine.slope == slope) {
        throw operaException("Line: error: no intersection point; input line is parallel. ", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);                
    }
    float xintersect = (inputLine.intercept - intercept)/(slope - inputLine.slope);
    return operaPoint(xintersect,getYcoord(xintersect));
}

Line Line::getPerpendicularLine() const {
    float perpendicularSlope = -(1.0/slope);
    return Line(perpendicularSlope, width, length, midPoint);
}

Line Line::getTopLine() const {
    float ShiftInYDirection = getLineYWidth()/2.0;
    Line topline(slope, intercept+ShiftInYDirection, width, length);
    Line perpendicularline = getPerpendicularLine();
    operaPoint topMidPoint = topline.getIntersectionPoint(perpendicularline);
    topline.setMidPoint(topMidPoint);
	return topline;
}

Line Line::getBottomLine() const {
    float ShiftInYDirection = getLineYWidth()/2.0;
    Line bottomline (slope, intercept-ShiftInYDirection, width, length);
    Line perpendicularline = getPerpendicularLine();
    operaPoint bottomMidPoint = bottomline.getIntersectionPoint(perpendicularline);
    bottomline.setMidPoint(bottomMidPoint);
	return bottomline;
}

Line Line::getLeftLine() const {
    float midPointXshift = getLineXLength()/2.0;
    operaPoint leftMidPoint(midPoint.getXcoord() - midPointXshift,getYcoord(midPoint.getXcoord() - midPointXshift));
    return Line(getPerpendicularLine().slope, length, width, leftMidPoint);
}

Line Line::getRightLine() const {
    float midPointXshift = getLineXLength()/2.0;
	operaPoint rightMidPoint(midPoint.getXcoord() + midPointXshift,getYcoord(midPoint.getXcoord() + midPointXshift));
	return Line(getPerpendicularLine().slope, length, width, rightMidPoint);
}

bool Line::pointOnLine(const operaPoint &TestPoint) const {
    float ShiftInXDirection = getLineXLength()/2.0;
	float xmin = midPoint.getXcoord() - ShiftInXDirection;
    float xmax = midPoint.getXcoord() + ShiftInXDirection;
	return TestPoint.getXcoord() >= xmin && TestPoint.getXcoord() <= xmax && TestPoint.getYcoord() == getYcoord(TestPoint.getXcoord());
}

bool Line::pointInShape(const operaPoint &TestPoint) const {
    return TestPoint.getYcoord() < getTopLine().getYcoord(TestPoint.getXcoord()) &&
       TestPoint.getYcoord() >= getBottomLine().getYcoord(TestPoint.getXcoord()) &&
       TestPoint.getXcoord() >= getLeftLine().getXcoord(TestPoint.getYcoord()) &&
       TestPoint.getXcoord() < getRightLine().getXcoord(TestPoint.getYcoord());
}

BoundingBox Line::getBoundingBox(void) const {
	float maxx = midPoint.getXcoord() + length/2;
    float minx = midPoint.getXcoord() - length/2;
    float maxy = midPoint.getYcoord() + width/2;
    float miny = midPoint.getYcoord() - width/2;
    Line topline = getTopLine();
    Line bottomline = getBottomLine();
    Line leftline = getLeftLine();
    Line rightline = getRightLine();
    operaPoint Corner1 = rightline.getIntersectionPoint(bottomline);
    operaPoint Corner2 = topline.getIntersectionPoint(rightline);
    operaPoint Corner3 = leftline.getIntersectionPoint(topline);
    operaPoint Corner4 = bottomline.getIntersectionPoint(leftline);
    if(slope > 0) {
        maxx = Corner1.getXcoord();
        maxy = Corner2.getYcoord();
        minx = Corner3.getXcoord();
        miny = Corner4.getYcoord();
    } else if (slope < 0) {
        maxx = Corner2.getXcoord();
        maxy = Corner3.getYcoord();
        minx = Corner4.getXcoord();
        miny = Corner1.getYcoord();        
    }
    return BoundingBox(minx, maxx, miny, maxy);
}
