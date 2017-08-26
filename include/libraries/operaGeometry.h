/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaGeometry
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

#ifndef OPERAGEOMETRY_H
#define OPERAGEOMETRY_H

#include "libraries/Polynomial.h"	

/*! 
 * \sa class operaGeometry
 * \brief Encapsulation of the Geom object.
 * \return none
 * \file operaGeometry.h
 * \ingroup libraries
 */

class operaSpectralOrder;

enum axis_t {Unknown_axis, cols, rows};
enum dispersiondirection_t {Unknown_direction, up, down};

class operaGeometry {
	
private:
	
	unsigned Ndatapoints; 
	doubleValueVector_t OrderCenters;			// vectors of ndatapoints long of x,y,values, errors along the center of each order
	doubleValueVector_t pixelValueVector;		// vectors of coords, values, errors
	Polynomial *geometryPolynomial;				// Polynomial describing the order geometry
	double ymin;
	double ymax;
	double centerX;
	double centerY;
	double apertureWidth; 	
	double orderLength;
	unsigned NumberofPointsToBinInYDirection;
	axis_t dispersionAxis;						// (axis_t cols rows)
	dispersiondirection_t dispersionDirection;	// (dispersiondirection_t up = wl-increase up, down = wl-increase down)  
	
	double orderSeparation;						// peak-to-peak separation between current and previous order (in pixel units)
	
public:
	
	/*
	 * Constructor
	 */
	operaGeometry();
	operaGeometry(unsigned maxdatapoints, unsigned maxValues);
	
	/*
	 * Destructor
	 */
	~operaGeometry();
	
	/*
	 * Methods
	 */
	
	double CalculateAndSetOrderLength(void);	
	
	double getOrderLength(void) const;
	
	double CalculateDistance(double yMin, double yMax) const;
	
	void setNumberofPointsToBinInYDirection(unsigned Points);
	
	unsigned getNumberofPointsToBinInYDirection() const;
	
	void traceOrder(unsigned coeffs, double &chisqr, bool witherrors);
	
	Polynomial *getCenterPolynomial();
	
	const Polynomial *getCenterPolynomial() const;
	
	void setNdatapoints(unsigned dp);
	
	unsigned getNdatapoints(void) const;
	
	double getYmin(void) const;
	
	void setYmin(double min);
	
	double getCenterX(unsigned index) const;
	
	double getCenterY(unsigned index) const;
	
	double getCenterV(unsigned index) const;
	
	void setYmax(double max);
	
	double getYmax(void) const;
	
	double getapertureWidth(void) const;
	
	void setapertureWidth(double width);
	
	axis_t getdispersionAxis(void) const;
	
	void setdispersionAxis(axis_t axis);
	
	dispersiondirection_t getdispersionDirection(void) const;
	
	void setdispersionDirection(dispersiondirection_t direction);
	
	void createPixelValueVector(unsigned NValues);
	
	void deletePixelValueVector(void);
	
	void addPixelValue(double x, double y, double value);	
	
	void addPixelValue(double x, double y, double value, double error);	
	
	void createOrderCentersVector(unsigned NValues);
	
	void deleteOrderCentersVector(void);
	
	void addOrderCenterValue(double x, double y, double value);
	
	void addOrderCenterValue(double x, double y, double value, double error);
	
	doubleValueVector_t getOrderCenters(void);
	
	doubleValueVector_t getPixelValueVector(void);
	
	void setorderSeparation(double OrderSeparation);	
	
	double getorderSeparation(void) const; 
    
    double CalculateMinimumYBinSize(double y) const;
    
    void resetCenter(double x, double y, double xerror, unsigned index);
	
};
#endif
