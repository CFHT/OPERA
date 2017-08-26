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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaGeometry.h"
#include "libraries/operaFit.h"
#include "libraries/operaMath.h"

/*!
 * \brief operaGeometry
 * \details Encapsulation of the Geom object..
 * \author Doug Teeple
 * \file operaGeometry.cpp
 * \ingroup libraries
 */

using namespace std;

/*
 * Constructor
 */
operaGeometry::operaGeometry() :
Ndatapoints(0), 
ymin(0.0),
ymax(0.0),
apertureWidth(0.0), 
orderLength(0.0),	
NumberofPointsToBinInYDirection(0),
dispersionAxis(cols),			// (axis_t 0=columns 1=rows)
dispersionDirection(up)		// (dispersondirection_t up=wl-increase up, down wl-increase down)
{
	geometryPolynomial = new Polynomial();
	OrderCenters.vectorlength = 0;
	pixelValueVector.vectorlength = 0;
}

operaGeometry::operaGeometry(unsigned maxdatapoints, unsigned maxValues) :
Ndatapoints(0), 
ymin(0.0),
ymax(0.0),
apertureWidth(0.0), 
orderLength(0.0),	
NumberofPointsToBinInYDirection(0),
dispersionAxis(cols),			// (axis_t 0=columns 1=rows)
dispersionDirection(up)		// (dispersondirection_t up=wl-increase up, down wl-increase down)
{
	if (maxdatapoints == 0) {
		throw operaException("operaGeometry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (maxValues == 0) {
		throw operaException("operaGeometry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	Ndatapoints = maxdatapoints;
	OrderCenters.vectorlength = 0;
	pixelValueVector.vectorlength = 0;
	createOrderCentersVector(maxdatapoints);
	createPixelValueVector(maxValues);
	geometryPolynomial = new Polynomial();
}

operaGeometry::~operaGeometry() {
	deleteOrderCentersVector();
	deletePixelValueVector();
	if (geometryPolynomial)
		delete geometryPolynomial;
	geometryPolynomial = NULL;
	Ndatapoints = 0;
	OrderCenters.vectorlength = 0;
	pixelValueVector.vectorlength = 0;
}

/*
 * Methods
 */

double operaGeometry::getOrderLength(void) const {
	return orderLength;
}

double operaGeometry::CalculateAndSetOrderLength() {
    orderLength = CalculateDistance(getYmin(), getYmax());
    return orderLength;
}

/*
 * Does not set order length -- used for temporary calculations of distance
 */
double operaGeometry::CalculateDistance(double yMin, double yMax) const {
    
    int npar = getCenterPolynomial()->getOrderOfPolynomial();
    double *par = (double *)getCenterPolynomial()->getVector();
	return LengthofPolynomial(yMin, yMax, par, npar);
}

double operaGeometry::CalculateMinimumYBinSize(double y) const {
	int npar = geometryPolynomial->getOrderOfPolynomial();
	double *par = (double *)geometryPolynomial->getVector();    
    
    double dx = fabs(DiffPolynomialFunction(y,par,npar)); // = dx/dy, assuming dy=1, unit is x-pixel/y-pixel
    
    // in order to stay below Nyquist sampling, i.e. binsize*dx > 0.5 x-pixel
    double minimumYBinSize = (0.5/dx); // binsize unit is in y-pixels

	return minimumYBinSize;    
}

void operaGeometry::setNumberofPointsToBinInYDirection(unsigned Points) {
	NumberofPointsToBinInYDirection = Points;
}

unsigned operaGeometry::getNumberofPointsToBinInYDirection() const {
	return NumberofPointsToBinInYDirection;
}

/*
 * creates a polynomial which traces the entire order. It tries from 2 to coeffs to try to get
 * the best fit by testing chisqr closest to 1.0. Rough first passes do not use errors
 * where refinements do.
 */

void operaGeometry::traceOrder(unsigned coeffs, double &chisqr, bool witherrors) {
	double par[MAXPOLYNOMIAL];
	double errs[MAXPOLYNOMIAL];
	unsigned nparbestfit = coeffs;
	double bestchisqr = BIG;
    
	if (coeffs > MAXPOLYNOMIAL) {
		throw operaException("operaGeometry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (unsigned currentfit=2; currentfit<=coeffs; currentfit++) {
        for	(unsigned i=0; i<currentfit; i++) {
            par[i] = 1.0;
            errs[i] = 0.0;
        }
#ifdef PRINT_DEBUG         
        cerr << "Fitting polynomial of degree = " << currentfit << endl;
#endif
        
		if (witherrors) {
			int errorcode = operaMPFitPolynomial(OrderCenters.vectorlength, OrderCenters.ys, OrderCenters.xs, OrderCenters.errors, currentfit, par, errs, &chisqr);
			if (errorcode <= 0) {
				throw operaException("operaGeometry: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);	
			}
#ifdef PRINT_DEBUG        
			for	(unsigned i=0; i<currentfit; i++) {
				cerr << "p[" << i << "]=" << par[i] << " +/- " << errs[i] << " chisqr = " << chisqr << endl;
			}    
#endif
			
		} else {
			operaLMFitPolynomial(OrderCenters.vectorlength, OrderCenters.ys, OrderCenters.xs, currentfit, par, &chisqr);
#ifdef PRINT_DEBUG        
			for	(unsigned i=0; i<currentfit; i++) {
				cerr << "p[" << i << "]=" << par[i] << " chisqr = " << chisqr << endl;
			}    
#endif    
		}
		if (chisqr < bestchisqr) {
			bestchisqr = chisqr;
			nparbestfit = currentfit;
		}
	}
	for	(unsigned i=0; i<coeffs; i++) {
		par[i] = 1.0;
		errs[i] = 0.0;
	}
	if (witherrors) {
		int errorcode = operaMPFitPolynomial(OrderCenters.vectorlength, OrderCenters.ys, OrderCenters.xs, OrderCenters.errors, nparbestfit, par, errs, &chisqr);
		if (errorcode <= 0) {
			throw operaException("operaGeometry: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
		operaLMFitPolynomial(OrderCenters.vectorlength, OrderCenters.ys, OrderCenters.xs, nparbestfit, par, &chisqr);
	}
	/* DT Apr 25 2013 -- bad interface, the polynomial* is changed in a call
	 * and may well surprise the caller.
	if (geometryPolynomial)
		delete geometryPolynomial;
	geometryPolynomial = new Polynomial(nparbestfit, par, errs);
	geometryPolynomial->setChisqr(chisqr);
	 */
	geometryPolynomial->resize(nparbestfit);
	for	(unsigned i=0; i<nparbestfit; i++) {
		geometryPolynomial->setCoefficient(i, par[i]);
		geometryPolynomial->setCoefficientError(i, errs[i]);
	}
	geometryPolynomial->setChisqr(chisqr);
	/* New Apr 25 2013, make sure order length is set */
	CalculateAndSetOrderLength();
}

Polynomial *operaGeometry::getCenterPolynomial() {
	return geometryPolynomial;
}

const Polynomial *operaGeometry::getCenterPolynomial() const {
	return geometryPolynomial;
}

void operaGeometry::setNdatapoints(unsigned np) {
	if (np > Ndatapoints) {
		throw operaException("operaGeometry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	Ndatapoints = np; 
}

unsigned operaGeometry::getNdatapoints(void) const {
	return Ndatapoints; 
}

double operaGeometry::getYmin(void) const {
	return ymin; 
}

void operaGeometry::setYmin(double min) {
	ymin = min; 
}

double operaGeometry::getYmax(void) const {
	return ymax; 
}

void operaGeometry::setYmax(double max) {
	ymax = max; 
}

double operaGeometry::getCenterX(unsigned index) const {
	if (index >= Ndatapoints) {
		throw operaException("operaGeometry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return OrderCenters.xs[index];
}

double operaGeometry::getCenterY(unsigned index) const {
	if (index >= Ndatapoints) {
		throw operaException("operaGeometry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return OrderCenters.ys[index];
}

double operaGeometry::getCenterV(unsigned index) const {
	if (index >= Ndatapoints) {
		throw operaException("operaGeometry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return OrderCenters.values[index];
}


void operaGeometry::resetCenter(double x, double y, double xerror, unsigned index) {
	if (index >= Ndatapoints) {
		throw operaException("operaGeometry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    OrderCenters.xs[index] = x;
    OrderCenters.errors[index] = xerror;
    OrderCenters.ys[index] = y;    
}

double operaGeometry::getapertureWidth(void) const {
	return apertureWidth; 
}
void operaGeometry::setapertureWidth(double width) {
	apertureWidth = width; 
}

axis_t operaGeometry::getdispersionAxis(void) const {
	return dispersionAxis;
}

void operaGeometry::setdispersionAxis(axis_t axis) {	
	dispersionAxis = axis;
}

dispersiondirection_t operaGeometry::getdispersionDirection(void) const {  
	return dispersionDirection;
}

void operaGeometry::setdispersionDirection(dispersiondirection_t direction) {  
	dispersionDirection = direction;
}

void operaGeometry::createPixelValueVector(unsigned NValues) {
	if (NValues == 0) {
		throw operaException("operaGeometry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	pixelValueVector.maxlength = NValues;
	pixelValueVector.vectorlength = 0;
	pixelValueVector.xs = (double *)malloc(sizeof(double)*NValues);
	pixelValueVector.ys = (double *)malloc(sizeof(double)*NValues);
	pixelValueVector.values = (double *)malloc(sizeof(double)*NValues);
	pixelValueVector.errors = (double *)malloc(sizeof(double)*NValues);
	if (!pixelValueVector.errors) {
		throw operaException("operaGeometry: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	memset(pixelValueVector.xs, 0, sizeof(double)*NValues);
	memset(pixelValueVector.ys, 0, sizeof(double)*NValues);
	memset(pixelValueVector.values, 0, sizeof(double)*NValues);
	memset(pixelValueVector.errors, 0, sizeof(double)*NValues);
}

void operaGeometry::deletePixelValueVector(void) {
	if (pixelValueVector.vectorlength) {
		if (pixelValueVector.xs) 
			free(pixelValueVector.xs);
		pixelValueVector.xs = NULL;
		if (pixelValueVector.ys) 
			free(pixelValueVector.ys);
		pixelValueVector.ys = NULL;
		if (pixelValueVector.values) 
			free(pixelValueVector.values);
		pixelValueVector.values = NULL;
		if (pixelValueVector.errors) 
			free(pixelValueVector.errors);
		pixelValueVector.errors = NULL;
		pixelValueVector.maxlength = 0;
		pixelValueVector.vectorlength = 0;		
	}
}

void operaGeometry::addPixelValue(double x, double y, double value) {
	if (pixelValueVector.vectorlength > pixelValueVector.maxlength) {
		throw operaException("operaGeometry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (pixelValueVector.vectorlength <= pixelValueVector.maxlength) {
		pixelValueVector.xs[pixelValueVector.vectorlength] = x;
		pixelValueVector.ys[pixelValueVector.vectorlength] = y;
		pixelValueVector.values[pixelValueVector.vectorlength] = value;
		pixelValueVector.errors[pixelValueVector.vectorlength] = 0.0;
		pixelValueVector.vectorlength++;
	}
}

void operaGeometry::addPixelValue(double x, double y, double value, double error) {
	if (pixelValueVector.vectorlength > pixelValueVector.maxlength) {
		throw operaException("operaGeometry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (pixelValueVector.vectorlength < pixelValueVector.maxlength) {
		pixelValueVector.xs[pixelValueVector.vectorlength] = x;
		pixelValueVector.ys[pixelValueVector.vectorlength] = y;
		pixelValueVector.values[pixelValueVector.vectorlength] = value;
		pixelValueVector.errors[pixelValueVector.vectorlength] = error;
		pixelValueVector.vectorlength++;
	}
}

void operaGeometry::createOrderCentersVector(unsigned NValues) {
	if (NValues == 0) {
		throw operaException("operaGeometry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	OrderCenters.maxlength = NValues;
	OrderCenters.vectorlength = 0;
	OrderCenters.xs = (double *)malloc(sizeof(double)*NValues);
	OrderCenters.ys = (double *)malloc(sizeof(double)*NValues);
	OrderCenters.values = (double *)malloc(sizeof(double)*NValues);
	OrderCenters.errors = (double *)malloc(sizeof(double)*NValues);
	if (!OrderCenters.errors) {
		throw operaException("operaGeometry: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	memset(OrderCenters.xs, 0, sizeof(double)*NValues);
	memset(OrderCenters.ys, 0, sizeof(double)*NValues);
	memset(OrderCenters.values, 0, sizeof(double)*NValues);
	memset(OrderCenters.errors, 0, sizeof(double)*NValues);
}

void operaGeometry::deleteOrderCentersVector(void) {
	if (OrderCenters.vectorlength) {
		if (OrderCenters.xs)
			free(OrderCenters.xs);
		OrderCenters.xs = NULL;
		if (OrderCenters.ys)
			free(OrderCenters.ys);
		OrderCenters.ys = NULL;
		if (OrderCenters.values)
			free(OrderCenters.values);
		OrderCenters.values = NULL;
		if (OrderCenters.errors)
			free(OrderCenters.errors);
		OrderCenters.errors = NULL;
		OrderCenters.maxlength = 0;
		OrderCenters.vectorlength = 0;		
	}
}

void operaGeometry::addOrderCenterValue(double x, double y, double value) {
	if (OrderCenters.vectorlength < OrderCenters.maxlength) {
		OrderCenters.xs[OrderCenters.vectorlength] = x;
		OrderCenters.ys[OrderCenters.vectorlength] = y;
		OrderCenters.values[OrderCenters.vectorlength] = value;
		OrderCenters.errors[OrderCenters.vectorlength] = 0.0;
		OrderCenters.vectorlength++;
	}
}

void operaGeometry::addOrderCenterValue(double x, double y, double value, double error) {
	if (OrderCenters.vectorlength < OrderCenters.maxlength) {
		OrderCenters.xs[OrderCenters.vectorlength] = x;
		OrderCenters.ys[OrderCenters.vectorlength] = y;
		OrderCenters.values[OrderCenters.vectorlength] = value;
		OrderCenters.errors[OrderCenters.vectorlength] = error;
		OrderCenters.vectorlength++;
	}
}

doubleValueVector_t operaGeometry::getOrderCenters(void) {
	return OrderCenters;
}

doubleValueVector_t operaGeometry::getPixelValueVector(void) {
	return pixelValueVector;
}

void operaGeometry::setorderSeparation(double OrderSeparation) {
	orderSeparation = OrderSeparation;
}

double operaGeometry::getorderSeparation(void) const {
	return orderSeparation;	
}
