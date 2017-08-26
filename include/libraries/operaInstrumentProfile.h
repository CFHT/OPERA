#ifndef OPERAINSTRUMENTPROFILE_H
#define OPERAINSTRUMENTPROFILE_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaInstrumentProfile
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

#include "libraries/operaFITSImage.h"
#include "libraries/Polynomial.h"
#include "libraries/operaVector.h"
#include "libraries/operaMatrix.h"

typedef Matrix<Polynomial> PolynomialMatrix;

/*! 
 * \brief class encapsulating the operaInstrumentProfile.
 * \details The instrument profile (IP) consists of a data set (pixelized image) 
 * \details representing the distribution of the fraction of flux for a uniformly illuminated 
 * \details monochromatic image of the entrance slit projected on the spectrograph focal plane. 
 * \details The fiducial set of coordinates chosen is such that the ordinate is oriented along 
 * \details with the dispersion direction and the abscissa with the spatial direction.
 * \return none
 * \file operaInstrumentProfile.h
 * \ingroup libraries
 */

class operaInstrumentProfile {
	
private:
	unsigned NXPoints;
	unsigned NYPoints;
	unsigned NTotalPoints;
	unsigned Xsampling;		// (fraction of pixel to sample PSF in spatial direction)
	unsigned Ysampling;		// (fraction of pixel to sample PSF in dispersion direction)
	unsigned xsize;			// (size along the x direction in pix units) 
	unsigned ysize;			// (size along the y direction in pix units)
	double geometricCenterX; // geometric x center in pixel units 
	double geometricCenterY; // geometric y center in pixel units 
	
	unsigned nDataPoints;
	unsigned maxnDataPoints;
	DCube dataCube;
	
	operaVector distd;
	
	PolynomialMatrix ipPolyModel;
    DMatrix chisqrMatrix;
	
public:
	
	/*
	 * Constructor
	 */
	operaInstrumentProfile(void);
	
	operaInstrumentProfile(unsigned ipxsize,unsigned ipxsampling,unsigned ipysize,unsigned ipysampling);
	
	operaInstrumentProfile(unsigned ipxsize,unsigned ipxsampling,unsigned ipysize,unsigned ipysampling,unsigned NDataPoints);
	
	/*
	 * Methods
	 */
	
	unsigned getNXPoints(void) const;
	unsigned getNYPoints(void) const;
	unsigned getNTotalPoints(void) const;
	unsigned getXsampling(void) const;
	unsigned getYsampling(void) const;
	unsigned getxsize(void) const;
	unsigned getysize(void) const;
	double getGeometricCenterX(void) const;
	double getGeometricCenterY(void) const;
	
	void setNXPoints(unsigned npx);
	void setNYPoints(unsigned npy);
	void setNTotalPoints(unsigned np);
	void setXsampling(unsigned xsamp);
	void setYsampling(unsigned ysamp);
	void setsampling(unsigned xsamp, unsigned ysamp);
	void setxsize(unsigned xs);
	void setysize(unsigned ys);
	void setsize(unsigned xs, unsigned ys);
	void setGeometricCenter(void); // set x,y coordinates of the geometric center (requires x,y size and x,y sampling)
	
	double getIPixXCoordinate(unsigned i) const;
	double getIPixYCoordinate(unsigned j) const;
	unsigned getIPixiIndex(double xcoord) const;
	unsigned getIPixjIndex(double ycoord) const;
	
	unsigned getnDataPoints(void) const;
	void setnDataPoints(unsigned NDataPoints);
	
	DMatrix getdataCubeValues(unsigned index) const;
	double getdataCubeValues(unsigned i, unsigned j, unsigned index) const;
    void setdataCubeValues(DMatrix DataMatrix, unsigned index);
    void setdataCubeValues(double DataValue, unsigned i, unsigned j, unsigned index);
	
	double getdataGivenCoords(double xcoord, double ycoord, unsigned index) const;
	
	void setdataCubeValues(operaFITSImage &image, operaFITSImage &badpix, double xcenter, double ycenter, unsigned index);
	
	void normalizeCubeData(void); 			// normalize data to Sum = 1
    
	void normalizeCubeData(unsigned index);
    
	void subtractOuterFrame(unsigned index);    
    
	void printModel(double DistanceInPixels, unsigned ordernumber, ostream *pout);
    
	const operaVector& getdistdVector(void) const;
	
	double getdistd(unsigned index) const;
	
	void setdistdVector(const operaVector& Distd);
	
	void setdistd(double DistdValue, unsigned index);
	
	void setipPolyModel(const PolynomialMatrix& IPPolyModel);
	
	const PolynomialMatrix& getipPolyModel(void);
	
	const Polynomial& getipPolyModelCoefficients(unsigned i,unsigned j) const;
    
	void setipPolyModelCoefficients(const Polynomial& PolyModelCoeffs,unsigned i,unsigned j);
    
	DMatrix getipDataFromPolyModel(double d) const;
    
    double getipDataFromPolyModel(double d, unsigned i, unsigned j) const;
	
	void FitPolyMatrixtoIPDataVector(unsigned coeffs, bool witherrors);
    
    void FitMediantoIPDataVector(void);
    
    double getchisqrMatrixValue(unsigned i, unsigned j) const;
	
	void setchisqrMatrixValue(double Value, unsigned i, unsigned j);
    
	void setdataCubeFromPolyModel(void);
    
    void setdataCubeFromPolyModel(unsigned index);
    
	void deleteDataCubes(void);
    
    double getIPphotoCenterX(double d) const;
    
    double getIPphotoCenterY(double d) const;
};

#endif
