#ifndef OPERASPECTRALFEATURE_H
#define OPERASPECTRALFEATURE_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralFeature
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

#include "libraries/Gaussian.h"

#ifndef MAXLINESPERFEATURE
#define MAXLINESPERFEATURE 20
#endif

/*!
 * \file operaSpectralFeature.h
 */

/*! 
 * \sa class operaSpectralFeature
 * \brief This class is used to manipulate a set of spectral lines within a feature.
 * \ingroup libraries
 */
class operaSpectralFeature {
	
private:

	unsigned nLines;            // number of spectral lines in the feature    
	unsigned maxnLines;         // maximum number of spectral lines in the feature    
    
	double BackgroundSlope;     // slope of a linear fit to background
	double BackgroundIntercept; // intercept of a linear fit to background
    
    Gaussian *gaussianFit;      // gaussian fit to lines
    
    unsigned nDataPoints;
    unsigned maxnDataPoints;
    unsigned originalIndex[2];  // initial and final index to original data
    double *xdata;              // x data
    double *ydata;              // y data
    double *yerrors;            // y errors    
    
public:
	
	/*
	 * Constructor
	 */
	operaSpectralFeature(void);
	
	operaSpectralFeature(unsigned NLines);
	
	operaSpectralFeature(unsigned NLines, unsigned MaxNDataPoints);
    
	/*
	 * Destructor
	 */
	~operaSpectralFeature(void);
	
	/*
	 * Setters/Getters
	 */
    
    void setnLines(unsigned NLines);
    
    unsigned getnLines(void) const;
    
    void setBackgroundSlope(double Slope);
    
    double getBackgroundSlope(void) const;
    
    void setBackgroundIntercept(double Intercept);
    
    double getBackgroundIntercept(void) const;
    
    void setGaussianFit(Gaussian *GaussianFit);
    
    const Gaussian *getGaussianFit(void) const;
    
    Gaussian *getGaussianFit(void);

    unsigned getNDataPoints(void) const;
    
    void setOriginalIndex(unsigned i0,unsigned i1);
    
    unsigned getOriginalInitialIndex(void) const;
 
    unsigned getOriginalFinalIndex(void) const;
    
    void CreateDataVectors(unsigned NDataPoints);    
    
	void DeleteDataVectors(void);
	
    void setdataVector(unsigned NDataPoints, double *Xdata, double *Ydata, double *Yerrors);
    
    void setDataPoint(double Xdata, double Ydata, double Yerrors, unsigned index);    

    /*
	 * Methods
	 */
    
    void fitGaussianModel(void);

    void fitBackground(void);

};

#endif
