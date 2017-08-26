/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralFeature
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
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaException.h"

#include "libraries/operaFit.h"
#include "libraries/operaMath.h"
#include "libraries/ladfit.h"

/*!
 * operaSpectralFeature
 * \author Eder Martioli
 * \brief This class encapsulates a spectral feature.
 * \file operaSpectralFeature.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaSpectralFeature
 * \brief TThis class encapsulates a specral feature.
 * \return none
 */

/*
 * Constructor
 */

operaSpectralFeature::operaSpectralFeature(void) :
nLines(0),
maxnLines(0),
BackgroundSlope(0),
BackgroundIntercept(0),
gaussianFit(NULL),
nDataPoints(0),
maxnDataPoints(0),
xdata(NULL),
ydata(NULL),
yerrors(NULL)
{    
}

operaSpectralFeature::operaSpectralFeature(unsigned NLines) :
nLines(0),
maxnLines(0),
BackgroundSlope(0),
BackgroundIntercept(0),
gaussianFit(NULL),
nDataPoints(0),
maxnDataPoints(0),
xdata(NULL),
ydata(NULL),
yerrors(NULL)
{
	if (NLines == 0) {
		throw operaException("operaSpectralFeature: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    nLines = maxnLines = NLines;
    gaussianFit = new Gaussian(getnLines());     
}

operaSpectralFeature::operaSpectralFeature(unsigned NLines, unsigned MaxNDataPoints) :
nLines(0),
maxnLines(0),
BackgroundSlope(0),
BackgroundIntercept(0),
gaussianFit(NULL),
nDataPoints(0),
maxnDataPoints(0),
xdata(NULL),
ydata(NULL),
yerrors(NULL)
{
	if (NLines == 0) {
		throw operaException("operaSpectralFeature: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    nLines = maxnLines = NLines;
        
    CreateDataVectors(MaxNDataPoints);
    
    gaussianFit = new Gaussian(getnLines()); 
}

/*
 * Destructor
 */
operaSpectralFeature::~operaSpectralFeature(void) {
	DeleteDataVectors();
	if (gaussianFit) {
		delete gaussianFit;
	}
	gaussianFit = NULL;
}

/*
 * Setter/Getters
 */

void operaSpectralFeature::setnLines(unsigned NLines) {
	if (NLines > maxnLines || NLines == 0) {
		throw operaException("operaSpectralFeature: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    nLines = maxnLines = NLines;
}

unsigned operaSpectralFeature::getnLines(void) const {
    return nLines;
}

double operaSpectralFeature::getBackgroundSlope(void) const {
    return BackgroundSlope;
}

void operaSpectralFeature::setBackgroundSlope(double Slope) {
    BackgroundSlope = Slope;
}

double operaSpectralFeature::getBackgroundIntercept(void) const {
    return BackgroundIntercept;
}

void operaSpectralFeature::setBackgroundIntercept(double Intercept) {
    BackgroundIntercept = Intercept;
}

void operaSpectralFeature::setGaussianFit(Gaussian *GaussianFit) {
    gaussianFit = GaussianFit;
}

const Gaussian *operaSpectralFeature::getGaussianFit(void) const {
    return gaussianFit;
}

Gaussian *operaSpectralFeature::getGaussianFit(void) {
    return gaussianFit;
}

unsigned operaSpectralFeature::getNDataPoints(void) const {
    return nDataPoints;
}

void operaSpectralFeature::setOriginalIndex(unsigned i0,unsigned i1) {
    originalIndex[0] = i0;
    originalIndex[1] = i1;
}

unsigned operaSpectralFeature::getOriginalInitialIndex(void) const {
    return originalIndex[0];    
}

unsigned operaSpectralFeature::getOriginalFinalIndex(void) const {
    return originalIndex[1];   
}

void operaSpectralFeature::CreateDataVectors(unsigned NDataPoints) {
	if (NDataPoints == 0) {
		throw operaException("operaSpectralFeature: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (NDataPoints <= maxnDataPoints && maxnDataPoints != 0) {
		nDataPoints = NDataPoints;
	} else {
		nDataPoints = maxnDataPoints = NDataPoints;
		if (xdata) {
			free(xdata);
		}
		xdata = (double *)malloc(NDataPoints*sizeof(double));	
		if (ydata) {
			free(ydata);
		}
		ydata = (double *)malloc(NDataPoints*sizeof(double));	
		if (yerrors) {
			free(yerrors);
		}
		yerrors = (double *)malloc(NDataPoints*sizeof(double));    
		if (!yerrors) {
			throw operaException("operaSpectralFeature: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
		}
	}
}

void operaSpectralFeature::DeleteDataVectors(void) {
    if(xdata)
        free(xdata);
	xdata = NULL;
    if(ydata)
        free(ydata); 
	ydata = NULL;
    if(yerrors)
        free(yerrors); 
	yerrors = NULL;
}

void operaSpectralFeature::setdataVector(unsigned NDataPoints, double *Xdata, double *Ydata, double *Yerrors) {
	if (NDataPoints == 0) {
		throw operaException("operaSpectralFeature: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (NDataPoints > nDataPoints) {
		DeleteDataVectors();
		CreateDataVectors(NDataPoints);	
	}
	nDataPoints = NDataPoints;
    for(unsigned i=0;i<nDataPoints;i++) {
        xdata[i] = Xdata[i];
        ydata[i] = Ydata[i];
        yerrors[i] = Yerrors[i];
    }    
}

void operaSpectralFeature::setDataPoint(double Xdata, double Ydata, double Yerrors, unsigned index) {
	if (index > nDataPoints) {
		throw operaException("operaSpectralFeature: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    xdata[index] = Xdata;  
    ydata[index] = Ydata;    
    yerrors[index] = Yerrors;  
}
/*
 * Methods
 */

void operaSpectralFeature::fitBackground(void) {
    unsigned backgroundAperture = 2+nDataPoints/(4*nLines); // -- bug fix for backgroundAperture == 0 -- DT May 13 2013
	if (nDataPoints == 0) {
		throw operaException("operaSpectralFeature: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (backgroundAperture == 0) {
		throw operaException("operaSpectralFeature: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    float *bkgleftsideYData = (float *)malloc(backgroundAperture*sizeof(float));    
    float *bkgrightsideYData = (float *)malloc(backgroundAperture*sizeof(float));  
	if (!bkgrightsideYData) {
		throw operaException("operaSpectralFeature: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    unsigned nBkgDataPoints = 0;
    double leftsideX = 0;  
    double leftsideY = 0;    
    
    for (unsigned i=0;i<backgroundAperture;i++) {
        leftsideX += xdata[i];
        bkgleftsideYData[nBkgDataPoints++] = (float)ydata[i];
    }
    if(nBkgDataPoints) { 
        leftsideX /= (double)nBkgDataPoints;
        leftsideY = (double)operaArrayMedianQuick(nBkgDataPoints,bkgleftsideYData);
    }
    
    nBkgDataPoints = 0;
    double rightsideX = 0;  
    double rightsideY = 0;    

    for (unsigned i=nDataPoints-backgroundAperture;i<nDataPoints;i++) {       
        rightsideX += xdata[i];
        bkgrightsideYData[nBkgDataPoints++] = (float)ydata[i];
    }
    if(nBkgDataPoints) {     
        rightsideX /= (double)nBkgDataPoints;
        rightsideY = (double)operaArrayMedianQuick(nBkgDataPoints,bkgrightsideYData);      
    }

    double tmp_intercept;
    
    if(rightsideX != leftsideX) {
        BackgroundSlope = (rightsideY - leftsideY)/(rightsideX - leftsideX); 
    }
    tmp_intercept = leftsideY - leftsideX*BackgroundSlope;    
  
    double minX=0, minY=BIG;
    
    for (unsigned i=backgroundAperture;i<(nDataPoints-backgroundAperture);i++) {
        if(ydata[i] < minY) {
            minY = ydata[i];
            minX = xdata[i];
        }
    }    
    
    double interceptShift = minY - (tmp_intercept+BackgroundSlope*minX);
    if(interceptShift < 0 && minX) {
        tmp_intercept += interceptShift;
    }
    
    BackgroundIntercept = tmp_intercept;
    
#ifdef PRINT_DEBUG      
    for (unsigned i=0;i<nDataPoints;i++) {
        cout << nLines << " " << xdata[i] << " " << ydata[i] << " " << yerrors[i] << " " << (ydata[i] - BackgroundIntercept - BackgroundSlope*xdata[i]) << " " << (BackgroundIntercept + BackgroundSlope*xdata[i]) << endl; 
    }
#endif          
    
    free(bkgleftsideYData);
    free(bkgrightsideYData);    
}

void operaSpectralFeature::fitGaussianModel(void) {
   
    gaussianFit->MPFitModeltoData(nLines,BackgroundIntercept,BackgroundSlope,nDataPoints,xdata,ydata,yerrors);
    
    BackgroundIntercept = gaussianFit->getBaselineIntercept();
    
    BackgroundSlope = gaussianFit->getBaselineSlope();    
}
