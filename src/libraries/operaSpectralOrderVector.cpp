/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralOrderVector
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

#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMEFFITSProduct.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaGeometry.h"		// for calculate order length
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaWavelength.h" // MAXORDEROFWAVELENGTHPOLYNOMIAL and MAXREFWAVELENGTHSPERORDER
#include "libraries/operaSpectralLines.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/GainBiasNoise.h"
#include "libraries/Polynomial.h"
#include "libraries/LaurentPolynomial.h"  

#include "libraries/operastringstream.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaFit.h"
#include "libraries/operaStats.h"
#include "libraries/gzstream.h"
#include "libraries/operaFFT.h"    
#include "libraries/ladfit.h" // for ladfit_d

#include "libraries/operaSpectralTools.h"			// void calculateUniformSample, getFluxAtWavelength

/*!
 * operaSpectralOrderVector
 * \author Doug Teeple
 * \brief spectral order vector.
 * \details {This library contains serialization and deserializatin of spectral orders.}
 * \file operaSpectralOrderVector.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaSpectralOrderVector();
 * \brief Base constructor.
 */
operaSpectralOrderVector::operaSpectralOrderVector() :
vector(MAXORDERS+1),
orderSpacingPolynomial(NULL),
length(0),
minorder(0),
maxorder(0),
sequence(0),
instrumentmode(MODE_UNKNOWN),
count(0),
gainBiasNoise(NULL)
{
	unsigned order = 0;
	for (order=minorder; order<MAXORDERS; order++) {
		vector[order] = new operaSpectralOrder(order);
	}
	vector[order] = NULL;
	length = MAXORDERS;
	orderSpacingPolynomial = new Polynomial();
    for (unsigned dispIndex=0; dispIndex<MAXORDEROFWAVELENGTHPOLYNOMIAL; dispIndex++) {
        dispersionPolynomial[dispIndex] = new LaurentPolynomial();
    }
	gainBiasNoise = new GainBiasNoise();
}
/* 
 * \class operaSpectralOrderVector(unsigned length, unsigned maxdatapoints, unsigned maxValues, unsigned nElements);
 * \brief Create a NULL-terminated SpectralOrderVector of spectralorders of type "None".
 */
operaSpectralOrderVector::operaSpectralOrderVector(unsigned Length, unsigned maxdatapoints, unsigned maxValues, unsigned nElements) :
vector(MAXORDERS+1),
orderSpacingPolynomial(NULL),
length(0),
minorder(0),
maxorder(0),
sequence(0),
instrumentmode(MODE_UNKNOWN),
count(0),
gainBiasNoise(NULL)
{
	if (Length == 0) {
		throw operaException("operaSpectralOrderVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (Length > MAXORDERS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned order = 0;
	length = Length;
	for (order=minorder; order<MAXORDERS; order++) {
		vector[order] = new operaSpectralOrder(order, maxdatapoints, maxValues, nElements, None);
	}
	vector[order] = NULL;
	orderSpacingPolynomial = new Polynomial();
    for (unsigned dispIndex=0; dispIndex<MAXORDEROFWAVELENGTHPOLYNOMIAL; dispIndex++) {
        dispersionPolynomial[dispIndex] = new LaurentPolynomial();
    }
	gainBiasNoise = new GainBiasNoise();
}
/*
 * Destructor
 */
operaSpectralOrderVector::~operaSpectralOrderVector() {
	freeSpectralOrderVector();
	delete orderSpacingPolynomial;
	orderSpacingPolynomial = NULL;
    for (unsigned dispIndex=0; dispIndex<MAXORDEROFWAVELENGTHPOLYNOMIAL; dispIndex++) {
        delete dispersionPolynomial[dispIndex];
        dispersionPolynomial[dispIndex] = NULL;
    }
	
	delete gainBiasNoise;
	gainBiasNoise = NULL;
}

/*
 * Methods
 */

/*!
 * unsigned getnumberOfDispersionPolynomials(void);
 * \brief returns the number of dispersion polynomials.
 * \return unsigned - numberOfDispersionPolynomials.
 */
unsigned operaSpectralOrderVector::getnumberOfDispersionPolynomials(void) const {
    return numberOfDispersionPolynomials;
}

/*!
 * void setnumberOfDispersionPolynomials(unsigned NumberOfDispersionPolynomials);
 * \brief sets the number of dispersion polynomials.
 * \return none.
 */
void operaSpectralOrderVector::setnumberOfDispersionPolynomials(unsigned NumberOfDispersionPolynomials) {\
	if (NumberOfDispersionPolynomials > MAXORDEROFWAVELENGTHPOLYNOMIAL) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    numberOfDispersionPolynomials = NumberOfDispersionPolynomials;
}

void operaSpectralOrderVector::freeSpectralOrderVector() {
	for (unsigned order = 0; order < MAXORDERS; order++) {
		delete vector[order];
	}
}
/* 
 * unsigned getGainBiasNoise();
 * \brief returns a pointer to the GainBiasNoise class instance.
 */
GainBiasNoise *operaSpectralOrderVector::getGainBiasNoise() {
	return gainBiasNoise;
}
const GainBiasNoise *operaSpectralOrderVector::getGainBiasNoise() const {
	return gainBiasNoise;
}

void operaSpectralOrderVector::setWavelengthsFromCalibration(int Minorder, int Maxorder) {
	for(unsigned order=(unsigned)Minorder; order<=(unsigned)Maxorder; order++) {
		GetSpectralOrder(order)->setWavelengthsFromCalibration();
	}
}

/* 
 * unsigned getCount();
 * \brief returns the count of spectral orders that have content.
 */
unsigned operaSpectralOrderVector::getCount() const {
	return count;
}

/* 
 * unsigned getMinorder();
 * \brief returns the least order number in the vector
 */
unsigned operaSpectralOrderVector::getMinorder() const {
	return minorder;
}

/* 
 * unsigned getMaxorder();
 * \brief returns the maximal order number in the vector
 */
unsigned operaSpectralOrderVector::getMaxorder() const {
	return maxorder;
}

/* 
 * unsigned setCount();
 * \brief sets the count of spectral orders that have content.
 */
void operaSpectralOrderVector::setCount(unsigned Count) {
#ifdef RANGE_CHECK
	if (Count > length) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	count = Count;
}

/* 
 * unsigned setMinorder();
 * \brief sets the least order number in the vector
 */
void operaSpectralOrderVector::setMinorder(unsigned Minorder) {
#ifdef RANGE_CHECK
	if (Minorder > length) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	minorder = Minorder;
}

/* 
 * void getMaxorder(unsigned Maxorder);
 * \brief sets the maximal order number in the vector
 */
void operaSpectralOrderVector::setMaxorder(unsigned Maxorder) {
#ifdef RANGE_CHECK
	if (Maxorder > length) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	maxorder = Maxorder;
}


void operaSpectralOrderVector::shiftOrdersDown(unsigned shift) {
	std::vector<operaSpectralOrder*> temp(MAXORDERS);
	for(unsigned order = 0; order < MAXORDERS-shift; order++) {
		temp[order] = vector[order+shift];
		temp[order]->setorder(order);
	}
	for(unsigned order = 0; order < shift; order++) {
		delete vector[order];
	}
	for(unsigned order = MAXORDERS-shift; order < MAXORDERS; order++) {
		temp[order] = new operaSpectralOrder(order);
	}
	vector = temp;
}

void operaSpectralOrderVector::shiftOrdersUp(unsigned shift) {
	std::vector<operaSpectralOrder*> temp(MAXORDERS);
	for(unsigned order = shift; order < MAXORDERS; order++) {
		temp[order] = vector[order-shift];
		temp[order]->setorder(order);
	}
	for(unsigned order = 0; order < shift; order++) {
		temp[order] = new operaSpectralOrder(order);
	}
	for(unsigned order = MAXORDERS-shift; order < MAXORDERS; order++) {
		delete vector[order];
	}
	vector = temp;
}

void operaSpectralOrderVector::shiftOrders(int shift) {
	if(shift < 0) shiftOrdersDown(abs(shift));
	else shiftOrdersUp(shift);
	setMinorder(minorder+shift);
	setMaxorder(maxorder+shift);
}

/*
 * \sa void setObject(string object)
 * \brief sets the object name
 * \return none.
 */
void operaSpectralOrderVector::setObject(string Object) {
	object = Object;
}

/*
 * \sa string getObject(void);
 * \brief get the object name
 * \return none.
 */
string operaSpectralOrderVector::getObject(void) const {
	return object;
}

/* 
 * \sa setSequence(unsigned sequence)
 * \brief sets the sequence number
 * \return none.
 */
void operaSpectralOrderVector::setSequence(unsigned Sequence) {
	sequence = Sequence;
}

/* 
 * \sa unsigned getSequence(void);
 * \brief get the sequence number
 * \return none.
 */
unsigned operaSpectralOrderVector::getSequence(void) const {
	return sequence;
}
/* 
 * Polynomial *getOrderSpacingPolynomial(void);
 * \brief gets the order spacing polynomial
 */
Polynomial *operaSpectralOrderVector::getOrderSpacingPolynomial(void) {
	return orderSpacingPolynomial;
}
const Polynomial *operaSpectralOrderVector::getOrderSpacingPolynomial(void) const {
	return orderSpacingPolynomial;
}

/* 
 * setOrderSpacingPolynomial(PolynomialCoeffs_t *pc);
 * \brief sets the order spacing polynomial
 */
void operaSpectralOrderVector::setOrderSpacingPolynomial(PolynomialCoeffs_t *pc) {
	delete orderSpacingPolynomial;
	orderSpacingPolynomial = new Polynomial(pc);
}


/*
 * \sa method Polynomial *getDispersionPolynomial(unsigned index);
 * \brief gets the dispersion polynomial
 */
LaurentPolynomial *operaSpectralOrderVector::getDispersionPolynomial(unsigned index) {
    if (index >= numberOfDispersionPolynomials) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	return dispersionPolynomial[index];
}
const LaurentPolynomial *operaSpectralOrderVector::getDispersionPolynomial(unsigned index) const {
    if (index >= numberOfDispersionPolynomials) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	return dispersionPolynomial[index];
}

/*
 * setDispersionPolynomial(unsigned index, const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, PolynomialCoeffs_t *pc);
 * \brief sets the dispersion polynomial
 */
void operaSpectralOrderVector::setDispersionPolynomial(unsigned index, const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, PolynomialCoeffs_t *pc) {
    if (index >= MAXORDEROFWAVELENGTHPOLYNOMIAL) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	delete dispersionPolynomial[index];
	dispersionPolynomial[index] = new LaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial,pc);
}

/* 
 * operaSpectralOrder* operaSpectralOrderVector::GetSpectralOrder(unsigned order);
 * \brief Gets an operaSpectralOrder* to a given order, else NULL
 */
const operaSpectralOrder* operaSpectralOrderVector::GetSpectralOrder(unsigned order) const {
	return vector[order]; //Why don't we just do this instead? Should always work?
}

operaSpectralOrder* operaSpectralOrderVector::GetSpectralOrder(unsigned order) {
	return vector[order];
}


void operaSpectralOrderVector::fitOrderSpacingPolynomial(operaFITSImage &masterFlatImage, operaFITSImage &badpixImage, float slit, unsigned nsamples, unsigned sampleCenterPosition, unsigned referenceOrderNumber, float referenceOrderSeparation, int detectionMethod, bool FFTfilter, float gain, float noise, unsigned x1, unsigned x2, unsigned y1, unsigned y2, unsigned cleanbinsize, float nsigcut, ostream *pout) {
    
    unsigned nx = x2 - x1;
    unsigned ny = y2 - y1;
    
    float *fx = (float *) malloc (nx * sizeof(float));
    float *fy = (float *) malloc (ny * sizeof(float));
    float *fytmp = (float *) malloc (ny * sizeof(float));
    float *fyerr = (float *) malloc (ny * sizeof(float));
	if (!fyerr) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);
	}
    float xmean[MAXORDERS],xmeanerr[MAXORDERS],ymean[MAXORDERS];
    
    float TemporaryOrderPosition[MAXORDERS], TemporaryOrderSeparation[MAXORDERS],TemporaryOrderSeparationError[MAXORDERS];
    
    unsigned nords=0;
    float threshold = DETECTTHRESHOLD;
    float sigma = slit/4;
    
    float *fys = (float *) malloc ((nsamples+1)* sizeof(float));
	if (!fys) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);
	}
    
    // Sample nsamples (default 3) equally spaced regions of the detector at e.g. 1/4, 1/2, 3/4 in the vertical direction
    // what we are doing here is taking three(NSAMPLE) samples to try to determine:
    // what is the polynomial that descibes the function of spacing between orders in the x direction
    unsigned y0,yf;
    if(sampleCenterPosition - nsamples/2 < y1) {
        y0 = y1;
        yf = y1 + nsamples;
    } else {
        y0 = sampleCenterPosition - nsamples/2;
        yf = y0 + nsamples;
    }
    if(yf > y2) {
        yf=y2;
    }
    
    unsigned np = 0;
    
    for (unsigned x=x1; x<nx; x++) {
        fx[np] = (float)x + 0.5;
        
        unsigned ns=0;
        for (unsigned y=y0; y<yf; y++) {
            if(badpixImage[y][x] == 1) {
                fys[ns++] = masterFlatImage[y][x];
            }
        }
        
        if (ns == 0) {
            if(FFTfilter){
                fytmp[np] = 0.0;	// SHOULD BE NAN, but library doesn't handle it...
                fyerr[np] = 0.0;	// SHOULD BE NAN, but library doesn't handle it...
            } else {
                fy[np] = 0.0;	// SHOULD BE NAN, but library doesn't handle it...
                fyerr[np] = 0.0;	// SHOULD BE NAN, but library doesn't handle it...
            }
        } else {
            if(FFTfilter){
                fytmp[np] = operaArrayMedian(ns,fys);
                fyerr[np] = operaArrayMedianSigma(ns, fys, fytmp[np]);
            } else {
                fy[np] = operaArrayMedian(ns,fys);
                fyerr[np] = operaArrayMedianSigma(ns, fys, fy[np]);
            }
        }
        
#ifdef PRINT_DEBUG
        cerr << np << ' ' << fx[np] << ' ' << fy[np] << ' ' << endl;
#endif
        np++;
    }
    
    if(FFTfilter){
        operaFFTLowPass(np,fytmp,fy, 0.1);
    }
    
    if(detectionMethod == 1) {
#ifdef WITH_ERRORS
        nords = operaCCDDetectPeaksWithErrorsUsingGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean,xmeanerr);
#else
        nords = operaCCDDetectPeaksWithGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean);
#endif
    } else if (detectionMethod == 2) {
#ifdef WITH_ERRORS
        nords = operaCCDDetectPeaksWithErrorsUsingGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean,xmeanerr);
#else
        nords = operaCCDDetectPeaksWithGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean);
#endif
    } else if (detectionMethod == 3) {
#ifdef WITH_ERRORS
        nords = operaCCDDetectPeaksWithErrorsUsingTopHat(np,fx,fy,(unsigned)slit,noise,gain,threshold,xmean,ymean,xmeanerr);
#else
        nords = operaCCDDetectPeaksWithTopHat(np,fx,fy,(unsigned)slit,noise,gain,threshold,xmean,ymean);
#endif
    }
    
    unsigned npspc = 0;
    
    for(unsigned i=1;i<nords;i++) {
        TemporaryOrderSeparation[npspc] = fabs(xmean[i] - xmean[i-1]);
        TemporaryOrderSeparationError[npspc] = sqrt(xmeanerr[i]*xmeanerr[i] + xmeanerr[i-1]*xmeanerr[i-1]);
        TemporaryOrderPosition[npspc] = xmean[i-1] + (xmean[i] + xmean[i-1])/2;
        
        if(i > 1) {
            float ratio = TemporaryOrderSeparation[npspc]/TemporaryOrderSeparation[npspc - 1];
            int roundedRatio = (int)round(ratio);
            
            // If ratio>1 it means the current separation is at least twice larger than the size of
            // previous separation. Therefore the current has skipped one (or more) order(s).
            if(roundedRatio > 1
               && TemporaryOrderSeparationError[npspc] < 0.5*TemporaryOrderSeparation[npspc]
               && TemporaryOrderSeparationError[npspc] < 0.5*TemporaryOrderSeparation[npspc - 1]) {
                TemporaryOrderSeparation[npspc] /= float(roundedRatio);
                
                // Else if ratio<1 it means the current separation is at least twice smaller than the size
                // of previous separations. In this case we need to fix all previous separations.
            } else if (roundedRatio < 1
                       && TemporaryOrderSeparationError[npspc] < 0.5*TemporaryOrderSeparation[npspc]
                       && TemporaryOrderSeparationError[npspc] < 0.5*TemporaryOrderSeparation[npspc - 1]) {
                unsigned j = npspc - 1;
                
                // Fix previous separation
                TemporaryOrderSeparation[j] /= 1.0/float(roundedRatio);
                
                //Check for separations before previous one
                for (unsigned ii = 0; ii < i-2; ii++) {
                    ratio = TemporaryOrderSeparation[j]/TemporaryOrderSeparation[j-1];
                    roundedRatio = (int)round(ratio);
                    if(roundedRatio > 1
                       && TemporaryOrderSeparationError[j] < 0.5*TemporaryOrderSeparation[j]
                       && TemporaryOrderSeparationError[j] < 0.5*TemporaryOrderSeparation[j - 1]) {
                        // Fix separation value.
                        TemporaryOrderSeparation[j] /= float(roundedRatio);
                    }
                    j--;
                }
            }
        }
        npspc++;
    }
    
    int SortIndex[MAXORDERS];
    operaArrayIndexSort((int)npspc,TemporaryOrderPosition,SortIndex);
    
    double OrderPosition[MAXORDERS];
    double OrderSeparation[MAXORDERS],OrderSeparationError[MAXORDERS];
    
    // identify reference order, which will be the order for which OrderSeparation ~ referenceOrderSeparation
    float minResidualSeparation = BIG;
    unsigned minj = 0;
    
    for(unsigned j=0; j<npspc; j++) {
        if(fabs(TemporaryOrderSeparation[SortIndex[j]] - referenceOrderSeparation) < minResidualSeparation) {
            minResidualSeparation = fabs(TemporaryOrderSeparation[SortIndex[j]] - referenceOrderSeparation);
            minj = j;
        }
        
        //OrderPosition[j] = (double)TemporaryOrderPosition[SortIndex[j]];
        OrderSeparation[j] = (double)TemporaryOrderSeparation[SortIndex[j]];
        OrderSeparationError[j] = (double)TemporaryOrderSeparationError[SortIndex[j]];
    }
    
    // Below we are replacing the variable detector position to order number, which
    // is based on the reference order provided by user -> referenceOrderNumber
    
    for(unsigned j=0; j<npspc; j++) {
        OrderPosition[j] = double(referenceOrderNumber + j - minj);
    }
	
    // here is the polynomial fit of spacing between orders across the array in the x direction
    // npar is different in this case than order tracing
    unsigned npars = 3;
    double par[3] = {1,1,1};
    double ecoeffs[3] = {0,0,0};
    double chisqr = 0.0;
    
#ifdef WITH_ERRORS
    int errorcode = operaMPFitPolynomial(npspc, OrderPosition, OrderSeparation, OrderSeparationError, npars, par, ecoeffs, &chisqr);
    if (errorcode <= 0) {
        throw operaException("operaSpectralOrderVector: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);
    }
    
#else
    operaLMFitPolynomial(npspc, OrderPosition, OrderSeparation, npars, par, &chisqr);
#endif
    
    
    /*
     *  Below it applies a sigma clip cleaning
     */
    unsigned binsize = cleanbinsize;
    float nsig = nsigcut;
    
    if (binsize==0 || npspc==0) {
        throw operaException("operaSpectralOrderVector: binsize=0 or nDataPoints=0", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
    }
    unsigned nDataPoints = npspc;
    
    double *cleanxdataVector = new double[nDataPoints];
    double *cleanydataVector = new double[nDataPoints];
    double *cleanyerrorVector = new double[nDataPoints];

    unsigned numberOfCleanPoints = 0;
    
    float *xtmp = new float[nDataPoints];
    float *ytmp = new float[nDataPoints];
    
    for(unsigned i=0; i<nDataPoints; i++) {
        
        int firstPoint = (int)i - (int)binsize;
        int lastPoint = (int)i + (int)binsize + 1;
        
        if(firstPoint < 0) {
            firstPoint = 0;
            lastPoint = 2*(int)binsize + 1;
        }
        if(lastPoint > (int)nDataPoints) {
            lastPoint = (int)nDataPoints;
            firstPoint = (int)nDataPoints - 2*(int)binsize - 1;
            if(firstPoint < 0) {
                firstPoint = 0;
            }
        }
        
        unsigned np = 0;
        for(unsigned ii=(unsigned)firstPoint; ii<(unsigned)lastPoint; ii++) {
            xtmp[np] = (float)OrderPosition[ii];
            ytmp[np] = (float)OrderSeparation[ii];
            np++;
        }
        
        float am,bm,abdevm;
        
        //--- Robust linear fit
        ladfit(xtmp,ytmp,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
        
        //--- Clean up
        float fitMedianSlope = (bm*(float)OrderPosition[i] + am);
        
        if(fabs((float)OrderSeparation[i] - fitMedianSlope) < nsig*abdevm) {
            cleanxdataVector[numberOfCleanPoints] = OrderPosition[i];
            cleanydataVector[numberOfCleanPoints] = OrderSeparation[i];
#ifdef WITH_ERRORS

            cleanyerrorVector[numberOfCleanPoints] = OrderSeparationError[i];  
#else
            cleanyerrorVector[numberOfCleanPoints] = 0.0;
#endif
            numberOfCleanPoints++;
        }
    }
    
    for(unsigned i=0; i<numberOfCleanPoints; i++) {
        OrderPosition[i] = cleanxdataVector[i];
        OrderSeparation[i] = cleanydataVector[i];
        OrderSeparationError[i] = cleanyerrorVector[i];
    }
    nDataPoints = numberOfCleanPoints;
    
    delete[] cleanxdataVector;
    delete[] cleanydataVector;
    delete[] cleanyerrorVector;
    
    npspc = nDataPoints;
    /*
     *  End of sigma clip cleaning
     */
    
    
    /*
     * Perform polynomial fit on clean data
     */
#ifdef WITH_ERRORS
    errorcode = operaMPFitPolynomial(npspc, OrderPosition, OrderSeparation, OrderSeparationError, npars, par, ecoeffs, &chisqr);
    if (errorcode <= 0) {
        throw operaException("operaGeometryCalibration: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);
    }
#else
    operaLMFitPolynomial(npspc, OrderPosition, OrderSeparation, npars, par, &chisqr);
#endif
    
    if (pout != NULL) {
        for(unsigned i=0;i<npspc;i++) {
            *pout <<  OrderPosition[i] << ' ' << OrderSeparation[i] << ' ' << OrderSeparationError[i] << ' ' << PolynomialFunction(OrderPosition[i],par,npars) << endl;
        }
    }
    
    // set in the order spacing polynomial
    PolynomialCoeffs_t orderSpacing;
    orderSpacing.orderofPolynomial = npars;
    orderSpacing.polychisqr = chisqr;
    for (unsigned i=0; i<npars; i++) {
        orderSpacing.p[i] = par[i];
        orderSpacing.e[i] = ecoeffs[i];
    }
    setOrderSpacingPolynomial(&orderSpacing);
    free(fx);
    free(fy);
    free(fytmp);
    free(fyerr);
}

void operaSpectralOrderVector::measureIPAlongRowsFromSamples(operaFITSImage &masterFlatImage, operaFITSImage &badpixImage, float slit, unsigned nsamples, bool FFTfilter, float gain, float noise, unsigned x1, unsigned x2, unsigned y1, unsigned y2,float *ipfunc, float *ipx, float *iperr) {
	
    unsigned nx = x2 - x1;
    unsigned ny = y2 - y1;    
    
    float *fx = (float *) malloc (nx * sizeof(float));
    float *fy = (float *) malloc (ny * sizeof(float));
    float *fytmp = (float *) malloc (ny * sizeof(float));		
	if (!fytmp) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    float xmean[MAXORDERS],ymean[MAXORDERS];	
    
    unsigned nords=0;
    float threshold = DETECTTHRESHOLD;
    float sigma = slit/4;
    unsigned uslit = (unsigned)slit;
    
    unsigned numberofpointsinydirectiontomerge = NPINSAMPLE;
    float *fys = (float *) malloc (NPINSAMPLE * sizeof(float)); 
	if (!fys) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	
    float *ipsample = new float[uslit];
    float *ipxsample = new float[uslit];
    float *ipvar = new float[uslit];
	
    memset(ipfunc, 0, sizeof(float)*uslit);
    memset(ipx, 0, sizeof(float)*uslit);
    memset(ipvar, 0, sizeof(float)*uslit);
    
    for(unsigned k=0; k<nsamples; k++) {
        unsigned np=0;
        for (unsigned x=x1; x<nx; x++) {			
            unsigned ns=0;
            for (unsigned y=(ny-y1)*(k+1)/(nsamples + 1) - numberofpointsinydirectiontomerge/2; y < (ny-y1)*(k+1)/(nsamples + 1) + numberofpointsinydirectiontomerge/2; y++) {
                if(badpixImage[y][x] == 1)
                    fys[ns++] = masterFlatImage[y][x];
            }
            
            fx[np] = (float)x + 0.5;
            
			if (ns == 0) {
				if(FFTfilter){
					fytmp[np] = 0.0;
				} else {
					fy[np] = 0.0;	
				}				
			} else {
				if(FFTfilter){
					fytmp[np] = operaArrayMedian(ns,fys);
				} else {
					fy[np] = operaArrayMedian(ns,fys);	
				}				
			}
			
#ifdef PRINT_DEBUG
            cerr << k << ' ' << fx[np] << ' ' << fy[np] << ' ' << endl;
#endif
            
            np++;
        }
        
        if(FFTfilter){
            operaFFTLowPass(np,fytmp,fy, 0.1);
        }
        
        nords = operaCCDDetectPeaksWithGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean);
        operaCCDFitIP(np,fx,fy,nords,xmean,ymean,ipsample,ipxsample,iperr,uslit);
        
        for(unsigned i=0;i<uslit;i++) {
            ipfunc[i] += ipsample[i];
            ipx[i] += ipxsample[i];
            ipvar[i] += (iperr[i]*iperr[i]);
        } 
    }
    float IPNormalizationFactor = 0;
    for(unsigned i=0;i<uslit;i++) {
        IPNormalizationFactor += ipfunc[i];
        ipx[i] /= (double)nsamples;
    }     
    for(unsigned i=0;i<uslit;i++) {
        ipfunc[i] /= IPNormalizationFactor;
        iperr[i] = sqrt(ipvar[i]/IPNormalizationFactor);
    }      
    delete[] ipsample;
    delete[] ipxsample;
    delete[] ipvar;
    free(fx);
    free(fy);
    free(fytmp);
}



unsigned operaSpectralOrderVector::getElemIndexAndOrdersByWavelength(int *orderForWavelength, unsigned *elemIndexForWavelength, double wavelength) {
    
    /* E.Martioli -- Jul 1 2013
     * It probably needs a range check for the size of *orderForWavelength and *elemIndexForWavelength, 
     * which must be greater than (maxorder - minorder + 1)
     * But I don't know what is the best way to do that.
     */
    unsigned nOrdersSelected = 0;
    
    for(unsigned order=minorder; order<=maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());
            
            double wl0 = SpectralElements->getwavelength(0);
            double wlf = SpectralElements->getwavelength(SpectralElements->getnSpectralElements()-1);

#ifdef PRINT_DEBUG
            cout << wl0 << ' ' << wlf << ' ' << wavelength << ' ' << order << endl;
#endif

            if((wl0 <= wlf && wavelength >= wl0 && wavelength < wlf) ||
               (wl0 >  wlf && wavelength > wlf && wavelength <= wl0) ) { // for ascending wavelengths
                
                orderForWavelength[nOrdersSelected] = order;
                
                for (unsigned elemIndex=0; elemIndex<SpectralElements->getnSpectralElements()-1; elemIndex++) {
                    double elem_wl0 = SpectralElements->getwavelength(elemIndex);
                    double elem_wlf = SpectralElements->getwavelength(elemIndex+1);
                    
                    if((elem_wl0 <= elem_wlf && wavelength >= elem_wl0 && wavelength < elem_wlf) ||
                       (elem_wl0 >  elem_wlf && wavelength > elem_wlf && wavelength <= elem_wl0) ) {
                        if(fabs(wavelength - elem_wl0) <= fabs(wavelength - elem_wlf)) {
                            elemIndexForWavelength[nOrdersSelected] = elemIndex;
                            break;
                        } else if (fabs(wavelength - elem_wl0) > fabs(wavelength - elem_wlf)) {
                            elemIndexForWavelength[nOrdersSelected] = elemIndex+1;
                            break;
                        }
                    }
                }
                
                nOrdersSelected++;
            }
        }
    }
    return nOrdersSelected;
}


unsigned operaSpectralOrderVector::getOrdersByWavelengthRange(int *orderForWavelengthRange, double Range_wl0, double Range_wlf) {
    
    unsigned nOrdersSelected = 0;
    
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());
            
            double wl0 = SpectralElements->getwavelength(0);
            double wlf = SpectralElements->getwavelength(SpectralElements->getnSpectralElements()-1);
                        
#ifdef PRINT_DEBUG
            cout << wl0 << ' ' << wlf << ' ' << Range_wl0 << ' ' << Range_wlf << ' ' << order << endl;
#endif
            
            if((wl0 <= wlf && Range_wl0 >= wl0 && Range_wl0 < wlf) ||
               (wl0 >  wlf && Range_wl0 > wlf && Range_wl0 <= wl0) ||
               (wl0 <= wlf && Range_wlf >= wl0 && Range_wlf < wlf) ||
               (wl0 >  wlf && Range_wlf > wlf && Range_wlf <= wl0)) { // for ascending wavelengths
                
                orderForWavelengthRange[nOrdersSelected] = order;
                nOrdersSelected++;
            }
        }
    }
    return nOrdersSelected;
}

void operaSpectralOrderVector::getOrdersByWavelengthRange(operaWavelengthRange wavelengthRange, int& Minorder, int& Maxorder) {
    unsigned tempmin = 0, tempmax = 0;
    for(unsigned order=(unsigned)Minorder; order<=(unsigned)Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            const operaVector& wl = spectralOrder->getSpectralElements()->getWavelength();
            if(wavelengthRange.contains(wl.first()) || wavelengthRange.contains(wl.last())) {
				if(!tempmin) tempmin = order;
				tempmax = order;
			} else if (tempmin) break;
        }
    }
    Minorder = tempmin;
    Maxorder = tempmax;
}

void operaSpectralOrderVector::measureContinuumAcrossOrders(unsigned binsize, int orderBin, unsigned nsigcut) {
    
    unsigned maxNDataPoints = 0;
    unsigned numberOfBeams = getNumberOfBeams(minorder, maxorder);
    
    int usefulMinorder = (int)minorder;
    int usefulMaxorder = (int)maxorder;
    
    bool hasUsefulMinorder = false;
    
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            
            if(!hasUsefulMinorder) {
                usefulMinorder = order;
                hasUsefulMinorder = true;
            }
            usefulMaxorder = order;
            
			spectralOrder->setWavelengthsFromCalibration();
			
            spectralOrder->calculateContinuum(binsize,nsigcut);
            
            if(spectralOrder->gethasSpectralEnergyDistribution()) {
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                
                if(maxNDataPoints < spectralEnergyDistribution->getnDataPoints()) {
                    maxNDataPoints = spectralEnergyDistribution->getnDataPoints();
                }
            }
        }
    }
    
    operaFluxVector *FluxCalibrationFitVector[MAXORDERS];
    operaFluxVector *BeamFluxCalibrationFitVector[MAXNUMBEROFBEAMS][MAXORDERS];
    
    float *fluxData = new float[MAXORDERS];
    float *orderData = new float[MAXORDERS];
        
    float *fluxBeamData[MAXNUMBEROFBEAMS];
    float *orderBeamData[MAXNUMBEROFBEAMS];
    float orderBeamFlux[MAXNUMBEROFBEAMS];
    unsigned npBeam[MAXNUMBEROFBEAMS];
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        npBeam[beam] = 0;
        fluxBeamData[beam] = new float[MAXORDERS];
        orderBeamData[beam] = new float[MAXORDERS];
    }
    
    float fluxvariance = 0;
    float beamfluxvariance = 0;

    unsigned order_count = 0;
    
    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        int loword = order - orderBin;
        int hiord = order + orderBin;
        if(loword < usefulMinorder) {
            loword = usefulMinorder;
        }
        if(hiord > usefulMaxorder) {
            hiord = usefulMaxorder;
        }
        
        int nord = (hiord - loword + 1);
        
        if(!(nord%2)) {
            if(hiord == usefulMaxorder) {
                loword--;
            } else {
                hiord++;
            }
        }
        
        nord = (hiord - loword + 1);
        
        if(GetSpectralOrder(order)->gethasSpectralElements() && GetSpectralOrder(order)->gethasSpectralEnergyDistribution() && GetSpectralOrder(order)->gethasWavelength()) {
            
            unsigned nDataPointsOfCurrentOrder = GetSpectralOrder(order)->getSpectralEnergyDistribution()->getnDataPoints();
            
            FluxCalibrationFitVector[order_count]  = new operaFluxVector(nDataPointsOfCurrentOrder);
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                BeamFluxCalibrationFitVector[beam][order_count] = new operaFluxVector(nDataPointsOfCurrentOrder);
            }
            
            for(unsigned dataIndex=0;dataIndex < maxNDataPoints; dataIndex++) {
                if(dataIndex>=nDataPointsOfCurrentOrder) {
                    break;
                }
                
                unsigned np = 0;
                float orderFlux = NAN;
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    npBeam[beam] = 0;
                    orderBeamFlux[beam] = NAN;
                }
                
                for(int o=loword; o<=hiord; o++) {
                    operaSpectralOrder *spectralOrder = GetSpectralOrder(o);
                    
                    if(spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution() && spectralOrder->gethasWavelength()) {
                        
                        operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                        
                        if(dataIndex < spectralEnergyDistribution->getnDataPoints() && dataIndex < nDataPointsOfCurrentOrder) {
                            float flux = (float)spectralEnergyDistribution->getfluxData(dataIndex);
                            fluxvariance = flux;
                            
                            if(!isnan(flux)) {
                                orderData[np] = (float)o;
                                fluxData[np] = flux;
                                if(o==order) {
                                    orderFlux = flux;
                                }
                                np++;
                            }
                            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                                operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                                float beamflux = beamSED->getfluxData(dataIndex);
                                beamfluxvariance = beamflux;
                                
                                if(!isnan(beamflux)) {
                                    orderBeamData[beam][npBeam[beam]] = (float)o;
                                    fluxBeamData[beam][npBeam[beam]] = beamflux;
                                    if(o==order) {
                                        orderBeamFlux[beam] = beamflux;
                                    }
                                    npBeam[beam]++;
                                }
                            }
                        }
                        
                    }
                }
                
                if(np && !(np%2)) {
                    np--;
                }
                if(np>=3) {
                    float am,bm,abdevm;

                    ladfit(orderData,fluxData,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
                    
                    float dytop = 0;
                    
                    for(unsigned i=0;i<np;i++){
                        float fitMedianSlope = (bm*orderData[i] + am);
                        
                        if(fabs(fluxData[i] - fitMedianSlope) < abdevm &&
                           fabs(fluxData[i] - fitMedianSlope) > dytop) {
                            dytop = fabs(fluxData[i] - fitMedianSlope);
                        }
                    }
                    if(dytop == 0) {
                        dytop = abdevm;
                    }
                    
                    FluxCalibrationFitVector[order_count]->setflux(double(bm*(float)order + am + dytop),dataIndex);
                    FluxCalibrationFitVector[order_count]->setvariance(double(0.674433*abdevm)*double(0.674433*abdevm),dataIndex);
                } else {
                    FluxCalibrationFitVector[order_count]->setflux((double)orderFlux,dataIndex);
                    FluxCalibrationFitVector[order_count]->setvariance((double)fluxvariance,dataIndex);
                }
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam] && !(npBeam[beam]%2)) {
                        npBeam[beam]--;
                    }
                }
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam]>=3) {
                        float amBeam,bmBeam,abdevmBeam;
                        ladfit(orderBeamData[beam],fluxBeamData[beam],npBeam[beam],&amBeam,&bmBeam,&abdevmBeam); /* robust linear fit: f(x) =  a + b*x */
                        float dytop = 0;
                        
                        for(unsigned i=0;i<npBeam[beam];i++){
                            float fitMedianSlope = (bmBeam*orderBeamData[beam][i] + amBeam);
                            
                            if(fabs(fluxBeamData[beam][i] - fitMedianSlope) < abdevmBeam &&
                               (fluxBeamData[beam][i] - fitMedianSlope) > dytop) {
                                dytop = fluxBeamData[beam][i] - fitMedianSlope;
                            }
                        }
                        
                        if(dytop == 0) {
                            dytop = abdevmBeam;
                        }
                        BeamFluxCalibrationFitVector[beam][order_count]->setflux(double(bmBeam*(float)order + amBeam + dytop),dataIndex);
                        BeamFluxCalibrationFitVector[beam][order_count]->setvariance(double(0.674433*abdevmBeam)*double(0.674433*abdevmBeam),dataIndex);
                    } else {
                        BeamFluxCalibrationFitVector[beam][order_count]->setflux(double(orderBeamFlux[beam]),dataIndex);
                        BeamFluxCalibrationFitVector[beam][order_count]->setvariance(double(beamfluxvariance),dataIndex);
                    }
                } // for(unsigned beam=0; beam < numberOfBeams; beam++) {
            } // for(unsigned dataIndex=0;dataIndex < maxNDataPoints; dataIndex++) {
            order_count++;
        } // if(GetSpectralOrder(order)->gethasSpectralElements() && GetSpectralOrder(order)->gethasSpectralEnergyDistribution()) {
    } // for(int order=(int)minorder; order<=(int)maxorder; order++) {

    order_count = 0;
    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if(spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution() && spectralOrder->gethasWavelength()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            unsigned nDataPointsOfCurrentOrder = spectralEnergyDistribution->getnDataPoints();
            
            for(unsigned dataIndex = 0; dataIndex < nDataPointsOfCurrentOrder; dataIndex++) {
                spectralEnergyDistribution->setfluxData(FluxCalibrationFitVector[order_count]->getflux(dataIndex),dataIndex);
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    beamSED->setfluxData(BeamFluxCalibrationFitVector[beam][order_count]->getflux(dataIndex),dataIndex);
                }
#ifdef PRINT_DEBUG
                cout << order << ' ' << dataIndex << ' '
                << spectralEnergyDistribution->getdistanceData(dataIndex) << ' '
                << spectralEnergyDistribution->getfluxData(dataIndex) << ' '
                << FluxCalibrationFitVector[order_count]->getflux(dataIndex) << ' ';
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    cout << beam << ' ' << beamSED->getfluxData(dataIndex) << ' ' << BeamFluxCalibrationFitVector[beam][order_count]->getflux(dataIndex) << ' ';
                }
                cout << endl;
#endif
            }
#ifdef PRINT_DEBUG
            cout << endl;
#endif
            order_count++;
        }
    }

    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution() && spectralOrder->gethasWavelength()) {
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            spectralEnergyDistribution->populateUncalibratedFluxFromContinuumData();
            
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                beamSED->populateUncalibratedFluxFromContinuumData();
            }
        }
    } //for(int order=(int)minorder; order<=(int)maxorder; order++) {
}


void operaSpectralOrderVector::measureContinuumAcrossOrders(unsigned binsize, int orderBin, unsigned nsigcut, unsigned nOrdersPicked, int *orderForWavelength) {
    
    int usefulMinorder = (int)minorder;
    int usefulMaxorder = (int)maxorder;
    
    bool hasUsefulMinorder = false;

    unsigned maxNDataPoints = 0;
	unsigned numberOfBeams = getNumberOfBeams(minorder, maxorder);
    
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            
            if(!hasUsefulMinorder) {
                usefulMinorder = order;
                hasUsefulMinorder = true;
            }
            usefulMaxorder = order;
            
            spectralOrder->setWavelengthsFromCalibration();
            
            spectralOrder->calculateContinuum(binsize,nsigcut);
            
            if(spectralOrder->gethasSpectralEnergyDistribution()) {
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                
                if(maxNDataPoints < spectralEnergyDistribution->getnDataPoints()) {
                    maxNDataPoints = spectralEnergyDistribution->getnDataPoints();
                }
            }
        }
    }
    
    int actualMinorder = usefulMinorder;
    int actualMaxorder = usefulMaxorder;
    
    if(orderForWavelength[0] < usefulMinorder) {
        actualMinorder = orderForWavelength[0];
    }
    
    if(orderForWavelength[nOrdersPicked-1] > usefulMaxorder) {
        actualMaxorder = orderForWavelength[nOrdersPicked-1];
    }
    
    operaFluxVector *FluxCalibrationFitVector[MAXORDERS];
    operaFluxVector *BeamFluxCalibrationFitVector[MAXNUMBEROFBEAMS][MAXORDERS];
    
    float *fluxData = new float[MAXORDERS];
    float *orderData = new float[MAXORDERS];
    
    float *fluxBeamData[MAXNUMBEROFBEAMS];
    float *orderBeamData[MAXNUMBEROFBEAMS];
    float orderBeamFlux[MAXNUMBEROFBEAMS];
    unsigned npBeam[MAXNUMBEROFBEAMS];
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        npBeam[beam] = 0;
        fluxBeamData[beam] = new float[MAXORDERS];
        orderBeamData[beam] = new float[MAXORDERS];
    }
    
    float fluxvariance = 0;
    float beamfluxvariance = 0;
    
    unsigned order_count = 0;
    
    for(int order=actualMinorder; order<=actualMaxorder; order++) {
        
        int loword = order - orderBin;
        int hiord = order + orderBin;
        if(loword < actualMinorder) {
            loword = actualMinorder;
        }
        if(hiord > actualMaxorder) {
            hiord = actualMaxorder;
        }
        
        int nord = (hiord - loword + 1);
        
        if(!(nord%2)) {
            if(hiord == actualMaxorder) {
                loword--;
            } else {
                hiord++;
            }
        }
        
        nord = (hiord - loword + 1);
                
        if(GetSpectralOrder(order)->gethasSpectralElements() && GetSpectralOrder(order)->gethasSpectralEnergyDistribution()) {
            
            unsigned nDataPointsOfCurrentOrder = GetSpectralOrder(order)->getSpectralEnergyDistribution()->getnDataPoints();
            
            FluxCalibrationFitVector[order_count]  = new operaFluxVector(nDataPointsOfCurrentOrder);
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                BeamFluxCalibrationFitVector[beam][order_count] = new operaFluxVector(nDataPointsOfCurrentOrder);
            }
            
            for(unsigned dataIndex=0;dataIndex < maxNDataPoints; dataIndex++) {
                if(dataIndex>=nDataPointsOfCurrentOrder) {
                    break;
                }
                
                unsigned np = 0;
                float orderFlux = NAN;
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    npBeam[beam] = 0;
                    orderBeamFlux[beam] = NAN;
                }
                
                for(int o=loword; o<=hiord; o++) {
                    operaSpectralOrder *spectralOrder = GetSpectralOrder(o);
                    
                    if(spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
                        
                        operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                        
                        if(dataIndex < spectralEnergyDistribution->getnDataPoints() && dataIndex < nDataPointsOfCurrentOrder) {
                            float flux = spectralEnergyDistribution->getfluxData(dataIndex);
                            fluxvariance = flux;
                            
                            if(!isnan(flux)) {
                                orderData[np] = (float)o;
                                fluxData[np] = flux;
                                if(o==order) {
                                    orderFlux = flux;
                                }
                                np++;
                            }
                            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                                operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                                float beamflux = beamSED->getfluxData(dataIndex);
                                beamfluxvariance = beamflux;
                                
                                if(!isnan(beamflux)) {
                                    orderBeamData[beam][npBeam[beam]] = (float)o;
                                    fluxBeamData[beam][npBeam[beam]] = beamflux;
                                    if(o==order) {
                                        orderBeamFlux[beam] = beamflux;
                                    }
                                    npBeam[beam]++;
                                }
                            }
                        }
                        
                    }
                }
                
                if(np && !(np%2)) {
                    np--;
                }
                
                if(np>=3) {
                    float am,bm,abdevm;
                    ladfit(orderData,fluxData,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
                    
                    float dytop = 0;
                    
                    for(unsigned i=0;i<np;i++){
                        float fitMedianSlope = (bm*orderData[i] + am);
                        
                        if(fabs(fluxData[i] - fitMedianSlope) < abdevm &&
                           fabs(fluxData[i] - fitMedianSlope) > dytop) {
                            dytop = fabs(fluxData[i] - fitMedianSlope);
                        }
                    }
                    
                    if(dytop == 0) {
                        dytop = abdevm;
                    }
                    
                    FluxCalibrationFitVector[order_count]->setflux(double(bm*(float)order + am + dytop),dataIndex);
                    FluxCalibrationFitVector[order_count]->setvariance(double(0.674433*abdevm)*double(0.674433*abdevm),dataIndex);
                } else {
                    FluxCalibrationFitVector[order_count]->setflux((double)orderFlux,dataIndex);
                    FluxCalibrationFitVector[order_count]->setvariance((double)fluxvariance,dataIndex);
                }
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam] && !(npBeam[beam]%2)) {
                        npBeam[beam]--;
                    }
                }
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam]>=3) {
                        float amBeam,bmBeam,abdevmBeam;
                        ladfit(orderBeamData[beam],fluxBeamData[beam],npBeam[beam],&amBeam,&bmBeam,&abdevmBeam); /* robust linear fit: f(x) =  a + b*x */
                        
                        float dytop = 0;
                        
                        for(unsigned i=0;i<npBeam[beam];i++){
                            double fitMedianSlope = (bmBeam*orderBeamData[beam][i] + amBeam);
                            
                            if(fabs(fluxBeamData[beam][i] - fitMedianSlope) < abdevmBeam &&
                               (fluxBeamData[beam][i] - fitMedianSlope) > dytop) {
                                dytop = fluxBeamData[beam][i] - fitMedianSlope;
                            }
                        }
                        
                        if(dytop == 0) {
                            dytop = abdevmBeam;
                        }
                        BeamFluxCalibrationFitVector[beam][order_count]->setflux(double(bmBeam*(float)order + amBeam + dytop),dataIndex);
                        BeamFluxCalibrationFitVector[beam][order_count]->setvariance(double(0.674433*abdevmBeam)*double(0.674433*abdevmBeam),dataIndex);
                    } else {
                        BeamFluxCalibrationFitVector[beam][order_count]->setflux((double)orderBeamFlux[beam],dataIndex);
                        BeamFluxCalibrationFitVector[beam][order_count]->setvariance((double)beamfluxvariance,dataIndex);
                    }
                }
            } // for(unsigned dataIndex=0;dataIndex < maxNDataPoints; dataIndex++) {
			order_count++;
        } // if(GetSpectralOrder(order)->gethasSpectralElements() && GetSpectralOrder(order)->gethasSpectralEnergyDistribution()) {
    } // for(int order=(int)minorder; order<=(int)maxorder; order++) {
    
    order_count = 0;
                       
    for(int order=(int)actualMinorder; order<=(int)actualMaxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if(spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            unsigned nDataPointsOfCurrentOrder = spectralEnergyDistribution->getnDataPoints();
            
            for(unsigned dataIndex = 0; dataIndex < nDataPointsOfCurrentOrder; dataIndex++) {
                spectralEnergyDistribution->setfluxData(FluxCalibrationFitVector[order_count]->getflux(dataIndex),dataIndex);
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    beamSED->setfluxData(BeamFluxCalibrationFitVector[beam][order_count]->getflux(dataIndex),dataIndex);
                }
#ifdef PRINT_DEBUG
                cout << order << ' ' << dataIndex << ' '
                << spectralEnergyDistribution->getdistanceData(dataIndex) << ' '
                << spectralEnergyDistribution->getfluxData(dataIndex) << ' '
                << FluxCalibrationFitVector[order_count]->getflux(dataIndex) << ' ';
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    cout << beam << ' ' << beamSED->getfluxData(dataIndex) << ' ' << BeamFluxCalibrationFitVector[beam][order_count]->getflux(dataIndex) << ' ';
                }
                cout << endl;
#endif
            }
#ifdef PRINT_DEBUG
            cout << endl;
#endif
            order_count++;
        }
    }
    
    for(int order=(int)actualMinorder; order<=(int)actualMaxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            spectralEnergyDistribution->populateUncalibratedFluxFromContinuumData();
            
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                beamSED->populateUncalibratedFluxFromContinuumData();
            }
        }
    } //for(int order=(int)minorder; order<=(int)maxorder; order++) {

}

void operaSpectralOrderVector::FitFluxCalibrationAcrossOrders(int lowOrderToClip, int highOrderToClip, int orderBin, bool throughput) {
    
    int usefulMinorder = (int)minorder;
    int usefulMaxorder = (int)maxorder;
    
    bool hasUsefulMinorder = false;
    
    unsigned maxNElements = 0;
    unsigned numberOfBeams = 0;
    
    // Step 1. Figure out maximum number of elements and number of beams
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralEnergyDistribution() &&
            spectralOrder->gethasSpectralElements() &&
            spectralOrder->gethasWavelength()) {
            
            if(!hasUsefulMinorder) {
                usefulMinorder = order;
                hasUsefulMinorder = true;
            }
            
            usefulMaxorder = order;
            
            numberOfBeams = spectralOrder->getnumberOfBeams();
            unsigned nCalibrationElements = spectralOrder->getSpectralEnergyDistribution()->getCalibrationWavelength().size();
            
            if(maxNElements < nCalibrationElements) {
                maxNElements = nCalibrationElements;
            }
        }
    }

    operaFluxVector *FluxCalibrationFitVector[MAXORDERS];
    float *fluxData = new float[MAXORDERS];
    float *orderData = new float[MAXORDERS];
    float *fluxBeamData[MAXNUMBEROFBEAMS];
    
    operaFluxVector *BeamFluxCalibrationFitVector[MAXNUMBEROFBEAMS][MAXORDERS];
    float *orderBeamData[MAXNUMBEROFBEAMS];
    float orderBeamFlux[MAXNUMBEROFBEAMS];
    unsigned npBeam[MAXNUMBEROFBEAMS];
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        npBeam[beam] = 0;
        fluxBeamData[beam] = new float[MAXORDERS];
        orderBeamData[beam] = new float[MAXORDERS];
    }
    
    float fluxvariance = 0;
    float beamfluxvariance = 0;
    
    unsigned order_count = 0;
    
    // Step 2. Calculate robust fit between neighbor orders
    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        
        int loword = order - orderBin;
        int hiord = order + orderBin;
        if(loword < usefulMinorder) {
            loword = usefulMinorder;
        }
        if(hiord > usefulMaxorder) {
            hiord = usefulMaxorder;
        }
        
        int nord = (hiord - loword + 1);
        
        if(!(nord%2)) {
            if(hiord == usefulMaxorder) {
                loword--;
            } else {
                hiord++;
            }
        }
        
        nord = (hiord - loword + 1);
        
        unsigned nElementsOfCurrentOrder = GetSpectralOrder(order)->getSpectralElements()->getnSpectralElements();
        
        FluxCalibrationFitVector[order_count]  = new operaFluxVector(nElementsOfCurrentOrder);
        
        for(unsigned beam=0; beam < numberOfBeams; beam++) {
            BeamFluxCalibrationFitVector[beam][order_count] = new operaFluxVector(nElementsOfCurrentOrder);
        }

        for(unsigned elemIndex=0;elemIndex < maxNElements; elemIndex++) {
            
            if(elemIndex>=nElementsOfCurrentOrder){break;}
            
            unsigned np = 0;
            float orderFlux = NAN;
            
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                npBeam[beam] = 0;
                orderBeamFlux[beam] = NAN;
            }
            
            if (order > lowOrderToClip && order < highOrderToClip) {
                for(int o=loword; o<=hiord; o++) {
                    operaSpectralOrder *spectralOrder = GetSpectralOrder(o);
                    
                    if(spectralOrder->gethasSpectralElements() &&
                       spectralOrder->gethasSpectralEnergyDistribution() &&
                       spectralOrder->gethasWavelength()) {
                        
                        operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                        
                        const operaFluxVector& fluxCalibration = spectralEnergyDistribution->getCalibration(throughput);
                        
                        if(elemIndex < fluxCalibration.getlength() && elemIndex < nElementsOfCurrentOrder) {
                            float flux = (float)fluxCalibration.getflux(elemIndex);
                            fluxvariance = (float)fluxCalibration.getvariance(elemIndex);
                            
                            if ((o > lowOrderToClip && o < highOrderToClip) || o==order) {
                                
                                if(!isnan(flux)) {
                                    orderData[np] = (float)o;
                                    fluxData[np] = flux;
                                    if(o==order) {
                                        orderFlux = flux;
                                    }
                                    np++;
                                }
                                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                                    
                                    const operaFluxVector& BeamFluxCalibration = beamSED->getCalibration(throughput);
                                    
                                    float beamflux = (float)BeamFluxCalibration.getflux(elemIndex);
                                    beamfluxvariance = (float)BeamFluxCalibration.getvariance(elemIndex);
                                    
                                    if(!isnan(beamflux)) {
                                        orderBeamData[beam][npBeam[beam]] = (float)o;
                                        fluxBeamData[beam][npBeam[beam]] = beamflux;
                                        if(o==order) {
                                            orderBeamFlux[beam] = beamflux;
                                        }
                                        npBeam[beam]++;
                                    }
                                }
                            }
                        }
                        
                    }
                }
                if(np && !(np%2)) {
                    np--;
                }
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam] && !(npBeam[beam]%2)) {
                        npBeam[beam]--;
                    }
                }
            } else {
                int o = order;
                operaSpectralOrder *spectralOrder = GetSpectralOrder(o);
                
                if(spectralOrder->gethasSpectralElements() &&
                   spectralOrder->gethasSpectralEnergyDistribution() &&
                   spectralOrder->gethasWavelength()) {
                    
                    operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                    
                    const operaFluxVector& fluxCalibration = spectralEnergyDistribution->getCalibration(throughput);
                    
                    if(elemIndex < fluxCalibration.getlength() && elemIndex < nElementsOfCurrentOrder) {
                        float flux = (float)fluxCalibration.getflux(elemIndex);
                        fluxvariance = (float)fluxCalibration.getvariance(elemIndex);
                        
                        if(!isnan(flux)) {
                            orderData[np] = (float)o;
                            fluxData[np] = flux;
                            orderFlux = flux;
                            np++;
                        }
                        for(unsigned beam=0; beam < numberOfBeams; beam++) {
                            operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                            const operaFluxVector& BeamFluxCalibration = beamSED->getCalibration(throughput);
                            float beamflux = (float)BeamFluxCalibration.getflux(elemIndex);
                            beamfluxvariance = (float)BeamFluxCalibration.getvariance(elemIndex);
                            if(!isnan(beamflux)) {
                                orderBeamData[beam][npBeam[beam]] = (float)o;
                                fluxBeamData[beam][npBeam[beam]] = beamflux;
                                orderBeamFlux[beam] = beamflux;
                                npBeam[beam]++;
                            }
                        }
                    }
                }
            }
            
            if(np>=3) {
                float am,bm,abdevm;
                ladfit(orderData,fluxData,np,&am,&bm,&abdevm); // robust linear fit: f(x) =  a + b*x 
                
                FluxCalibrationFitVector[order_count]->setflux(double(bm*(double)order + am),elemIndex);
                FluxCalibrationFitVector[order_count]->setvariance(double(0.674433*abdevm)*double(0.674433*abdevm),elemIndex);
            } else {
                FluxCalibrationFitVector[order_count]->setflux((double)orderFlux,elemIndex);
                FluxCalibrationFitVector[order_count]->setvariance((double)fluxvariance,elemIndex);
            }
            
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                if(npBeam[beam]>=3) {
                    float amBeam,bmBeam,abdevmBeam;
                    ladfit(orderBeamData[beam],fluxBeamData[beam],npBeam[beam],&amBeam,&bmBeam,&abdevmBeam); // robust linear fit: f(x) =  a + b*x 
                    BeamFluxCalibrationFitVector[beam][order_count]->setflux(double(bmBeam*(double)order + amBeam),elemIndex);
                    BeamFluxCalibrationFitVector[beam][order_count]->setvariance(double(0.674433*abdevmBeam)*double(0.674433*abdevmBeam),elemIndex);
                } else {
                    BeamFluxCalibrationFitVector[beam][order_count]->setflux((double)orderBeamFlux[beam],elemIndex);
                    BeamFluxCalibrationFitVector[beam][order_count]->setvariance((double)beamfluxvariance,elemIndex);
                }
            }
        } //  for(unsigned elemIndex=0;elemIndex < maxNElements; elemIndex++) {
        order_count++;
    } // for(int order=(int)minorder; order<=(int)maxorder; order++) {
    
    order_count = 0;
    // Step 3. Feed back fit quantities to order SED
    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if(spectralOrder->gethasSpectralElements() &&
           spectralOrder->gethasSpectralEnergyDistribution() && 
           spectralOrder->gethasWavelength()) {
            
            for(unsigned b=0; b < spectralOrder->MainAndBeamCount(); b++) {
				if(b == 0) spectralOrder->MainAndBeamSED(b).setFluxCalibration(*FluxCalibrationFitVector[order_count]);
				else spectralOrder->MainAndBeamSED(b).setFluxCalibration(*BeamFluxCalibrationFitVector[b-1][order_count]);
            }
            
            /*operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            
            operaSpectralElements& fluxCalibrationElements = spectralEnergyDistribution->getCalibrationElements(throughput);
            
            for(unsigned elemIndex=0;elemIndex < fluxCalibrationElements.getnSpectralElements(); elemIndex++) {
                
                fluxCalibrationElements.setFlux(FluxCalibrationFitVector[order_count]->getflux(elemIndex),elemIndex);
                fluxCalibrationElements.setFluxVariance(FluxCalibrationFitVector[order_count]->getvariance(elemIndex),elemIndex);
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    operaSpectralElements& BeamFluxCalibrationElements = beamSED->getCalibrationElements(throughput);
                    BeamFluxCalibrationElements.setFlux(BeamFluxCalibrationFitVector[beam][order_count]->getflux(elemIndex),elemIndex);
                    BeamFluxCalibrationElements.setFluxVariance(BeamFluxCalibrationFitVector[beam][order_count]->getvariance(elemIndex),elemIndex);
                }
            }*/
            order_count++;
        }
    }
}

void operaSpectralOrderVector::getContinuumFluxesForNormalization(double *uncalibratedContinuumFluxForNormalization, double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS],unsigned binsize, int orderBin, unsigned nsigcut) {
    
    /*
     * Note: This function returns the continuum flux measured for a given wavelength obtained from 
     *       spectral energy distribution class.  If the "wavelengthfornormalization" given in the SED  
     *       is not in the range covevered by the orders, then it returns NaNs.
     *
     */
    
    unsigned NumberofBeams = 0;
    double wavelengthForNormalization=0;
    
    for (int order=(int)minorder; order<=(int)maxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralEnergyDistribution() &&
            spectralOrder->gethasSpectralElements() &&
            spectralOrder->gethasWavelength()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            wavelengthForNormalization = spectralEnergyDistribution->getwavelengthForNormalization();
            NumberofBeams = spectralOrder->getnumberOfBeams();
            break;
        }
    }
    
    if(wavelengthForNormalization == 0 || NumberofBeams == 0) {
        //throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		return;	// ths case of a constant DT Aug 9 2013
	}
    
    int *orderWithReferenceFluxForNormalization = new int[MAXORDERS];
    unsigned *elemIndexWithReferenceFluxForNormalization = new unsigned[MAXORDERS];
    
    unsigned nOrdersPicked = getElemIndexAndOrdersByWavelength(orderWithReferenceFluxForNormalization,elemIndexWithReferenceFluxForNormalization,wavelengthForNormalization);
    
    measureContinuumAcrossOrders(binsize,orderBin,nsigcut,nOrdersPicked,orderWithReferenceFluxForNormalization);
    
    // DT May 20 2014 -- notused --int orderpicked = 0;
    unsigned elemIndexPicked = 0;
    
    *uncalibratedContinuumFluxForNormalization = -BIG;
    bool hasfluxforNormalization = false;
    
    for(unsigned i=0;i<nOrdersPicked;i++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(orderWithReferenceFluxForNormalization[i]);
        
        if (spectralOrder->gethasSpectralEnergyDistribution()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            double tmp_uncalibratedContinuumFluxForNormalization = spectralEnergyDistribution->getUncalibratedFlux().getflux(elemIndexWithReferenceFluxForNormalization[i]);
            
            if(*uncalibratedContinuumFluxForNormalization < tmp_uncalibratedContinuumFluxForNormalization) {
                *uncalibratedContinuumFluxForNormalization = tmp_uncalibratedContinuumFluxForNormalization;
                //orderpicked = orderWithReferenceFluxForNormalization[i];
                elemIndexPicked = elemIndexWithReferenceFluxForNormalization[i];
                for(unsigned beam=0; beam < NumberofBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    uncalibratedContinuumBeamFluxForNormalization[beam] = beamSED->getUncalibratedFlux().getflux(elemIndexPicked);
                }
                hasfluxforNormalization = true;
            }
        }
    }
    
    if(!hasfluxforNormalization) {
        *uncalibratedContinuumFluxForNormalization = NAN;
        for(unsigned beam=0; beam < NumberofBeams; beam++) {
            uncalibratedContinuumBeamFluxForNormalization[beam] = NAN;
        }
    }
    
    delete[] orderWithReferenceFluxForNormalization;
    delete[] elemIndexWithReferenceFluxForNormalization;
}

unsigned operaSpectralOrderVector::getLEElementCount(string LEfluxCalibration) {
	ifstream astream;
	string dataline;
	unsigned length = 0;
    
	astream.open(LEfluxCalibration.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
                    length++;  
                }
            }
		} // while (astream.good())
		astream.close();
	}	// if (astream.open()
	return length;
}

void operaSpectralOrderVector::readLEFluxCalibration(string LEfluxCalibration, operaSpectralElements *fluxCalibrationElements) {
	ifstream astream;
	string dataline;
	unsigned length = 0;
	double wl, fcal;
    
	astream.open(LEfluxCalibration.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &wl, &fcal);
					fluxCalibrationElements->setwavelength(wl, length);
					fluxCalibrationElements->setFlux(fcal, length);
					fluxCalibrationElements->setHasWavelength(true);
                    length++;  
                }
            }
		} // while (astream.good())
		astream.close();
	}	// if (astream.open()
}


/* Notes: This function returns a clean uniform sample of the continuum flux for both the main elements and the beam elements of all orders.
 *        The spectral regions for the continuum is obtained from an input mask inputWavelengthMaskForUncalContinuum.
 *        binsize defines how many points should be binned to get rid of noise/features
 *        delta_wl defines which is the minimum wavelength range to merge data
 *        The size of final sample is given by numberOfPointsInUniformSample
 *        The wavelength vector and flux vectors of the sample are returned
 *        This function uses a number of routines in the operaSpectralTools library */
operaSpectrum operaSpectralOrderVector::calculateCleanUniformSampleOfContinuum(int Minorder, int Maxorder, unsigned binsize, double delta_wl, string inputWavelengthMaskForUncalContinuum, unsigned numberOfPointsInUniformSample, bool useBeams) {

    // Loop over orders to set maximum number of beams
    unsigned numberOfBeams = getNumberOfBeams(Minorder, Maxorder);
    if(!useBeams) numberOfBeams = 0;
    
    // Collect continuum sample using input mask
    operaSpectrum uncal(numberOfBeams+1, 0);
    double minwl = BIG;
    double maxwl = -BIG;
    
    operaWavelengthRanges rangesInWLMask = readContinuumWavelengthMask(inputWavelengthMaskForUncalContinuum);
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            for(unsigned samplestart=0; samplestart<spectralElements->getnSpectralElements(); samplestart+=binsize){
                unsigned sampleend = samplestart + binsize;
                if (sampleend > spectralElements->getnSpectralElements()) {
                    sampleend = spectralElements->getnSpectralElements();
                    if(binsize <= sampleend) samplestart = sampleend - binsize;
                }
                
                operaSpectrum sample(uncal.fluxcount(), 0);
                for(unsigned elemIndex=samplestart;elemIndex < sampleend; elemIndex++) {
					if(rangesInWLMask.contains(spectralElements->getwavelength(elemIndex))) {
						if(spectralElements->getFlux(elemIndex) && !isnan((float)spectralElements->getFlux(elemIndex))) {
							if(useBeams) sample.insert(spectralElements->getwavelength(elemIndex), spectralOrder->getMainAndBeamFluxes(elemIndex));
							else sample.insert(spectralElements->getwavelength(elemIndex), spectralElements->getFlux(elemIndex));
							if(sample.lastwl() < minwl) minwl = sample.lastwl();
							if(sample.lastwl() > maxwl) maxwl = sample.lastwl();
						}
					}
                }
                
                if(sample.size() > MINNUMBEROFPOINTSINSIDEBIN) uncal.insert(sample.meanwl(), sample.medianfluxes());
            }
        }
    }
    uncal.sort(); // Sort sample
    
    // Merge data from all orders
    operaSpectrum merged(uncal.fluxcount(), 0);
    
    operaSpectrum temp(uncal.fluxcount(), 0);
    for(unsigned index=0; index<uncal.size(); index++) {
        temp.insert(uncal.getwavelength(index), uncal.getfluxes(index));
        
        if(index >= uncal.size()-1 || fabs(uncal.getwavelength(index+1) - uncal.getwavelength(index)) > delta_wl) {
			operaVector meanfluxes = temp.meanfluxes();
			if(merged.size() == 0) merged.insert(minwl, meanfluxes); //insert minwl at beginning
			merged.insert(temp.meanwl(), meanfluxes);
			if(merged.size() == uncal.size()+1) merged.insert(maxwl, meanfluxes); //insert maxwl at end
			temp.clear();
        }
    }
    
    // Calculate an uniform sample. This is a necessary step in order to make further spline interpolations to work.
    return calculateUniformSample(merged, numberOfPointsInUniformSample);
}

// Version of the function supporting float arrays. Still called like this in a few places, so easier to just keep this overload for now.
void operaSpectralOrderVector::calculateCleanUniformSampleOfContinuum(int Minorder, int Maxorder, unsigned binsize, double delta_wl, string inputWavelengthMaskForUncalContinuum, unsigned numberOfPointsInUniformSample, float *uniform_wl, float *uniform_flux,float *uniform_Beamflux[MAXNUMBEROFBEAMS], bool useBeams) {
	// Call the real version of the function
    operaSpectrum uniform = calculateCleanUniformSampleOfContinuum(Minorder, Maxorder, binsize, delta_wl, inputWavelengthMaskForUncalContinuum, numberOfPointsInUniformSample, useBeams);
    
    // Then copy to out parameters
    for(unsigned i = 0; i < numberOfPointsInUniformSample; i++) {
		uniform_wl[i] = (float) uniform.getwavelength(i);
		uniform_flux[i] = (float) uniform.getflux(i);
		for(unsigned b = 1; b < uniform.fluxcount(); b++) uniform_Beamflux[b-1][i] = (float) uniform.fluxvector(b)[i];
	}
}

void operaSpectralOrderVector::trimOrdersByWavelengthRanges(int Minorder, int Maxorder) {
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelengthRange()) {
			spectralOrder->TrimOrderToWavelengthRange();
        }
    }
}

void operaSpectralOrderVector::applyTelluricRVCorrectionINTOExtendendSpectra(double tellCorrection, int Minorder, int Maxorder) {
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralOrder->applyRvelWavelengthCorrection(tellCorrection);
            spectralElements->copyTOtell();	// Save the tell
            operaWavelength *wavelength = spectralOrder->getWavelength();
            spectralElements->setwavelengthsFromCalibration(wavelength); // Reload the non-shifted wavelengths
            spectralElements->setHasWavelength(true);
        }
    }
}

void operaSpectralOrderVector::setRVCorrectionINTOExtendendSpectra(double rvelCorrection, int Minorder, int Maxorder) {
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralOrder->setExtendedRvelWavelengthCorrection(rvelCorrection);
            spectralElements->setHasWavelength(true);
        }
    }
}

void operaSpectralOrderVector::correctFlatField(int Minorder, int Maxorder, bool StarPlusSky, bool starplusskyInvertSkyFiber, double skyOverStarFiberAreaRatio) {
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            spectralOrder->divideSpectralElementsBySEDElements(true);
			if (StarPlusSky) spectralOrder->calculateStarAndSkyElements(starplusskyInvertSkyFiber, skyOverStarFiberAreaRatio);
        } else {
            spectralOrder->sethasSpectralElements(false);
        }
    }
}

// Saves the flux into both the normalized flux and fcal flux of the extended spectra
void operaSpectralOrderVector::saveFluxINTOExtendedSpectra(int Minorder, int Maxorder) {
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
			spectralOrder->CopyFluxVectorIntoFcalFlux();
			spectralOrder->CopyFluxVectorIntoNormalizedFlux();
        }
    }
}

// Normalizes the flux order-by-order and saves the result into the normalized flux of the extended spectra
void operaSpectralOrderVector::normalizeFluxOrderbyOrderINTOExtendendSpectra(unsigned normalizationBinsize, int Minorder, int Maxorder, bool normalizeBeams) {
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
			spectralOrder->CopyNormalizedFluxIntoFluxVector();
            spectralOrder->applyNormalization(normalizationBinsize, 0, false, true, normalizeBeams);
            spectralOrder->CopyFluxVectorIntoNormalizedFlux();
        }
    }
}

// Normalizes the flux to the continuum by sampling the continuum across all orders, and saves the result into the normalized flux of the extended spectra
void operaSpectralOrderVector::normalizeFluxINTOExtendendSpectra(string inputWavelengthMaskForUncalContinuum, unsigned numberOfPointsInUniformSample, unsigned normalizationBinsize, double delta_wl, int Minorder, int Maxorder, bool normalizeBeams) {
    
    if (inputWavelengthMaskForUncalContinuum.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    for (int order=Minorder; order<=Maxorder; order++) GetSpectralOrder(order)->CopyNormalizedFluxIntoFluxVector(); //Start off by copying in normalized flux
    
    operaSpectrum continuumSample = calculateCleanUniformSampleOfContinuum(Minorder, Maxorder, normalizationBinsize, delta_wl, inputWavelengthMaskForUncalContinuum, numberOfPointsInUniformSample, normalizeBeams);
    
	for (int order=Minorder; order<=Maxorder; order++) {
		operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
		if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
			spectralOrder->fitSEDUncalibratedFluxToSample(continuumSample); //Fit the SED uncalibrated flux to our uniform continuum sample
			spectralOrder->applyNormalizationFromExistingContinuum(normalizeBeams); //Divide flux vector by the SED uncalibrated flux
			spectralOrder->CopyFluxVectorIntoNormalizedFlux(); //Don't forget to copy back to normalized flux once we're done
		}
    }
}

void operaSpectralOrderVector::normalizeAndCalibrateFluxINTOExtendendSpectra(string inputWavelengthMaskForUncalContinuum, double exposureTime, bool AbsoluteCalibration, unsigned numberOfPointsInUniformSample, unsigned normalizationBinsize, double delta_wl, int Minorder, int Maxorder, bool normalizeBeams, double SkyOverStarFiberAreaRatio, bool StarPlusSky) {
    
    if (inputWavelengthMaskForUncalContinuum.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    unsigned NumberofBeams = getNumberOfBeams(Minorder, Maxorder);
    
    // For StarPlusSky, we copy flux calibration and throughput (well actually flux calibration twice) from SED into BeamSED only of sky fibers?
    if(StarPlusSky) {
        for (int order=Minorder; order<=Maxorder; order++) {
            operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
            if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
                
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                operaFluxVector fluxCalibrationFluxVector = spectralEnergyDistribution->getFluxCalibration();
                operaFluxVector thruputFluxVector = spectralEnergyDistribution->getFluxCalibration(); //These are the same?
                
                for(unsigned beam=0; beam < NumberofBeams; beam++) {
                    if(beam*2 >= NumberofBeams) {// Sky Fiber
                        operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);
                        BeamSED->setFluxCalibration(fluxCalibrationFluxVector);
                        BeamSED->setThroughput(thruputFluxVector);
                    }
                }
            }
        }
    }
    
    // Don't normalize beams for StarPlusSky?
    if(StarPlusSky) normalizeBeams = false;
    
    //Below is same as normalize function, but we reuse the continuum sample
    for (int order=Minorder; order<=Maxorder; order++) GetSpectralOrder(order)->CopyNormalizedFluxIntoFluxVector(); //Start off by copying in normalized flux
    
    operaSpectrum continuumSample = calculateCleanUniformSampleOfContinuum(Minorder, Maxorder, normalizationBinsize, delta_wl, inputWavelengthMaskForUncalContinuum, numberOfPointsInUniformSample, normalizeBeams);
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
			spectralOrder->fitSEDUncalibratedFluxToSample(continuumSample); //Fit the SED uncalibrated flux to our uniform continuum sample
			spectralOrder->applyNormalizationFromExistingContinuum(normalizeBeams); //Divide flux vector by the SED uncalibrated flux
			spectralOrder->CopyFluxVectorIntoNormalizedFlux(); //Don't forget to copy back to normalized flux once we're done
		}
    }
    
    //Below we apply flux calibrations
    for (int order=Minorder; order<=Maxorder; order++) GetSpectralOrder(order)->CopyFcalFluxIntoFluxVector(); //Start off by copying in fcal flux
    
    double wavelengthForNormalization=0;
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            if(spectralEnergyDistribution->getwavelengthForNormalization() > 0) {
                wavelengthForNormalization = spectralEnergyDistribution->getwavelengthForNormalization();
                break;
            }
        }
    }
    
    operaVector spectralBinConstants(NumberofBeams+1);
    if(AbsoluteCalibration) spectralBinConstants = exposureTime;
    else spectralBinConstants = getFluxesAtWavelength(continuumSample, wavelengthForNormalization);
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
		if(StarPlusSky && beam*2 >= NumberofBeams) spectralBinConstants[beam+1] /= SkyOverStarFiberAreaRatio; // Sky Fiber
	}

    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            if(AbsoluteCalibration) { // use fluxcalibation
                spectralOrder->multiplySpectralElementsBySEDElements(false, spectralBinConstants);
            } else { // use throughput
                spectralOrder->multiplySpectralElementsBySEDElements(true, spectralBinConstants);
            }
            spectralOrder->CopyFluxVectorIntoFcalFlux();
        }
    }
}

// Applies a flat response calibration to the flux and saves the result into the fcal flux of the extended spectra
void operaSpectralOrderVector::applyFlatResponseINTOExtendendSpectra(string flatResponse, bool FITSformat, int Minorder, int Maxorder) {
    
    if (flatResponse.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }

    unsigned NumberofBeams = getNumberOfBeams(Minorder, Maxorder);
    
    operaVector spectralBinConstants(NumberofBeams+1);
    spectralBinConstants = 1.0;
    
    for (int order=Minorder; order<=Maxorder; order++) GetSpectralOrder(order)->CopyFcalFluxIntoFluxVector(); //Start off by copying in fcal flux
    
    readFlatResponseIntoSED(flatResponse,Minorder,Maxorder,FITSformat); // Read flat response data into Spectral Energy Distribution class

    // The flux calibration is effectively calculated below, but it expects that SED elements are in there.
    bool AbsoluteCalibration = false;
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            if(AbsoluteCalibration) { // use fluxcalibation
                spectralOrder->multiplySpectralElementsBySEDElements(false, spectralBinConstants);
            } else { // use throughput
                spectralOrder->multiplySpectralElementsBySEDElements(true, spectralBinConstants);
            }
            spectralOrder->CopyFluxVectorIntoFcalFlux();
        }
    }
}

// Applies the removal of the continuum polarization and saves the result into the continuum removed polarization
void operaSpectralOrderVector::removeContinuumPolarization(int Minorder, int Maxorder) {
	for (int order=Minorder; order<=Maxorder; order++) {
		operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
		operaPolarimetry *polarimetry = spectralOrder->getPolarimetry();
		for(unsigned i = 0; i < 4; i++) {
			stokes_parameter_t param = static_cast<stokes_parameter_t> (i);
			if(polarimetry->getHasDegreeOfStokes(param)) {
				const operaVector& wl = spectralOrder->getSpectralElements()->getWavelength();
				operaFluxVector temp = polarimetry->getDegreeOfPolarization(param) + 1.0;
				operaFluxVector yfit = fitSpectrumToPolynomial(wl, temp.getflux(), wl, 2);
				temp /= yfit;
				temp -= 1.0;
				polarimetry->setContinuumRemoved(param, temp);
			}
		}
	}
}

unsigned operaSpectralOrderVector::getMaxNumberOfElementsInOrder(int Minorder, int Maxorder) const {
    unsigned maxNElements = 0;
    for (int order=Minorder; order<=Maxorder; order++) {
        const operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            const operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            if(maxNElements < spectralElements->getnSpectralElements()) {
                maxNElements = spectralElements->getnSpectralElements();
            }
        }
    }
    return maxNElements;
}

unsigned operaSpectralOrderVector::getNumberOfBeams(int Minorder, int Maxorder) const {
    unsigned NumberofBeams = 0;
    for (int order=Minorder; order<=Maxorder; order++) {
        const operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            if(NumberofBeams==0 && spectralOrder->getnumberOfBeams()) {
                NumberofBeams = spectralOrder->getnumberOfBeams();
                break;
            }
        }
    }
    return NumberofBeams;
}

operaSpectrum operaSpectralOrderVector::getSpectrumWithinTelluricMask(string inputWavelengthMaskForTelluric, int Minorder, int Maxorder, bool normalized, unsigned normalizationBinsize) {
    operaWavelengthRanges telluricwlmask = readContinuumWavelengthMask(inputWavelengthMaskForTelluric);
    operaSpectrum spectrum;
    for (int order=Maxorder; order>=Minorder; order--) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            //weird way to do this...
            spectralOrder->getSpectralElements()->createExtendedVectors();
            spectralOrder->getSpectralElements()->copyTOrawFlux();
            if(normalized && normalizationBinsize>0) {
                spectralOrder->applyNormalization(normalizationBinsize,0,false,true,false);
            }
            
            operaWavelength *wavelength =  spectralOrder->getWavelength();
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralElements->setwavelengthsFromCalibration(wavelength);
            
            for (unsigned elemIndex=0; elemIndex<spectralElements->getnSpectralElements(); elemIndex++) {
				double wl = spectralElements->getwavelength(elemIndex);
				double flux = spectralElements->getFlux(elemIndex);
				if (telluricwlmask.contains(wl) && !isnan(flux)) {
					spectrum.insert(wl, flux, spectralElements->getFluxVariance(elemIndex));
                }
            }
        }
    }
    spectrum.sort();
    return spectrum;
}

operaSpectralLineList operaSpectralOrderVector::detectSpectralLinesWithinWavelengthMask(string inputWavelengthMaskForTelluric, int Minorder, int Maxorder, bool normalize, unsigned normalizationBinsize, double spectralResolution, bool emissionSpectrum,double LocalMaxFilterWidth,double MinPeakDepth,double DetectionThreshold,double nsigclip) {
    operaWavelengthRanges telluricwlmask = readContinuumWavelengthMask(inputWavelengthMaskForTelluric);
    operaSpectralLineList spectrallines;
    for (int order=Maxorder; order>=Minorder; order--) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            if(normalize && normalizationBinsize>0) {
                spectralOrder->applyNormalization(normalizationBinsize,0,false,true,false);
            }
            
            operaSpectralLineList orderlines = spectralOrder->detectSpectralLines(LocalMaxFilterWidth, MinPeakDepth, DetectionThreshold, nsigclip, spectralResolution, emissionSpectrum);
            appendSpectrumWithinRange(telluricwlmask, orderlines, spectrallines);
        }
    }
    spectrallines.sort();
    return spectrallines;
}

/*
 *  Function to extract normalized spectrum around lines in a given list of lines in sourceLines
 *  E. Martioli Jan 13, 2016
 */
operaSpectrum operaSpectralOrderVector::getSpectrumAroundLines(operaSpectrum sourceLines, int Minorder, int Maxorder, bool normalized, unsigned normalizationBinsize, double spectralResolution,double nsig, double snrClip, unsigned numberOfPointsToCutInOrderEnds) {

    operaWavelengthRanges linesmask = getWavelengthMaskAroundLines(sourceLines,spectralResolution,nsig);
    
    operaSpectrum spectrum;
    
    for (int order=Maxorder; order>=Minorder; order--) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            
            if(normalized && normalizationBinsize>0) {
                spectralOrder->applyNormalization(normalizationBinsize,0,false,true,false);
            }
            
            operaWavelength *wavelength =  spectralOrder->getWavelength();
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralElements->setwavelengthsFromCalibration(wavelength);
            
            for (unsigned elemIndex=numberOfPointsToCutInOrderEnds; elemIndex<spectralElements->getnSpectralElements()-numberOfPointsToCutInOrderEnds; elemIndex++) {
                double wl = spectralElements->getwavelength(elemIndex);
                double flux = spectralElements->getFlux(elemIndex);
                double fluxvar = spectralElements->getFluxVariance(elemIndex);
                
                if (linesmask.contains(wl) && !isnan(flux) && flux/sqrt(fluxvar) > snrClip) {
                    spectrum.insert(wl, flux,fluxvar);
                }
            }
        }
    }
    spectrum.sort();
    return spectrum;
}


/*
 *  Function to extract normalized spectrum within wavelength range
 *  E. Martioli Jan 13, 2016
 */
operaSpectrum operaSpectralOrderVector::getSpectrumWithinWavelengthRange(operaWavelengthRange range, int Minorder, int Maxorder, bool normalized, unsigned normalizationBinsize, double snrClip, unsigned numberOfPointsToCutInOrderEnds) {

    operaSpectrum spectrum;

    for (int order=Maxorder; order>=Minorder; order--) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            
            if(normalized && normalizationBinsize>0) {
                spectralOrder->applyNormalization(normalizationBinsize,0,false,true,false);
            }
            
            operaWavelength *wavelength =  spectralOrder->getWavelength();
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralElements->setwavelengthsFromCalibration(wavelength);
            
            for (unsigned elemIndex=numberOfPointsToCutInOrderEnds; elemIndex<spectralElements->getnSpectralElements()-numberOfPointsToCutInOrderEnds; elemIndex++) {
                double wl = spectralElements->getwavelength(elemIndex);
                double flux = spectralElements->getFlux(elemIndex);
                double fluxvar = spectralElements->getFluxVariance(elemIndex);
                
                if (range.contains(wl) && !isnan(flux) && flux/sqrt(fluxvar) > snrClip) {
                    spectrum.insert(wl, flux,fluxvar);
                }
            }
        }
    }
    spectrum.sort();
    return spectrum;
}

void operaSpectralOrderVector::calculateRawFluxQuantities(int Minorder, int Maxorder, double *integratedFlux, double *meanFlux, double *maxSNR, double *meanSNR) {
    
    unsigned totalNelem = 0;
    *integratedFlux = 0;
    *maxSNR = -BIG;
    double sumSNR = 0;
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralElements->calculateFluxSNR();
            for (unsigned elemIndex=0; elemIndex<spectralElements->getnSpectralElements(); elemIndex++) {
                if(!isnan(spectralElements->getFlux(elemIndex))) {
                    
                    *integratedFlux += spectralElements->getFlux(elemIndex);
                    
                    sumSNR += spectralElements->getFluxSNR(elemIndex);
                    
                    if(spectralElements->getFluxSNR(elemIndex) > *maxSNR) {
                        *maxSNR = spectralElements->getFluxSNR(elemIndex);
                    }
                    
                    totalNelem++;
                }
            }
        }
    }
    if(totalNelem) {
        *meanFlux = *integratedFlux/(double)totalNelem;
        *meanSNR = sumSNR/(double)totalNelem;
    } else {
        *meanFlux = 0;
        *meanSNR = 0;
    }
}

/*
 * unsigned readFITSFlatResponse(string filename,float *frwavelength,float *flatresp)
 * \brief Read OPERA FITS flat response from the .fits.gz file.
 */
operaSpectrum operaSpectralOrderVector::readFITSFlatResponse(string filename) {
    operaFITSImage inFlatResp(filename, tfloat, READONLY);
    unsigned np = inFlatResp.getXDimension();
    if(np > MAXNUMBEROFPOINTSINFLATRESPONSE) {
        throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    operaSpectrum fr;
    for (unsigned i=0; i<np; i++) fr.insert(inFlatResp[unsigned(0)][i], inFlatResp[1][i]);
    inFlatResp.operaFITSImageClose();
    return fr;
}

/*
 * unsigned readLibreEspritFlatResponse(string filename,float *frwavelength,float *flatresp)
 * \brief Read Libre-Esprit flat response from the .s file.
 */
operaSpectrum operaSpectralOrderVector::readLibreEspritFlatResponse(string filename) {
    operaSpectrum fr;
    operaistream fspectrum(filename.c_str());
    if (fspectrum.is_open()) {
        string dataline;
        unsigned line = 0;
        float w, fresp;
        while (fspectrum.good()) {
            getline(fspectrum, dataline);
            if (strlen(dataline.c_str())) {
                if (dataline.c_str()[0] == '*') {
                    // skip comments
                } else if (line == 0) {
					unsigned temp;
                    sscanf(dataline.c_str(), "%u", &temp);
                    line++;
                } else {
                    sscanf(dataline.c_str(), "%f %f", &w, &fresp);
                    fr.insert(w, fresp);
                    if(fr.size() > MAXNUMBEROFPOINTSINFLATRESPONSE) {
                        throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
                    }
                    line++;
                }
            }
        }
        fspectrum.close();
    }
    return fr;
}

/*
 * void readLibreEspritFlatResponseIntoSED(string filename,int Minorder, int Maxorder, bool FITSformat)
 * \brief Read flat response from FITS file or from LE .s file (when FITSformat=false).
 */
void operaSpectralOrderVector::readFlatResponseIntoSED(string filename,int Minorder, int Maxorder, bool FITSformat) {
    // Read in flat response file
    operaSpectrum rawflatresp;
    if (FITSformat) rawflatresp = readFITSFlatResponse(filename);
    else rawflatresp = readLibreEspritFlatResponse(filename);
    
    // Resample to a uniform wavelength distribution with the same number of elements
    operaSpectrum flatresp = calculateUniformSample(rawflatresp, rawflatresp.size());
    
    // Get wavelength where flat response is closest to 1.0.
    double wavelengthForNormalization = flatresp.getwavelength(MinIndex(Abs(flatresp.fluxvector() - 1.0)));
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
			// Set wavelengthForNormalization
			spectralOrder->getSpectralEnergyDistribution()->setwavelengthForNormalization(wavelengthForNormalization);
			// Use the flat response to calculate the FluxCalibration and InstrumentThroughput of the SEDs
            spectralOrder->fitSEDFcalAndThroughputToFlatResp(flatresp);
        }
    }
}

/*
 *  Returns an operaSpectrum containing the specified parts of the extended spectrum orders based on input parameters.
 *
 *  E. Martioli Jan 26, 2016
 *
 *  This function is useful when one wants to quickly grab a spectrum out of a .spc product. 
 *  One can select the wavelength type: wavelength calibration, telluric  corrected,
 *  and/or heliocentric corrected. One can also input an additional source radial velocity
 *  correction (sourceRV_KMS), which will subtracted from the wavelength scale.
 *  The flux type may also be selected through the variable RawNormalizedOrCalibrated,
 *  where one can input a signal-to-noise cut threshold (snrClip) and also clip the
 *  order ends through the variable numberOfPointsToCutInOrderEnds.
 */
operaSpectrum operaSpectralOrderVector::getExtendedSpectrum(int Minorder, int Maxorder, unsigned RawNormalizedOrCalibrated, bool wavelengthTelluricCorrected, bool wavelengthHelioCorrected, double snrClip, unsigned numberOfPointsToCutInOrderEnds, double sourceRV_KMS) {
    
    operaSpectrum spectrum;
    
    for (int order=Maxorder; order>=Minorder; order--) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            
            // In case it doesn't have extended data then set values to avoid using the extended data
            if(!spectralElements->getHasExtendedBeamFlux()) {
                RawNormalizedOrCalibrated = 0;
                wavelengthTelluricCorrected = false;
                wavelengthHelioCorrected = false;
            }
            
            for (unsigned elemIndex=numberOfPointsToCutInOrderEnds; elemIndex<spectralElements->getnSpectralElements()-numberOfPointsToCutInOrderEnds; elemIndex++) {
                double wl = 0;
                double flux = 0;
                double fluxvar = 0;
    
                if(wavelengthTelluricCorrected) {
                    wl = spectralElements->gettell(elemIndex);
                } else {
                    wl = spectralElements->getwavelength(elemIndex);
                }
                
                if(wavelengthHelioCorrected) {
                    wl += (spectralElements->getrvel(elemIndex) - sourceRV_KMS*wl/SPEED_OF_LIGHT_KMS);
                }
                
                switch (RawNormalizedOrCalibrated) {
                    case 1:
                        flux = spectralElements->getnormalizedFlux(elemIndex);
                        fluxvar = spectralElements->getnormalizedFluxVariance(elemIndex);
                        break;
                    case 2:
                        flux = spectralElements->getfcalFlux(elemIndex);
                        fluxvar = spectralElements->getfcalFluxVariance(elemIndex);
                        break;
                    default:
                        flux = spectralElements->getFlux(elemIndex);
                        fluxvar = spectralElements->getFluxVariance(elemIndex);
                        break;
                }
                
                if (flux/sqrt(fluxvar) > snrClip) {
                    spectrum.insert(wl, flux,fluxvar);
                }
            }
        }
    }
    spectrum.sort();
    return spectrum;
}
