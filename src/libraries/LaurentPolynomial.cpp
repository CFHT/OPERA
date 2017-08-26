/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: LaurentPolynomial
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

/*!
 * LaurentPolynomial
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the LaurentPolynomial object.
 * \file LaurentPolynomial.cpp
 * \ingroup libraries
 */

#include <math.h>						// for pow

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/LaurentPolynomial.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaFit.h"

/*
 * Constructors / Destructors
 */
LaurentPolynomial::LaurentPolynomial() :
maxNDataPoints(0),
nDataPoints(0),
xdataVector(NULL),
ydataVector(NULL),
yerrorVector(NULL),
minorderOfLaurentPolynomial(0),
maxorderOfLaurentPolynomial(MAXPOLYNOMIAL-1),
numberOfCoefficients(MAXPOLYNOMIAL),
polychisqr(0.0)
{
	for (unsigned i=0; i<numberOfCoefficients; i++) {
		LaurentPolynomialVector[i] = 0.0;
	}
}

LaurentPolynomial::LaurentPolynomial(const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial) :
maxNDataPoints(0),
nDataPoints(0),
xdataVector(NULL),
ydataVector(NULL),
yerrorVector(NULL),
polychisqr(0.0)
{
    setMinMaxOrderOfLaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial);

	for (unsigned i=0; i<numberOfCoefficients; i++) {
		LaurentPolynomialVector[i] = 0.0;
		LaurentPolynomialErrors[i] = 0.0;
	}
}

LaurentPolynomial::LaurentPolynomial(const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, const double* CoefficientVector) :
maxNDataPoints(0),
nDataPoints(0),
xdataVector(NULL),
ydataVector(NULL),
yerrorVector(NULL),
polychisqr(0.0)
{
    setMinMaxOrderOfLaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial);
    
	for (unsigned i=0; i<numberOfCoefficients; i++) {
		LaurentPolynomialVector[i] = CoefficientVector[i];
	}
}

LaurentPolynomial::LaurentPolynomial(const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, const double* CoefficientVector, const double* CoefficientErrorVector) :
maxNDataPoints(0),
nDataPoints(0),
xdataVector(NULL),
ydataVector(NULL),
yerrorVector(NULL),
polychisqr(0.0)
{
    setMinMaxOrderOfLaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial);

	for (unsigned i=0; i<numberOfCoefficients; i++) {
		LaurentPolynomialVector[i] = CoefficientVector[i];
		LaurentPolynomialErrors[i] = CoefficientErrorVector[i];
	}
}

LaurentPolynomial::LaurentPolynomial(const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, const PolynomialCoeffs_t* Coefficients) :
maxNDataPoints(0),
nDataPoints(0),
xdataVector(NULL),
ydataVector(NULL),
yerrorVector(NULL),
polychisqr(0.0)
{
    setMinMaxOrderOfLaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial);
	if(numberOfCoefficients != Coefficients->orderofPolynomial) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
	polychisqr = Coefficients->polychisqr;
	for (unsigned i=0; i<numberOfCoefficients; i++) {
		LaurentPolynomialVector[i] = Coefficients->p[i];
		LaurentPolynomialErrors[i] = Coefficients->e[i];
	}
}

LaurentPolynomial::LaurentPolynomial(const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, const doublePolynomialCoeffs_t* Coefficients) :
maxNDataPoints(0),
nDataPoints(0),
xdataVector(NULL),
ydataVector(NULL),
yerrorVector(NULL),
polychisqr(0.0)
{
    setMinMaxOrderOfLaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial);    
	if(numberOfCoefficients != Coefficients->orderofPolynomial) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
	polychisqr = Coefficients->polychisqr;
	for (unsigned i=0; i<numberOfCoefficients; i++) {
		LaurentPolynomialVector[i] = Coefficients->p[i];
		LaurentPolynomialErrors[i] = Coefficients->e[i];
	}
}

LaurentPolynomial::~LaurentPolynomial() {
    deleteDataVectors();
}

/*!
 * void setnDataPoints(unsigned NDataPoints)
 * \brief This function sets the number of data points
 * \note usage:  setnDataPoints(10);
 * \param NDataPoints is a const unsigned number of data points
 * \return void
 */
void LaurentPolynomial::setnDataPoints(unsigned NDataPoints) {
	// this needs to be checked to prevent possible memory overwrite...
	if (NDataPoints > maxNDataPoints) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	if (NDataPoints == 0) {
		throw operaException("LaurentPolynomial: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
	}
    nDataPoints = NDataPoints;
}

/*!
 * unsigned getnDataPoints(void)
 * \brief This function gets the number of data points
 * \note usage:  n = getnDataPoints();
 * \return unsigned
 */
unsigned LaurentPolynomial::getnDataPoints(void) const {
    return nDataPoints;
}

/*!
 * void createDataVectors(unsigned NDataPoints)
 * \brief This function creates the data vectors with length NDataPoints
 * \note usage:  createDataVectors(10);
 * \param NDataPoints is a const unsigned number of data points
 * \return void
 */
void LaurentPolynomial::createDataVectors(unsigned NDataPoints) {
	if (NDataPoints == 0) {
		throw operaException("LaurentPolynomial: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
	}
    maxNDataPoints = NDataPoints;
    xdataVector = (double *)malloc(sizeof(double)*maxNDataPoints);
    ydataVector = (double *)malloc(sizeof(double)*maxNDataPoints);
    yerrorVector = (double *)malloc(sizeof(double)*maxNDataPoints);
}

/*!
 * void deleteDataVectors()
 * \brief This function deletes and resets all data vectors
 * \note usage:  deleteDataVectors();
 * \return void
 */
void LaurentPolynomial::deleteDataVectors(void) {
    if(xdataVector) {
        free(xdataVector);
        xdataVector = NULL;
    }
    if(ydataVector) {
        free(ydataVector);
        ydataVector = NULL;
    }
    if(yerrorVector) {
        free(yerrorVector);
        yerrorVector = NULL;
    }
    maxNDataPoints = 0;
    nDataPoints = 0;
}

/*!
 * void setDataValues(unsigned index, double Xvalue, double Yvalue, double Yerror)
 * \brief This function sets data values for a given index position in the data vectors
 * \note usage:  setDataValues(1, 10.54, 0.33, 0.01);
 * \return void
 */
void LaurentPolynomial::setDataValues(unsigned index, double Xvalue, double Yvalue, double Yerror) {
	if (index >= maxNDataPoints) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    xdataVector[index] = Xvalue;
    ydataVector[index] = Yvalue;
    yerrorVector[index] = Yerror;
}

/*!
 * void setDataVectors(unsigned NDataPoints, double *Xvector, double *Yvector, double *Yerrorvector)
 * \brief This function sets the data vectors
 * \note usage:  setDataVectors(80, x, y, err);
 * \return void
 */
void LaurentPolynomial::setDataVectors(unsigned NDataPoints, double *Xvector, double *Yvector, double *Yerrorvector) {
	if (NDataPoints > maxNDataPoints) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    nDataPoints = NDataPoints;    
    for(unsigned index=0;index<nDataPoints;index++) {
        xdataVector[index] = Xvector[index];
        ydataVector[index] = Yvector[index];
        yerrorVector[index] = Yerrorvector[index];
    }
}

double LaurentPolynomial::getXdataValue(unsigned index) const {
	if (index >= maxNDataPoints) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    return xdataVector[index];
}

double LaurentPolynomial::getYdataValue(unsigned index) const {
	if (index >= maxNDataPoints) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    return ydataVector[index];
}

double LaurentPolynomial::getYerrorValue(unsigned index) const {
	if (index >= maxNDataPoints) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    return yerrorVector[index];
}

/*
 * double Get(const unsigned index)
 * \brief This function gets the LaurentPolynomial value at index.
 * \brief usage: float p3 = Get<float>(3);
 * \brief usage: double p3 = Get3);
 * \param index is a const unsigned index of the coefficient to get
 * \return void
 */
double LaurentPolynomial::Get(const unsigned index) const {
	if (index >= MAXPOLYNOMIAL) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return LaurentPolynomialVector[index];
}

/*
 * void Set(const double x, const unsigned index)
 * \brief This function sets the LaurentPolynomial value at index.
 * \brief usage: Set<float>((float)x, 3);
 * \brief usage: Set((double)x, 3);
 * \param x is a double value of the LaurentPolynomial coefficient at index
 * \return void
 */
void LaurentPolynomial::Set(const double x, const unsigned index) {
	if (index >= MAXPOLYNOMIAL) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	LaurentPolynomialVector[index] = x;
}

/*
 * SimpletType Evaluate(double x)
 * \brief This function returns the value of a given LaurentPolynomial function.
 * \brief usage: float value = Evaluate<float>((float)x);
 * \brief usage: double value = Evaluate((double)x);
 * \param x is a double input value for which the given LaurentPolynomial is evaluated
 * \return double value of evaluation of the LaurentPolynomial
 */
double LaurentPolynomial::Evaluate(const double x) const {
	double fpoly = 0;
	for (int i=minorderOfLaurentPolynomial; i<=maxorderOfLaurentPolynomial; i++) {
		fpoly += LaurentPolynomialVector[(unsigned)(i-minorderOfLaurentPolynomial)]*pow(x, i);
	}
	return fpoly;
}

/*
 * SimpletType double** getVector();
 * \brief This function returns the double *LaurentPolynomialVector.
 * \brief usage: float *vec = getvector<float>();
 * \brief usage: double *vec = getvector();
 * \return double * - the vector of LaurentPolynomial order coefficients
 */
double* LaurentPolynomial::getVector() {
	return LaurentPolynomialVector;
}

/*
 * double double* getErrorVector();
 * \brief This function returns the double *LaurentPolynomialErrorVector.
 * \brief usage: double *errs = getvector();
 * \return double * - the vector of LaurentPolynomial coefficient errors
 */
double* LaurentPolynomial::getErrorVector() {
	return LaurentPolynomialErrors;
}


/*!
 * void setMinMaxOrderOfLaurentPolynomial(int MinOrder,int MaxOrder);
 * \brief This function sets the int minimum and maximum orders of LaurentPolynomial.
 * \note usage: setMinMaxOrderOfLaurentPolynomial(3);
 * \return void
 */
void LaurentPolynomial::setMinMaxOrderOfLaurentPolynomial(int MinOrder,int MaxOrder) {
	
    minorderOfLaurentPolynomial = MinOrder;
	maxorderOfLaurentPolynomial = MaxOrder;
    
    if(minorderOfLaurentPolynomial <= maxorderOfLaurentPolynomial) {
        numberOfCoefficients = (unsigned)(maxorderOfLaurentPolynomial - minorderOfLaurentPolynomial) + 1;
    } else {
		throw operaException("LaurentPolynomial: maxorder must be greater than minorder", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
}
/*!
 * unsigned getMinorderOfLaurentPolynomial(void);
 * \brief This function returns the int minimum order of LaurentPolynomial.
 * \note usage: int npar = getMinorderOfLaurentPolynomial<float>();
 * \note usage: int npar = getMinorderOfLaurentPolynomial();
 */
int LaurentPolynomial::getMinorderOfLaurentPolynomial(void) const {
    return minorderOfLaurentPolynomial;
}

/*!
 * unsigned getMaxorderOfLaurentPolynomial(void);
 * \brief This function returns the int maximum order of LaurentPolynomial.
 * \note usage: int npar = getMaxorderOfLaurentPolynomial<float>();
 * \note usage: int npar = getMaxorderOfLaurentPolynomial();
 */
int LaurentPolynomial::getMaxorderOfLaurentPolynomial(void) const {
    return maxorderOfLaurentPolynomial;
}

/*!
 * void getNumberOfCoefficients(void);
 * \brief This function returns the number of coefficients.
 * \note usage: getNumberOfCoefficients(3);
 * \return void
 */
unsigned LaurentPolynomial::getNumberOfCoefficients(void) const {
    return numberOfCoefficients;
}

/*
 * void double getCoefficient(unsigned index);
 * \brief This function sets a coefficent at the index.
 * \brief usage: float coeff3 = getCoefficient(3);
 * \return void
 */
double LaurentPolynomial::getCoefficient(unsigned index) const {
	if (index >= MAXPOLYNOMIAL) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return LaurentPolynomialVector[index];
}

/*
 * void setCoefficient(unsigned index, double value);
 * \brief This function sets a coefficent at the index.
 * \brief usage: setOrderOfLaurentPolynomial(3);
 * \return void
 */
void LaurentPolynomial::setCoefficient(unsigned index, double value) {
	if (index >= MAXPOLYNOMIAL) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	LaurentPolynomialVector[index] = value;
}

/*
 * void double getCoefficientError(unsigned index);
 * \brief This function gets a coefficent error at the index.
 * \brief usage: double err = getCoefficientError(3);
 * \return double
 */
double LaurentPolynomial::getCoefficientError(unsigned index) const {
	if (index >= MAXPOLYNOMIAL) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return LaurentPolynomialErrors[index];
}

/*
 * void setCoefficientError(unsigned index, double value);
 * \brief This function sets a coefficent error at the index.
 * \brief usage: setCoefficientError(3, 0.001);
 * \return void
 */
void LaurentPolynomial::setCoefficientError(unsigned index, double value) {
	if (index >= MAXPOLYNOMIAL) {
		throw operaException("LaurentPolynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	LaurentPolynomialErrors[index] = value;
}

/*
 * struct LaurentPolynomialCoeffs_t* getLaurentPolynomialCoeffs();
 * \brief This function returns a PolynomialCoeffs_t struct.
 * \brief usage: PolynomialCoeffs_t *p = getLaurentPolynomialCoeffs<float>();
 * \note allocates storage that must be freed
 * \return PolynomialCoeffs_t  * - the PolynomialCoeffs_t struct *
 */
PolynomialCoeffs_t* LaurentPolynomial::getLaurentPolynomialCoeffs() {
	PolynomialCoeffs_t *pcoefficients = (PolynomialCoeffs_t *)malloc(sizeof(PolynomialCoeffs_t));
	pcoefficients->orderofPolynomial = numberOfCoefficients;
	pcoefficients->polychisqr = polychisqr;
	for (unsigned i=0; i<numberOfCoefficients; i++) {
		pcoefficients->p[i] = LaurentPolynomialVector[i];
		pcoefficients->e[i] = LaurentPolynomialErrors[i];
	}
	return pcoefficients;
}	

/*
 * void setLaurentPolynomialCoeffs();
 * \brief This function returns a PolynomialCoeffs_t struct.
 * \brief usage: LaurentPolynomialCoeffs_t *p = getLaurentPolynomialCoeffs<float>();
 * \note allocates storage that must be freed
 * \return PolynomialCoeffs_t  * - the PolynomialCoeffs_t struct *
 */
void  LaurentPolynomial::setLaurentPolynomialCoeffs(PolynomialCoeffs_t* pcoefficients) {
	numberOfCoefficients = pcoefficients->orderofPolynomial;
	polychisqr = pcoefficients->polychisqr;
	for (unsigned i=0; i<numberOfCoefficients; i++) {
		LaurentPolynomialVector[i] = pcoefficients->p[i];
		LaurentPolynomialErrors[i] = pcoefficients->e[i];
	}
}	

/*
 * void setChisqr(double Chisqr)
 * \brief This function sets chisqr.
 * \brief usage: setChisqr(0.98);
 * \return void
 */

void LaurentPolynomial::setChisqr(double Chisqr) {
	polychisqr = Chisqr;
}

/*
 * double getChisqr(void)
 * \brief This function gets chisqr.
 * \brief usage: double c = getChisqr();
 * \return double
 */
double LaurentPolynomial::getChisqr(void) const {
	return polychisqr;
}

/*
 * double calculateRMSofResiduals(void)
 * \brief This function calculates the root mean square of residuals.
 * \brief usage: double rms = calculateRMSofResiduals();
 * \return double
 */
double LaurentPolynomial::calculateRMSofResiduals(void) const {

	double WSSR = 0;
	
	for(unsigned i=0;i<nDataPoints;i++)
		WSSR += pow((Evaluate(xdataVector[i]) - ydataVector[i]),2.0);
	
	double rms = sqrt(WSSR/(double)(nDataPoints - numberOfCoefficients));
    
	return rms;
}

/*!
 * void printEquation(ostream *pout)
 * \brief This function prints the LaurentPolynomial equation in Gnuplot format.
 * \note usage: printEquation(&cout);
 * \return void
 */
void LaurentPolynomial::printEquation(ostream *pout) {
    if (pout != NULL) {
        *pout << "f(x) =";
        for(int i=minorderOfLaurentPolynomial;i<=maxorderOfLaurentPolynomial;i++) {
            if(i==0) {
                *pout << " + " << LaurentPolynomialVector[(unsigned)(i-minorderOfLaurentPolynomial)];
            } else if (i==1) {
                *pout << " + " << LaurentPolynomialVector[(unsigned)(i-minorderOfLaurentPolynomial)] << "*x";
            } else {
                *pout << " + " << LaurentPolynomialVector[(unsigned)(i-minorderOfLaurentPolynomial)] << "*x**" << i;
            }
        }
        *pout << endl;
    }
}


/*!
 * void FitModeltoData(ostream *pout)
 * \brief This function performs a least-squares fit of the Laurent Polynomial model to the data.
 * \brief usage: FitModeltoData();
 * \return void
 */
void LaurentPolynomial::FitModeltoData(void) {
#ifdef PRINT_DEBUG
    cout <<" Before::polychisqr=" << polychisqr << " ";    
    for(unsigned i=0;i<numberOfCoefficients;i++) {
        cout << " par=" << LaurentPolynomialVector[i] << " " ;
    }
    cout << endl;
#endif
    // The fitting routine only works for negative power, which means
    // minorderOfLaurentPolynomial must be <= 0 and
    // maxorderOfLaurentPolynomial must be = 0, therefore:
    
    if(minorderOfLaurentPolynomial >0 || maxorderOfLaurentPolynomial != 0) {
        throw operaException("LaurentPolynomial: maxorderOfLaurentPolynomial must be 0 for fitting. ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaLMFitLaurentPolynomial(nDataPoints,xdataVector,ydataVector,minorderOfLaurentPolynomial,maxorderOfLaurentPolynomial,LaurentPolynomialVector,&polychisqr);
#ifdef PRINT_DEBUG
    cout <<" After::polychisqr=" << polychisqr << " ";
    for(unsigned i=0;i<numberOfCoefficients;i++) {
        cout << " par=" << LaurentPolynomialVector[i] << " " ;
    }
    cout << endl;
#endif
}

/*!
 * void removeOutLiersFromDataSet(unsigned binsize, float nsig)
 * \brief This function eliminate outliers
 * \brief usage: cleanOutLiers();
 * \return void
 */
void LaurentPolynomial::removeOutLiersFromDataSet(unsigned binsize, float nsig) {
    binsize = (unsigned)floor((float)binsize/2);
    
	if (binsize==0) {
		throw operaException("LaurentPolynomial: binsize=0", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
	}
	if (nDataPoints==0) {
		throw operaException("LaurentPolynomial: nDataPoints=0", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
	}
    
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
            xtmp[np] = (float)xdataVector[ii];
            ytmp[np] = (float)ydataVector[ii];
            np++;
        }

        float am,bm,abdevm;
        
        //--- Robust linear fit
        ladfit(xtmp,ytmp,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
        
        //--- Clean up
        float fitMedianSlope = (bm*(float)xdataVector[i] + am);
        
        if(fabs((float)ydataVector[i] - fitMedianSlope) < nsig*abdevm) {
            cleanxdataVector[numberOfCleanPoints] = xdataVector[i];
            cleanydataVector[numberOfCleanPoints] = ydataVector[i];
            cleanyerrorVector[numberOfCleanPoints] = yerrorVector[i];
            
            numberOfCleanPoints++;
        }
    }

    for(unsigned i=0; i<numberOfCleanPoints; i++) {
        xdataVector[i] = cleanxdataVector[i];
        ydataVector[i] = cleanydataVector[i];
        yerrorVector[i] = cleanyerrorVector[i]; 
    }
    nDataPoints = numberOfCleanPoints;
    
    delete[] cleanxdataVector;
    delete[] cleanydataVector;
    delete[] cleanyerrorVector; 
    
    delete[] xtmp;
    delete[] ytmp;
}


