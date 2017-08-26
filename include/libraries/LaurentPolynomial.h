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

#ifndef LAURENTPOLYNOMIAL_H
#define LAURENTPOLYNOMIAL_H

#include "libraries/operaLibCommon.h"	// for MAXPOLYNOMIAL

/*! 
 * \sa class LaurentPolynomial
 * \brief Encapsulates a LaurentPolynomial. 
 * \file LaurentPolynomial.h
 * \ingroup libraries
 */
class LaurentPolynomial {
	
private:
    unsigned maxNDataPoints;
    unsigned nDataPoints;
    double *xdataVector;
    double *ydataVector;
    double *yerrorVector;
    
	int minorderOfLaurentPolynomial;    
	int maxorderOfLaurentPolynomial;
    unsigned numberOfCoefficients;
    
	double LaurentPolynomialVector[MAXPOLYNOMIAL];
	double LaurentPolynomialErrors[MAXPOLYNOMIAL];
	double polychisqr;
	
public:
	/*
	 * Constructors / Destructors
	 */
	LaurentPolynomial();	
	LaurentPolynomial(const int minorderOfLaurentPolynomial,const int maxorderOfLaurentPolynomial);
	LaurentPolynomial(const int minorderOfLaurentPolynomial,const int maxorderOfLaurentPolynomial, const double* CoefficientVector);	
	LaurentPolynomial(const int minorderOfLaurentPolynomial,const int maxorderOfLaurentPolynomial, const double* CoefficientVector, const double* CoefficientErrorVector);	
    LaurentPolynomial(const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, const PolynomialCoeffs_t* Coefficients);
    LaurentPolynomial(const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, const doublePolynomialCoeffs_t* Coefficients);
	
	~LaurentPolynomial();
    
	/*!
	 * void setnDataPoints(unsigned NDataPoints)
	 * \brief This function sets the number of data points
	 * \note usage:  setnDataPoints(10);
	 * \param NDataPoints is a const unsigned number of data points
	 * \return void
	 */
    void setnDataPoints(unsigned NDataPoints);
    
	/*!
	 * unsigned getnDataPoints(void)
	 * \brief This function gets the number of data points
	 * \note usage:  n = getnDataPoints();
	 * \return unsigned
	 */
    unsigned getnDataPoints(void) const;
	
	/*!
	 * void createDataVectors(unsigned NDataPoints)
	 * \brief This function creates the data vectors with length NDataPoints
	 * \note usage:  createDataVectors(10);
	 * \param NDataPoints is a const unsigned number of data points
	 * \return void
	 */
	void createDataVectors(unsigned NDataPoints);
    
   	/*!
	 * void deleteDataVectors()
	 * \brief This function deletes and resets all data vectors
	 * \note usage:  deleteDataVectors();
	 * \return void
	 */
	void deleteDataVectors(void);
    
   	/*!
	 * void setDataValues(unsigned index, double Xvalue, double Yvalue, double Yerror)
	 * \brief This function sets data values for a given index position in the data vectors
	 * \note usage:  setDataValues(1, 10.54, 0.33, 0.01);
	 * \return void
	 */
	void setDataValues(unsigned index, double Xvalue, double Yvalue, double Yerror);

    /*!
     * void setDataVectors(unsigned NDataPoints, double *Xvector, double *Yvector, double *Yerrorvector)
     * \brief This function sets the data vectors
     * \note usage:  setDataVectors(80, x, y, err);
     * \return void
     */
    void setDataVectors(unsigned NDataPoints, double *Xvector, double *Yvector, double *Yerrorvector);
    
    double getYerrorValue(unsigned index) const;
    double getYdataValue(unsigned index) const;
    double getXdataValue(unsigned index) const;
    
	/*!
	 * double Get(const unsigned index)
	 * \brief This function gets the LaurentPolynomial value at index.
	 * \note usage: float p3 = Get(3);
	 * \param index is a const unsigned index of the coefficient to get
	 * \return double
	 */
	double Get(const unsigned index) const;
	/*!
	 * void Set(const double x, const unsigned index)
	 * \brief This function sets the LaurentPolynomial value at index.
	 * \note usage: Set((double)x, 3);
	 * \param x is a double value of the LaurentPolynomial coefficient at index
	 * \return void
	 */
	void Set(const double x, const unsigned index);	
	/*!
	 * double Evaluate(double x)
	 * \brief This function returns the value of a given LaurentPolynomial function.
	 * \note usage: double value = Evaluate((double)x);
	 * \param x is a double input value for which the given LaurentPolynomial is evaluated
	 * \return double value of evaluation of the LaurentPolynomial
	 */
	double Evaluate(const double x) const;
	/*!
	 * double double* getVector();
	 * \brief This function returns the double *LaurentPolynomialVector.
	 * \note usage: double *vec = getvector();
	 * \return double * - the vector of LaurentPolynomial order coefficients
	 */
	double* getVector();	
	/*!
	 * double double* getErrorVector();
	 * \brief This function returns the double *LaurentPolynomialErrorVector.
	 * \note usage: double *errs = getvector();
	 * \return double * - the vector of LaurentPolynomial coefficient errors
	 */
	double* getErrorVector();	
	/*!
	 * void setMinMaxOrderOfLaurentPolynomial(int MinOrder,int MaxOrder);
	 * \brief This function sets the int minimum and maximum orders of LaurentPolynomial.
	 * \note usage: setMinMaxOrderOfLaurentPolynomial(3);
	 * \return void
	 */
	void setMinMaxOrderOfLaurentPolynomial(int MinOrder,int MaxOrder);
    /*!
	 * unsigned getMinorderOfLaurentPolynomial(void);
	 * \brief This function returns the int minimum order of LaurentPolynomial.
	 * \note usage: int npar = getMinorderOfLaurentPolynomial<float>();
	 * \note usage: int npar = getMinorderOfLaurentPolynomial();
	 */
	int getMinorderOfLaurentPolynomial(void) const;
	/*!
	 * unsigned getMaxorderOfLaurentPolynomial(void);
	 * \brief This function returns the int maximum order of LaurentPolynomial.
	 * \note usage: int npar = getMaxorderOfLaurentPolynomial<float>();
	 * \note usage: int npar = getMaxorderOfLaurentPolynomial();
	 */
	int getMaxorderOfLaurentPolynomial(void) const;
    
	/*!
	 * void getNumberOfCoefficients(void);
	 * \brief This function returns the number of coefficients.
	 * \note usage: getNumberOfCoefficients();
	 * \return void
	 */
	unsigned getNumberOfCoefficients(void) const;
    
	double getCoefficient(unsigned index) const;
    
	/*!
	 * void setCoefficient(unsigned index, double value);
	 * \brief This function sets a coefficent at the index.
	 * \note usage: setOrderOfLaurentPolynomial(3);
	 * \return void
	 */
	void setCoefficient(unsigned index, double value);	
	/*!
	 * void double getCoefficientError(unsigned index);
	 * \brief This function gets a coefficent error at the index.
	 * \note usage: double err = getCoefficientError(3);
	 * \return double
	 */
	double getCoefficientError(unsigned index) const;
	/*!
	 * void setCoefficientError(unsigned index, double value);
	 * \brief This function sets a coefficent error at the index.
	 * \note usage: setCoefficientError(3, 0.001);
	 * \return void
	 */
	void setCoefficientError(unsigned index, double value);	
	/*!
	 * struct LaurentPolynomialCoeffs_t* getLaurentPolynomialCoeffs();
	 * \brief This function returns a PolynomialCoeffs_t struct.
	 * \note usage: PolynomialCoeffs_t *p = getLaurentPolynomialCoeffs<float>();
	 * \note allocates storage that must be freed
	 * \return PolynomialCoeffs_t  * - the PolynomialCoeffs_t struct *
	 */
	PolynomialCoeffs_t* getLaurentPolynomialCoeffs();
	/*!
	 * void setLaurentPolynomialCoeffs(PolynomialCoeffs_t* coeffs);
	 * \brief This function sets a PolynomialCoeffs_t struct.
	 * \return void
	 */
	 void setLaurentPolynomialCoeffs(PolynomialCoeffs_t* coeffs);	
	/*!
	 * void setChisqr(double Chisqr)
	 * \brief This function sets chisqr.
	 * \note usage: setChisqr(0.98);
	 * \return void
	*/
	void setChisqr(double Chisqr);	
	/*!
	 * double getChisqr(void)
	 * \brief This function gets chisqr.
	 * \note usage: double c = getChisqr();
	 * \return double
	 */
	double getChisqr(void) const;
    /*
     * double calculateRMSofResiduals(void)
     * \brief This function calculates the root mean square of residuals.
     * \brief usage: double rms = calculateRMSofResiduals();
     * \return double
     */
    double calculateRMSofResiduals(void) const;
    
	/*!
	 * void printEquation(ostream *pout)
	 * \brief This function prints the LaurentPolynomial equation in Gnuplot format.
	 * \note usage: printEquation(&cout);
	 * \return void
	 */
	void printEquation(ostream *pout);
    
	/*!
	 * void FitModeltoData(ostream *pout)
     * \brief This function performs a least-squares fit of the Laurent Polynomial model to the data.
     * \brief usage: FitModeltoData();
	 * \return void
	 */
    void FitModeltoData(void);
 
    /*!
     * void removeOutLiersFromDataSet(unsigned binsize, unsigned nsig)
     * \brief This function eliminate outliers
     * \brief usage: cleanOutLiers();
     * \return void
     */
    void removeOutLiersFromDataSet(unsigned binsize, float nsig);
    
};

/*!
 * SimpletType EvaluateLaurentPolynomialQuick(double x, LaurentPolynomial &pol)
 * \brief This function returns the value of a given LaurentPolynomial function.
 * \note usage: float value = Evaluate<float>((float)x);
 * \note usage: double value = Evaluate((double)x);
 * \note x is a double input value for which the given LaurentPolynomial is evaluated
 * \return double value of evaluation of the LaurentPolynomial
 */
static inline double EvaluateLaurentPolynomialQuick(const double x, LaurentPolynomial &poly) {
	double *LaurentPolynomialVector = poly.getVector();
	double fpoly = LaurentPolynomialVector[0];
	for (int i=poly.getMinorderOfLaurentPolynomial(); i<=poly.getMaxorderOfLaurentPolynomial(); i++) {
		fpoly += LaurentPolynomialVector[i]*pow(x,i);
	}
	return fpoly;
}

#endif
