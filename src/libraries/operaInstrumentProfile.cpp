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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaInstrumentProfile.h"
#include "libraries/operaVectorOperations.h"
#include "libraries/operaException.h"

#include "libraries/operaFit.h"
#include "libraries/operaMath.h"
#include "libraries/operaStats.h"

#ifndef SATURATIONLIMIT
#define SATURATIONLIMIT 65535  // this should be retrieved from the config/param file
#endif

/*! 
 * operaInstrumentProfile
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the Instrument Profile.
 * \file operaInstrumentProfile.cpp
 * \ingroup libraries
 */

using namespace std;

/*! 
 * \brief operaInstrumentProfile
 * \details The instrument profile (IP) consists of a data set (pixelized image) 
 * \details representing the distribution of the fraction of flux for a uniformly illuminated 
 * \details monochromatic image of the entrance slit projected on the spectrograph focal plane. 
 * \details The fiducial set of coordinates chosen is such that the ordinate is oriented along 
 * \details with the dispersion direction and the abscissa with the spatial direction.
 */

/*
 * Constructor
 */

operaInstrumentProfile::operaInstrumentProfile(void) :
NXPoints(1), // number of points in x direction
NYPoints(1), // number of points in y directon
NTotalPoints(1), // total number of points
Xsampling(1),		// (sampling factor for the PSF in x direction)
Ysampling(1),		// (sampling factor for the  PSF in y direction)
xsize(1),			// (size along the x direction in pix units) 
ysize(1),			// (size along the y direction in pix units)
geometricCenterX(0.5),
geometricCenterY(0.5),
nDataPoints(0),
maxnDataPoints(0),
dataCube(0, 0, 0),
ipPolyModel(NYPoints, NXPoints),
chisqrMatrix(NYPoints, NXPoints)
{
	
}

operaInstrumentProfile::operaInstrumentProfile(unsigned ipxsize,unsigned ipxsampling,unsigned ipysize,unsigned ipysampling) :
NXPoints(ipxsize*ipxsampling),
NYPoints(ipysize*ipysampling),
NTotalPoints(ipxsize*ipxsampling*ipysize*ipysampling),
Xsampling(ipxsampling),
Ysampling(ipysampling),
xsize(ipxsize),
ysize(ipysize),
geometricCenterX(0.5),
geometricCenterY(0.5),
nDataPoints(1),
maxnDataPoints(1),
dataCube(1, NYPoints, NXPoints),
distd(1),
ipPolyModel(NYPoints, NXPoints),
chisqrMatrix(NYPoints, NXPoints)
{
	if (NXPoints == 0 || NYPoints == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	setGeometricCenter();
}

operaInstrumentProfile::operaInstrumentProfile(unsigned ipxsize,unsigned ipxsampling,unsigned ipysize,unsigned ipysampling,unsigned NDataPoints) :
NXPoints(ipxsize*ipxsampling),
NYPoints(ipysize*ipysampling),
NTotalPoints(ipxsize*ipxsampling*ipysize*ipysampling),
Xsampling(ipxsampling),
Ysampling(ipysampling),
xsize(ipxsize),
ysize(ipysize),
geometricCenterX(0.5),
geometricCenterY(0.5),
nDataPoints(NDataPoints),
maxnDataPoints(NDataPoints),
dataCube(NDataPoints, NYPoints, NXPoints),
distd(NDataPoints),
ipPolyModel(NYPoints, NXPoints),
chisqrMatrix(NYPoints, NXPoints)
{
	if (NXPoints == 0 || NYPoints == 0 || NDataPoints == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	setGeometricCenter();
}

/*
 * Setter/Getters
 */

unsigned operaInstrumentProfile::getNXPoints(void) const {
	return NXPoints;
}

unsigned operaInstrumentProfile::getNYPoints(void) const {
	return NYPoints;
}

unsigned operaInstrumentProfile::getNTotalPoints(void) const {
	return NTotalPoints;
}

unsigned operaInstrumentProfile::getXsampling(void) const {
	return Xsampling;
}

unsigned operaInstrumentProfile::getYsampling(void) const {
	return Ysampling;
}

unsigned operaInstrumentProfile::getxsize(void) const {
	return xsize;
}

unsigned operaInstrumentProfile::getysize(void) const {
	return ysize;
}

double operaInstrumentProfile::getGeometricCenterX(void) const {
	return geometricCenterX;
}

double operaInstrumentProfile::getGeometricCenterY(void) const {
	return geometricCenterY;
}

void operaInstrumentProfile::setNXPoints(unsigned npx) {
	NXPoints = npx;
}

void operaInstrumentProfile::setNYPoints(unsigned npy) {
	NYPoints = npy;
}

void operaInstrumentProfile::setNTotalPoints(unsigned np) {
	NTotalPoints = np;
}

void operaInstrumentProfile::setXsampling(unsigned xsamp) {
	Xsampling = xsamp;
}

void operaInstrumentProfile::setYsampling(unsigned ysamp) {
	Ysampling = ysamp;
}

void operaInstrumentProfile::setsampling(unsigned xsamp, unsigned ysamp) {
	Xsampling = xsamp;
	Ysampling = ysamp;
}

void operaInstrumentProfile::setxsize(unsigned xs) {
	xsize = xs;
}

void operaInstrumentProfile::setysize(unsigned ys) {
	ysize = ys;
}

void operaInstrumentProfile::setsize(unsigned xs, unsigned ys) {
	xsize = xs;
	ysize = ys;
}

// set x,y coordinates of geometric center (requires x,y size and x,y sampling)
void operaInstrumentProfile::setGeometricCenter(void){
	geometricCenterX = xsize/2.0;
	geometricCenterY = ysize/2.0;	
}

double operaInstrumentProfile::getIPixXCoordinate(unsigned i) const {
    return (i + 0.5)/Xsampling - geometricCenterX;
}

double operaInstrumentProfile::getIPixYCoordinate(unsigned j) const {
    return (j + 0.5)/Ysampling - geometricCenterY;
}

unsigned operaInstrumentProfile::getIPixiIndex(double xcoord) const {
	int i = (int)((xcoord + geometricCenterX)*Xsampling);
	if(i < 0) return 0;
	if(i >= (int)NXPoints) return NXPoints;
	return (unsigned)i;
}

unsigned operaInstrumentProfile::getIPixjIndex(double ycoord) const {
	int j = (int)((ycoord + geometricCenterY)*Ysampling);
	if(j < 0) return 0;
	if(j >= (int)NYPoints) return NYPoints;
	return (unsigned)j;
}

unsigned operaInstrumentProfile::getnDataPoints(void) const {
	return nDataPoints;
}

void operaInstrumentProfile::setnDataPoints(unsigned NDataPoints) {
	if (NDataPoints > maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	nDataPoints = NDataPoints;
}

DMatrix operaInstrumentProfile::getdataCubeValues(unsigned index) const {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return dataCube[index];
}

double operaInstrumentProfile::getdataCubeValues(unsigned i, unsigned j, unsigned index) const {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return dataCube[index][j][i];
}

void operaInstrumentProfile::setdataCubeValues(DMatrix DataMatrix, unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	dataCube[index] = DataMatrix;
}

void operaInstrumentProfile::setdataCubeValues(double DataValue, unsigned i, unsigned j, unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if(index > dataCube.slices()) cout << "slice=" << index << " vs " << dataCube.slices() << endl; 
	else if(j > dataCube.rows()) cout << "j=" << j << " vs " << dataCube.rows() << endl; 
	else if(i > dataCube.cols()) cout << "i=" << i << " vs " << dataCube.cols() << endl;
	dataCube[index][j][i] = DataValue;
}

double operaInstrumentProfile::getdataGivenCoords(double xcoord, double ycoord, unsigned index) const {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return dataCube[index][getIPixjIndex(ycoord)][getIPixiIndex(xcoord)];
}

void operaInstrumentProfile::setdataCubeValues(operaFITSImage &image, operaFITSImage &badpix, double xcenter, double ycenter, unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    for (unsigned j=0; j<NYPoints; j++) {
        unsigned yy = (unsigned)floor(ycenter + getIPixYCoordinate(j));
        for (unsigned i=0; i<NXPoints; i++) {
            unsigned xx = (unsigned)floor(xcenter + getIPixXCoordinate(i));
            if (xx > 0 && xx < image.getnaxis1() && yy > 0 && yy < image.getnaxis2()){
                if (image[yy][xx] < SATURATIONLIMIT && badpix[yy][xx] == 1 && image[yy][xx] > 0) {
                    dataCube[index][j][i] = image[yy][xx];
                } else {
                    dataCube[index][j][i] = NAN;
                }
            }
        }
    }
}

void operaInstrumentProfile::subtractOuterFrame(unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    operaVector xx;
	operaVector yy;
	operaVector z;
	operaVector fz;
	
    unsigned numberofNANs = 0;
    
    for (unsigned j=0; j<NYPoints; j++) {	
        for (unsigned i=0; i<NXPoints; i++) {
            if((j==0 || j==NYPoints-1) && (i==0 || i==NXPoints-1)) {
                if(!isnan(getdataCubeValues(i,j,index))) {
                    xx.insert((double)i);
                    yy.insert((double)j);
                    fz.insert(getdataCubeValues(i,j,index));
                    z.insert(getdataCubeValues(i,j,index));
                } else {
                    numberofNANs++;
                }
            }
        }
    }
    unsigned k = xx.size();
    
    if(k > numberofNANs) {
        unsigned npar = 4;
        double pars[4] = {1,1,1,1};
        double chisqr;
        
        operaLMFit2DPolynomial(k,xx.datapointer(),yy.datapointer(),z.datapointer(),npar,pars,&chisqr);
        
        for (unsigned j=0; j<NYPoints; j++) {
            for (unsigned i=0; i<NXPoints; i++) {
                double value = getdataCubeValues(i,j,index) - Polynomial2DFunction((double)i,(double)j,pars,npar);
                setdataCubeValues(value,i,j,index);
            }
        }
    } else if (k > 3) {
        double minValue = Min(fz);
        for (unsigned j=0; j<NYPoints; j++) {
            for (unsigned i=0; i<NXPoints; i++) {
                double value = getdataCubeValues(i,j,index) - minValue;
                setdataCubeValues(value,i,j,index);
            }
        }
    }
}


void operaInstrumentProfile::normalizeCubeData(unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
#endif
    double sum=0;
    for (unsigned j=0; j<NYPoints; j++) {
        for (unsigned i=0; i<NXPoints; i++) {
            double temp = getdataCubeValues(i,j,index);
            if(!isnan(temp)) sum += temp;
        }
    }
    for (unsigned j=0; j<NYPoints; j++) {
        for (unsigned i=0; i<NXPoints; i++) {
            if(sum) {
                setdataCubeValues(getdataCubeValues(i,j,index)/sum,i,j,index);
            }
        }
    }
}

void operaInstrumentProfile::normalizeCubeData(void) {
	for(unsigned index = 0; index < getnDataPoints(); index++) {
		normalizeCubeData(index);
	}
}

const operaVector& operaInstrumentProfile::getdistdVector(void) const {
	return distd;
}

double operaInstrumentProfile::getdistd(unsigned index) const {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
#endif
	return distd[index];
}

void operaInstrumentProfile::setdistdVector(const operaVector& Distd) {
	distd = Distd;
}

void operaInstrumentProfile::setdistd(double DistdValue, unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	distd[index] = DistdValue;
}

void operaInstrumentProfile::printModel(double DistanceInPixels, unsigned ordernumber, ostream *pout) {
	if (pout) {
		*pout << "#o\ti\tj\tx\ty\td\tz\tzerr\n" << endl;
		for (unsigned j=0; j<NYPoints; j++) {
			for (unsigned i=0; i<NXPoints; i++) {
				*pout << ordernumber << "\t"
                << i << "\t"
                << j << "\t";
				*pout << getIPixXCoordinate(i) <<
				"\t" << getIPixYCoordinate(j) <<
				"\t" << DistanceInPixels <<
                "\t" << getipDataFromPolyModel(DistanceInPixels,i,j) <<
				"\t" << sqrt(getchisqrMatrixValue(i,j)) << endl;
			}
			*pout << endl;
		}
	}
}

void operaInstrumentProfile::setipPolyModel(const PolynomialMatrix& IPPolyModel) {
	ipPolyModel = IPPolyModel;
}

const PolynomialMatrix& operaInstrumentProfile::getipPolyModel(void) {
	return ipPolyModel;
}

const Polynomial& operaInstrumentProfile::getipPolyModelCoefficients(unsigned i,unsigned j) const {
#ifdef RANGE_CHECK
	if (j >= NYPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (i >= NXPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return ipPolyModel[j][i];
}

void operaInstrumentProfile::setipPolyModelCoefficients(const Polynomial& PolyModelCoeffs,unsigned i,unsigned j) {
#ifdef RANGE_CHECK
	if (j >= NYPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (i >= NXPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	ipPolyModel[j][i] = PolyModelCoeffs;
}

DMatrix operaInstrumentProfile::getipDataFromPolyModel(double d) const {
	DMatrix ipmatrix(NYPoints, NXPoints);
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			ipmatrix[j][i] = getipDataFromPolyModel(d,i,j);
		}
	}
	return ipmatrix;
}

double operaInstrumentProfile::getipDataFromPolyModel(double d, unsigned i, unsigned j) const {
	return ipPolyModel[j][i].Evaluate(d);
	//Below might be slightly faster?
	/*unsigned order = ipPolyModel[j][i].getOrderOfPolynomial(); 
	const double* coeffs = ipPolyModel[j][i].getVector();
	double total = 0;
	for(unsigned i = order; i > 0; i--) total = d*total + coeffs[i-1];
	return total;*/
}

void operaInstrumentProfile::setdataCubeFromPolyModel(void) {
	for(unsigned index = 0; index < getnDataPoints(); index++) {
		setdataCubeFromPolyModel(index);
	}
}

void operaInstrumentProfile::setdataCubeFromPolyModel(unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
#endif
    for (unsigned j=0; j<NYPoints; j++) {
        for (unsigned i=0; i<NXPoints; i++) {
            if(!isnan(getdataCubeValues(i,j,index))) {
				setdataCubeValues(getipDataFromPolyModel(distd[index],i,j),i,j,index);
			}
        }
    }			
}

void operaInstrumentProfile::FitPolyMatrixtoIPDataVector(unsigned coeffs, bool witherrors) {
    ipPolyModel = PolynomialMatrix(NYPoints, NXPoints);
	chisqrMatrix = DMatrix(NYPoints, NXPoints);
	
	operaVector par(coeffs);
	operaVector parError(coeffs);
	operaVector ytmp(nDataPoints);
	operaVector ytmpError(nDataPoints);
	operaVector xtmp(nDataPoints);
	
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			
			int nparbestfit = coeffs;
			double bestchisqr = BIG;
			double chisqr;
			
            unsigned npts = 0;
			for(unsigned index=0;index<nDataPoints;index++){	
                if(!isnan(getdataCubeValues(i,j,index))) {
                    ytmp[npts] = getdataCubeValues(i,j,index);				
                    xtmp[npts] = getdistd(index);
                    ytmpError[npts] = 0;
                    npts++;
                }
			}
			
            if(npts < coeffs) {
                coeffs = npts;
            }
            
            if(coeffs) {
                for (unsigned currentfit=1; currentfit<=coeffs; currentfit++) {
                    for	(unsigned k=0; k<currentfit; k++) {
                        par[k] = 1.0;
                    }
                    if (witherrors) {
                        operaMPFitPolynomial(npts, xtmp.datapointer(), ytmp.datapointer(), ytmpError.datapointer(), currentfit, par.datapointer(), parError.datapointer(), &chisqr);
                        if (fabs(chisqr-1.0) < bestchisqr) {
                            bestchisqr = chisqr;
                            nparbestfit = currentfit;
                        }
                        
                    } else {
                        operaLMFitPolynomial(npts, xtmp.datapointer(), ytmp.datapointer(), currentfit, par.datapointer(), &chisqr);
                        if (chisqr < bestchisqr) {                
                            bestchisqr = chisqr;
                            nparbestfit = currentfit;
                        }
                    }
                }
                for	(int k=0; k<nparbestfit; k++) {
                    par[k] = 1.0;
                }
                
                if (witherrors) {
                    operaMPFitPolynomial(npts, xtmp.datapointer(), ytmp.datapointer(), ytmpError.datapointer(), nparbestfit, par.datapointer(), parError.datapointer(), &chisqr);                        
                } else {
                    operaLMFitPolynomial(npts, xtmp.datapointer(), ytmp.datapointer(), nparbestfit, par.datapointer(), &chisqr);
                }
                
                Polynomial pp(nparbestfit, par.datapointer());
                setipPolyModelCoefficients(pp, i,j);
                setchisqrMatrixValue(chisqr,i,j);
                
#ifdef PRINT_DEBUG
                cout << "i="<< i << " j=" << j << " np=" << nparbestfit;
                for	(unsigned k=0; k<nparbestfit; k++) {
                    cout << " p["<<k<<"]=" << par[k];
                }
                cout << " chi2=" << chisqr << endl;	
                
                if(i==(0*NXPoints/8) || i==(1*NXPoints/8) || i==(2*NXPoints/8) || i==(3*NXPoints/8) || i==(4*NXPoints/8) || 
                   i==(5*NXPoints/8) || i==(6*NXPoints/8) || i==(7*NXPoints/8) || i==(8*NXPoints/8)){
                    for(unsigned jj=1;jj<npts;jj++){
                        cout << i << " " << xtmp[jj] << " " << ytmp[jj] << " " << ytmpError[jj] << " " << getipDataFromPolyModel(xtmp[jj],i,j) << " "  << fabs(xtmp[jj] - xtmp[jj-1]) << endl;
                    }
                }  
#endif      
            }
        }
	}	
}

void operaInstrumentProfile::FitMediantoIPDataVector(void) {
    ipPolyModel = PolynomialMatrix(NYPoints, NXPoints);
	chisqrMatrix = DMatrix(NYPoints, NXPoints);
	
	operaVector ytmp(nDataPoints);
	
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			unsigned nelem=0;
			for(unsigned index=0;index<nDataPoints;index++) {
				if(!isnan(getdataCubeValues(i,j,index))) {
					ytmp[nelem] = getdataCubeValues(i,j,index);
					nelem++;
				}
			}
			
			Polynomial pp(1);
			if(nelem) {
				ytmp.resize(nelem);
				pp.setCoefficient(0, Median(ytmp));
			} else {
				pp.setCoefficient(0, NAN);
			}
            double chisqr = ChiSqr(ytmp, pp.getCoefficient(0), nelem - 1);
            
            setipPolyModelCoefficients(pp, i,j);
            setchisqrMatrixValue(chisqr,i,j);
        }
	}
}

double operaInstrumentProfile::getchisqrMatrixValue(unsigned i, unsigned j) const {
#ifdef RANGE_CHECK
	if (j >= NYPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (i >= NXPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return chisqrMatrix[j][i];
}

void operaInstrumentProfile::setchisqrMatrixValue(double Value, unsigned i, unsigned j){
#ifdef RANGE_CHECK
	if (j >= NYPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (i >= NXPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	chisqrMatrix[j][i] = Value;
}

void operaInstrumentProfile::deleteDataCubes(void) {
	dataCube.clear();
	distd.clear();
	nDataPoints = 0;    
}

double operaInstrumentProfile::getIPphotoCenterX(double d) const {
    double sumflux = 0;
    double xc = 0;
    for (unsigned j=0; j<NYPoints; j++) {		
        for (unsigned i=0; i<NXPoints; i++) {				
            sumflux += getipDataFromPolyModel(d,i,j);
            xc += getipDataFromPolyModel(d,i,j)*getIPixXCoordinate(i);
        }
    }	    
    return xc/sumflux;
}

double operaInstrumentProfile::getIPphotoCenterY(double d) const {
    double sumflux = 0;
    double yc = 0;
    for (unsigned i=0; i<NXPoints; i++) {    
        for (unsigned j=0; j<NYPoints; j++) {
            sumflux += getipDataFromPolyModel(d,i,j);
            yc += getipDataFromPolyModel(d,i,j)*getIPixYCoordinate(j);
        }
    }	    
    return yc/sumflux;    
}
