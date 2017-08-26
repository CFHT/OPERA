/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralLines
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

#include <iostream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaSpectralLines.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaException.h"

#include "libraries/operaFit.h"
#include "libraries/operaMath.h"
#include "libraries/operaStats.h"

/*!
 * operaSpectralLines
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the set of Spectral Lines.
 * \file operaSpectralLines.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaSpectralLines
 * \brief This class is a container to manipulate a set of spectral lines
 * \return none
 */

/*
 * Constructor
 */

operaSpectralLines::operaSpectralLines(void) :
nLines(0),
comparisonSpectrum(NULL),
nFeatures(0),
referenceLineWidth(0.0),
dispersiontype(distance_disp)
{
	for (unsigned f=0; f<MAXNUMBEROFFEATURES; f++) {
		spectralFeatures[f] = NULL;
	}
}

operaSpectralLines::operaSpectralLines(operaSpectralElements *ComparisonSpectrum, double ReferenceLineWidth) :
nLines(0),
comparisonSpectrum(NULL),
nFeatures(0),
referenceLineWidth(0.0),
dispersiontype(distance_disp)
{
	if (ComparisonSpectrum->getnSpectralElements() == 0) {
		throw operaException("operaSpectralLines: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	comparisonSpectrum = new operaSpectralElements(ComparisonSpectrum->getnSpectralElements(), ComparisonSpectrum->getSpectrumType());
	*comparisonSpectrum = *ComparisonSpectrum; // sets element by element
    referenceLineWidth = ReferenceLineWidth;
	for (unsigned f=0; f<MAXNUMBEROFFEATURES; f++) {
		spectralFeatures[f] = NULL;
	}
    
#ifdef PRINT_DEBUG    
    for(unsigned indexElem=0; indexElem < comparisonSpectrum->getnSpectralElements(); indexElem++) {
		cout << comparisonSpectrum->getdistd(indexElem) <<
        "\t" << comparisonSpectrum->getrawsumflux(indexElem) <<
        "\t" << sqrt(comparisonSpectrum->getrawsumfluxVariance(indexElem)) << endl;
    }
#endif    
}

operaSpectralLines::operaSpectralLines(operaSpectralElements *ComparisonSpectrum, double ReferenceLineWidth, dispersionaxis_t Dispersiontype) :
nLines(0),
comparisonSpectrum(NULL),
nFeatures(0),
referenceLineWidth(0.0),
dispersiontype(distance_disp)
{
    
    dispersiontype = Dispersiontype;
    
	if (ComparisonSpectrum->getnSpectralElements() == 0) {
		throw operaException("operaSpectralLines: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	comparisonSpectrum = new operaSpectralElements(ComparisonSpectrum->getnSpectralElements(), ComparisonSpectrum->getSpectrumType());
	*comparisonSpectrum = *ComparisonSpectrum; // sets element by element
    referenceLineWidth = ReferenceLineWidth;
	for (unsigned f=0; f<MAXNUMBEROFFEATURES; f++) {
		spectralFeatures[f] = NULL;
	}
}

operaSpectralLines::operaSpectralLines(unsigned Features, unsigned LinesPerFeature) :
nLines(0),
comparisonSpectrum(NULL),
nFeatures(0),
referenceLineWidth(0.0),
dispersiontype(distance_disp)
{
	if (Features > MAXNUMBEROFFEATURES) {
		throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (Features == 0) {
		throw operaException("operaSpectralLines: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (LinesPerFeature == 0) {
		throw operaException("operaSpectralLines: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	nLines = Features * LinesPerFeature;
	nFeatures = Features;
	for (unsigned f=0; f<Features; f++) {
		spectralFeatures[f] = new operaSpectralFeature(LinesPerFeature);
	}
}

operaSpectralLines::operaSpectralLines(unsigned Features, unsigned LinesPerFeature, dispersionaxis_t Dispersiontype) :
nLines(0),
comparisonSpectrum(NULL),
nFeatures(0),
referenceLineWidth(0.0),
dispersiontype(distance_disp)
{
    dispersiontype = Dispersiontype;
    
	if (Features > MAXNUMBEROFFEATURES) {
		throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (Features == 0) {
		throw operaException("operaSpectralLines: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (LinesPerFeature == 0) {
		throw operaException("operaSpectralLines: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	nLines = Features * LinesPerFeature;
	nFeatures = Features;
	for (unsigned f=0; f<Features; f++) {
		spectralFeatures[f] = new operaSpectralFeature(LinesPerFeature);
	}
}

/*
 * Destructor
 */
operaSpectralLines::~operaSpectralLines(void) {
	for (unsigned f=0; f<nFeatures; f++) {
		if (spectralFeatures[f])
			delete spectralFeatures[f];
		spectralFeatures[f] = NULL;
	}
	if (comparisonSpectrum) {
		delete comparisonSpectrum;
		comparisonSpectrum = NULL;
	}
	nLines = 0;
	nFeatures = 0;
}

/*
 * Setter/Getters
 */

void operaSpectralLines::setDispersiontype(dispersionaxis_t Dispersiontype) {
    dispersiontype = Dispersiontype;
}

dispersionaxis_t operaSpectralLines::getDispersiontype(void) const {
    return dispersiontype;
}

void operaSpectralLines::setnLines(unsigned NLines) {
      nLines = NLines;  
}

unsigned operaSpectralLines::getnLines(void) const {
      return nLines;  
}

void operaSpectralLines::setNFeatures(unsigned NFeatures) {
	if (NFeatures > MAXNUMBEROFFEATURES) {
		throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    nFeatures = NFeatures;
}

unsigned operaSpectralLines::getNFeatures(void) const {
    return nFeatures; 
}


void operaSpectralLines::setReferenceLineWidth(double ReferenceLineWidth) {
    referenceLineWidth = ReferenceLineWidth;
}

double operaSpectralLines::getReferenceLineWidth(void) const {
    return referenceLineWidth;
}    

const operaSpectralFeature *operaSpectralLines::getSpectralFeature(unsigned indexFeature) const{
	if (indexFeature >= nFeatures) {
		throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return spectralFeatures[indexFeature];  
}

operaSpectralFeature *operaSpectralLines::getSpectralFeature(unsigned indexFeature){
	if (indexFeature >= nFeatures) {
		throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return spectralFeatures[indexFeature];  
}

void operaSpectralLines::detectSpectralFeatures(double DetectionThreshold, double LocalMaxFilterWidth, double MinPeakDepth) {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (!comparisonSpectrum->getHasXCorrelation()) {
		throw operaException("operaSpectralLines::detectSpectralFeatures: ", operaErrorHasNoXCorrelation, __FILE__, __FUNCTION__, __LINE__);	
	}
    unsigned NPeaks = 0;
    unsigned NLines = 0;    
    unsigned np = comparisonSpectrum->getnSpectralElements();
	if (np == 0) {
		throw operaException("operaSpectralLines::detectSpectralFeatures: no spectrum: ", operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);	
	}    
    
	// DT NOTE FIX ME!!!, just a workaround, should be np-1
	// horrible hack!!! Feb 21 2013
	unsigned endofspectralvaluesindex = 1;
	while (getMX(np-endofspectralvaluesindex) == 0) {
		endofspectralvaluesindex++;
		if (endofspectralvaluesindex == np) {
			throw operaException("operaSpectralLines::detectSpectralFeatures: no spectrum: ", operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);	
		}    
	}
    double xscale = fabs(getMX(np-endofspectralvaluesindex)-getMX(0))/(double)(np-endofspectralvaluesindex); // note: this could easily cause an error if either 1st or last point in array is screwed up
    
    unsigned slit = 4*(unsigned)round(referenceLineWidth/xscale);
    unsigned subslit = (unsigned)round(LocalMaxFilterWidth/xscale);

#ifdef PRINT_DEBUG
    cout << "operaSpectralLines: slit=" << slit << " subslit=" << subslit << " xscale=" << xscale << endl;
#endif
    
    if(subslit > slit)
        subslit=slit;
    
    double *xpeak = (double *) malloc (np * sizeof(double));
    double *ypeak = (double *) malloc (np * sizeof(double));
    double *peakAmplitude = (double *) malloc (np * sizeof(double));
    double *peakBackground = (double *) malloc (np * sizeof(double));     
    unsigned *peakIndex = (unsigned *) malloc (np * sizeof(unsigned));
	if (!peakIndex) {
		throw operaException("operaSpectralLines: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	
    for(unsigned i=slit;i<(np-slit);i++) {
                
        double x0 = getMX(i-slit/2);
        double y0 = getMY(i-slit/2);
        double xf = getMX(i+slit/2);
        double yf = getMY(i+slit/2);
        
        unsigned npts=0;
        
        double netYdisplacement=0;
        double totalYincrement=0;
        
        double itismaxprob = 1;
        
        for(unsigned j=0;j<slit;j++) {
            
            double ybkg0 = y0; 
            double ybkg1 = y0;
            double yincrement;
            
            if(j<slit/2) {
                ybkg0 += ((yf-y0)/(xf-x0))*(getMX(i-slit/2+j) - x0); 
                ybkg1 += ((yf-y0)/(xf-x0))*(getMX(i-slit/2+j+1) - x0);
                
                yincrement = (getMY(i-slit/2+j+1)-ybkg1) - (getMY(i-slit/2+j)-ybkg0);
                if(yincrement > 0) {
                    npts++;     //count points whenever it grows before center or decreases after center
                }
                netYdisplacement += yincrement; 
                totalYincrement += fabs(yincrement);
            } else if (j>slit/2) {
                ybkg0 += ((yf-y0)/(xf-x0))*(getMX(i-slit/2+j-1) - x0); 
                ybkg1 += ((yf-y0)/(xf-x0))*(getMX(i-slit/2+j) - x0);   
                
                yincrement = (getMY(i-slit/2+j)-ybkg1) - (getMY(i-slit/2+j-1)-ybkg0);
                if(yincrement < 0) {
                    npts++;     //count points whenever it grows before center or decreases after center
                } 
                netYdisplacement -= yincrement;
                totalYincrement += fabs(yincrement);                
            } 

            if(j != slit/2 && j >= (slit/2 - subslit/2) && j < (slit/2 + subslit/2)) {
                // make sure my[i] is a local maximum within the subslit
                if(getMY(i-slit/2+j) >= getMY(i)) {
                    itismaxprob = 0; 
                }
            }
            
            
        }	 
        
        double detectprob = 0;
        if(netYdisplacement > MinPeakDepth) {
            detectprob = netYdisplacement/totalYincrement; 
        }
        
        double xcorrelation = 0;
        if(comparisonSpectrum->getXCorrelation(i) > 0) {
            xcorrelation = comparisonSpectrum->getXCorrelation(i);  
        }
        double peakprob = itismaxprob*detectprob*xcorrelation;
        
#ifdef PRINT_DEBUG
        cout <<  getMX(i) << " " << getMY(i) << " " << itismaxprob << " " << detectprob << " " << xcorrelation << " " << peakprob << " " << netYdisplacement << " " << totalYincrement << " " << MinPeakDepth << " " << DetectionThreshold << endl;
#endif
        if(peakprob > DetectionThreshold) {
            xpeak[NPeaks] = getMX(i);
            ypeak[NPeaks] = getMY(i);
            
            double ybkgFromSlitBorders = y0 + ((yf-y0)/(xf-x0))*(getMX(i) - x0);
            double ybkgFromSubSlitBorders = ((getMY(i - subslit/2) + getMY(i + subslit/2 - 1))/2);
            
            double ybkg = 0;
            if(ybkgFromSlitBorders < ybkgFromSubSlitBorders) {
                ybkg = ybkgFromSlitBorders;
            } else if (ybkgFromSlitBorders >= ybkgFromSubSlitBorders && ybkgFromSubSlitBorders < getMY(i)) {
                ybkg = ybkgFromSubSlitBorders;
            }
            
            peakAmplitude[NPeaks] = getMY(i) - ybkg; 
            peakBackground[NPeaks] = ybkg;           
            peakIndex[NPeaks] = i;
            
#ifdef PRINT_DEBUG             
            cout << "peak #" << NPeaks << " " << xpeak[NPeaks] << " " << ypeak[NPeaks] << endl;            
#endif            
            NPeaks++;
        }
    }	
    
    double minPeakSeparation = (double)slit + (double)slit/4.0 + 2; // half slit for each line plus space for background
    
    unsigned NFeatures = 0;

    for(unsigned k=0;k<NPeaks;k++) {
        unsigned firstK = k;
        unsigned multiplicity = 1;

        while (k < (NPeaks - 1) && ((double)peakIndex[k+1] - (double)peakIndex[k]) <= minPeakSeparation) {
            multiplicity++;
            k++;
        }
		if (NFeatures > MAXNUMBEROFFEATURES) {
			throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		}
        
        unsigned i1 = peakIndex[firstK]-(slit/2+slit/4+1);
        unsigned i2 = peakIndex[k]+(slit/2+slit/4+1);
        
        if(i1 < 0) {i1 = 0;}
		// DT NOTE FIX ME, just a workaround, should be np-1
        //if(i2 >= np) {i2 = np-1;}
        if(i2 >= np-(endofspectralvaluesindex-1)) {i2 = np-endofspectralvaluesindex;}
        
        // Below only add feature if number of points is greater than the number of parameters to fit
        // Since we are fitting gaussians with background, the # of parameters is given by (4 x npeaks + 2)
        
        unsigned NPointsInFeature = unsigned(i2-i1);

        if(unsigned((NPointsInFeature-2) / 4*multiplicity) >= 1) {
            
            if (spectralFeatures[NFeatures]) {
                delete spectralFeatures[NFeatures];
            }
            spectralFeatures[NFeatures] = new operaSpectralFeature(multiplicity);
            for (unsigned kk=0; kk<multiplicity; kk++) {
                spectralFeatures[NFeatures]->getGaussianFit()->setAmplitude(peakAmplitude[firstK+kk],kk);
                spectralFeatures[NFeatures]->getGaussianFit()->setSigma(referenceLineWidth,kk);
                spectralFeatures[NFeatures]->getGaussianFit()->setCenter(xpeak[firstK+kk],kk);
            }
            
            spectralFeatures[NFeatures]->setOriginalIndex(i1,i2);
            
            spectralFeatures[NFeatures]->CreateDataVectors(NPointsInFeature);
            
            unsigned dataIndex = 0;
            for(unsigned i=i1;i<i2;i++) {
                spectralFeatures[NFeatures]->setDataPoint(getMX(i),getMY(i),sqrt(getVar(i)),dataIndex);
                dataIndex++;
            }
            
            /*
             * At this point the data, the number of peaks, and the inital guess for the
             * gaussian parameters are already loaded into the feature class. Next step is
             * to perform a least square fit to a multiple gaussian plus background model.
             */
            
            spectralFeatures[NFeatures]->fitBackground();
            
            spectralFeatures[NFeatures]->fitGaussianModel();
            
            NLines += multiplicity;
            
            NFeatures++;
        }
    }

    nFeatures = NFeatures;
    
    setnLines(NLines);
    
    free(xpeak);
    free(ypeak);
    free(peakIndex);
    free(peakAmplitude);	// DT two more to free...
    free(peakBackground);

}

void operaSpectralLines::detectAbsorptionSpectralFeatures(double DetectionThreshold, double LocalMaxFilterWidth, double MinPeakDepth) {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (!comparisonSpectrum->getHasXCorrelation()) {
		throw operaException("operaSpectralLines::detectSpectralFeatures: ", operaErrorHasNoXCorrelation, __FILE__, __FUNCTION__, __LINE__);
	}
    unsigned NPeaks = 0;
    unsigned NLines = 0;
    unsigned np = comparisonSpectrum->getnSpectralElements();
	if (np == 0) {
		throw operaException("operaSpectralLines::detectSpectralFeatures: no spectrum: ", operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
	}
    
	// DT NOTE FIX ME, just a workaround, should be np-1
	unsigned endofspectralvaluesindex = 1;
	while (getMX(np-endofspectralvaluesindex) == 0) {
		endofspectralvaluesindex++;
		if (endofspectralvaluesindex == np) {
			throw operaException("operaSpectralLines::detectSpectralFeatures: no spectrum: ", operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);	
		}    
	}
    double xscale = fabs(getMX(np-endofspectralvaluesindex)-getMX(0))/(double)(np-endofspectralvaluesindex); // note: this could easily cause an error if either 1st or last point in array is screwed up
    
    unsigned slit = 4*(unsigned)round(referenceLineWidth/xscale);
    unsigned subslit = (unsigned)round(LocalMaxFilterWidth/xscale);
    if(subslit > slit)
        subslit=slit;
    
    double *xpeak = (double *) malloc (np * sizeof(double));
    double *ypeak = (double *) malloc (np * sizeof(double));
    double *peakAmplitude = (double *) malloc (np * sizeof(double));
    double *peakBackground = (double *) malloc (np * sizeof(double));
    unsigned *peakIndex = (unsigned *) malloc (np * sizeof(unsigned));
	if (!peakIndex) {
		throw operaException("operaSpectralLines: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);
	}
	// hack DT Mar 15 2013 -- get out if the slit doesn't make sense
	if (slit > 16) {
		nFeatures = 0;
		setnLines(0);
		return;
	}
	
    for(unsigned i=slit;i<(np-slit);i++) {
        
        double x0 = getMX(i-slit/2);
        double y0 = getAbsorptionMY(i-slit/2);
        double xf = getMX(i+slit/2);
        double yf = getAbsorptionMY(i+slit/2);
        
        unsigned npts=0;
        
        double netYdisplacement=0;
        double totalYincrement=0;
        
        double itismaxprob = 1;
        
        for(unsigned j=0;j<slit;j++) {
            
            double ybkg0 = y0;
            double ybkg1 = y0;
            double yincrement;
            
            if(j<slit/2) {
                ybkg0 += ((yf-y0)/(xf-x0))*(getMX(i-slit/2+j) - x0);
                ybkg1 += ((yf-y0)/(xf-x0))*(getMX(i-slit/2+j+1) - x0);
                
                yincrement = (getAbsorptionMY(i-slit/2+j+1)-ybkg1) - (getAbsorptionMY(i-slit/2+j)-ybkg0);
                if(yincrement > 0) {
                    npts++;     //count points whenever it grows before center or decreases after center
                }
                netYdisplacement += yincrement;
                totalYincrement += fabs(yincrement);
            } else if (j>slit/2) {
                ybkg0 += ((yf-y0)/(xf-x0))*(getMX(i-slit/2+j-1) - x0);
                ybkg1 += ((yf-y0)/(xf-x0))*(getMX(i-slit/2+j) - x0);
                
                yincrement = (getAbsorptionMY(i-slit/2+j)-ybkg1) - (getAbsorptionMY(i-slit/2+j-1)-ybkg0);
                if(yincrement < 0) {
                    npts++;     //count points whenever it grows before center or decreases after center
                }
                netYdisplacement -= yincrement;
                totalYincrement += fabs(yincrement);
            }
            
            if(j != slit/2 && j >= (slit/2 - subslit/2) && j < (slit/2 + subslit/2)) {
                // make sure my[i] is a local maximum within the subslit
                if(getAbsorptionMY(i-slit/2+j) >= getAbsorptionMY(i)) {
                    itismaxprob = 0;
                }
            }
            
            
        }
        
        double detectprob = 0;
        if(netYdisplacement > MinPeakDepth) {
            detectprob = netYdisplacement/totalYincrement;
        }
        
        double xcorrelation = 0;
        if(comparisonSpectrum->getXCorrelation(i) > 0) {
            xcorrelation = 1 - comparisonSpectrum->getXCorrelation(i);
        }
        double peakprob = itismaxprob*detectprob*xcorrelation;
        
#ifdef PRINT_DEBUG
        cout <<  getMX(i) << " " << getAbsorptionMY(i) << " " << itismaxprob << " " << detectprob << " " << xcorrelation << " " << peakprob << " " << netYdisplacement << " " << totalYincrement << " " << MinPeakDepth << " " << DetectionThreshold << endl;
#endif
        if(peakprob > DetectionThreshold) {
            xpeak[NPeaks] = getMX(i);
            ypeak[NPeaks] = getAbsorptionMY(i);
            
            double ybkgFromSlitBorders = y0 + ((yf-y0)/(xf-x0))*(getMX(i) - x0);
            double ybkgFromSubSlitBorders = ((getAbsorptionMY(i - subslit/2) + getAbsorptionMY(i + subslit/2 - 1))/2);
            
            double ybkg = 0;
            if(ybkgFromSlitBorders < ybkgFromSubSlitBorders) {
                ybkg = ybkgFromSlitBorders;
            } else if (ybkgFromSlitBorders >= ybkgFromSubSlitBorders && ybkgFromSubSlitBorders < getAbsorptionMY(i)) {
                ybkg = ybkgFromSubSlitBorders;
            }
            
            peakAmplitude[NPeaks] = getAbsorptionMY(i) - ybkg;
            peakBackground[NPeaks] = ybkg;
            peakIndex[NPeaks] = i;
            
#ifdef PRINT_DEBUG
            cout << "peak #" << NPeaks << " " << xpeak[NPeaks] << " " << ypeak[NPeaks] << endl;
#endif
            NPeaks++;
        }
    }
    double minPeakSeparation = (double)slit + (double)slit/4.0 + 2; // half slit for each line plus space for background
        
    unsigned NFeatures = 0;
    
    for(unsigned k=0;k<NPeaks;k++) {
        unsigned firstK = k;
        unsigned multiplicity = 1;
        
        while (k < (NPeaks - 1) && ((double)peakIndex[k+1] - (double)peakIndex[k]) <= minPeakSeparation) {
            multiplicity++;
            k++;
        }
		if (NFeatures > MAXNUMBEROFFEATURES) {
			throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
		}
        
        unsigned i1 = peakIndex[firstK]-(slit/2+slit/4+1);
        unsigned i2 = peakIndex[k]+(slit/2+slit/4+1);
        
        if(i1 < 0) {i1 = 0;}
		// DT NOTE FIX ME, just a workaround, should be np-1
        //if(i2 >= np) {i2 = np-1;}
        if(i2 >= np-(endofspectralvaluesindex-1)) {i2 = np-endofspectralvaluesindex;}
        
        // Below only add feature if number of points is greater than the number of parameters to fit
        // Since we are fitting gaussians with background, the # of parameters is given by (4 x npeaks + 2)
        unsigned NPointsInFeature = unsigned(i2-i1);
        
        if(unsigned((NPointsInFeature-2) / 4*multiplicity) >= 1) {
            
            if (spectralFeatures[NFeatures]) {
                delete spectralFeatures[NFeatures];
            }
            spectralFeatures[NFeatures] = new operaSpectralFeature(multiplicity);
            for (unsigned kk=0; kk<multiplicity; kk++) {
                spectralFeatures[NFeatures]->getGaussianFit()->setAmplitude(peakAmplitude[firstK+kk],kk);
                spectralFeatures[NFeatures]->getGaussianFit()->setSigma(referenceLineWidth,kk);
                spectralFeatures[NFeatures]->getGaussianFit()->setCenter(xpeak[firstK+kk],kk);
            }
            
            spectralFeatures[NFeatures]->setOriginalIndex(i1,i2);
                    
            spectralFeatures[NFeatures]->CreateDataVectors(NPointsInFeature);
            
            unsigned dataIndex = 0;
            for(unsigned i=i1;i<i2;i++) {
                spectralFeatures[NFeatures]->setDataPoint(getMX(i),getMY(i),sqrt(getVar(i)),dataIndex);
                dataIndex++;
            }
            
            /*
             * At this point the data, the number of peaks, and the inital guess for the
             * gaussian parameters are already loaded into the feature class. Next step is
             * to perform a least square fit to a multiple gaussian plus background model.
             */
            
            spectralFeatures[NFeatures]->fitBackground();
            
            spectralFeatures[NFeatures]->fitGaussianModel();
            
            NLines += multiplicity;
            
            NFeatures++;
        }
    }
    
    nFeatures = NFeatures;
    
    setnLines(NLines);
    
    free(xpeak);
    free(ypeak);
    free(peakIndex);
    free(peakAmplitude);	// DT two more to free...
    free(peakBackground);
    
}


void operaSpectralLines::subtractFeatureModel(void) {

	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned k=0;k<nFeatures;k++) {
                
        unsigned i1 = spectralFeatures[k]->getOriginalInitialIndex();
        unsigned i2 = spectralFeatures[k]->getOriginalFinalIndex();        
        
        for(unsigned i=i1;i<i2;i++) {
            setMY(getMY(i)-spectralFeatures[k]->getGaussianFit()->EvaluateGaussian(getMX(i)), i);
        }     
    } 
#ifdef PRINT_DEBUG
    for(unsigned i=0;i<comparisonSpectrum->getnSpectralElements();i++) {
        cout << getMX(i) << " " << getMY(i) << endl;
    }      
#endif    
}

void operaSpectralLines::printLines(ostream *pout) const {
	if (pout) {    
        unsigned totalNlines = 0;
        for(unsigned k=0;k<nFeatures;k++) {
            double *center = spectralFeatures[k]->getGaussianFit()->getCenterVector();
            double *sigma = spectralFeatures[k]->getGaussianFit()->getSigmaVector();
            double *amplitude = spectralFeatures[k]->getGaussianFit()->getAmplitudeVector();        
            double *centerError = spectralFeatures[k]->getGaussianFit()->getCenterErrorVector();
            double *sigmaError = spectralFeatures[k]->getGaussianFit()->getSigmaErrorVector();
            double *amplitudeError = spectralFeatures[k]->getGaussianFit()->getAmplitudeErrorVector(); 
            
            for(unsigned line=0; line<spectralFeatures[k]->getnLines(); line++) {
                *pout << totalNlines <<
                " " << k <<             
                " " << spectralFeatures[k]->getnLines() <<
                " " << spectralFeatures[k]->getGaussianFit()->getBaselineIntercept() << 
                " " << spectralFeatures[k]->getGaussianFit()->getBaselineInterceptError() <<            
                " " << spectralFeatures[k]->getGaussianFit()->getBaselineSlope() << 
                " " << spectralFeatures[k]->getGaussianFit()->getBaselineSlopeError() << 
                " " << center[line] << 
                " " << centerError[line] <<             
                " " << sigma[line] <<
                " " << sigmaError[line] <<            
                " " << amplitude[line] << 
                " " << amplitudeError[line] << 
                " " << spectralFeatures[k]->getGaussianFit()->getGaussianChisqr() << endl;
                totalNlines++;
            }        
        }   
    }
}

void operaSpectralLines::printReferenceSpectrum(ostream *pout) const {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (pout) {   
        unsigned np = comparisonSpectrum->getnSpectralElements();
        for(unsigned i=0;i<np;i++) {        
                *pout << getMX(i) <<
                " " << getMY(i) << 
                " " << sqrt(getVar(i)) << endl;
        }   
    }
}

unsigned operaSpectralLines::selectLines(double MaxContamination, unsigned nSig, double amplitudeCutOff, double *LinePositionVector,double *LineSigmaVector,double *LineAmplitudeVector) {
    unsigned nSelectedLines=0;  
    
    double *center = (double *) malloc (nLines * sizeof(double));
    float *sigma = (float *) malloc (nLines * sizeof(float));
    double *amplitude = (double *) malloc (nLines * sizeof(double)); 
    
    double *featureArea = (double *) malloc (nFeatures * sizeof(double));
	if (!featureArea) {
		throw operaException("operaSpectralLines: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    
    unsigned totalNlines = 0;    
    for(unsigned k=0;k<nFeatures;k++) {
        
        double *center_tmp = spectralFeatures[k]->getGaussianFit()->getCenterVector();        
        double *sigma_tmp = spectralFeatures[k]->getGaussianFit()->getSigmaVector();
        double *amplitude_tmp = spectralFeatures[k]->getGaussianFit()->getAmplitudeVector(); 
        
        featureArea[k] = 0;
        if (spectralFeatures[k]->getGaussianFit()->getNumberOfPeaks() < spectralFeatures[k]->getnLines()) {
			throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		}
        for(unsigned line=0; line<spectralFeatures[k]->getnLines(); line++) {
			if (totalNlines >= nLines) {
				throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
			}
            center[totalNlines] = center_tmp[line];            
            sigma[totalNlines] = (float)sigma_tmp[line];
            amplitude[totalNlines] = amplitude_tmp[line];
            
            featureArea[k] += sigma_tmp[line]*amplitude_tmp[line];
            
            totalNlines++;
        }        
    }             
    
    float MedianSigma = operaArrayMedian(totalNlines,sigma);  
    float MedSigSigma = operaArrayMedianSigma(totalNlines,sigma,MedianSigma);       
    
	unsigned j = 0;
	//
	// danger - how do we know that nSelectedLines length is >= sizeof(the vectors)
	//
    for(unsigned k=0;k<nFeatures;k++) {  
		if (nLines < spectralFeatures[k]->getnLines()) {
			throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		}
        for(unsigned line=0; line<spectralFeatures[k]->getnLines(); line++) {
			if (j >= nLines) {
				throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
			}
            if(sigma[j]*amplitude[j]/featureArea[k] >= (1.0 - MaxContamination) &&
               sigma[j] > MedianSigma - float(nSig)*MedSigSigma &&
               sigma[j] < MedianSigma + float(nSig)*MedSigSigma &&
               amplitude[j] > amplitudeCutOff) {
                
                LinePositionVector[nSelectedLines] = center[j];
                LineSigmaVector[nSelectedLines] = (double)sigma[j];                
                LineAmplitudeVector[nSelectedLines] = amplitude[j];                
                            
                nSelectedLines++;
            }
            j++;
        }        
    }  
    
    free(center);
    free(sigma);
    free(amplitude);
    free(featureArea);
    
    return nSelectedLines;
}

// Note that these are length protected...
inline double operaSpectralLines::getMX(unsigned k) const {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#ifdef RANGE_CHECK
	if (k >= comparisonSpectrum->getnSpectralElements()) {
		throw operaException("operaSpectralLines::detectSpectralFeatures:getMX "+itos(k)+">="+itos(comparisonSpectrum->getnSpectralElements())+' ', operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	switch (dispersiontype) {
		case wavelength_disp: {
			return comparisonSpectrum->getwavelength(k);
		}
			break;
		case distance_disp: {
			return comparisonSpectrum->getdistd(k);        
		}
			break;
		case y_distance_disp: {
			return comparisonSpectrum->getphotoCenterY(k);
		}
			break;            
		default:
			break;
	}
	return NAN;
}

inline double operaSpectralLines::getMY(unsigned k) const {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return comparisonSpectrum->getFlux(k);
}

inline double operaSpectralLines::getAbsorptionMY(unsigned k) const {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return (1.0 - comparisonSpectrum->getFlux(k));
}

inline void operaSpectralLines::setMY(double value, unsigned k) {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	comparisonSpectrum->setFlux(value, k);
}

inline double operaSpectralLines::getVar(unsigned k) const {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return comparisonSpectrum->getFluxVariance(k);
}

inline void operaSpectralLines::setVar(double value, unsigned k) {
	if (comparisonSpectrum == NULL) {
		throw operaException("operaSpectralLines ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	comparisonSpectrum->setFluxVariance(value, k);
}

