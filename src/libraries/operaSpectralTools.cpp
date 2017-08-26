/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
Library name: operaSpectralTools - common C++ library functions
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

/*!
 * operaSpectralTools
 * \author Eder Martioli
 * \brief operaSpectralTools - common C++ library functions.
 * \file operaSpectralTools.cpp
 * \ingroup libraries
 */


#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaSpectralTools.h"

#include "libraries/operaSpectralLines.h"       // for operaSpectralLines
#include "libraries/operaSpectralFeature.h"    // for operaSpectralFeature
#include "libraries/operaLibCommon.h"           // for SPEED_OF_LIGHT_M
#include "libraries/operaStats.h"               // for operaCrossCorrelation
#include "libraries/operaFFT.h"                 // for operaXCorrelation
#include "libraries/operaException.h"
#include "libraries/operaFit.h"

using namespace std;

class Always {
public:
	bool operator()(double center, double sigma) { return true; }
};

class CenterInRange {
	operaWavelengthRange range;
public:
	CenterInRange(double wl0, double wlf) : range(wl0, wlf) { }
	bool operator()(double center, double sigma) { return range.contains(center); }
};

class WidthInRange {
	operaWavelengthRange range;
public:
	WidthInRange(double min, double max) : range(min, max) { }
	bool operator()(double center, double width) { return range.contains(width); }
};

template <class T>
operaSpectralLineList getSpectralLinesMatchingCondition(const operaSpectralLines& spectralLines, T condition) {
	operaSpectralLineList linelist;
	for (unsigned featurenumber = 0; featurenumber < spectralLines.getNFeatures(); featurenumber++) {
		const operaSpectralFeature *spectralFeature = spectralLines.getSpectralFeature(featurenumber);
		const Gaussian* gaussianFit = spectralFeature->getGaussianFit();
		const double *center = gaussianFit->getCenterVector();
		const double *centerError = gaussianFit->getCenterErrorVector();
		const double *sigma = gaussianFit->getSigmaVector();
		const double *amplitude = gaussianFit->getAmplitudeVector();
		for (unsigned line=0; line<spectralFeature->getnLines(); line++) {
			if (condition(center[line], sigma[line])) {
				linelist.center.insert(center[line]);
				linelist.centerError.insert(centerError[line]);
				linelist.sigma.insert(sigma[line]);
				linelist.amplitude.insert(amplitude[line]);
			}
		}
	}
	if(linelist.size() > 0) {
		linelist.medianWidth = Median(linelist.sigma);
		linelist.medianWidthError = MedianStdDev(linelist.sigma, linelist.medianWidth);
	} else {
		linelist.medianWidth = NAN;
		linelist.medianWidthError = NAN;
	}
	return linelist;
}

operaSpectralLineList getAllSpectralLines(const operaSpectralLines& spectralLines) {
	return getSpectralLinesMatchingCondition(spectralLines, Always());
}

operaSpectralLineList getSpectralLinesInWavelengthRange(const operaSpectralLines& spectralLines, operaWavelengthRange wlrange) {
	return getSpectralLinesMatchingCondition(spectralLines, CenterInRange(wlrange.getwl0(), wlrange.getwlf()));
}

operaSpectralLineList getSpectralLinesNearMedianWidth(const operaSpectralLines& spectralLines, double nsigclip) {
	operaVector linesigmas;
	for (unsigned feature = 0; feature<spectralLines.getNFeatures(); feature++) {
		const operaSpectralFeature *currentFeature = spectralLines.getSpectralFeature(feature);
		linesigmas.append(operaVector(currentFeature->getGaussianFit()->getSigmaVector(), currentFeature->getnLines()));
	}
	double medianWidth = NAN;
	double medianWidthError = NAN;
	if(linesigmas.size() > 0) {
		medianWidth = MedianQuick(linesigmas);
		medianWidthError = MedianStdDev(linesigmas, medianWidth);
	}
	double minwidth = medianWidth - nsigclip * medianWidthError;
	double maxwidth = medianWidth + nsigclip * medianWidthError;
	return getSpectralLinesMatchingCondition(spectralLines, WidthInRange(minwidth, maxwidth));
}

operaSpectralLines DetectSpectralLines(operaSpectralElements* spectralElements, double linewidth, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, dispersionaxis_t dispersiontype) {
	operaSpectralLines compLines(spectralElements, linewidth, dispersiontype);
	double CompLocalMaxFilterWidth = LocalMaxFilterWidth*linewidth;
	double meanVariance = Mean(FilterNans(spectralElements->getFluxVector().getvariance()));
	double CompMinPeakDepth = MinPeakDepth*sqrt(meanVariance);
	compLines.detectSpectralFeatures(DetectionThreshold, CompLocalMaxFilterWidth, CompMinPeakDepth);
	return compLines;
}


class GaussianFunc {
	double a, b, c;
public:
	GaussianFunc(double height, double center, double sigma) : a(height), b(center), c(sigma) { }
	double operator()(double x) { return a*exp(-(x-b)*(x-b)/(2*c*c)); }
};

class NormalizedGaussianFunc : public GaussianFunc {
public:
	NormalizedGaussianFunc(double center, double sigma) : GaussianFunc(1.0/(sigma * sqrt(2.0*M_PI)), center, sigma) { }
};

double operaCrossCorrelation(operaVector a, operaVector b) {
	a -= Mean(a);
	b -= Mean(b);
    return InnerProduct(a, b) / sqrt(InnerProduct(a, a) * InnerProduct(b, b));
}

operaVector calculateXCorrWithGaussian(const operaVector& wavelength, const operaVector& flux, double sigma) {
    unsigned np = wavelength.size();
    operaVector outputXcorr;
    for(unsigned i=0; i<np; i++) {
		double wlstep;
        if(i==0) {
            wlstep = fabs(wavelength[i+1] - wavelength[i]);
        } else if (i==np-1) {
            wlstep = fabs(wavelength[i] - wavelength[i-1]);
        } else {
            wlstep = fabs(wavelength[i+1] - wavelength[i-1])/2.0;
        }
        
        unsigned window = (unsigned)ceil(2*sigma/wlstep);
        if (window > np/2) window = np/2;
        
        operaVector windowFunc;
        operaVector mainFunc;
        
        NormalizedGaussianFunc g(i, window/2.0);
        unsigned minj = i > window ? i-window : 0;
        unsigned maxj = i+window < np ? i+window : np;
        for(unsigned j=minj; j<maxj; j++) {
			windowFunc.insert(g(j));
            mainFunc.insert(flux[j]);
        }
        outputXcorr.insert(operaCrossCorrelation(windowFunc, mainFunc));
    }
    return outputXcorr;
}

operaVector convolveSpectrumWithGaussian(const operaVector& wavelength, const operaVector& flux, double sigma) {
    unsigned np = wavelength.size();
	operaVector convolvedSpectrum(np);
    for(unsigned i=0;i<np;i++) {
		double wlstep;
        if(i==0) {
            wlstep = fabs(wavelength[i+1] - wavelength[i]);
        } else if (i==np-1) {
            wlstep = fabs(wavelength[i] - wavelength[i-1]);
        } else {
            wlstep = fabs(wavelength[i+1] - wavelength[i-1])/2.0;
        }
        
        unsigned window = (unsigned)ceil(2*sigma/wlstep);
        if (window > np/2) window = np/2;
        
        double weighSum = 0;
        NormalizedGaussianFunc g(wavelength[i], sigma);
        unsigned minj = i > window ? i-window : 0;
        unsigned maxj = i+window < np ? i+window : np;
        for(unsigned j=minj; j<maxj; j++) {
            double temp = g(wavelength[j]);
            convolvedSpectrum[i] += flux[j] * temp;
            weighSum += temp;
        }
        if (weighSum) {
            convolvedSpectrum[i] /= weighSum;
        }
    }
    return convolvedSpectrum;
}

operaVector convolveSpectrumWithGaussianByResolution(const operaVector& wavelength, const operaVector& flux, double spectralResolution) {
	unsigned np = wavelength.size();
	operaVector convolvedSpectrum(np);
    for(unsigned i=0;i<np;i++) {
		double wlstep;
        if(i==0) {
            wlstep = fabs(wavelength[i+1] - wavelength[i]);
        } else if (i==np-1) {
            wlstep = fabs(wavelength[i] - wavelength[i-1]);
        } else {
            wlstep = fabs(wavelength[i+1] - wavelength[i-1])/2.0;
        }
        
        double sigma = wavelength[i] / spectralResolution;
        unsigned window = (unsigned)ceil(2*sigma/wlstep);
        if (window > np/2) window = np/2;
        
        double weighSum = 0;
        NormalizedGaussianFunc g(wavelength[i], sigma);
        unsigned minj = i > window ? i-window : 0;
        unsigned maxj = i+window < np ? i+window : np;
        for(unsigned j=minj; j<maxj; j++) {
            if(j >= 0 && j < np) {
                double temp = g(wavelength[j]);
                convolvedSpectrum[i] += flux[j] * temp;
                weighSum += temp;
            }
        }
        if (weighSum) {
            convolvedSpectrum[i] /= weighSum;
        }
    }
    return convolvedSpectrum;
}

/*
 * normalizeSpectrum(unsigned nLines, double *lineflux)
 * \brief This function normalize an input vector of fluxes by the maximum value in the array
 * \param unsigned nLines
 * \param double *lineflux
 * \return void
 */
void normalizeSpectrum(unsigned nLines, double *lineflux) {
	double maxflux = 0;
    
    for (unsigned i=0; i<nLines; i++) {
		if(lineflux[i] > maxflux) {
            maxflux = lineflux[i];
        }
	}
	for (unsigned i=0; i<nLines; i++) {
		lineflux[i] /= maxflux/100;
	}
}

/*
 * normalizeSpectrum(unsigned nLines, double *lineflux)
 * \brief This function normalize an input vector of fluxes and respective variances by the maximum value in the array
 * \param unsigned nLines
 * \param double *lineflux
 * \param double *linevariance
 * \return void
 */
void normalizeSpectrum(unsigned nLines, double *lineflux, double *linevariance) {
	double maxflux = 0;
    
    for (unsigned i=0; i<nLines; i++) {
		if(lineflux[i] > maxflux) {
            maxflux = lineflux[i];
        }
	}
	for (unsigned i=0; i<nLines; i++) {
		lineflux[i] /= maxflux/100;
        linevariance[i] /= (maxflux/100)*(maxflux/100);
	}
}

/*
 * double convertVacuumToAirWavelength(double vac_wl)
 * \brief This function converts from vacuum to air wavelength using the IAU standard for conversion
 * \brief  from air to vacuum wavelengths as given in Morton (1991, ApJS, 77, 119)
 * \param double vac_wl (in Angstrom)
 * \return double air_wl (in Angstrom)
 */
double convertVacuumToAirWavelength(double vac_wl) {
    return vac_wl/(1.0 + 0.0002735182 + (131.4182/(vac_wl*vac_wl)) + (276249000.0/(vac_wl*vac_wl*vac_wl*vac_wl)));
}


/*
 * double calculateBlackBodyVFlux(double Temperature)
 * \brief This function calculates the flux (in ph/(m^2 s m)) for a Black Body
 * \param double Temperature (in Kelvin)
 * \return double flux (in ph/(m^2 s m))
 */
double calculateBlackBodyVFlux(double Temperature) {
    
    double lamb[24], Vresp[24];
    double dlambda = 1e-8;
    
    for(unsigned i=0;i<24;i++) {
        lamb[i] = 470e-9 + (double)i*dlambda;
    }
    // Below is the spectral transmission for the V-band Johnson filter
    Vresp[0] = 0.000;
    Vresp[1] = 0.030;
    Vresp[2] = 0.163;
    Vresp[3] = 0.458;
    Vresp[4] = 0.780;
    Vresp[5] = 0.967;
    Vresp[6] = 1.000;
    Vresp[7] = 0.973;
    Vresp[8] = 0.898;
    Vresp[9] = 0.792;
    Vresp[10] = 0.684;
    Vresp[11] = 0.574;
    Vresp[12] = 0.461;
    Vresp[13] = 0.359;
    Vresp[14] = 0.270;
    Vresp[15] = 0.197;
    Vresp[16] = 0.135;
    Vresp[17] = 0.081;
    Vresp[18] = 0.045;
    Vresp[19] = 0.025;
    Vresp[20] = 0.017;
    Vresp[21] = 0.013;
    Vresp[22] = 0.009;
    Vresp[23] = 0.000;
    
    double fluxV = 0;
    
    for(unsigned i=0;i<24;i++) {
        fluxV += PlanckFunction(Temperature,lamb[i])*Vresp[i]*dlambda;
    }
    return fluxV;
}

/*
 * double planck(double T, double wl)
 * \brief This function returns the spectral radiance of a black body in photons/(s m^2 dlambda).
 * \param Temperature is a double input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a double input that represents the wavelength at which the black body is observed in nanometers
 * \return double value for the spectral radiance of the black body in ph/(s m^2 dlambda)
 */
double PlanckFunction(double T, double wl) {
    // http://spiff.rit.edu/classes/phys317/lectures/planck.html
    double K_B = 1.380658e-23;
    double H_PLANCK = 6.6260755e-34;
    
    // in units of W/(m^2 dlambda)
    //  flux = (2*M_PI*H_PLANCK*SPEED_OF_LIGHT_M*SPEED_OF_LIGHT_M)/(pow(wl,5)*(exp(H_PLANCK*SPEED_OF_LIGHT_M/(wl*K_B*T)) - 1.));
    
    // in units of photons/(s . m^2  dlambda)
    double flux = (2*M_PI*SPEED_OF_LIGHT_M)/(pow(wl,4)*(exp(H_PLANCK*SPEED_OF_LIGHT_M/(wl*K_B*T)) - 1.));
    
    return flux;
}


/*
 * double IntegrateSpectralElementOfBlackBody(double wl0, double wlf, double T)
 * \brief This function returns the integrated spectral flux of a black body in photons/(s m^2) for a given
 * \brief wavelength range using the simpson method.
 * \param Temperature is a double input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a double input that represents the wavelength at which the black body is observed in nanometers
 * \return double value for the spectral radiance of the black body in ph/(s m^2 dlambda)
 */
double IntegrateSpectralElementOfBlackBody(double wl0, double wlf, double T) {
    double f;
    unsigned N = 1000;
    double h = fabs(wlf - wl0)/((double)N - 1.);
    double xi = wl0;
    
    double sum = 0;
    for(unsigned i=0;i<=N;i++)
    {
        if(i==0) {
            f = PlanckFunction(T, wl0);
        } else if(i==N) {
            f = PlanckFunction(T, wlf);
        } else {
            xi += h;
            double aux = fmodf((double)(i),2.0);
            if(aux != 0) {
                f = 4.*PlanckFunction(T,xi);
            } else {
                f = 2.*PlanckFunction(T,xi);
            }
        }
        sum += (h/3.)*f;
    }
    
    return sum;
}


double getFactorToMatchFluxesBetweenElements(operaSpectralElements *refElements,operaSpectralElements *elementsToMatch,double delta_wl) {
    bool debug = false;
    double ref_wl0 = refElements->getwavelength(0);
    double ref_wlf = refElements->getwavelength(refElements->getnSpectralElements()-1);
    double elem2match_wl0 = elementsToMatch->getwavelength(0);
    double elem2match_wlf = elementsToMatch->getwavelength(elementsToMatch->getnSpectralElements()-1);
    
    if(ref_wl0 >= ref_wlf || elem2match_wl0 >= elem2match_wlf) {
        throw operaException("getOverlappingWLRange: initial wl must not be greater than final wl. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    // find intersecting range:
    bool elementsIntersect = false;
    double intersect_wl0 = 0;
    double intersect_wlf = 0;
    
    if(elem2match_wl0 >= ref_wl0 && elem2match_wl0 <= ref_wlf) {
        intersect_wl0 = elem2match_wl0;
        intersect_wlf = ref_wlf;
        elementsIntersect = true;
    }
    
    if (elem2match_wlf >= ref_wl0 && elem2match_wlf <= ref_wlf) {
        intersect_wlf = elem2match_wlf;
        if(intersect_wl0 == 0) {
            intersect_wl0 = ref_wl0;
        }
        elementsIntersect = true;
    }
    
    if(debug) {
        cout << "refElements     = " << ref_wl0 << " -> " << ref_wlf << endl;
        cout << "elementsToMatch = " << elem2match_wl0 << " -> " << elem2match_wlf << endl;
        if (elementsIntersect == true) {
            cout << "Intersection    = " << intersect_wl0 << " -> " << intersect_wlf << endl;
        } else {
            cout << " Elements Do Not Intersect" << endl;
        }
    }
    
    
    // first collect data:
    unsigned refNElements = refElements->getnSpectralElements();
    float *refFlux = new float[refNElements];
    unsigned nref = 0;
    
    unsigned nElements = elementsToMatch->getnSpectralElements();
    float *elemFlux = new float[nElements];
    unsigned nelem = 0;
    double factorToMatch = 1.0;
    
    if (elementsIntersect == true) {
        
        for(unsigned elemIndex=0; elemIndex<refNElements;elemIndex++) {
            if(refElements->getwavelength(elemIndex) >= intersect_wl0 &&
               refElements->getwavelength(elemIndex) <= intersect_wlf) {
                refFlux[nref] = refElements->getFlux(elemIndex);
                nref++;
            }
        }
        
        double medianRefFlux = (double)operaArrayMedian(nref,refFlux);
        
        for(unsigned elemIndex=0; elemIndex<nElements;elemIndex++) {
            if(elementsToMatch->getwavelength(elemIndex) >= intersect_wl0 &&
               elementsToMatch->getwavelength(elemIndex) <= intersect_wlf) {
                elemFlux[nelem] = elementsToMatch->getFlux(elemIndex);
                nelem++;
            }
        }
        double medianElemFlux = (double)operaArrayMedian(nelem,elemFlux);
        
        factorToMatch = medianRefFlux/medianElemFlux;
        if(debug) {
            cout << medianRefFlux << " " << medianElemFlux << " " << factorToMatch << endl;
        }
    } else {
        // find out which element set comes first
        // to find the wavelength between the two orders

        double refinflimit_wl, refsuplimit_wl;
        double eleminflimit_wl, elemsuplimit_wl;
        
        if(ref_wlf < elem2match_wl0) { // ref comes first
            refinflimit_wl = ref_wlf - delta_wl;
            refsuplimit_wl = ref_wlf;
            eleminflimit_wl = elem2match_wl0;
            elemsuplimit_wl = elem2match_wl0 + delta_wl;
            
        } else {            
            refinflimit_wl = ref_wl0;
            refsuplimit_wl = ref_wl0 + delta_wl;
            eleminflimit_wl = elem2match_wlf - delta_wl;
            elemsuplimit_wl = elem2match_wlf;
        }
        if(debug) {
            cout << "refinflimit_wl=" << refinflimit_wl << " refsuplimit_wl=" << refsuplimit_wl << endl;
            cout << "eleminflimit_wl=" << eleminflimit_wl << " elemsuplimit_wl=" << elemsuplimit_wl << endl;
        }
        
        for(unsigned elemIndex=0; elemIndex<refNElements;elemIndex++) {
            if(refElements->getwavelength(elemIndex) >= refinflimit_wl &&
               refElements->getwavelength(elemIndex) <= refsuplimit_wl) {
                refFlux[nref] = refElements->getFlux(elemIndex);
                nref++;
            }
        }
        
        double medianRefFlux = (double)operaArrayMedian(nref,refFlux);
        
        for(unsigned elemIndex=0; elemIndex<nElements;elemIndex++) {
            if(elementsToMatch->getwavelength(elemIndex) >= eleminflimit_wl &&
               elementsToMatch->getwavelength(elemIndex) <= elemsuplimit_wl) {
                elemFlux[nelem] = elementsToMatch->getFlux(elemIndex);
                nelem++;
            }
        }
        double medianElemFlux = (double)operaArrayMedian(nelem,elemFlux);
        
        factorToMatch = medianRefFlux/medianElemFlux;
        if(debug) {
            cout << medianRefFlux << " " << medianElemFlux << " " << factorToMatch << endl;
        }
    }
    delete[] refFlux;
    delete[] elemFlux;
    
    return factorToMatch;
}

operaVector interpolateFlux(const operaVector& fluxA, const operaVector& fluxB, double wlA, double wlB, double wlout) {
	operaVector slopes = (fluxB - fluxA) / (wlB - wlA);
	return fluxB + (wlout-wlB)*slopes; //intercept = fluxin-(slope*wlin), fluxout = intercept+(slope*wlout)
}

double interpolateFlux(double fluxA, double fluxB, double wlA, double wlB, double wlout) {
	double slope = (fluxB - fluxA) / (wlB - wlA);
	return fluxB + (wlout-wlB)*slope;
}

operaSpectrum calculateUniformSample(const operaSpectrum& input, unsigned npout) {
	if(npout < 2) {
        throw operaException("operaSpectralTools: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
    }
    operaSpectrum output(input.fluxcount(), 0);
    const double wl0 = input.firstwl(), wlf = input.lastwl();
    const double wlstep = fabs(wlf - wl0)/(npout-1);
    unsigned k=0;
    for(unsigned i=0; i<npout-1; i++) {
        double wlout = wl0 + i*wlstep;
        while(input.getwavelength(k) < wlout) k++;
        if(input.getwavelength(k) == wlout) output.insert(input.getwavelength(k), input.getfluxes(k));
        else output.insert(wlout, interpolateFlux(input.getfluxes(k-1), input.getfluxes(k), input.getwavelength(k-1), input.getwavelength(k), wlout));
    }
    output.insert(wlf, input.getfluxes(input.size()-1)); //insert last point manually to avoid missing it due to rounding errors
    return output;
}

float getFluxAtWavelength(unsigned np,float *wl,float *flux,float wavelengthForNormalization) {
    float wl0 = wl[0];
    float wlf = wl[np-1];
    
    if(wavelengthForNormalization < wl0 || wavelengthForNormalization > wlf) {
        cout << wl0 << " " << wavelengthForNormalization << " " << wlf << endl;

        throw operaException("operaSpectralTools: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    float fluxAtwavelengthForNormalization = 0;
    
    for (unsigned k=0; k<np;k++) {
        if(wl[k] >= wavelengthForNormalization && k>0) {
            float slope = (flux[k] - flux[k - 1])/ (wl[k] - wl[k - 1]);
            float intercept = flux[k] - slope*wl[k];
            fluxAtwavelengthForNormalization = intercept + slope*wavelengthForNormalization;
            break;
        }
    }
    return  fluxAtwavelengthForNormalization;
}

double getFluxAtWavelength(const operaVector& wavelength, const operaVector& flux, double target) {
    if(target < wavelength.first() || target > wavelength.last()) {
        throw operaException("operaSpectralTools: ", operaErrorPickOutofRange, __FILE__, __FUNCTION__, __LINE__);
    }
    for (unsigned i=0; i<wavelength.size(); i++) {
        if(wavelength[i] >= target) {
			if(wavelength[i] == target) return wavelength[i];
            return interpolateFlux(flux[i-1], flux[i], wavelength[i-1], wavelength[i], target);
        }
    }
    throw operaException("operaErrorCodeNOTIMPLEMENTED: ", operaErrorCodeBracketingError, __FILE__, __FUNCTION__, __LINE__); // This should never happen.
}

operaVector getFluxesAtWavelength(const operaSpectrum& spectrum, double target) {
    if(target < spectrum.firstwl() || target > spectrum.lastwl()) {
        throw operaException("operaSpectralTools: ", operaErrorPickOutofRange, __FILE__, __FUNCTION__, __LINE__);
    }
    for (unsigned i=0; i<spectrum.size(); i++) {
        if(spectrum.getwavelength(i) >= target) {
			if(spectrum.getwavelength(i) == target) return spectrum.getfluxes(i);
            return interpolateFlux(spectrum.getfluxes(i-1), spectrum.getfluxes(i), spectrum.getwavelength(i-1), spectrum.getwavelength(i), target);
        }
    }
    throw operaException("operaErrorCodeNOTIMPLEMENTED: ", operaErrorCodeBracketingError, __FILE__, __FUNCTION__, __LINE__); // This should never happen.
}

/*
 * Read mask
 */
operaWavelengthRanges readContinuumWavelengthMask(string wavelength_mask) {
	operaWavelengthRanges wlranges;
	ifstream astream(wavelength_mask.c_str());
	if (astream.is_open()) {
		string dataline;
		while (getline(astream, dataline)) {
			if (!dataline.empty() && dataline[0] != '#') { // skip comments
				double tmpwl0, tmpwlf;
				istringstream ss(dataline);
				if(ss >> tmpwl0 >> tmpwlf) wlranges.addWavelengthRange(tmpwl0, tmpwlf);
			}
		}
		astream.close();
	}
	return wlranges;
}


operaWavelengthRanges getWavelengthMaskAroundLines(const operaSpectrum sourceLines, double spectralResolution, double nsig) {
    
    operaWavelengthRanges wlranges;
    
    double wl0_ref  = sourceLines.firstwl()*(1.0 - nsig/spectralResolution);
    double wlf_ref  = sourceLines.firstwl()*(1.0 + nsig/spectralResolution);
    
    for (unsigned l=0; l<sourceLines.size(); l++) {
        double lineWidth = sourceLines.getwavelength(l)/spectralResolution;
        
        double tmpwl0 = sourceLines.getwavelength(l) - nsig*lineWidth;
        double tmpwlf = sourceLines.getwavelength(l) + nsig*lineWidth;
        
        if (tmpwl0 < wlf_ref && tmpwlf > wlf_ref) {
            wlf_ref = tmpwlf;
        } else if (tmpwl0 > wlf_ref) {
            wlranges.addWavelengthRange(wl0_ref, wlf_ref);
            wl0_ref = tmpwl0;
            wlf_ref = tmpwlf;
        }
    
        if(l==sourceLines.size()-1) {
            wlranges.addWavelengthRange(wl0_ref, wlf_ref);
        }
    }

    return wlranges;
}


unsigned readContinuumWavelengthMask(string wavelength_mask, double *wl0, double *wlf) {
	ifstream astream;
	string dataline;
    
	double tmpwl0 = -1.0;
	double tmpwlf = -1.0;
	unsigned np = 0;
	
	astream.open(wavelength_mask.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl0, &tmpwlf);
                    wl0[np] = tmpwl0;
                    wlf[np] = tmpwlf;
                    np++;
                }	// skip comments
            }
		} // while (astream.good())
		astream.close();
	}	// if (astream.open())
	return np;
}

unsigned getSpectrumWithinWLRange(operaSpectralElements *inputSpectrum, double wl0, double wlf, double *outputFlux, double *outputWavelength) {
    if(!inputSpectrum->getHasWavelength()) {
        throw operaException("getSpectrumWithinWLRange: no wavelength in input SpectralElements. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    unsigned np = 0;
    for (unsigned i=0; i<inputSpectrum->getnSpectralElements(); i++) {
        if(inputSpectrum->getwavelength(i) >= wl0 &&
           inputSpectrum->getwavelength(i) <= wlf ) {
            
            outputFlux[np] = inputSpectrum->getFlux(i);
            outputWavelength[np] = inputSpectrum->getwavelength(i);
            np++;
        }
    }
    return np;
}

bool getOverlappingWLRange(operaSpectralElements *refElements, operaSpectralElements *elementsToMatch, double &wl0, double &wlf) {
    bool debug = false;
    double ref_wl0 = refElements->getwavelength(0);
    double ref_wlf = refElements->getwavelength(refElements->getnSpectralElements()-1);
    double elem2match_wl0 = elementsToMatch->getwavelength(0);
    double elem2match_wlf = elementsToMatch->getwavelength(elementsToMatch->getnSpectralElements()-1);
    
    if(ref_wl0 > ref_wlf || elem2match_wl0 > elem2match_wlf) {
        throw operaException("getOverlappingWLRange: initial wl must not be greater than final wl. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    // find overlapping range:
    bool overlap = false;
    double intersect_wl0 = 0;
    double intersect_wlf = 0;
    
    if(elem2match_wl0 >= ref_wl0 && elem2match_wl0 <= ref_wlf) {
        intersect_wl0 = elem2match_wl0;
        intersect_wlf = ref_wlf;
        overlap = true;
    }
    
    if (elem2match_wlf >= ref_wl0 && elem2match_wlf <= ref_wlf) {
        intersect_wlf = elem2match_wlf;
        if(intersect_wl0 == 0) {
            intersect_wl0 = ref_wl0;
        }
        overlap = true;
    }
    
    wl0 = intersect_wl0;
    wlf = intersect_wlf;
    
    if(debug) {
        cout << "refElements     = " << ref_wl0 << " -> " << ref_wlf << endl;
        cout << "elementsToMatch = " << elem2match_wl0 << " -> " << elem2match_wlf << endl;
        if (overlap == true) {
            cout << "Intersection    = " << intersect_wl0 << " -> " << intersect_wlf << endl;
        } else {
            cout << " Elements do not overlap" << endl;
        }
    }
    
    return overlap;
}

double calculateDeltaRadialVelocityInKPS(double telluricWL, double observedWL) {
    return (telluricWL - observedWL)*SPEED_OF_LIGHT_KMS/telluricWL;
}

operaVector fitSpectrum(const operaVector& inputWavelength, const operaVector& inputFlux, const operaVector& outputWavelength) {
	operaVector outputFlux(outputWavelength.size());
	operaFitSplineDouble(inputWavelength.size(), inputWavelength.datapointer(), inputFlux.datapointer(), outputWavelength.size(), outputWavelength.datapointer(), outputFlux.datapointer());
	return outputFlux;
}

operaSpectrum fitSpectrum(const operaSpectrum& inputSpectrum, const operaVector& outputWavelength) {
    operaSpectrum outputSpectrum(inputSpectrum.fluxcount(), outputWavelength);
    for(unsigned v=0; v<inputSpectrum.fluxcount(); v++) {
		operaFitSplineDouble(inputSpectrum.size(), inputSpectrum.wavelength_ptr(), inputSpectrum.flux_ptr(v), outputSpectrum.size(), outputSpectrum.wavelength_ptr(), outputSpectrum.flux_ptr(v));
		operaFitSplineDouble(inputSpectrum.size(), inputSpectrum.wavelength_ptr(), inputSpectrum.variance_ptr(v), outputSpectrum.size(), outputSpectrum.wavelength_ptr(), outputSpectrum.variance_ptr(v));
	}
    return outputSpectrum;
}

operaSpectrum getSpectrumWithinMask(string wavelength_mask, const operaSpectrum& inputSpectrum) {
    operaWavelengthRanges wlranges = readContinuumWavelengthMask(wavelength_mask);
    return getSpectrumWithinRange(wlranges, inputSpectrum);
}

/*
 * Remove spectral points around lines -- this is useful for masking telluric lines
 */
operaSpectrum maskSpectrumAroundLines(const operaSpectrum inputSpectrum, const operaSpectrum telluricLines, double spectralResolution, double nsig) {
    
    operaSpectrum outputSpectrum;
    
    operaWavelengthRanges wlranges = getWavelengthMaskAroundLines(telluricLines,spectralResolution,nsig);
    
    operaWavelengthRanges invertedRanges;

    for (unsigned r=0; r<wlranges.size(); r++) {
        if (r==0) {
            invertedRanges.addWavelengthRange(-BIG,wlranges.wl0(r));
        } else if (r==wlranges.size()-1) {
            invertedRanges.addWavelengthRange(wlranges.wlf(r),BIG);
        } else {
            invertedRanges.addWavelengthRange(wlranges.wlf(r-1),wlranges.wl0(r));
        }
    }
    
    unsigned firstr = 0;
    
    for (unsigned i=0; i<inputSpectrum.size(); i++) {
        double wl = inputSpectrum.getwavelength(i);
        
        for (unsigned r=firstr; r<wlranges.size(); r++) {
            if(invertedRanges.contains(wl,r)) {
                outputSpectrum.insert(wl, inputSpectrum.getflux(i), inputSpectrum.getvariance(i));
                firstr = r;
                break;
            }
        }
    }

    return outputSpectrum;
}


operaVector convolveSpectrum(operaSpectrum inputSpectrum, double spectralResolution) {
    return convolveSpectrumWithGaussianByResolution(inputSpectrum.wavelengthvector(), inputSpectrum.fluxvector(), spectralResolution);
}

operaVector fitSpectrumToPolynomial(const operaVector& inputWavelength, const operaVector& inputFlux, const operaVector& outputWavelength, unsigned order) {
    operaSpectrum filteredinput;
    for(unsigned i = 0; i < inputWavelength.size(); i++) {
		if(!isnan(inputFlux[i])) filteredinput.insert(inputWavelength[i], inputFlux[i], 1.0);
	}
    operaVector coeffs(order+1);
    coeffs = 1.0;
    operaMPFitPolynomial(filteredinput.size(), filteredinput.wavelength_ptr(), filteredinput.flux_ptr(), filteredinput.variance_ptr(), coeffs.size(), coeffs.datapointer(), 0, 0);
    return Operation(Polynomial(coeffs), outputWavelength);
}

// Helper function for LinearFit
double LinearMedianDelta(const operaVector& x, const operaVector& y, double b, double& a, double& absdev, double eps) {
	operaVector d = y - (b * x);
	a = Median(d);
	d -= a; // d = y - (b * x + a)
	absdev = Sum(Abs(d));
	double sum = 0.0;
	for (unsigned i=0; i<d.size(); i++) {
		if (y[i] != 0.0) d[i] /= fabs(y[i]);
		if (fabs(d[i]) > eps) sum += (d[i] >= 0.0 ? x[i] : -x[i]);
	}
	return sum;
}

void LinearFit(const operaVector& x, const operaVector& y, double& a, double& b, double& absdev) {
	float eps = EPS;
	
	double sx = Sum(x);
	double sy = Sum(y);
	double sxy = InnerProduct(x, y);
	double sxx = InnerProduct(x, x);
	double del = x.size() * sxx - sx * sx;
	
	if (del == 0.0) { // All X's are the same
		b = Median(y);
		a = 0.0; // Bisect the range w/ a flat line
		return;
	}
	double aa = (sxx * sy - sx * sxy) / del; // Least squares solution y = aa + bb * x
	double bb = (x.size() * sxy - sx * sy) / del;
	
	const operaVector& t = y - (aa + bb * x);
	double sigb = Magnitude(t)/sqrt(del); // Standard deviation sig = sqrt(chisqr/del)
	
	double b1 = bb;
	double f1 = LinearMedianDelta(x, y, b1, aa, absdev, eps);
	
	//  Quick return. The initial least squares gradient is the LAD solution.
	if (f1 != 0.0) {
		double delb = (f1 >= 0 ? 3.0 : -3.0) * sigb;
		double b2 = b1 + delb;
		double f2 = LinearMedianDelta(x, y, b2, aa, absdev, eps);
		while (f1*f2 > 0.0) {     // Bracket the zero of the function
			b1 = b2;
			f1 = f2;
			b2 = b1 + delb;
			f2 = LinearMedianDelta(x, y, b2, aa, absdev, eps);
		}
		//  In case we finish early.
		bb = b2;
		double f = f2;
		
		// Narrow tolerance to refine 0 of fcn.
		sigb = 0.01 * sigb;
		while (fabs(b2-b1) > sigb && f != 0.0) { // bisection of interval b1,b2.
			bb = 0.5 * (b1 + b2);
			if (bb == b1 || bb == b2)
				break;
			f = LinearMedianDelta(x, y, bb, aa, absdev, eps);
			if (f*f1 >= 0.0) {
				f1 = f;
				b1 = bb;
			} else {
				f2 = f;
				b2 = bb;
			}
		}
	}
	absdev /= x.size();
	b = bb;
	a = aa;
}

void LinearFit(const operaVector& x, const operaVector& y, double& a, double& aError, double& b, double& bError, double& absdev) {
	LinearFit(x, y, a, b, absdev);
	double sigy = MedianStdDev(y - a - b*x, 0);
	double sx = Sum(x);
	double sxx = InnerProduct(x, x);
	double del = x.size() * sxx - sx * sx;
	aError = sigy * sqrt(sxx / del);
	bError = sigy * sqrt(x.size() / del);
}

Polynomial PolynomialFit(const operaVector& x, const operaVector& y, operaVector initcoeffs) {
	double chisqr;
	operaLMFitPolynomial(x.size(), x.datapointer(), y.datapointer(), initcoeffs.size(), initcoeffs.datapointer(), &chisqr);
	Polynomial poly(initcoeffs.size(), initcoeffs.datapointer());
	poly.setChisqr(chisqr);
	return poly;
}

Polynomial PolynomialFit(const operaVector& x, const operaVector& y, const operaVector& yerr, operaVector initcoeffs, operaVector initcoefferrors) {
	double chisqr;
	int errorcode = operaMPFitPolynomial(x.size(), x.datapointer(), y.datapointer(), yerr.datapointer(), initcoeffs.size(), initcoeffs.datapointer(), initcoefferrors.datapointer(), &chisqr);
	if (errorcode <= 0) {
		throw operaException("PolynomialFit: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);	
	}
	Polynomial poly(initcoeffs.size(), initcoeffs.datapointer(), initcoefferrors.datapointer());
	poly.setChisqr(chisqr);
	return poly;
}

unsigned findClosestInSortedRange(const operaVector& vector, double target, unsigned startindex, unsigned endindex) {
	if(startindex >= endindex || vector.begin()+endindex > vector.end()) return vector.size();
	std::vector<double>::const_iterator iter = std::lower_bound(vector.begin()+startindex, vector.begin()+endindex, target); //get the first element >= target
	if(iter == vector.end()) return vector.size() - 1;
	if(iter > vector.begin()+startindex && fabs(*iter - target) >= fabs(*(iter-1) - target)) return iter - 1 - vector.begin();
	return iter - vector.begin();
}

void measureFluxContinuum(const operaFluxVector& uncalibratedFlux, unsigned binsize, unsigned nsigcut, operaVector& continuumElemSamples, operaVector& continuumFluxSamples) {
	const unsigned NumberofPoints = uncalibratedFlux.getlength();
	if(NumberofPoints < 2*binsize) throw operaException("NumberofPoints < 2*binsize", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    
    for(unsigned start=0; start<NumberofPoints; start+=binsize){
        // If our bin runs past the end, shift it back
        unsigned end = start+binsize;
        if(end > NumberofPoints) {
            start = NumberofPoints - binsize;
            end = NumberofPoints;
        }
        // Add an extra binsize worth of points to both sides for our first calcuations, but if we are near the edge only add one extra binsize
		unsigned firstPoint, lastPoint;
		if(start < binsize) {
            firstPoint = 0;
            lastPoint = 2*binsize;
        } else if(end == NumberofPoints) {
            firstPoint = NumberofPoints - 2*binsize;
            lastPoint = NumberofPoints;
        } else {
            firstPoint = start - binsize;
            lastPoint = end + binsize;
            if (lastPoint > NumberofPoints) lastPoint = NumberofPoints; // Othewise we can run past the end on the 2nd to last bin
        }
		
		operaVector uncalflux_tmp;
        operaVector elemindex_tmp;
        for(unsigned i=firstPoint;i<lastPoint;i++) {
            uncalflux_tmp.insert(uncalibratedFlux.getflux(i));
            elemindex_tmp.insert((double)i);
        }
        
        double am, bm, abdevm;
		// 1st pass: calculate continuum slope
		LinearFit(elemindex_tmp, uncalflux_tmp, am, bm, abdevm); // robust linear fit: f(x) =  a + b*x
        // Filter out points that deviate more than abdev from the robust linear fit
        uncalflux_tmp.clear();
        elemindex_tmp.clear();
        for(unsigned i=firstPoint;i<lastPoint;i++) {
            if(fabs(uncalibratedFlux.getflux(i) - (bm*i + am)) < abdevm) {
                uncalflux_tmp.insert(uncalibratedFlux.getflux(i));
                elemindex_tmp.insert((double)i);
            }
        }
        // 2nd pass: calculate continuum slope    
        if(elemindex_tmp.size() > 0) {
            LinearFit(elemindex_tmp, uncalflux_tmp, am, bm, abdevm); // robust linear fit: f(x) =  a + b*x
        }
        
        // Calculate residuals in our original binsize and sort
        operaVector residuals_tmp;
        for(unsigned i=start; i<end; i++) residuals_tmp.insert(uncalibratedFlux.getflux(i) - (bm*i + am));
        residuals_tmp.sort();
        
        // Select largest residual which does not exceed nsigcut*(absolute deviation)
        double dytop = nsigcut*abdevm; // if none are found then use nsigcut*abdev
        double continuumBinFraction = 0.5;
        for(unsigned i=binsize-1; i>(unsigned)(binsize*(1.0 - continuumBinFraction)); i--) {
            if(residuals_tmp[i] < nsigcut*abdevm) {
                dytop = residuals_tmp[i];
                break;
            }
        }
        
        // If this is the first or last bin, sample at the edge of the bin, otherwise sample at the center of the bin
        if(start==0) {
            continuumElemSamples.insert(start);
            continuumFluxSamples.insert(bm*start + am + dytop);
        }
        continuumElemSamples.insert(start + binsize/2.0);
        continuumFluxSamples.insert(bm*(start + binsize/2.0) + am + dytop);
        if (end==NumberofPoints) {
			continuumElemSamples.insert(end-1);
			continuumFluxSamples.insert(bm*(end-1) + am + dytop);
        }
	}
}
