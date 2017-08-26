#ifndef OPERASPECTRALLINES_H
#define OPERASPECTRALLINES_H

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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#ifndef MAXNUMBEROFFEATURES
#define MAXNUMBEROFFEATURES 1000
#endif

class operaSpectralOrder;
class operaSpectralElements;
class operaSpectralFeature;

enum dispersionaxis_t {wavelength_disp, distance_disp, y_distance_disp};

/*! 
 * \class operaSpectralLines
 * \brief This class is used to manipulate a set of spectral lines within an order.
 * \return none
 * \file operaSpectralLines.h
 * \ingroup libraries
 */
class operaSpectralLines {
	
private:
    unsigned nLines; // total number of lines detected
    
    operaSpectralElements *comparisonSpectrum; // pointer to the operaSpectralElements containing comparison extracted spectrum

	unsigned nFeatures; // number of Features in the spectrum
    
    operaSpectralFeature *spectralFeatures[MAXNUMBEROFFEATURES]; // vector of pointers to operaSpectralFeature
    
    double referenceLineWidth; // width to be used as reference
    
    dispersionaxis_t dispersiontype;
    
public:
	
	/*
	 * Constructor
	 */
	operaSpectralLines(void);
	
    operaSpectralLines(operaSpectralElements *ComparisonSpectrum, double ReferenceLineWidth);
    
    operaSpectralLines(operaSpectralElements *ComparisonSpectrum, double ReferenceLineWidth, dispersionaxis_t Dispersiontype);
    
    operaSpectralLines(unsigned Features, unsigned Linesinfeature);
    
    operaSpectralLines(unsigned Features, unsigned Linesinfeature, dispersionaxis_t Dispersiontype);    

	/*
	 * Destructor
	 */
	~operaSpectralLines(void);
	
	/*
	 * Setters/Getters	
     */

    void setDispersiontype(dispersionaxis_t Dispersiontype);   

    dispersionaxis_t getDispersiontype(void) const;
    
    void setnLines(unsigned NLines);
    
    unsigned getnLines(void) const;

    unsigned getNFeatures(void) const;
    
    void setNFeatures(unsigned NFeatures);
    
    void setReferenceLineWidth(double ReferenceLineWidth);
    
    double getReferenceLineWidth(void) const;
    
#if 0    
    operaSpectralElements *getcomparisonSpectrum(void);
      
    void setSpectralFeature(operaSpectralFeature &SpectralFeature, unsigned indexFeature);
#endif    
    const operaSpectralFeature *getSpectralFeature(unsigned indexFeature) const;
    
    operaSpectralFeature *getSpectralFeature(unsigned indexFeature);
	
    /*
	 * Methods	
     */
         
    void detectSpectralFeatures(double DetectionThreshold, double LocalMaxFilterWidth, double MinPeakDepth);
    
    void detectAbsorptionSpectralFeatures(double DetectionThreshold, double LocalMaxFilterWidth, double MinPeakDepth);
    
    void subtractFeatureModel(void); 
    
    void printLines(ostream *pout) const;
    
    void printReferenceSpectrum(ostream *pout) const;
    
    unsigned selectLines(double MaxContamination, unsigned nSig, double amplitudeCutOff, double *LinePositionVector,double *LineSigmaVector,double *LineAmplitudeVector);
    
	// Note that these are length protected...
	inline double getMX(unsigned k) const;
	
	inline double getMY(unsigned k) const;
	
    inline double getAbsorptionMY(unsigned k) const;
    
	inline void setMY(double value, unsigned k);
	
	inline double getVar(unsigned k) const;
	
	inline void setVar(double value, unsigned k);
	
};

#endif
