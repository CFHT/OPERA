#ifndef OPERASPECTRALELEMENTS_H
#define OPERASPECTRALELEMENTS_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralElements
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

#include "libraries/operaLibCommon.h" //for operaSpectralOrder_t
#include "libraries/operaFluxVector.h"
#include "libraries/operaWavelength.h"

/*! 
 * \sa class operaSpectralElements
 * \brief operaSpectralElements
 * \details A spectral element (SE) consists of a spectrally resolved 
 * \details subdivision of a spectral order taken in the dispersion direction i.e. along the trace polynomial
 * \details more or less in the y direction for espadons, but is actually curved. 
 * \details The minimal width of the spectral element is defined by the resolution 
 * \details of the spectrograph (RE). The actual width could be defined as any size 
 * \details smaller than the order size and as long as RE < RS, where RE is the 
 * \details resolution of the element.
 * \return none
 * \file operaSpectralElements.h
 * \ingroup libraries
 */

class operaSpectralElements {
	
private:
	unsigned nSpectralElements;			// a count of the spectral elements stored
	double elementHeight;				// height of spectral element in pixel units measured along the dispersion direction (CenterPolynomial)
	operaSpectralOrder_t SpectrumType;	// what kind of spectrum this is

	operaFluxVector fluxvector;			// (dbl in counts)
	operaVector XCorrelation;			// cross-correlation vector, used for detection of spectral lines
	operaVector photoCenterX;			// x coordinates of the photocenter in the image reference frame (pixel units) 
	operaVector photoCenterY;			// y coordinates of the photocenter in the image reference frame (pixel units) 
	operaVector distd;					// distance in pixel units measured along the dispersion, starting at the first element d[0] = 0 
	operaVector wavelength;				// wavelength in nm 
	operaVector fluxSNR;				// SNR at each spectralelement
	
	bool hasRawFlux;
	bool hasStandardFlux;
	bool hasOptimalFlux;
	bool hasOperaOptimalFlux;
	
	bool hasXCorrelation;    
	bool hasWavelength;    
	bool hasDistance;    
	bool hasFluxSNR;    
	
	// Extended spectrum includes variations of wl and flux data...
	operaVector tell;					// telluric wl information
	operaVector rvel;					// heliocentric wl information
	
	operaFluxVector rawFlux;			// raw flux
	operaFluxVector normalizedFlux;		// normalized flux
	operaFluxVector fcalFlux;			// flux calibrated flux
	
	bool hasExtendedBeamFlux;
    
public:
	
	/*!
     * \brief Creates empty set of operaSpectralElements.
     */
	operaSpectralElements();
	
	/*!
     * \brief Creates operaSpectralElements with the specfied number of elements.
     * \param nElements The number of spectral elements.
     */
	operaSpectralElements(unsigned nElements);
	
	/*!
     * \brief Creates operaSpectralElements as specfied by parameters.
     * \param nElements The number of spectral elements.
     * \param format A flag specifying The format of the spectrum.
     * \param extended Whether or not to create extended wavelength and flux vectors (default is false).
     */
	operaSpectralElements(unsigned nElements, operaSpectralOrder_t format, bool extended = false);
	
	/*!
     * \brief Creates extended wavelength and flux vectors.
     * \details Creates vectors containing wavelength corrections for telluric and heliocentric radial velocity
     * \details as well as individual flux vectors for the raw, flux calibrated, and normalized flux.
     */
	void createExtendedVectors(); 
	
	/*!
     * \brief Resizes to the specified number of spectral elements.
     * \param nElements The new number of spectral elements.
     * \details Existing elements are preserved up to nElements.
     */
	void resize(unsigned nElements);
	
	/*!
     * \brief Removes all spectral elments outside of the specified range.
     * \param range The index range of elements to keep.
     */
	void trim(operaIndexRange range);
	
	/*!
     * \brief Gets the number of spectral elements.
     * \return The number of spectral elements.
     */
	unsigned getnSpectralElements(void) const;
	
	/*!
     * \brief Gets the main flux vector.
     * \return An immutable reference to the flux vector.
     */
	const operaFluxVector& getFluxVector(void) const;
	
	/*!
     * \brief Gets the cross correlation vector.
     * \return An immutable reference to the XCorrelation vector.
     */
	const operaVector& getXCorrelation() const;
	
	/*!
     * \brief Gets the vector of photocenter x-coordinates.
     * \return An immutable reference to the photoCenterX vector.
     */
	const operaVector& getPhotoCenterX() const;
	
	/*!
     * \brief Gets the vector of photocenter y-coordinates.
     * \return An immutable reference to the photoCenterY vector.
     */
	const operaVector& getPhotoCenterY() const;
	
	/*!
     * \brief Gets the distance vector.
     * \return An immutable reference to the distd vector.
     */
	const operaVector& getDistd() const;
	
	/*!
     * \brief Gets the wavelength vector.
     * \return An immutable reference to the wavelength vector.
     */
	const operaVector& getWavelength(void) const;
	
	/*!
     * \brief Gets the fluxSNR vector.
     * \return An immutable reference to the fluxSNR vector.
     */
	const operaVector& getFluxSNR() const;
	
	/*!
     * \brief Gets the wavelength correction vector.
     * \return An immutable reference to the rvel vector.
     */
	const operaVector& getRvel(void) const;
	
	/*!
     * \brief Sets the main flux vector.
     * \param FluxVector A flux vector with same number of spectral elements.
     */
	void setFluxVector(const operaFluxVector &FluxVector);
	
	/*!
     * \brief Sets the cross correlation vector.
     * \param XCorr A cross correlation vector with same number of spectral elements.
     */
	void setXCorrelation(const operaVector& XCorr);
	
	/*!
     * \brief Sets the photo center vectors.
     * \param PhotoCenterX A photocenter x-coordinate vector with same number of spectral elements.
     * \param PhotoCenterY A photocenter y-coordinate vector with same number of spectral elements.
     */
	void setPhotoCenter(const operaVector& PhotoCenterX, const operaVector& PhotoCenterY);
	
	/*!
     * \brief Sets the distance vector.
     * \param Distd A distance vector with same number of spectral elements.
     */
	void setDistd(const operaVector& Distd);
	
	/*!
     * \brief Sets the wavelength vector.
     * \param Wavelength A wavelength vector with same number of spectral elements.
     */
	void setWavelength(const operaVector &Wavelength);
	
	/*!
     * \brief Sets the fluxSNR vector.
     * \param FluxSNR A fluxSNR vector with same number of spectral elements.
     */
	void setFluxSNR(const operaVector& FluxSNR);
	
	/*!
     * \brief Sets the wavelength correction vector.
     * \param Rvel A wavelength correction vector with same number of spectral elements.
     */
	void setRvel(const operaVector &Rvel);
	
	/*!
     * \brief Calculates the flux SNR from the flux vector.
     * \details Updates the flux SNR vector with the calculated values.
     */
	void calculateFluxSNR();
	
	/*!
     * \brief Set the wavelengths using the given wavelength calibration.
     * \param Wavelength The wavelength polynomial to calculate wavelengths with.
     * \details Updates the wavelength vector with the calculated value.
     */
	void setwavelengthsFromCalibration(const operaWavelength *Wavelength);
	
	/*!
     * \brief Gets the index range of spectral elements in the specified wavelength range.
     * \param wl0 The lower bound of the wavelength range.
     * \param wlf The upper bound of the wavelength range.
     * \details The wavelength vector must be in ascending order, and wl0 must not be greater than wlf.
     * \return The index range of wavelengths between the bounds.
     */
	operaIndexRange getContainedWavelengthSubrangeIndexes(double wl0, double wlf) const;
    
    /*!
     * \brief Gets the flux of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The flux at that index.
     */
	double getFlux(unsigned indexElem) const;
	
	/*!
     * \brief Sets the flux of a spectral element.
     * \param Flux The new flux value.
     * \param indexElem The index of the spectral element.
     */
	void setFlux(double Flux, unsigned indexElem);
	
	/*!
     * \brief Gets the flux variance of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The flux variance at that index.
     */
	double getFluxVariance(unsigned indexElem) const;
	
	/*!
     * \brief Sets the flux variance of a spectral element.
     * \param FluxVariance The new flux variance value.
     * \param indexElem The index of the spectral element.
     */
	void setFluxVariance(double FluxVariance, unsigned indexElem);
	
	/*!
     * \brief Gets the flux SNR of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The flux SNR at that index.
     */
	double getFluxSNR(unsigned indexElem) const;
	
	/*!
     * \brief Sets the flux SNR of a spectral element.
     * \param FluxSNR The new flux SNR value.
     * \param indexElem The index of the spectral element.
     */
	void setFluxSNR(double FluxSNR, unsigned indexElem);
	
	/*!
     * \brief Gets the x-coordinate of the photocenter a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The x-coordinate (in pixels) of the photocenter at that index.
     */
	double getphotoCenterX(unsigned indexElem) const;
	
	/*!
     * \brief Gets the y-coordinate of the photocenter a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The y-coordinate (in pixels) of the photocenter at that index.
     */
	double getphotoCenterY(unsigned indexElem) const;
	
	/*!
     * \brief Sets the coordinates of the photocenter of a spectral element.
     * \param x The new x-coordinate (in pixels).
     * \param y The new y-coordinate (in pixels).
     * \param indexElem The index of the spectral element.
     */
	void setphotoCenter(double x, double y, unsigned indexElem);
	
	/*!
     * \brief Gets the distance along the dispersion direction of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The distance (in pixels) at that index.
     */
	double getdistd(unsigned indexElem) const;
	
	/*!
     * \brief Sets the distance along the dispersion direction of a spectral element.
     * \param Distd The new distance (in pixels).
     * \param indexElem The index of the spectral element.
     */
	void setdistd(double Distd, unsigned indexElem);
	
	/*!
     * \brief Gets the wavelength of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The wavelength at that index.
     */
	double getwavelength(unsigned indexElem) const;
	
	/*!
     * \brief Sets the wavelength of a spectral element.
     * \param Wavelength The new wavelength.
     * \param indexElem The index of the spectral element.
     */
	void setwavelength(double Wavelength, unsigned indexElem);
    
	/*!
     * \brief Gets the cross-correlation of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The cross-correlation at that index.
     */
	double getXCorrelation(unsigned indexElem) const;  
	
	/*!
     * \brief Sets the cross-correlation of a spectral element.
     * \param Xcorr The new cross-correlation.
     * \param indexElem The index of the spectral element.
     */
	void setXCorrelation(double Xcorr, unsigned indexElem);    

	/*!
     * \brief Gets the telluric wavelength correction of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The telluric correction at that index.
     */
	double gettell(unsigned indexElem) const;  
	/*!
     * \brief Sets the telluric wavelength correction of a spectral element.
     * \param value The new telluric correction.
     * \param indexElem The index of the spectral element.
     */
	void settell(double value, unsigned indexElem);
	/*!
     * \brief Copies the wavelength vector into the telluric wavelength vector.
     * \details Replaces the telluric wavelength correction of each spectral element with the wavelength of that spectral element.
     */
	void copyTOtell(void);
	/*!
     * \brief Copies the telluric wavelength vector into the wavelength vector.
     * \details Replaces the wavelength of each spectral element with the telluric wavelength correction of that spectral element.
     */
	void copyFROMtell(void);
	
	/*!
     * \brief Gets the radial velocity wavelength correction of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The radial velocity correction at that index.
     */
	double getrvel(unsigned indexElem) const;  
	/*!
     * \brief Sets the radial velocity wavelength correction of a spectral element.
     * \param value The new radial velocity correction.
     * \param indexElem The index of the spectral element.
     */
	void setrvel(double value, unsigned indexElem);    
	
	/*!
     * \brief Gets the normalized flux of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The normalized flux at that index.
     */
	double getnormalizedFlux(unsigned indexElem) const;
	/*!
     * \brief Sets the normalized flux of a spectral element.
     * \param Flux The new normalized flux value.
     * \param indexElem The index of the spectral element.
     */
	void setnormalizedFlux(double value, unsigned indexElem);
	/*!
     * \brief Gets the normalized flux variance of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The normalized flux variance at that index.
     */
	double getnormalizedFluxVariance(unsigned indexElem) const;
	/*!
     * \brief Sets the normalized flux variance of a spectral element.
     * \param FluxVariance The new normalized flux variance value.
     * \param indexElem The index of the spectral element.
     */
	void setnormalizedFluxVariance(double value, unsigned indexElem);
	/*!
     * \brief Copies the flux vector into the normalized flux vector.
     */
	void copyTOnormalizedFlux(void);
	/*!
     * \brief Copies the normalized flux vector into the flux vector.
     */
	void copyFROMnormalizedFlux(void);
	
	/*!
     * \brief Gets the calibrated flux of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The calibrated flux at that index.
     */
	double getfcalFlux(unsigned indexElem) const;
	/*!
     * \brief Sets the calibrated flux of a spectral element.
     * \param Flux The new calibrated flux value.
     * \param indexElem The index of the spectral element.
     */
	void setfcalFlux(double value, unsigned indexElem);
	/*!
     * \brief Gets the calibrated flux variance of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The calibrated flux variance at that index.
     */
	double getfcalFluxVariance(unsigned indexElem) const;    
    /*!
     * \brief Sets the calibrated flux variance of a spectral element.
     * \param FluxVariance The new calibrated flux variance value.
     * \param indexElem The index of the spectral element.
     */
	void setfcalFluxVariance(double value, unsigned indexElem);    
	/*!
     * \brief Copies the flux vector into the calibrated flux vector.
     */
	void copyTOfcalFlux(void);
	/*!
     * \brief Copies the calibrated flux vector into the flux vector.
     */
	void copyFROMfcalFlux(void);
	
	/*!
     * \brief Gets the raw flux of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The raw flux at that index.
     */
	double getrawFlux(unsigned indexElem) const;  
	/*!
     * \brief Sets the raw flux of a spectral element.
     * \param Flux The new raw flux value.
     * \param indexElem The index of the spectral element.
     */
	void setrawFlux(double value, unsigned indexElem);
	/*!
     * \brief Gets the raw flux variance of a spectral element.
     * \param indexElem The index of the spectral element.
     * \return The raw flux variance at that index.
     */
	double getrawFluxVariance(unsigned indexElem) const;
    /*!
     * \brief Sets the raw flux variance of a spectral element.
     * \param FluxVariance The new raw flux variance value.
     * \param indexElem The index of the spectral element.
     */
	void setrawFluxVariance(double value, unsigned indexElem);
	/*!
     * \brief Copies the flux vector into the raw flux vector.
     */
	void copyTOrawFlux(void);
	/*!
     * \brief Copies the raw flux vector into the flux vector.
     */
	void copyFROMrawFlux(void);
	
	operaSpectralOrder_t getSpectrumType(void) const { return SpectrumType; }
	void setSpectrumType(operaSpectralOrder_t format) { SpectrumType = format; }
	
	double getelementHeight() const { return elementHeight; }
	void setelementHeight(double Height) { elementHeight = Height; }
	
	bool getHasRawSpectrum() const { return hasRawFlux; };
	void setHasRawSpectrum(bool HasRawFlux) { hasRawFlux = HasRawFlux; };
	
	bool getHasStandardSpectrum() const { return hasStandardFlux; };
	void setHasStandardSpectrum(bool HasStandardFlux) { hasStandardFlux = HasStandardFlux; };
	
	bool getHasOptimalSpectrum() const { return hasOptimalFlux; };
	void setHasOptimalSpectrum(bool HasOptimalFlux) { hasOptimalFlux = HasOptimalFlux; };
	
	bool getHasOperaOptimalSpectrum() const { return hasOperaOptimalFlux; };
	void setHasOperaOptimalSpectrum(bool HasOperaOptimalFlux) { hasOperaOptimalFlux = HasOperaOptimalFlux; };
	
	bool getHasXCorrelation() const { return hasXCorrelation; };
	void setHasXCorrelation(bool HasXCorrelation) { hasXCorrelation = HasXCorrelation; };
    
	bool getHasWavelength() const { return hasWavelength; };
	void setHasWavelength(bool HasWavelength) { hasWavelength = HasWavelength; };
    
	bool getHasDistance() const { return hasDistance; };
	void setHasDistance(bool HasDistance) { hasDistance = HasDistance; };
    
	bool getHasFluxSNR() const { return hasFluxSNR; };
	void setHasFluxSNR(bool HasFluxSNR) { hasFluxSNR = HasFluxSNR; };
    
	bool getHasExtendedBeamFlux() const { return hasExtendedBeamFlux; };
	void setHasExtendedBeamFlux(bool HasExtendedBeamFlux) { hasExtendedBeamFlux = HasExtendedBeamFlux; };
};

#endif
