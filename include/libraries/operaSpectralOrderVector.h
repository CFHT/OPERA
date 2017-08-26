#ifndef OPERASPECTRALORDERVECTOR_H
#define OPERASPECTRALORDERVECTOR_H
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

#include <vector>
#include "operaError.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralElements
#include "libraries/operaSpectralTools.h"			// for operaSpectrum
#include "libraries/operaSpectralOrder.h"
#include "libraries/GainBiasNoise.h"
#include "libraries/Polynomial.h"	
#include "libraries/LaurentPolynomial.h"	
#include "libraries/operaWavelength.h" // for MAXORDEROFWAVELENGTHPOLYNOMIAL
#include "libraries/operaEspadonsImage.h" // for instrumentmode_t

#define NEWLINES_BETWEEN_ORDERS false

#define MAXNUMBEROFPOINTSINFLATRESPONSE 1000

using namespace std;

/*! 
 * \sa class operaSpectralOrderVector
 * \brief Stores an ordered vector of SpectralOrder class instances. 
 * \return none
 * \file operaSpectralOrderVector.h
 * \ingroup libraries
 */
class operaSpectralOrderVector {
	
private:
	std::vector<operaSpectralOrder*> vector; // a vector of spectral order pointers
	Polynomial *orderSpacingPolynomial; // captures the spacing between orders
	
	string object;					// for Libre-Esprit output
    
    unsigned numberOfDispersionPolynomials;
	LaurentPolynomial *dispersionPolynomial[MAXORDEROFWAVELENGTHPOLYNOMIAL]; // captures the dispersion polynomials
	
    unsigned length;				// total length of vector
	unsigned minorder;				// lowest order number with a value
	unsigned maxorder;				// highest order number with a value
	unsigned sequence;				// for master fluxcalibrations to store sequence #
	instrumentmode_t instrumentmode;	// instrumentmode
	unsigned count;					// Note that this count is the count of actual orders, 
									// not the length of the vector (which will be MAXORDERS)
	GainBiasNoise *gainBiasNoise;	// flat stats
	
public:
	/*
	 * Constructors
	 */
	
	/*! 
	 * \sa class operaSpectralOrderVector();
	 * \brief Base constructor.
	 * \return void
	 */
	operaSpectralOrderVector();
	/*! 
	 * \sa class operaSpectralOrderVector(unsigned length, operaSpectralOrder_t OrderType, unsigned maxdatapoints, unsigned maxValues, unsigned nElements);
	 * \brief Create a SpectralOrderVector of a given type.
	 * \return void
	 */
	operaSpectralOrderVector(unsigned length, unsigned maxdatapoints, unsigned maxValues, unsigned nElements);
	/*
	 * Destructor
	 */
	~operaSpectralOrderVector();

	/*
	 * Methods
	 */
	
	/*!
	 * \sa method unsigned getnumberOfDispersionPolynomials(void);
	 * \brief returns the number of dispersion polynomials.
	 * \return unsigned - numberOfDispersionPolynomials.
	 */
	unsigned getnumberOfDispersionPolynomials(void) const;
	
	/*!
	 * \sa method void setnumberOfDispersionPolynomials(unsigned NumberOfDispersionPolynomials);
	 * \brief sets the number of dispersion polynomials.
	 * \return none.
	 */
	void setnumberOfDispersionPolynomials(unsigned NumberOfDispersionPolynomials);
    
	/*! 
	 * \sa method void freeoperaSpectralOrderVector();
	 * \brief frees the vector contents and the vector.
	 * \return none.
	 */
	void freeSpectralOrderVector();
	/*! 
	 * \sa method unsigned getGainBiasNoise();
	 * \brief returns a pointer to the GainBiasNoise class instance.
	 * \return GainBiasNoise pointer.
	 */
	GainBiasNoise *getGainBiasNoise();
	const GainBiasNoise *getGainBiasNoise() const;
	
	void setWavelengthsFromCalibration(int Minorder, int Maxorder);
	
	/*! 
	 * \sa method unsigned getCount();
	 * \brief returns the count of spectral orders that have content.
	 * \return unsigned - count.
	 */
	unsigned getCount() const;
	
	/*! 
	 * \sa method unsigned getMinorder();
	 * \brief returns the least order number in the vector
	 * \return unsigned - minorder.
	 */
	unsigned getMinorder() const;
	
	/*! 
	 * \sa method unsigned getMaxorder();
	 * \brief returns the maximal order number in the vector
	 * \return unsigned - minorder.
	 */
	unsigned getMaxorder() const;
	
	/*! 
	 * \sa method unsigned setCount();
	 * \brief sets the count of spectral orders that have content.
	 * \return none.
	 */
	void setCount(unsigned Count);
	
	/*! 
	 * \sa method void shiftOrdersDown(unsigned shift);
	 * \brief decreases the index and label of each order by the number specified
	 * \return none.
	 */
	void shiftOrdersDown(unsigned shift);
	
	/*! 
	 * \sa method void shiftOrdersUp(unsigned shift);
	 * \brief increases the index and label of each order by the number specified
	 * \return none.
	 */
	void shiftOrdersUp(unsigned shift);
	
	/*! 
	 * \sa method void shiftOrders();
	 * \brief moves the index and label of each order by the number specified.
	 * \return none.
	 */
	void shiftOrders(int shift);
	
	/*! 
	 * \sa void setObject(string object)
	 * \brief sets the object name
	 * \return none.
	 */
	void setObject(string object);
	
	/*! 
	 * \sa string getObject(void);
	 * \brief get the object name
	 * \return none.
	 */
	string getObject(void) const;
	
	/*! 
	 * \sa void setInstrumentmode(instrumentmode_t Instrumentmode)
	 * \brief sets the Instrumentmode
	 * \return none.
	 */
	void setInstrumentmode(instrumentmode_t Instrumentmode) { instrumentmode = Instrumentmode; }
	
	/*! 
	 * \sa string getInstrumentmode(void);
	 * \brief get the getInstrumentmode
	 */
	instrumentmode_t getInstrumentmode(void) const { return instrumentmode; }
	
	/*! 
	 * \sa setSequence(unsigned sequence)
	 * \brief sets the sequence number
	 * \return none.
	 */
	void setSequence(unsigned Sequence);
	
	/*! 
	 * \sa unsigned getSequence(void);
	 * \brief get the sequence number
	 * \return none.
	 */
	unsigned getSequence(void) const;
	
	/*! 
	 * \sa method unsigned setMinorder();
	 * \brief sets the least order number in the vector
	 * \return none.
	 */
	void setMinorder(unsigned Minorder);
	
	/*! 
	 * \sa method void setMaxorder(unsigned Maxorder);
	 * \brief sets the maximal order number in the vector
	 * \return none.
	 */
	void setMaxorder(unsigned Maxorder);
	
	/*! 
	 * \sa method Polynomial *getOrderSpacingPolynomial(void);
	 * \brief gets the order spacing polynomial
	 * \return Polynomial *.
	 */
	Polynomial *getOrderSpacingPolynomial(void);
	const Polynomial *getOrderSpacingPolynomial(void) const;
	
	/*! 
	 * \sa method setOrderSpacingPolynomial(PolynomialCoeffs_t *pc);
	 * \brief sets the order spacing polynomial
	 * \return void.
	 */
	void setOrderSpacingPolynomial(PolynomialCoeffs_t *pc);
	
	/*!
	 * \sa method Polynomial *getDispersionPolynomial(unsigned index);
	 * \brief gets the dispersion polynomial
	 */
	LaurentPolynomial *getDispersionPolynomial(unsigned index);
	const LaurentPolynomial *getDispersionPolynomial(unsigned index) const;
	
	/*!
	 * \sa method setDispersionPolynomial(unsigned index, const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, PolynomialCoeffs_t *pc);
	 * \brief sets the dispersion polynomial
	 */
    void setDispersionPolynomial(unsigned index, const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, PolynomialCoeffs_t *pc);
	
	/*! 
	 * \sa method operaSpectralOrder* operaSpectralOrderVector::GetSpectralOrder(unsigned order);
	 * \brief Gets an operaSpectralOrder* to a given order, else NULL
	 * \return operaSpectralOrder - pointer to the operaSpectralOrder.
	 */
	operaSpectralOrder* GetSpectralOrder(unsigned order);
	const operaSpectralOrder* GetSpectralOrder(unsigned order) const;

    void fitOrderSpacingPolynomial(operaFITSImage &masterFlatImage, operaFITSImage &badpixImage, float slit, unsigned nsamples, unsigned sampleCenterPosition, unsigned referenceOrderNumber, float referenceOrderSeparation, int detectionMethod, bool FFTfilter, float gain, float noise, unsigned x1, unsigned x2, unsigned y1, unsigned y2, unsigned cleanbinsize, float nsigcut, ostream *pout);
    
    void measureIPAlongRowsFromSamples(operaFITSImage &masterFlatImage, operaFITSImage &badpixImage, float slit, unsigned nsamples, bool FFTfilter, float gain, float noise, unsigned x1, unsigned x2, unsigned y1, unsigned y2,float *ipfunc, float *ipx, float *iperr);
    
    unsigned getElemIndexAndOrdersByWavelength(int *orderWithReferenceFluxForNormalization, unsigned *elemIndexWithReferenceFluxForNormalization, double wavelength);

    unsigned getOrdersByWavelengthRange(int *orderForWavelengthRange, double Range_wl0, double Range_wlf);
    void getOrdersByWavelengthRange(operaWavelengthRange wavelengthRange, int& Minorder, int& Maxorder);

    void measureContinuumAcrossOrders(unsigned binsize, int orderBin, unsigned nsigcut);

    void measureContinuumAcrossOrders(unsigned binsize, int orderBin, unsigned nsigcut, unsigned nOrdersPicked, int *orderForWavelength);
    
    void FitFluxCalibrationAcrossOrders(int lowOrderToClip, int highOrderToClip, int orderBin, bool throughput);
    
    void getContinuumFluxesForNormalization(double *uncalibratedContinuumFluxForNormalization, double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS],unsigned binsize, int orderBin, unsigned nsigcut);

	unsigned getLEElementCount(string LEfluxCalibration);
	
	void readLEFluxCalibration(string LEfluxCalibration, operaSpectralElements *fluxCalibrationElements);
    
    /// Returns a clean uniform sample of the continuum flux for the spectral elements (and optionally the beam elements) of all orders
    operaSpectrum calculateCleanUniformSampleOfContinuum(int Minorder, int Maxorder, unsigned binsize, double delta_wl, string inputWavelengthMaskForUncalContinuum, unsigned numberOfPointsInUniformSample, bool useBeams);
    // Overload using float arrays instead of operaSpectrum
    void calculateCleanUniformSampleOfContinuum(int Minorder, int Maxorder, unsigned binsize, double delta_wl, string inputWavelengthMaskForUncalContinuum, unsigned numberOfPointsInUniformSample, float *uniform_wl, float *uniform_flux,float *uniform_Beamflux[MAXNUMBEROFBEAMS], bool useBeams);
    
    void trimOrdersByWavelengthRanges(int Minorder, int Maxorder);
    void applyTelluricRVCorrectionINTOExtendendSpectra(double tellCorrection, int Minorder, int Maxorder);
    void setRVCorrectionINTOExtendendSpectra(double rvelCorrection, int Minorder, int Maxorder);
    void correctFlatField(int Minorder, int Maxorder, bool StarPlusSky, bool starplusskyInvertSkyFiber=false, double skyOverStarFiberAreaRatio=1.0);
    
    /// Saves the flux into both the normalized flux and fcal flux of the extended spectra
    void saveFluxINTOExtendedSpectra(int Minorder, int Maxorder);
    
    /// Normalizes the flux order-by-order and saves the result into the normalized flux of the extended spectra
    void normalizeFluxOrderbyOrderINTOExtendendSpectra(unsigned normalizationBinsize, int Minorder, int Maxorder, bool normalizeBeams);
    
    /// Normalizes the flux to the continuum by sampling the continuum across all orders, and saves the result into the normalized flux of the extended spectra
    void normalizeFluxINTOExtendendSpectra(string inputWavelengthMaskForUncalContinuum, unsigned numberOfPointsInUniformSample, unsigned normalizationBinsize, double delta_wl, int Minorder, int Maxorder, bool normalizeBeams);
    
    /// Normalizes the flux to the continuum and applies flux calibrations (does some weird stuff currently) and saves the result into the normalized flux and the fcal flux of the extended spectra, respectively
    void normalizeAndCalibrateFluxINTOExtendendSpectra(string inputWavelengthMaskForUncalContinuum, double exposureTime, bool AbsoluteCalibration, unsigned numberOfPointsInUniformSample, unsigned normalizationBinsize, double delta_wl, int Minorder, int Maxorder, bool normalizeBeams, double SkyOverStarFiberAreaRatio, bool StarPlusSky);
    
    /// Applies a flat response calibration to the flux and saves the result into the fcal flux of the extended spectra
    void applyFlatResponseINTOExtendendSpectra(string flatResponse, bool FITSformat, int Minorder, int Maxorder);
    
    /// Applies the removal of the continuum polarization and saves the result into the continuum removed polarization
    void removeContinuumPolarization(int Minorder, int Maxorder);
    
    unsigned getMaxNumberOfElementsInOrder(int Minorder, int Maxorder) const;
    unsigned getNumberOfBeams(int Minorder, int Maxorder) const;
    
    operaSpectrum getSpectrumWithinTelluricMask(string inputWavelengthMaskForTelluric, int Minorder, int Maxorder, bool normalized, unsigned normalizationBinsize);
    operaSpectralLineList detectSpectralLinesWithinWavelengthMask(string inputWavelengthMaskForTelluric, int Minorder, int Maxorder, bool normalized, unsigned normalizationBinsize, double spectralResolution, bool emissionSpectrum,double LocalMaxFilterWidth,double MinPeakDepth,double DetectionThreshold,double nsigclip);

    operaSpectrum getSpectrumAroundLines(operaSpectrum sourceLines, int Minorder, int Maxorder, bool normalized, unsigned normalizationBinsize, double spectralResolution,double nsig, double snrClip, unsigned numberOfPointsToCutInOrderEnds);
    operaSpectrum getSpectrumWithinWavelengthRange(operaWavelengthRange range, int Minorder, int Maxorder, bool normalized, unsigned normalizationBinsize, double snrClip, unsigned numberOfPointsToCutInOrderEnds);

    void calculateRawFluxQuantities(int Minorder, int Maxorder, double *integratedFlux, double *meanFlux, double *maxSNR, double *meanSNR);

    void readFlatResponseIntoSED(string filename,int Minorder, int Maxorder, bool FITSformat);
    operaSpectrum readFITSFlatResponse(string filename);
	operaSpectrum readLibreEspritFlatResponse(string filename);
    
    /// Returns an operaSpectrum containing the specified parts of the extended spectrum orders based on input parameters.
    operaSpectrum getExtendedSpectrum(int Minorder, int Maxorder, unsigned RawNormalizedOrCalibrated, bool wavelengthTelluricCorrected, bool wavelengthHelioCorrected, double snrClip, unsigned numberOfPointsToCutInOrderEnds, double sourceRV_KMS);
};

#endif
