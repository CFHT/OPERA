#ifndef OPERASPECTRALORDER_H
#define OPERASPECTRALORDER_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralOrder
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

#include "operaError.h"

#include "libraries/operaLibCommon.h"           // for operaSpectralOrder_t
#include "libraries/Polynomial.h"				// for Polynomial
#include "libraries/operaGeometry.h"			// for operaGeometry
#include "libraries/operaPolarimetry.h"			// for Polarimetry
#include "libraries/operaWavelength.h"			// for operaWavelength
#include "libraries/operaSpectralElements.h"	// for operaSpectralElements
#include "libraries/operaInstrumentProfile.h"	// for operaInstrumentProfile
#include "libraries/operaSpectralLines.h"       // for operaSpectralLines
#include "libraries/operaFITSImage.h"			// for Matrix
#include "libraries/operaExtractionAperture.h"  // for operaExtractionAperture
#include "libraries/operaSpectralEnergyDistribution.h" // for operaSpectralEnergyDistribution
#include "libraries/GainBiasNoise.h" // for GainBiasNoise
#include "libraries/operaSpectralTools.h"

using namespace std;

/*! 
 * \sa class operaSpectralOrder
 * \brief operaSpectralOrder
 * \details A spectral order (SO) consists of a data set containing 
 * \details the information concerned with a full spectral order.
 * \return none
 * \file operaSpectralOrder.h
 * \ingroup libraries
 */
class operaSpectralOrder {
	
private:
	
	unsigned orderNumber;
	
	operaSpectralOrder_t SpectrumType;				// what type of spectral order data is stored
	
	double SNR;										// SNR at wlc for this order
    
	doubleValue_t tiltInDegrees;					// tiltInDegrees    
	
	operaSpectralElements *SpectralElements;		// pointer to the operaSpectralElements
	
    operaSpectralElements *SkyElements;             // pointer to the sky elements        
    
	operaGeometry *Geometry;						// pointer to the operaGeometry class instance
	
	operaWavelength *Wavelength;					// pointer to the operaWavelength class instance
	
	operaInstrumentProfile *InstrumentProfile;		// pointer to the operaInstrumentProfile 
    
	operaSpectralLines *SpectralLines;              // pointer to the operaSpectralLines     
	
    operaPolarimetry *Polarimetry;					// pointer to the polarimetry    
    
	operaSpectralEnergyDistribution *SpectralEnergyDistribution;	// pointer to the operaSpectralEnergyDistribution class instance    
    
    unsigned numberOfBeams;
    
    operaSpectralElements *BeamElements[MAXNUMBEROFBEAMS]; // pointer to spectral elements for beams
    
    operaInstrumentProfile *BeamProfiles[MAXNUMBEROFBEAMS];    // pointer to instrument profiles for beams         
    
    operaSpectralElements *BackgroundElements[LEFTANDRIGHT]; // pointer to spectral elements for background 
    
	operaExtractionAperture<Line> *ExtractionApertures[MAXNUMBEROFBEAMS];    // pointer to the apertures for extraction
	
    operaExtractionAperture<Line> *BackgroundApertures[LEFTANDRIGHT];    // pointer to the apertures for background    
	
	operaSpectralEnergyDistribution *BeamSED[MAXNUMBEROFBEAMS];	// pointer to the operaSpectralEnergyDistribution class instance        
    
	bool hasSpectralElements;						// true if we have information of this type about this order
	bool hasSkyElements;    
	bool hasGeometry;
	bool hasWavelength;
	bool hasInstrumentProfile;
	bool hasSpectralLines;
	bool hasExtractionApertures;
	bool hasPolarimetry;    
	bool hasSNR;
	bool hasCenterSNROnly;
	bool hasSpectralEnergyDistribution;
	
	double minwavelength;
	double maxwavelength;
	bool hasWavelengthRange;
	double snrSpectralBinSize;
	
public:
	
	/*
	 * Constructors
	 */
	operaSpectralOrder();
	
	operaSpectralOrder(unsigned order) ;
	
	operaSpectralOrder(unsigned order, unsigned maxDatapoints, unsigned maxValues, unsigned maxElements, operaSpectralOrder_t format);
	
	/*
	 * Destructor
	 */
	~operaSpectralOrder();
	
	/*
	 * Common Methods
	 */
	
	void deleteAll();
	
	void deleteInstrumentProfile(void);	
	
	void deleteBeamProfiles(void);
    
	void deleteApertures(void);    
    
	void deleteBeamsAndBackgroundElements(void);    
	
	void createBeamsAndBackgrounds(unsigned nElements, unsigned nBeams, operaSpectralOrder_t format, bool extendedbeams = false);
    
    void deletePolarimetry(void);
    
    void createPolarimetry(unsigned nElements);
    
	unsigned getorder(void) const;
	
	void setorder(unsigned ordernumber);
	
	void sethasSpectralElements(bool HasSpectralElement);
	
	void sethasSkyElements(bool HasSkyElements);
	
	void sethasGeometry(bool HasGeometry);
	
	void sethasWavelength(bool HasWavelength);
	
	void sethasInstrumentProfile(bool hasInstrumentProfile);
	
	void sethasSpectralLines(bool HasSpectralLines);	
    
	void sethasExtractionApertures(bool HasSpectralLines); 
    
	void sethasPolarimetry(bool HasPolarimetry);
	
	void sethasSNR(bool HasSNR);
    
	void sethasCenterSNROnly(bool HasCenterSNROnly);
    
	void sethasSpectralEnergyDistribution(bool HasSpectralEnergyDistribution);     
	
	bool gethasSpectralElements(void) const;
	
	bool gethasSkyElements(void) const;
    
	bool gethasGeometry(void) const;
	
	bool gethasWavelength(void) const;
	
	bool gethasInstrumentProfile(void) const;
	
	bool gethasSpectralLines(void) const;
    
	bool gethasExtractionApertures(void) const;
	
	bool gethasPolarimetry(void) const;
    
	bool gethasSNR(void) const;
    
	bool gethasCenterSNROnly(void) const;
    
	bool gethasSpectralEnergyDistribution(void) const;
	
	void createGeometry(unsigned maxdatapoints, unsigned maxValues);
	
	void createWavelength(unsigned maxnumberofcoefficients);
	
	void createSpectralElements(unsigned maxdatapoints, operaSpectralOrder_t SpectrumType, bool extended = false);
	
    void createSkyElements(unsigned maxdatapoints, operaSpectralOrder_t SpectrumType);

	operaSpectralElements *getSpectralElements(void);
	const operaSpectralElements *getSpectralElements(void) const;
	
	operaSpectralElements *getSkyElements(void);
    const operaSpectralElements *getSkyElements(void) const;
    
	operaGeometry *getGeometry(void);
	const operaGeometry *getGeometry(void) const;
	
	operaWavelength *getWavelength(void);
	const operaWavelength *getWavelength(void) const;
	
	operaInstrumentProfile *getInstrumentProfile(void);
	const operaInstrumentProfile *getInstrumentProfile(void) const;
	
	void setSpectralLines(operaSpectralLines *spectralLines);
	
	operaSpectralLines *getSpectralLines(void);
	const operaSpectralLines *getSpectralLines(void) const;
	
    operaPolarimetry *getPolarimetry(void);
    const operaPolarimetry *getPolarimetry(void) const;
    
    operaSpectralEnergyDistribution *getSpectralEnergyDistribution(void);
	const operaSpectralEnergyDistribution *getSpectralEnergyDistribution(void) const;
	
    operaSpectralElements *getBeamElements(unsigned beam);
	const operaSpectralElements *getBeamElements(unsigned beam) const;
	
    operaInstrumentProfile *getBeamProfiles(unsigned beam);
    const operaInstrumentProfile *getBeamProfiles(unsigned beam) const;
    
    operaSpectralEnergyDistribution *getBeamSED(unsigned beam);
	const operaSpectralEnergyDistribution *getBeamSED(unsigned beam) const;
	
    operaSpectralElements *getBackgroundElements(unsigned LeftOrRight);
    const operaSpectralElements *getBackgroundElements(unsigned LeftOrRight) const;
    
	operaExtractionAperture<Line> *getExtractionApertures(unsigned beam);
	const operaExtractionAperture<Line> *getExtractionApertures(unsigned beam) const;
	
	operaExtractionAperture<Line> *calculateMainApertureFromExtractionBeams(bool useIP);
    
    operaExtractionAperture<Line> *getBackgroundApertures(unsigned LeftOrRight);
    const operaExtractionAperture<Line> *getBackgroundApertures(unsigned LeftOrRight) const;
    
    void setBeamElements(unsigned beam, operaSpectralElements *beamElements);
    
    void setBeamProfiles(unsigned beam, operaInstrumentProfile *beamProfiles);    
    
    void setBackgroundElements(unsigned LeftOrRight, operaSpectralElements *backgroundElements);
    
	void setExtractionApertures(unsigned beam, operaExtractionAperture<Line> *extractionApertures);
	
    void setBackgroundApertures(unsigned LeftOrRight, operaExtractionAperture<Line> *backgroundApertures);
	
	double getCenterSNR(void) const;
	
	void setCenterSNR(double Snr);
    
	unsigned getnumberOfBeams(void) const;
    
	void setnumberOfBeams(unsigned NumberOfBeams);    
    
	doubleValue_t getTiltInDegrees(void) const;
    
	void setTiltInDegrees(doubleValue_t TiltInDegrees); 
    
    void setTiltInDegrees(double tilt, double error);
    
    double getTiltInDegreesValue(void) const;
    
    double getTiltInDegreesError(void) const;
	
    double *getSNRVector(void);    
	
	void calculateSNR(void);
	
	float getCentralSmoothedSNR(int upperlowerbound) const;
	
	float getPeakSmoothedSNR(int upperlowerbound) const;

	void NormalizeFlat(operaFITSImage &flatMatrix, operaFITSImage &outputMatrix, unsigned nx, unsigned ny, unsigned binsize);
	
	void extractRawSum(operaFITSImage &inputImage, ofstream &sout);
	
	void extractRawSum(operaFITSImage &inputImage, float noise, float gain);
	
    void measureInstrumentProfileAlongRows(operaFITSImage &masterFlatImage, unsigned binsize, unsigned sampleElementForPlot, ostream *pout);
	
    void measureInstrumentProfileAlongRowsInto2DWithGaussian(operaFITSImage &masterFlatImage, operaFITSImage &badpix, unsigned binsize, float gaussSig, float tiltInDegrees, bool witherrors, unsigned sampleElementForPlot, ostream *pout, unsigned minimumLines);
	
	void CalculateWavelengthSolution(void);
	
	void setSpectralElementsByHeight(double Height);

    void setSpectralElementsByStitchingApertures(double effectiveApertureFraction);
    
	void setInstrumentProfileVector(unsigned IPxsize, unsigned IPxsampling, unsigned IPysize, unsigned IPysampling, unsigned NDataPoints);
	
    void setSpectralLines(operaFITSImage &masterCompImage, operaFITSImage &badpix, operaFITSImage &bias, float noise, float gain, float ReferenceLineWidth,float DetectionThreshold, float LocalMaxFilterWidth, float MinPeakDepth);        
	
    void measureInstrumentProfile(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, unsigned minimumLines);
	
    void measureInstrumentProfileWithBinning(operaFITSImage &masterCompImage, operaFITSImage &badpix, double binsize, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, unsigned minimumLines);
    
    void measureInstrumentProfileUsingMedian(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, unsigned minimumLines);
    
    void measureInstrumentProfileUsingWeightedMean(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout);
    
    void recenterOrderPosition(void);
    
    void setApertureElements(operaSpectralOrder_t format);
    
    operaFluxVector extractFluxElement(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, unsigned indexElem, const PixelSet *aperturePixels);
    operaFluxVector extractSubpixelFlux(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, unsigned badpixValue, GainBiasNoise &gainBiasNoise, unsigned indexElem, const PixelSet *aperturePixels, Vector<unsigned>* pixcol=0, Vector<unsigned> *pixrow=0);
    
    void extractSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, const operaFluxVector& backgroundModelFlux);
    
    operaFluxVector extractBackground(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, unsigned NumberofElementsToBin);
    
    void normalizeFluxToAperture();
    
    void extractRawSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction);
	
    void extractRawSpectrumRejectingBadpixels(operaFITSImage &objectImage, operaFITSImage &flatImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize, double minSigmaClip, unsigned iterations, bool onTargetProfile, bool usePolynomialFit, bool removeBackground, bool verbose, bool calculateXCorrelation, ostream *pout);
    
    void extractStandardSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize);
    
    void extractStandardSpectrumNoBackground(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction);
	
    void extractOptimalSpectrum(operaFITSImage &objectImage, operaFITSImage &flatImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize, double minSigmaClip, double sigmaClipRange, unsigned iterations, bool onTargetProfile, bool usePolynomialFit, bool removeBackground, bool verbose, bool calculateXCorrelation, ostream *pout);
	
    void measureBeamSpatialProfiles(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, bool usePolynomialFit);
    
    void updateBadPixelsToRejectCosmicRays(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double minSigmaClip);
    
    void measureRawSpectrumRejectingBadpixels(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise);
    
    void measureOptimalSpectrum(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double minSigmaClip, double sigmaClipRange);
    
    void calculateXCorrBetweenIPandImage(operaFITSImage &Image, operaFITSImage &badpix, ostream *pout);
	
	operaSpectralOrder_t getSpectrumType(void) const;
	
	void setSpectrumType(operaSpectralOrder_t format);
	
    void printBeamSpectrum(ostream *pout);
	
	void printBeamSpectrum(string addFirstColumnEntry, ostream *pout);
	
	double getminwavelength() const;
	double getmaxwavelength() const;
	bool gethasWavelengthRange() const;
	void setminwavelength(double wl);
	void setmaxwavelength(double wl);
	void sethasWavelengthRange(bool hasRange);
	double getsnrSpectralBinSize() const;
	void setsnrSpectralBinSize(double spectralbinsize);
	
	void setWavelengthsFromCalibration();

	/*
	 * Normalization/Flux Calibration...
	 */
	operaSpectralElements& MainAndBeamElements(unsigned index);

	operaSpectralEnergyDistribution& MainAndBeamSED(unsigned index);

	unsigned MainAndBeamCount();

	operaVector getMainAndBeamFluxes(unsigned index);

	void fitSEDUncalibratedFluxToSample(const operaSpectrum& uniformSample);
	
	void fitSEDFcalAndThroughputToFlatResp(const operaSpectrum& flatresp);
	
	void applyNormalization(unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial, bool overwriteUncalFlux, bool normalizeBeams);
	
    void applyNormalizationForEmissionSpectrum(unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial, bool overwriteUncalFlux, bool normalizeBeams);

    void applyNormalizationFromExistingContinuum(bool normalizeBeams);

    void normalizeSpectrum(const operaFluxVector &uncalibratedFlux, operaFluxVector &normalizedFlux, operaFluxVector &outputContinuum, unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial);
	
    void measureContinuum(const operaFluxVector &uncalibratedFlux,operaFluxVector &outputContinuum, unsigned binsize, unsigned nsigcut, unsigned orderOfPolynomial, bool usePolynomial);
	
    void deleteSpectralEnergyDistribution(void);
    
    void createSpectralEnergyDistribution();
    
    void calculateContinuum(unsigned binsize, unsigned nsigcut);
    
    void normalizeSpectralElementsByConstant(const operaVector& maxFluxForNormalization);

	void divideSpectralElementsBySEDElements(bool useThroughput);

    void multiplySpectralElementsBySEDElements(bool useThroughput, const operaVector& spectralBinConstants);
    
    /*
     * Star+Sky Mode
     */
    void calculateStarAndSkyElements(bool starplusskyInvertSkyFiber, double skyOverStarFiberAreaRatio);
    
	/*
	 * Radial Velocity Wavelength Correction
	 */
    void applyRvelWavelengthCorrection(double RVcorrectionInKmPerSecond);
	void setExtendedRvelWavelengthCorrection(double RVcorrectionInKmPerSecond);
    void applyWavelengthCorrectionFromExtendedRvel(void);
    
    void CreateExtendedVectors();
	void CopyFluxVectorIntoRawFlux();
	void CopyRawFluxIntoFluxVector();
	void CopyFluxVectorIntoFcalFlux();
	void CopyFcalFluxIntoFluxVector();
	void CopyFluxVectorIntoNormalizedFlux();
	void CopyNormalizedFluxIntoFluxVector();
	
    void TrimOrderToWavelengthRange();

	operaSpectralLineList getRawLinesFromUncalSpectrum(double linewidth, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, double nsigclip, dispersionaxis_t dispersiontype);
	operaSpectralLineList detectSpectralLines(double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, double nsigclip, double spectralResolution, bool emissionSpectrum);
};

#endif

