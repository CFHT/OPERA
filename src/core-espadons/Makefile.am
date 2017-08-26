#
# espadons build makefile Generator
#
# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/ -I../../include/ -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -I../../include/ -I/usr/include/ -I/usr/local/include/ -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/
AM_LDFLAGS = -loperaCommonModuleElements -loperaIOFormats -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaDateTime -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -loperaVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaImageVector -loperaGeometricShapes -loperaMatrix -loperaSpectralEnergyDistribution -lPixelSet -loperaStats -loperaLib -loperaArgumentHandler -lArgumentHandler -lfftw3 -lgzstream -lcfitsio -lz -lpthread -lm
# This is for Linux...
LIBS = -loperaCommonModuleElements -loperaIOFormats -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaDateTime -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -loperaVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaStats -loperaLib -loperaArgumentHandler -lArgumentHandler -lfftw3 -lgzstream -lcfitsio -lz -lpthread -lm

# this lists the binaries to produce
bin_PROGRAMS = operaSNR operaWavelengthCalibration \
			operaMasterFlat operaPolar operaBias operaExtraction \
			operaGain operaReductionSet operaGeometryCalibration \
			operaMasterComparison operaMasterBias operaMasterFabPerot \
			operaInstrumentProfileCalibration operaPixelSensitivityMap \
			operaCreateProduct operaExtractionApertureCalibration \
			operaStarPlusSky operaStarOnly operaPolarIntensity \
			operaCreateFluxCalibration operaTelluricWavelengthCorrection \
			operaOrderSpacingCalibration operaPolarimetryCorrection \
			operaCalculateSpectralResolution operaMasterFluxCalibration \
			operaCreateFlatFieldFluxCalibration \
			operaGenerateLEFormats operaCreateFlatResponse \
			operaHeliocentricWavelengthCorrection

#
# wcs support
#				
if LIBWCS
AM_LDFLAGS += -lwcs
LIBS += -lwcs
endif

operaCreateFlatResponse_SOURCES = operaCreateFlatResponse.cpp

operaGenerateLEFormats_SOURCES = operaGenerateLEFormats.cpp

operaCreateFlatFieldFluxCalibration_SOURCES = operaCreateFlatFieldFluxCalibration.cpp

operaMasterFluxCalibration_SOURCES = operaMasterFluxCalibration.cpp

operaCalculateSpectralResolution_SOURCES = operaCalculateSpectralResolution.cpp

operaPolarimetryCorrection_SOURCES = operaPolarimetryCorrection.cpp operaExtendedSpectrumCreation.cpp

operaOrderSpacingCalibration_SOURCES = operaOrderSpacingCalibration.cpp

operaCreateFluxCalibration_SOURCES = operaCreateFluxCalibration.cpp

operaStarPlusSky_SOURCES = operaStarPlusSky.cpp operaExtendedSpectrumCreation.cpp

operaPolarIntensity_SOURCES = operaPolarIntensity.cpp operaExtendedSpectrumCreation.cpp

operaStarOnly_SOURCES = operaStarOnly.cpp operaExtendedSpectrumCreation.cpp

operaExtraction_SOURCES = operaExtraction.cpp

operaCreateProduct_SOURCES = operaCreateProduct.cpp

operaGeometryCalibration_SOURCES = operaGeometryCalibration.cpp

operaSNR_SOURCES = operaSNR.cpp

operaBias_SOURCES = operaBias.cpp

operaMasterFlat_SOURCES = operaMasterFlat.cpp operaMasterCalibration.cpp
			
operaMasterComparison_SOURCES = operaMasterComparison.cpp
			
operaMasterBias_SOURCES = operaMasterBias.cpp operaMasterCalibration.cpp
			
operaMasterFabPerot_SOURCES = operaMasterFabPerot.cpp operaMasterCalibration.cpp
			
operaPolar_SOURCES = operaPolar.cpp

operaGain_SOURCES = operaGain.cpp
			
operaReductionSet_SOURCES = operaReductionSet.cpp
			
operaWavelengthCalibration_SOURCES = operaWavelengthCalibration.cpp
			
operaInstrumentProfileCalibration_SOURCES = operaInstrumentProfileCalibration.cpp

operaExtractionApertureCalibration_SOURCES = operaExtractionApertureCalibration.cpp

operaPixelSensitivityMap_SOURCES = operaPixelSensitivityMap.cpp

operaTelluricWavelengthCorrection_SOURCES = operaTelluricWavelengthCorrection.cpp

operaHeliocentricWavelengthCorrection_SOURCES = operaHeliocentricWavelengthCorrection.cpp
