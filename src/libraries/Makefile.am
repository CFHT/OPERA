# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -g -O1 -L$(PREFIX)/src/libraries/  -L../../lib -L/usr/lib/ -L/usr/local/lib/ -I../../include/ -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -g -O1 -L$(PREFIX)/src/libraries/  -L../../lib -I../../include/ -I/usr/include/ -I/usr/local/include/ -L/usr/lib/ -L/usr/local/lib/
#AM_LDFLAGS = -loperaConfigurationAccess -loperaParameterAccess -loperaImage -lcfitsio -lm
# This is for Linux...
#LIBS = -loperaConfigurationAccess -loperaParameterAccess -loperaImage -lcfitsio -lm
LIBS = -lm

# this lists the libraries to produce, static and dynamic

lib_LIBRARIES = liboperaFFT.a 

lib_LTLIBRARIES = liboperaException.la libPolynomial.la liboperaVector.la liboperaFluxVector.la liboperaMath.la libGaussian.la liboperaWavelength.la \
	liboperaSpectralOrderVector.la liboperaSpectralElements.la liboperaSpectralOrder.la liboperaInstrumentProfile.la \
	liboperaGeometry.la liboperaCCD.la liboperaFITSSubImage.la liboperaEspadonsSubImage.la liboperaImage.la liboperaStats.la \
	liboperaLMFit.la liboperaFit.la liboperaConfigurationAccess.la liboperaParameterAccess.la liboperaFITSImage.la liboperaPolarimetry.la \
	liboperaEspadonsImage.la liboperaFITSProduct.la liboperaMEFFITSProduct.la liboperaMatrix.la liboperaSpectralLines.la \
	liboperaLib.la liboperaSpectralFeature.la liboperaGeometricShapes.la liboperaExtractionAperture.la libGainBiasNoise.la \
	liboperaMultiExtensionFITSImage.la liboperaMultiExtensionFITSCube.la liboperaFITSCube.la liboperaMuellerMatrix.la \
	liboperaImageVector.la liboperaStokesVector.la libPixelSet.la liboperaSpectralEnergyDistribution.la \
	libgzstream.la libLaurentPolynomial.la liboperaSpectralTools.la liboperaDateTime.la\
	liboperaObservingConditions.la liboperaObjectInTheSky.la liboperaSpectrograph.la liboperaTelescope.la liboperaInstrumentEnvironmentSetup.la \
	libArgumentHandler.la liboperaArgumentHandler.la liboperaCommonModuleElements.la liboperaIOFormats.la

liboperaInstrumentEnvironmentSetup_la_SOURCES = operaInstrumentEnvironmentSetup.cpp operaInstrumentEnvironmentSetup.h
liboperaInstrumentEnvironmentSetup_la_LDFLAGS = -version-info 1:0:0

liboperaTelescope_la_SOURCES = operaTelescope.cpp operaTelescope.h
liboperaTelescope_la_LDFLAGS = -version-info 1:0:0

liboperaSpectrograph_la_SOURCES = operaSpectrograph.cpp operaSpectrograph.h
liboperaSpectrograph_la_LDFLAGS = -version-info 1:0:0

liboperaObjectInTheSky_la_SOURCES = operaObjectInTheSky.cpp operaObjectInTheSky.h
liboperaObjectInTheSky_la_LDFLAGS = -version-info 1:0:0

liboperaObservingConditions_la_SOURCES = operaObservingConditions.cpp operaObservingConditions.h
liboperaObservingConditions_la_LDFLAGS = -version-info 1:0:0

liboperaDateTime_la_SOURCES = operaDateTime.cpp operaDateTime.h
liboperaDateTime_la_LDFLAGS = -version-info 1:0:0

liboperaSpectralTools_la_SOURCES = operaSpectralTools.cpp operaSpectralTools.h
liboperaSpectralTools_la_LDFLAGS = -version-info 1:0:0
liboperaSpectralTools_la_LIBADD = liboperaSpectralElements.la

libLaurentPolynomial_la_SOURCES = LaurentPolynomial.cpp LaurentPolynomial.h
libLaurentPolynomial_la_LDFLAGS = -version-info 1:0:0

libgzstream_la_SOURCES = gzstream.cpp gzstream.h
libgzstream_la_LDFLAGS = -version-info 1:0:0

liboperaSpectralEnergyDistribution_la_SOURCES = operaSpectralEnergyDistribution.cpp operaSpectralEnergyDistribution.h
liboperaSpectralEnergyDistribution_la_LDFLAGS = -version-info 1:0:0

libPixelSet_la_SOURCES = PixelSet.cpp PixelSet.h
libPixelSet_la_LDFLAGS = -version-info 1:0:0

liboperaImageVector_la_SOURCES = operaImageVector.cpp operaImageVector.h
liboperaImageVector_la_LDFLAGS = -version-info 1:0:0

liboperaMuellerMatrix_la_SOURCES = operaMuellerMatrix.cpp operaMuellerMatrix.h
liboperaMuellerMatrix_la_LDFLAGS = -version-info 1:0:0

liboperaStokesVector_la_SOURCES = operaStokesVector.cpp operaStokesVector.h
liboperaStokesVector_la_LDFLAGS = -version-info 1:0:0

liboperaFITSCube_la_SOURCES = operaFITSCube.cpp operaFITSCube.h
liboperaFITSCube_la_LDFLAGS = -version-info 1:0:0

liboperaMultiExtensionFITSCube_la_SOURCES = operaMultiExtensionFITSCube.cpp operaMultiExtensionFITSCube.h
liboperaMultiExtensionFITSCube_la_LDFLAGS = -version-info 1:0:0

liboperaMultiExtensionFITSImage_la_SOURCES = operaMultiExtensionFITSImage.cpp operaMultiExtensionFITSImage.h
liboperaMultiExtensionFITSImage_la_LDFLAGS = -version-info 1:0:0

liboperaPolarimetry_la_SOURCES = operaPolarimetry.cpp operaPolarimetry.h
liboperaPolarimetry_la_LDFLAGS = -version-info 1:0:0
liboperaPolarimetry_la_LIBADD = liboperaFluxVector.la

liboperaExtractionAperture_la_SOURCES = operaExtractionAperture.cpp operaExtractionAperture.h
liboperaExtractionAperture_la_LDFLAGS = -version-info 1:0:0

libGainBiasNoise_la_SOURCES = GainBiasNoise.cpp GainBiasNoise.h
libGainBiasNoise_la_LDFLAGS = -version-info 1:0:0

liboperaFluxVector_la_SOURCES = operaFluxVector.cpp operaFluxVector.h
liboperaFluxVector_la_LDFLAGS = -version-info 1:0:0

liboperaVector_la_SOURCES = operaVector.cpp operaVectorOperations.tpp operaVector.h operaVectorOperations.h
liboperaVector_la_LDFLAGS = -version-info 1:0:0

liboperaGeometricShapes_la_SOURCES = operaGeometricShapes.cpp operaGeometricShapes.h
liboperaGeometricShapes_la_LDFLAGS = -version-info 1:0:0

libGaussian_la_SOURCES = Gaussian.cpp Gaussian.h
libGaussian_la_LDFLAGS = -version-info 1:0:0

liboperaSpectralFeature_la_SOURCES = operaSpectralFeature.cpp operaSpectralFeature.h
liboperaSpectralFeature_la_LDFLAGS = -version-info 1:0:0

liboperaLib_la_SOURCES = operaLib.cpp operaLib.h
liboperaLib_la_LDFLAGS = -version-info 1:0:0

liboperaException_la_SOURCES = operaException.cpp operaException.h
liboperaException_la_LDFLAGS = -version-info 1:0:0

liboperaImage_la_SOURCES = operaImage.c operaImage.h operaLibCommon.h
liboperaImage_la_LDFLAGS = -version-info 1:0:0

liboperaMatrix_la_SOURCES = operaMatrix.cpp operaMatrix.h 
liboperaMatrix_la_LDFLAGS = -version-info 1:0:0

libPolynomial_la_SOURCES = Polynomial.cpp Polynomial.h
libPolynomial_la_LDFLAGS = -version-info 1:0:0

liboperaSpectralLines_la_SOURCES = operaSpectralLines.cpp operaSpectralLines.h
liboperaSpectralLines_la_LDFLAGS = -version-info 1:0:0

liboperaCCD_la_SOURCES = operaCCD.cpp operaCCD.h operaLibCommon.h
liboperaCCD_la_LDFLAGS = -version-info 1:0:0

liboperaMath_la_SOURCES = operaMath.c operaMath.h operaLibCommon.h
liboperaMath_la_LDFLAGS = -version-info 1:0:0

liboperaStats_la_SOURCES = operaStats.c operaStats.h operaLibCommon.h
liboperaStats_la_LDFLAGS = -version-info 1:0:0

liboperaLMFit_la_SOURCES = operaLMFit.c operaLMFit.h operaLibCommon.h
liboperaLMFit_la_LDFLAGS = -version-info 1:0:0

liboperaFit_la_SOURCES = operaFit.c operaFit.h ladfit.c mpfit.c mpfit.h operaLibCommon.h
liboperaFit_la_LDFLAGS = -version-info 1:0:0

liboperaFFT_a_SOURCES = operaFFT.c operaFFT.h operaLibCommon.h

liboperaFITSImage_la_SOURCES = operaFITSImage.cpp operaFITSImage.h operaLibCommon.h
liboperaFITSImage_la_LDFLAGS = -version-info 1:0:0

liboperaEspadonsImage_la_SOURCES = operaEspadonsImage.cpp operaEspadonsImage.h operaLibCommon.h
liboperaEspadonsImage_la_LDFLAGS = -version-info 1:0:0

liboperaFITSSubImage_la_SOURCES = operaFITSSubImage.cpp operaFITSSubImage.h operaLibCommon.h
liboperaFITSSubImage_la_LDFLAGS = -version-info 1:0:0

liboperaEspadonsSubImage_la_SOURCES = operaEspadonsSubImage.cpp operaEspadonsSubImage.h operaLibCommon.h
liboperaEspadonsSubImage_la_LDFLAGS = -version-info 1:0:0

liboperaFITSProduct_la_SOURCES = operaFITSProduct.cpp operaFITSProduct.h operaLibCommon.h
liboperaFITSProduct_la_LDFLAGS = -version-info 1:0:0

liboperaMEFFITSProduct_la_SOURCES = operaMEFFITSProduct.cpp operaMEFFITSProduct.h operaLibCommon.h
liboperaMEFFITSProduct_la_LDFLAGS = -version-info 1:0:0

liboperaConfigurationAccess_la_SOURCES = operaConfigurationAccess.c operaConfigurationAccess.h
liboperaConfigurationAccess_la_LDFLAGS = -version-info 1:0:0

liboperaParameterAccess_la_SOURCES = operaParameterAccess.c operaParameterAccess.h
liboperaParameterAccess_la_LDFLAGS = -version-info 1:0:0

liboperaSpectralElements_la_SOURCES = operaSpectralElements.cpp operaSpectralElements.h
liboperaSpectralElements_la_LDFLAGS = -version-info 1:0:0
liboperaSpectralElements_la_LIBADD = liboperaFluxVector.la liboperaWavelength.la

liboperaGeometry_la_SOURCES = operaGeometry.cpp operaGeometry.h
liboperaGeometry_la_LDFLAGS = -version-info 1:0:0

liboperaWavelength_la_SOURCES = operaWavelength.cpp operaWavelength.h
liboperaWavelength_la_LDFLAGS = -version-info 1:0:0
liboperaWavelength_la_LIBADD = libPolynomial.la libGaussian.la liboperaMath.la

liboperaSpectralOrder_la_SOURCES = operaSpectralOrder.cpp operaSpectralOrder.h
liboperaSpectralOrder_la_LDFLAGS = -version-info 1:0:0

liboperaInstrumentProfile_la_SOURCES = operaInstrumentProfile.cpp operaInstrumentProfile.h
liboperaInstrumentProfile_la_LDFLAGS = -version-info 1:0:0
liboperaInstrumentProfile_la_LIBADD = libPolynomial.la

liboperaSpectralOrderVector_la_SOURCES = operaSpectralOrderVector.cpp operaSpectralOrderVector.h
liboperaSpectralOrderVector_la_LDFLAGS = -version-info 1:0:0

libArgumentHandler_la_SOURCES = ArgumentHandler.cpp ArgumentHandler.h
libArgumentHandler_la_LDFLAGS = -version-info 1:0:0

liboperaArgumentHandler_la_SOURCES = operaArgumentHandler.cpp operaArgumentHandler.h
liboperaArgumentHandler_la_LDFLAGS = -version-info 1:0:0

liboperaCommonModuleElements_la_SOURCES = operaCommonModuleElements.cpp operaCommonModuleElements.h
liboperaCommonModuleElements_la_LDFLAGS = -version-info 1:0:0

liboperaIOFormats_la_SOURCES = operaIOFormats.cpp operaIOFormats.h
liboperaIOFormats_la_LDFLAGS = -version-info 1:0:0

#
# if we want png plotting support, bring in the freetype includes
# you may need to change the freetype-config path
#				
if LIBPNG
FT_ARG := `$(FREE_TYPE_CONFIG_PATH)freetype-config --cflags`
AM_CXXFLAGS += $(FT_ARG)
lib_LIBRARIES += liboperaPNG.a 
liboperaPNG_a_SOURCES = operaPNG.cpp operaPNG.h
endif

#
# wcs support
#				
if LIBWCS
lib_LIBRARIES += libwcs.a
libwcs_a_SOURCES =
endif


