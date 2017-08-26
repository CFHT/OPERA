#
# WIRCam build makefile Generator
#
# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/ -I../../include/ -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -I../../include/ -I/usr/include/ -I/usr/local/include/ -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/
AM_LDFLAGS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaEspadonsImage -loperaWIRCamImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaFITSImage -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaImageVector -loperaGeometricShapes -loperaMatrix -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -lPixelSet -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lz -lsofa_c -lpthread -lm
# This is for Linux...
LIBS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaWIRCamImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaEspadonsImage -loperaFITSImage -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lsofa_c -lz -lpthread -lm

# this lists the binaries to produce
bin_PROGRAMS = wirDetrend wirAstrometry wirPhotometry wirMasterDark wirMasterTwilightFlat wirTwilightFlat wirPickSkies wirSubtractSky wirPickTwilightFlats wirCreateZeroPoints wirSkyFlat

wirDetrend_SOURCES = wirDetrend.cpp wirDetrend.h

wirAstrometry_SOURCES = wirAstrometry.cpp wirAstrometry.h

wirPhotometry_SOURCES = wirPhotometry.cpp wirPhotometry.h

wirMasterDark_SOURCES = wirMasterDark.cpp wirMasterDark.h

wirSkyFlat_SOURCES = wirSkyFlat.cpp wirSkyFlat.h

wirTwilightFlat_SOURCES = wirTwilightFlat.cpp wirTwilightFlat.h

wirMasterTwilightFlat_SOURCES = wirMasterTwilightFlat.cpp wirMasterTwilightFlat.h

wirPickSkies_SOURCES = wirPickSkies.cpp wirPickSkies.h

wirSubtractSky_SOURCES = wirSubtractSky.cpp wirSubtractSky.h

wirPickTwilightFlats_SOURCES = wirPickTwilightFlats.cpp wirPickTwilightFlats.h

wirCreateZeroPoints_SOURCES = wirCreateZeroPoints.cpp wirCreateZeroPoints.h
