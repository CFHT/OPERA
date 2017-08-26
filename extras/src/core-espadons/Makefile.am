#
# espadons build makefile Generator
#
# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/ -I../../include/ -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -I../../include/ -I/usr/include/ -I/usr/local/include/ -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/
AM_LDFLAGS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaImageVector -loperaGeometricShapes -loperaMatrix -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -lPixelSet -loperaStats -loperaLib -lArgumentHandler -loperaArgumentHandler -loperaCommonModuleElements -lfftw3 -lgzstream -lcfitsio -lz -lsofa_c -lpthread -lm
# This is for Linux...
LIBS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaJD -loperaHelio -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -loperaStats -loperaLib -lArgumentHandler -loperaArgumentHandler -loperaCommonModuleElements -lfftw3 -lgzstream -lcfitsio -lsofa_c -lz -lpthread -lm

# this lists the binaries to produce
bin_PROGRAMS = operaNormalize operaFluxCalibration operaEchelleDispersionCalibration operaNormalizeAcrossOrders operaExtractSpectralLines \
			operaObjectInTheSkySetup operaObservingConditionsSetup operaSpectrographSetup operaTelescopeSetup

#
# wcs support
#				
if LIBWCS
AM_LDFLAGS += -lwcs
LIBS += -lwcs
endif

operaTelescopeSetup_SOURCES = operaTelescopeSetup.cpp operaTelescopeSetup.h

operaSpectrographSetup_SOURCES = operaSpectrographSetup.cpp operaSpectrographSetup.h

operaObservingConditionsSetup_SOURCES = operaObservingConditionsSetup.cpp operaObservingConditionsSetup.h

operaObjectInTheSkySetup_SOURCES = operaObjectInTheSkySetup.cpp operaObjectInTheSkySetup.h

operaExtractSpectralLines_SOURCES = operaExtractSpectralLines.cpp operaExtractSpectralLines.h

operaNormalizeAcrossOrders_SOURCES = operaNormalizeAcrossOrders.cpp operaNormalizeAcrossOrders.h

operaEchelleDispersionCalibration_SOURCES = operaEchelleDispersionCalibration.cpp operaEchelleDispersionCalibration.h

operaFluxCalibration_SOURCES = operaFluxCalibration.cpp operaFluxCalibration.h

operaNormalize_SOURCES = operaNormalize.cpp operaNormalize.h
