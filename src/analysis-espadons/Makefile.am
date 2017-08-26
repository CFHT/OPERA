#
# espadons build makefile Generator for analysis espadons
#
# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -L/usr/local/lib/ -L/usr/lib/ -I../../include/ -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -ggdb3 -O1 -L$(PREFIX)/src/libraries/ -L../../lib -I../../include/ -I/usr/local/include/ -L/usr/local/lib/ -L/usr/lib/
AM_LDFLAGS = -loperaCommonModuleElements -loperaIOFormats -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaDateTime -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -loperaVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaImageVector -loperaGeometricShapes -loperaMatrix -loperaSpectralEnergyDistribution -lPixelSet -loperaStats -loperaLib -loperaArgumentHandler -lArgumentHandler -lfftw3 -lgzstream -lcfitsio -lz -lpthread -lm
# This is for Linux...
LIBS = -loperaCommonModuleElements  -loperaIOFormats -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaDateTime -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -loperaVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaStats -loperaLib -loperaArgumentHandler -lArgumentHandler -lfftw3 -lgzstream -lcfitsio -lz -lpthread -lm

# this lists the binaries to produce
bin_PROGRAMS = operaBinPolarData operaBinFluxData operaRadialVelocity operaStackObjectSpectra operaRadialVelocityFromSelectedLines

operaBinPolarData_SOURCES = operaBinPolarData.cpp

operaBinFluxData_SOURCES = operaBinFluxData.cpp operaBinFluxData.h

operaRadialVelocity_SOURCES = operaRadialVelocity.cpp operaRadialVelocity.h

operaRadialVelocityFromSelectedLines_SOURCES = operaRadialVelocityFromSelectedLines.cpp operaRadialVelocity.h

operaStackObjectSpectra_SOURCES = operaStackObjectSpectra.cpp

