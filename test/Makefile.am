# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -g -O1 -L$(PREFIX)/src/libraries/  -L../lib -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib/ -L/usr/local/lib/ -I../include/ -I. -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -g -O1 -L$(PREFIX)/src/libraries/ -L../lib -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib/ -L/usr/local/lib/ -I../include/ -I. -I/usr/include/ -I/usr/local/include/
AM_LDFLAGS = -loperaImageVector -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectrumSimulation -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaFITSImage -loperaFITSSubImage  -loperaImageVector -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector  -loperaFluxVector -loperaGeometricShapes -loperaMatrix -lPixelSet -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -loperaMuellerMatrix -loperaFit -loperaLMFit -lPolynomial -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lz -lpthread -lm
#AM_LDFLAGS = -Wl,--no-as-needed
# This is for Linux...
LIBS = -loperaImageVector -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaSpectrumSimulation -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaFITSProduct -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaFITSImage -loperaFITSSubImage  -loperaImageVector -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector  -loperaFluxVector -loperaGeometricShapes -loperaMatrix -lPixelSet -loperaSpectralEnergyDistribution -loperaSpectrumSimulation -loperaMuellerMatrix -loperaFit -loperaLMFit -lPolynomial -loperaStats -loperaLib -lfftw3 -lgzstream -lcfitsio -lz -lpthread -lm
# this lists the binaries to produce
bin_PROGRAMS = operaAsmTest operaMatrixLibTest operaMathLibTest operaJDTest testmpfit operaFITSProductTest \
	operaMPFitLibTest operaFitLibTest operaImageOperatorTest operaFITSSubImageTest operaConfigurationAccesstest \
	operaFITSImageTest operaEspadonsImageTest operaStatsLibTest operaGeometricShapesTest operaExtractionApertureTest \
	operaFluxVectorTest operaPolarimetryTest operaCubeTest \
	operaPolarTest basicFITSImageTest gzstreamtest operaSextractorTest sitelletest SBIGtest FITSImageVectorTest \
	operaAOBImageTest operaNICIImageTest operaNIFSImageTest operaCreateInstrumentEnvironmentSetup nancheck \
	operastringstreamtest operaOESTest

#
# wcs support
#				
if LIBWCS
AM_LDFLAGS += -lwcs
LIBS += -lwcs
endif

operaOESTest_SOURCES = operaOESTest.cpp

operastringstreamtest_SOURCES = operastringstreamtest.cpp

nancheck_SOURCES = nancheck.cpp

operaCreateInstrumentEnvironmentSetup_SOURCES = operaCreateInstrumentEnvironmentSetup.cpp operaCreateInstrumentEnvironmentSetup.h

FITSImageVectorTest_SOURCES = FITSImageVectorTest.cpp

SBIGtest_SOURCES = SBIGtest.cpp

operaAOBImageTest_SOURCES = operaAOBImageTest.cpp

operaNICIImageTest_SOURCES = operaNICIImageTest.cpp

operaNIFSImageTest_SOURCES = operaNIFSImageTest.cpp

sitelletest_SOURCES = sitelletest.cpp

operaSextractorTest_SOURCES = operaSextractorTest.cpp

gzstreamtest_SOURCES = gzstreamtest.cpp

basicFITSImageTest_SOURCES = basicFITSImageTest.cpp

operaPolarTest_SOURCES = operaPolarTest.cpp

operaCubeTest_SOURCES = operaCubeTest.cpp

operaPolarimetryTest_SOURCES = operaPolarimetryTest.cpp

operaFluxVectorTest_SOURCES = operaFluxVectorTest.cpp

operaExtractionApertureTest_SOURCES = operaExtractionApertureTest.cpp

operaGeometricShapesTest_SOURCES = operaGeometricShapesTest.cpp

operaAsmTest_SOURCES = operaAsmTest.cpp
			
operaFITSImageTest_SOURCES = operaFITSImageTest.cpp operaFITSImageTest.h
			
operaEspadonsImageTest_SOURCES = operaEspadonsImageTest.cpp operaEspadonsImageTest.h
			
operaStatsLibTest_SOURCES = operaStatsLibTest.c operaStats.h

operaFitLibTest_SOURCES = operaFitLibTest.c operaFit.h

operaMPFitLibTest_SOURCES = operaMPFitLibTest.c

testmpfit_SOURCES = testmpfit.c

operaConfigurationAccesstest_SOURCES = operaConfigurationAccesstest.cpp operaConfigurationAccesstest.h
			
operaImageOperatorTest_SOURCES = operaImageOperatorTest.cpp operaImageOperatorTest.h
			
operaFITSSubImageTest_SOURCES = operaFITSSubImageTest.cpp operaFITSSubImageTest.h

operaFITSProductTest_SOURCES = operaFITSProductTest.cpp

operaJDTest_SOURCES = operaJDTest.c

operaMathLibTest_SOURCES = operaMathLibTest.c operaMath.h

operaMatrixLibTest_SOURCES = operaMatrixLibTest.c operaMatrix.h

