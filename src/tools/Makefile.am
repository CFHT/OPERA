# what flags you want to pass to the C compiler & linker
AM_CFLAGS = --pedantic -Wall -std=c99 -g -O1 -L$(PREFIX)/src/libraries/ -L../../lib -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/ -I../../include/ -I/usr/local/include/
AM_CXXFLAGS = --pedantic -Wall -g -O1 -L$(PREFIX)/src/libraries/ -L../../lib -I../../include/ -I/usr/include/ -I/usr/local/include/ -L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -L/usr/local/lib/
AM_LDFLAGS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaDateTime -loperaIOFormats -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaFluxVector -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -lPolynomial -lLaurentPolynomial -loperaMath -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -loperaFITSSubImage -loperaEspadonsSubImage -loperaPolarimetry -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaException -lGainBiasNoise -loperaMuellerMatrix -loperaStokesVector -loperaVector -loperaFFT -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPixelSet -loperaSpectralEnergyDistribution -loperaStats -loperaLib -loperaArgumentHandler -lArgumentHandler -lfftw3 -lgzstream -lcfitsio -lz -lpthread -lm
# This is for Linux...
LIBS = -loperaInstrumentEnvironmentSetup -loperaObservingConditions -loperaObjectInTheSky -loperaSpectrograph -loperaTelescope -loperaDateTime -loperaIOFormats -loperaSpectralOrderVector -loperaSpectralOrder -loperaSpectralElements -loperaSpectralLines -loperaSpectralFeature -loperaExtractionAperture -loperaSpectralTools -loperaFluxVector -loperaVector -lGaussian -loperaInstrumentProfile -loperaGeometry -loperaWavelength -loperaEspadonsSubImage -loperaFITSProduct -loperaMEFFITSProduct -loperaFITSImage -loperaEspadonsImage -loperaMultiExtensionFITSCube -loperaMultiExtensionFITSImage -loperaFITSCube -loperaPolarimetry -loperaException -lGainBiasNoise -loperaFFT -loperaMuellerMatrix -loperaStokesVector -loperaFluxVector -loperaImageVector -loperaFITSSubImage -loperaImageVector -loperaGeometricShapes -loperaMatrix -lPolynomial -lLaurentPolynomial -loperaMath -loperaCCD -loperaFit -loperaImage -loperaLMFit -loperaConfigurationAccess -loperaParameterAccess -lPixelSet -loperaSpectralEnergyDistribution -loperaStats -loperaLib -loperaArgumentHandler -lArgumentHandler -lfftw3 -lgzstream -lcfitsio -lz -lpthread -lm

# this lists the binaries to produce
bin_PROGRAMS =  operaConfigurationAccess operaParameterAccess \
				opera operaImage operaFITS2txt operaCreateBadpixMask \
				operagetheader operasaturated operaConvert2ampTo1amp operaExtractRawSum \
				operaimarith operaQueryImageInfo operasetheader operacompress \
				operabiasinjector operaPlotInstrumentProfile \
				operaStatistics espqlh catz operaFITSDisplayImage operaimagestats \
				operads9thumbs operaRotate \
				operaEspadonsETC operaExtractImage operaPlotOut \
				operaRotateMirrorCrop operaMedianCombine operaMJD
#
# if we want png plotting support, bring in the png libs and freetype
#				
if LIBPNG
FT_ARG := `$(FREE_TYPE_CONFIG_PATH)freetype-config --cflags` `$(FREE_TYPE_CONFIG_PATH)freetype-config --libs`
AM_CXXFLAGS += $(FT_ARG) 
AM_CFLAGS += $(FT_ARG) 
bin_PROGRAMS += operaFITStoPNG operaPlotGeom operaPlotWave operaPlot
AM_LDFLAGS += -loperaPNG -lpng -lz -lfreetype
LIBS +=  -loperaPNG -lpng -lz -lfreetype
operaFITStoPNG_SOURCES = operaFITStoPNG.c operaFITStoPNG.h
operaPlotGeom_SOURCES = operaPlotGeom.cpp operaPlotGeom.h
operaPlotWave_SOURCES = operaPlotWave.cpp operaPlotWave.h
operaPlot_SOURCES = operaPlot.cpp operaPlot.h
endif

#
# wcs support
#				
if LIBWCS
AM_LDFLAGS += -lwcs
LIBS += -lwcs
endif

operaMJD_SOURCES = operaMJD.cpp

operaRotateMirrorCrop_SOURCES = operaRotateMirrorCrop.cpp

operaMedianCombine_SOURCES = operaMedianCombine.cpp

operaPlotOut_SOURCES = operaPlotOut.cpp

operaExtractImage_SOURCES = operaExtractImage.cpp

operaEspadonsETC_SOURCES = operaEspadonsETC.cpp operaEspadonsETC.h

operaRotate_SOURCES = operaRotate.cpp

operads9thumbs_SOURCES = operads9thumbs.cpp

operaimagestats_SOURCES = operaimagestats.cpp

operaFITSDisplayImage_SOURCES = operaFITSDisplayImage.cpp

espqlh_SOURCES = espqlh.cpp

catz_SOURCES = catz.cpp

operaStatistics_SOURCES = operaStatistics.cpp

operaPlotInstrumentProfile_SOURCES = operaPlotInstrumentProfile.cpp

operabiasinjector_SOURCES = operabiasinjector.cpp

operaExtractRawSum_SOURCES = operaExtractRawSum.cpp operaExtractRawSum.h

operaConfigurationAccess_SOURCES = operaConfigurationAccess.cpp operaConfigurationAccess.h

operaParameterAccess_SOURCES = operaParameterAccess.cpp operaParameterAccess.h

opera_SOURCES = opera.cpp opera.h

operaImage_SOURCES = operaImage.cpp operaImage.h

operagetheader_SOURCES = operagetheader.cpp operagetheader.h

operasetheader_SOURCES = operasetheader.cpp operasetheader.h

operacompress_SOURCES = operacompress.cpp operacompress.h

operaimarith_SOURCES = operaimarith.cpp operaimarith.h

operasaturated_SOURCES = operasaturated.cpp operasaturated.h

operaCreateBadpixMask_SOURCES = operaCreateBadpixMask.cpp operaCreateBadpixMask.h

operaFITS2txt_SOURCES = operaFITS2txt.cpp operaFITS2txt.h

operaConvert2ampTo1amp_SOURCES = operaConvert2ampTo1amp.cpp operaConvert2ampTo1amp.h

operaQueryImageInfo_SOURCES = operaQueryImageInfo.cpp operaQueryImageInfo.h
