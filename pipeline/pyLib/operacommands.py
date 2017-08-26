# -*- coding: utf-8 -*-
"""
    Created on Feb 20 2017
    Description: Library of OPERA commands.
    @author: Eder Martioli
    Laboratorio Nacional de Astrofisica, Brazil
    """

#### Function to generate a command line for BIAS list: ####
def BiasListCommand(Dirs, Instmode, Readmode, keywords,listfilename) :
    commandline = Dirs.EXE + 'operaQueryImageInfo -r ' + Dirs.DATADIR + \
    ' -q "'+keywords.INSTRUMEKEY+' '+keywords.READMODEKEY+' '+keywords.OBSTYPEKEY+'"' + \
    ' '+keywords.INSTRUMEKEY+'="'+ Instmode.INSTRUME +'" '+keywords.READMODEKEY+'="'+Readmode.READOUTSPEED+'" '+keywords.OBSTYPEKEY+'='+keywords.BIASKEYWORD + \
    ' > ' + listfilename
    return commandline
###########################################

#### Function to test BIAS list: ####
def testBiasListCommand(Dirs, Instmode, Readmode, keywords) :
    commandline = Dirs.EXE + 'operaQueryImageInfo -r ' + Dirs.DATADIR + \
    ' -q "'+keywords.INSTRUMEKEY+' '+keywords.READMODEKEY+' '+keywords.OBSTYPEKEY+'"' + \
    ' '+keywords.INSTRUMEKEY+'="'+ Instmode.INSTRUME +'" '+keywords.READMODEKEY+'="'+Readmode.READOUTSPEED+'" '+keywords.OBSTYPEKEY+'='+keywords.BIASKEYWORD
    return commandline
###########################################

#### Function to generate a command line for FLAT, COMP lists: ####
def CalibrationListCommand(Dirs, Instmode, Readmode, obstypekeyvalue, keywords, listfilename, allowanyreadout) :
    commandline = ""
    if allowanyreadout :
        commandline = Dirs.EXE + 'operaQueryImageInfo -r ' + Dirs.DATADIR + \
        ' -q "'+keywords.INSTRUMEKEY+' '+keywords.INSTMODEKEY+' '+keywords.OBSTYPEKEY+'"' + \
        ' '+keywords.INSTRUMEKEY+'="'+ Instmode.INSTRUME +'" '+keywords.INSTMODEKEY+'="'+Instmode.INSTRUMENTMODEKEY+'" '+keywords.OBSTYPEKEY+'='+obstypekeyvalue + \
        ' > ' + listfilename
    else :
        commandline = Dirs.EXE + 'operaQueryImageInfo -r ' + Dirs.DATADIR + \
        ' -q "'+ keywords.INSTRUMEKEY+' '+keywords.INSTMODEKEY+' '+keywords.READMODEKEY+' '+keywords.OBSTYPEKEY +'"' + \
        ' '+ keywords.INSTRUMEKEY+'="'+ Instmode.INSTRUME +'" '+keywords.INSTMODEKEY+'="'+Instmode.INSTRUMENTMODEKEY+'" '+keywords.READMODEKEY+'="'+Readmode.READOUTSPEED+'" '+keywords.OBSTYPEKEY +'='+obstypekeyvalue + \
        ' > ' + listfilename
    return commandline
###########################################

#### Function to test FLAT, COMP lists: ####
def testCalibrationListCommand(Dirs, Instmode, Readmode, obstypekeyvalue, keywords, allowanyreadout) :
    commandline = ""
    if allowanyreadout :
        commandline = Dirs.EXE + 'operaQueryImageInfo -r ' + Dirs.DATADIR + \
        ' -q "'+keywords.INSTRUMEKEY+' '+keywords.INSTMODEKEY+' '+keywords.OBSTYPEKEY+'"' + \
        ' '+keywords.INSTRUMEKEY+'="'+ Instmode.INSTRUME +'" '+keywords.INSTMODEKEY+'="'+Instmode.INSTRUMENTMODEKEY+'" '+keywords.OBSTYPEKEY+'='+obstypekeyvalue
    else :
        commandline = Dirs.EXE + 'operaQueryImageInfo -r ' + Dirs.DATADIR + \
        ' -q "'+ keywords.INSTRUMEKEY+' '+keywords.INSTMODEKEY+' '+keywords.READMODEKEY+' '+keywords.OBSTYPEKEY +'"' + \
        ' '+ keywords.INSTRUMEKEY+'="'+ Instmode.INSTRUME +'" '+keywords.INSTMODEKEY+'="'+Instmode.INSTRUMENTMODEKEY+'" '+keywords.READMODEKEY+'="'+Readmode.READOUTSPEED+'" '+keywords.OBSTYPEKEY +'='+obstypekeyvalue
    return commandline
###########################################

#### Function to generate a command line for OBJECT list on screen: ####
def ObjectListCommand(Dirs, Instmode, Readmode, obstypekeyvalue, keywords) :
    commandline = Dirs.EXE + 'operaQueryImageInfo -r ' + Dirs.DATADIR + \
    ' -q "'+ keywords.INSTRUMEKEY+' '+keywords.INSTMODEKEY+' '+keywords.READMODEKEY+' '+keywords.OBSTYPEKEY +'"' + \
    ' '+ keywords.INSTRUMEKEY+'="'+ Instmode.INSTRUME +'" '+keywords.INSTMODEKEY+'="'+Instmode.INSTRUMENTMODEKEY+'" '+keywords.READMODEKEY+'="'+Readmode.READOUTSPEED+'" '+keywords.OBSTYPEKEY +'='+obstypekeyvalue
    
    return commandline
###########################################

#### Function to generate a command line to create OBJECT list into a file: ####
def ObjectListCommandToFile(Dirs, Instmode, Readmode, obstypekeyvalue, keywords, listfilename) :
    commandline = ObjectListCommand(Dirs, Instmode, Readmode, obstypekeyvalue, keywords) + \
    ' > ' + listfilename
    return commandline
###########################################

#### Function to generate a command line for mastercalibrations: ####
def MasterCalibrationCommand(Dirs, bin, output, list) :
    
    inputlistpar = ''
    if bin == 'operaMedianCombine' :
        inputlistpar += ' --list=' + list
    else :
        inputlistpar += ' --imagelistfile=' + list
    
    commandline = Dirs.EXE + bin + ' --output=' + output + inputlistpar
    return commandline
###########################################

#### Function to generate a command line for master flat response: ####
def MasterFlatResponseCommand(Dirs, bin, output, pythonlist) :
    list = ''
    for file in pythonlist :
        list += ' --images=' + file
    
    commandline = Dirs.EXE + bin + ' --output=' + output + list
    return commandline
###########################################

#### Function to generate a command line for mastercomparison: ####
def MasterComparisonCommand(Dirs, product, list, badpix, masterbias, Instmode) :
    
    commandline = Dirs.EXE + 'operaMasterComparison --output=' + product + \
    ' --imagelistfile=' + list + ' --badpixelmask=' + badpix + ' --masterbias=' + masterbias + \
    ' --combineMethod=' + Instmode.MASTERCOMP_COMBINEMETHOD + ' --saturationLimit='+Instmode.MASTERCOMP_SATURATIONLIMIT + \
    ' --outputExposureTime=' + Instmode.MASTERCOMP_OUTPUTEXPOSURETIME + \
    ' --truncateOuputFluxToSaturation=' + Instmode.MASTERCOMP_TRUNCATEOUTPUTFLUXTOSATURATION + \
    ' --biasConstant='+Instmode.MASTERCOMP_BIASCONSTANT + \
    ' --expTimeFITSKeyword='+Instmode.MASTERCOMP_EXPTIMEFITSKEYWORD

    return commandline
###########################################

#### Function to generate a command line for operaGain: ####
def GainCommand(Dirs, output, biaslist, flatlist, badpix, Readmode, Instmode) :

    datasec_str = ''
    if Readmode.GAIN_NUMBEROFAMPLIFIERS == '1' :
        datasec_str += ' --DATASEC="' + Instmode.GAIN_DATASEC + '"'
    elif Readmode.GAIN_NUMBEROFAMPLIFIERS == '2' :
        datasec_str += ' --DSECA="'+Instmode.GAIN_DSECA + '"' + ' --DSECB="'+Instmode.GAIN_DSECB + '"'

    commandline = Dirs.EXE + 'operaGain --output=' + output + \
    ' --biaslistfile=' + biaslist + ' --flatlistfile=' + flatlist + \
    ' --badpixelmask=' + badpix + ' --defaultgain=' +  Readmode.DEFAULTGAIN + \
    ' --defaultnoise=' + Readmode.DEFAULTNOISE + datasec_str + \
    ' --numberofamplifiers='+Readmode.GAIN_NUMBEROFAMPLIFIERS+' --subwindow="' + Instmode.GAIN_SUBWINDOW + '"' + \
    ' --gainMinPixPerBin=' + Instmode.GAIN_MINPIXPERBIN + ' --gainMaxNBins=' + Instmode.GAIN_MAXNBINS + \
    ' --gainLowestCount=' + Instmode.GAIN_LOWESTCOUNT + ' --gainHighestCount=' + Instmode.GAIN_HIGHESTCOUNT
    
    return commandline
###########################################

###########################################

#### Function to generate a command line for operaOrderSpacingCalibration: ####
def OrderSpacingCommand(Dirs, output, gainproduct, masterbias, masterflat, badpix, Instmode, plots) :
    
    plotstring = ' --plotfilename=' + plots["ORDSPCPLOTFILE"] + ' --datafilename=' + plots["ORDSPCDATAFILE"] + ' --scriptfilename=' + plots["ORDSPCSCRIPTFILE"]
    
    commandline = Dirs.EXE + 'operaOrderSpacingCalibration --orderspacingoutput=' + output + \
    ' --inputGainFile=' + gainproduct + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + \
    ' --badpixelmask=' + badpix + ' --aperture=' + Instmode.SPACING_APERTURE + \
    ' --referenceOrderNumber=' + Instmode.SPACING_REFERENCEORDERNUMBER + \
    ' --referenceOrderSeparation=' + Instmode.SPACING_REFERENCEORDERSEPARATION + \
    ' --numberOfsamples=' + Instmode.SPACING_NUMBEROFSAMPLES + ' --sampleCenterPosition='+ Instmode.SPACING_SAMPECENTERPOSITION + \
    ' --subformat="' + Instmode.SPACING_SUBFORMAT + '" --detectionMethod=' + Instmode.SPACING_DETECTIONMETHOD + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for operaGeometryCalibration: ####
def GeometryCommand(Dirs, output, gain, masterbias, masterflat, badpix, orderspacing, Instmode, plots) :
    
    plotstring = ' --plotfilename=' + plots["GEOMPLOTFILE"] + ' --datafilename=' + plots["GEOMDATAFILE"] + ' --scriptfilename=' + plots["GEOMSCRIPTFILE"]
    
    commandline = Dirs.EXE + 'operaGeometryCalibration --outputGeomFile=' + output + \
    ' --inputGainFile=' + gain + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + \
    ' --badpixelmask=' + badpix + ' --inputOrderSpacing=' + orderspacing + ' --aperture=' + Instmode.GEOM_APERTURE + \
    ' --maxorders=' + Instmode.GEOM_MAXNORDERS + ' --minordertouse=' + Instmode.GEOM_MINORDERTOUSE + \
    ' --recenterIPUsingSliceSymmetry=' + Instmode.GEOM_RECENTERIPUSINGSLICESYMMETRY + \
    ' --totalNumberOfSlices=' + Instmode.GEOM_NUMBEROFSLICES + ' --subformat="' + Instmode.GEOM_SUBFORMAT + '"' + \
    ' --detectionMethod=' + Instmode.GEOM_DETECTIONMETHOD + ' --FFTfilter=' + Instmode.GEOM_FFTFILTER + ' --nsamples=' + Instmode.GEOM_NSAMPLES + \
    '  --orderOfTracingPolynomial=' + Instmode.GEOM_ORDEROFTRACINGPOLYNOMIAL + ' --binsize=' + Instmode.GEOM_BINSIZE + \
    ' --colDispersion=' + Instmode.GEOM_COLDISPERSION + ' --invertOrders=' +  Instmode.GEOM_INVERTORDERS + \
    ' --referenceOrderSamplePosition=' + Instmode.GEOM_REFERENCEORDERSAMPLEPOSITION + ' --graces=' + Instmode.GEOM_GRACES + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for operaInstrumentProfileCalibration: ####
def InstrumentProfileCommand(Dirs, output, geom, gain, masterbias, masterflat, mastercomp, badpix, Instmode, plots) :
    
    plotstring = ' --plotfilename=' + plots["PROFPLOTFILE"] + ' --datafilename=' + plots["PROFDATAFILE"] + ' --scriptfilename=' + plots["PROFSCRIPTFILE"]
    
    commandline = Dirs.EXE + "operaInstrumentProfileCalibration --outputProf=" + output + \
    ' --geometryfilename=' + geom + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + \
    ' --badpixelmask=' + badpix + ' --mastercomparison=' +  mastercomp + ' --gainfilename=' + gain + \
    ' --xSize=' + Instmode.PROF_IPXSIZE + ' --ySize=' + Instmode.PROF_IPYSIZE + ' --xSampling=' + Instmode.PROF_IPXSAMPLING + \
    ' --ySampling=' + Instmode.PROF_IPYSAMPLING + ' --referenceLineWidth=' + Instmode.PROF_REFERENCELINEWIDTH + \
    ' --binsize=' + Instmode.PROF_BINSIZE + ' --method=' + Instmode.PROF_IPMETHOD + ' --tilt=' + Instmode.PROF_TILTANGLE + \
    ' --spectralElementHeight=' + Instmode.PROF_SPECTRALELEMENTHEIGHT + ' --maxthreads=' + Instmode.PROF_MAXTHREADS + \
    ' --minimumlines=' + Instmode.PROF_MINIMUMLINES + ' --LocalMaxFilterWidth=' + Instmode.PROF_LOCALMAXFILTERWIDTH + \
    ' --DetectionThreshold=' + Instmode.PROF_DETECTIONTHRESHOLD + ' --MinPeakDepth=' + Instmode.PROF_MINPEAKDEPTH + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for operaExtractionApertureCalibration: ####
def ApertureCommand(Dirs, output, geom, prof, orderspacing, Instmode, plots) :
    
    plotstring = ' --plotfilename=' + plots["APERPLOTFILE"] + \
        ' --datafilename=' + plots["APERDATAFILE"] + \
            ' --scriptfilename=' + plots["APERSCRIPTFILE"] + \
                ' --tiltplotfilename=' + plots["APERTILTPLOTFILE"] + \
                    ' --tiltdata1filename=' + plots["APERTILTDATA1FILE"] + \
                        ' --tiltdata2filename=' + plots["APERTILTDATA2FILE"] + \
                            ' --tiltscriptfilename=' + plots["APERTILTSCRIPTFILE"]

    commandline = Dirs.EXE + "operaExtractionApertureCalibration --outputApertureFile=" + output + \
    ' --inputgeom=' + geom + ' --inputprof=' + prof + ' --inputorderspacing=' + orderspacing + \
    ' --numberOfBeams=' + Instmode.APER_NUMBEROFBEAMS + ' --gapBetweenBeams=' + Instmode.APER_GAP + \
    ' --apertureHeight=' + Instmode.APER_APERTUREHEIGHT + ' --apertureWidth=' + Instmode.APER_APERTURE + \
    ' --constantTilt=' + Instmode.APER_CONSTANTTILTFLAG + \
    ' --backgroundAperture=' + Instmode.APER_SKYAPERAPERTURE + ' --pickImageRow=' + Instmode.APER_PICKIMAGEROW + \
    ' --nRowSamples=' + Instmode.APER_NROWSAMPLES + ' --xbin=' + Instmode.APER_XBIN + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT + \
    plotstring
    
    return commandline
###########################################

#### Function to generate a command line for Raw Extraction of Comparison spectra: ####
def compRawExtractionCommand(Dirs, output, inputImage, masterbias, masterflat, badpix, gainproduct, geomproduct, profproduct, aperproduct, Instmode) :

    commandline = Dirs.EXE + 'operaExtraction --outputSpectraFile=' + output + \
    ' --inputImage=' + inputImage + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + ' --badpixelmask=' + badpix + \
    ' --inputGainFile=' + gainproduct + ' --inputGeometryFile=' + geomproduct + \
    ' --inputInstrumentProfileFile=' + profproduct + ' --inputApertureFile=' + aperproduct + \
    ' --spectrumtype=5 --spectrumtypename=RawBeamSpectrum --maxthreads=' + Instmode.EXTRACTION_MAXTHREADS + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT
    
    return commandline
###########################################

#### Function to generate a command line for Optimal Extraction of spectra: ####
def calibrationExtractionCommand(Dirs, output, inputImage, masterbias, masterflat, badpix, gainproduct, geomproduct, profproduct, aperproduct, Instmode) :
    
    commandline = Dirs.EXE + 'operaExtraction --outputSpectraFile=' + output + \
    ' --inputImage=' + inputImage + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + ' --badpixelmask=' + badpix + \
    ' --inputGainFile=' + gainproduct + ' --inputGeometryFile=' + geomproduct + \
    ' --inputInstrumentProfileFile=' + profproduct + ' --inputApertureFile=' + aperproduct + \
    ' --spectrumtype=' + Instmode.CALEXTRACTION_SPECTRUMTYPE + ' --spectrumtypename=' + Instmode.CALEXTRACTION_SPECTRUMTYPENAME + \
    ' --backgroundBinsize='+Instmode.CALEXTRACTION_BACKGROUNDBINSIZE + \
    ' --minsigmaclip='+Instmode.CALEXTRACTION_MINSIGMACLIP + ' --sigmacliprange='+Instmode.CALEXTRACTION_SIGMACLIPRANGE + \
    '  --removeBackground='+ Instmode.CALEXTRACTION_REMOVEBACKGROUND + ' --iterations=' + Instmode.CALEXTRACTION_ITERATIONS + \
    ' --rejectBadpixInRawExtraction=' + Instmode.CALEXTRACTION_REJECTBADPIXINRAWEXTRACTION + ' --onTargetProfile='+Instmode.CALEXTRACTION_ONTARGETPROFILE + \
    ' --usePolynomialFit='+Instmode.CALEXTRACTION_USEPOLYNOMIALFIT+' --maxthreads=' + Instmode.EXTRACTION_MAXTHREADS + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT
    
    return commandline
###########################################

#### Function to generate a command line for operaWavelengthCalibration: ####
def WavelengthCommand(Dirs, output, outputresolution, geomproduct, compspectrum, Instmode, config, lineSetFilename, plots) :
    
    plotstring = ' --ordersplotfilename=' + plots["WAVEORDSPLOTFILE"] + ' --specplotfilename=' + plots["WAVESPECPLOTFILE"] + \
    ' --ordersscriptfilename=' + plots["WAVEORDSCRIPTFILE"] + ' --specscriptfilename=' + plots["WAVESPECSCRIPTFILE"] + \
    ' --ordersdatafilename=' + plots["WAVEORDSDATAFILE"] + ' --atlasdatafilename=' + plots["WAVEATLASDATAFILE"] + \
    ' --compdatafilename=' + plots["WAVECOMPDATAFILE"] + ' --linesdatafilename=' + plots["WAVELINESDATAFILE"]
    
    commandline = Dirs.EXE + 'operaWavelengthCalibration --outputWaveFile=' + output + \
    ' --outputResolutionFile=' + outputresolution + \
    ' --inputGeomFile=' + geomproduct + ' --uncalibrated_spectrum=' + compspectrum +  \
    ' --uncalibrated_linewidth=' + Instmode.WAVE_UNCALLINEWIDTH + \
    ' --atlas_lines=' + config.THARATLASLINES + ' --atlas_spectrum=' + config.THARATLASSPECTRUM + \
    ' --parseSolution='+Instmode.WAVE_PARSESOLUTION+' --ParRangeSizeInPerCent='+Instmode.WAVE_PARRANGESIZEINPERCENT+\
    ' --NpointsPerPar='+Instmode.WAVE_NPOINTSPERPAR+' --maxNIter='+Instmode.WAVE_MAXNITER + \
    ' --minNumberOfLines='+Instmode.WAVE_MINNUMBEROFLINES+' --maxorderofpolynomial='+Instmode.WAVE_MAXORDEROFPOLYNOMIAL + \
    ' --dampingFactor='+Instmode.WAVE_DAMPINGFACTOR+' --initialAcceptableMismatch='+Instmode.WAVE_INITIALACCEPTABLEMISMATCH + \
    ' --nsigclip='+Instmode.WAVE_NSIGCLIP+' --normalizeUncalibratedSpectrum='+Instmode.WAVE_NORMALIZEUNCALIBRATEDSPECTRUM + \
    ' --normalizationBinSize='+Instmode.WAVE_NORMALIZATIONBINSIZE+' --LocalMaxFilterWidth='+Instmode.WAVE_LOCALMAXFILTERWIDTH + \
    ' --DetectionThreshold=' + Instmode.WAVE_DETECTIONTHRESHOLD + ' --MinPeakDepth=' + Instmode.WAVE_MINPEAKDEPTH + \
    ' --referenceOrder=' + Instmode.WAVE_REFERENCEORDER + ' --nOrdersToSearchAround=' + Instmode.WAVE_NORDERSTOSEARCHAROUND + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT + \
    lineSetFilename + plotstring
    
    return commandline
###########################################

#### Function to generate a command line for operaCreateFlatFieldFluxCalibration: ####
def FlatFluxCalibrationCommand(Dirs, output, flatspectrum, Instmode, wave) :
    
    commandline = Dirs.EXE + 'operaCreateFlatFieldFluxCalibration --outputFluxCalibrationFile=' + output + \
    ' --inputMasterFlatSpectrum=' + flatspectrum + ' --wavelengthCalibration=' + wave + \
    ' --wavelengthForNormalization=' + Instmode.FCAL_WAVELENGTHFORNORMALIZATION + ' --binsize='+Instmode.FLATFLUXCAL_BINSIZE+ \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT
    
    return commandline
###########################################

#### Function to generate a command line to create OBJECT list into a file: ####
def RotateMirrorCropCommand(Dirs, Instmode, file) :
    
    commandline = Dirs.EXE + 'operaRotateMirrorCrop --images=' + file + \
    ' --outputdir=' + Dirs.PRODUCTDIR + ' --cropsubwindow="' + Instmode.ROTMIRRCROP_CROPSUBWINDOW + '"' + \
    ' --rotate=' + Instmode.ROTMIRRCROP_ROTATE + ' --mirrorcols=' + Instmode.ROTMIRRCROP_MIRRORCOLS + \
    ' --mirrorrows=' + Instmode.ROTMIRRCROP_MIRRORROWS + ' --sufix=' + Instmode.ROTMIRRCROP_SUFFIX
    
    return commandline
###########################################

#### Function to generate a command line for Extraction of object spectra: ####
def objectExtractionCommand(Dirs, product, inputImage, masterbias, masterflat, badpix, gainproduct, geomproduct, profproduct, aperproduct, Instmode) :

    commandline = Dirs.EXE + 'operaExtraction --outputSpectraFile=' + product + \
    ' --inputImage=' + inputImage + ' --masterbias=' + masterbias + ' --masterflat=' + masterflat + ' --badpixelmask=' + badpix + \
    ' --inputGainFile=' + gainproduct + ' --inputGeometryFile=' + geomproduct + \
    ' --inputInstrumentProfileFile=' + profproduct + ' --inputApertureFile=' + aperproduct + \
    ' --starplusskymode=' + str(Instmode.STARPLUSKYMODEFLAG) + ' ' + Instmode.EXTRACTION_INVERTSKYFIBERFLAG + \
    ' --spectrumtype=' + Instmode.EXTRACTION_SPECTRUMTYPE + ' --spectrumtypename='+Instmode.EXTRACTION_SPECTRUMTYPENAME+\
    ' --backgroundBinsize='+Instmode.EXTRACTION_BACKGROUNDBINSIZE + \
    ' --minsigmaclip='+Instmode.EXTRACTION_MINSIGMACLIP + ' --sigmacliprange='+Instmode.EXTRACTION_SIGMACLIPRANGE + \
    ' --removeBackground='+ Instmode.EXTRACTION_REMOVEBACKGROUND + ' --iterations=' + Instmode.EXTRACTION_ITERATIONS + \
    ' --rejectBadpixInRawExtraction=' + Instmode.EXTRACTION_REJECTBADPIXINRAWEXTRACTION + ' --onTargetProfile='+Instmode.EXTRACTION_ONTARGETPROFILE + \
    ' --usePolynomialFit='+Instmode.EXTRACTION_USEPOLYNOMIALFIT+' --maxthreads=' + Instmode.EXTRACTION_MAXTHREADS + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT
    
    return commandline
###########################################

#### Function to generate a command line for Telluric Wavelength Correction: ####
def TelluricWaveCommand(Dirs, product, inputSpectrum, wave, flatSpectrum, config, Instmode, plotbool) :
    
    if plotbool :
        plotstring = ' --rvcorrsplotfilename=' + "rvcorr.eps" + ' --specplotfilename=' + "tellspec.eps" + \
    ' --rvcorrscriptfilename=' + "rvcorr.gnu" + ' --specscriptfilename=' + "spec_tmp.gnu" + \
    ' --rvcorrdatafilename=' + "rvcorr.dat" + ' --rvcorrfitdatafilename=' + "rvcorr-fit.dat" + \
    ' --specdatafilename=' + "tellspec.dat"
    else :
        plotstring = ''
    
    if Instmode.STARPLUSKYMODEFLAG != 0 :
        flagstring = ' --StarPlusSky'
    else :
        flagstring = ''

    commandline = Dirs.EXE + 'operaTelluricWavelengthCorrection --outputWaveFile=' + product + \
    ' --inputObjectSpectrum=' + inputSpectrum + ' --inputWaveFile=' + wave + \
    ' --telluric_lines=' + config.TELLURICLINES + \
    ' --inputWavelengthMaskForTelluric=' + config.TELLURICWAVELENGTHMASK + \
    ' --spectralResolution=' + Instmode.SPECTRALRESOLUTION + \
    ' --radialVelocityRange=' + Instmode.TELL_RADIALVELOCITYRANGE + \
    ' --radialVelocityStep=' + Instmode.TELL_RADIALVELOCITYSTEP + \
    ' --XCorrelationThreshold='+Instmode.TELL_XCORRELATIONTHRESHOLD+' --normalizationBinsize='+Instmode.TELL_NORMALIZATIONBINSIZE + \
    ' --inputFlatFluxCalibration=' + flatSpectrum + ' --useFitToFindMaximum' + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT + \
    ' --RVCorrectionMethod=' + Instmode.TELL_RVCORRECTIONMETHOD + \
    ' --LocalMaxFilterWidth=' + Instmode.TELL_LOCALMAXFILTERWIDTH + \
    ' --MinPeakDepth=' + Instmode.TELL_MINPEAKDEPTH + \
    ' --DetectionThreshold=' + Instmode.TELL_DETECTIONTHRESHOLD + \
    ' --nsigclip=' + Instmode.TELL_NSIGCLIP + \
    ' --minNumberOfMatchedLines=' + Instmode.TELL_MINNUMBEROFMATCHEDLINES + \
    ' --duplicateLineThreshold=' + Instmode.TELL_DUPLICATELINETHRESHOLD + \
    Instmode.TELL_INVERTSKYFIBERFLAG + flagstring + plotstring

    return commandline
##########################################

#### Function to generate a command line for Radial Velocity module: ####
def RadialVelocityCommand(Dirs, output, outputdata, inputSpectrum, rvel, config, Instmode, mjdate) :
    
    commandline = Dirs.EXE + 'operaRadialVelocity --outputRVFile=' + output + \
    ' --outputIndividualMeasurements=' + outputdata + \
    ' --inputObjectSpectrum=' + inputSpectrum + \
    ' --telluric_lines=' + config.TELLURICLINES + \
    ' --template_spectrum=' + config.SOMEREFERENCESPECTRUM + \
    ' --inputWavelengthRangesForRVMeasurements=' + config.WAVELENGTHRANGESFORRVMEASUREMENTS + \
    ' --inputHeliocentricCorrection=' + rvel + \
    ' --mjdate=' + str(mjdate) + \
    ' --useFitToFindMaximum' + \
    ' --radialVelocityRange=' + Instmode.RV_RADIALVELOCITYRANGE + \
    ' --radialVelocityStep=' + Instmode.RV_RADIALVELOCITYSTEP + \
    ' --threshold=' + Instmode.RV_THRESHOLD + \
    ' --spectralResolution=' + Instmode.SPECTRALRESOLUTION + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT
    
    return commandline
##########################################

#### Function to generate a command line for Radial Velocity module: ####
def RadialVelocity2Command(Dirs, output, inputSpectrum, wave, rvel, tell, flatSpectrum, config, Instmode, mjdate, plotbool) :
    
    if Instmode.STARPLUSKYMODEFLAG != 0 :
        flagstring = ' --StarPlusSky'
    else :
        flagstring = ''
    
    commandline = Dirs.EXE + 'operaRadialVelocityFromSelectedLines --outputRVFile=' + output + \
    ' --inputObjectSpectrum=' + inputSpectrum + ' --inputWaveFile=' + wave + \
    ' --telluric_lines=' + config.TELLURICLINES + \
    ' --source_lines=' + config.SOURCELINES + \
    ' --mjdate=' + str(mjdate) + \
    ' --inputWavelengthMaskForTelluric=' + config.TELLURICWAVELENGTHMASK + \
    ' --spectralResolution=' + Instmode.SPECTRALRESOLUTION + \
    ' --normalizationBinsize=' + Instmode.RV2_NORMALIZATIONBINSIZE + \
    ' --initialRVguess=' + Instmode.RV2_INITIALRVGUESS +\
    ' --inputHeliocentricCorrection=' + rvel + \
    ' --inputTelluricCorrection=' + tell + \
    ' --inputFlatFluxCalibration=' + flatSpectrum  + \
    ' --robustFit' + \
    ' --sourceLineWidthResolution=' + Instmode.RV2_SOURCELINERESOLUTION + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT + \
    ' --gnuScriptFileName=' + \
    flagstring

    return commandline
##########################################

#### Function to generate a command line for Barycentric Wavelength Correction: ####
def HeliocentricWaveCommand(Dirs, output, ra, dec, mjdate, exptime, startHA, config, Instmode) :

    commandline = Dirs.EXE + 'operaHeliocentricWavelengthCorrection --outputRVelFile=' + output + \
    ' --observatory_coords="' + Instmode.observatory_coords + '" --observatory_elevation=' + Instmode.observatory_elevation + \
    ' --object_coords="' + str(ra) + ' ' + str(dec) + '" --MJDTime=' + str(mjdate) + ' --etime=' + str(exptime) + \
    ' --leapseconds=' + config.LEAPSECONDSFILE
    
    return commandline
##########################################

#### Function to generate a command line to create a flux calibration spectrum *.fcal: ####
def CreateFcalCommand(Dirs, output, inputSpectrum, stdcaldata, flatSpectrum, config, Instmode, aperture, wave, exptime):
    
    commandline = Dirs.EXE + 'operaCreateFluxCalibration --outputFluxCalibrationFile=' + output + \
    ' --inputUncalibratedSpectrum=' + inputSpectrum + ' --inputCalibratedSpectrum=' + stdcaldata  + \
    ' --inputFlatFluxCalibration=' + flatSpectrum + ' --inputWavelengthMaskForRefContinuum=' + config.ATLASWAVELENGTHMASK + \
    ' --inputWavelengthMaskForUncalContinuum=' + config.ATYPEWAVELENGTHMASK + ' --inputWaveFile=' + wave + \
    ' --inputApertureFile=' + aperture  + ' --wavelengthForNormalization=' + Instmode.FCAL_WAVELENGTHFORNORMALIZATION + \
    ' --exposureTime=' + str(exptime) + ' --numberOfPointsInUniformSample=' + Instmode.FCAL_NUMBEROFPOINTSINUNIFORMSAMPLE + \
    ' --numberOfPointsInUniformRefSample=' + Instmode.FCAL_NUMBEROFPOINTSINUNIFORMREFSAMPLE + ' --binsize=' + Instmode.FCAL_BINSIZE + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT
    
    return commandline
##########################################

#### Function to generate a command line to create a LE flat response file *.s: ####
def CreateFlatResponseCommand(Dirs, output, inputSpectrum, inputImage, stdcaldata, flatSpectrum, config, Instmode, aperture, wave):

    commandline = Dirs.EXE + 'operaCreateFlatResponse --outputFlatResponseFile=' + output + \
    ' --inputUncalibratedSpectrum=' + inputSpectrum + ' --inputSpectrumFITSImage=' + inputImage + \
    ' --outputFITS' + ' --inputCalibratedSpectrum=' + stdcaldata  + \
    ' --inputFlatFluxCalibration=' + flatSpectrum + ' --inputWavelengthMaskForRefContinuum=' + config.ATLASWAVELENGTHMASK + \
    ' --inputWavelengthMaskForUncalContinuum=' + config.SOLARTYPEWAVELENGTHMASK + ' --inputWaveFile=' + wave + \
    ' --wavelengthForNormalization=' + Instmode.FCAL_WAVELENGTHFORNORMALIZATION + \
    ' --numberOfPointsInUniformSample=' + Instmode.FLATRESP_NUMBEROFPOINTSINUNIFORMSAMPLE + \
    ' --numberOfPointsInUniformRefSample=' + Instmode.FLATRESP_NUMBEROFPOINTSINUNIFORMREFSAMPLE + ' --binsize=' + Instmode.FLATRESP_BINSIZE + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT
    
    return commandline
##########################################

#### Function to generate a command line for calibrated spectrum *.spc: ####
def SpcModuleCommand(Dirs, output, Instmode, config, inputSpectrum, flatSpectrum, inputFcal, rvelwave, tellwave, wave, objectname, exptime) :
    
    commandline = Dirs.EXE + Instmode.SPC_MODULENAME + ' --outputCalibratedSpectrum=' + output + \
    ' --inputUncalibratedSpectrum=' + inputSpectrum + ' --inputFlatFluxCalibration=' + flatSpectrum  + \
    ' --fluxCalibration=' + inputFcal + \
    ' --radialvelocitycorrection=' + rvelwave + ' --telluriccorrection=' + tellwave + ' --wavelengthCalibration=' + wave + \
    ' --object="' + objectname + '" --etime=' + str(exptime) + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT + \
    ' --spectrumtype=' + Instmode.SPC_SPECTRUMTYPE + ' --numberOfPointsInUniformSample=' + Instmode.SPC_NUMBEROFPOINTSINUNIFORMSAMPLE + \
    ' --normalizationBinsize=' + Instmode.SPC_NORMALIZATIONBINSIZE + ' --AbsoluteCalibration=' + Instmode.SPC_ABSOLUTECALIBRATION + \
    Instmode.SPC_SKYOVERSTARFIBERAREARATIO + Instmode.SPC_INVERTSKYFIBERFLAG
    
    return commandline
##########################################

#### Function to generate a command line to stack *.spc spectra: ####
def StackObjectSpectra(Dirs, output, objectname, Instmode, inputSpectralFiles, inputRVFiles, combineMethod, spectrumTypeToExtract) :
    
    inputfileentries = ' --inputspectrum="'
    for file in inputSpectralFiles :
        inputfileentries += file + ' '
    inputfileentries += '"'

    inputrvfileentries = ' --inputRVfile="'
    for file in inputRVFiles :
        inputrvfileentries += file + ' '
    inputrvfileentries += '"'
    
    commandline = Dirs.EXE + 'operaStackObjectSpectra --outputspectrum=' + output + \
    inputfileentries + inputrvfileentries +  \
    ' --combineMethod=' + str(combineMethod) + \
    ' --spectrumTypeToExtract=' + str(spectrumTypeToExtract) + \
    ' --applyTelluricWaveCorrection=' + Instmode.STACKSPC_APPLYTELL + \
    ' --applyHeliocentricRVCorrection=' + Instmode.STACKSPC_APPLYHELIO + \
    ' --object=' + objectname + \
    ' --RadialVelocityBin=' + Instmode.STACKSPC_RVBIN + \
    ' --firstWavelength=' + Instmode.STACKSPC_MINWL + \
    ' --lastWavelength=' + Instmode.STACKSPC_MAXWL + \
    ' --snrClip=' + Instmode.STACKSPC_SNRCLIP + \
    ' --numberOfPointsToCutInOrderEnds=' + Instmode.STACKSPC_NCUTATORDERENDS
    
    return commandline
##########################################

#### Function to generate SNR file: ####
def GenerateSNR(Dirs, input, output, wave, object, Instmode) :
    commandline = Dirs.EXE + 'operaSNR' + \
    ' --input=' + input + \
    ' --output=' + output + \
    ' --wavelengthCalibration=' + wave + \
    ' --object=' + object + \
    ' --spectrumtype=' + Instmode.SNR_SPECTRUMTYPE + \
    ' --centralsnr=' + Instmode.SNR_CENTRALSNR + \
    ' --spectralbinsize=' + Instmode.SNR_SPECTRALBINSIZE
    return commandline
##########################################

#### Function to generate a command line to export *.spc to fits: ####
def exportToFits(Dirs, output, spectrumfile, originalfitsfile, wavePixfile) :
    
    commandline = Dirs.PYTHONMODULES + 'exportToFits.py' + \
    ' --outputfitsfilename=' + output + \
    ' --spectrumfile=' + spectrumfile + \
    ' --originalfitsfile=' + originalfitsfile + \
    ' --wavePixfile=' + wavePixfile
    
    return commandline
##########################################

#### Function to generate a command line for converting from *.spc to Libre-Esprit format: ####
def GenLEFormatsCommand(Dirs, config, product, inputSpectrum, Instmode, objectname, fluxtype, wavetype) :
    
    commandline = Dirs.EXE + 'operaGenerateLEFormats --outputLEfilename=' + product + \
    ' --inputOperaSpectrum=' + inputSpectrum + ' --LibreEspritSpectrumType=' + Instmode.LESPECTRUMTYPE + \
    ' --object="' + objectname + '" --fluxType=' + str(fluxtype) + ' --wavelengthType=' + str(wavetype) + \
    ' --LEorderwavelength=' + config.LEORDERWAVELENGTH
    
    return commandline
##########################################

#### Function to generate a command line for raw polarimetry *.p.gz: ####
def PolarCommand(Dirs, Instmode, product, INPUT1, INPUT2, INPUT3, INPUT4, wave, stokes) :
    
    commandline = Dirs.EXE + 'operaPolar --output=' + product + \
    ' --input1=' + INPUT1 + ' --input2=' + INPUT2 + ' --input3=' + INPUT3 + ' --input4=' + INPUT4 + \
    ' --inputWaveFile=' + wave + ' --stokesparameter=' + str(stokes) + \
    ' --numberofexposures=4 --method=2 --ordernumber=-999 ' + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT
    
    return commandline
##########################################

#### Function to generate a command line for polarimetry calibrated product *.pol.gz : ####
def CalibratedPolarCommand(Dirs, product, Instmode, config, polar, flatSpectrum, inputFcal, rvelwave, tellwave, wave, objectname, exptime) :
    
    commandline = Dirs.EXE + 'operaPolarimetryCorrection --outputCalibratedSpectrum=' + product + \
    ' --polar=' + polar + ' --inputFlatFluxCalibration=' + flatSpectrum + \
    ' --fluxCalibration=' + inputFcal + ' --flatResponse=' + config.OLAPAFLATRESPONSE + \
    ' --radialvelocitycorrection=' + rvelwave + ' --telluriccorrection=' + tellwave + ' --wavelengthCalibration=' + wave + \
    ' --object="' + objectname + '" --etime=' + str(exptime) + ' --inputWavelengthMaskForUncalContinuum=' + config.ATYPEWAVELENGTHMASK + \
    ' --spectrumtype=18 --numberOfPointsInUniformSample=150 --normalizationBinsize=750 --AbsoluteCalibration=0' + \
    ' --minorder=' + Instmode.MINORDERTOEXTRACT + ' --maxorder=' + Instmode.MAXORDERTOEXTRACT + \
    ' --wlrangefile=' + config.LEORDERWAVELENGTH
    
    return commandline
##########################################

#### Function to generate the final product in FITS format out of LE files: ####
def operaCreateProduct_LE_Intensity(Dirs, objectfile, output, inputUS, inputUN, inputUW, inputNW, rvelwave, tellwave, snrfile, sresfile, operaversion, date) :
    
    commandline = Dirs.EXE + 'operaCreateProduct' + \
    ' --input=' + objectfile + \
    ' --output=' + output + \
    ' --ufile=' + inputUS + \
    ' --nfile=' + inputUN + \
    ' --uwfile=' + inputUW + \
    ' --nwfile=' + inputNW + \
    ' --rvel=' + rvelwave + \
    ' --tell=' + tellwave + \
    ' --snr=' + snrfile + \
    ' --sres=' + sresfile + \
    ' --version="' + operaversion + '"' + \
    ' --date="' + date + '"' + \
    ' --spectrumtype=21 --compressiontype=21'
    
    return commandline
##########################################

#### Function to generate the final product in FITS format out of SPC files: ####
def operaCreateProduct_SPC_Intensity(Dirs, objectfile, spectrumfile, output, rvelwave, tellwave, snrfile, sresfile, operaversion, date) :
    
    commandline = Dirs.EXE + 'operaCreateProduct' + \
    ' --input=' + objectfile + \
    ' --spectrumfile=' + spectrumfile + \
    ' --output=' + output + \
    ' --rvel=' + rvelwave + \
    ' --tell=' + tellwave + \
    ' --snr=' + snrfile + \
    ' --sres=' + sresfile + \
    ' --version="' + operaversion + '"' + \
    ' --date="' + date + '"' + \
    ' --spectrumtype=21 --compressiontype=21'

    return commandline
##########################################

#### Function to generate the final Polar product in FITS format out of LE files: ####
def operaCreateProduct_LE_Polar(Dirs, objectfile, output, inputUS, inputUN, inputUW, inputNW, rvelwave, tellwave, snrfile, sresfile, operaversion, date) :
    
    commandline = Dirs.EXE + 'operaCreateProduct' + \
    ' --input=' + objectfile + \
    ' --output=' + output + \
    ' --ufile=' + inputUS + \
    ' --nfile=' + inputUN + \
    ' --uwfile=' + inputUW + \
    ' --nwfile=' + inputNW + \
    ' --rvel=' + rvelwave + \
    ' --tell=' + tellwave + \
    ' --snr=' + snrfile + \
    ' --sres=' + sresfile + \
    ' --version="' + operaversion + '"' + \
    ' --date="' + date + '"' + \
    ' --spectrumtype=23 --compressiontype=21'
    
    return commandline
##########################################

#### Function to generate the final Polar product in FITS format out of POL files: ####
def operaCreateProduct_SPC_Polar(Dirs, objectfile, polfile, output, rvelwave, tellwave, snrfile, sresfile, operaversion, date) :
    
    commandline = Dirs.EXE + 'operaCreateProduct' + \
    ' --input=' + objectfile + \
    ' --spectrumfile=' + polfile + \
    ' --output=' + output + \
    ' --rvel=' + rvelwave + \
    ' --tell=' + tellwave + \
    ' --snr=' + snrfile + \
    ' --sres=' + sresfile + \
    ' --version="' + operaversion + '"' + \
    ' --date="' + date + '"' + \
    ' --spectrumtype=23 --compressiontype=21'
    
    return commandline
##########################################
