#!/bin/csh
# Description: simple script to run calibrations on ESPaDOnS data (for Star-only and Star+Sky)
# Author: Eder Martioli
# Laboratorio Nacional de Astrofisica, Brazil
# Apr 2014
#

set EXECUTECALIBRATION="OK"
set EXECUTEREDUCTION=""

####### Set variables ######
# Change variables below to reduce a new data set
# Note that:
#     It assumes all calibrations and object raw images are placed in the DATADIR directory
set OP_HOME="/Users/edermartioli/opera-1.0/"
set CONFIGDIR="$OP_HOME/config/"
set EXE="$OP_HOME/bin/"
#set DATAROOTDIR="/data/espadons/"
set DATAROOTDIR="/data/GRACES/"
##########################

####### NIGHT dir ########
set NIGHT="2011-07-14" # this is an example
##########################

####### Options for instrument mode ######
set STARONLY=1
set STARPLUSSKY=2
####### Options for read out speed ######
set FASTREADOUT=1
set NORMALREADOUT=2
set SLOWREADOUT=3
##########################

####### Set instrument configuration ######
set INSTRUMENTMODE="$STARONLY"  # set this manually to either 1=STARONLY or 2=STARPLUSSKY
set READOUTSPEEDCHOICE="$NORMALREADOUT" # set this manually to either FASTREADOUT, NORMALREADOUT, or SLOWREADOUT
set INSTRUMENT="ESPaDOnS"
set DETECTOR="OLAPA"
##########################

####### Below the $argv[1] overwrites night dir, $argv[2] overwrites instmode,
####### and $argv[3] overwrites readout speed from command line inputs
if($#argv)  then
    if ($argv[1] != "") then
        set NIGHT="$argv[1]"
        set INSTRUMENTMODE="$argv[2]"
        set READOUTSPEEDCHOICE="$argv[3]"
        echo "Input Night: $NIGHT"
    endif
endif
##########################

####### Set directories ######
set DATADIR="$DATAROOTDIR/$NIGHT/"
set DEFAULTDIR="/Users/edermartioli/Reductions/default/"
#set PRODUCTDIR="/Users/edermartioli/Reductions/Espadons/$NIGHT/"
set PRODUCTDIR="/Users/edermartioli/Reductions/GRACES/$NIGHT/"
##########################

####### Set read out speed ######
if ($READOUTSPEEDCHOICE == 1) then
    echo "Read out speed is Fast"
    set READOUTSPEED="Fast: 4.70e noise, 1.60e/ADU, 32s"
    set DEFAULTGAIN="1.6"
    set DEFAULTNOISE="4.14"
else if ($READOUTSPEEDCHOICE == 2) then
    echo "Read out speed is Normal"
    set READOUTSPEED="Normal: 4.20e noise, 1.30e/ADU, 38s"
    set DEFAULTGAIN="1.3"
    set DEFAULTNOISE="3.8"
else if ($READOUTSPEEDCHOICE == 3) then
    echo "Read out speed is Slow"
    set READOUTSPEED="Slow: 2.90e noise, 1.20e/ADU, 60s"
    set DEFAULTGAIN="1.1"
    set DEFAULTNOISE="2.98"
else
    echo "No read out speed has been set"
    exit
endif
##########################

####### Set instrument mode specific configuration ######
if ($INSTRUMENTMODE == 1) then
    echo "Instrument mode is Star-Only"
    set INSTRUMENTMODEKEY="Spectroscopy, star only, R=80,000"
    set STARPLUSKYMODEFLAG=0
    set SPCMODULE="operaStarOnly"
    set recenterIPUsingSliceSymmetry="0"
    set NUMBEROFSLICES="6"
    set NUMBEROFBEAMS="1"
    set SPECTRALRESOLUTION="80000"
    set LESPECTRUMTYPE=21
else if ($INSTRUMENTMODE == 2) then
    echo "Instrument mode is Star+Sky"
    set INSTRUMENTMODEKEY="Spectroscopy, star+sky, R=65,000"
    set STARPLUSKYMODEFLAG=1
    set SPCMODULE="operaStarPlusSky"
    set recenterIPUsingSliceSymmetry="1"
    set NUMBEROFSLICES="6"
    set NUMBEROFBEAMS="2"
    set SPECTRALRESOLUTION="65000"
    set LESPECTRUMTYPE=20
else
    echo "No instrument mode has been set."
    exit
endif
##########################

####### Set OBSTYPE keywords ######
set BIASKEYWORD="BIAS"
set FLATKEYWORD="FLAT"
set COMPKEYWORD="COMPARISON"
set OBJECTKEYWORD="OBJECT"
##########################

####### Set config files ######
set BADPIXELMASK=$CONFIGDIR"badpix_olapa-a.fits.gz"
set THARATLASLINES=$CONFIGDIR"thar_MM201006.dat.gz"
set THARATLASSPECTRUM=$CONFIGDIR"LovisPepe_ThArAtlas.dat.gz"
set WAVEFIRSTGUESS=$CONFIGDIR"wcal_ref.dat.gz"

set SOLARTYPEWAVELENGTHMASK=$CONFIGDIR"wavelengthMaskForUncalContinuumDetection_SolarTypeStars.txt"
set TELLURICLINES=$CONFIGDIR"opera_HITRAN08-extracted.par.gz"
set TELLURICSPECTRUM=$CONFIGDIR"KPNO_atmtrans.dat.gz"
##################################

###############################################
###############################################
######### SET CALIBRATION PRODUCTS ############
###############################################
set BIASLIST=$NIGHT"_bias.list"
set FLATLIST=$NIGHT"_flat.list"
set COMPLIST=$NIGHT"_comp.list"
set MASTERBIAS=$NIGHT"_masterbias.fits.gz"
set MASTERFLAT=$NIGHT"_masterflat.fits.gz"
set MASTERCOMP=$NIGHT"_mastercomp.fits.gz"
##
set GAINPRODUCT=$NIGHT".gain.gz"
set ORDERSPACINGPRODUCT=$NIGHT".ordp.gz"
set GEOMETRYPRODUCT=$NIGHT".geom.gz"
set INSTRUMENTPROFILEPRODUCT=$NIGHT".prof.gz"
set APERTUREPRODUCT=$NIGHT".aper.gz"
#
set COMPEXTRACTEDSPECTRUM=$NIGHT"_comp.e.gz"
set FLATEXTRACTEDSPECTRUM=$NIGHT"_flat.e.gz"
set FLATFLUXCALIBRATIONSPECTRUM=$NIGHT"_flat.fcal.gz"
#
set FIRSTWAVELENGTHPRODUCT=$NIGHT".wcar.gz"
set WAVELENGTHPRODUCT=$NIGHT".wcal.gz"
###############################################

###### Print out parameters ######
echo "Running OPERA-1.0 pipeline: Calibrations"
echo " "
echo "NIGHT = $NIGHT"
echo ""
##################################

###############################################
###############################################
####### START C A L I B R A T I O N ###########
###############################################
###############################################
if ($EXECUTECALIBRATION == "OK") then
###############################################
echo "--------"
echo "STARTING CALIBRATION"
echo ""

####### Create file lists for bias, flat, and comp ######
echo "--------"
echo "Creating bias list: $BIASLIST"
echo "Creating flat list: $FLATLIST"
echo "Creating comp list: $COMPLIST"
echo ""
$EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME INSTMODE EREADSPD DETECTOR OBSTYPE" INSTRUME=$INSTRUMENT INSTMODE="$INSTRUMENTMODEKEY" EREADSPD="$READOUTSPEED" DETECTOR=$DETECTOR OBSTYPE=$BIASKEYWORD > $BIASLIST
$EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME INSTMODE EREADSPD DETECTOR OBSTYPE" INSTRUME=$INSTRUMENT INSTMODE="$INSTRUMENTMODEKEY" EREADSPD="$READOUTSPEED" DETECTOR=$DETECTOR OBSTYPE=$FLATKEYWORD > $FLATLIST
$EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME INSTMODE EREADSPD DETECTOR OBSTYPE" INSTRUME=$INSTRUMENT INSTMODE="$INSTRUMENTMODEKEY" EREADSPD="$READOUTSPEED" DETECTOR=$DETECTOR OBSTYPE=$COMPKEYWORD > $COMPLIST
###############################

####### Create masterimages for bias, flat, and comp ######
echo "--------"
echo "Creating master bias : $MASTERBIAS"
echo "Creating master flat : $MASTERFLAT"
echo "Creating master comp : $MASTERCOMP"
echo ""
$EXE/operaMasterBias --output=$MASTERBIAS --list=$BIASLIST
$EXE/operaMasterFlat --output=$MASTERFLAT --list=$FLATLIST
$EXE/operaMasterComparison --output=$MASTERCOMP --list=$COMPLIST --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --combineMethod=1 --saturationLimit=65535 --outputExposureTime=60 --truncateOuputFluxToSaturation=1 --expTimeFITSKeyword=EXPTIME
###############################

####### Calculate detector gain and noise ######
echo "--------"
echo "Creating gain and noise calibration product : $GAINPRODUCT"
echo ""
$EXE/operaGain --output=$GAINPRODUCT --listofbiasimgs=$BIASLIST --listofflatimgs=$FLATLIST --DATASEC="1 2048 1 4608" --badpixelmask=$BADPIXELMASK --defaultgain=$DEFAULTGAIN --defaultnoise=$DEFAULTNOISE --numberofamplifiers=1 --maximages=12 --subwindow="100 800 500 3000" --gainMinPixPerBin=1000 --gainMaxNBins=100 --gainLowestCount=1000 --gainHighestCount=30000
###############################

####### Calculate order spacing calibration ######
echo "--------"
echo "Creating order spacing calibration product: $ORDERSPACINGPRODUCT"
echo ""
$EXE/operaOrderSpacingCalibration --orderspacingoutput=$ORDERSPACINGPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --subformat="8 2040 3 4600" --aperture=30 --numberOfsamples=30 --sampleCenterPosition=2300
###############################

####### Calculate geometry calibration ######
echo "--------"
echo "Creating geometry calibration product: $GEOMETRYPRODUCT"
echo ""
$EXE/operaGeometryCalibration --outputGeomFile=$GEOMETRYPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --subformat="8 2040 3 4600" --aperture=30 --detectionMethod=2 --FFTfilter=0 --nsamples=3 --maxorders=44 --minordertouse=18 --orderOfTracingPolynomial=3 --binsize=25 --colDispersion=1 --invertOrders=1 --recenterIPUsingSliceSymmetry=$recenterIPUsingSliceSymmetry --totalNumberOfSlices=$NUMBEROFSLICES --inputOrderSpacing=$ORDERSPACINGPRODUCT --referenceOrderNumber=55 --referenceOrderSeparation=66.7 --referenceOrderSamplePosition=2300
###############################

####### Calculate instrument profile calibration ######
echo "--------"
echo "Creating instrument profile calibration product: $INSTRUMENTPROFILEPRODUCT"
echo ""
$EXE/operaInstrumentProfileCalibration --outputProf=$INSTRUMENTPROFILEPRODUCT --geometryfilename=$GEOMETRYPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --mastercomparison=$MASTERCOMP --badpixelmask=$BADPIXELMASK --ipDimensions="30 5 6 5" --binsize=110 --ordernumber=-999 --method=2 --tilt=-2.0 --gain=$GAINPRODUCT --referenceLineWidth=1.3 --spectralElementHeight=1.0 --maxthreads=4
###############################

####### Calculate aperture calibration ######
echo "--------"
echo "Creating aperture calibration product: $APERTUREPRODUCT"
echo ""
$EXE/operaExtractionApertureCalibration  --outputApertureFile=$APERTUREPRODUCT --inputgeom=$GEOMETRYPRODUCT --inputprof=$INSTRUMENTPROFILEPRODUCT --inputorderspacing=$ORDERSPACINGPRODUCT --numberOfBeams=$NUMBEROFBEAMS --gapBetweenBeams=0 --apertureWidth=28 --apertureHeight=0.6923 --backgroundAperture=1.0 --pickImageRow=0 --nRowSamples=10 --xbin=10 -v
###############################

####### Extract comparison and flat-field spectra ######
echo "--------"
echo "Creating comparison spectrum: $COMPEXTRACTEDSPECTRUM"
echo "Creating flat-field spectrum: $FLATEXTRACTEDSPECTRUM"
echo ""
$EXE/operaExtraction --inputImage=$MASTERCOMP --outputSpectraFile=$COMPEXTRACTEDSPECTRUM  --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --spectrumtype=5 --spectrumtypename=RawBeamSpectrum --starplusskymode=0  --maxthreads=4 -v
$EXE/operaExtraction --outputSpectraFile=$FLATEXTRACTEDSPECTRUM --inputImage=$MASTERFLAT --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --backgroundBinsize=300 --sigmaclip=6 --removeBackground=0 --iterations=3 --onTargetProfile=1 --usePolynomialFit=0 --starplusskymode=0 --maxthreads=4 -v
###############################

####### Wavelength calibration ######
echo "--------"
echo "Creating 1st wavelength calibration product: $FIRSTWAVELENGTHPRODUCT"
echo ""
$EXE/operaWavelengthCalibration --outputWaveFile=$FIRSTWAVELENGTHPRODUCT --atlas_lines=$THARATLASLINES --atlas_spectrum=$THARATLASSPECTRUM --uncalibrated_spectrum=$COMPEXTRACTEDSPECTRUM --uncalibrated_linewidth=1.5 --inputGeomFile=$GEOMETRYPRODUCT --inputWaveFile=$WAVEFIRSTGUESS --parseSolution=0 --ParRangeSizeInPerCent=0.5 --NpointsPerPar=1000 --maxNIter=40 --minNumberOfLines=40 --maxorderofpolynomial=4 --dampingFactor=0.85 --initialAcceptableMismatch=2.0 --nsigclip=2.5 --normalizeUncalibratedSpectrum=0 --normalizationBinSize=180 --LocalMaxFilterWidth=6
echo ""
echo "Creating wavelength calibration product after stitching orders together: $WAVELENGTHPRODUCT"
echo ""
$EXE/operaStitchOrders --outputWaveFile=$WAVELENGTHPRODUCT --inputSpectrum=$COMPEXTRACTEDSPECTRUM --inputWaveFile=$FIRSTWAVELENGTHPRODUCT --orderOfReference=51 --DWavelengthRange=0.1 --DWavelengthStep=0.00005 --XCorrelationThreshold=0.1 --sigmaThreshold=2.0
###############################

####### Create flat-field flux calibration spectrum ######
echo ""
echo "Creating flat-field flux calibration spectrum: $FLATFLUXCALIBRATIONSPECTRUM"
echo ""
$EXE/operaCreateFlatFieldFluxCalibration --outputFluxCalibrationFile=$FLATFLUXCALIBRATIONSPECTRUM --inputMasterFlatSpectrum=$FLATEXTRACTEDSPECTRUM --wavelengthCalibration=$WAVELENGTHPRODUCT --binsize=500 --wavelengthForNormalization=548 -v
###############################

echo ""
echo "END CALIBRATION"
echo "--------"
###############################################
########### END C A L I B R A T I O N #########
###############################################
###############################################
endif
###############################################


###############################################
###############################################
########## SET REDUCTION PRODUCTS #############
###############################################
set OBJECTLIST=$NIGHT"_object.list"
################################################

###### Print out parameters ######
echo "Running OPERA-1.0 pipeline: Reduction"
echo " "
echo "NIGHT = $NIGHT"
echo ""
##################################

###############################################
###############################################
######### START R E D U C T I O N #############
###############################################
###############################################
if ($EXECUTEREDUCTION == "OK") then
###############################################
echo "--------"
echo "STARTING REDUCTION"
echo ""

####### Create list of Object files ######
# input: directory and qualifiers
# This module creates a list of object files for reduction
echo "--------"
echo "Creating object list: $OBJECTLIST"
echo ""
$EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME INSTMODE EREADSPD DETECTOR OBSTYPE" INSTRUME=$INSTRUMENT INSTMODE="$INSTRUMENTMODEKEY" EREADSPD="$READOUTSPEED" DETECTOR=$DETECTOR OBSTYPE=$OBJECTKEYWORD > $OBJECTLIST
##########################################

foreach OBJIMAGE (`cat $OBJECTLIST`)

    ####### Figure out image base name ######
    set DATADIRLEN=`echo $DATADIR | awk '{print length($1)}'`
    set OBJIMGBASENAME=`echo $OBJIMAGE | awk '{print substr($0,'$DATADIRLEN'+1,8)}'`
    ##########################################

    ####### Figure out info from object image header ######
    set OBJECTNAME=`$EXE/operagetheader --keyword=OBJECT $OBJIMAGE`
    set MJDATE=`$EXE/operagetheader --keyword=MJDATE $OBJIMAGE`
    set absra_center=`$EXE/operagetheader --keyword=RA_DEG $OBJIMAGE`
    set absdec_center=`$EXE/operagetheader --keyword=DEC_DEG $OBJIMAGE`
    echo "Object=$OBJECTNAME JD=$MJDATE RA=$absra_center Dec=$absdec_center"
    ##########################################

    ####### Extract Object Spectrum ######
    set OBJECTSPECTRUM=$OBJIMGBASENAME".e.gz"
    echo "--------"
    echo "Extracting object spectrum product: $OBJECTSPECTRUM "
    echo ""
    $EXE/operaExtraction --outputSpectraFile=$OBJECTSPECTRUM --inputImage=$OBJIMAGE --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --backgroundBinsize=300 --sigmaclip=6 --onTargetProfile=1 --starplusskymode=$STARPLUSKYMODEFLAG --usePolynomialFit=0 --removeBackground=0 --iterations=3 --maxthreads=4 -v

    ####### Calculate telluric wavelength correction ######
    set TELLWAVECAL=$OBJIMGBASENAME".tell.gz"
    echo "--------"
    echo "Calculating telluric wavelength correction product: $TELLWAVECAL "
    echo ""
    $EXE/operaTelluricWavelengthCorrection --outputWaveFile=$TELLWAVECAL --inputObjectSpectrum=$OBJECTSPECTRUM --inputWaveFile=$WAVELENGTHPRODUCT --telluric_lines=$TELLURICLINES --telluric_spectrum=$TELLURICSPECTRUM --spectralResolution=$SPECTRALRESOLUTION --initialWavelengthRange=0.1 --initialWavelengthStep=0.002 --XCorrelationThreshold=0.1 --subtractCentralWavelength=1 --normalizationBinsize=110 --sigmaThreshold=1.25 -v

    ####### Calculate barycentric wavelength correction ######
    set BARYWAVECAL=$OBJIMGBASENAME".rvel.gz"
    echo "--------"
    echo "Calculating barycentric wavelength correction product: $BARYWAVECAL "
    echo ""
    $EXE/operaBarycentricWavelengthCorrection --outputRVelFile=$BARYWAVECAL --inputWaveFile=$WAVELENGTHPRODUCT --observatory_coords="19:49:36 -155:28:18" --object_coords="$absra_center $absdec_center" --observatory_elevation=4207 --MJDTime=$MJDATE -v

    ####### Calculate final calibrated spectrum *.spc ######
    set CALIBRATEDSPECTRUM=$OBJIMGBASENAME".spc.gz"
    echo "--------"
    echo "Calculating final calibrated spectrum product: $CALIBRATEDSPECTRUM "
    echo ""
    $EXE/$SPCMODULE --outputCalibratedSpectrum=$CALIBRATEDSPECTRUM --inputUncalibratedSpectrum=$OBJECTSPECTRUM --radialvelocitycorrection=$BARYWAVECAL --telluriccorrection=$TELLWAVECAL --spectrumtype=17 --wavelengthCalibration=$WAVELENGTHPRODUCT --inputFlatFluxCalibration=$FLATFLUXCALIBRATIONSPECTRUM --inputWavelengthMaskForUncalContinuum=$SOLARTYPEWAVELENGTHMASK --object="$OBJECTNAME" --numberOfPointsInUniformSample=150 --normalizationBinsize=750 --AbsoluteCalibration=0 --etime=1.0 -v

    ####### Generate LE formats ######
    set LESPCNW=$OBJIMGBASENAME".inw.s.gz"
    set LESPCU=$OBJIMGBASENAME".iu.s.gz"
    set LESPCUW=$OBJIMGBASENAME".iuw.s.gz"
    set LESPCN=$OBJIMGBASENAME".in.s.gz"
    echo "--------"
    echo "Calculating telluric wavelength correction: $OBJECTSPECTRUM "
    echo ""
    $EXE/operaGenerateLEFormats --outputLEfilename=$LESPCNW --inputOperaSpectrum=$CALIBRATEDSPECTRUM --LibreEspritSpectrumType=$LESPECTRUMTYPE --object="$OBJECTNAME" --fluxType=2 --wavelengthType=3 -v
    $EXE/operaGenerateLEFormats --outputLEfilename=$LESPCU --inputOperaSpectrum=$CALIBRATEDSPECTRUM --LibreEspritSpectrumType=$LESPECTRUMTYPE --object="$OBJECTNAME" --fluxType=3 --wavelengthType=3 -v
    $EXE/operaGenerateLEFormats --outputLEfilename=$LESPCUW --inputOperaSpectrum=$CALIBRATEDSPECTRUM --LibreEspritSpectrumType=$LESPECTRUMTYPE --object="$OBJECTNAME" --fluxType=3 --wavelengthType=4 -v
    $EXE/operaGenerateLEFormats --outputLEfilename=$LESPCN --inputOperaSpectrum=$CALIBRATEDSPECTRUM --LibreEspritSpectrumType=$LESPECTRUMTYPE --object="$OBJECTNAME" --fluxType=2 --wavelengthType=3 -v
    ###############################################

end

###############################################
echo ""
echo "END REDUCTION"
echo "--------"
###############################################
########### END R E D U C T I O N #############
###############################################
###############################################
endif
###############################################


echo " "
echo "The pipeline ran successfully!"
echo " "

exit


