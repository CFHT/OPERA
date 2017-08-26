#!/bin/csh
# Description: script to reduce GRACES data (for Star-only and Star+Sky modes)
# Author: Eder Martioli
# Laboratorio Nacional de Astrofisica, Brazil
# Last update on July 28 2014
#

set EXECUTECALIBRATION=0
set EXECUTEREDUCTION=0

####### Set variables ######
# Change variables below to reduce a new data set
# Note that:
#     It assumes all calibrations and object raw images are placed in the DATADIR directory
set OP_HOME="/Users/edermartioli/opera-1.0/"
set DATAROOTDIR="/data/GRACES/"
set PRODUCTROOTDIR="/Users/edermartioli/Reductions/GRACES/"
##########################

####### NIGHT dir ########
set NIGHT="20140506" # this is an example
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
set INSTRUMENT="GRACES"
set DETECTOR="OLAPA"
set ALLOWANYREADOUTSPEED=0
set ALLOWUSEDEFAULTCAL=0
##########################

####### Below the $argv[1] overwrites night dir, $argv[2] overwrites instmode,
####### and $argv[3] overwrites readout speed from command line inputs
if($#argv)  then
    if ($argv[1] != "") then
        set OP_HOME="$argv[1]"
        set DATAROOTDIR="$argv[2]"
        set PRODUCTROOTDIR="$argv[3]"
        set NIGHT="$argv[4]"
        set INSTRUMENTMODE="$argv[5]"
        set READOUTSPEEDCHOICE="$argv[6]"
        set EXECUTECALIBRATION="$argv[7]"
        set EXECUTEREDUCTION="$argv[8]"
        set ALLOWANYREADOUT="$argv[9]"
        set ALLOWUSEDEFAULTCAL="$argv[10]"
    endif
endif
##########################

####### Set dir variables ######
set CONFIGDIR="$OP_HOME/config/"
set EXE="$OP_HOME/bin/"
##########################

####### Set directories ######
set DATADIR="$DATAROOTDIR/$NIGHT/"
set PRODUCTDIR="$PRODUCTROOTDIR/$NIGHT/"
##########################

####### Set read out speed ######
if ($READOUTSPEEDCHOICE == 1) then
    echo "Read out speed is Fast"
    set READOUTSPEEDSHORTNAME="fast"
    set READOUTSPEED="Fast: 4.70e noise, 1.60e/ADU, 32s"
    set DEFAULTGAIN="1.6"
    set DEFAULTNOISE="4.14"
else if ($READOUTSPEEDCHOICE == 2) then
    echo "Read out speed is Normal"
    set READOUTSPEEDSHORTNAME="normal"
    set READOUTSPEED="Normal: 4.20e noise, 1.30e/ADU, 38s"
    set DEFAULTGAIN="1.3"
    set DEFAULTNOISE="3.8"
else if ($READOUTSPEEDCHOICE == 3) then
    echo "Read out speed is Slow"
    set READOUTSPEEDSHORTNAME="slow"
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
    set INSTRUMENTMODESHORTNAME="StarOnly"
    set INSTRUMENTMODEKEY="FOURSLICE"
    set STARPLUSKYMODEFLAG=0
    set INVERTSKYFIBERFLAG=""
    set SPCMODULE="operaStarOnly"
    set recenterIPUsingSliceSymmetry="0"
    set NUMBEROFSLICES="4"
    set NUMBEROFBEAMS="1"
    set SPECTRALRESOLUTION="60000"
    set LESPECTRUMTYPE=21
    set ORDSPCAPERTURE=30
    set SPACINGREFERENCEORDERNUMBER=45
    set SPACINGREFERENCEORDERSEPARATION=58.0
    set GEOMAPERTURE=30
    set GEOMMAXNORDERS=42
    set GEOMMINORDERTOUSE=17
    set REFERENCELINEWIDTH=1.8
    set IPDIMENSIONS="32 5 7 5"
    set APERAPERTURE=32
    set APERGAP=0
    set WAVEUNCALLINEWIDTH=1.6

else if ($INSTRUMENTMODE == 2) then
    echo "Instrument mode is Star+Sky"
    set INSTRUMENTMODESHORTNAME="StarPlusSky"
    set INSTRUMENTMODEKEY="TWOSLICE"
    set STARPLUSKYMODEFLAG=1
    set INVERTSKYFIBERFLAG="--starplusskyInvertSkyFiber=1"
    set SPCMODULE="operaStarPlusSky"
    set recenterIPUsingSliceSymmetry="1"
    set NUMBEROFSLICES="4"
    set NUMBEROFBEAMS="2"
    set SPECTRALRESOLUTION="40000"
    set LESPECTRUMTYPE=20
    set ORDSPCAPERTURE=36
    set SPACINGREFERENCEORDERNUMBER=45
    set SPACINGREFERENCEORDERSEPARATION=58.0
    set GEOMAPERTURE=31
    set GEOMMAXNORDERS=42
    set GEOMMINORDERTOUSE=17
    set REFERENCELINEWIDTH=2.5
    set IPDIMENSIONS="34 5 7 5"
    set APERAPERTURE=34
    set APERGAP=1
    set WAVEUNCALLINEWIDTH=1.8

else
    echo "No instrument mode has been set."
    exit
endif
##########################

set INSTCONFIGPREFIX=$NIGHT"_"$INSTRUMENTMODESHORTNAME"_"$READOUTSPEEDSHORTNAME

####### Set OBSTYPE keywords ######
set BIASKEYWORD="BIAS"
set FLATKEYWORD="FLAT"
set COMPKEYWORD="COMPARISON"
set OBJECTKEYWORD="OBJECT"
##########################

####### Set config files ######
set BADPIXELMASK=$CONFIGDIR"badpix_olapa-a.fits.fz"
set THARATLASLINES=$CONFIGDIR"thar_MM201006.dat.gz"
set THARATLASSPECTRUM=$CONFIGDIR"LovisPepe_ThArAtlas.dat.gz"
set WAVEFIRSTGUESS=$CONFIGDIR"wcal_ref.dat.gz"

set SOLARTYPEWAVELENGTHMASK=$CONFIGDIR"wavelengthMaskForUncalContinuumDetection_SolarTypeStars.txt"
set ATYPEWAVELENGTHMASK=$CONFIGDIR"wavelengthMaskForUncalContinuumDetection.txt"
set ATLASWAVELENGTHMASK=$CONFIGDIR"wavelengthMaskForRefContinuumDetection.txt"
set TELLURICLINES=$CONFIGDIR"opera_HITRAN08-extracted.par.gz"
set TELLURICSPECTRUM=$CONFIGDIR"KPNO_atmtrans.dat.gz"
##################################

###############################################
###############################################
######### SET CALIBRATION PRODUCTS ############
###############################################
set BIASLIST=$PRODUCTDIR$INSTCONFIGPREFIX"_bias.list"
set FLATLIST=$PRODUCTDIR$INSTCONFIGPREFIX"_flat.list"
set COMPLIST=$PRODUCTDIR$INSTCONFIGPREFIX"_comp.list"
set MASTERBIAS=$PRODUCTDIR$INSTCONFIGPREFIX"_masterbias.fits.gz"
set MASTERFLAT=$PRODUCTDIR$INSTCONFIGPREFIX"_masterflat.fits.gz"
set MASTERCOMP=$PRODUCTDIR$INSTCONFIGPREFIX"_mastercomp.fits.gz"
##
set GAINPRODUCT=$PRODUCTDIR$INSTCONFIGPREFIX".gain.gz"
set ORDERSPACINGPRODUCT=$PRODUCTDIR$INSTCONFIGPREFIX".ordp.gz"
set GEOMETRYPRODUCT=$PRODUCTDIR$INSTCONFIGPREFIX".geom.gz"
set INSTRUMENTPROFILEPRODUCT=$PRODUCTDIR$INSTCONFIGPREFIX".prof.gz"
set APERTUREPRODUCT=$PRODUCTDIR$INSTCONFIGPREFIX".aper.gz"
#
set COMPEXTRACTEDSPECTRUM=$PRODUCTDIR$INSTCONFIGPREFIX"_comp.e.gz"
set FLATEXTRACTEDSPECTRUM=$PRODUCTDIR$INSTCONFIGPREFIX"_flat.e.gz"
set FLATFLUXCALIBRATIONSPECTRUM=$PRODUCTDIR$INSTCONFIGPREFIX"_flat.fcal.gz"
#
set FIRSTWAVELENGTHPRODUCT=$PRODUCTDIR$INSTCONFIGPREFIX".wcar.gz"
set WAVELENGTHPRODUCT=$PRODUCTDIR$INSTCONFIGPREFIX".wcal.gz"
###############################################

########## Reduction Parameters ##########
set WAVELENGTHFORNORMALIZATION=548
##########################################

###### Print out parameters ######
echo "Running GRACES Pipeline 1.0: Calibrations"
echo " "
echo "Processing NIGHT = $NIGHT"
echo "For INSTRUMENT MODE = $INSTRUMENTMODESHORTNAME"
echo "For READOUT SPEED  = $READOUTSPEEDSHORTNAME"
echo ""
##################################

###############################################
###############################################
####### START C A L I B R A T I O N ###########
###############################################
###############################################
if ($EXECUTECALIBRATION == 1) then
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

echo "$EXE/operaQueryImageInfo --directory=$DATADIR -q "'"INSTRUME EREADSPD OBSTYPE"'" INSTRUME=$INSTRUMENT EREADSPD="'"'"$READOUTSPEED"'"'" OBSTYPE=$BIASKEYWORD > $BIASLIST"
$EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME EREADSPD OBSTYPE" INSTRUME=$INSTRUMENT EREADSPD="$READOUTSPEED" OBSTYPE=$BIASKEYWORD > $BIASLIST

if ($ALLOWANYREADOUTSPEED == 0) then
    echo "$EXE/operaQueryImageInfo --directory=$DATADIR -q "'"INSTRUME GSLICER EREADSPD OBSTYPE"'" INSTRUME=$INSTRUMENT GSLICER="'"'"$INSTRUMENTMODEKEY"'"'" EREADSPD="'"'"$READOUTSPEED"'"'" OBSTYPE=$FLATKEYWORD > $FLATLIST"
    $EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME GSLICER EREADSPD OBSTYPE" INSTRUME=$INSTRUMENT GSLICER="$INSTRUMENTMODEKEY" EREADSPD="$READOUTSPEED" OBSTYPE=$FLATKEYWORD > $FLATLIST
    echo "$EXE/operaQueryImageInfo --directory=$DATADIR -q "'"INSTRUME GSLICER EREADSPD OBSTYPE"'" INSTRUME=$INSTRUMENT GSLICER="'"'"$INSTRUMENTMODEKEY"'"'" EREADSPD="'"'"$READOUTSPEED"'"'" OBSTYPE=$COMPKEYWORD > $COMPLIST"
    $EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME GSLICER EREADSPD OBSTYPE" INSTRUME=$INSTRUMENT GSLICER="$INSTRUMENTMODEKEY" EREADSPD="$READOUTSPEED" OBSTYPE=$COMPKEYWORD > $COMPLIST
else
    echo "$EXE/operaQueryImageInfo --directory=$DATADIR -q "'"INSTRUME GSLICER OBSTYPE"'" INSTRUME=$INSTRUMENT GSLICER="'"'"$INSTRUMENTMODEKEY"'"'" OBSTYPE=$FLATKEYWORD > $FLATLIST"
    $EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME GSLICER OBSTYPE" INSTRUME=$INSTRUMENT GSLICER="$INSTRUMENTMODEKEY" OBSTYPE=$FLATKEYWORD > $FLATLIST
    echo "$EXE/operaQueryImageInfo --directory=$DATADIR -q "'"INSTRUME GSLICER OBSTYPE"'" INSTRUME=$INSTRUMENT GSLICER="'"'"$INSTRUMENTMODEKEY"'"'" OBSTYPE=$COMPKEYWORD > $COMPLIST"
    $EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME GSLICER OBSTYPE" INSTRUME=$INSTRUMENT GSLICER="$INSTRUMENTMODEKEY" OBSTYPE=$COMPKEYWORD > $COMPLIST
endif

###############################

####### Create masterimages for bias, flat, and comp ######
echo "--------"
echo ""

if ( -s $BIASLIST ) then
    echo "Creating master bias : $MASTERBIAS"
    echo "$EXE/operaMasterBias --output=$MASTERBIAS --list=$BIASLIST -v"
    $EXE/operaMasterBias --output=$MASTERBIAS --list=$BIASLIST
else
    if (ALLOWUSEDEFAULTCAL == 0) then
        echo "No calibration available - exiting program!"
        exit
    else
        set MASTERBIAS=$OP_HOME"/DefaultCalibration/"$READOUTSPEEDSHORTNAME"_masterbias.fits.gz"
        echo "Using default bias: $MASTERBIAS"
    endif
endif


if ( -s $FLATLIST ) then
    echo "Creating master flat : $MASTERFLAT"
    echo "$EXE/operaMasterFlat --output=$MASTERFLAT --list=$FLATLIST -v"
    $EXE/operaMasterFlat --output=$MASTERFLAT --list=$FLATLIST
else
    if (ALLOWUSEDEFAULTCAL == 0) then
        echo "No calibration available - exiting program!"
        exit
    else
        set MASTERFLAT=$OP_HOME"/DefaultCalibration/"$INSTRUMENTMODESHORTNAME"_masterflat.fits.gz"
        echo "Using default flat: $MASTERFLAT"
    endif
endif


if ( -s $COMPLIST ) then
    echo "Creating master comp : $MASTERCOMP"
    echo "$EXE/operaMasterComparison --output=$MASTERCOMP --listofimages=$COMPLIST --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --combineMethod=1 --saturationLimit=65535 --outputExposureTime=40 --truncateOuputFluxToSaturation=1 --expTimeFITSKeyword=EXPTIME -v"
    $EXE/operaMasterComparison --output=$MASTERCOMP --listofimages=$COMPLIST --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --combineMethod=1 --saturationLimit=65535 --outputExposureTime=40 --truncateOuputFluxToSaturation=1 --expTimeFITSKeyword=EXPTIME -v
else
    if (ALLOWUSEDEFAULTCAL == 0) then
        echo "No calibration available - exiting program!"
        exit
    else
        set MASTERCOMP=$OP_HOME"/DefaultCalibration/"$INSTRUMENTMODESHORTNAME"_mastercomp.fits.gz"
        echo "Using default comparison: $MASTERCOMP"
    endif
endif
###############################

####### Figure out gain and noise calibration ######
if ( -s $FLATLIST && -s $BIASLIST) then

    ####### Calculate detector gain and noise ######
    echo "--------"
    echo "Creating gain and noise calibration product : $GAINPRODUCT"
    echo ""
    echo "$EXE/operaGain --output=$GAINPRODUCT --listofbiasimgs=$BIASLIST --listofflatimgs=$FLATLIST --DATASEC="'"1 2048 1 4608"'" --badpixelmask=$BADPIXELMASK --defaultgain=$DEFAULTGAIN --defaultnoise=$DEFAULTNOISE --numberofamplifiers=1 --maximages=12 --subwindow="'"100 800 500 4000"'" --gainMinPixPerBin=1000 --gainMaxNBins=100 --gainLowestCount=1000 --gainHighestCount=30000 -v"
    $EXE/operaGain --output=$GAINPRODUCT --listofbiasimgs=$BIASLIST --listofflatimgs=$FLATLIST --DATASEC="1 2048 1 4608" --badpixelmask=$BADPIXELMASK --defaultgain=$DEFAULTGAIN --defaultnoise=$DEFAULTNOISE --numberofamplifiers=1 --maximages=12 --subwindow="100 800 500 4000" --gainMinPixPerBin=1000 --gainMaxNBins=100 --gainLowestCount=1000 --gainHighestCount=30000 -v
    ###############################
else
    if (ALLOWUSEDEFAULTCAL == 0) then
        echo "No calibration available - exiting program!"
        exit
    else
        set GAINPRODUCT=$OP_HOME"/DefaultCalibration/"$READOUTSPEEDSHORTNAME".gain.gz"
        echo "Using default gain: $GAINPRODUCT"
    endif
endif
###############################

####### Calculate order spacing calibration ######
echo "--------"
set ORDSPCPLOTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_spcplot.eps"
set ORDSPCDATAFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_spcplot.dat"
set ORDSPCSCRIPTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_spcplot.gnu"
echo "Creating order spacing calibration product: $ORDERSPACINGPRODUCT"
echo ""
echo "$EXE/operaOrderSpacingCalibration --inputGainFile=$GAINPRODUCT --plotfilename=$ORDSPCPLOTFILE --datafilename=$ORDSPCDATAFILE --scriptfilename=$ORDSPCSCRIPTFILE --orderspacingoutput=$ORDERSPACINGPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --subformat="'"8 2040 3 4600"'" --aperture=$ORDSPCAPERTURE --numberOfsamples=30 --sampleCenterPosition=2300 --referenceOrderNumber=$SPACINGREFERENCEORDERNUMBER --referenceOrderSeparation=$SPACINGREFERENCEORDERSEPARATION -v"
$EXE/operaOrderSpacingCalibration --inputGainFile=$GAINPRODUCT --plotfilename=$ORDSPCPLOTFILE --datafilename=$ORDSPCDATAFILE --scriptfilename=$ORDSPCSCRIPTFILE --orderspacingoutput=$ORDERSPACINGPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --subformat="8 2040 3 4600" --aperture=$ORDSPCAPERTURE --numberOfsamples=10 --sampleCenterPosition=2300 --referenceOrderNumber=$SPACINGREFERENCEORDERNUMBER --referenceOrderSeparation=$SPACINGREFERENCEORDERSEPARATION -v
###############################

####### Calculate geometry calibration ######
echo "--------"
set GEOMPLOTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_geomplot.eps"
set GEOMDATAFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_geomplot.dat"
set GEOMSCRIPTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_geomplot.gnu"
echo "Creating geometry calibration product: $GEOMETRYPRODUCT"
echo ""
echo "$EXE/operaGeometryCalibration --inputGainFile=$GAINPRODUCT --plotfilename=$GEOMPLOTFILE --datafilename=$GEOMDATAFILE --scriptfilename=$GEOMSCRIPTFILE --outputGeomFile=$GEOMETRYPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --subformat="'"8 2040 3 4600"'" --aperture=$GEOMAPERTURE --detectionMethod=2 --FFTfilter=0 --nsamples=30 --maxorders=$GEOMMAXNORDERS --minordertouse=$GEOMMINORDERTOUSE --orderOfTracingPolynomial=3 --binsize=10 --colDispersion=1 --invertOrders=1 --recenterIPUsingSliceSymmetry=$recenterIPUsingSliceSymmetry --totalNumberOfSlices=$NUMBEROFSLICES --inputOrderSpacing=$ORDERSPACINGPRODUCT --referenceOrderSamplePosition=2300 -v"
$EXE/operaGeometryCalibration --inputGainFile=$GAINPRODUCT --plotfilename=$GEOMPLOTFILE --datafilename=$GEOMDATAFILE --scriptfilename=$GEOMSCRIPTFILE --outputGeomFile=$GEOMETRYPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --subformat="8 2040 3 4600" --aperture=$GEOMAPERTURE --detectionMethod=2 --FFTfilter=0 --nsamples=30 --maxorders=$GEOMMAXNORDERS --minordertouse=$GEOMMINORDERTOUSE --orderOfTracingPolynomial=3 --binsize=10 --colDispersion=1 --invertOrders=1 --recenterIPUsingSliceSymmetry=$recenterIPUsingSliceSymmetry --totalNumberOfSlices=$NUMBEROFSLICES --inputOrderSpacing=$ORDERSPACINGPRODUCT --referenceOrderSamplePosition=1000 -v
###############################

####### Calculate instrument profile calibration ######
echo "--------"
set PROFPLOTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_profplot.eps"
set PROFDATAFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_profplot.dat"
set PROFSCRIPTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_profplot.gnu"
echo "Creating instrument profile calibration product: $INSTRUMENTPROFILEPRODUCT"
echo ""
echo "$EXE/operaInstrumentProfileCalibration --plotfilename=$PROFPLOTFILE --datafilename=$PROFDATAFILE --scriptfilename=$PROFSCRIPTFILE --outputProf=$INSTRUMENTPROFILEPRODUCT --geometryfilename=$GEOMETRYPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --mastercomparison=$MASTERCOMP --badpixelmask=$BADPIXELMASK --ipDimensions="'"'"$IPDIMENSIONS"'"'" --binsize=100 --ordernumber=-999 --method=2 --tilt=-3.0 --gain=$GAINPRODUCT --referenceLineWidth=$REFERENCELINEWIDTH --spectralElementHeight=1.0 --maxthreads=4 -v"
$EXE/operaInstrumentProfileCalibration --plotfilename=$PROFPLOTFILE --datafilename=$PROFDATAFILE --scriptfilename=$PROFSCRIPTFILE --outputProf=$INSTRUMENTPROFILEPRODUCT --geometryfilename=$GEOMETRYPRODUCT --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --mastercomparison=$MASTERCOMP --badpixelmask=$BADPIXELMASK --ipDimensions="$IPDIMENSIONS" --binsize=100 --ordernumber=-999 --method=2 --tilt=-3.0 --gain=$GAINPRODUCT --referenceLineWidth=$REFERENCELINEWIDTH --spectralElementHeight=1.0 --maxthreads=4 -v
###############################

####### Calculate aperture calibration ######
echo "--------"
set APERPLOTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_aperplot.eps"
set APERDATAFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_aperplot.dat"
set APERSCRIPTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_aperplot.gnu"
echo "Creating aperture calibration product: $APERTUREPRODUCT"
echo ""
echo "$EXE/operaExtractionApertureCalibration  --plotfilename=$APERPLOTFILE --datafilename=$APERDATAFILE --scriptfilename=$APERSCRIPTFILE --outputApertureFile=$APERTUREPRODUCT --inputgeom=$GEOMETRYPRODUCT --inputprof=$INSTRUMENTPROFILEPRODUCT --inputorderspacing=$ORDERSPACINGPRODUCT --numberOfBeams=$NUMBEROFBEAMS --gapBetweenBeams=$APERGAP --apertureWidth=$APERAPERTURE --apertureHeight=0.6923 --backgroundAperture=1.0 --pickImageRow=0 --nRowSamples=10 --xbin=10 -v"
$EXE/operaExtractionApertureCalibration  --plotfilename=$APERPLOTFILE --datafilename=$APERDATAFILE --scriptfilename=$APERSCRIPTFILE --outputApertureFile=$APERTUREPRODUCT --inputgeom=$GEOMETRYPRODUCT --inputprof=$INSTRUMENTPROFILEPRODUCT --inputorderspacing=$ORDERSPACINGPRODUCT --numberOfBeams=$NUMBEROFBEAMS --gapBetweenBeams=$APERGAP --apertureWidth=$APERAPERTURE --apertureHeight=0.6923 --backgroundAperture=1.0 --pickImageRow=0 --nRowSamples=10 --xbin=10 -v
###############################

####### Extract comparison and flat-field spectra ######
echo "--------"
echo "Creating comparison spectrum: $COMPEXTRACTEDSPECTRUM"
echo "Creating flat-field spectrum: $FLATEXTRACTEDSPECTRUM"
echo ""
echo "$EXE/operaExtraction --inputImage=$MASTERCOMP --outputSpectraFile=$COMPEXTRACTEDSPECTRUM  --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --spectrumtype=5 --spectrumtypename=RawBeamSpectrum --starplusskymode=0  --maxthreads=4 -v"
$EXE/operaExtraction --inputImage=$MASTERCOMP --outputSpectraFile=$COMPEXTRACTEDSPECTRUM  --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --spectrumtype=5 --spectrumtypename=RawBeamSpectrum --starplusskymode=0  --maxthreads=4 -v
echo "$EXE/operaExtraction --outputSpectraFile=$FLATEXTRACTEDSPECTRUM --inputImage=$MASTERFLAT --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --backgroundBinsize=300 --sigmaclip=6 --removeBackground=0 --iterations=3 --onTargetProfile=1 --usePolynomialFit=0 --starplusskymode=0 --maxthreads=4 -v"
$EXE/operaExtraction --outputSpectraFile=$FLATEXTRACTEDSPECTRUM --inputImage=$MASTERFLAT --masterflat=$MASTERFLAT --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --backgroundBinsize=300 --sigmaclip=6 --removeBackground=0 --iterations=3 --onTargetProfile=1 --usePolynomialFit=0 --starplusskymode=0 --maxthreads=4 -v
###############################

####### Wavelength calibration ######
echo "--------"
echo "Creating 1st wavelength calibration product: $FIRSTWAVELENGTHPRODUCT"
set WAVEORDSPLOTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_waveordsplot.eps"
set WAVESPECPLOTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_wavespecplot.eps"
set WAVESPECSCRIPTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_wavespecplot.gnu"
set WAVEORDSCRIPTFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_waveordplot.gnu"
set WAVEORDSDATAFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_waveordsplot.dat"
set WAVEATLASDATAFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_waveatlasplot.dat"
set WAVECOMPDATAFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_wavecompplot.dat"
set WAVELINESDATAFILE=$PRODUCTDIR$INSTCONFIGPREFIX"_wavelinesplot.dat"
echo ""
echo "$EXE/operaWavelengthCalibration --ordersplotfilename=$WAVEORDSPLOTFILE --specplotfilename=$WAVESPECPLOTFILE --ordersscriptfilename=$WAVEORDSCRIPTFILE --specscriptfilename=$WAVESPECSCRIPTFILE --ordersdatafilename=$WAVEORDSDATAFILE --atlasdatafilename=$WAVEATLASDATAFILE --compdatafilename=$WAVECOMPDATAFILE --linesdatafilename=$WAVELINESDATAFILE --outputWaveFile=$FIRSTWAVELENGTHPRODUCT --atlas_lines=$THARATLASLINES --atlas_spectrum=$THARATLASSPECTRUM --uncalibrated_spectrum=$COMPEXTRACTEDSPECTRUM --uncalibrated_linewidth=$WAVEUNCALLINEWIDTH --inputGeomFile=$GEOMETRYPRODUCT --inputWaveFile=$WAVEFIRSTGUESS --parseSolution=0 --ParRangeSizeInPerCent=1.0 --NpointsPerPar=1000 --maxNIter=40 --minNumberOfLines=40 --maxorderofpolynomial=4 --dampingFactor=0.85 --initialAcceptableMismatch=1.5 --nsigclip=2.5 --normalizeUncalibratedSpectrum=0 --normalizationBinSize=180 --LocalMaxFilterWidth=6 -v"
$EXE/operaWavelengthCalibration --ordersplotfilename=$WAVEORDSPLOTFILE --specplotfilename=$WAVESPECPLOTFILE --ordersscriptfilename=$WAVEORDSCRIPTFILE --specscriptfilename=$WAVESPECSCRIPTFILE --ordersdatafilename=$WAVEORDSDATAFILE --atlasdatafilename=$WAVEATLASDATAFILE --compdatafilename=$WAVECOMPDATAFILE --linesdatafilename=$WAVELINESDATAFILE --outputWaveFile=$FIRSTWAVELENGTHPRODUCT --atlas_lines=$THARATLASLINES --atlas_spectrum=$THARATLASSPECTRUM --uncalibrated_spectrum=$COMPEXTRACTEDSPECTRUM --uncalibrated_linewidth=$WAVEUNCALLINEWIDTH --inputGeomFile=$GEOMETRYPRODUCT --inputWaveFile=$WAVEFIRSTGUESS --parseSolution=0 --ParRangeSizeInPerCent=1.0 --NpointsPerPar=3000 --maxNIter=40 --minNumberOfLines=40 --maxorderofpolynomial=4 --dampingFactor=0.85 --initialAcceptableMismatch=1.5 --nsigclip=2.5 --normalizeUncalibratedSpectrum=0 --normalizationBinSize=180 --LocalMaxFilterWidth=6 -v

echo ""
echo "Creating wavelength calibration product after stitching orders together: $WAVELENGTHPRODUCT"
echo ""
echo "$EXE/operaStitchOrders --outputWaveFile=$WAVELENGTHPRODUCT --inputSpectrum=$COMPEXTRACTEDSPECTRUM --inputWaveFile=$FIRSTWAVELENGTHPRODUCT --orderOfReference=37 --DWavelengthRange=0.1 --DWavelengthStep=0.00005 --XCorrelationThreshold=0.1 --sigmaThreshold=2.0 -v"
$EXE/operaStitchOrders --outputWaveFile=$WAVELENGTHPRODUCT --inputSpectrum=$COMPEXTRACTEDSPECTRUM --inputWaveFile=$FIRSTWAVELENGTHPRODUCT --orderOfReference=37 --DWavelengthRange=0.1 --DWavelengthStep=0.00005 --XCorrelationThreshold=0.1 --sigmaThreshold=2.0 -v
###############################

####### Create flat-field flux calibration spectrum ######
echo ""
echo "Creating flat-field flux calibration spectrum: $FLATFLUXCALIBRATIONSPECTRUM"
echo ""
echo "$EXE/operaCreateFlatFieldFluxCalibration --outputFluxCalibrationFile=$FLATFLUXCALIBRATIONSPECTRUM --inputMasterFlatSpectrum=$FLATEXTRACTEDSPECTRUM --wavelengthCalibration=$WAVELENGTHPRODUCT --binsize=500 --wavelengthForNormalization=$WAVELENGTHFORNORMALIZATION -v"
$EXE/operaCreateFlatFieldFluxCalibration --outputFluxCalibrationFile=$FLATFLUXCALIBRATIONSPECTRUM --inputMasterFlatSpectrum=$FLATEXTRACTEDSPECTRUM --wavelengthCalibration=$WAVELENGTHPRODUCT --binsize=500 --wavelengthForNormalization=$WAVELENGTHFORNORMALIZATION -v
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
set OBJECTLIST=$PRODUCTDIR$INSTCONFIGPREFIX"_object.list"
set STANDARDSLIST=$PRODUCTDIR$INSTCONFIGPREFIX"_standards.list"
################################################

###############################################
###############################################
######### START R E D U C T I O N #############
###############################################
###############################################
if ($EXECUTEREDUCTION == 1) then
###############################################
###### Print out parameters ######
echo "Running GRACES Pipeline 1.0: Reduction"
echo " "
echo "NIGHT = $NIGHT"
echo ""
##################################

echo "--------"
echo "STARTING REDUCTION"
echo ""

####### Create list of Object files ######
# input: directory and qualifiers
# This module creates a list of object files for reduction
echo "--------"
echo "Creating object list: $OBJECTLIST"
echo ""
if ($ALLOWANYREADOUTSPEED == 0) then
    echo "$EXE/operaQueryImageInfo --directory=$DATADIR -q "'"INSTRUME GSLICER EREADSPD OBSTYPE"'" INSTRUME=$INSTRUMENT GSLICER="'"'"$INSTRUMENTMODEKEY"'"'" EREADSPD="'"'"$READOUTSPEED"'"'" OBSTYPE=$OBJECTKEYWORD > $OBJECTLIST"
    $EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME GSLICER EREADSPD OBSTYPE" INSTRUME=$INSTRUMENT GSLICER="$INSTRUMENTMODEKEY" EREADSPD="$READOUTSPEED" OBSTYPE=$OBJECTKEYWORD > $OBJECTLIST
else
    echo "$EXE/operaQueryImageInfo --directory=$DATADIR -q "'"INSTRUME GSLICER OBSTYPE"'" INSTRUME=$INSTRUMENT GSLICER="'"'"$INSTRUMENTMODEKEY"'"'" OBSTYPE=$OBJECTKEYWORD > $OBJECTLIST"
    $EXE/operaQueryImageInfo --directory=$DATADIR -q "INSTRUME GSLICER OBSTYPE" INSTRUME=$INSTRUMENT GSLICER="$INSTRUMENTMODEKEY" OBSTYPE=$OBJECTKEYWORD > $OBJECTLIST
endif
##########################################

foreach OBJIMAGE (`cat $OBJECTLIST`)

    ####### Figure out image base name ######
    set DATADIRLEN=`echo $DATADIR | awk '{print length($1)}'`
    set OBJIMGBASENAME=`echo $OBJIMAGE | awk '{print substr($0,'$DATADIRLEN'+1,13)}'`
    ##########################################

    ####### Figure out info from object image header ######
    set OBJECTNAME=`$EXE/operagetheader --keyword=OBJECT $OBJIMAGE`
    set MJDATE=`$EXE/operagetheader --keyword=MJDATE $OBJIMAGE`
    set EXPTIME=`$EXE/operagetheader --keyword=EXPTIME $OBJIMAGE`
    echo "$EXPTIME"

    set absra_center=`$EXE/operagetheader --keyword=RA_DEG $OBJIMAGE`
    set absdec_center=`$EXE/operagetheader --keyword=DEC_DEG $OBJIMAGE`
    echo "Object=$OBJECTNAME JD=$MJDATE RA=$absra_center Dec=$absdec_center"
    ##########################################

    ####### Extract Object Spectrum ######
    set OBJECTSPECTRUM=$PRODUCTDIR$OBJIMGBASENAME".e.gz"
    echo "--------"
    echo "Extracting object spectrum product: $OBJECTSPECTRUM "
    echo ""
    echo "$EXE/operaExtraction --outputSpectraFile=$OBJECTSPECTRUM --inputImage=$OBJIMAGE --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --backgroundBinsize=300 --sigmaclip=6 --onTargetProfile=1 --starplusskymode=$STARPLUSKYMODEFLAG $INVERTSKYFIBERFLAG --usePolynomialFit=0 --removeBackground=0 --iterations=3 --maxthreads=4 -v"
    $EXE/operaExtraction --outputSpectraFile=$OBJECTSPECTRUM --inputImage=$OBJIMAGE --badpixelmask=$BADPIXELMASK --masterbias=$MASTERBIAS --masterflat=$MASTERFLAT --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum --inputInstrumentProfileFile=$INSTRUMENTPROFILEPRODUCT --inputGeometryFile=$GEOMETRYPRODUCT --inputApertureFile=$APERTUREPRODUCT --inputGainFile=$GAINPRODUCT --backgroundBinsize=300 --sigmaclip=6 --onTargetProfile=1 --starplusskymode=$STARPLUSKYMODEFLAG $INVERTSKYFIBERFLAG --usePolynomialFit=0 --removeBackground=0 --iterations=3 --maxthreads=4 -v

    ####### Calculate telluric wavelength correction ######
    set TELLWAVECAL=$PRODUCTDIR$OBJIMGBASENAME".tell.gz"
    echo "--------"
    echo "Calculating telluric wavelength correction product: $TELLWAVECAL "
    echo ""
#   E. Martioli 14 May 2014:
#   Telluric correction turned-off because it hasn't been tested for GRACES.
#   $EXE/operaTelluricWavelengthCorrection --outputWaveFile=$TELLWAVECAL --inputObjectSpectrum=$OBJECTSPECTRUM --inputWaveFile=$WAVELENGTHPRODUCT --telluric_lines=$TELLURICLINES --telluric_spectrum=$TELLURICSPECTRUM --spectralResolution=$SPECTRALRESOLUTION --initialWavelengthRange=0.1 --initialWavelengthStep=0.002 --XCorrelationThreshold=0.1 --subtractCentralWavelength=1 --normalizationBinsize=110 --sigmaThreshold=1.25 -v

    ####### Calculate barycentric wavelength correction ######
    set BARYWAVECAL=$PRODUCTDIR$OBJIMGBASENAME".rvel.gz"
    echo "--------"
    echo "Calculating barycentric wavelength correction product: $BARYWAVECAL "
    echo ""
#   E. Martioli 14 May 2014:
#   Barycentric correction not possible because target coordinates are not in the header yet
#   $EXE/operaBarycentricWavelengthCorrection --outputRVelFile=$BARYWAVECAL --inputWaveFile=$WAVELENGTHPRODUCT --observatory_coords="19:49:36 -155:28:18" --object_coords="$absra_center $absdec_center" --observatory_elevation=4207 --MJDTime=$MJDATE -v

    ####### Calculate final calibrated spectrum *.spc ######
    set CALIBRATEDSPECTRUM=$PRODUCTDIR$OBJIMGBASENAME".spc.gz"
    echo "--------"
    echo "Calculating final calibrated spectrum product: $CALIBRATEDSPECTRUM "
    echo ""
    echo "$EXE/$SPCMODULE --outputCalibratedSpectrum=$CALIBRATEDSPECTRUM --inputUncalibratedSpectrum=$OBJECTSPECTRUM --spectrumtype=17 --wavelengthCalibration=$WAVELENGTHPRODUCT --inputFlatFluxCalibration=$FLATFLUXCALIBRATIONSPECTRUM --inputWavelengthMaskForUncalContinuum=$ATYPEWAVELENGTHMASK --object="'"'"$OBJECTNAME"'"'" --numberOfPointsInUniformSample=150 --normalizationBinsize=750 --AbsoluteCalibration=0 --etime=1.0 -v $INVERTSKYFIBERFLAG"
    $EXE/$SPCMODULE --outputCalibratedSpectrum=$CALIBRATEDSPECTRUM --inputUncalibratedSpectrum=$OBJECTSPECTRUM --spectrumtype=17 --wavelengthCalibration=$WAVELENGTHPRODUCT --inputFlatFluxCalibration=$FLATFLUXCALIBRATIONSPECTRUM --inputWavelengthMaskForUncalContinuum=$ATYPEWAVELENGTHMASK --object="$OBJECTNAME" --numberOfPointsInUniformSample=150 --normalizationBinsize=750 --AbsoluteCalibration=0 --etime=1.0 -v $INVERTSKYFIBERFLAG

    ####### Generate LE formats ######
    set LESPCNW=$PRODUCTDIR$OBJIMGBASENAME".inw.s.gz"
    set LESPCU=$PRODUCTDIR$OBJIMGBASENAME".iu.s.gz"
    set LESPCUW=$PRODUCTDIR$OBJIMGBASENAME".iuw.s.gz"
    set LESPCN=$PRODUCTDIR$OBJIMGBASENAME".in.s.gz"
    echo "--------"
    echo "Calculating telluric wavelength correction: $OBJECTSPECTRUM "
    echo ""
#   E. Martioli 14 May 2014:
#   The modules below generate files in Libre-Esprit format.
#    $EXE/operaGenerateLEFormats --outputLEfilename=$LESPCNW --inputOperaSpectrum=$CALIBRATEDSPECTRUM --LibreEspritSpectrumType=$LESPECTRUMTYPE --object="$OBJECTNAME" --fluxType=2 --wavelengthType=3 -v
#    $EXE/operaGenerateLEFormats --outputLEfilename=$LESPCU --inputOperaSpectrum=$CALIBRATEDSPECTRUM --LibreEspritSpectrumType=$LESPECTRUMTYPE --object="$OBJECTNAME" --fluxType=3 --wavelengthType=3 -v
#    $EXE/operaGenerateLEFormats --outputLEfilename=$LESPCUW --inputOperaSpectrum=$CALIBRATEDSPECTRUM --LibreEspritSpectrumType=$LESPECTRUMTYPE --object="$OBJECTNAME" --fluxType=3 --wavelengthType=4 -v
#    $EXE/operaGenerateLEFormats --outputLEfilename=$LESPCN --inputOperaSpectrum=$CALIBRATEDSPECTRUM --LibreEspritSpectrumType=$LESPECTRUMTYPE --object="$OBJECTNAME" --fluxType=2 --wavelengthType=3 -v
    ###############################################

#   E. Martioli 14 May 2014:
#   Flux calibration creation is crashing. I haven't looked into this problem yet.
#   The call below should be ok, but it requires a list of standards done by hand.
#   The list of standards consist of a list of objects with flux calibrated spectra in
#   the pipeline config directory. Both the list and data files only work if written in
#   a particular format.
#
#   foreach STDIMAGE (`cat $STANDARDSLIST`)
#        if ($STDIMAGE == $OBJIMAGE) then
            ####### Generate flux calibration from standard ######
#            set FLUXCALIBRATIONPRODUCT=$OBJIMGBASENAME".fcal.gz"
#            set STDCALIBRATEDSPECTRUM=$CONFIGDIR/standardstars/"feige66_operaFluxCal.dat"
#            echo "--------"
#            echo "Generating flux calibration product: $FLUXCALIBRATIONPRODUCT "
#           echo ""
#            echo "$EXE/operaCreateFluxCalibration --inputCalibratedSpectrum=$STDCALIBRATEDSPECTRUM --inputUncalibratedSpectrum=$OBJECTSPECTRUM --inputFlatFluxCalibration=$FLATFLUXCALIBRATIONSPECTRUM --outputFluxCalibrationFile=$FLUXCALIBRATIONPRODUCT --inputWavelengthMaskForRefContinuum=$ATLASWAVELENGTHMASK --inputWavelengthMaskForUncalContinuum=$ATYPEWAVELENGTHMASK --inputWaveFile=$WAVELENGTHPRODUCT --inputApertureFile=$APERTUREPRODUCT --wavelengthForNormalization=$WAVELENGTHFORNORMALIZATION --exposureTime=$EXPTIME --numberOfPointsInUniformSample=200 --numberOfPointsInUniformRefSample=70 --binsize=100 -v"
#           $EXE/operaCreateFluxCalibration --inputCalibratedSpectrum=$STDCALIBRATEDSPECTRUM --inputUncalibratedSpectrum=$OBJECTSPECTRUM --inputFlatFluxCalibration=$FLATFLUXCALIBRATIONSPECTRUM --outputFluxCalibrationFile=$FLUXCALIBRATIONPRODUCT --inputWavelengthMaskForRefContinuum=$ATLASWAVELENGTHMASK --inputWavelengthMaskForUncalContinuum=$ATYPEWAVELENGTHMASK --inputWaveFile=$WAVELENGTHPRODUCT --inputApertureFile=$APERTUREPRODUCT --wavelengthForNormalization=$WAVELENGTHFORNORMALIZATION --exposureTime=$EXPTIME --numberOfPointsInUniformSample=100 --numberOfPointsInUniformRefSample=150 --binsize=200 -v
            ###############################################
#        endif
#    end

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


