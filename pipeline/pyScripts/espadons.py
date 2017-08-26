# -*- coding: utf-8 -*-
"""
Created on Oct 21 2014
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on Mar 8, 2017
"""
import sys, os
import subprocess
from time import strftime

pylibdir = os.path.split(os.path.dirname(__file__))[0] + '/pyLib'
sys.path.insert(0,pylibdir)

from operacommands import *

########## DIRECTORIES ############
class Directories :
     
    'Common base class for directories'
    def __init__(self, pipelinehomedir,datarootdir,productrootdir,night):
        self.OP_HOME=pipelinehomedir
        self.DATAROOTDIR=datarootdir
        self.PRODUCTROOTDIR=productrootdir
        self.PIPELINEDIR=self.OP_HOME+"/pipeline/"
        self.DEFAULTCALIBRATIONDIR=self.OP_HOME+"/DefaultCalibration/"
        self.CONFIGDIR=self.OP_HOME+"/config/"
        self.STANDARDSDIR=self.CONFIGDIR+"/standardstars/"
        self.EXE=self.OP_HOME+"/bin/"
        self.DATADIR=self.DATAROOTDIR+"/"+night+"/"
        self.PRODUCTDIR=self.PRODUCTROOTDIR+"/"+night+"/"
        
        self.createNightDir(self.PRODUCTDIR)

    def createNightDir(self, inputdir):	
        if os.path.exists(inputdir) :
            print "Product directory: ",inputdir," (already exists)"
        else :
            print "Product directory: ",inputdir," (created new directory)"
            createdircommand = "mkdir " + inputdir
            os.system(createdircommand)

##################################

######### CONFIG FILES ###########
class ConfigFiles :
    'Common base class for config files'
    standarddata = {}
    
    def __init__(self, Dirs):
        self.BADPIXELMASK = Dirs.CONFIGDIR + "badpix_olapa-a.fits.fz"
        self.THARATLASLINES = Dirs.CONFIGDIR + "thar_MM201006.dat.gz"
        self.THARATLASSPECTRUM = Dirs.CONFIGDIR + "LovisPepe_ThArAtlas.dat.gz"
        self.WAVEFIRSTGUESS = Dirs.CONFIGDIR + "wcal_ref.dat.gz"
        self.SOLARTYPEWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForUncalContinuumDetection_SolarTypeStars.txt"
        #self.ATYPEWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForUncalContinuumDetection_SolarTypeStars.txt"
        self.ATYPEWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForUncalContinuumDetection.txt"
        self.ATLASWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForRefContinuumDetection.txt"
        self.TELLURICWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForTelluricAbsorption.txt"
        self.TELLURICLINES = Dirs.CONFIGDIR + "opera_HITRAN08-extracted.par.gz"
        #self.TELLURICLINES = Dirs.CONFIGDIR + "skyline_skycal-Eso_mod.txt"
        self.TELLURICSPECTRUM = Dirs.CONFIGDIR + "KPNO_atmtrans.dat.gz"
        self.LEORDERWAVELENGTH = Dirs.CONFIGDIR + "LE_order_wavelength.dat"

        self.STANDARDLISTFILE = Dirs.STANDARDSDIR + "operaStandardStars.dat"
        self.readStandards(Dirs,self.STANDARDLISTFILE)
    
        self.SYNTHETICSPECTRUM = Dirs.CONFIGDIR + "/syntethicSpectra/Teff5500.spec"
        self.RVXCORRWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForRVxcorr.txt"
        
        #self.OLAPAFLATRESPONSE = Dirs.CONFIGDIR + "flat_resp_olapa.s"
        self.OLAPAFLATRESPONSE = Dirs.CONFIGDIR + "flat_resp_olapa.fits.gz"
        self.EEV1FLATRESPONSE = Dirs.CONFIGDIR + "flat_resp_eev1.s"

        self.SOURCELINES = Dirs.CONFIGDIR + "sourcelinesForRVMeasurements.txt"
        #self.SOURCELINES = Dirs.CONFIGDIR + "sourcelinesForRVMeasurements_SolarType.txt"
        self.INTERSTELLARLINES = Dirs.CONFIGDIR + "/spectralLinesLibraries/interstellarLinesForRVMeasurements.txt"
        self.WAVELENGTHRANGESFORRVMEASUREMENTS = Dirs.CONFIGDIR + "/spectralLinesLibraries/wavelengthRangesForRVMeasurements.txt"
        
        self.SOMEREFERENCESPECTRUM = "/Users/edermartioli/Reductions/GRACES/20150807/KIC09472174_20150807_master.norm.spc"

        self.LEAPSECONDSFILE = Dirs.CONFIGDIR + "/leapseconds.dat"
    
        self.TIMESTAMPFILE =Dirs.OP_HOME + "/TIMESTAMP"
        self.OPERAVERSION = self.readfirstline(self.TIMESTAMPFILE)
    
    def readfirstline(self, file) :
        with open(file, 'r') as f:
            first_line = f.readline()
            return first_line.rstrip()

# Function below reads list of standard stars for which there is available calibration data
    def readStandards(self, Dirs, stdListFilename) :
        stdlistfile = open(stdListFilename)
        for line in stdlistfile :
            if(line[0] != '#') :
                name = line.split(' ')[0]
                datafilename = Dirs.STANDARDSDIR + name + "_operaFluxCal.dat"
                ra = line.split(' ')[1]
                dec = line.split(' ')[2]
                pmra = line.split(' ')[3]
                pmdec = line.split(' ')[4]
                vmag = line.split(' ')[5]
                teff = line.split(' ')[6]
                rv = line.split(' ')[7]
                stype = line.split(' ')[8]
                if os.path.exists(datafilename) :
                    self.standarddata[name] = [datafilename,ra,dec,pmra,pmdec,vmag,teff,rv,stype]
        
        stdlistfile.close()

# Function to check if a given standard star exist
    def hasStandard(self, StandardName) :
        if StandardName in self.standarddata :
            return True
        else :
            return False

# Function below returns the corresponding data for a standard star in the list
    def getStandardDataFile(self, key) :
        return self.standarddata.get(key, 0)[0]

    def getStandardRA(self, key) :
        return self.standarddata.get(key, 0)[1]

    def getStandardDec(self, key) :
        return self.standarddata.get(key, 0)[2]

    def getStandardPMRA(self, key) :
        return self.standarddata.get(key, 0)[3]

    def getStandardPMDec(self, key) :
        return self.standarddata.get(key, 0)[4]
            
    def getStandardVmag(self, key) :
        return self.standarddata.get(key, 0)[5]
            
    def getStandardTeff(self, key) :
        return self.standarddata.get(key, 0)[6]

    def getStandardRV(self, key) :
        return self.standarddata.get(key, 0)[7]
    
    def getStandardSType(self, key) :
        return self.standarddata.get(key, 0)[8]
#-------------------------------------------
##################################

######### KEYWORDS ###########
class Keywords :
    'Common base class for keywords'
    def __init__(self):
        self.BIASKEYWORD="BIAS"
        self.FLATKEYWORD="FLAT"
        self.COMPKEYWORD="COMPARISON"
        self.ALIGNKEYWORD="ALIGN"
        self.OBJECTKEYWORD="OBJECT"

        self.INSTRUMEKEY = 'INSTRUME'
        self.OBSTYPEKEY = 'OBSTYPE'
        self.READMODEKEY = 'EREADSPD'
        self.INSTMODEKEY = 'INSTMODE'
##############################

######### DEFAULT CALIBRATION ###########
class DefaultCalibration :    
    'Common base class for default calibration'
    def __init__(self, Dirs, Instmode, Readmode):
        self.DEFAULTMASTERBIAS = Dirs.DEFAULTCALIBRATIONDIR + Readmode.READOUTSPEEDSHORTNAME + "_masterbias.fits.gz"
        self.DEFAULTGAIN = Dirs.DEFAULTCALIBRATIONDIR + Readmode.READOUTSPEEDSHORTNAME + ".gain.gz"
        self.DEFAULTMASTERALIGN = Dirs.DEFAULTCALIBRATIONDIR + Instmode.INSTRUMENTMODESHORTNAME + "_masteralign.fits.gz"
        self.DEFAULTMASTERCOMP = Dirs.DEFAULTCALIBRATIONDIR + Instmode.INSTRUMENTMODESHORTNAME + "_mastercomp.fits.gz"
        self.DEFAULTMASTERFLAT = Dirs.DEFAULTCALIBRATIONDIR + Instmode.INSTRUMENTMODESHORTNAME + "_masterflat.fits.gz"
        self.DEFAULTMASTERFLUXCAL = Dirs.DEFAULTCALIBRATIONDIR + Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME + ".fcal.gz"

##################################


######### INSTRUMENT MODE ###########
class InstMode :
    'Common base class for all instrument modes'
    mode = 0
    def __init__(self, intrumentmode):
        self.mode = int(intrumentmode)

        self.INSTRUME="ESPaDOnS"

        self.observatory_coords = '19:49:36 -155:28:18'
        self.observatory_elevation = '4207'

        self.MINORDERTOEXTRACT = '22'
        self.MAXORDERTOEXTRACT = '61'

        #self.MINORDERTOEXTRACT = '36'
        #self.MAXORDERTOEXTRACT = '59'

        self.MASTERCOMP_COMBINEMETHOD='2'
        self.MASTERCOMP_SATURATIONLIMIT='65535'
        self.MASTERCOMP_OUTPUTEXPOSURETIME='40'
        self.MASTERCOMP_TRUNCATEOUTPUTFLUXTOSATURATION='1'
        self.MASTERCOMP_BIASCONSTANT='0'
        self.MASTERCOMP_EXPTIMEFITSKEYWORD='EXPTIME'
        
        self.GAIN_DATASEC='1 2048 1 4608'
        self.GAIN_DSECA='21 1044 1 4608'
        self.GAIN_DSECB='1045 2068 1 4608'
        self.GAIN_SUBWINDOW='1 1000 1 2000'
        self.GAIN_MINPIXPERBIN='1000'
        self.GAIN_MAXNBINS='100'
        self.GAIN_LOWESTCOUNT='1000'
        self.GAIN_HIGHESTCOUNT='30000'
        
        self.SPACING_NUMBEROFSAMPLES = '30'
        self.SPACING_SAMPECENTERPOSITION = '2300'
        self.SPACING_DETECTIONMETHOD = '2'
        self.SPACING_SUBFORMAT = '300 2040 3 4600'
        
        self.GEOM_SUBFORMAT='21 2040 3 4100'
        self.GEOM_DETECTIONMETHOD = '2'
        self.GEOM_FFTFILTER = '0'
        self.GEOM_NSAMPLES = '5'
        self.GEOM_ORDEROFTRACINGPOLYNOMIAL = '3'
        self.GEOM_BINSIZE = '25'
        self.GEOM_COLDISPERSION = '1'
        self.GEOM_INVERTORDERS = '1'
        self.GEOM_REFERENCEORDERSAMPLEPOSITION = '2300'
        self.GEOM_GRACES = '1'
        
        self.PROF_BINSIZE = '100'
        self.PROF_SPECTRALELEMENTHEIGHT = '1.0'
        self.PROF_MAXTHREADS = '4'
        self.PROF_MINIMUMLINES = '3'
        self.PROF_LOCALMAXFILTERWIDTH = '3.0'
        self.PROF_DETECTIONTHRESHOLD = '0.03'
        self.PROF_MINPEAKDEPTH = '0.2'
        self.PROF_IPMETHOD='2'
        
        self.APER_APERTUREHEIGHT = '0.6923'
        self.APER_PICKIMAGEROW = '0'
        self.APER_NROWSAMPLES = '10'
        self.APER_XBIN = '8'
        
        self.CALEXTRACTION_SPECTRUMTYPE='5'
        self.CALEXTRACTION_SPECTRUMTYPENAME='RawBeamSpectrum'
        self.CALEXTRACTION_BACKGROUNDBINSIZE='300'
        self.CALEXTRACTION_REMOVEBACKGROUND='0'
        self.CALEXTRACTION_ITERATIONS='3'
        self.CALEXTRACTION_ONTARGETPROFILE='1'
        self.CALEXTRACTION_USEPOLYNOMIALFIT='0'
        self.CALEXTRACTION_REJECTBADPIXINRAWEXTRACTION='1'
        self.CALEXTRACTION_MINSIGMACLIP='15'
        self.CALEXTRACTION_SIGMACLIPRANGE='6'
        
        self.WAVE_PARSESOLUTION='0'
        self.WAVE_PARRANGESIZEINPERCENT='2'
        self.WAVE_NPOINTSPERPAR='1000'
        self.WAVE_MAXNITER='40'
        self.WAVE_MINNUMBEROFLINES='40'
        self.WAVE_MAXORDEROFPOLYNOMIAL='5'
        self.WAVE_DAMPINGFACTOR='0.9'
        self.WAVE_INITIALACCEPTABLEMISMATCH='1.5'
        self.WAVE_NSIGCLIP='2.25'
        self.WAVE_NORMALIZEUNCALIBRATEDSPECTRUM='0'
        self.WAVE_NORMALIZATIONBINSIZE='100'
        self.WAVE_LOCALMAXFILTERWIDTH='5'
        self.WAVE_DETECTIONTHRESHOLD='0.05'
        self.WAVE_MINPEAKDEPTH='1.0'
        self.WAVE_NORDERSTOSEARCHAROUND = '3'
        self.WAVE_REFERENCEORDER='50'
        
        self.FLATFLUXCAL_BINSIZE='500'
        
        self.EXTRACTION_MAXTHREADS='4'
        self.EXTRACTION_SPECTRUMTYPE='7'
        self.EXTRACTION_SPECTRUMTYPENAME='OptimalBeamSpectrum'
        self.EXTRACTION_BACKGROUNDBINSIZE='300'
        self.EXTRACTION_REMOVEBACKGROUND='0'
        self.EXTRACTION_ITERATIONS='2'
        self.EXTRACTION_ONTARGETPROFILE='1'
        self.EXTRACTION_USEPOLYNOMIALFIT='0'
        self.EXTRACTION_REJECTBADPIXINRAWEXTRACTION='0'
        self.EXTRACTION_MINSIGMACLIP='25'
        self.EXTRACTION_SIGMACLIPRANGE='6'
        
        self.TELL_XCORRELATIONTHRESHOLD='0.1'
        self.TELL_NORMALIZATIONBINSIZE='110'
        self.TELL_RVCORRECTIONMETHOD='1'
        self.TELL_LOCALMAXFILTERWIDTH='3.0'
        self.TELL_MINPEAKDEPTH='1.0'
        self.TELL_DETECTIONTHRESHOLD='0.01'
        self.TELL_NSIGCLIP='3.0'
        self.TELL_MINNUMBEROFMATCHEDLINES='10'
        self.TELL_DUPLICATELINETHRESHOLD='0.001'
        self.TELL_RADIALVELOCITYRANGE='0.8'
        self.TELL_RADIALVELOCITYSTEP='0.05'
        self.TELL_INVERTSKYFIBERFLAG=''

        self.SNR_SPECTRUMTYPE='31'
        self.SNR_CENTRALSNR='1'
        self.SNR_SPECTRALBINSIZE=self.APER_APERTUREHEIGHT
        
        self.RV2_NORMALIZATIONBINSIZE='400'
        self.RV2_INITIALRVGUESS='0.0'
        self.RV2_SOURCELINERESOLUTION = '2500'
        
        self.FCAL_NUMBEROFPOINTSINUNIFORMSAMPLE='150'
        self.FCAL_NUMBEROFPOINTSINUNIFORMREFSAMPLE='70'
        self.FCAL_BINSIZE='500'
        self.FCAL_WAVELENGTHFORNORMALIZATION='548'
        
        self.FLATRESP_NUMBEROFPOINTSINUNIFORMSAMPLE='300'
        self.FLATRESP_NUMBEROFPOINTSINUNIFORMREFSAMPLE='70'
        self.FLATRESP_BINSIZE='750'
        
        self.SPC_SPECTRUMTYPE='17'
        self.SPC_NUMBEROFPOINTSINUNIFORMSAMPLE='150'
        self.SPC_NORMALIZATIONBINSIZE='110'
        self.SPC_ABSOLUTECALIBRATION='0'
        
        self.STACKSPC_RVBIN = '1.6'    # spectral bin in km/s
        self.STACKSPC_MINWL = '0.0'    # min wavelength in nm (0 for full range)
        self.STACKSPC_MAXWL = '0.0'    # max wavelength in nm (0 for full range)
        self.STACKSPC_SNRCLIP = '10'    # Min SNR to accept a given spectral point
        self.STACKSPC_NCUTATORDERENDS = '100'    # Number of points to exclude at order edges
        self.STACKSPC_APPLYTELL = '1'
        self.STACKSPC_APPLYHELIO = '1'

        if (self.mode == 1) :
            
            self.INSTRUMENTMODESHORTNAME='StarOnly'
            self.INSTRUMENTMODEKEY='Spectroscopy, star only, R=80,000'
            self.STARPLUSKYMODEFLAG=0
            
            self.EXTRACTION_INVERTSKYFIBERFLAG=''
            
            self.SPC_MODULENAME='operaStarOnly'
            self.SPC_SKYOVERSTARFIBERAREARATIO=' --SkyOverStarFiberAreaRatio=1.0'
            self.SPC_INVERTSKYFIBERFLAG=''

            self.SPECTRALRESOLUTION='80000'
            self.LESPECTRUMTYPE='21'
            self.OUTPUTFITSSPECTRUMTYPE=21

            self.SPACING_APERTURE='30'
            self.SPACING_REFERENCEORDERNUMBER='48'
            self.SPACING_REFERENCEORDERSEPARATION='57.5'
            
            self.GEOM_APERTURE='30'
            self.GEOM_MAXNORDERS='40'
            self.GEOM_MINORDERTOUSE='22'
            self.GEOM_RECENTERIPUSINGSLICESYMMETRY='0'
            self.GEOM_NUMBEROFSLICES='6'
            
            self.PROF_REFERENCELINEWIDTH='2.0'
            self.PROF_IPXSIZE='30'
            self.PROF_IPYSIZE='4'
            self.PROF_IPXSAMPLING='3'
            self.PROF_IPYSAMPLING='5'
            self.PROF_TILTANGLE='-2.0'
            
            self.APER_NUMBEROFBEAMS='1'
            self.APER_CONSTANTTILTFLAG='0'
            self.APER_APERTURE='30'
            self.APER_SKYAPERAPERTURE='1'
            self.APER_GAP='0'
    
            self.WAVE_UNCALLINEWIDTH='1.6'
            
            self.SPC_MODULENAME='operaStarOnly'
            self.SPC_SKYOVERSTARFIBERAREARATIO=''
            self.SPC_INVERTSKYFIBERFLAG=''
            
            self.RV_RADIALVELOCITYRANGE='0.8'
            self.RV_RADIALVELOCITYSTEP='0.05'
            self.RV_THRESHOLD='0.05'
            
            self.RADIALVELOCITYRANGE=0.8
            self.RADIALVELOCITYSTEP=0.05
            self.RADIALVELOCITYSEARCHRANGE=200
            self.RADIALVELOCITYSEARCHSTEP=0.5

        elif (self.mode == 2) :
            
            self.INSTRUMENTMODESHORTNAME='StarPlusSky'
            self.INSTRUMENTMODEKEY='Spectroscopy, star+sky, R=65,000'
            self.STARPLUSKYMODEFLAG=1
            
            self.EXTRACTION_INVERTSKYFIBERFLAG=''
            
            self.SPECTRALRESOLUTION='65000'
            self.LESPECTRUMTYPE='20'
            self.OUTPUTFITSSPECTRUMTYPE=20
            
            self.SPACING_APERTURE='32'
            self.SPACING_REFERENCEORDERNUMBER='48'
            self.SPACING_REFERENCEORDERSEPARATION='57.5'
            
            self.GEOM_APERTURE='32'
            self.GEOM_MAXNORDERS='40'
            self.GEOM_MINORDERTOUSE='22'
            self.GEOM_RECENTERIPUSINGSLICESYMMETRY='1'
            self.GEOM_NUMBEROFSLICES='3'
            
            self.PROF_REFERENCELINEWIDTH='2.0'
            self.PROF_IPXSIZE='32'
            self.PROF_IPYSIZE='4'
            self.PROF_IPXSAMPLING='3'
            self.PROF_IPYSAMPLING='5'
            self.PROF_TILTANGLE='-1.0'
            
            self.APER_NUMBEROFBEAMS='2'
            self.APER_CONSTANTTILTFLAG='1'
            self.APER_APERTURE='32'
            self.APER_SKYAPERAPERTURE='1'
            self.APER_GAP='0'
            
            self.WAVE_UNCALLINEWIDTH='1.8'
            
            self.SPC_MODULENAME='operaStarPlusSky'
            self.SPC_SKYOVERSTARFIBERAREARATIO=' --SkyOverStarFiberAreaRatio=1.890625'
            self.SPC_INVERTSKYFIBERFLAG=' --starplusskyInvertSkyFiber=0'
            
            self.RV_RADIALVELOCITYRANGE='0.8'
            self.RV_RADIALVELOCITYSTEP='0.05'
            self.RV_THRESHOLD='0.05'
            
            self.RADIALVELOCITYRANGE=0.8
            self.RADIALVELOCITYSTEP=0.05
            self.RADIALVELOCITYSEARCHRANGE=200
            self.RADIALVELOCITYSEARCHSTEP=0.5

        elif (self.mode == 3) :
            
            self.INSTRUMENTMODESHORTNAME='Polar'
            self.INSTRUMENTMODEKEY='Polarimetry, R=65,000'
            self.STARPLUSKYMODEFLAG=0
            
            self.EXTRACTION_INVERTSKYFIBERFLAG=''
            
            self.SPECTRALRESOLUTION='65000'
            self.LESPECTRUMTYPE='22'
            self.OUTPUTFITSSPECTRUMTYPE=22
            
            self.SPACING_APERTURE='30'
            self.SPACING_REFERENCEORDERNUMBER='48'
            self.SPACING_REFERENCEORDERSEPARATION='57.2'
            
            self.GEOM_APERTURE='30'
            self.GEOM_MAXNORDERS='40'
            self.GEOM_MINORDERTOUSE='22'
            self.GEOM_RECENTERIPUSINGSLICESYMMETRY='1'
            self.GEOM_NUMBEROFSLICES='3'
            
            self.PROF_REFERENCELINEWIDTH='2.0'
            self.PROF_IPXSIZE='32'
            self.PROF_IPYSIZE='4'
            self.PROF_IPXSAMPLING='3'
            self.PROF_IPYSAMPLING='5'
            self.PROF_TILTANGLE='-2.0'
            
            self.APER_NUMBEROFBEAMS='2'
            self.APER_CONSTANTTILTFLAG='0'
            self.APER_APERTURE='32'
            self.APER_SKYAPERAPERTURE='1'
            self.APER_GAP='0'

            self.WAVE_UNCALLINEWIDTH='1.8'
            
            self.SPC_MODULENAME='operaPolarIntensity'
            self.SPC_SKYOVERSTARFIBERAREARATIO=''
            self.SPC_INVERTSKYFIBERFLAG=''
            
            self.RV_RADIALVELOCITYRANGE='0.8'
            self.RV_RADIALVELOCITYSTEP='0.05'
            self.RV_THRESHOLD='0.05'
            
            self.RADIALVELOCITYRANGE=0.8
            self.RADIALVELOCITYSTEP=0.05
            self.RADIALVELOCITYSEARCHRANGE=200
            self.RADIALVELOCITYSEARCHSTEP=0.5

##################################

######### READOUT MODE ###########
class ReadoutMode :
    'Common base class for all readout modes'
    mode = 0
    def __init__(self, readoutmode):
        self.mode = int(readoutmode)
        if (self.mode == 1) :
            self.READOUTSPEEDSHORTNAME="fast"
            self.READOUTSPEED="Fast: 4.70e noise, 1.60e/ADU, 32s"
            self.DEFAULTGAIN="1.6"
            self.DEFAULTNOISE="4.14"
            self.GAIN_NUMBEROFAMPLIFIERS='1'
        elif (self.mode == 2) :
            self.READOUTSPEEDSHORTNAME="normal"
            self.READOUTSPEED="Normal: 4.20e noise, 1.30e/ADU, 38s"
            self.DEFAULTGAIN="1.3"
            self.DEFAULTNOISE="3.8"
            self.GAIN_NUMBEROFAMPLIFIERS='1'
        elif (self.mode == 3) :
            self.READOUTSPEEDSHORTNAME="slow"
            self.READOUTSPEED="Slow: 2.90e noise, 1.20e/ADU, 60s"
            self.DEFAULTGAIN="1.1"
            self.DEFAULTNOISE="2.98"
            self.GAIN_NUMBEROFAMPLIFIERS='1'

        if (self.mode == 4) :
            self.READOUTSPEEDSHORTNAME="fast_2amp"
            self.READOUTSPEED="Fast: ??e noise, ??e/ADU, ??s"
            self.DEFAULTGAIN="1.6"
            self.DEFAULTNOISE="4.14"
            self.GAIN_NUMBEROFAMPLIFIERS='2'
        elif (self.mode == 5) :
            self.READOUTSPEEDSHORTNAME="normal_2amp"
            self.READOUTSPEED="Normal: 4.15e noise, 1.30e/ADU, 19s"
            self.DEFAULTGAIN="1.30"
            self.DEFAULTNOISE="4.15"
            self.GAIN_NUMBEROFAMPLIFIERS='2'
        elif (self.mode == 6) :
            self.READOUTSPEEDSHORTNAME="slow_2amp"
            self.READOUTSPEED="Slow: 2.90e noise, 1.15e/ADU, 30s"
            self.DEFAULTGAIN="1.15"
            self.DEFAULTNOISE="2.90"
            self.GAIN_NUMBEROFAMPLIFIERS='2'
##################################

######### OBJECT PRODUCT FILE NAMES, DEPENDENCIES, AND  COMMAND LINES ###########
def setObjectProducts(products, Dirs, night, Instmode, Readmode, DefaultCal, Keywords, config, verbose) :
    
    emptystring = ""
    
    verstr = ""
    if (verbose) :
        verstr = " --verbose"
    
    objproducts = {}
    objdependencies = {}
    objcommands = {}
    objcommandsType = {}
    
    commandline = ObjectListCommand(Dirs, Instmode, Readmode, Keywords.OBJECTKEYWORD, Keywords)
    objectfiles = subprocess.check_output(commandline,stderr=subprocess.STDOUT,shell=True).rstrip('\n') 
    objfilelist = sorted(objectfiles.split())
    
    polquad = []
    polextkey = []
    polcount = 0
    pkey = ""
    ptarget = ""

    listofobjects = []
    listofobjectstackkeys = []
    listofobjecttargets = []
    
    listofobjectspcs = {}
    listofobjectspcsDependencies = {}
    
    listofobjectrvs = {}
    listofobjectrvsDependencies = {}
    
    listofstdfcalkey = []
    listofstdfcaltargets = []

    masterfcalkeytarget = config.OLAPAFLATRESPONSE
    inputMaskForSPCModule = ""
    masterfcalkey = "MASTERFLATRESP"+Instmode.INSTRUMENTMODESHORTNAME

    for file in objfilelist:
        
        base=os.path.basename(file)
        
        if ".fits.gz" in base :
            basename = os.path.splitext(os.path.splitext(base)[0])[0]
        elif ".fits" in file :
            basename = os.path.splitext(base)[0]
        elif ".fit" in file :
            basename = os.path.splitext(base)[0]
        else :
            print "Error: unknown extension, exiting ..."
            exit()
        
        ### Get header data info #########
        objectcommand = Dirs.EXE +"operagetheader --keyword=OBJECT " + file
        objectname = subprocess.check_output(objectcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        objectNameWithoutWhites = objectname.replace(" ", "")
        
        currentTime = strftime("%c")

        if verbose :
            print "Object filename: ",file, ", basename = ", basename, ", object=",objectNameWithoutWhites

        if objectNameWithoutWhites not in listofobjects :
            listofobjects.append(objectNameWithoutWhites)
            listofobjectspcs[objectNameWithoutWhites] = []
            listofobjectrvs[objectNameWithoutWhites] = []
            listofobjectspcsDependencies[objectNameWithoutWhites] = []
            listofobjectrvsDependencies[objectNameWithoutWhites] = []

        exptimecommand = Dirs.EXE +"operagetheader --keyword=EXPTIME " + file
        exptime = subprocess.check_output(exptimecommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        racommand = Dirs.EXE +"operagetheader --keyword=RA_DEG " + file
        absra_center = subprocess.check_output(racommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        deccommand = Dirs.EXE +"operagetheader --keyword=DEC_DEG " + file
        absdec_center = subprocess.check_output(deccommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        mjcommand = Dirs.EXE +"operagetheader --keyword=MJDATE " + file
        mjdate = subprocess.check_output(mjcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        airmasscommand = Dirs.EXE +"operagetheader --keyword=AIRMASS " + file
        airmass = subprocess.check_output(airmasscommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        expnumcommand = Dirs.EXE +"operagetheader --keyword=EXPNUM " + file
        expnum = subprocess.check_output(expnumcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        outsideTemperaturecommand = Dirs.EXE +"operagetheader --keyword=TEMPERAT " + file
        outsideTemperature = subprocess.check_output(outsideTemperaturecommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        windspeedcommand = Dirs.EXE +"operagetheader --keyword=WINDSPED " + file
        windspeed = subprocess.check_output(windspeedcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        relativeHumididycommand = Dirs.EXE +"operagetheader --keyword=RELHUMID " + file
        relativeHumididy = subprocess.check_output(relativeHumididycommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        atmoPressurecommand = Dirs.EXE +"operagetheader --keyword=PRESSURE " + file
        atmoPressure = subprocess.check_output(atmoPressurecommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        startHAcommand = Dirs.EXE +"operagetheader --keyword=HA " + file
        startHA = subprocess.check_output(startHAcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        ########################################

        ### Create TARGETS for Extraction #########
        extkey = "EXTRACT" + basename
        exttarget = Dirs.PRODUCTDIR + basename + ".e.gz"
        objproducts[extkey] = exttarget
        objdependencies[extkey] = ["MASTERBIAS","MASTERFLAT","GAINPRODUCT","GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","APERTUREPRODUCT"]
        objcommands[extkey] = objectExtractionCommand(Dirs, exttarget, file, products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["GAINPRODUCT"], products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["APERTUREPRODUCT"], Instmode) + verstr
        objcommandsType[extkey] = True
        ########################################

        ### Create TARGETS for telluric wavelength correction #########
        tellkey = "TELLWAVE" + basename
        telltarget = Dirs.PRODUCTDIR + basename + ".tell.gz"
        objproducts[tellkey] = telltarget
        objdependencies[tellkey] = [extkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
        objcommands[tellkey] = TelluricWaveCommand(Dirs, telltarget, exttarget, products["WAVELENGTHPRODUCT"], products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode, False) + verstr
        objcommandsType[tellkey] = True
        ###############################################################
    
        ### Create TARGETS for heliocentric wavelength correction ######
        rvelkey = "HELIORVEL" + basename
        rveltarget = Dirs.PRODUCTDIR + basename + ".rvel.gz"
        objproducts[rvelkey] = rveltarget
        objdependencies[rvelkey] = []
        objcommands[rvelkey] = HeliocentricWaveCommand(Dirs, rveltarget, absra_center, absdec_center, mjdate, exptime, startHA, config, Instmode) + verstr
        objcommandsType[rvelkey] = True
        ###############################################################
    
        ###############################################################
        radvel2key = "RADIALVELOCITY2" + basename
        radvel2target = Dirs.PRODUCTDIR + basename + ".rv2.gz"
        objproducts[radvel2key] = radvel2target
        objdependencies[radvel2key] = [extkey,"WAVELENGTHPRODUCT",rvelkey, tellkey,"FLATFLUXCALIBRATIONSPECTRUM"]
        objcommands[radvel2key] = RadialVelocity2Command(Dirs, radvel2target, exttarget, products["WAVELENGTHPRODUCT"], rveltarget,telltarget,"", config, Instmode, mjdate , False) + verstr
        objcommandsType[radvel2key] = True
        ###############################################################

        ###############################################################
        snrkey = "SNR" + basename
        snrtarget = Dirs.PRODUCTDIR + basename + ".sn.gz"
        objproducts[snrkey] = snrtarget
        objdependencies[snrkey] = [extkey,"WAVELENGTHPRODUCT"]
        objcommands[snrkey] = GenerateSNR(Dirs, exttarget, snrtarget, products["WAVELENGTHPRODUCT"], objectNameWithoutWhites, Instmode)
        objcommandsType[snrkey] = True
        ###############################################################

        ### Create TARGETS for calibrated spectrum *.spc ##############
        spckey = "OPSPC" + basename
        spctarget = Dirs.PRODUCTDIR + basename + ".spc.gz"
        objproducts[spckey] = spctarget
        objdependencies[spckey] = [extkey,tellkey,rvelkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM",masterfcalkey]
        #objcommands[spckey] = SpcModuleCommand(Dirs, spctarget, Instmode, config, exttarget, products["FLATFLUXCALIBRATIONSPECTRUM"],objproducts["MASTERFLUXCALIBRATION"], rveltarget, telltarget, products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
        objcommands[spckey] = SpcModuleCommand(Dirs, spctarget, Instmode, config, exttarget, products["FLATFLUXCALIBRATIONSPECTRUM"],emptystring, rveltarget, telltarget,products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
        objcommandsType[spckey] = True

        if objectNameWithoutWhites in listofobjects :
            listofobjectspcs[objectNameWithoutWhites].append(spctarget)
            listofobjectspcsDependencies[objectNameWithoutWhites].append(spckey)
            # -- E Martioli Feb 05 2016 - The lines below are commented because
            #   we still can't rely on either one of the RadialVelocity modules
            #   to calculate the source radial velocities for all kinds of sources.
            #   By uncommenting these one sets dependency of rv measurements on
            #   the stack module. The source RVs are important for the stack module
            #   only if the source's radial velocity is large and variable.
            listofobjectrvs[objectNameWithoutWhites].append(radvel2target)
            listofobjectrvsDependencies[objectNameWithoutWhites].append(radvel2key)
        ###############################################################

        ### Create TARGETS for Radial Velocity #########
        radvelkey = "RADIALVELOCITY1" + basename
        radveltarget = Dirs.PRODUCTDIR + basename + ".rv.gz"
        radvelexttarget = Dirs.PRODUCTDIR + basename + ".rv_ext.gz"
        objproducts[radvelkey] = radveltarget
        objdependencies[radvelkey] = [spckey, rvelkey]
        objcommands[radvelkey] = RadialVelocityCommand(Dirs, radveltarget, radvelexttarget, spctarget, rveltarget, config, Instmode, mjdate ) + verstr
        objcommandsType[radvelkey] = True
        ###############################################################

        ###############################################################
        onedspeckey = "ONEDSPEC" + basename
        onedspectarget = Dirs.PRODUCTDIR + basename + ".1d.spc"
        objproducts[onedspeckey] = onedspectarget
        objdependencies[onedspeckey] = [spckey,rvelkey]
        dummylist = [spctarget]
        objcommands[onedspeckey] = StackObjectSpectra(Dirs, onedspectarget, objectNameWithoutWhites, Instmode, dummylist, "", 0, 2) + verstr
        objcommandsType[onedspeckey] = True
        ###############################################################

        ### Create TARGETS for final FITS spectrum *.fits.gz ##############
        spcfitskey = "OPFITSSPC" + basename
        spcfitstarget = Dirs.PRODUCTDIR + basename + ".m.fits.gz"
        objproducts[spcfitskey] = spcfitstarget
        objdependencies[spcfitskey] = ["WAVELENGTHPRODUCT",spckey,rvelkey,tellkey,snrkey]
        objcommands[spcfitskey] = operaCreateProduct_SPC_Intensity(Dirs, file, spctarget, spcfitstarget, rveltarget, telltarget, snrtarget, products["WAVELENGTHRESOLUTIONPRODUCT"], config.OPERAVERSION, currentTime) + verstr
        objcommandsType[spcfitskey] = True
        ###############################################################

        ### Create TARGETS for flux calibration Standards ######
        standardname = ""

        if config.hasStandard(objectname.replace(" ", "")) :
            standardname = objectname.replace(" ", "")
        elif config.hasStandard(objectname) :
            standardname = objectname
    
        if standardname :
            stdfcalkey = "STDFCAL" + basename
            stdfcaltarget = Dirs.PRODUCTDIR + basename + ".fcal.gz"

            listofstdfcalkey.append(stdfcalkey)
            listofstdfcaltargets.append(stdfcaltarget)

            objproducts[stdfcalkey] = stdfcaltarget
            objdependencies[stdfcalkey] = [extkey,"APERTUREPRODUCT","WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
            objcommands[stdfcalkey] = CreateFcalCommand(Dirs, stdfcaltarget, exttarget, config.getStandardDataFile(standardname), products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode, products["APERTUREPRODUCT"], products["WAVELENGTHPRODUCT"], exptime) + verstr
            objcommandsType[stdfcalkey] = True
            
            # -- The target below produces a LE-style flat response for each observation of
            #    standard stars. However it could be limited to the Moon spectrum only.
            #    E. Martioli -- Feb 22 2015
            flatrespkey = "FLATRESPONSE" + basename
            flatresptarget = Dirs.PRODUCTDIR + basename + "_flat_resp.fits.gz"
            objproducts[flatrespkey] = flatresptarget
            objdependencies[flatrespkey] = [extkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
            objcommands[flatrespkey] = CreateFlatResponseCommand(Dirs, flatresptarget, exttarget, file, config.getStandardDataFile(standardname), products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode, products["APERTUREPRODUCT"], products["WAVELENGTHPRODUCT"]) + verstr
            objcommandsType[flatrespkey] = True
        ###############################################################
        
        ### Create TARGETS for LE formats #############################
        lespcnwkey = "LESPCNW" + basename
        lespcnwtarget = Dirs.PRODUCTDIR + basename + ".inw.s.gz"
        objproducts[lespcnwkey] = lespcnwtarget
        objdependencies[lespcnwkey] = [spckey]
        objcommands[lespcnwkey] = GenLEFormatsCommand(Dirs, config, lespcnwtarget, spctarget, Instmode, objectname, 2, 3) + verstr
        objcommandsType[lespcnwkey] = True
        #
        lespcukey = "LESPCU" + basename
        lespcutarget = Dirs.PRODUCTDIR + basename + ".iu.s.gz"
        objproducts[lespcukey] = lespcutarget
        objdependencies[lespcukey] = [spckey]
        objcommands[lespcukey] = GenLEFormatsCommand(Dirs, config, lespcutarget, spctarget, Instmode, objectname, 3, 4) + verstr
        objcommandsType[lespcukey] = True
        #
        lespcuwkey = "LESPCUW" + basename
        lespcuwtarget = Dirs.PRODUCTDIR + basename + ".iuw.s.gz"
        objproducts[lespcuwkey] = lespcuwtarget
        objdependencies[lespcuwkey] = [spckey]
        objcommands[lespcuwkey] = GenLEFormatsCommand(Dirs, config, lespcuwtarget, spctarget, Instmode, objectname, 3, 3) + verstr
        objcommandsType[lespcuwkey] = True
        #
        lespcnkey = "LESPCN" + basename
        lespcntarget = Dirs.PRODUCTDIR + basename + ".in.s.gz"
        objproducts[lespcnkey] = lespcntarget
        objdependencies[lespcnkey] = [spckey]
        objcommands[lespcnkey] = GenLEFormatsCommand(Dirs, config, lespcntarget, spctarget, Instmode, objectname, 2, 4) + verstr
        objcommandsType[lespcnkey] = True
        #
        fitsLEproductkey = "LEFITSPRODUCT" + basename
        fitsLEproducttarget = Dirs.PRODUCTDIR + basename + ".i.fits.gz"
        objproducts[fitsLEproductkey] = fitsLEproducttarget
        objdependencies[fitsLEproductkey] = ["WAVELENGTHPRODUCT",lespcnwkey,lespcukey,lespcuwkey,lespcnkey,rvelkey,tellkey,snrkey]
        objcommands[fitsLEproductkey] = operaCreateProduct_LE_Intensity(Dirs, file, fitsLEproducttarget, lespcutarget, lespcntarget, lespcuwtarget, lespcnwtarget, rveltarget, telltarget, snrtarget, products["WAVELENGTHRESOLUTIONPRODUCT"], config.OPERAVERSION, currentTime) + verstr
        objcommandsType[fitsLEproductkey] = True
        ###############################################################

        if (Instmode.mode == 3) :
            command = Dirs.EXE +"operagetheader " + file + " --keyword=CMMTSEQ | awk '{print $3}' | sed -e s/,//"

            POLAREXPOSURE = subprocess.check_output(command,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

            if (int(POLAREXPOSURE) == 1 and polcount == 0) :
                ### Create TARGETS for raw polarimetry spectrum *.p.gz ########
                pkey = "POLAR" + basename
                ptarget = Dirs.PRODUCTDIR + basename + ".p.gz"
                objdependencies[pkey] = []
                objcommands[pkey] = ""
                objproducts[pkey] = ptarget
                objcommandsType[pkey] = True
                ###############################################################

                ### Create TARGETS for final polarimetry spectrum *.pol.gz ####
                polkey = "OPPOL" + basename
                poltarget = Dirs.PRODUCTDIR + basename + ".pol.gz"
                objproducts[polkey] = poltarget
                objdependencies[polkey] = [pkey,tellkey,rvelkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
                #objcommands[polkey] = CalibratedPolarCommand(Dirs, poltarget, Instmode, config, ptarget, products["FLATFLUXCALIBRATIONSPECTRUM"], objproducts["MASTERFLUXCALIBRATION"], rveltarget, telltarget, products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
                objcommands[polkey] = CalibratedPolarCommand(Dirs, poltarget, Instmode, config, ptarget, products["FLATFLUXCALIBRATIONSPECTRUM"], emptystring, rveltarget, telltarget, products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
                objcommandsType[polkey] = True
                ###############################################################

                ### Create TARGETS for final FITS polarimetry *pol.fits.gz ##############
                polfitskey = "OPFITSPOL" + basename
                polfitstarget = Dirs.PRODUCTDIR + basename + ".pol.fits.gz"
                objproducts[polfitskey] = polfitstarget
                objdependencies[polfitskey] = ["WAVELENGTHPRODUCT",polkey,rvelkey,tellkey,snrkey]
                objcommands[polfitskey] = operaCreateProduct_SPC_Polar(Dirs, file, poltarget, polfitstarget, rveltarget, telltarget, snrtarget, products["WAVELENGTHRESOLUTIONPRODUCT"], config.OPERAVERSION, currentTime) + verstr
                objcommandsType[polfitskey] = True
                ###############################################################

                ### Create TARGETS for Polarimetric LE formats #############################
                lepolnwkey = "LEPOLNW" + basename
                lepolnwtarget = Dirs.PRODUCTDIR + basename + ".pnw.s.gz"
                objproducts[lepolnwkey] = lepolnwtarget
                objdependencies[lepolnwkey] = [polkey]
                objcommands[lepolnwkey] = GenLEFormatsCommand(Dirs, config, lepolnwtarget, poltarget, Instmode, objectname, 2, 3) + verstr
                objcommandsType[lepolnwkey] = True
                #
                lepolukey = "LEPOLU" + basename
                lepolutarget = Dirs.PRODUCTDIR + basename + ".pu.s.gz"
                objproducts[lepolukey] = lepolutarget
                objdependencies[lepolukey] = [polkey]
                objcommands[lepolukey] = GenLEFormatsCommand(Dirs, config, lepolutarget, poltarget, Instmode, objectname, 3, 4) + verstr
                objcommandsType[lepolukey] = True
                #
                lepoluwkey = "LEPOLUW" + basename
                lepoluwtarget = Dirs.PRODUCTDIR + basename + ".puw.s.gz"
                objproducts[lepoluwkey] = lespcuwtarget
                objdependencies[lepoluwkey] = [polkey]
                objcommands[lepoluwkey] = GenLEFormatsCommand(Dirs, config, lepoluwtarget, poltarget, Instmode, objectname, 3, 3) + verstr
                objcommandsType[lepoluwkey] = True
                #
                lepolnkey = "LEPOLN" + basename
                lepolntarget = Dirs.PRODUCTDIR + basename + ".pn.s.gz"
                objproducts[lepolnkey] = lepolntarget
                objdependencies[lepolnkey] = [polkey]
                objcommands[lepolnkey] = GenLEFormatsCommand(Dirs, config, lepolntarget, poltarget, Instmode, objectname, 2, 4) + verstr
                objcommandsType[lepolnkey] = True
                ###############################################################
                
                ### Create TARGETS for final LE FITS polarimetry *p.fits.gz ##############
                lepolfitskey = "LEFITSPOL" + basename
                lepolfitstarget = Dirs.PRODUCTDIR + basename + ".p.fits.gz"
                objproducts[lepolfitskey] = polfitstarget
                objdependencies[lepolfitskey] = ["WAVELENGTHPRODUCT",lepolnwkey,lepolukey,lepoluwkey,lepolnkey,polkey,rvelkey,tellkey,snrkey]
                objcommands[lepolfitskey] = operaCreateProduct_LE_Polar(Dirs, file, lepolfitstarget, lepolutarget, lepolntarget, lepoluwtarget, lepolnwtarget, rveltarget, telltarget, snrtarget, products["WAVELENGTHRESOLUTIONPRODUCT"], config.OPERAVERSION, currentTime) + verstr
                objcommandsType[lepolfitskey] = True
                ###############################################################

                polquad.append(exttarget)
                polextkey.append(extkey)
                polcount += 1
            
            elif( (int(POLAREXPOSURE) == 2 or int(POLAREXPOSURE) == 3) and polcount > 0 and polcount < 3) :
                polquad.append(exttarget)
                polextkey.append(extkey)                
                polcount +=1
            elif( int(POLAREXPOSURE) == 4 and polcount == 3) :
                polquad.append(exttarget)
                polextkey.append(extkey)
                
                stokescommand = Dirs.EXE +"operagetstokes --numeric " + file
                stokes = subprocess.check_output(stokescommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
                   
                objdependencies[pkey] = [polextkey[0], polextkey[1], polextkey[2], polextkey[3],"WAVELENGTHPRODUCT"]
                objcommands[pkey] = PolarCommand(Dirs, Instmode, ptarget, polquad[0], polquad[1], polquad[2], polquad[3], products["WAVELENGTHPRODUCT"], stokes) + verstr
                objcommandsType[pkey] = True

                pkey = ""
                ptarget = ""
                polquad = []
                polextkey = []
                polcount = 0

    for object in listofobjects :
        stacknormobjectkey = "STACKNORM" + object
        stacknormobjecttarget = Dirs.PRODUCTDIR + object + "_" + night + "_master.norm.spc"
        objproducts[stacknormobjectkey] = stacknormobjecttarget
        objdependencies[stacknormobjectkey] = listofobjectspcsDependencies[object] + listofobjectrvsDependencies[object]
        objcommands[stacknormobjectkey] = StackObjectSpectra(Dirs, stacknormobjecttarget, object, Instmode, listofobjectspcs[object], listofobjectrvs[object], 0, 1) + verstr
        objcommandsType[stacknormobjectkey] = True
        
        stackrawobjectkey = "STACKRAW" + object
        stackrawobjecttarget = Dirs.PRODUCTDIR + object + "_" + night + "_master.raw.spc"
        objproducts[stackrawobjectkey] = stackrawobjecttarget
        objdependencies[stackrawobjectkey] = listofobjectspcsDependencies[object] + listofobjectrvsDependencies[object]
        objcommands[stackrawobjectkey] = StackObjectSpectra(Dirs, stackrawobjecttarget, object, Instmode, listofobjectspcs[object], listofobjectrvs[object], 0, 0) + verstr
        objcommandsType[stackrawobjectkey] = True
        
        stackfcalobjectkey = "STACKFCAL" + object
        stackfcalobjecttarget = Dirs.PRODUCTDIR + object + "_" + night + "_master.fcal.spc"
        objproducts[stackfcalobjectkey] = stackfcalobjecttarget
        objdependencies[stackfcalobjectkey] = listofobjectspcsDependencies[object] + listofobjectrvsDependencies[object]
        objcommands[stackfcalobjectkey] = StackObjectSpectra(Dirs, stackfcalobjecttarget, object, Instmode, listofobjectspcs[object], listofobjectrvs[object], 0, 2) + verstr
        objcommandsType[stackfcalobjectkey] = True



    # E. Martioli Jun 22 2016 -- The part below is to produce a master flat response
    # out of standard observations and include this into the spc file.
    ### Create DEPENDENCIES and COMMANDLINE for master flat response ######
    hasStandard = False
    if listofstdfcalkey :
        masterfcalkeytarget = Dirs.PRODUCTDIR + night + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_flat_resp.fits.gz"
        objproducts[masterfcalkey] = masterfcalkeytarget
        objdependencies[masterfcalkey] = listofstdfcalkey
        objcommands[masterfcalkey] = MasterFlatResponseCommand(Dirs, "operaMedianCombine", masterfcalkeytarget, listofstdfcaltargets) + verstr
        objcommandsType[masterfcalkey] = True
        hasStandard = True

    for file in objfilelist:
        base=os.path.basename(file)
        if ".fits.gz" in base :
            basename = os.path.splitext(os.path.splitext(base)[0])[0]
        elif ".fits" in base :
            basename = os.path.splitext(base)[0]
        elif ".fit" in base :
            basename = os.path.splitext(base)[0]
        else :
            print "Error: unknown extension, exiting ..."
            exit()
        #### Append flux calib stuff to command lines for calibrated spectrum *.spc ########
        spckey = "OPSPC" + basename
        if not hasStandard :
            objdependencies[spckey].remove(masterfcalkey);

        objcommands[spckey] += ' --flatResponse=' + masterfcalkeytarget  + ' --inputWavelengthMaskForUncalContinuum=' + config.ATYPEWAVELENGTHMASK
        objcommandsType[spckey] = True
    ###############################################################
    
    return objproducts, objdependencies, objcommands,objcommandsType
#-------------------------------------------

######### PRODUCT FILE NAMES ###########
def setCalibrationProducts(Dirs,night,Instmode,Readmode,DefaultCal,Keywords,allowanyreadout,config, plots,verbose) :

    INSTCONFIGPREFIX = Dirs.PRODUCTDIR + night + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME

    products = {}
    dependencies = {}
    commands = {}
    commandsType = {}
    
    verstr = ""
    if (verbose) :
        verstr = " --verbose"
    
    command = ""

    biascheck = testBiasListCommand(Dirs, Instmode, Readmode, Keywords)
    if subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["BIASLIST"] = INSTCONFIGPREFIX + "_bias.list"
        dependencies["BIASLIST"] = []
        commands["BIASLIST"] = BiasListCommand(Dirs, Instmode, Readmode, Keywords, products["BIASLIST"])
        commandsType["BIASLIST"] = True
        
        products["MASTERBIAS"] = INSTCONFIGPREFIX + "_masterbias.fits.gz"
        dependencies["MASTERBIAS"] = ["BIASLIST"]
        commands["MASTERBIAS"] = MasterCalibrationCommand(Dirs, "operaMasterBias", products["MASTERBIAS"], products["BIASLIST"]) + verstr
        commandsType["MASTERBIAS"] = True
    else :
        products["MASTERBIAS"] = DefaultCal.DEFAULTMASTERBIAS
        dependencies["MASTERBIAS"] = []

    flatcheck = testCalibrationListCommand(Dirs, Instmode, Readmode, Keywords.FLATKEYWORD, Keywords, allowanyreadout)
    if subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["FLATLIST"] = INSTCONFIGPREFIX + "_flat.list"
        dependencies["FLATLIST"] = []
        commands["FLATLIST"] = CalibrationListCommand(Dirs, Instmode, Readmode, Keywords.FLATKEYWORD, Keywords, products["FLATLIST"],allowanyreadout)
        commandsType["FLATLIST"] = True
        
        products["MASTERFLAT"] = INSTCONFIGPREFIX + "_masterflat.fits.gz"
        dependencies["MASTERFLAT"] = ["FLATLIST"]
        commands["MASTERFLAT"] = MasterCalibrationCommand(Dirs, "operaMasterFlat", products["MASTERFLAT"], products["FLATLIST"]) + verstr
        commandsType["MASTERFLAT"] = True
    else :
        products["MASTERFLAT"] = DefaultCal.DEFAULTMASTERFLAT
        dependencies["MASTERFLAT"] = []

    compcheck = testCalibrationListCommand(Dirs, Instmode, Readmode, Keywords.COMPKEYWORD, Keywords, allowanyreadout)
    if subprocess.check_output(compcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["COMPLIST"] = INSTCONFIGPREFIX + "_comp.list"
        dependencies["COMPLIST"] = []
        commands["COMPLIST"] = CalibrationListCommand(Dirs, Instmode, Readmode, Keywords.COMPKEYWORD, Keywords, products["COMPLIST"],allowanyreadout)
        commandsType["COMPLIST"] = True
        
        products["MASTERCOMP"] = INSTCONFIGPREFIX + "_mastercomp.fits.gz"
        dependencies["MASTERCOMP"] = ["COMPLIST"]
        commands["MASTERCOMP"] = MasterComparisonCommand(Dirs, products["MASTERCOMP"], products["COMPLIST"], config.BADPIXELMASK, products["MASTERBIAS"], Instmode) + verstr
        commandsType["MASTERCOMP"] = True
    else :
        products["MASTERCOMP"] = DefaultCal.DEFAULTMASTERCOMP
        dependencies["MASTERCOMP"] = []

    if (Instmode.mode == 1 or Instmode.mode == 3) :
        aligncheck = testCalibrationListCommand(Dirs, Instmode, Readmode, Keywords.ALIGNKEYWORD, Keywords, allowanyreadout)
        if subprocess.check_output(aligncheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
            products["ALIGNLIST"] = INSTCONFIGPREFIX + "_align.list"
            dependencies["ALIGNLIST"] =[]
            commands["ALIGNLIST"] = CalibrationListCommand(Dirs, Instmode, Readmode, Keywords.ALIGNKEYWORD, Keywords, products["ALIGNLIST"],allowanyreadout)
            commandsType["ALIGNLIST"] = True
            
            products["MASTERALIGN"] = INSTCONFIGPREFIX + "_masteralign.fits.gz"
            dependencies["MASTERALIGN"] = ["ALIGNLIST"]
            commands["MASTERALIGN"] = MasterCalibrationCommand(Dirs, "operaMasterFlat", products["MASTERALIGN"], products["ALIGNLIST"]) + verstr
            commandsType["MASTERALIGN"] = True
        else :
            if "COMPLIST" in products :
                products["MASTERALIGN"] = products["MASTERCOMP"]
            else :
                products["MASTERALIGN"] = DefaultCal.DEFAULTMASTERALIGN
                dependencies["MASTERALIGN"] = []

    if subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') and \
        subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["GAINPRODUCT"] = INSTCONFIGPREFIX + ".gain.gz"
        dependencies["GAINPRODUCT"] = ["BIASLIST","FLATLIST"]
        commands["GAINPRODUCT"] = GainCommand(Dirs,products["GAINPRODUCT"], products["BIASLIST"], products["FLATLIST"], config.BADPIXELMASK, Readmode, Instmode) + verstr
        commandsType["GAINPRODUCT"] = True
    else :
        products["GAINPRODUCT"] = DefaultCal.DEFAULTGAIN
        dependencies["GAINPRODUCT"] = []

    products["ORDERSPACINGPRODUCT"] = INSTCONFIGPREFIX + ".ordp.gz"
    dependencies["ORDERSPACINGPRODUCT"] = ["GAINPRODUCT","MASTERBIAS","MASTERFLAT"]
    commands["ORDERSPACINGPRODUCT"] = OrderSpacingCommand(Dirs,products["ORDERSPACINGPRODUCT"],products["GAINPRODUCT"],products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,Instmode,plots["ORDERSPACINGPRODUCT"]) + verstr
    commandsType["ORDERSPACINGPRODUCT"] = True
    
    products["GEOMETRYPRODUCT"] = INSTCONFIGPREFIX + ".geom.gz"
    dependencies["GEOMETRYPRODUCT"] = ["GAINPRODUCT","MASTERBIAS","MASTERFLAT","ORDERSPACINGPRODUCT"]
    commands["GEOMETRYPRODUCT"] = GeometryCommand(Dirs,products["GEOMETRYPRODUCT"],products["GAINPRODUCT"],products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["ORDERSPACINGPRODUCT"],Instmode,plots["GEOMETRYPRODUCT"]) + verstr
    commandsType["GEOMETRYPRODUCT"] = True
    
    products["INSTRUMENTPROFILEPRODUCT"] = INSTCONFIGPREFIX + ".prof.gz"
    commandsType["INSTRUMENTPROFILEPRODUCT"] = True
    if (Instmode.mode == 1 or Instmode.mode == 3) :
        dependencies["INSTRUMENTPROFILEPRODUCT"] = ["GEOMETRYPRODUCT","GAINPRODUCT","MASTERBIAS","MASTERFLAT","MASTERALIGN"]
        commands["INSTRUMENTPROFILEPRODUCT"] = InstrumentProfileCommand(Dirs,products["INSTRUMENTPROFILEPRODUCT"], products["GEOMETRYPRODUCT"],products["GAINPRODUCT"],products["MASTERBIAS"],products["MASTERFLAT"],products["MASTERALIGN"],config.BADPIXELMASK,Instmode,plots["INSTRUMENTPROFILEPRODUCT"]) + verstr
    elif (Instmode.mode == 2) :
        dependencies["INSTRUMENTPROFILEPRODUCT"] = ["GEOMETRYPRODUCT","GAINPRODUCT","MASTERBIAS","MASTERFLAT","MASTERCOMP"]
        commands["INSTRUMENTPROFILEPRODUCT"] = InstrumentProfileCommand(Dirs,products["INSTRUMENTPROFILEPRODUCT"], products["GEOMETRYPRODUCT"],products["GAINPRODUCT"],products["MASTERBIAS"],products["MASTERFLAT"],products["MASTERCOMP"],config.BADPIXELMASK,Instmode,plots["INSTRUMENTPROFILEPRODUCT"]) + verstr

    products["APERTUREPRODUCT"] = INSTCONFIGPREFIX + ".aper.gz"
    dependencies["APERTUREPRODUCT"] = ["GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","ORDERSPACINGPRODUCT"]
    commands["APERTUREPRODUCT"] = ApertureCommand(Dirs,products["APERTUREPRODUCT"],products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["ORDERSPACINGPRODUCT"],Instmode,plots["APERTUREPRODUCT"]) + verstr
    commandsType["APERTUREPRODUCT"] = True
    
    products["COMPEXTRACTEDSPECTRUM"] = INSTCONFIGPREFIX + "_comp.e.gz"
    dependencies["COMPEXTRACTEDSPECTRUM"] = ["MASTERCOMP","MASTERBIAS","MASTERFLAT","GAINPRODUCT","GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","APERTUREPRODUCT"]
    commands["COMPEXTRACTEDSPECTRUM"] = compRawExtractionCommand(Dirs,products["COMPEXTRACTEDSPECTRUM"],products["MASTERCOMP"],products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["GAINPRODUCT"], products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["APERTUREPRODUCT"],Instmode) + verstr
    commandsType["COMPEXTRACTEDSPECTRUM"] = True
    
    products["FLATEXTRACTEDSPECTRUM"] = INSTCONFIGPREFIX + "_flat.e.gz"
    dependencies["FLATEXTRACTEDSPECTRUM"] = ["MASTERFLAT","MASTERBIAS","MASTERFLAT","GAINPRODUCT","GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","APERTUREPRODUCT"]
    commands["FLATEXTRACTEDSPECTRUM"] = calibrationExtractionCommand(Dirs,products["FLATEXTRACTEDSPECTRUM"],products["MASTERFLAT"],products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["GAINPRODUCT"], products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["APERTUREPRODUCT"],Instmode) + verstr
    commandsType["FLATEXTRACTEDSPECTRUM"] = True
    
    products["WAVELENGTHPRODUCT"] = INSTCONFIGPREFIX + ".wcal.gz"
    dependencies["WAVELENGTHPRODUCT"] = ["GEOMETRYPRODUCT","COMPEXTRACTEDSPECTRUM"]
    products["WAVELENGTHRESOLUTIONPRODUCT"] = INSTCONFIGPREFIX + ".sres.gz"
    dependencies["WAVELENGTHRESOLUTIONPRODUCT"] = ["WAVELENGTHPRODUCT"]
    inputwaveguess = ' --inputWaveFile=' + config.WAVEFIRSTGUESS
    commands["WAVELENGTHPRODUCT"] = WavelengthCommand(Dirs, products["WAVELENGTHPRODUCT"], products["WAVELENGTHRESOLUTIONPRODUCT"], products["GEOMETRYPRODUCT"], products["COMPEXTRACTEDSPECTRUM"], Instmode, config ,inputwaveguess, plots["WAVELENGTHPRODUCT"]) + verstr
    commands["WAVELENGTHRESOLUTIONPRODUCT"] = ""
    commandsType["WAVELENGTHRESOLUTIONPRODUCT"] = True
    commandsType["WAVELENGTHPRODUCT"] = True
    
    products["FLATFLUXCALIBRATIONSPECTRUM"] = INSTCONFIGPREFIX + "_flat.fcal.gz"
    dependencies["FLATFLUXCALIBRATIONSPECTRUM"] = ["FLATEXTRACTEDSPECTRUM","WAVELENGTHPRODUCT"]
    commands["FLATFLUXCALIBRATIONSPECTRUM"] = FlatFluxCalibrationCommand(Dirs, products["FLATFLUXCALIBRATIONSPECTRUM"],products["FLATEXTRACTEDSPECTRUM"], Instmode, products["WAVELENGTHPRODUCT"]) + verstr
    commandsType["FLATFLUXCALIBRATIONSPECTRUM"] = True
    
    return products,dependencies,commands,commandsType
##################################

######### PRODUCT FILE NAMES ###########
def setPlotFilenames(Dirs,night,Instmode,Readmode, plotbool) :
    
    INSTCONFIGPREFIX = Dirs.PRODUCTDIR + night + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME
    
    plots = {}

    if plotbool : 
        plots["ORDERSPACINGPRODUCT"] = {'ORDSPCPLOTFILE': INSTCONFIGPREFIX + "_spcplot.eps", \
                                        'ORDSPCDATAFILE': INSTCONFIGPREFIX + "_spcplot.dat", \
                                        'ORDSPCSCRIPTFILE': INSTCONFIGPREFIX + "_spcplot.gnu"}
    
        plots["GEOMETRYPRODUCT"] =  {'GEOMPLOTFILE': INSTCONFIGPREFIX + "_geomplot.eps", \
                                     'GEOMDATAFILE': INSTCONFIGPREFIX + "_geomplot.dat", \
                                     'GEOMSCRIPTFILE': INSTCONFIGPREFIX + "_geomplot.gnu"}
    
        plots["INSTRUMENTPROFILEPRODUCT"] =  {'PROFPLOTFILE': INSTCONFIGPREFIX + "_profplot.eps", \
                                              'PROFDATAFILE': INSTCONFIGPREFIX + "_profplot.dat", \
                                              'PROFSCRIPTFILE': INSTCONFIGPREFIX + "_profplot.gnu"}
                                              
        plots["APERTUREPRODUCT"] =  {'APERPLOTFILE': INSTCONFIGPREFIX + "_aperplot.eps", \
                                     'APERDATAFILE': INSTCONFIGPREFIX + "_aperplot.dat", \
                                     'APERSCRIPTFILE': INSTCONFIGPREFIX + "_aperplot.gnu", \
                                     'APERTILTPLOTFILE': INSTCONFIGPREFIX + "_tiltplot.eps", \
                                     'APERTILTDATA1FILE': INSTCONFIGPREFIX + "_tiltplot1.dat", \
                                     'APERTILTDATA2FILE': INSTCONFIGPREFIX + "_tiltplot2.dat", \
                                     'APERTILTSCRIPTFILE': INSTCONFIGPREFIX + "_tiltplot.gnu"}

        plots["WAVELENGTHPRODUCT"] =  {'WAVEORDSPLOTFILE': INSTCONFIGPREFIX + "_waveordsplot.eps", \
                                            'WAVESPECPLOTFILE': INSTCONFIGPREFIX + "_wavespecplot.eps", \
                                            'WAVESPECSCRIPTFILE': INSTCONFIGPREFIX + "_wavespecplot.gnu", \
                                            'WAVEORDSCRIPTFILE': INSTCONFIGPREFIX + "_waveordplot.gnu", \
                                            'WAVEORDSDATAFILE': INSTCONFIGPREFIX + "_waveordsplot.dat", \
                                            'WAVEATLASDATAFILE': INSTCONFIGPREFIX + "_waveatlasplot.dat", \
                                            'WAVECOMPDATAFILE': INSTCONFIGPREFIX + "_wavecompplot.dat", \
                                            'WAVELINESDATAFILE': INSTCONFIGPREFIX + "_wavelinesplot.dat"}

# Uncomment the lines below and remove the one above to go back using StitchOrders
        """        plots["FIRSTWAVELENGTHPRODUCT"] =  {'WAVEORDSPLOTFILE': INSTCONFIGPREFIX + "_waveordsplot.eps", \
                                            'WAVESPECPLOTFILE': INSTCONFIGPREFIX + "_wavespecplot.eps", \
                                            'WAVESPECSCRIPTFILE': INSTCONFIGPREFIX + "_wavespecplot.gnu", \
                                            'WAVEORDSCRIPTFILE': INSTCONFIGPREFIX + "_waveordplot.gnu", \
                                            'WAVEORDSDATAFILE': INSTCONFIGPREFIX + "_waveordsplot.dat", \
                                            'WAVEATLASDATAFILE': INSTCONFIGPREFIX + "_waveatlasplot.dat", \
                                            'WAVECOMPDATAFILE': INSTCONFIGPREFIX + "_wavecompplot.dat", \
                                            'WAVELINESDATAFILE': INSTCONFIGPREFIX + "_wavelinesplot.dat"}

        plots["WAVELENGTHPRODUCT"] = {'STITCHPLOTFILE': INSTCONFIGPREFIX + "_stitchplot.eps", \
                                      'STITCHDATAFILE': INSTCONFIGPREFIX + "_stitchplot.dat", \
                                      'STITCHSCRIPTFILE': INSTCONFIGPREFIX + "_stitchplot.gnu"}

        """

    else :
        plots["ORDERSPACINGPRODUCT"] = {'ORDSPCPLOTFILE': "",'ORDSPCDATAFILE': "",'ORDSPCSCRIPTFILE': ""}
        plots["GEOMETRYPRODUCT"] = {'GEOMPLOTFILE': "",'GEOMDATAFILE': "",'GEOMSCRIPTFILE': ""}
        plots["INSTRUMENTPROFILEPRODUCT"] = {'PROFPLOTFILE': "",'PROFDATAFILE': "",'PROFSCRIPTFILE': ""}
        plots["APERTUREPRODUCT"] = {'APERPLOTFILE': "",'APERDATAFILE': "",'APERSCRIPTFILE': "",\
                                    'APERTILTPLOTFILE': "",'APERTILTDATA1FILE': "",'APERTILTDATA2FILE': "",'APERTILTSCRIPTFILE': ""}
        plots["WAVELENGTHPRODUCT"] = {'WAVEORDSPLOTFILE': "",'WAVESPECPLOTFILE': "",'WAVESPECSCRIPTFILE': "", \
                                           'WAVEORDSCRIPTFILE': "",'WAVEORDSDATAFILE': "",'WAVEATLASDATAFILE': "", \
                                           'WAVECOMPDATAFILE': "", 'WAVELINESDATAFILE': ""}
        # Uncomment the lines below and remove the one above to go back using StitchOrders
        """        plots["FIRSTWAVELENGTHPRODUCT"] = {'WAVEORDSPLOTFILE': "",'WAVESPECPLOTFILE': "",'WAVESPECSCRIPTFILE': "", \
                                           'WAVEORDSCRIPTFILE': "",'WAVEORDSDATAFILE': "",'WAVEATLASDATAFILE': "", \
                                           'WAVECOMPDATAFILE': "", 'WAVELINESDATAFILE': ""}
        """
    return plots
##################################


############# Class to encapsulate modes for reduction ####################
# This class contains the information on all available 
class ReductionModes :
    'Common base class for reduction modes'

    # modes = {[modekey]: [modeName, InstMode, ReadMode, Nobjects, Nbiases, Nflats, Ncomps, Naligns] }
    modes = {}
    files = {}
    
    def __init__(self, dirs, keywords, allowanyreadout, forcecalibration):
        if os.path.exists(dirs.DATADIR) :
            self.displayStats = "---\n"
            self.displayStats += "STATISTICS for DATA in DIR: " + dirs.DATADIR + "\n"
            self.displayStats += "Allow any readout? " + str(bool(allowanyreadout)) + "\n"
            self.displayStats += "Include mode for calibration alone? " + str(forcecalibration) + "\n"
            self.displayStats += "---\n"
            self.displayStats += "InstMode\tReadout\tobject\tbias\tflat\tcomp\talign\tSELECTED?\n"
            self.displayStats += "-------------------------------------------------------------------\n"

            for mode in range(1,4) :
                for read in range(1,7) :

                    Instmode = InstMode(mode)
                    Readmode = ReadoutMode(read)
                    
                    modeKey = Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME
                    
                    nobjects = 0
                    nbiases = 0
                    nflats = 0
                    ncomps = 0
                    naligns = 0
                    
                    objectcheck = ObjectListCommand(dirs, Instmode, Readmode, keywords.OBJECTKEYWORD,keywords)
                    objects = (subprocess.check_output(objectcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(objects) :
                        nobjects = len(objects)
                    
                    biascheck = testBiasListCommand(dirs, Instmode, Readmode, keywords)
                    biases = (subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(biases) :
                        nbiases = len(biases)
                
                    flatcheck = testCalibrationListCommand(dirs, Instmode, Readmode, keywords.FLATKEYWORD, keywords, allowanyreadout)
                    flats = (subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(flats) :
                        nflats = len(flats)
            
                    compcheck = testCalibrationListCommand(dirs, Instmode, Readmode, keywords.COMPKEYWORD, keywords, allowanyreadout)
                    comps = (subprocess.check_output(compcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(comps) :
                        ncomps = len(comps)
            
                    aligncheck = testCalibrationListCommand(dirs, Instmode, Readmode, keywords.ALIGNKEYWORD, keywords, allowanyreadout)
                    aligns = (subprocess.check_output(aligncheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(aligns) :
                        naligns = len(aligns)
                
                    selected = "NO"
                
                    if forcecalibration :
                        if nbiases and nflats and ncomps :
                            self.modes[modeKey] = [modeKey, mode, read, nobjects, nbiases, nflats, ncomps, naligns]
                            self.files[modeKey] = [objects, biases, flats, comps, aligns]
                            selected = "YES"
                    else :
                        if nobjects and nbiases and nflats and ncomps :
                            self.modes[modeKey] = [modeKey, mode, read, nobjects, nbiases, nflats, ncomps, naligns]
                            self.files[modeKey] = [objects, biases, flats, comps, aligns]
                            selected = "YES"
    
                    xtab = "\t"
                    if mode == 3 :
                        xtab = "\t\t"
                    self.displayStats += Instmode.INSTRUMENTMODESHORTNAME + xtab + Readmode.READOUTSPEEDSHORTNAME + "\t" + str(nobjects) + "\t" + str(nbiases) + "\t" + str(nflats) + "\t" + str(ncomps) + "\t" + str(naligns) + "\t" + selected + "\n"
    
        else :
            print "Error: can\'t find data directory: ", dirs.DATADIR

        self.displayStats += "-------------------------------------------------------------------\n"

    def getNmodes(self) :
        return len(self.modes)
                
    def displayAllModeStats(self) :
        modeItems = self.modes.items()
        for modeItem in modeItems:
            print "MODE SELECTED: " + modeItem[1][0] + "\n" +\
                str(modeItem[1][3]) + " objects  " +  \
                str(modeItem[1][4]) + " biases  " +  \
                str(modeItem[1][5]) + " flats  " + \
                str(modeItem[1][6]) + " comps  " + \
                str(modeItem[1][7]) + " aligns"

    def displayModeStats(self,intrumentmode,readoutspeed) :
        modeItems = self.modes.items()
        for modeItem in modeItems:
            if(modeItem[1][1] == intrumentmode and modeItem[1][2] == readoutspeed) :
                print "MODE SELECTED: " + modeItem[1][0] + "\n" +\
                        str(modeItem[1][3]) + " objects  " +  \
                        str(modeItem[1][4]) + " biases  " +  \
                        str(modeItem[1][5]) + " flats  " + \
                        str(modeItem[1][6]) + " comps  " + \
                        str(modeItem[1][7]) + " aligns"

    def displayOverallStats(self) :
        print self.displayStats

    def getInstReadModes(self) :
        listmodes = []
        modeItems = self.modes.items()
        for modeItem in modeItems:
            listmodes.append([modeItem[1][1],modeItem[1][2]])
        return listmodes

    def displayObjectData(self) :
        filesItems = self.files.items()
        for filesItem in filesItems:
            print "-- "
            print "OBJECTS for MODE: " + filesItem[0]
            for file in filesItem[1][0] :
                print file

    def displayCalibrationData(self) :
        filesItems = self.files.items()
        for filesItem in filesItems:
            print "-- "            
            print "CALIBRATION DATA for MODE: " + filesItem[0]
            print "-- "            
            print "BIAS: "
            for file in filesItem[1][1] :
                print file
            print "-- "            
            print "FLAT: "
            for file in filesItem[1][2] :
                print file
            print "-- "        
            print "COMPARISON: "
            for file in filesItem[1][3] :
                print file
            print "-- "                    
            print "ALIGN: "
            for file in filesItem[1][4] :
                print file

    def cleanModes(self) :
        for key, item in self.modes.items():
            del self.modes[key]
            del self.files[key]


###############################################################
