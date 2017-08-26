# -*- coding: utf-8 -*-
"""
Created on Oct 21 2014
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on Mar 08, 2017
"""
import sys, os
import subprocess
from time import strftime

pylibdir = os.path.split(os.path.dirname(__file__))[0] + '/pyLib'
sys.path.insert(0,pylibdir)

from operacommands import *

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

########## DIRECTORIES ############
class Directories :
     
    'Common base class for directories'
    def __init__(self, pipelinehomedir,datarootdir,productrootdir,night):
        self.OP_HOME=pipelinehomedir
        self.DATAROOTDIR=datarootdir
        self.PRODUCTROOTDIR=productrootdir
        self.PIPELINEDIR=self.OP_HOME+"/pipeline/"
        #self.DEFAULTCALIBRATIONDIR=self.OP_HOME+"/DefaultCalibration/OES/"
        self.OESCONFIGDIR=self.PIPELINEDIR+"/pyOES/config/"
        self.OESMASKSDIR=self.OESCONFIGDIR+"/masks/"
        self.CONFIGDIR=self.OP_HOME+"/config/"
        self.STANDARDSDIR=self.OESCONFIGDIR+"/standardstars/"
        self.EXE=self.OP_HOME+"/bin/"
        self.DATADIR=self.DATAROOTDIR+"/"+night+"/"
        self.PRODUCTDIR=self.PRODUCTROOTDIR+"/"+night+"/"
        self.PYTHONMODULES=self.OP_HOME+"/pipeline/pyOES/"

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
        
        self.BADPIXELMASK = Dirs.OESCONFIGDIR + "OES_badpixelmask.fits.gz"
        self.WAVEFIRSTGUESS = Dirs.OESCONFIGDIR + "OES_wcal_ref.dat.gz"

        self.THARATLASLINES = Dirs.CONFIGDIR + "thar_MM201006.dat.gz"
        self.THARATLASSPECTRUM = Dirs.CONFIGDIR + "LovisPepe_ThArAtlas.dat.gz"
        self.SOLARTYPEWAVELENGTHMASK = Dirs.OESCONFIGDIR + "wavelengthMaskForUncalContinuumDetection_SolarTypeStars.txt"
        self.ATYPEWAVELENGTHMASK = Dirs.OESCONFIGDIR + "wavelengthMaskForUncalContinuumDetection_SolarTypeStars.txt"
        self.ATYPEWAVELENGTHMASK = Dirs.OESCONFIGDIR + "wavelengthMaskForUncalContinuumDetection.txt"
        self.ATLASWAVELENGTHMASK = Dirs.OESCONFIGDIR + "wavelengthMaskForRefContinuumDetection.txt"
        self.TELLURICWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForTelluricAbsorption.txt"
        self.TELLURICLINES = Dirs.CONFIGDIR + "opera_HITRAN08-extracted.par.gz"
        self.TELLURICSPECTRUM = Dirs.CONFIGDIR + "KPNO_atmtrans.dat.gz"

        self.STANDARDLISTFILE = Dirs.STANDARDSDIR + "operaStandardStarsForOES.dat"
        self.readStandards(Dirs,self.STANDARDLISTFILE)
    
        self.SYNTHETICSPECTRUM = Dirs.CONFIGDIR + "/syntethicSpectra/Teff5500.spec"
        self.RVXCORRWAVELENGTHMASK = Dirs.CONFIGDIR + "wavelengthMaskForRVxcorr.txt"
        
        self.OESFLATRESPONSE = ""
        self.LEAPSECONDSFILE = Dirs.CONFIGDIR + "/leapseconds.dat"
    
        self.SOURCELINES = Dirs.CONFIGDIR + "/spectralLinesLibraries/sourcelinesForRVMeasurements.txt"
        
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
        self.BIASKEYWORD="zero"
        self.FLATKEYWORD="flat"
        self.COMPKEYWORD="comp"
        self.OBJECTKEYWORD="object"

        self.INSTRUMEKEY = 'INSTRUME'
        self.OBSTYPEKEY = 'IMAGETYP'
        self.READMODEKEY = 'READSPD'
        self.INSTMODEKEY = 'SLITTYPE'
##############################

######### DEFAULT CALIBRATION ###########
class DefaultCalibration :    
    'Common base class for default calibration'
    def __init__(self, Dirs, Instmode, Readmode):
        self.DEFAULTGAINCAL = Dirs.OESCONFIGDIR + Instmode.INSTRUME + "_" + Readmode.READOUTSPEEDSHORTNAME + ".gain.gz"

        self.DEFAULTMASTERBIAS = Dirs.OESCONFIGDIR + Instmode.INSTRUME + "_" + Readmode.READOUTSPEEDSHORTNAME + "_masterbias.fits.gz"
        self.DEFAULTMASTERCOMP = Dirs.OESCONFIGDIR + Instmode.INSTRUME + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_mastercomp.fits.gz"
        self.DEFAULTMASTERFLAT = Dirs.OESCONFIGDIR + Instmode.INSTRUME + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_masterflat.fits.gz"
        self.DEFAULTMASTERFLUXCAL = Dirs.OESCONFIGDIR + Instmode.INSTRUME + "_" + Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME + ".fcal.gz"
##################################


######### INSTRUMENT MODE ###########
class InstMode :
    'Common base class for all instrument modes'
    mode = 0
    def __init__(self, intrumentmode):
        self.mode = int(intrumentmode)

        self.INSTRUME="OES"
        
        self.observatory_coords = '49:54:55 14:46:52'
        self.observatory_elevation = '528'
            
        self.ROTMIRRCROP_CROPSUBWINDOW='420 1450 0 2048'
        self.ROTMIRRCROP_ROTATE='1'
        self.ROTMIRRCROP_MIRRORCOLS='1'
        self.ROTMIRRCROP_MIRRORROWS='0'
        self.ROTMIRRCROP_SUFFIX='.gz'
        
        self.MINORDERTOEXTRACT = '26'
        self.MAXORDERTOEXTRACT = '81'
        
        self.SPACING_NUMBEROFSAMPLES = '1'
        self.SPACING_SAMPECENTERPOSITION = '1020'
        self.SPACING_DETECTIONMETHOD = '2'
        self.SPACING_APERTURE='12'
        self.SPACING_SUBFORMAT='100 930 1 2048'
        
        self.GEOM_SUBFORMAT='1 1030 1 2048'
        self.GEOM_DETECTIONMETHOD = '2'
        self.GEOM_FFTFILTER = '0'
        self.GEOM_NSAMPLES = '5'
        self.GEOM_ORDEROFTRACINGPOLYNOMIAL = '4'
        self.GEOM_BINSIZE = '15'
        self.GEOM_COLDISPERSION = '1'
        self.GEOM_INVERTORDERS = '1'
        self.GEOM_REFERENCEORDERSAMPLEPOSITION = '1020'
        self.GEOM_GRACES = '1'

        self.PROF_BINSIZE = '100'
        self.PROF_SPECTRALELEMENTHEIGHT = '1.0'
        self.PROF_MAXTHREADS = '4'
        self.PROF_MINIMUMLINES = '5'
        self.PROF_LOCALMAXFILTERWIDTH = '3.0'
        self.PROF_DETECTIONTHRESHOLD = '0.05'
        self.PROF_MINPEAKDEPTH = '0.5'
        self.PROF_IPMETHOD='2'

        self.APER_APERTUREHEIGHT = '0.8'
        self.APER_PICKIMAGEROW = '0'
        self.APER_NROWSAMPLES = '10'
        self.APER_XBIN = '8'
        self.APER_GAP='0'
        self.APER_NUMBEROFBEAMS='1'

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
        self.WAVE_INITIALACCEPTABLEMISMATCH='2.0'
        self.WAVE_NSIGCLIP='2.25'
        self.WAVE_NORMALIZEUNCALIBRATEDSPECTRUM='0'
        self.WAVE_NORMALIZATIONBINSIZE='100'
        self.WAVE_LOCALMAXFILTERWIDTH='5'
        self.WAVE_DETECTIONTHRESHOLD='0.01'
        self.WAVE_MINPEAKDEPTH='0.01'
        self.WAVE_NORDERSTOSEARCHAROUND='3'

        self.FLATFLUXCAL_BINSIZE='500'

        self.EXTRACTION_MAXTHREADS='4'
        self.EXTRACTION_SPECTRUMTYPE='7'
        self.EXTRACTION_SPECTRUMTYPENAME='OptimalBeamSpectrum'
        self.EXTRACTION_BACKGROUNDBINSIZE='300'
        self.EXTRACTION_REMOVEBACKGROUND='0'
        self.EXTRACTION_ITERATIONS='2'
        self.EXTRACTION_ONTARGETPROFILE='1'
        self.EXTRACTION_USEPOLYNOMIALFIT='0'
        self.EXTRACTION_REJECTBADPIXINRAWEXTRACTION='1'
        self.EXTRACTION_MINSIGMACLIP='15'
        self.EXTRACTION_SIGMACLIPRANGE='6'
        self.EXTRACTION_INVERTSKYFIBERFLAG=''

        self.TELL_XCORRELATIONTHRESHOLD='0.1'
        self.TELL_NORMALIZATIONBINSIZE='100'
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

        self.RV_RADIALVELOCITYRANGE='5.0'
        self.RV_RADIALVELOCITYSTEP='0.15'
        self.RV_THRESHOLD='0.05'

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
        self.SPC_NORMALIZATIONBINSIZE='100'
        self.SPC_ABSOLUTECALIBRATION='0'
        self.SPC_MODULENAME='operaStarOnly'
        self.SPC_SKYOVERSTARFIBERAREARATIO=''
        self.SPC_INVERTSKYFIBERFLAG=''
        
        self.STACKSPC_RVBIN = '2.4'    # spectral bin in km/s
        self.STACKSPC_MINWL = '0.0'    # min wavelength in nm (0 for full range)
        self.STACKSPC_MAXWL = '0.0'    # max wavelength in nm (0 for full range)
        self.STACKSPC_SNRCLIP = '2'    # Min SNR to accept a given spectral point
        self.STACKSPC_NCUTATORDERENDS = '50'    # Number of points to exclude at order edges
        self.STACKSPC_APPLYTELL = '1'
        self.STACKSPC_APPLYHELIO = '1'

        if (self.mode == 1) :
            self.INSTRUMENTMODESHORTNAME='OESBLADE'
            self.INSTRUMENTMODEKEY='BLADE'
            self.STARPLUSKYMODEFLAG=0
            self.SPECTRALRESOLUTION='50000'
        
            self.SPACING_REFERENCEORDERNUMBER='63'
            self.SPACING_REFERENCEORDERSEPARATION='18.5'

            self.GEOM_RECENTERIPUSINGSLICESYMMETRY='0'
            self.GEOM_NUMBEROFSLICES='1'
            self.GEOM_APERTURE='10'
            self.GEOM_MAXNORDERS='61'
            self.GEOM_MINORDERTOUSE='23'
                        
            self.PROF_REFERENCELINEWIDTH='1.25'
            self.PROF_IPXSIZE='10'
            self.PROF_IPYSIZE='5'
            self.PROF_IPXSAMPLING='5'
            self.PROF_IPYSAMPLING='5'
            self.PROF_TILTANGLE='5.0'

            self.APER_APERTURE='10'
            self.APER_CONSTANTTILTFLAG='0'
            self.APER_SKYAPERAPERTURE='1'

            self.WAVE_UNCALLINEWIDTH='1.6'
            self.WAVE_THARSELECTEDLINES = 'OES_ThArLines.dat'
            self.WAVE_REFERENCEORDER = '60'

##################################

######### READOUT MODE ###########
class ReadoutMode :
    'Common base class for all readout modes'
    mode = 0
    def __init__(self, readoutmode):
        self.mode = int(readoutmode)
        if (self.mode == 1) :
            self.READOUTSPEEDSHORTNAME="100kHz"
            self.READOUTSPEED="100kHz"
            self.DEFAULTGAIN="2.0"
            self.DEFAULTNOISE="10.0"

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
    objfilelist = objectfiles.split()
    
    listofobjects = []
    listofobjectstackkeys = []
    listofobjecttargets = []
    
    listofobjectspcs = {}
    listofobjectspcsDependencies = {}
    
    listofobjectrvs = {}
    listofobjectrvsDependencies = {}

    listofstdfcalkey = []
    listofstdfcaltargets = []
    masterfcalkeytarget = config.OESFLATRESPONSE
    inputMaskForSPCModule = ""
    masterfcalkey = "MASTERFLATRESP"+Instmode.INSTRUMENTMODESHORTNAME

    for file in objfilelist:
        
        base=os.path.basename(file)
        
        ## Caution: the action below strongly depends on the file name format!!
        # assuming the following possible file names: c201307040004.fits.gz, c201307040004.fits, or c201307040004.fit
        if ".fits.gz" in base :
            basename = os.path.splitext(os.path.splitext(base)[0])[0]
        elif ".fits" in file :
            basename = os.path.splitext(base)[0]
        elif ".fit" in file :
            Instmode.ROTMIRRCROP_SUFFIX='s.gz'
            basename = os.path.splitext(base)[0]
        else :
            print "Error: unknown extension, exiting ..."
            exit()
        
        if verbose :
            print "Object filename: ",file, ", basename = ", basename
        
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
        
        # E. Martioli Mar 24 2015 - the three entries below are not correct, since RA and DEC are
        # in hexa format in OES data and DATE is in UT. These need to be converted to degrees and MJD.
        # E. Martioli Jan 05 2016 - below I'm using SkyCoord and Time from astropy to convert units,
        # but it adds extra dependency, so one needs astropy package installed.

        racommand = Dirs.EXE +"operagetheader --keyword=RA " + file
        absra_center = subprocess.check_output(racommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        #
        deccommand = Dirs.EXE +"operagetheader --keyword=DEC " + file
        absdec_center = subprocess.check_output(deccommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        #
        lstcommand = Dirs.EXE +"operagetheader --keyword=ST " + file
        lst = subprocess.check_output(lstcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')

        coordstr = absra_center + ' ' + absdec_center
        objCoord = SkyCoord(coordstr,unit=(u.hourangle,u.deg),frame='icrs')
        ra_deg = objCoord.ra.deg
        dec_deg = objCoord.dec.deg

        mjcommand = Dirs.EXE +"operagetheader --keyword=DATE " + file
        utFromHeader = subprocess.check_output(mjcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        uttime =  utFromHeader.rstrip(' ')[:10] + ' ' + utFromHeader.rstrip(' ')[11:]
        tt = Time(uttime, format='iso', scale='utc')
        mjdate = tt.mjd

        startHA = '0.0'
        ##################################

        # The entries below are not important
        #airmasscommand = Dirs.EXE +"operagetheader --keyword=AIRMASS " + file
        #airmass = subprocess.check_output(airmasscommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        airmass = "none"
        #windspeedcommand = Dirs.EXE +"operagetheader --keyword=WINDSPED " + file
        #windspeed = subprocess.check_output(windspeedcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        windspeed = "none"

        expnumcommand = Dirs.EXE +"operagetheader --keyword=FILENAME " + file
        expnum = subprocess.check_output(expnumcommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        outsideTemperaturecommand = Dirs.EXE +"operagetheader --keyword=OUTTEMP " + file
        outsideTemperature = subprocess.check_output(outsideTemperaturecommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        relativeHumididycommand = Dirs.EXE +"operagetheader --keyword=AIRHUMEX " + file
        relativeHumididy = subprocess.check_output(relativeHumididycommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        
        atmoPressurecommand = Dirs.EXE +"operagetheader --keyword=AIRPRESS " + file
        atmoPressure = subprocess.check_output(atmoPressurecommand,stderr=subprocess.STDOUT,shell=True).rstrip('\n')
        ########################################

        ### Create TARGETS for Rotate Mirror Crop ###
        newfilekey = "ROTATEMIRRORCROP" + basename
        newfiletarget = Dirs.PRODUCTDIR + basename + ".fits.gz"
        objproducts[newfilekey] = newfiletarget
        objdependencies[newfilekey] = []
        objcommands[newfilekey] = RotateMirrorCropCommand(Dirs, Instmode,file) + verstr
        objcommandsType[newfilekey] = True
        #############################################

        ### Create TARGETS for Extraction #########
        extkey = "EXTRACT" + basename
        exttarget = Dirs.PRODUCTDIR + basename + ".e.gz"
        objproducts[extkey] = exttarget
        objdependencies[extkey] = [newfilekey,"MASTERBIAS","MASTERFLAT","GAINPRODUCT","GEOMETRYPRODUCT","INSTRUMENTPROFILEPRODUCT","APERTUREPRODUCT"]
        objcommands[extkey] = objectExtractionCommand(Dirs, exttarget, newfiletarget, products["MASTERBIAS"],products["MASTERFLAT"],config.BADPIXELMASK,products["GAINPRODUCT"], products["GEOMETRYPRODUCT"],products["INSTRUMENTPROFILEPRODUCT"],products["APERTUREPRODUCT"], Instmode) + verstr
        objcommandsType[extkey] = True
        ########################################

        ### Create TARGETS for telluric wavelength correction #########
        tellkey = "TELLWAVE" + basename
        telltarget = Dirs.PRODUCTDIR + basename + ".tell.gz"
        objproducts[tellkey] = telltarget
        objdependencies[tellkey] = [extkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
        objcommands[tellkey] = TelluricWaveCommand(Dirs, telltarget, exttarget, products["WAVELENGTHPRODUCT"], products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode,False) + verstr
        objcommandsType[tellkey] = True
        ###############################################################
    
        ### Create TARGETS for heliocentric wavelength correction ######
        rvelkey = "RVEL" + basename
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
        objcommands[spckey] = SpcModuleCommand(Dirs, spctarget, Instmode, config, exttarget, products["FLATFLUXCALIBRATIONSPECTRUM"],emptystring, rveltarget, telltarget,products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
        #objdependencies[spckey] = [extkey,rvelkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM",masterfcalkey]
        #objcommands[spckey] = SpcModuleCommand(Dirs, spctarget, Instmode, config, exttarget, products["FLATFLUXCALIBRATIONSPECTRUM"],emptystring, rveltarget,emptystring,products["WAVELENGTHPRODUCT"], objectname, exptime) + verstr
        objcommandsType[spckey] = True

        if objectNameWithoutWhites in listofobjects :
            listofobjectspcs[objectNameWithoutWhites].append(spctarget)
            listofobjectspcsDependencies[objectNameWithoutWhites].append(spckey)
            # -- E Martioli Feb 05 2016 - One should comment lines below if
            #   either one of the RadialVelocity modules can't provide reliable
            #   calculations of source radial velocities for all kinds of sources.
            #   By uncommenting these two lines, it sets dependency of rv measurements
            #   on the stack module. The source RVs are important for the stack module
            #   only if the source's radial velocities are large and variable.
            listofobjectrvs[objectNameWithoutWhites].append(radvel2target)
            listofobjectrvsDependencies[objectNameWithoutWhites].append(radvel2key)
        ###############################################################

        ### Create TARGETS for final FITS spectrum *.fits.gz ##############
        spcfitskey = "OPFITSSPC" + basename
        spcfitstarget = Dirs.PRODUCTDIR + basename + ".m.fits.gz"
        objproducts[spcfitskey] = spcfitstarget
        objdependencies[spcfitskey] = ["WAVELENGTHPRODUCT",spckey,rvelkey,tellkey,snrkey]
        objcommands[spcfitskey] = operaCreateProduct_SPC_Intensity(Dirs, file, spctarget, spcfitstarget, rveltarget, telltarget, snrtarget, products["WAVELENGTHRESOLUTIONPRODUCT"], config.OPERAVERSION, currentTime) + verstr
        objcommandsType[spcfitskey] = True
        ###############################################################

        ###############################################################
        onedspeckey = "ONEDSPEC" + basename
        onedspectarget = Dirs.PRODUCTDIR + basename + ".1d.spc"
        objproducts[onedspeckey] = onedspectarget
        objdependencies[onedspeckey] = [spckey,rvelkey]
        dummylist = [spctarget]
        objcommands[onedspeckey] = StackObjectSpectra(Dirs, onedspectarget, objectNameWithoutWhites, Instmode, dummylist, "", 0, 1) + verstr
        objcommandsType[onedspeckey] = True
        ###############################################################

        ### Create TARGETS for *.fits spectrum ##############
        spcfitskey = "SPCFITS" + basename
        spcfitstarget = Dirs.PRODUCTDIR + basename + "_spc.fits"
        spcfitsbyproduct = Dirs.PRODUCTDIR + basename + "_wavpix_spcfits.txt"
        objproducts[spcfitskey] = spcfitstarget
        objdependencies[spcfitskey] = [spckey, newfilekey ]
        objcommands[spcfitskey] = exportToFits(Dirs, spcfitstarget, spctarget, newfiletarget, spcfitsbyproduct)
        objcommandsType[spcfitskey] = True
        ###############################################################

        ### Create TARGETS for flux calibration Standards ######
        standardname = ""
        if config.hasStandard(objectname.replace(" ", "")) :
            standardname = objectname.replace(" ", "")
        elif config.hasStandard(objectname) :
            standardname = objectname

        if standardname :
            print "Generating FLATRESP for STANDARD:", standardname
            # -- The target below produces a LE-style flat response for each observation of
            #    standard stars. However it could be limited to the Moon spectrum only.
            #    E. Martioli -- Feb 22 2015
            flatrespkey = "FLATRESPONSE" + basename
            flatresptarget = Dirs.PRODUCTDIR + basename + "_flat_resp.fits.gz"
            
            listofstdfcalkey.append(flatrespkey)
            listofstdfcaltargets.append(flatresptarget)
            
            objproducts[flatrespkey] = flatresptarget
            objdependencies[flatrespkey] = [extkey,"WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"]
            objcommands[flatrespkey] = CreateFlatResponseCommand(Dirs, flatresptarget, exttarget, newfiletarget, config.getStandardDataFile(standardname), products["FLATFLUXCALIBRATIONSPECTRUM"], config, Instmode, products["APERTUREPRODUCT"], products["WAVELENGTHPRODUCT"]) + verstr
            objcommandsType[flatrespkey] = True
        ###############################################################

    for object in listofobjects :
        stacknormobjectkey = "STACKNORM" + object
        stacknormobjecttarget = Dirs.PRODUCTDIR + object + "_" + night + "_master.norm.spc"
        objproducts[stacknormobjectkey] = stacknormobjecttarget
        #objdependencies[stacknormobjectkey] = listofobjectspcsDependencies[object] + listofobjectrvsDependencies[object]
        objdependencies[stacknormobjectkey] = listofobjectspcsDependencies[object]
        objcommands[stacknormobjectkey] = StackObjectSpectra(Dirs, stacknormobjecttarget, object, Instmode, listofobjectspcs[object], "", 0, 1) + verstr
        objcommandsType[stacknormobjectkey] = True
        
        stackrawobjectkey = "STACKRAW" + object
        stackrawobjecttarget = Dirs.PRODUCTDIR + object + "_" + night + "_master.raw.spc"
        objproducts[stackrawobjectkey] = stackrawobjecttarget
        #objdependencies[stackrawobjectkey] = listofobjectspcsDependencies[object] + listofobjectrvsDependencies[object]
        objdependencies[stackrawobjectkey] = listofobjectspcsDependencies[object]
        objcommands[stackrawobjectkey] = StackObjectSpectra(Dirs, stackrawobjecttarget, object, Instmode, listofobjectspcs[object], "", 0, 0) + verstr
        objcommandsType[stackrawobjectkey] = True
        
        stackfcalobjectkey = "STACKFCAL" + object
        stackfcalobjecttarget = Dirs.PRODUCTDIR + object + "_" + night + "_master.fcal.spc"
        objproducts[stackfcalobjectkey] = stackfcalobjecttarget
        #objdependencies[stackfcalobjectkey] = listofobjectspcsDependencies[object] + listofobjectrvsDependencies[object]
        objdependencies[stackfcalobjectkey] = listofobjectspcsDependencies[object]
        objcommands[stackfcalobjectkey] = StackObjectSpectra(Dirs, stackfcalobjecttarget, object, Instmode, listofobjectspcs[object], "", 0, 2) + verstr
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
        inputMaskForSPCModule = config.ATYPEWAVELENGTHMASK
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

        objcommands[spckey] += ' --flatResponse=' + masterfcalkeytarget  + ' --inputWavelengthMaskForUncalContinuum=' + inputMaskForSPCModule
    ###############################################################

    return objproducts, objdependencies, objcommands, objcommandsType
#-------------------------------------------

######### PRODUCT FILE NAMES ###########
def setCalibrationProducts(Dirs,night,Instmode,Readmode,DefaultCal,Keywords,allowanyreadout,config,plots,verbose) :

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
        
        products["FITMASTERBIAS"] = INSTCONFIGPREFIX + "_masterbias.fits"
        dependencies["FITMASTERBIAS"] = ["BIASLIST"]
        commands["FITMASTERBIAS"] = MasterCalibrationCommand(Dirs, "operaMedianCombine", products["FITMASTERBIAS"], products["BIASLIST"]) + verstr
        commandsType["FITMASTERBIAS"] = True
        
        products["MASTERBIAS"] = INSTCONFIGPREFIX + "_masterbias.fits.gz"
        dependencies["MASTERBIAS"] = ["FITMASTERBIAS"]
        commands["MASTERBIAS"] = RotateMirrorCropCommand(Dirs,Instmode,products["FITMASTERBIAS"]) + verstr
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
        
        products["FITMASTERFLAT"] = INSTCONFIGPREFIX + "_masterflat.fits"
        dependencies["FITMASTERFLAT"] = ["FLATLIST"]
        commands["FITMASTERFLAT"] = MasterCalibrationCommand(Dirs, "operaMedianCombine", products["FITMASTERFLAT"], products["FLATLIST"]) + verstr
        commandsType["FITMASTERFLAT"] = True
        
        products["MASTERFLAT"] = INSTCONFIGPREFIX + "_masterflat.fits.gz"
        dependencies["MASTERFLAT"] = ["FITMASTERFLAT"]
        commands["MASTERFLAT"] = RotateMirrorCropCommand(Dirs,Instmode,products["FITMASTERFLAT"]) + verstr
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
        
        products["FITMASTERCOMP"] = INSTCONFIGPREFIX + "_mastercomp.fits"
        dependencies["FITMASTERCOMP"] = ["COMPLIST"]
        commands["FITMASTERCOMP"] = MasterCalibrationCommand(Dirs, "operaMedianCombine", products["FITMASTERCOMP"], products["COMPLIST"]) + verstr
        commandsType["FITMASTERCOMP"] = True
        
        products["MASTERCOMP"] = INSTCONFIGPREFIX + "_mastercomp.fits.gz"
        dependencies["MASTERCOMP"] = ["FITMASTERCOMP"]
        commands["MASTERCOMP"] = RotateMirrorCropCommand(Dirs,Instmode,products["FITMASTERCOMP"]) + verstr
        commandsType["MASTERCOMP"] = True
    else :
        products["MASTERCOMP"] = DefaultCal.DEFAULTMASTERCOMP
        dependencies["MASTERCOMP"] = []

    """
    # E. Martioli Mar 24 2015 - Below is commented to use only default gain
    if subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') and \
        subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n') :
        products["GAINPRODUCT"] = INSTCONFIGPREFIX + ".gain.gz"
        dependencies["GAINPRODUCT"] = ["BIASLIST","FLATLIST"]
        commands["GAINPRODUCT"] = GainCommand(Dirs,products["GAINPRODUCT"], products["BIASLIST"], products["FLATLIST"], config.BADPIXELMASK, Readmode) + verstr
        commandsType["GAINPRODUCT"] = True
    else :
    """
    products["GAINPRODUCT"] = DefaultCal.DEFAULTGAINCAL
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
    dependencies["INSTRUMENTPROFILEPRODUCT"] = ["GEOMETRYPRODUCT","GAINPRODUCT","MASTERBIAS","MASTERFLAT","MASTERCOMP"]
    commands["INSTRUMENTPROFILEPRODUCT"] = InstrumentProfileCommand(Dirs,products["INSTRUMENTPROFILEPRODUCT"], products["GEOMETRYPRODUCT"],products["GAINPRODUCT"],products["MASTERBIAS"],products["MASTERFLAT"],products["MASTERCOMP"],config.BADPIXELMASK,Instmode,plots["INSTRUMENTPROFILEPRODUCT"]) + verstr
    commandsType["INSTRUMENTPROFILEPRODUCT"] = True
    
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
    inputwaveguess = ' --inputLineSetFilename=' + Dirs.OESCONFIGDIR+Instmode.WAVE_THARSELECTEDLINES
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

    else :
        plots["ORDERSPACINGPRODUCT"] = {'ORDSPCPLOTFILE': "",'ORDSPCDATAFILE': "",'ORDSPCSCRIPTFILE': ""}
        plots["GEOMETRYPRODUCT"] = {'GEOMPLOTFILE': "",'GEOMDATAFILE': "",'GEOMSCRIPTFILE': ""}
        plots["INSTRUMENTPROFILEPRODUCT"] = {'PROFPLOTFILE': "",'PROFDATAFILE': "",'PROFSCRIPTFILE': ""}
        plots["APERTUREPRODUCT"] = {'APERPLOTFILE': "",'APERDATAFILE': "",'APERSCRIPTFILE': ""}
        plots["WAVELENGTHPRODUCT"] = {'WAVEORDSPLOTFILE': "",'WAVESPECPLOTFILE': "",'WAVESPECSCRIPTFILE': "", \
                                           'WAVEORDSCRIPTFILE': "",'WAVEORDSDATAFILE': "",'WAVEATLASDATAFILE': "", \
                                           'WAVECOMPDATAFILE': "", 'WAVELINESDATAFILE': ""}

    return plots
##################################

############# Class to encapsulate modes for reduction ####################
# This class contains the information on all available 
class ReductionModes :
    'Common base class for reduction modes'

    # modes = {[modekey]: [modeName, InstMode, ReadMode, Nobjects, Nbiases, Nflats, Ncomps] }
    modes = {}
    files = {}
    
    def __init__(self, dirs, keywords, allowanyreadout, forcecalibration):
        if os.path.exists(dirs.DATADIR) :
            self.displayStats = "---\n"
            self.displayStats += "STATISTICS for DATA in DIR: " + dirs.DATADIR + "\n"
            self.displayStats += "Allow any readout? " + str(bool(allowanyreadout)) + "\n"
            self.displayStats += "Include mode for calibration alone? " + str(forcecalibration) + "\n"
            self.displayStats += "---\n"
            self.displayStats += "InstMode\tReadout\tobject\tbias\tflat\tcomp\tSELECTED?\n"
            self.displayStats += "-------------------------------------------------------------------\n"

            for mode in range(1,2) :
                for read in range(1,2) :
                    
                    Instmode = InstMode(mode)
                    Readmode = ReadoutMode(read)
                
                    modeKey = Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME
                    
                    nobjects = 0
                    nbiases = 0
                    nflats = 0
                    ncomps = 0
                    
                    objectcheck = ObjectListCommand(dirs, Instmode, Readmode, keywords.OBJECTKEYWORD, keywords)
                    objects = (subprocess.check_output(objectcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(objects) :
                        nobjects = len(objects)
                    
                    biascheck = testBiasListCommand(dirs, Instmode, Readmode, keywords)
                    biases = (subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(biases) :
                        nbiases = len(biases)
                
                    flatcheck = testCalibrationListCommand(dirs, Instmode, Readmode,  keywords.FLATKEYWORD, keywords, allowanyreadout)
                    flats = (subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(flats) :
                        nflats = len(flats)
            
                    compcheck = testCalibrationListCommand(dirs, Instmode, Readmode,  keywords.COMPKEYWORD, keywords, allowanyreadout)
                    comps = (subprocess.check_output(compcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(comps) :
                        ncomps = len(comps)
                
                    selected = "NO"
                
                    if forcecalibration :
                        if nbiases and nflats and ncomps :
                            self.modes[modeKey] = [modeKey, mode, read, nobjects, nbiases, nflats, ncomps]
                            self.files[modeKey] = [objects, biases, flats, comps]
                            selected = "YES"
                    else :
                        if nobjects and nbiases and nflats and ncomps :
                            self.modes[modeKey] = [modeKey, mode, read, nobjects, nbiases, nflats, ncomps]
                            self.files[modeKey] = [objects, biases, flats, comps]
                            selected = "YES"
    
                    xtab = "\t"
                    self.displayStats += Instmode.INSTRUMENTMODESHORTNAME + xtab + Readmode.READOUTSPEEDSHORTNAME + "\t" + str(nobjects) + "\t" + str(nbiases) + "\t" + str(nflats) + "\t" + str(ncomps) + "\t" + selected + "\n"
    
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
                str(modeItem[1][6]) + " comps  "

    def displayModeStats(self,intrumentmode,readoutspeed) :
        modeItems = self.modes.items()
        for modeItem in modeItems:
            if(modeItem[1][1] == intrumentmode and modeItem[1][2] == readoutspeed) :
                print "MODE SELECTED: " + modeItem[1][0] + "\n" +\
                        str(modeItem[1][3]) + " objects  " +  \
                        str(modeItem[1][4]) + " biases  " +  \
                        str(modeItem[1][5]) + " flats  " + \
                        str(modeItem[1][6]) + " comps  "

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

    def cleanModes(self) :
        for key, item in self.modes.items():
            del self.modes[key]
            del self.files[key]

###############################################################
