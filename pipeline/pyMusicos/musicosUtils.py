# -*- coding: utf-8 -*-
"""
Created on Mar 26 2015
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on Mar 26 2015
"""

import os
from numpy import *
from astropy.io.fits import getheader
import astropy.io.fits as fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle

######################
def get_fitsfzfilepaths(directory):
	"""
        This function generates a list of file names in a directory
        tree by walking the tree either top-down or bottom-up. For each
        directory in the tree rooted at directory top (including top itself),
        it yields a 3-tuple (dirpath, dirnames, filenames).
        """
	file_paths = []  # List which will store all of the full filepaths.
    
	# Walk the tree.
	for root, directories, files in os.walk(directory):
		for filename in files:
			# Merge strings to form full filepaths
			filepath = os.path.join(root, filename)
			if filename.endswith(".fits.fz") :
				file_paths.append(filepath)
    
	return file_paths  # Self-explanatory.
######################

######################
def get_fitsgzfilepaths(directory):
	"""
        This function generates a list of file names in a directory
        tree by walking the tree either top-down or bottom-up. For each
        directory in the tree rooted at directory top (including top itself),
        it yields a 3-tuple (dirpath, dirnames, filenames).
        """
	file_paths = []  # List which will store all of the full filepaths.
    
	# Walk the tree.
	for root, directories, files in os.walk(directory):
		for filename in files:
			# Merge strings to form full filepaths
			filepath = os.path.join(root, filename)
			if filename.endswith(".fits.gz") :
				file_paths.append(filepath)
    
	return file_paths  # Self-explanatory.
######################

######################
def get_fitsfilepaths(directory):
	"""
	This function generates a list of file names in a directory
	tree by walking the tree either top-down or bottom-up. For each
	directory in the tree rooted at directory top (including top itself),
	it yields a 3-tuple (dirpath, dirnames, filenames).
	"""
	file_paths = []  # List which will store all of the full filepaths.

	# Walk the tree.
	for root, directories, files in os.walk(directory):
		for filename in files:
			# Merge strings to form full filepaths
			filepath = os.path.join(root, filename)
			if filename.endswith(".fits") :
				file_paths.append(filepath)

	return file_paths  # Self-explanatory.
######################

######################
def selectObjects(file_paths) :
    """
    This function selects images from a list for which OBSTYPE='OBJECT'. 
    Output is an array of file paths.
    """
    obstypekeyvalue = 'OBJECT'
    objectlist = []

    for file in file_paths :
        header = getheader(file.rstrip('\n'), 0)
        if(header['OBSTYPE'] == obstypekeyvalue) :
            objectlist.append(file.rstrip('\n'))
    return objectlist
######################

###############################################
def selectObjectsToFile(file_paths, objectname, objectlistfilename) :

    """
        This function generates a list of object file names.
        If option savetofile is true, then it will create an file with a list of object file paths
        and it will return the list file name. If option savetofile is false then it returns an
        array of object file paths.
    """

    obstypekeyvalue = 'OBJECT'
    objectlist = []
    
    for file in file_paths:
        header = getheader(file.rstrip('\n'), 0)
        if(header['OBSTYPE'] == obstypekeyvalue and header['OBJECT'] == objectname) :
            objectlist.append(file.rstrip('\n'))

    if (objectlistfilename) :
        listfile = open(objectlistfilename, 'w')
        for object in objectlist:
            listfile.write(object+'\n')
        listfile.close()

    return objectlist
###############################################

########### EXPORT PRIMARY HEADER EXTENSION ##############
def	exportPrimaryHeaderExtension(inputImage, outputHeaderExtension, verbose) :
    if inputImage.endswith(".fits") and os.path.exists(inputImage) :
        if verbose :
            print "Exporting primary header from image: " + inputImage

        hdulist = fits.open(inputImage)
        hdulist[0].writeto(outputHeaderExtension)
        hdulist.close()
    else :
        print "Input image: ", inputImage, " is not FITS or does not exist!"
        exit()
###########################################

########### PRINT IMAGE INFO ##############
def	printImageInfo(inputImage) :
    if inputImage.endswith(".fits") and os.path.exists(inputImage) :
                
        hdulist = fits.open(inputImage)
        hdulist.info()
        print "\n\n"
        i = 0
        for entry in hdulist[0].header :
            print entry, " = ", hdulist[0].header[i]
            i = i+1
###########################################

############# GUNZIP OR FUNPACK FILES ###############
def UnpackFilesIfNeeded(datadir, funpackEXE, gunzipEXE, verbose) :
    
    for file in get_fitsgzfilepaths(datadir) :
        commandgunzip = gunzipEXE + " " + file
        try :
            if verbose :
                print "Executing command: ",commandgunzip
            os.system(commandgunzip)
        except :
            print "Error: can\'t execute command: ",commandgunzip
            exit()

    for file in get_fitsfzfilepaths(datadir) :
        commandfunpack = funpackEXE + " " + file
        try :
            if verbose :
                print "Executing command: ",commandfunpack
            os.system(commandfunpack)
        except :
            print "Error: can\'t execute command: ",commandfunpack
            exit()
###########################################


############# FIX HEADERS OF DATA FILES ###############
def fixMusicosDataHeaders(datadir, outputdir, verbose, simulate) :

    update = False
    if outputdir == datadir or outputdir == "":
        update = True

    filelist = get_fitsgzfilepaths(datadir)
    filelist += get_fitsfilepaths(datadir)
    
    for file in filelist :
        
        if update==True :
            hdulist = fits.open(file, mode='update')
        else :
            hdulist = fits.open(file)

        prihdr = hdulist[0].header

        # Check whether it is MUSICOS data
        isitMusicos = False
        instrumenkeyvalue = prihdr['INSTRUME']

        if "MUSICOS" in instrumenkeyvalue or "Musicos" in instrumenkeyvalue or "musicos" in instrumenkeyvalue:
            isitMusicos = True
            instrumekeyvalue = "MUSICOS"
            prihdr['INSTRUME'] = (instrumekeyvalue, 'instrument name')

        # ---

        if ('MODDATA' not in prihdr) and isitMusicos:

            # Figure out OBSTYPE (bias, flat, comp, or object)
            obstypekeyvalue = "INDEF"
        
            imagekeyvalue = hdulist[0].header['IMAGE']
            objectkeyvalue = hdulist[0].header['OBJECT']
            exptimekeyvalue = hdulist[0].header['EXPTIME']
            expvalue = float(exptimekeyvalue.replace(',', '.'))
            prihdr['EXPTIME2'] = (expvalue, 'Exposure time inserted by OPERA')

            readtimekeyvalue = '%.1e' % float(hdulist[0].header['READTIME'])
            prihdr['READTIME'] = (readtimekeyvalue, 'Read time inserted by OPERA')

            if ("bias" in imagekeyvalue or "Bias" in imagekeyvalue or "BIAS" in imagekeyvalue or \
                "zero" in imagekeyvalue or "Zero" in imagekeyvalue or "ZERO" in imagekeyvalue or \
                "bias" in objectkeyvalue or "Bias" in objectkeyvalue or "BIAS" in objectkeyvalue or \
                "zero" in objectkeyvalue or "Zero" in objectkeyvalue or "ZERO" in objectkeyvalue) and \
                expvalue <= 0.0001 :
            
                obstypekeyvalue = "BIAS"

            elif ("flat" in imagekeyvalue or "Flat" in imagekeyvalue or "FLAT" in imagekeyvalue or \
                  "flat" in objectkeyvalue or "Flat" in objectkeyvalue or "FLAT" in objectkeyvalue) and \
                    expvalue > 0.0001 :
            
                obstypekeyvalue = "FLAT"

            elif ("thar" in imagekeyvalue or "ThAr" in imagekeyvalue or "THAR" in imagekeyvalue or \
                  "thar" in objectkeyvalue or "ThAr" in objectkeyvalue or "THAR" in objectkeyvalue or \
                  "comp" in imagekeyvalue or "Comp" in imagekeyvalue or "COMP" in imagekeyvalue or \
                  "comp" in objectkeyvalue or "Comp" in objectkeyvalue or "COMP" in objectkeyvalue or \
                  "arc" in imagekeyvalue or "Arc" in imagekeyvalue or "ARC" in imagekeyvalue or \
                  "arc" in objectkeyvalue or "Arc" in objectkeyvalue or "ARC" in objectkeyvalue) and \
                    expvalue > 0.0001 :
            
                obstypekeyvalue = "COMP"

            else :
                obstypekeyvalue = "OBJECT"

            #--

            prihdr['OBSTYPE'] = (obstypekeyvalue, 'BIAS, FLAT, COMP, or OBJECT')

            # Figure out INSTMODE (red or blue)
            instmodekeyvalue = "INDEF"

            if ("RED" in imagekeyvalue or "Red" in imagekeyvalue or "red" in imagekeyvalue or \
                "_R_" in imagekeyvalue or "_r_" in imagekeyvalue or "_R" in imagekeyvalue) or \
                ("RED" in objectkeyvalue or "Red" in objectkeyvalue or "red" in objectkeyvalue or \
                 "_R" in objectkeyvalue or "_r_" in objectkeyvalue or "_R" in objectkeyvalue) :
        
                instmodekeyvalue = "RED"

            elif ("BLUE" in imagekeyvalue or "Blue" in imagekeyvalue or "blue" in imagekeyvalue or \
                  "_B_" in imagekeyvalue or "_b_" in imagekeyvalue or "_B" in imagekeyvalue) or \
                ("BLUE" in objectkeyvalue or "Blue" in objectkeyvalue or "blue" in objectkeyvalue or \
                 "_B" in objectkeyvalue or "_b_" in objectkeyvalue or "_B" in objectkeyvalue) :

                instmodekeyvalue = "BLUE"
                    #--

            prihdr['INSTMODE'] = (instmodekeyvalue, 'prism configuration: RED or BLUE')

            if obstypekeyvalue == 'OBJECT' :
                ra_hex = hdulist[0].header['RA']
                dec_hex = hdulist[0].header['DEC']
                c = SkyCoord(ra_hex+' '+dec_hex, unit=(u.hourangle, u.deg))

                prihdr['RA_DEG'] = (round(c.ra.degree,5), 'RA in degrees')
                prihdr['DEC_DEG'] = (round(c.dec.degree,5), 'Dec in degrees')

            julianDate = 2400000.5
            if obstypekeyvalue == 'OBJECT' :
                julianDate = hdulist[0].header['JD']
            mjd = 0
            if julianDate :
                mjd = float(julianDate) - 2400000.5

            prihdr['MJD'] = (round(mjd,5), 'MJD = JD - 2400000.5')

            prihdr['MODDATA'] = (True, 'Header modified by OPERA pipeline')

            if update==True :
                if simulate :
                    print "Modifying original file:",file, " OBSTYPE="+obstypekeyvalue, "INSTMODE="+instmodekeyvalue
                else :
                    if verbose :
                        print "Modifying original file:",file, " OBSTYPE="+obstypekeyvalue, "INSTMODE="+instmodekeyvalue
#                    os.remove(file)
#                    hdulist.writeto(file)
                    hdulist.flush()
            else :
                outputFile = outputdir + os.path.basename(file)
                if simulate:
                    print "Creating a new file:",outputFile, " OBSTYPE="+obstypekeyvalue, "INSTMODE="+instmodekeyvalue
                else :
                    if os.path.exists(outputFile) :
                        if verbose :
                            print "Skipping.. : ",outputFile, "already exists"
                    else :
                        if verbose :
                            print "Creating a new file:",outputFile, " OBSTYPE="+obstypekeyvalue, "INSTMODE="+instmodekeyvalue
                        hdulist.writeto(outputFile)
    
        else:
            if 'MODDATA' in prihdr and verbose :
                print "Skipping header modification: MODDATA keyword already exists"
            if not isitMusicos and verbose:
                print "Skipping header modification: INSTRUMEN keyword is not MUSICOS"


        hdulist.close()
###########################################

############# FIX HEADERS OF DATA FILES ###############
def createMusicosBadPixelMap(inputReffile, outputfile) :

    #inputReffile = "/data/MUSICOS/14set05_mod/flat_R_000.fits"
    #outputfile = "/Users/edermartioli/Reductions/MUSICOS/14set05_mod/MUSICOS_badpixelmask.fits"

    hdulist = fits.open(inputReffile)
    data = hdulist[0].data
    for pixIndex in range(len(data)) :
        data[pixIndex] = 1.0

    prihdr = hdulist[0].header

    prihdr['OBSTYPE'] = ("BADPIX", 'Badpixel mask')
    prihdr['MODDATA'] = ("TRUE", 'Header modified by OPERA pipeline')
    prihdr['INSTMODE'] = ("INDEF", 'prism configuration: RED or BLUE')

    prihdr['IMAGE'] = ("MUSICOS_badpixelmask.fits", 'Header modified by OPERA pipeline')
    prihdr['OBJECT'] = ("BADPIX", 'Header modified by OPERA pipeline')
    prihdr['EXPTIME'] = ("INDEF", 'Header modified by OPERA pipeline')

    hdulist.writeto(outputfile)
    hdulist.close()

###########################################

