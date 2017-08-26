# -*- coding: utf-8 -*-
"""
Created on May 28 2015
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on May 28 2015
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
def get_allfitsfilepaths(directory):
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
            if filename.endswith(".fits") or filename.endswith(".fits.gz"):
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

###############################################
def getcompkey(Dirs, Keywords, Instmode, Readmode) :
    
    file_paths = get_fitsfilepaths(Dirs.DATADIR)
    
    selectedkey = ""
    
    for key in Keywords.COMPKEYWORD :
        complist = []
    
        for file in file_paths:
            header = getheader(file.rstrip('\n'), 0)
            if header[Keywords.OBSTYPEKEY] == key and header[Keywords.INSTRUMEKEY] == Instmode.INSTRUME and header[Keywords.INSTMODEKEY] == Instmode.INSTRUMENTMODEKEY and header[Keywords.READMODEKEY] == Readmode.READOUTSPEED :
                complist.append(file.rstrip('\n'))

        if len(complist) :
            selectedkey = key

    if selectedkey == "" :
        selectedkey="COMP"
    
    return selectedkey
###############################################

