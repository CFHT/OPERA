# -*- coding: utf-8 -*-
"""
Created on Mon Jan 09, 2017
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
"""

import os
from astropy.io.fits import getheader
import astropy.io.fits as fits

############# Class to encapsulate instrument modes for reduction ###############
# This class contains the information on all available reduction modes
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
            
            for grating in range(1,3) :
                for occmask in range(1,3) :
                    
                    Instmode = InstMode(mode)
                    Readmode = ReadoutMode(read)
                    
                    modeKey = Instmode.INSTRUMENTMODESHORTNAME + "_" + Readmode.READOUTSPEEDSHORTNAME
                    
                    nobjects = 0
                    ndarks = 0
                    nflats = 0
                    nronchis = 0
                    narcs = 0
                    
                    objectcheck = ObjectListCommand(dirs, Instmode, Readmode, keywords.OBJECTKEYWORD)
                    objects = (subprocess.check_output(objectcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(objects) :
                        nobjects = len(objects)
                    
                    biascheck = testBiasListCommand(dirs, Instmode, Readmode, keywords)
                    biases = (subprocess.check_output(biascheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(biases) :
                        nbiases = len(biases)
                    
                    flatcheck = testCalibrationListCommand(dirs, Instmode, Readmode,  keywords.FLATKEYWORD, allowanyreadout)
                    flats = (subprocess.check_output(flatcheck,stderr=subprocess.STDOUT,shell=True).rstrip('\n')).split()
                    if len(flats) :
                        nflats = len(flats)
                    
                    compcheck = testCalibrationListCommand(dirs, Instmode, Readmode,  keywords.COMPKEYWORD, allowanyreadout)
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
            exit()

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
