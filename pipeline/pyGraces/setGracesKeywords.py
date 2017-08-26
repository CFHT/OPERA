#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
"""
Created on Jul 27 2014
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
"""
import sys, getopt
from astropy.io import fits
from datetime import datetime

def main(argv):
    inputlist = ''
    instmode = ''            # options available are: staronly, starsky 
    readspeed = ''           # options available are: fast, normal, slow
    obstype = ''             # options available are: bias, flat, comp, object
    
    verbose = 0
    try:
        opts, args = getopt.getopt(argv,"hvi:m:r:t:",["inputlist=","instmode=","readspeed=", "obstype="])
    except getopt.GetoptError:
        print 'setGracesKeywords.py -i <inputlist> -m <instmode> -r <readspeed> -t <obstype>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'setGracesKeywords.py -i <inputlist> -m <instmode> -r <readspeed> -t <obstype>\n'
            print 'Available options are:'            
            print '-i , --inputlist=: "listoffiles.txt" '
            print '-m , --instmode=: "staronly" or "starsky" '
            print '-r , --readspeed=: "fast" or "normal" or "slow" '
            print '-t , --obstype=: "bias" or "flat" or "comp" or "object" \n'
            sys.exit()
        elif opt in ("-i", "--inputlist"):
            inputlist = arg
        elif opt in ("-m", "--instmode"):
            instmode = arg
        elif opt in ("-r", "--readspeed"):
            readspeed = arg            
        elif opt in ("-t", "--obstype"):
            obstype = arg
        elif opt in ("-v", "--verbose"):
            verbose = 1

    if verbose:
        print 'inputlist: "', inputlist
        print 'instmode: "', instmode
        print 'readspeed: "', readspeed
        print 'obstype: "', obstype
    
    
    instmodekeyvalue = ''
    readspeedkeyvalue = ''
    obstypekeyvalue = ''
    instrumentkeyvalue = 'GRACES'    
    
    # Below it sets the header keyword value for selected instrument mode    
    if (instmode == "staronly") :
        instmodekeyvalue = "FOURSLICE"
    elif (instmode == "starsky") :
        instmodekeyvalue = "TWOSLICE"
         
    # Below it sets the header keyword value for selected readout mode             
    if (readspeed == "fast") :
        readspeedkeyvalue = "Fast: 4.70e noise, 1.60e/ADU, 32s"
    elif (readspeed == "normal") :
        readspeedkeyvalue = "Normal: 4.20e noise, 1.30e/ADU, 38s"       
    elif (readspeed == "slow") :
        readspeedkeyvalue = "Slow: 2.90e noise, 1.20e/ADU, 60s"         
    
    # Below it sets the header keyword value for selected obstype             
    if (obstype == "bias") :
        obstypekeyvalue = "BIAS"
    elif (obstype == "flat") :
        obstypekeyvalue = "FLAT"
    elif (obstype == "comp") :
        obstypekeyvalue = "COMP"
    elif (obstype == "object") :
        obstypekeyvalue = "OBJECT"
    
    if verbose:    
        print "The following keywords are being updated: "
        print "INSTRUME: " + instrumentkeyvalue
        print "GSLICER: " + instmodekeyvalue
        print "EREADSPD: " + readspeedkeyvalue
        print "OBSTYPE: " + obstypekeyvalue
    
    updatecomment = 'by setGracesKeywords ' + datetime.now().strftime('%Y-%m-%d')   
    
    listfile = open(inputlist, 'r')
    
    for image in listfile:
        print "Processing image: " + image.rstrip('\n')
        hdulist = fits.open(image.rstrip('\n'), mode='update')
        prihdr = hdulist[0].header
                      
        if (instmode) :
            try:
                prihdr['INSTRUME'] = (instrumentkeyvalue, updatecomment)
            except:
                prihdr.set('INSTRUME', instrumentkeyvalue)            
                
            try:
                prihdr['GSLICER'] = (instmodekeyvalue, updatecomment)
            except:
                prihdr.set('GSLICER', instmodekeyvalue)
            
        if (readspeed) : 
            try:
                prihdr['EREADSPD'] = (readspeedkeyvalue, updatecomment)
            except:
                prihdr.set('EREADSPD', readspeedkeyvalue)

        if (obstype) :
            try:
                prihdr['OBSTYPE'] = (obstypekeyvalue, updatecomment)
            except:
                prihdr.set('OBSTYPE', obstypekeyvalue)  
                
        hdulist.flush()
        hdulist.close()
        
    listfile.close()
        
if __name__ == "__main__":
    main(sys.argv[1:])

    
