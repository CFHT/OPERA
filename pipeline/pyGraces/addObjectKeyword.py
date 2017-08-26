#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
"""
Created on Jul 23 2015
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
"""
import sys, getopt
from astropy.io import fits
from datetime import datetime

def main(argv):
    inputlist = ''
    objectname = ''
    
    verbose = 0
    try:
        opts, args = getopt.getopt(argv,"hvi:n:",["inputlist=","objectname="])
    except getopt.GetoptError:
        print 'addObjectKeyword.py -i <inputlist> -n <objectname>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'addObjectKeyword.py -i <inputlist> -n <objectname>\n'
            print 'Available options are:'            
            print '-i , --inputlist=: "listoffiles.txt" '
            print '-n , --objectname=: "objectname" \n'
            sys.exit()
        elif opt in ("-i", "--inputlist"):
            inputlist = arg
        elif opt in ("-n", "--objectname"):
            objectname = arg
        elif opt in ("-v", "--verbose"):
            verbose = 1

    if verbose:
        print 'inputlist: "', inputlist
        print 'objectname: "', objectname

    updatecomment = 'by setGracesKeywords ' + datetime.now().strftime('%Y-%m-%d')   
    
    listfile = open(inputlist, 'r')
    
    for image in listfile:
        print "Processing image: " + image.rstrip('\n')
        hdulist = fits.open(image.rstrip('\n'), mode='update')
        prihdr = hdulist[0].header
        
        prihdr.set('OBJECT', objectname)

        hdulist.flush()
        hdulist.close()
        
    listfile.close()
        
if __name__ == "__main__":
    main(sys.argv[1:])

    
