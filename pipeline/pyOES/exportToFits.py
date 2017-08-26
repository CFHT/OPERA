#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

"""
    --> #!/usr/bin/python
    Created on Nov 17 2014

    Description: A module to plot an opera product spectrum.
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python exportToFits.py --spectrumfile=inputfilename.spc.gz --originalfitsfile=inputfilename.fits.gz --wavePixfile=inputfile.txt --outputfitsfilename=outputfilename.fits 
"""

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import sys,os
import spectralclasses_polyn
import numpy as np

parser = OptionParser()
parser.add_option("-S", "--spectrumfile", dest="spectrumfile", help='spectrum file',type='string',default="")
#parser.add_option("-W", "--wcalfile", dest="wcalfile", help='wcal file', type='string', default="")
parser.add_option("-I", "--originalfitsfile", dest="originalfitsfile", help='original fits file', type='string',default="")
parser.add_option("-O", "--outputfitsfilename", dest="outputfitsfilename", help='output fits filename', type='string',default="")
parser.add_option("-W", "--wavePixfile", dest="wavePixfile", help='wavelength-pixel filename', type='string',default="")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with plotSpectrum.py -h ";sys.exit(1);

if options.verbose:
    print 'Spectrum file: ', options.spectrumfile
    print 'Original spectrum: ', options.originalfitsfile
    print 'Output spectrum: ', options.outputfitsfilename
    print 'Wavelength-pixel textfile: ', options.wavePixfile


spc = spectralclasses_polyn.Spectrum(options.spectrumfile)

spectrumsize=2048

spc.exportToFits(options.originalfitsfile, options.outputfitsfilename, spectrumsize,options.wavePixfile)
