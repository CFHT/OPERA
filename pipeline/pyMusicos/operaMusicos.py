#!/opt/anaconda/bin/python
# -*- coding: iso-8859-1 -*-
"""
    *** IMPORTANT NOTE ***
    Use line below as shebang for default python location
    #!/usr/bin/python
    
    Use a line similar to the ones below as shebang for custom python location
    #!/$HOME/Ureka/variants/common/bin/python
    #!/opt/anaconda/bin/python
   
    Created on Mar 26 2015
	
    Description: A wrapper to run Opera Musicos reduction pipeline.
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    /Users/edermartioli/opera-1.0/pipeline/pyMusicos/operaMusicos.py --datarootdir=/data/MUSICOS/ 
    --pipelinehomedir=/Users/edermartioli/opera-1.0 --productrootdir=/Users/edermartioli/Reductions/MUSICOS/ 
    --night=14set05 --product="CALIB" -pvts
"""

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import sys,os
import musicospipeline
import musicos
import musicosUtils

parser = OptionParser()
parser.add_option("-N", "--night", dest="night", help="night directory",type='string',default="")
parser.add_option("-D", "--datarootdir", dest="datarootdir", help="data root directory",type='string',default="/data/MUSICOS/")
parser.add_option("-O", "--pipelinehomedir", dest="pipelinehomedir", help="pipeline directory",type='string',default="/Users/edermartioli/opera-1.0/")
parser.add_option("-P", "--productrootdir", dest="productrootdir", help="data product root directory",type='string',default="/Users/edermartioli/Reductions/MUSICOS/")
parser.add_option("-T", "--product", dest="product", help='target product: "CALIBRATIONS", "OBJECTS" (default)"',type='string',default="OBJECTS")
parser.add_option("-a", action="store_true", dest="cleanall", help="JUST clean all products",default=False)
parser.add_option("-c", action="store_true", dest="clean", help="clean products",default=False)
parser.add_option("-s", action="store_true", dest="simulate", help="simulate",default=False)
parser.add_option("-p", "--plot", action="store_true", dest="plot", help="plots",default=False)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="verbose",default=False)
parser.add_option("-t", "--trace", action="store_true", dest="trace", help="trace",default=False)
parser.add_option("-m", "--moddatadir", action="store_true", dest="moddatadir", help="moddatadir",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with operaMusicos.py -h ";sys.exit(1);

if options.verbose:
    print 'PIPELINE HOME DIR: ', options.pipelinehomedir
    print 'DATA ROOT DIR: ', options.datarootdir
    print 'PRODUCT ROOT DIR: ', options.productrootdir
    print 'NIGHT: ', options.night

"""
Set up directories:
"""
Dirs = musicos.Directories(options.pipelinehomedir,options.datarootdir,options.productrootdir,options.night)

"""
Fix header keywords:
"""
if(options.moddatadir) :
    Dirs.createModDataDir(options.simulate)
    musicosUtils.fixMusicosDataHeaders(Dirs.DATADIR,Dirs.DATADIRMOD,options.verbose,options.simulate)
    Dirs.DATADIR = Dirs.DATADIRMOD
else :
    musicosUtils.fixMusicosDataHeaders(Dirs.DATADIR,"",options.verbose,options.simulate)

"""
Set up config files:
"""
config = musicos.ConfigFiles(Dirs)

"""
Set up MUSICOS keywords:
"""
keywords = musicos.Keywords()

"""
Set up available modes for reduction:
"""
allowanyreadout = False

forcecalibration = False # This is set to "True" for calibrations even when there is no object files
if(options.product == "CALIBRATIONS") :
    forcecalibration = True

modes = musicos.ReductionModes(Dirs, keywords, allowanyreadout, forcecalibration)

modes.displayOverallStats()

for mode in modes.getInstReadModes() :
    intrumentmode = mode[0]
    readoutspeed = mode[1]
    modes.displayModeStats(intrumentmode,readoutspeed)
    input = [options.night,intrumentmode,readoutspeed,options.clean,options.simulate,options.plot,options.verbose,options.trace,allowanyreadout,options.product,options.cleanall]
    musicospipeline.executePipeline(input, Dirs, config, keywords)

modes.cleanModes()



