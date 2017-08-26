# -*- coding: iso-8859-1 -*-
"""
Created on May 28 2015
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
Last update on Mar 07 2017
"""
import sys,os
import graces

pylibdir = os.path.split(os.path.dirname(__file__))[0] + '/pyLib'
sys.path.insert(0,pylibdir)

from pipelinebrain import Products

############# MAIN FUNCTION TO EXECUTE PIPELINE ###############
def executePipeline(input, Dirs, config, keywords):

    """ 
    First assign input variables
    """
    night = input[0]
    intrumentmode = input[1]
    readoutspeed = input[2]
    clean = input[3]
    simulate = input[4]
    plot = input[5]
    verbose = input[6]
    trace = input[7]
    allowanyreadout = input[8]
    specificProduct = input[9]
    cleanall = input[10]
    
    """
    Set up instrument mode and readout mode
    """
    Instmode = graces.InstMode(intrumentmode)
    Readmode = graces.ReadoutMode(readoutspeed)
    
    """
    Set up default calibrations
    """
    DefaultCal = graces.DefaultCalibration(Dirs, Instmode, Readmode)
    
    """
    CALIBRATION: Create calibration product file names, dependencies and command lines: 
    """
    plotFilenames = graces.setPlotFilenames(Dirs,night,Instmode,Readmode,plot)
    productFilenames,productDependencies,productCommands,productCommandTypes = graces.setCalibrationProducts(Dirs,night,Instmode,Readmode,DefaultCal,keywords,allowanyreadout,config,plotFilenames,verbose)
    
    """
    REDUCTION: Create object product file names, dependencies and command lines:
    """
    objectproducts = graces.setObjectProducts(productFilenames, Dirs, night, Instmode, Readmode, DefaultCal, keywords, config, verbose)

    """
    Configure calibration and object product targets:
    """
    products = Products(productFilenames,plotFilenames,productDependencies,productCommands,productCommandTypes,trace)
    products.addTargets(objectproducts)

    """
    Activate simulation mode
    """    
    if (simulate) :
        print "\n--- START SIMULATION ---\n"
        products.setSimulation()
    
    
    """
    ****** RUN PIPELINE ********
    """
    if (cleanall) :
        if(verbose and bool(simulate) == False) :
            print "Removing the following products: "
            products.displayTargets()
        
        products.cleanAll()
    else :
        #products.executeTarget("INSTRUMENTPROFILEPRODUCT") # This is an example how to execute a given target
        # Note that any target will trigger a novel process until the final product can be produced
            
        if (specificProduct in "OBJECTS") :
            products.executeTargetsWithSubstrings(["OPSPC"])
        elif (specificProduct in "CALIBRATIONS") :
            products.executeTargetsWithSubstrings(["GEOMETRY","WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"])
        elif (specificProduct == "") :
            products.executeAllTargets()
        else :
            products.executeTargetsWithSubstrings([specificProduct])

        if clean :
            if (specificProduct in "LIBRE-ESPRIT") :
                products.removeTargets(["LESPC","LEPOL"])
            elif (specificProduct in "OBJECTS") :
                products.removeTargets(["OPSPC","OPPOL"])
            elif (specificProduct in "CALIBRATIONS") :
                products.removeTargets(["GEOMETRY","WAVELENGTHPRODUCT","FLATFLUXCALIBRATIONSPECTRUM"])
            elif (specificProduct != "") :
                products.removeTargets([specificProduct])
            
            if(verbose and bool(simulate) == False) :
                print "Removing the following products: "
                products.displayTargets()

            products.cleanAll()
    #---------------------------
    
    """
    Reset simulation mode
    """
    if (simulate) :
        products.resetSimulation()
        print "\n --- END SIMULATION --- \n"

############# END ####################


