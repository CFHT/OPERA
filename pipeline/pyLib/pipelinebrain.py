# -*- coding: utf-8 -*-
"""
    Created on Dec 01 2014
    @author: Eder Martioli
    Laboratorio Nacional de Astrofisica, Brazil
    Last update on Dec 01, 2014
    """
import sys, os
import subprocess

######### PRODUCTS ###########
class Products :
    'Common base class for pipeline products'
    targets = {}
    dependencies = {}
    commands = {}
    commandstype = {}
    
    plots = {}
    existenceStatus = {}
    
    trace = False
    simulation = False
    
    def __init__(self,Targets, Plots, Dependencies, Commands, CommandsType, Trace):
        self.targets = Targets
        self.plots = Plots
        self.dependencies = Dependencies
        self.commands = Commands
        self.commandstype = CommandsType
        self.setExistenceStatus()
        self.trace = Trace
    
    # Function below sets processing in simulation mode
    def setSimulation(self) :
        self.simulation = True
    #-------------------------------------------
    
    # Function below resets simulation back to non-simulation mode
    def resetSimulation(self) :
        self.simulation = False
        self.setExistenceStatus()
    #-------------------------------------------
    
    # Function below returns the product file name <string> associated to a target
    def getTarget(self, key) :
        return self.targets[key]
    #-------------------------------------------
    
    # Function to display targets
    def displayTargets(self) :
        for targetItem in self.targets.items():
            print '"',targetItem[0],'" :', targetItem[1]
    #-------------------------------------------
    
    # Function to add targets from input dictionary.
    # Input is a vector of four dictionaries: target{}, dependencies{}, commands{}, commandstype{}
    def addTargets(self, inputDictVec) :
        # input: inputDictVec = [target{}, dependencies{}, commands{}, commandstype{}]
        
        self.targets.update(inputDictVec[0])
        
        newkeys = []
        for item in inputDictVec[0].items():
            newkeys.append(item[0])
    
        dictbool = dict.fromkeys(newkeys, False)
        self.existenceStatus.update(dictbool)
    
        self.dependencies.update(inputDictVec[1])
        self.commands.update(inputDictVec[2])
        self.commandstype.update(inputDictVec[3])
    
        self.setExistenceStatus()
    #-------------------------------------------
    
    # Function below returns the list of dependencies to produce product file
    def getDependencies(self, key) :
        return self.dependencies[key]
    #-------------------------------------------
    
    # Function below returns the command line to produce product file
    def getCommandLine(self, key) :
        return self.commands[key]
    #-------------------------------------------
    
    # Function below returns the command line to produce product file
    def getCommandType(self, key) :
        return self.commandstype[key]
    #-------------------------------------------
        
    # Function below executes a single target
    def executeTarget(self, key) :
        if (self.simulation == False) :
            self.setExistenceStatus()
        # If target exists then set target existenceStatus to True and exit
        if self.existenceStatus[key] == True :
            return True
        # Otherwise execute command and recheck target existence and exit
        else :
            if self.resolveDependencies(key) == True :
                if (self.simulation == True) :
                    if key in self.commands :
                        if self.getCommandType(key) :
                            print self.getCommandLine(key)
                        else :
                            print "python: " + self.getCommandLine(key)
                        self.existenceStatus[key] = True
                else :
                    try :
                        if(self.trace == True) :
                            print self.getCommandLine(key)
                    
                        if key in self.commands :
                            if self.getCommandType(key) :
                                #subprocess.check_output(self.getCommandLine(key),stderr=subprocess.STDOUT,shell=True)
                                os.system(self.getCommandLine(key))
                            else :
                                exec self.getCommandLine(key)
                    
                        if os.path.exists(self.targets[key]) :
                            self.existenceStatus[key] = True
                    except :
                        print "Error: can\'t execute command: ",self.getCommandLine(key)

        return self.getExistenceStatus(key)
    #-------------------------------------------
    
    # Function below executes all targets matching any substring given in the input list of substrings
    def executeTargetsWithSubstrings(self,substrings) :
        self.setExistenceStatus()
        targetItems = self.targets.items()
        for targetItem in targetItems :
            for str in substrings :
                if (str in targetItem[0]) or (str in targetItem[1]) :
                    self.executeTarget(targetItem[0])
        
        if (self.simulation == True) :
            self.setExistenceStatus()

    # Function below returns the status of a target <bool>
    def getExistenceStatus(self, key) :
        return self.existenceStatus[key]
    #-------------------------------------------

    # Function below sets the status of each target based on existence check
    def setExistenceStatus(self) :
        # Convert to list of tuples
        targetItems = self.targets.items()
        for targetItem in targetItems:
            if os.path.exists(targetItem[1]) :
                self.existenceStatus[targetItem[0]] = True
            else :
                self.existenceStatus[targetItem[0]] = False
    #-------------------------------------------

    # Function below execute every step where ExistenceStatus is FALSE
    def executeAllTargets(self) :
        self.setExistenceStatus()
        targetItems = self.targets.items()
        for targetItem in targetItems:
            if (self.existenceStatus[targetItem[0]] == False) :
                self.executeTarget(targetItem[0])
        if (self.simulation == True) :
            self.setExistenceStatus()
    #-------------------------------------------

    # Function to check dependencies.
    # Return False if any dependency is missing.
    def checkDependencies(self, key) :
        dependenciesstatus = True
        for dkey in self.dependencies[key] :
            if self.existenceStatus[dkey] == False :
                dependenciesstatus = False
        return dependenciesstatus
    #-------------------------------------------

    # Function to resolve dependencies
    def resolveDependencies(self, key) :
        if (self.checkDependencies(key) == False) :
            for dkey in self.dependencies[key] :
                if self.existenceStatus[dkey] == False :
                    self.executeTarget(dkey)
        
        return self.checkDependencies(key)
    #-------------------------------------------

    # Function to delete all targets
    def cleanTargets(self) :
        targetItems = self.targets.items()
        for targetItem in targetItems:
            if os.path.exists(targetItem[1]) :
                try :
                    if(self.simulation == True) :
                        print "rm " + targetItem[1]
                    else :
                        os.remove(targetItem[1])
                except:
                    print "Error: can\'t remove file: ", targetItem[1]
                else:
                    self.existenceStatus[targetItem[0]] = False
        
        return
    #-------------------------------------------

    # Function to delete all plots
    def cleanPlots(self) :
        for plotItem in self.plots.items():
            for plotItemItem in plotItem[1].items() :
                if os.path.exists(plotItemItem[1]) :
                    try :
                        if(self.simulation == True) :
                            print "rm " + plotItemItem[1]
                        else :
                            os.remove(plotItemItem[1])
                    except:
                        print "Error: can\'t remove file: ", plotItemItem[1]
        return
    #-------------------------------------------

    # Function to delete all targets and plots
    def cleanAll(self) :
        self.cleanPlots()
        self.cleanTargets()
        return
    #-------------------------------------------

    # Function to remove targets matching any substring given in the input list of substrings
    def removeTargets(self,substrings) :
        targetItems = self.targets.items()
        for targetItem in targetItems :
            for str in substrings :
                if (str in targetItem[0]) or (str in targetItem[1]) :
                    self.targets.pop(targetItem[0], None)
        return
#-------------------------------------------

##################################
