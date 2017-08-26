#########################################################################################
#
# Makefile name: Makefile
# Version: 1.0
# Description: main harness entry point - includes all other packages
# Author(s): CFHT OPERA team
# Affiliation: Canada France Hawaii Telescope 
# Location: Hawaii USA
# Date: Aug/2011
# Contact: opera@cfht.hawaii.edu
# 
# Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see:
# http://software.cfht.hawaii.edu/licenses
# -or-
# http://www.gnu.org/licenses/gpl-3.0.html
#
#///
#/// \package Makefile
#/// \brief This is the main harness entry point - includes all other packages.
#/// \brief {It is the root of the harness, which is called by 
#///  the "opera" tool. All other specialized Makefiles are 
#///  included here.} 
#/// \brief {This version implements OLAPA as well as EEV1 modes.
#///  As of Sept 27 2010 there is some level of DANCE support.
#/// 	The new header keywords for DANCE:
#/// 	DETECTOR = OLAPA
#/// 	AMPLIST will be "a", "b", or "a,b".  
#/// 	If AMPLIST is "a,b", then there will be GAINA/GAINB and RDNOISEA/RDNOISEB
#/// 	instead of GAIN and READNOISE.}
#/// \defgroup harness Harness
#
#########################################################################################

#########################################################################################
#///
#/// \cond makefile
#///
#########################################################################################

.PHONY: 
	
.SUFFIXES: .fits .s .l .e .geom .wcal .fcal .aper .dist .geom .sn .cpu .cal .eps

.PRECIOUS:  
.SECONDARY:

# for Debian... 
export SHELL := $(shell which bash)

#########################################################################################
# IMPORTANT to define the instrument here, selects the set of Makefiles to use...
#########################################################################################

instrument	:= espadons

#########################################################################################
# basic configuration and parameter definitions
#########################################################################################

include $(opera)/harness/$(instrument)/Makefile.parameters
include $(opera)/harness/$(instrument)/Makefile.configuration

#########################################################################################
# CFHT-specific Makefiles for espadons
#########################################################################################

ifeq ($(observatory),CFHT)

ifeq ($(TRACE),1)
	pref	:= \#
	wpref	:= \#
	epref	:= \#
else
ifdef director
	pref	:= status: [`whoami`@$(MACHINE)]
	wpref	:= warning: [`whoami`@$(MACHINE)]
	epref	:= error: [`whoami`@$(MACHINE)]
else
	pref	:= $(shell date +'%D %H:%M:%S') [`whoami`@$(MACHINE)] 
	wpref	:= $(shell date +'%D %H:%M:%S') [`whoami`@$(MACHINE)] 
	epref	:= $(shell date +'%D %H:%M:%S') [`whoami`@$(MACHINE)] 
endif
endif

define deltat
	`echo "$$SECONDS - $$start" | bc | awk '{if ($$1 > 3600) printf("%d:%02d:%02d",$$1/3600,($$1/60)%60,$$1%60); else if ($$1 > 60) printf("%2d Min. %2d Sec.",($$1/60)%60,$$1%60); else printf("%2d Sec.",$$1%60);}'`
endef

all:
	@echo "$(pref) OPERA pipeline $(versionstr) configured for $(observatory)"
	
include $(opera)/harness/$(instrument)/Makefile.calibration
include $(opera)/harness/$(instrument)/Makefile.core

endif

#########################################################################################
# harness utilities
#########################################################################################

include $(opera)/harness/$(instrument)/Makefile.util

#########################################################################################
# CFHT and Libre-Esprit compatibility targets 
#########################################################################################

ifeq ($(observatory),CFHT)
	include $(opera)/harness/$(instrument)/Makefile.upena
endif

#########################################################################################
# Add Analysis targets here
# e.g. include Makefile.lsd.espadons
#########################################################################################

include $(opera)/harness/$(instrument)/Makefile.analysis

VPATH		=  $(reductiondir):$(specdir):$(byproductsdir):$(calibrationdir):$(tmpdir):$(visualsdir):$(analysesdir):

#########################################################################################
# general utilities
#########################################################################################

# eval ARG to instantiate values for parameter and configuration access
eval:
	@echo $($(ARG))

#########################################################################################
# last ditch - unknown command
#########################################################################################

%:
	@echo "$(epref) Invalid command: \"$*\"."
	@echo "$(epref) Issue the command \"opera --help\" for the list of commands."
	
#########################################################################################
#///
#/// \endcond
#///
#########################################################################################

