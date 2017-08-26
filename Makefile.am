#########################################################################################
# This configuration file lists all the subirectories to be traversed 
# to create Makefiles...
#########################################################################################

AUTOMAKE_OPTIONS = foreign

AM_MAKEFLAGS	= PREFIX=$(prefix)
ACLOCAL_AMFLAGS = -I m4

# add doc to this if you want to build the doxygen documentation
SUBDIRS = \
	src/libraries src/tools \
	src/core-espadons src/core-guest \
	src/analysis-espadons src/analysis-guest \
	scripts scripts/espadons scripts/guest
