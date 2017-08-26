#! /bin/bash
#########################################################################################
#
# Script name: a remote runner for opera
# Version: 1.0
# Description: a remote runner for opera
# Author(s): CFHT OPERA team
# Affiliation: Canada France Hawaii Telescope 
# Location: Hawaii USA
# Date: May/2013
# Contact: opera@cfht.hawaii.edu
# 
# Copyright (C) 2013  Opera Pipeline team, Canada France Hawaii Telescope
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
#/// \package opera.sh
#/// \brief a remote runner for opera
#/// \ingroup scripts
#
#########################################################################################
# 
#
function usage() {
	echo "$(basename $0): [--machine=...] <opera commands>"
	echo "		remote execution of opera commands."
   return 0
}

for arg in $@
do
   case ${arg} in
     *help|-h) 
         usage;
         exit 0;
         ;;
     --machine=*) 
        machine=${arg#--machine=};
		shift;
		  ;;
	   *)
		  ;;
	esac
done
home=$HOME/opera-1.0
if [ ! -d $home ]
then
	if [[ "`which opera`" == "" ]]
	then
		echo "$(basename $0): Can't find opera installation in $HOME/opera-1.0/"
		exit 1
	else
		home=$(dirname `which opera`)
	fi
fi
if [[ "`hostname`" == "$machine" || "" == "$machine" ]]
then
	export opera=$home; export PATH=$PATH:$opera/bin/:/usr/local/bin/; opera $@
else
	rsh $machine "export opera=$home; export PATH=$PATH:$opera/bin/:/usr/local/bin/; opera $@"
fi
exit $?

