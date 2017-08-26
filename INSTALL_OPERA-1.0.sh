#!/bin/bash
#
# Perform a default OPERA instalation in ~
# Usage : . ./$(basename $0)
#
cd ~
operaversion=opera-1.0
zip=opera-1.0.zip
delete=0
if [[ "$1" == "-h" || "$1" == "--help" ]]
then
	echo " Usage : cd ~; . ./$(basename $0) [--delete] [--png]"
	echo "	Installs opera in home directory. The --png option builds with png plotting support."
else
	if [[ "$1" == "--delete" ||  "$1" == "-d" ]]
	then
		delete=1
		shift
	fi
	if [[ "$1" != "" ]]
	then
		zip=$1
	else
		zip="`dirname $0`/$zip"
	fi
	if [[ -d $operaversion && delete == 0 ]]
	then
			echo "Please rm -rf $operaversion/"
	else
		if [[ ! -e $zip ]]
		then
			echo "$zip not found."
		else
			if [[ -d $operaversion && delete == 1 ]]
			then
				rm -rf $operaversion
			fi
			rm -rf __MACOSX/
			if [[ $zip =~ [a-z\.[0-9]*\.zip ]]
			then
				unzip $zip
			else
				if [[ $zip =~ [a-z\.[0-9]*\.tar\.gz ]]
				then
					tar xvzf $zip
				else
					echo "$zip must end in .zip or .tar.gz"
					exit 1
				fi
			fi
			rm -rf __MACOSX/
			export opera=$HOME/$operaversion
			export PATH=$PATH:$opera/bin:
			cd $opera
			#
			# llvm-gcc needs special cflags MacOSX Lion
			#
			if [[ "`ls -l $(which gcc) | grep llvm`" != "" ]]
			then
				if [[ "$2" == "--png" ]]
				then
					echo "./configure --prefix=$opera --enable-LIBPNG=true  CFLAGS=\"-arch x86_64\" CXXFLAGS=\"-arch x86_64 -Wno-long-long\""
					./configure --prefix=$opera --enable-LIBPNG=true CFLAGS="-arch x86_64" CXXFLAGS="-arch x86_64 -Wno-long-long"
				else
					echo "./configure --prefix=$opera CFLAGS=\"-arch x86_64\" CXXFLAGS=\"-arch x86_64 -Wno-long-long\""
					./configure --prefix=$opera CFLAGS="-arch x86_64" CXXFLAGS="-arch x86_64 -Wno-long-long"
				fi
			else
				if [[ "$2" == "--png" ]]
				then
					echo "./configure --prefix=$opera --enable-LIBPNG=true"
					./configure --prefix=$opera --enable-LIBPNG=true
				else
					echo "./configure --prefix=$opera"
					./configure --prefix=$opera
				fi
			fi
			make install --jobs=8
			if (( $? != 0 ))
			then
				echo "###############################################################################"
				echo "# The opera build failed."
				echo "###############################################################################"
			else
				#make mostlyclean
				echo "###############################################################################"
				echo "# opera is installed in $operaversion on `date`"
				echo "# `cat <TIMESTAMP`"
				echo "# You may want to put these lines in your .bashrc for convenience:"
				echo "# export opera=\$HOME/$operaversion"
				echo "# export PATH=\$PATH:$opera/bin:"
				echo "# "
				echo "# Enjoy! -- The OPERA Team"
				echo "###############################################################################"
			fi
		fi
	fi
fi
