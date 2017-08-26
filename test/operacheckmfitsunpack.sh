#!/bin/bash
#
# check the unpacking of the fits products
#
# Example checks:
#
# operaCreateProduct --beam=/data/uhane5/opera//spectra/11AQ12-Jun16/1310060i.e.gz --centralsnr --polar= --gain=/data/uhane5/opera//calibrations/11AQ12-Jun16/OLAPAab_pol_Normal.gain.gz --geom=/data/uhane5/opera//calibrations/11AQ12-Jun16/OLAPAab_pol_Normal.geom.gz --ordp=/data/uhane5/opera//calibrations/11AQ12-Jun16/OLAPAab_pol_Normal.ordp.gz --disp=/data/uhane5/opera//calibrations/11AQ12-Jun16/OLAPAab_pol_Normal.disp.gz --rvel=/data/uhane5/opera//calibrations/11AQ12-Jun16/1310060i.rvel.gz --prvel= --tell=/data/uhane5/opera//calibrations/11AQ12-Jun16/1310060i.tell.gz --ptell= --wave=/data/uhane5/opera//calibrations/11AQ12-Jun16/OLAPAab_pol_Normal.wcal.gz --fcal=/data/uhane5/opera//calibrations/11AQ12-Jun16/masterfluxcalibration_OLAPAab_pol_Normal3.fcal.gz --sequence=3 --snr=/data/uhane5/opera//spectra/11AQ12-Jun16/1310060ie.sn.gz --aper=/data/uhane5/opera//calibrations/11AQ12-Jun16/OLAPAab_pol_Normal.aper.gz --prof=/data/uhane5/opera//calibrations/11AQ12-Jun16/OLAPAab_pol_Normal.prof.gz --parameters=/tmp/11AQ12-Jun16/1310060m.parm --version="opera-1.0.652 build date Tue Mar 19 16:02:48 HST 2013" --date="Wed Mar 20 09:11:16 HST 2013" --spectrumtype=CalibratedOptimalBeamSpectrum --compressiontype=0 --input=/data/uhane5/opera/11AQ12-Jun16/1310060o.fits --output=/data/uhane5/opera//spectra/11AQ12-Jun16/1310060m.fits -v
# sh ./test/operacheckmfitsunpack.sh 11AQ12-Jun16 1310060m.fits --clean
# sh ./test/operacheckmfitsunpack.sh 11AQ12-Jun16 1310060m.fits
#
if [[ $# < 2 || "$1" == "--help" ]]
then
	echo "usage: $(basename $0) <NIGHT> [<odometer>m.fits | <odometer>p.fits.fz | <odometer>i.fits.fz] [--clean]"
	exit
fi
rootdir=/data/uhane5/opera/
spectradir=$rootdir/spectra/$1
calsdir=$rootdir/calibrations/$1
fits=$2
clean=0
if [[ "$3" == "--clean" ]]
then
	clean=1
fi
case $fits
in
	*p.fits.fz)
		fitstype="ip";
	;;
	*i.fits.fz)
		fitstype="ip";
	;;
	*m.fits)
		fitstype="m";
	;;
	*)
		echo "usage: $(basename $0) <NIGHT> [<odometer>m.fits | <odometer>p.fits.fz | <odometer>i.fits.fz] [--clean]"
		exit
	;;
esac
products=`ls *.gz *.s 2>/dev/null`
if (( clean ))
then
	for product in $products
	do
		unzipped=$(basename $product .gz)
		case $fitstype
		in
			m)
				case $unzipped
				in
				*.e|*.sn|*.tell|*.rvel)
					rm -f $unzipped
					rm -f $product
					rm -f $unzipped.difflog
					cd $spectradir >/dev/null
					rm -f $unzipped
					cd - >/dev/null
				;;
				*)
					rm -f $unzipped
					rm -f $product
					rm -f $unzipped.difflog
					cd $calsdir >/dev/null
					rm -f $unzipped
					cd - >/dev/null
				;;
				esac
			;;
			ip)
				rm -f $unzipped
				rm -f $unzipped.difflog
				cd $spectradir >/dev/null
				rm -f $unzipped
				cd - >/dev/null
			;;
		esac
	done
	exit
fi
echo "-------------------------------------------"
operaExtractProducts -v $spectradir/$fits
echo "-------------------------------------------"
products=`ls *.gz *.s 2>/dev/null`
for product in $products
do
	unzipped=$(basename $product .gz)
	case $fitstype
	in
		m)
			case $unzipped
			in
			*.e|*.sn|*.tell|*.rvel)
				cp -p $product $product.bak
				gunzip -f -d $product
				mv $product.bak $product
				cd $spectradir >/dev/null
				cp -p $product $product.bak
				gunzip -f -d $product
				mv $product.bak $product
				cd - >/dev/null
				if [[ "`diff -w -b -q $unzipped $spectradir/`" != "" ]]
				then
					echo "$unzipped differs"
					diff -w -b -s $unzipped $spectradir/ >$unzipped.difflog
				else
					echo "$unzipped OK"
				fi
			;;
			*)
				cp -p $product $product.bak
				gunzip -f -d $product
				mv $product.bak $product
				cd $calsdir >/dev/null
				cp -p $product $product.bak
				gunzip -f -d $product
				mv $product.bak $product
				cd - >/dev/null
				if [[ "`diff -w -b -q $unzipped $calsdir/`" != "" ]]
				then
					echo "$unzipped differs"
					diff -w -b $unzipped $calsdir/ >$unzipped.difflog
				else
					echo "$unzipped OK"
				fi
			;;
			esac
		;;
		ip)
			cd $spectradir >/dev/null
			cp -p ${product}.gz ${product}.gz.bak
			gunzip -f -d ${product}.gz
			mv ${product}.gz.bak ${product}.gz
			cd - >/dev/null
			if [[ "`diff -w -b -q $unzipped $spectradir/`" != "" ]]
			then
				echo "$unzipped differs"
				diff -w -b $unzipped $spectradir/ >$unzipped.difflog
			else
				echo "$unzipped OK"
			fi
		;;
	esac
done
exit
