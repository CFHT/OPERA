#!/bin/bash
#
# run some tests -- this is a CFHT-only script
# Other installations should write a test script for their site
#
trap "echo '$(basename $0) aborted...'; exit 1" SIGINT SIGTERM

OSTYPE=`uname -s | awk '{print tolower($0)}'`

function usage() {
	echo "This is a CFHT-only test script."
	echo "$(basename $0) commands:"
	echo "$(basename $0) NIGHT=<...> where <...> is like 11AQ02-Jan03 [-p]"
	echo "$(basename $0) <...> <...> ... where <...> is like 11AQ02-Jan03 [-p]"
	echo "$(basename $0) use the -p option to enable plot testing (requires libpng)"
}

optargs=""
NIGHT=11AQ12-Jun16
NIGHTS=""
((errors=0))
((plotting=0))
for arg in $@
do
   case "${arg}" in
     *h*) 
         usage;
         exit 0;
         ;;
     [0-9][0-9][AB][QE][0-9][0-9]-[A-Z][a-z][a-z][0-9][0-9]) 
         nights="$nights $arg";
			;;
     NIGHT=[0-9][0-9][AB][QE][0-9][0-9]-[A-Z][a-z][a-z][0-9][0-9]) 
         NIGHTS="$NIGHTS ${arg#NIGHT=}";
         NIGHT="${arg#NIGHT=}";
			;;
      -p) 
         ((plotting=1));
         ;;
      -v) 
         optargs += "${arg} ";
         ;;
      -t) 
         optargs += "${arg} ";
         ;;
      -d) 
         optargs += "${arg} ";
         ;;
      *) 
         usage;
         exit 0;
         ;;
	esac
done

if [[ $OSTYPE == darwin ]]
then
	datapath=/data/espadons/$NIGHT/
	spectrapath=/data/espadons/spectra/
	byproductspath=$HOME/opera/byproducts/
	calibrationspath=$HOME/opera/calibrations/
	configdir=$HOME/opera-1.0/config/
	visualspath=$HOME/opera/visuals/
fi

if [[ $OSTYPE == linux ]]
then
	datapath=/data/uhane5/opera/$NIGHT/
	byproductspath=/data/uhane5/opera/byproducts/
	calibrationspath=/data/uhane5/opera/calibrations/
	spectrapath=/data/uhane5/upena/spectra/
	configdir=$HOME/opera-1.0/config/
	visualspath=/data/uhane5/opera/visuals/
fi
if [[ $OSTYPE == darwin ]]
then
	otool -Vvt $HOME/opera-1.0/test/operaAsmTest.o
fi
if [[ $OSTYPE == linux ]]
then
	objdump -dS $HOME/opera-1.0/test/operaAsmTest.o
fi
if [[ -d $datapath ]]
then 
	#################################################################
	# Other image formats:
	#
	# ./bin/operaNIFSImageTest --inputImage=/data/NIFS-TESTDATA/N20120808S0704.fits
	# ./bin/operaAOBImageTest --inputImage=/data/AOB-TESTDATA/1426179o.fits
	# ./bin/operaNICIImageTest --inputImage=/data/NICI-TESTDATA/S20120913S0130.fits.gz
	#
	#################################################################

	#################################################################
	# Analysis:
	#
	# ./bin/opera DATADIRS="/data/niele/espadons/12BQ10-Nov27 /data/niele/espadons/12BQ10-Nov28 /data/niele/espadons/12BQ10-Nov29 /data/niele/espadons/12BQ10-Nov30 /data/niele/espadons/12BQ10-Dec01" OUT=rvel51Peg.seq OBJECT="51 Peg" radialvelocitysequence
	# ./bin/opera DATADIR=/data/uhane5/opera/12BQ11-Dec23 polarbinning
	#
	#################################################################

	#################################################################
	# Balmer Series:
	#
	# ./bin/balmer-series.sh 1515007.e.gz calibrations/GalileanMoons.OLAPAa_pol_Normal,wcar,gz 1515007
	#
	#################################################################

	#################################################################
	# Master Flux Calibrations:
	#
	# ./bin/opera DATADIR=/data/niele/espadons/12BQ10-Nov27 calibrations
	# ./bin/opera DATADIR=/data/niele/espadons/12BQ10-Nov27 intensity
	# ./bin/opera DATADIR=/data/niele/espadons/12BQ10-Nov27 masterfluxcalibrations
	#
	#################################################################

	#################################################################
	# Master Darks Analysis:
	#
	# opera instrument=wircam masterdarkanalysis DATADIR=/data/wircam/12AQ03-Mar01/ -t -p -v
	#
	#################################################################

	#################################################################
	echo "#################################################################"
	echo "operaPolarTest --stokesparameter=1 --output=testpolar_pu.s ${optargs}"
	operaPolarTest --stokesparameter=1 --output=testpolar_pu.s ${optargs}
	if (( $? != 0 ))
	then
		((errors++))
	fi
	#################################################################
	echo "#################################################################"
	echo "operaPolarimetryTest ${optargs}"
	operaPolarimetryTest ${optargs}
	if (( $? != 0 ))
	then
		((errors++))
	fi
	#################################################################
	echo "#################################################################"
	echo "operaFluxVectorTest ${optargs}"
	operaFluxVectorTest ${optargs}
	if (( $? != 0 ))
	then
		((errors++))
	fi
	#################################################################
	echo "#################################################################"
	echo "wircolorcomposite wircolorcomposite --R=/data/wircam/1317665p.fits --G=/data/wircam/1317671p.fits --B=/data/wircam/1317677p.fits --badpixelmask=/data/wircam/badpix_20061003.fits --output=color.fits ${optargs}"
	wircolorcomposite --R=/data/wircam/1317665p.fits --G=/data/wircam/1317671p.fits --B=/data/wircam/1317677p.fits --badpixelmask=/data/wircam/badpix_20061003.fits --output=color.fits ${optargs}
	if (( $? != 0 ))
	then
		((errors++))
	fi
	#################################################################
	echo "#################################################################"
	echo "operaExtractionApertureTest --input=$datapath/1310085o.fits --output=$visualspath/foo.fits ${optargs}"
	operaExtractionApertureTest --input=$datapath/1310085o.fits --output=$visualspath/foo.fits ${optargs}
	fitsverify $visualspath/foo.fits
	if (( $? == 0 ))
	then
		rm -f $visualspath/foo.fits
	else
		((errors++))
	fi
	#################################################################
	echo "#################################################################"
	if [[ -e $calibrationspath/$NIGHT/OLAPAa_sp2_Normal.prof && -e $calibrationspath/$NIGHT/OLAPAa_sp2_Normal.geom ]]
	then
		echo "operaPlotInstrumentProfile --inputprof=$calibrationspath/$NIGHT/OLAPAa_sp2_Normal.prof --inputgeom=$calibrationspath/$NIGHT/OLAPAa_sp2_Normal.geom -R 1000 --minorder=22 --maxorder=61  -P $visualspath/testipplot.eps -F $visualspath/testipplot.dat -S $visualspath/testipplot.gnu ${optargs}"
		operaPlotInstrumentProfile --inputprof=$calibrationspath/$NIGHT/OLAPAa_sp2_Normal.prof --inputgeom=$calibrationspath/$NIGHT/OLAPAa_sp2_Normal.geom -R 1000 --minorder=22 --maxorder=61  -P $visualspath/testipplot.eps -F $visualspath/testipplot.dat -S $visualspath/testipplot.gnu ${optargs}
		if (( $? == 0 ))
		then
			rm -f $visualspath/testipplot.dat
			rm -f $visualspath/testipplot.gnu
			echo "view $visualspath/testipplot.eps to see result"
		else
			((errors++))
		fi
	else
		echo "$calibrationspath/$NIGHT/OLAPAa_sp2_Normal.prof and $calibrationspath/$NIGHT/OLAPAa_sp2_Normal.geom are needed for this test"
	fi

	#################################################################
	echo "#################################################################"
	echo "operaImageOperatorTest --bias=$datapath/1310052b.fits --flat1=$datapath/1310041f.fits --flat2=$datapath/1310050f.fits --badpix=$opera/config/badpix_olapa-a.fits.fz --output=$visualspath/foo.fits ${optargs}"
	operaImageOperatorTest --bias=$datapath/1310052b.fits --flat1=$datapath/1310041f.fits --flat2=$datapath/1310050f.fits --badpix=$opera/config/badpix_olapa-a.fits.fz --output=$visualspath/foo.fits ${optargs}
	fitsverify $visualspath/foo.fits
	if (( $? == 0 ))
	then
		rm -f $visualspath/foo.fits
	else
		((errors++))
	fi
	#################################################################
	echo "#################################################################"
	echo "operaFITSImageTest --images=$datapath/1310169o.fits --images=$datapath/1310170o.fits --images=$datapath/1310171o.fits --output=$visualspath/foo.fits ${optargs}"
	operaFITSImageTest --images=$datapath/1310169o.fits --images=$datapath/1310170o.fits --images=$datapath/1310171o.fits --output=$visualspath/foo.fits --output=$visualspath/float.fits ${optargs}
	fitsverify $visualspath/foo.fits
	if (( $? == 0 ))
	then
		rm -f $visualspath/foo.fits
	else
		((errors++))
	fi
	fitsverify $visualspath/float.fits
	if (( $? == 0 ))
	then
		rm -f $visualspath/float.fits
	else
		((errors++))
	fi
	#################################################################
	echo "#################################################################"
	echo "operaFITSSubImageTest --bias=$datapath/1310052b.fits --flat1=$datapath/1310041f.fits --flat2=$datapath/1310050f.fits --badpix=$opera/config/badpix_olapa-a.fits.fz --output=$visualspath/foo.fits ${optargs}"
	operaFITSSubImageTest --bias=$datapath/1310052b.fits --flat1=$datapath/1310041f.fits --flat2=$datapath/1310050f.fits --badpix=$opera/config/badpix_olapa-a.fits.fz --output=$visualspath/foo.fits ${optargs}
	fitsverify $visualspath/foo.fits
	if [[ `which imlist` ]]
	then
		imlist $visualspath/foo.fits[1:8,1:8]
		echo "Note that the imlist is upside down..."
	fi
	#################################################################
	echo "#################################################################"
	echo "operaEspadonsImageTest --image=$datapath/1310171o.fits --output=$visualspath/foo.fits ${optargs}"
	operaEspadonsImageTest --image=$datapath/1310171o.fits --output=$visualspath/foo.fits ${optargs}
	fitsverify $visualspath/foo.fits
	#################################################################
	echo "#################################################################"
	echo "operagetheader --keyword=INSTMODE --keyword=EREADSPD $datapath/1310052b.fits"
	operagetheader --keyword=INSTMODE --keyword=EREADSPD $datapath/1310052b.fits
	echo "operagetheader --printfilename INSTMODE $datapath/1310052b.fits $datapath/1310041f.fits"
	operagetheader --printfilename INSTMODE $datapath/1310052b.fits $datapath/1310041f.fits
	if (( $? != 0 ))
	then
		((errors++))
	fi
	#################################################################
	echo "#################################################################"
	echo "operasaturated --saturated=35000.0 --count=500 $datapath/1310041f.fits"
	operasaturated --saturated=35000.0 --count=500 $datapath/1310041f.fits
	if (( $? != 0 ))
	then
		((errors++))
	fi
	if [[ -d $spectrapath/10BQ06-Oct20/ && -d $spectrapath/11AQ12-Jun16/ ]]
	then
		#################################################################
		echo "#################################################################"
		echo "Polar:"
		echo "operaFITSProductTest --product=$spectrapath/10BQ06-Oct20/1252398p.fits"
		operaFITSProductTest --product=$spectrapath/10BQ06-Oct20/1252398p.fits
		if (( $? != 0 ))
		then
			((errors++))
		fi
		echo "#################################################################"
		echo "Polar, based on o.fits:"
		echo "operaFITSProductTest --product=/data/espadons/spectra/11AQ12-Jun16/1310058i.fits --basedon=$datapath/1310058o.fits"
		operaFITSProductTest --product=$spectrapath/11AQ12-Jun16/1310058i.fits --basedon=$datapath/1310058o.fits
		if (( $? != 0 ))
		then
			((errors++))
		fi
		
		echo "Star Only:"
		echo "operaFITSProductTest --product=$spectrapath/10BQ06-Oct20/1252386i.fits"
		operaFITSProductTest --product=$spectrapath/10BQ06-Oct20/1252386i.fits
		if (( $? != 0 ))
		then
			((errors++))
		fi
		echo "Star+Sky:"
		echo "operaFITSProductTest --product=$spectrapath/09AQ02-Feb04/1057897i.fits"
		operaFITSProductTest --product=$spectrapath/09AQ02-Feb04/1057897i.fits
		if (( $? != 0 ))
		then
			((errors++))
		fi
	else
		echo "Couldnt find $spectrapath on this machine, skipping..."
	fi
	if (( plotting ))
	then
		# create the visuals dir if it doesn't exist
		if [[ ! -d $visualspath/11AQ12-Jun16/ ]]
		then
			mkdir -p $visualspath/11AQ12-Jun16/
		fi
		if [[ ! -d $visualspath/11AQ14-Jul07/ ]]
		then
			mkdir -p $visualspath/11AQ14-Jul07/
		fi
		if [[ ! -d$ visualspath/11AQ14-Jul14/ ]]
		then
			mkdir -p $visualspath/11AQ14-Jul14/
		fi
		if [[ ! -d $visualspath/11AQ12-Jun16/ ]]
		then
			mkdir -p $visualspath/11AQ12-Jun16/
		fi
		if [[ ! -d $visualspath/10BQ06-Oct20/ ]]
		then
			mkdir -p $visualspath/10BQ06-Oct20/
		fi
		if [[ ! -d $visualspath/09AQ02-Feb04/ ]]
		then
			mkdir -p $visualspath/09AQ02-Feb04/
		fi
		if [[ -d $spectrapath/11AQ14-Jul07/ && -d $spectrapath/11AQ14-Jul14/ && -d $datapath/11AQ12-Jun16/ ]]
		then
			#################################################################
			echo "Plots"
			echo "operaPlot --product=$spectrapath/11AQ14-Jul07/1314960p.fits --output=$visualspath/11AQ14-Jul07/1314960p.png -v"
			operaPlot --product=$spectrapath/11AQ14-Jul07/1314960p.fits --output=$visualspath/11AQ14-Jul07/1314960p.png -v
			if (( $? != 0 ))
			then
				((errors++))
			fi
			echo "operaPlot --product=$spectrapath/11AQ14-Jul14/1316722i.fits --output=$visualspath/11AQ14-Jul14/1316722i.png -v"
			operaPlot --product=operaPlot --product=$spectrapath/11AQ14-Jul14/1316722i.fits --output=$visualspath/11AQ14-Jul14/1316722i.png -v
			if (( $? != 0 ))
			then
				((errors++))
			fi
			echo "operaPlot --input=$datapath/11AQ12-Jun16/1310101o.fits --output=$visualspath/11AQ12-Jun16/1310101o.png -v"
			operaPlot --input=$datapath/11AQ12-Jun16/1310101o.fits --output=$visualspath/11AQ12-Jun16/1310101o.png -v
			if (( $? != 0 ))
			then
				((errors++))
			fi
			echo "operaPlot --input=$datapath/11AQ12-Jun16/1310101o.fits --output=$visualspath/11AQ12-Jun16/ortho-1310101o.png -v -m"
			operaPlot --input=$datapath/11AQ12-Jun16/1310101o.fits --output=$visualspath/11AQ12-Jun16/ortho-1310101o.png -v -m
			if (( $? != 0 ))
			then
				((errors++))
			fi
			echo "Polar (with datapoints):"
			echo "operaPlotGeom --geom=$calibrationspath/11AQ12-Jun16/OLAPAab_pola_Normal.geom --uncalibrated=$byproductspath/11AQ12-Jun16/OLAPAab_pola_Normal.s --masterflat=$calibrationspath/11AQ12-Jun16/masterflat_OLAPAab_pola_Normal.fits --output=$visualspath/11AQ12-Jun16/all-OLAPAab_pol_a_Normal.png -v"
			operaPlotGeom --geom=$calibrationspath/11AQ12-Jun16/OLAPAab_pola_Normal.geom --uncalibrated=$byproductspath/11AQ12-Jun16/OLAPAab_pola_Normal.s --masterflat=$calibrationspath/11AQ12-Jun16/masterflat_OLAPAab_pola_Normal.fits --output=$visualspath/11AQ12-Jun16/all-OLAPAab_pol_a_Normal.png -v
			if (( $? != 0 ))
			then
				((errors++))
			fi
			echo "#################################################################"
			echo "Polar:"
			echo "operaPlotGeom --geom=$calibrationspath/11AQ12-Jun16/OLAPAab_pola_Normal.geom --masterflat=$calibrationspath/11AQ12-Jun16/masterflat_OLAPAab_pola_Normal.fits --output=$visualspath/11AQ12-Jun16/OLAPAab_pol_a_Normal.png -v"
			operaPlotGeom --geom=$calibrationspath/11AQ12-Jun16/OLAPAab_pola_Normal.geom --masterflat=$calibrationspath/11AQ12-Jun16/masterflat_OLAPAab_pola_Normal.fits --output=$visualspath/11AQ12-Jun16/OLAPAab_pol_a_Normal.png -v
			if (( $? != 0 ))
			then
				((errors++))
			fi
			echo "Star Only:"
			echo "operaPlotGeom --geom=$calibrationspath/11AQ12-Jun16/OLAPAab_pola_Normal.geom --masterflat=$calibrationspath/10BQ06-Oct20/masterflat_OLAPAab_pola_Normal.fits --output=$visualspath/10BQ06-Oct20/OLAPAab_pol_a_Normal.png -v"
			operaPlotGeom --geom=$calibrationspath/10BQ06-Oct20/OLAPAab_pola_Normal.geom --masterflat=$calibrationspath/10BQ06-Oct20/masterflat_OLAPAab_pola_Normal.fits --output=$visualspath/10BQ06-Oct20/OLAPAab_pol_a_Normal.png -v
			if (( $? != 0 ))
			then
				((errors++))
			fi
			echo "Star+Sky:"
			echo "operaPlotGeom --geom=$calibrationspath/09AQ02-Feb04/OLAPAab_pola_Normal.geom --masterflat=$calibrationspath/09AQ02-Feb04/masterflat_OLAPAab_pola_Normal.fits --output=$visualspath/09AQ02-Feb04/OLAPAab_pol_a_Normal.png -v"
			operaPlotGeom --geom=$calibrationspath/09AQ02-Feb04/OLAPAab_pola_Normal.geom --masterflat=$calibrationspath/09AQ02-Feb04/masterflat_OLAPAab_pola_Normal.fits --output=$visualspath/09AQ02-Feb04/OLAPAab_pol_a_Normal.png -v
			if (( $? != 0 ))
			then
				((errors++))
			fi
		else
			echo "Couldnt find $spectrapath on this machine, skipping..."
		fi
	fi
	#################################################################
	echo "#################################################################"
	echo "opera NIGHT=$NIGHT clean ${optargs}"
	opera NIGHT=$NIGHT clean ${optargs}
else
	echo "Couldnt find $datapath on this machine, skipping..."
fi

#################################################################
# do library tests
#################################################################
echo "#################################################################"
echo "operaStatsLibTest"
operaStatsLibTest  ${optargs}
if (( $? != 0 ))
then
	((errors++))
fi
echo "#################################################################"
echo "operaMatrixLibTest"
operaMatrixLibTest  ${optargs}
if (( $? != 0 ))
then
	((errors++))
fi

echo "#################################################################"
echo "operaMathLibTest"
operaMathLibTest  ${optargs}
if (( $? != 0 ))
then
	((errors++))
fi

#################################################################
echo "#################################################################"
echo "operaFitLibTest"
operaFitLibTest  ${optargs}
if (( $? != 0 ))
then
	((errors++))
fi

#################################################################
echo "#################################################################"
echo "operaMPFitLibTest"
operaMPFitLibTest  ${optargs}
if (( $? != 0 ))
then
	((errors++))
fi

#################################################################
echo "#################################################################"
echo "testmpfit"
echo "should look like this:"
echo "*** testlinfit status = 1"
echo "  CHI-SQUARE = 2.756285    (8 DOF)"
echo "        NPAR = 2"
echo "       NFREE = 2"
echo "     NPEGGED = 0"
echo "     NITER = 3"
echo "      NFEV = 8"
echo ""
echo "  P[0] = 3.209966 +/- 0.022210     (ACTUAL 3.200000)"
echo "  P[1] = -1.770954 +/- 0.018938     (ACTUAL 1.780000)"
echo "*** testquadfit status = 1"
echo "  CHI-SQUARE = 5.679323    (7 DOF)"
echo "        NPAR = 3"
echo "       NFREE = 3"
echo "     NPEGGED = 0"
echo "     NITER = 3"
echo "      NFEV = 10"
echo ""
echo "  P[0] = 4.703829 +/- 0.097512     (ACTUAL 4.700000)"
echo "  P[1] = 0.062586 +/- 0.054802     (ACTUAL 0.000000)"
echo "  P[2] = 6.163087 +/- 0.054433     (ACTUAL 6.200000)"
echo "*** testquadfix status = 1"
echo "  CHI-SQUARE = 6.983588    (8 DOF)"
echo "        NPAR = 3"
echo "       NFREE = 2"
echo "     NPEGGED = 0"
echo "     NITER = 3"
echo "      NFEV = 8"
echo ""
echo "  P[0] = 4.696254 +/- 0.097286     (ACTUAL 4.700000)"
echo "  P[1] = 0.000000 +/- 0.000000     (ACTUAL 0.000000)"
echo "  P[2] = 6.172954 +/- 0.053743     (ACTUAL 6.200000)"
echo "*** testgaussfit status = 1"
echo "  CHI-SQUARE = 10.350032    (6 DOF)"
echo "        NPAR = 4"
echo "       NFREE = 4"
echo "     NPEGGED = 0"
echo "     NITER = 28"
echo "      NFEV = 139"
echo ""
echo "  P[0] = 0.480443 +/- 0.232235     (ACTUAL 0.000000)"
echo "  P[1] = 4.550752 +/- 0.395434     (ACTUAL 4.700000)"
echo "  P[2] = -0.062562 +/- 0.074715     (ACTUAL 0.000000)"
echo "  P[3] = 0.397472 +/- 0.089996     (ACTUAL 0.500000)"
echo "*** testgaussfix status = 1"
echo "  CHI-SQUARE = 15.516134    (8 DOF)"
echo "        NPAR = 4"
echo "       NFREE = 2"
echo "     NPEGGED = 0"
echo "     NITER = 12"
echo "      NFEV = 35"
echo ""
echo "  P[0] = 0.000000 +/- 0.000000     (ACTUAL 0.000000)"
echo "  P[1] = 5.059244 +/- 0.329307     (ACTUAL 4.700000)"
echo "  P[2] = 0.000000 +/- 0.000000     (ACTUAL 0.000000)"
echo "  P[3] = 0.479746 +/- 0.053804     (ACTUAL 0.500000)"

testmpfit  ${optargs}
if (( $? != 0 ))
then
	((errors++))
fi

#################################################################
# do some cals
#################################################################

if [[ $OSTYPE == linux ]]
then
	#################################################################
	echo "#################################################################"
	echo "opera ARCHIVE=2010/12/15-2010/12/16 NIGHT=00AQ00-Jan00 reductionset ${optargs}"
	opera ARCHIVE="2010/12/15-2010/12/16"  NIGHT=00AQ00-Jan00 reductionset ${optargs}
	#################################################################
	echo "#################################################################"
	echo "opera NIGHT=00AQ00-Jan00 clean ${optargs}"
	opera NIGHT=00AQ00-Jan00 clean ${optargs}
fi

#################################################################
if [[ $OSTYPE == darwin ]]
then
	NIGHTS="10BQ09-Nov26 10BQ11-Dec15 11AQ12-Jun16"
fi

if [[ $OSTYPE == linux ]]
then
	NIGHTS="09AQ02-Feb04 10BQ06-Oct20 10BQ11-Dec16 11AQ14-Jul04 11AQ14-Jul08 09AQ08-May01 10BQ09-Nov26 11AQ12-Jun16 10AQ02-Mar06 10BQ09-Nov27 11AQ14-Jul01 10AQ15-Jul27 10BQ11-Dec13 11AQ14-Jul02 11AQ14-Jul07 10BQ01-Aug05 10BQ11-Dec14 11AQ14-Jul03 11AQ14-Jul08"
fi

for NIGHT in $NIGHTS
do
	if [[ $OSTYPE == darwin ]]
	then
		datapath=/data/espadons/$NIGHT/
	fi
	if [[ $OSTYPE == linux ]]
	then
		datapath=/data/uhane5/opera/$NIGHT/
	fi
	if [[ -d $datapath ]]
	then
		opera datapath=$datapath clean ${optargs}
		#################################################################
		echo "#################################################################"
		echo "Here are the opera reductionlists"
		echo "opera datapath=$datapath reductionlist ${optargs}"
		opera datapath=$datapath reductionlist ${optargs}

		if [[ $OSTYPE == linux ]]
		then
			#################################################################
			echo "#################################################################"
			echo "Compare with the upena reductionlists"
			echo "upena NIGHT=$NIGHT reductionlist ${optargs}"
			opera datapath=$datapath reductionlist ${optargs}
		fi

		#################################################################
		echo "#################################################################"
		echo "opera datapath=$datapath masterflats ${optargs}"
		opera datapath=$datapath masterflats ${optargs}

		#################################################################
		echo "#################################################################"
		echo "opera datapath=$datapath masterbiases --pick=0 ${optargs}"
		opera datapath=$datapath masterbiases --pick=0 ${optargs}

		#################################################################
		echo "#################################################################"
		echo "opera datapath=$datapath mastercomparisons --pick=1 ${optargs}"
		opera datapath=$datapath mastercomparisons --pick=1 ${optargs}

		#################################################################
		echo "#################################################################"
		echo "opera datapath=$datapath masterfabperots --pick=2 ${optargs}"
		opera datapath=$datapath masterfabperots --pick=2 ${optargs}

		#################################################################
		echo "#################################################################"
		echo "opera datapath=$datapath gains ${optargs}"
		opera datapath=$datapath gains ${optargs}

		#################################################################
		echo "#################################################################"
		echo "opera datapath=$datapath geometries ${optargs}"
		opera datapath=$datapath geometries ${optargs}

		#################################################################
		echo "#################################################################"
		echo "opera datapath=$datapath wcals ${optargs}"
		opera datapath=$datapath wcals ${optargs}

		#################################################################
		echo "#################################################################"
		echo "opera NIGHT=$NIGHT clean ${optargs}"
		opera NIGHT=$NIGHT clean ${optargs}

	else
		echo "Couldnt find $datapath on this machine, skipping..."
	fi
done
#################################################################
echo "#################################################################"
echo "Tests done -- $errors error(s)"
echo "#################################################################"
