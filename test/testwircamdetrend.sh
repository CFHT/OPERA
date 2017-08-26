#! /bin/bash
###############################################################################
#
# test wircam detrend
#
###############################################################################
#name_raw=/data/iiwi/byproducts/12Aw04/12AD89//1557107o.fits
#name_flat=/data/iiwi/calibrations/mastertwilightflat_Ks_12Aw04_v100.fits
#name_dark=/data/iiwi/calibrations/masterdark_025s_12Aw04_v100.fits
#name_badpix=/data/iiwi/calibrations/badpix16_20120727HST192309_v100.fits
#name_weight=/data/ula/wircam/iiwi2/reductions/12Aw04/12AD89//1557107w.fits
#name_detrended=/data/ula/wircam/iiwi2/byproducts/12Aw04/12AD89//1557107s-uncalibrated.fits
#param_iiwiversion=2.1.100
#param_procdate=2012-08-17HST15:40:01
#param_linearize=       1
#param_subrefpix=       1
#param_needweight=       1
#param_gwinxtalk=MASK
#param_maskingthreshold=      5.00000
#grade=       1
###############################################################################
bindir=$HOME/opera-1.0/bin/
path=/data/WIRCamtest/

${bindir}/wirDetrend \
	--name_raw=${path}1557107o.fits \
	--name_flat=${path}mastertwilightflat_Ks_12Aw04_v100.fits \
	--name_dark=${path}masterdark_025s_12Aw04_v100.fits \
	--name_weight=${path}/1557107w.fits \
	--name_detrended=${path}/1557107s-uncalibrated.fits \
	--param_linearize=1 \
	--param_subrefpix=1 \
	--param_needweight=1 \
	--param_gwinxtalk=MASK \
	--param_maskingthreshold=5.0
	
exit


