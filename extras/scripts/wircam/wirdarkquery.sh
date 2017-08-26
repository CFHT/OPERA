#! /bin/bash
###############################################################################
#
# darks query - figure which nights have twilight flats and darks
# vs no twilight flats and darks
#
###############################################################################
bindir=$HOME/opera-1.0/bin/

if (( $# == 0 ))
then
	echo "usage: $(basename $0): <crunid> "
	exit 0
fi
crunid=$1

############################## GET THE DATES ##############################
datesWithAMFlats="`${bindir}wirenvdb \"select distinct convert(varchar,dateadd(hour,-10,mydatetime),105) from exposure where crunid='${crunid}' and object='TwilightFlats' and convert(varchar,dateadd(hour,-10,mydatetime),108) between '04:00:00' and '09:00:00' and _use>=0\"`"
datesWithAMDarks="`${bindir}wirenvdb \"select distinct convert(varchar,dateadd(hour,-10,mydatetime),105) from exposure where crunid='${crunid}' and object='Darks' and convert(varchar,dateadd(hour,-10,mydatetime),108) between '04:00:00' and '11:00:00' and _use>=0\"`"

# WITH FLATS dates for both darks and AM flats
both=""
for amdarkdate in $datesWithAMDarks 
do
	for amflatdate in $datesWithAMFlats 
	do
		if [[ "${amflatdate}" == "${amdarkdate}" ]]
		then
			#echo "${amflatdate} has both darks and AM flats."
			both="${both} ${amflatdate}"
		fi
	done
done
# WITHOUT FLATS dates for darks but no AM flats
onlydarks=""
for date in ${datesWithAMDarks}
do
	if [[ "`operafindword ${date} \"${datesWithAMFlats}\"`" == "" ]]
	then
		#echo "${date} has only darks."
		onlydarks="${onlydarks} ${date}"
	fi
done

############################## GET THE ODOMETERS ##############################
# WITH FLATS - get the odometers of the darks on dates that also have flats -- something wrong,  only getting odometers for April 30....
for night in ${both}
do
	${bindir}wirenvdb "select distinct odometer from exposure where crunid='${crunid}' and object='Darks' and convert(varchar,dateadd(hour,-10,mydatetime),105) = '${night}'"
done

#separate the output with a single blank space
echo

# WITHOUT FLATS - get the odometers for dates with only darks
for night in ${onlydarks}
do
	${bindir}wirenvdb "select distinct odometer from exposure where crunid='${crunid}' and object='Darks' and convert(varchar,dateadd(hour,-10,mydatetime),105) = '${night}'"
done

exit