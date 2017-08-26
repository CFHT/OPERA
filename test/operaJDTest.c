
/* 
 * Test Julian Date library
 * hjd  -  programme to calculate the JD, MJD, HJD dates
 *
 * (c) 2002 Richard Ogley.
 *
 * Version 0.1   06/02/2002  RNO.
 *
 * Does not rely on the Starlink heliocentric and barycentric position
 * and velocity routines.  This is system independent.
*/

#include <stdio.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaJD.h"
#include "sofa/sofa.h"
#include "sofa/sofam.h"

/*! \file operaJDTest.cpp */

/*! 
 * operaJDTest
 * \author Doug Teeple
 * \brief Perform various tests on Julian dates.
 * \arg argc
 * \arg argv
 * \note --keyword=...
 * \return EXIT_STATUS
 * \ingroup test
 */

int main(void)
{
	
	double jd;				/* Julian date*/
	double mjd;				/* Modified Julian date */
	double hjd;				/* Heliocentric Julian date */
	double hmjd;			/* Heliocentric modified Julian date */
	struct tm *ltime;		/* Time and timezone structure */
	time_t system_time;		/* The system time in seconds from 1970 */
	int early_date = 0;		/* A flag whether the date is before 1582 or not */
	
	double djm0 = 0.0, djm = 0.0; /* SOFA Julian Date */
	
	system_time = time(NULL);
	ltime = localtime(&system_time);
	if (ltime == NULL) {
		printf("Error ltime is NULL, aborting...\n");
		exit(EXIT_FAILURE);
	}
	
	/* Issue a warning if the date is before Oct 15 1582 */
	
	if (ltime->tm_year < -318) {
		early_date = 1;
	} else {
		if (ltime->tm_year == -318) {
			if (ltime->tm_mon < 9)
				early_date = 1;
			if (ltime->tm_mon == 9) {
				if (ltime->tm_mday < 16)
					early_date = 1;
			}
		}
    }
	
	if (early_date == 1) {
		printf("The Gregorian calender is inacurate before 15 Oct 1582.\n");
		exit(EXIT_FAILURE);
	}
	
	/* Output the time and date specified, to check it is what the user wanted.  */
	
    printf("Time: %02i:%02i:%02i %i/%i/%i\n", ltime->tm_hour,
		   ltime->tm_min, ltime->tm_sec, ltime->tm_mday, ltime->tm_mon+1,
		   ltime->tm_year+1900);
	
    /* Calculate the SOFA Julian date */
	
	int status = 0;
	if ((status=iauCal2jd(ltime->tm_year+1900, ltime->tm_mon+1, ltime->tm_mday, &djm0, &djm))!=0)
		printf("SOFA JD Failure %d...\n", status);
	
	/* Now add in the time to the SOFA Julian Date */
	djm += ltime->tm_hour / 24.0;
	djm += ltime->tm_min / 24.0 / 60.0;
	djm += ltime->tm_sec / 24.0 / 60.0 / 60.0;
	
    /* Calculate the Julian date */
    
    jd = Julian_date(ltime);
	
    /* Calculate the modified Julian date */
	
    mjd = modified_Julian_date(jd);
	
    /* Calculate the modified heliocentric Julian date */
	
	float MaunaKeaRA = 10.36478;
	float MaunaKeaDec =  19.8267;
    hmjd = heliocentric_modified_Julian_date(mjd, MaunaKeaRA, MaunaKeaDec);
	
    /* Calculate the heliocentric Julian date */
	
    hjd = hmjd + 2400000.5;
	
    /* Print the results */
	
    printf("SOFA                      Julian Date = %lf\n", djm0+djm);
    printf("SOFA             Modified Julian Date =   %lf\n", djm);
    printf("                          Julian Date = %lf\n", jd);
    printf("                 Modified Julian Date =   %f\n", mjd);
    printf("    Heliocentric Modified Julian Date =   %f\n", hmjd);
    printf("             Heliocentric Julian Date = %lf\n", hjd);
	
	exit(EXIT_SUCCESS);
	
}


