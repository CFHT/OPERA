/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaHeliocentricWavelengthCorrection
 Version: 1.0
 Description: Apply Heliocentric velocity wavelength correction 
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
 Contact: opera@cfht.hawaii.edu
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see:
 http://software.cfht.hawaii.edu/licenses
 -or-
 http://www.gnu.org/licenses/gpl-3.0.html
 ********************************************************************/

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include "libraries/operaDateTime.h"
#include "libraries/operaIOFormats.h"
#include "libraries/operaArgumentHandler.h"
#include <cmath>
#include <fstream>

/*! \file operaHeliocentricWavelengthCorrection.cpp */

/*! 
 * operaHeliocentricWavelengthCorrection
 * \author Eder Martioli & Lison Malo
 * \brief Calculate and apply Heliocentric velocity wavelength correction.
 * \details Uses algorithms from the IDL Astronomy User's Library at http://idlastro.gsfc.nasa.gov/
 * \details Uses algorithms from the IRAF at http://iraf.noao.edu/
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \note --wave=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

using namespace std;

inline double HrsToDeg(double hours) { return hours * 15.0; }
inline double DegToHrs(double degrees) { return degrees / 15.0; }
inline double DegToRad(double degrees) { return degrees * PI / 180.0; }
inline double HrsToRad(double hours) { return hours * PI / 12.0; }
inline double AsecToDeg(double arcseconds) { return arcseconds / 3600.0; }

/*
 * \brief Calculates the Heliocentric Modified Julian Date
 * \param mjdate Modified Julian Date
 * \param ra Right ascension in radians
 * \param dec Declination in radians
 * \return The Heliocentric Modified Julian Date
 */
double CalculateHMJD(double mjdate, double ra, double dec);

/*
 * \brief Calculates the radial velocity component of Earth's rotation
 * \param latitude_degrees Lattitude of the observer in degrees
 * \param longitude_degrees Longitude of the observer in degrees
 * \param elevation_meters Elevation of the observer in meters
 * \param ra_radians Right ascension in radians
 * \param dec_radians Declination in radians
 * \return The diurnal radial velocity
 */
double DiurnalVelocity(double latitude_degrees, double longitude_degrees, double elevation_meters, double ra_radians, double dec_radians, double JDTime);

/*
 * \brief Calculates the radial velocity component of the Earth-Moon barycenter's orbit around the Sun
 * \param ra_radians Right ascension in radians
 * \param dec_radians Declination in radians
 * \param JDTime The heliocentric Julian date
 * \return The orbital radial velocity
 */
double OrbitalVelocity(double ra_radians, double dec_radians, double JDTime);

/*
 * \brief Calculates the radial velocity component of Earth's orbit around the Earth-Moon barycenter
 * \param ra_radians Right ascension in radians
 * \param dec_radians Declination in radians
 * \param JDTime The heliocentric Julian date
 * \return The orbital radial velocity
 */
double LunarVelocity(double ra_radians, double dec_radians, double JDTime);

unsigned NumberOfLeapSecondsElapsed(double mjd, string leapSecondList) {
	ifstream fin(leapSecondList.c_str());
	unsigned i = 0;
	for(Date leapSecondDate; fin >> leapSecondDate; i++) {
		double leapSecondMJD = JDtoMJD(leapSecondDate.ToJulianDayNumber() + 0.5); //Add 0.5, since leap second happens at the end of that day rather than noon
		if(mjd < leapSecondMJD) break;
	}
	return i;
}

int main(int argc, char *argv[]) {
    string modulename = "operaHeliocentricWavelengthCorrection";
    operaArgumentHandler args;
    
    string outputRVelFile;
    double MJDTime = 0.0;
    string object_coords_s = "0.0 0.0";
    string observatory_coords_s = "19:49:41.86 -155:28:18.00";
    double observatory_elevation = 4200;
    double etime=0.0;
    string leapseconds;
    
    args.AddRequiredArgument("outputRVelFile", outputRVelFile, "output radial velocity correction file (.rvel)");
    args.AddRequiredArgument("MJDTime", MJDTime, "time at the exposure start in modified Julian date");
    args.AddRequiredArgument("object_coords", object_coords_s, "object sky coordinates \"RA Dec\"");
    args.AddRequiredArgument("observatory_coords", observatory_coords_s, "observatory geographic coordinates \"latitude longitude\"");
    args.AddRequiredArgument("observatory_elevation", observatory_elevation, "observatory elevation in meters");
    args.AddRequiredArgument("etime", etime, "exposure time (shutter open)");
    args.AddRequiredArgument("leapseconds", leapseconds, "file with a list of the dates of each leap second");
    
	try {
		args.Parse(argc, argv);
		
		// we need an outputRVelFile file...
		if (outputRVelFile.empty()) {
			throw operaException(modulename+": ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		
		double dec_degrees = 0.0, ra_degrees = 0.0;
		if (!object_coords_s.empty()) {
			istringstream ss(object_coords_s);
			ss >> ra_degrees >> dec_degrees;
		}
		
		double lat_degrees = 0.0, long_degrees = 0.0;
		if (!observatory_coords_s.empty()) {
			istringstream ss(observatory_coords_s);
			Sexigesimal observatory_latitude, observatory_longitude;
			ss >> observatory_latitude >> observatory_longitude;
			lat_degrees = observatory_latitude.ToDecimal();
			long_degrees = abs(observatory_longitude.ToDecimal());
		}
		
		if (args.verbose) {
			cout.precision(6);
			cout << fixed;
			cout << modulename << ": outputRVelFile = " << outputRVelFile << endl;
			cout << modulename << ": MJDTime = " << MJDTime << endl;
			cout << modulename << ": sky coordinates RA = " << ra_degrees << " deg, Dec = " << dec_degrees << " deg" << endl;
			cout << modulename << ": geographic coordinates Latitude = " << lat_degrees  << "deg, Longitude = " << long_degrees << " deg\n";
            cout << modulename << ": observatory_elevation = " << observatory_elevation << " m" << endl;
            cout << modulename << ": etime = " << etime << endl;
            cout << modulename << ": leapseconds file = " << leapseconds << endl;
		}
		
		MJDTime += SecToDays(etime/2.0); //Divide by 2 to get the middle of the exposure, and convert from seconds to days
        if (args.verbose) {
            cout << modulename << ": mid-exposure etime = " << etime/2.0 << "s" << endl;
            cout << modulename << ": mid-exposure MJD = " << MJDTime << endl;
        }
		
		double ra_radians = DegToRad(ra_degrees);
		double dec_radians = DegToRad(dec_degrees);
		
		double hmjd = CalculateHMJD(MJDTime, ra_radians, dec_radians);
		double hjdtime_utc = MJDtoJD(hmjd);
		unsigned nleapseconds = NumberOfLeapSecondsElapsed(MJDTime, leapseconds);
		double hjdtime_tt = hjdtime_utc + SecToDays(nleapseconds + 10 + 32.184); //Add leap seconds + 10 to convert UTC to TAI, and add 32.184 more to convert to TT.
		
		double vdiurnal = DiurnalVelocity(lat_degrees, long_degrees, observatory_elevation, ra_radians, dec_radians, hjdtime_tt);
		double vorbital = OrbitalVelocity(ra_radians, dec_radians, hjdtime_tt);
		double vlunar = LunarVelocity(ra_radians, dec_radians, hjdtime_tt);
		double rvCorrection = vlunar + vorbital + vdiurnal;
		
		if (args.verbose) {
			cout << modulename << ": HJD UTC = " << hjdtime_utc << endl;
            cout << modulename << ": HJD TT = " << hjdtime_tt << endl;
			cout << modulename << ": vlunar = " << vlunar << endl;
			cout << modulename << ": vorbital = " << vorbital << endl;
			cout << modulename << ": vdiurnal = " << vdiurnal << endl;
			cout << modulename << ": Total RV correction is " << rvCorrection << endl;
		}
		
		FormatHeader outputheader("Radial Velocity Correction");
		outputheader << "radialvelocity (km/s)" << "lunar rvel" << "orbital rvel" << "diurnal rvel" << "HJD UTC" << "HJD TT" << newline;
		FormatData outputdata;
		outputdata << rvCorrection << vlunar << vorbital << vdiurnal << fixed << hjdtime_utc << hjdtime_tt << endl;
		operaIOFormats::WriteCustomFormat("rvel", outputheader, outputdata, outputRVelFile);
	}
	catch (operaException e) {
		cerr << modulename << ": " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << modulename << ": " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/*
 * \brief Calculate positions at B1950.0 from positions at J2000.0
 * \param ra Right acension at J2000.0 in radians 
 * \param dec Declination at J2000.0 in radians 
 * \param ra_1950 Output right acension at B1950.0 in radians
 * \param dec_1950 Output declination at B1950.0 in radians 
 */
void bprecess(double ra, double dec, double &ra_1950, double &dec_1950);

/*
 * \brief Calculate geocentric X,Y,Z coordinates of the Sun
 * \param rjdate Reduced Julian Date
 * \param x Output x-coordinate of the Sun
 * \param y Output y-coordinate of the Sun
 * \param z Output z-coordinate of the Sun
 */
void xyz(double rjdate, double &x, double &y, double &z);

double CalculateHMJD(double mjdate, double ra, double dec) {
	double date = mjdate + 0.5; //Convert from MJD to RJD
	bprecess(ra, dec, ra, dec); //Because xyz uses default B1950 coordinates, we'll convert everything to B1950
	
	double delta_t = (date - 33282.42345905) / 36525.0;
	double epsilon_sec = 44.836 - 46.8495 * delta_t - 0.00429 * pow(delta_t, 2) + 0.00181 * pow(delta_t, 3);
	double epsilon = DegToRad(23.433333 + epsilon_sec / 3600.0);
	
	double x, y, z;
	xyz(date, x, y, z);
	
	// Find extra distance light must travel in AU, multiply by 1.49598e13 cm / AU, and divide by the speed of light
	double time = -499.00522*(cos(dec)*cos(ra)*x + (tan(epsilon)*sin(dec) + cos(dec)*sin(ra))*y);
	return date + SecToDays(time) - 0.5; //Subtract 0.5 to convert back from RJD to MJD
}

/*
 * \brief Converts angular coordiantes to new origin and pole
 * \param ao, bo Origin of new coordinates in radians
 * \param ap, bp Pole of new coordinates radians
 * \param a1, b1 Coordinates to be converted in radians
 * \param a2, b2 Converted coordinates in radians
 */
void ast_coord(double ao, double bo, double ap, double bp, double a1, double b1, double &a2, double &b2);

double DiurnalVelocity(double latitude_degrees, double longitude_degrees, double elevation_meters, double ra_radians, double dec_radians, double JDTime) {
	double latitude_radians = DegToRad(latitude_degrees);
	//Reduction of geodetic latitude to geocentric latitude.
	double dlat = -1.0*(11.0 * 60.0 + 32.743000)*sin(2.0 * latitude_radians) + 1.163300*sin(4.0 * latitude_radians) - 0.002600*sin(6.0 * latitude_radians); //In arcseconds
	latitude_radians += DegToRad(dlat/3600.0); //arcseconds to radians
	
	//r is the radius vector from the Earth's center to the observer (meters).
	double r = 6378160.0 * (0.998327073 + 0.00167643800 * cos(2.0 * latitude_radians) - 0.00000351 * cos(4.0 * latitude_radians) + 0.000000008 * cos(6.0 * latitude_radians)) + elevation_meters;
	const double sideral_day_in_hours = 23.934469591229; //(1986)
	double vc = TWOPI * (r / 1000.0)  / (sideral_day_in_hours * 3600.0); //The corresponding circular velocity (meters/sidereal day converted to km / sec).
	
	double dzero = JDTime-2451545.0; //Time elapsed since Noon UT on Jan 1, 2000
	
	double gmst = 18.697374558 + 24.06570982441908*dzero; //Greenwich mean sidereal time
	double lmst = fmod(gmst + DegToHrs(longitude_degrees), 24.0); //Local mean sidereal time
	
	return vc * cos(latitude_radians) * cos(dec_radians) * sin(ra_radians - HrsToRad(lmst));
}

double OrbitalVelocity(double ra_radians, double dec_radians, double JDTime) {
	double t = (JDTime - 2415020.0) / 36525.0; //The number of Julian centuries since J1900.
	
	double manom = 358.47583 + t * (35999.04975 - t * (0.000150 + t * 0.000003)); //The mean anomaly of the Earth's orbit (degrees)
	double lperi = 101.22083 + t * (1.7191733 + t * (0.000453 + t * 0.000003)); //The the mean longitude of perihelion (degrees)
	double oblq = 23.452294 - t * (0.0130125 + t * (0.00000164 - t * 0.000000503)); //The mean obliquity of the ecliptic (degrees)
	double eccen = 0.01675104 - t * (0.00004180 + t * 0.000000126); //The eccentricity of the Earth's orbit (dimensionless)

	//Convert to principle angles, in radians
	manom = DegToRad(fmod(manom, 360.0));
	lperi = DegToRad(fmod(lperi, 360.0));
	oblq = DegToRad(oblq);

	double tanom = manom + (2.0 * eccen - 0.25 * pow(eccen, 3)) * sin (manom) + 1.25 * pow(eccen, 2) * sin (2.0 * manom) + 13.0/12.0 * pow(eccen, 3) * sin (3.0 * manom); //The true anomaly (approximate formula) (radians)
	double slong = lperi + tanom + PI; //The true longitude of the Sun seen from the Earth (radians)

	double l, b; //The longitude and latitude of the star in the orbital plane of the Earth (radians)
	ast_coord(0.0, 0.0, -PI_OVER_2, PI_OVER_2 - oblq, ra_radians, dec_radians, l, b);
	
	const double semi_major_axis = 149598500.0; //The Earth's semi-major axis in km
	const double year_in_days = 365.2564;
	double vorb = SecToDays((TWOPI / year_in_days) * semi_major_axis / sqrt (1.0 - pow(eccen, 2))); //The component of the Earth's orbital velocity perpendicular to the radius vector (km/s)

	//Returns the projection onto the line of sight to the observation of the velocity of the Earth-Moon barycenter with respect to the Sun (km/s).
	return vorb * cos(b) * (sin(slong - l) - eccen * sin(lperi - l));
}

double LunarVelocity(double ra_radians, double dec_radians, double JDTime) {
	double t = (JDTime - 2415020.0) / 36525.0; //The number of Julian centuries since J1900.
	
	double oblq = DegToRad(23.452294 - t * (0.0130125 + t * (0.00000164 - t * 0.000000503))); //The mean obliquity of the ecliptic
	double omega = DegToRad(259.183275 - t * (1934.142008 + t * (0.002078 + t * 0.000002))); //The longitude of the mean ascending node
	double llong = DegToRad(270.434164 + t * (481267.88315 + t * (-0.001133 + t * 0.0000019))) - omega; //The mean lunar longitude (should be 13.1763965268 deg)
	double lperi = DegToRad(334.329556 + t * (4069.034029 - t * (0.010325 + t * 0.000012))) - omega; //The mean lunar longitude of perigee
	double inclin = DegToRad(5.1453964); //The inclination of the lunar orbit to the ecliptic
	double em = 0.054900489; //The eccentricity of the lunar orbit (dimensionless)
	
	//Determine true longitude.  Compute mean anomaly, convert to true anomaly (approximate formula), and convert back to longitude.
	//The mean anomaly is only approximate because lperi should be the true rather than the mean longitude of lunar perigee.
	double anom = llong - lperi;
	anom = anom + (2.0 * em - 0.25 * pow(em, 3)) * sin(anom) + 1.25 * pow(em, 2) * sin(2.0 * anom) + 13.0/12.0 * pow(em, 3) * sin(3.0 * anom);
	llong = anom + lperi;
	
	double l, b; //The ecliptic longitude and latitude of the observation.
	ast_coord(0.0, 0.0, -PI_OVER_2, PI_OVER_2 - oblq, ra_radians, dec_radians, l, b);
	
	double lm, bm; //The lunar longitude and latitude of the observation in the lunar orbital plane relative to the ascending node.
	ast_coord(omega, 0.0, omega - PI_OVER_2, PI_OVER_2 - inclin, l, b, lm, bm);
	
	double vmoon = SecToDays((TWOPI / 27.321661) * 384403.12040 / sqrt(1.0 - pow(em, 2))); //The component of the lunar velocity perpendicular to the radius vector.
	const double m_ratio = 81.53; //The ratio of the Earth's mass to the Moon's mass.
	
	//Returns the projection onto the line of sight to the observation of the velocity of the Earth's center with respect to the Earth-Moon barycenter.
	return vmoon * cos (bm) * (sin (llong - lm) - em * sin (lperi - lm)) / m_ratio;
}

void bprecess(double ra, double dec, double &ra_1950, double &dec_1950) {
	double A_arr[3] = { -1.62557e-6, -0.31919e-6, -0.13843e-6 }; //in radians
	double A_dot_arr[3] = { 1.244e-3, -1.579e-3, -0.660e-3 }; //in arc seconds per century
	operaVector A(A_arr, 3);
	operaVector A_dot(A_dot_arr, 3);
	A_dot *= DegToRad(AsecToDeg(1.0)); //convert to radians

	double r0_arr[3] = { cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec) };
	operaVector r0(r0_arr, 3);
	
	// Include the effects of the E-terms of aberration to form r and r_dot.
	double M[6][3] = {
		{0.9999256795, 0.0111814828, 0.0048590039},
		{-0.0111814828, 0.9999374849, -0.0000271771},
		{-0.0048590040, -0.0000271557, 0.9999881946},
		{-0.000551, 0.238509, -0.435614},
		{-0.238560, -0.002667, 0.012254},
		{0.435730, -0.008541, 0.002117}
	};
	operaVector r1(3);
	operaVector r1_dot(3);
	for (unsigned i = 0; i < 3; i++) {
		r1[i] = InnerProduct(r0, operaVector(M[i], 3));
		r1_dot[i] = InnerProduct(r0, operaVector(M[i + 3], 3));
	}
	r1_dot *= DegToRad(AsecToDeg(1.0)); //convert to radians
	
	double rmag = Magnitude(r1);

	r1 += r1_dot * 0.5;
	A += A_dot * 0.5;
	
	operaVector s1 = r1/rmag;
	operaVector s1_dot = r1_dot/rmag;

	operaVector s = s1;
	operaVector r;
	for (unsigned j = 0; j < 3; j++) {
		r = s1 + A - (s * InnerProduct(A, s));
		s = r/rmag;
	}
	rmag = Magnitude(r);

	double x = r[0], y = r[1], z = r[2];
	dec_1950 = asin(z / rmag);
	ra_1950 = atan2(y, x);

	if (ra_1950 < 0) ra_1950 = ra_1950 + TWOPI;
}

void xyz(double rjdate, double &x, double &y, double &z) {
	double t = (rjdate - 15020.0) / 36525.0; //Relative Julian century from 1900

	// NOTE: longitude arguments below are given in *equinox* of date.
	// Precess these to equinox 1950 to give everything an even footing.
	double pp = (1.396041 + 0.000308*(t + 0.5))*(t - 0.499998); //Compute argument of precession from equinox of date back to 1950
	double el = 279.696678 + 36000.76892*t + 0.000303*t*t - pp; // Compute mean solar longitude, precessed back to 1950
	double c = 270.434164 + 480960.0*t + 307.883142*t - 0.001133*t*t - pp; //Compute Mean longitude of the Moon
	double n = 259.183275 - 1800.0*t - 134.142008*t + 0.002078*t*t - pp; // Compute longitude of Moon's ascending node
	double g = 358.475833 + 35999.04975*t - 0.00015*t*t; // Compute mean solar anomaly
	double j = 225.444651 + 2880.0*t + 154.906654*t*t; //Compute the mean jupiter anomaly
	double v = 212.603219 + 58320.0*t + 197.803875*t + 0.001286*t*t; //Compute mean anomaly of Venus
	double m = 319.529425 + 19080.0*t + 59.8585*t + 0.000181*t*t; //Compute mean anomaly of Mars
	// Convert degrees to radians for trig functions
	el = DegToRad(el);
	c = DegToRad(c);
	v = DegToRad(v);
	g = DegToRad(g);
	j = DegToRad(j);
	n = DegToRad(n);
	m = DegToRad(m);

	// Calculate X, Y, Z using trigonometric series
	x = 0.999860*cos(el)
		- 0.025127*cos(g - el)
		+ 0.008374*cos(g + el)
		+ 0.000105*cos(g + g + el)
		+ 0.000063*t*cos(g - el)
		+ 0.000035*cos(g + g - el)
		- 0.000026*sin(g - el - j)
		- 0.000021*t*cos(g + el)
		+ 0.000018*sin(2.0*g + el - 2.0*v)
		+ 0.000017*cos(c)
		- 0.000014*cos(c - 2.0*el)
		+ 0.000012*cos(4.0*g + el - 8.0*m + 3.0*j)
		- 0.000012*cos(4.0*g - el - 8.0*m + 3.0*j)
		- 0.000012*cos(g + el - v)
		+ 0.000011*cos(2.0*g + el - 2.0*v)
		+ 0.000011*cos(2.0*g - el - 2.0*j);

	y = 0.917308*sin(el)
		+ 0.023053*sin(g - el)
		+ 0.007683*sin(g + el)
		+ 0.000097*sin(g + g + el)
		- 0.000057*t*sin(g - el)
		- 0.000032*sin(g + g - el)
		- 0.000024*cos(g - el - j)
		- 0.000019*t*sin(g + el)
		- 0.000017*cos(2.0*g + el - 2.0*v)
		+ 0.000016*sin(c)
		+ 0.000013*sin(c - 2.0*el)
		+ 0.000011*sin(4.0*g + el - 8.0*m + 3.0*j)
		+ 0.000011*sin(4.0*g - el - 8.0*m + 3.0*j)
		- 0.000011*sin(g + el - v)
		+ 0.000010*sin(2.0*g + el - 2.0*v)
		- 0.000010*sin(2.0*g - el - 2.0*j);


	z = 0.397825*sin(el)
		+ 0.009998*sin(g - el)
		+ 0.003332*sin(g + el)
		+ 0.000042*sin(g + g + el)
		- 0.000025*t*sin(g - el)
		- 0.000014*sin(g + g - el)
		- 0.000010*cos(g - el - j);
}

void ast_coord(double ao, double bo, double ap, double bp, double a1, double b1, double &a2, double &b2) {
	double x = cos (a1) * cos (b1);
	double y = sin (a1) * cos (b1);
	double z = sin (b1);
	double xp = cos (ap) * cos (bp);
	double yp = sin (ap) * cos (bp);
	double zp = sin (bp);

	//Rotate the origin about z.
	double sao = sin (ao);
	double cao = cos (ao);
	double sbo = sin (bo);
	double cbo = cos (bo);
	double temp = -1*xp * sao + yp * cao;
	xp = xp * cao + yp * sao;
	yp = temp;
	temp = -1*x * sao + y * cao;
	x = x * cao + y * sao;
	y = temp;

	//Rotate the origin about y.
	temp = -1*xp * sbo + zp * cbo;
	xp = xp * cbo + zp * sbo;
	zp = temp;
	temp = -1*x * sbo + z * cbo;
	x = x * cbo + z * sbo;
	z = temp;

	//Rotate pole around x.
	double sbp = zp;
	double cbp = yp;
	temp = y * cbp + z * sbp;
	y = y * sbp - z * cbp;
	z = temp;

	//Final angular coordinates.
	a2 = atan2(y, x);
	b2 = asin(z);
}
