#ifndef OPERADATETIME_H
#define OPERADATETIME_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaDateTime
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Mar/2012
 
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

#include <string>
#include <iostream>

/*!
 * \brief This class encapsulates a sexigesimal number, mainly for representing hours or degrees.
 * \ingroup libraries
 * \details
 * This class contains a whole number count, a count of 60ths and a count of 3600ths.
 * There is also a sign to determine positive or negative values.
 */
class Sexigesimal {
	public:
		/*!
		 * \brief Sets the value from a string.
		 * \param str In format "aa:bb:cc.ccc", with an optional + or - prefix.
		 * \return Whether the value was successfuly set.
		 */
		bool SetFromString(std::string str);
		
		/*!
		 * \brief Converts the fractional components to units of the most significant value and combines them.
		 * \return The decimal representation of the number.
		 */
		double ToDecimal() const;
		
		friend std::istream& operator>>(std::istream& in, Sexigesimal& var);
		friend std::ostream& operator<<(std::ostream& out, const Sexigesimal& var);
	protected:
		int sign;
		unsigned units;
		unsigned minutes;
		double seconds;
};

class Time : private Sexigesimal {
	public:
		/*!
		 * \brief Sets the time from a string.
		 * \param str In format "hh:mm:ss.sss", with an optional + or - prefix.
		 * \return Whether the value was successfuly set.
		 */
		bool SetFromString(std::string str);
		
		/*!
		 * \brief Converts the time to number of seconds.
		 * \return The time in days.
		 */
		double ToDays() const;
		
		/*!
		 * \brief Converts the time to number of hours.
		 * \return The time in days.
		 */
		double ToHours() const;
		
		/*!
		 * \brief Converts the time to number of minutes.
		 * \return The time in days.
		 */
		double ToMinutes() const;
		
		/*!
		 * \brief Converts the time to number of seconds.
		 * \return The time in days.
		 */
		double ToSeconds() const;
};

/*!
 * \brief This class encapsulates a date.
 * \ingroup libraries
 */
class Date {
	public:
		/*!
		 * \brief Sets the value from a string.
		 * \param str In format "yyyy-mm-dd".
		 * \return Whether the value was successfuly set.
		 */
		bool SetFromString(std::string str);
		
		/*!
		 * \brief Converts the Gregorian calendar date to the Julian day number.
		 * \details This is equal to the Julian date on that day at noon.
		 * \return The JDN.
		 */
		unsigned ToJulianDayNumber() const;
		
		friend std::istream& operator>>(std::istream& in, Date& var);
		friend std::ostream& operator<<(std::ostream& out, const Date& var);
	private:
		int year;
		int month;
		int day;
};

/*!
 * \brief This class encapsulates a date and time.
 * \ingroup libraries
 */
class DateTime {
	public:
		bool SetFromString(std::string str);
		bool SetFromStrings(std::string datestr, std::string timestr);
		double ToJulianDate() const;
	private:
		Date date;
		Time time;
};

double JDtoMJD(double jd);

double MJDtoJD(double mjd);

double SecToDays(double seconds);

#endif
