#include "globaldefines.h"
#include "libraries/operaDateTime.h"
#include <sstream>

bool Sexigesimal::SetFromString(std::string str) {
	if(str.empty()) return false;
	std::istringstream ss(str);
	if(str[0] == '-' || str[0] == '+') ss.get();
	unsigned h2, m2;
	double s2;
	if(!(ss >> h2) || ss.get() != ':' || !(ss >> m2) || ss.get() != ':' || !(ss >> s2)) return false;
	sign = (str[0] == '-' ? -1 : 1);
	units = h2;
	minutes = m2;
	seconds = s2;
	return true;
}

double Sexigesimal::ToDecimal() const {
	return sign*(units + minutes/60.0 + seconds/3600.0);
}

std::istream& operator>>(std::istream& in, Sexigesimal& var) {
	std::string temp;
	in >> temp;
	if(!var.SetFromString(temp)) in.setstate(std::ios_base::failbit);
	return in;
}

std::ostream& operator<<(std::ostream& out, const Sexigesimal& var) {
	if(var.sign < 0) out << '-';
	out << var.units << ":" << var.minutes << ":" << var.seconds;
	return out;
}

bool Time::SetFromString(std::string str) {
	return Sexigesimal::SetFromString(str);
}

double Time::ToDays() const {
	return sign*(units/24.0 + minutes/1440.0 + seconds/86400.0);
}

double Time::ToHours() const {
	return sign*(units + minutes/60.0 + seconds/3600.0);
}

double Time::ToMinutes() const {
	return sign*(units*60 + minutes + seconds/60.0);
}

double Time::ToSeconds() const {
	return sign*(units*3600.0 + minutes*60.0 + seconds);
}

bool Date::SetFromString(std::string str) {
	if(str.empty()) return false;
	std::istringstream ss(str);
	int y2, m2, d2;
	if(!(ss >> y2) || ss.get() != '-' || !(ss >> m2) || ss.get() != '-' || !(ss >> d2)) return false;
	year = y2;
	month = m2;
	day = d2;
	return true;
}

unsigned Date::ToJulianDayNumber() const {
	const int a = (month - 14) / 12; // a = -1 for Jan and Feb, otherwise 0
	return ((1461 * (year + 4800 + a)) / 4) +
      ((367 * (month - 2 - (12 * a))) / 12) -
      ((3 * ((year + 4900 + a) / 100)) / 4) + day - 32075;
}

std::istream& operator>>(std::istream& in, Date& var) {
	std::string temp;
	in >> temp;
	if(!var.SetFromString(temp)) in.setstate(std::ios_base::failbit);
	return in;
}

std::ostream& operator<<(std::ostream& out, const Date& var) {
	out << var.year << ":" << var.month << ":" << var.day;
	return out;
}

bool DateTime::SetFromString(std::string str) {
	std::istringstream ss(str);
	std::string datestr, timestr;
	return std::getline(ss, datestr, 'T') && std::getline(ss, timestr) && SetFromStrings(datestr, timestr);
}

bool DateTime::SetFromStrings(std::string datestr, std::string timestr) {
	return time.SetFromString(timestr) && date.SetFromString(datestr);
}

double DateTime::ToJulianDate() const {
	return date.ToJulianDayNumber() + time.ToDays()  - 0.5;
}

const double mjd0 = 2400000.5; //The Julian date when MJD = 0

double JDtoMJD(double jd) {
	return jd - mjd0;
}

double MJDtoJD(double mjd) {
	return mjd + mjd0;
}

double SecToDays(double seconds) {
	return seconds / 86400.0;
}
