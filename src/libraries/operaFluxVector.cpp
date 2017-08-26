/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: oepraFluxvector
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFluxVector.h"
#include <cmath>

/*!
 * \brief This file holds the implementation of the class operaFluxVector.
 * \file operaFluxVector.cpp
 * \ingroup libraries
 */

operaFluxVector::operaFluxVector(TendsTowards_t towards) : towards(towards) { }

operaFluxVector::operaFluxVector(unsigned length, TendsTowards_t towards) : towards(towards) {
	flux.resize(length);
	variance.resize(length);
}

operaFluxVector::operaFluxVector(double *fluxes, double *variances, unsigned length, TendsTowards_t towards) : towards(towards), flux(fluxes, length), variance(variances, length) { }

operaFluxVector::operaFluxVector(const operaVector& fluxes, const operaVector& variances, TendsTowards_t towards) : towards(towards), flux(fluxes), variance(variances) { }

operaFluxVector::operaFluxVector(const operaFluxVector &b, TendsTowards_t towards) : towards(towards), flux(b.flux), variance(b.variance) { }

operaFluxVector::operaFluxVector(const operaVector &b, TendsTowards_t towards) : towards(towards), flux(b), variance(b.size()) { }

void operaFluxVector::clear() {
	flux.clear();
	variance.clear();
}

void operaFluxVector::trim(operaIndexRange range) {
	flux.trim(range);
	variance.trim(range);
}

void operaFluxVector::resize(unsigned newlength) {
	flux.resize(newlength);
	variance.resize(newlength);
}

void operaFluxVector::insert(double newflux, double newvariance) {
	flux.insert(newflux);
	variance.insert(newvariance);
}

void operaFluxVector::reverse() {
	flux.reverse();
	variance.reverse();
}

void operaFluxVector::reorder(const operaIndexMap& indexmap) {
	flux.reorder(indexmap);
	variance.reorder(indexmap);
}

const operaVector &operaFluxVector::getflux() const {
	return flux;
}

const operaVector &operaFluxVector::getvariance() const {
	return variance;
}

void operaFluxVector::setflux(const operaVector& newflux) {
	if (flux.size() != newflux.size()) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	flux = newflux;
}

void operaFluxVector::setvariance(const operaVector& newvariance) {
	if (variance.size() != newvariance.size()) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	variance = newvariance;
}

void operaFluxVector::setflux(double newflux) {
	flux = newflux;
}

void operaFluxVector::setvariance(double newvariance) {
	variance = newvariance;
}

double* operaFluxVector::getfluxpointer() {
	return flux.datapointer();
}

const double* operaFluxVector::getfluxpointer() const {
	return flux.datapointer();
}

double* operaFluxVector::getvariancepointer() {
	return variance.datapointer();
}

const double* operaFluxVector::getvariancepointer() const {
	return variance.datapointer();
}

double operaFluxVector::getflux(unsigned index) const {
	return flux[index];
}

double operaFluxVector::getvariance(unsigned index) const {
	return variance[index];
}

void operaFluxVector::setflux(double newflux, unsigned index) {
	flux[index] = newflux;
}

void operaFluxVector::setvariance(double newvariance, unsigned index) {
	variance[index] = newvariance;
}

unsigned operaFluxVector::getlength() const {
	return flux.size();
}

TendsTowards_t operaFluxVector::gettowards(void) const {
	return towards;
}

double operaFluxVector::geterror(unsigned index) const {
	return sqrt(variance[index]);
}

std::pair<double,double> operaFluxVector::operator[](unsigned index) const {
	return std::pair<double,double>(flux[index], variance[index]);
}

operaFluxVector& operaFluxVector::operator+=(const operaFluxVector& b) {
	if (flux.size() != b.flux.size()) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	variance += b.variance;
	flux += b.flux;
	return *this;
}

operaFluxVector& operaFluxVector::operator-=(const operaFluxVector& b) {
	if (flux.size() != b.flux.size()) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	variance += b.variance;
	flux -= b.flux;
	return *this;
}

operaFluxVector& operaFluxVector::operator*=(const operaFluxVector& b) {
	if (flux.size() != b.flux.size()) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	variance = flux * flux * b.variance + b.flux * b.flux * variance;
	flux *= b.flux;
	return *this;
}

operaFluxVector& operaFluxVector::operator/=(const operaFluxVector& b) {
	if (flux.size() != b.flux.size()) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	for(unsigned i = 0; i < flux.size(); i++) {
		if (towards == ToDefault || !isinf(flux[i]) || !isinf(b.flux[i])) {
			const double b2 = b.flux[i]*b.flux[i]; //small optimization
			const double aoverb2 = flux[i]/b2; //small optimization
			variance[i] = variance[i]/b2 + aoverb2 * aoverb2 * b.variance[i]; // = (v(a)*b^2 + v(b)*a^2)/b^4
			flux[i] /= b.flux[i];
		} else if (towards == ToINF) {
			flux[i] = FP_INFINITE;
			variance[i] = 0.0;
		} else if (towards == ToNAN) {
			flux[i] = FP_NAN;
			variance[i] = FP_NAN;
		} else if (towards == ToZero) {
			flux[i] = 0.0;
			variance[i] = 0.0;
		} else if (towards == ToOne) {
			flux[i] = 1.0;
			variance[i] = 0.0;
		}
	}
	return *this;
}

operaFluxVector& operaFluxVector::operator=(double d) {
	flux.fill(d);
	variance.fill(0.0);
	return *this;
}

operaFluxVector& operaFluxVector::operator=(const operaVector& b) {
	flux = b;
	variance.resize(b.size());
	variance.fill(0.0);
	return *this;
}

operaFluxVector& operaFluxVector::operator+=(double d) {
	flux += d;
	return *this;
}

operaFluxVector& operaFluxVector::operator-=(double d) {
	flux -= d;
	return *this;
}

operaFluxVector& operaFluxVector::operator*=(double d) {
	flux *= d;
	variance *= d*d;
	return *this;
}

operaFluxVector& operaFluxVector::operator/=(double d) {
	flux /= d;
	variance /= d*d;
	return *this;
}

operaFluxVector operator-(double d, operaFluxVector a) {
	a.flux = d - a.flux;
	return a;
}

operaFluxVector operator/(double d, operaFluxVector a) {
	operaVector asqr = a.flux * a.flux;
	a.variance *= (d * d) / (asqr * asqr);
	a.flux = d / a.flux;
	return a;
}

operaFluxVector Sqrt(operaFluxVector b) {
	b.variance /= (b.flux * 4.0);
	b.flux = Sqrt(b.flux);
	return b;
}

operaFluxVector Pow(operaFluxVector b, double d) {
	b.variance *= Pow(Pow(b.flux, d-1.0)*d, 2.0);
	b.flux = Pow(b.flux, d);
	return b;
}

std::pair<double,double> Sum(const operaFluxVector& b) {
	return std::pair<double,double>(Sum(b.flux), Sum(b.variance));
}

std::pair<double,double> Mean(const operaFluxVector& b) {
	return std::pair<double,double>(Mean(b.flux), Mean(b.variance));
}
