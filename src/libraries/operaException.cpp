/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaException
 Version: 1.0
 Description: class implements opera exceptions..
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
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

/*
 * This class defines an opera Exception.
 *
 */
#include <exception>
#include <string>

#include "fitsio.h"
#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaLib.h"					// for itos

/*!
 * \brief This class defines an opera Exception
 * \file operaException.cpp
 */

using namespace std;

/*!
 * operaException
 * \author Doug Teeple
 * \brief The opera exception class.
 * \ingroup libraries
 */

/* 
 * \class operaException()
 * \brief Basic Exception class constructor.
 * \return none
 */
operaException::operaException() : 
errorcode(operaErrorCodeOK),
line(0)
{
}

/* 
 * \class operaException
 * \brief operaException(operaErrorCode Errorcode)
 * \brief Basic Exception class constructor with an errorcode.
 * \param Errorcode
 * \return none
 */
operaException::operaException(operaErrorCode Errorcode) :
line(0)
{
	errorcode = Errorcode;
}

/* 
 * \class operaException
 * \brief operaException(string Message, operaErrorCode Errorcode)
 * \brief Basic Exception class constructor with an additional informatonal message and error code.
 * \param Message
 * \param Errorcode
 * \return none
 */
operaException::operaException(string Message, operaErrorCode Errorcode) : 
line(0)
{
	errorcode = Errorcode;
	message = Message;
}

/* 
 * \class operaException
 * \brief operaException(string Message, operaErrorCode Errorcode, string Filename, string Function, int Line)
 * \brief Basic Exception class constructor with an additional informatonal message, error code, and file, function, line.
 * \param Message
 * \param Errorcode
 * \param Filename
 * \param Function
 * \param Line
 * \return none
 */
operaException::operaException(string Message, operaErrorCode Errorcode, string Filename, string Function, int Line) {
	errorcode = Errorcode;
	message = Message;
	file = Filename;
	function = Function;
	line = Line;
}

/*
 * getters and setters
 */

/* 
 * operaException::getFormattedMessage()
 * \brief Return a string containing a fully formatted error message.
 * \return string
 */
string operaException::getFormattedMessage() {
	string output;

	if (!file.empty()) {
		output += " File: " + file;
		
	}
	if (!function.empty()) {
		output += " Function: " + function;
		
	}
	if (line != 0) {
		output += " Line: " + itos(line);
	}
	output += ": " + message + " " + operaStrError(errorcode);
	return output;
}

/* 
 * operaException::setMessage(string m)
 * \brief Set the string part of an error message.
 * \param m - message string
 * \return void
 */
void operaException::setMessage(string m) {
	message=m;
}

/* 
 * string operaException::getMessage()
 * \brief Return the string part of an error message.
 * \return string
 */
string operaException::getMessage(){
	return message;
}

/* 
 * operaException::setErrorCode(const operaErrorCode errcode)
 * \brief Set the error code part of an error message.
 * \param errcode const operaErrorCode
 * \return void
 */
void operaException::setErrorCode(const operaErrorCode errcode) {
	errorcode=errcode;
}

/* 
 * operaErrorCode operaException::getErrorCode()
 * \brief Return the error code part of an error message.
 * \return operaErrorCode
 */
operaErrorCode operaException::getErrorCode() {
	return errorcode;
}

/* 
 * operaException::setLine((int l)
 * \brief Set the line number part of an error message.
 * \param l int
 * \return void
 */
void operaException::setLine(int l){
	line=l;
}

/* 
 * string operaException::getLine()
 * \brief Return the line number part of an error message.
 * \return int
 */
int operaException::getLine(){
	return line;
}

/* 
 * operaException::setFunction(string func)
 * \brief Set the function name part of an error message.
 * \param func string
 * \return void
 */
void operaException::setFunction(string func){
	function=func;
}

/* 
 * string operaException::getFunction()
 * \brief Return the function name part of an error message.
 * \return string
 */
string operaException::getFunction(){
	return function;
}

/* 
 * operaException::setFile(string func)
 * \brief Set the file name part of an error message.
 * \param func string
 * \return void
 */
void operaException::setFile(string func){
	function=func;
}

/* 
 * string operaException::getFile()
 * \brief Return the file name part of an error message.
 * \return string
 */
string operaException::getFile(){
	return function;
}



