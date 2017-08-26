#ifndef OPERAEXCEPTION_H_
#define OPERAEXCEPTION_H_

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
#include <iostream>
#include <exception>
#include <string>

#include "operaError.h"

using namespace std;

/*! 
 * \class operaException
 * \brief extends exception
 * \brief Basic Exception class constructor.
 * \return none
 * \file operaException.h
 * \ingroup libraries
 */
class operaException: public exception {

private:

	operaErrorCode errorcode;			// the error code itself
	string file;						// optional file of originating exception
	int line;							// optional linr number of originating exception
	string function;					// optional function name of originating exception
	string message;						// optional message
	
public:

	/*! 
	 * \class operaException()
	 * \brief Basic Exception class constructor.
	 * \return none
	 */
	operaException();
	/*! 
	 * \class operaException
	 * \brief operaException(operaErrorCode Errorcode)
	 * \brief Basic Exception class constructor with an errorcode.
	 * \param Errorcode
	 * \return none
	 */
	operaException(operaErrorCode Errorcode);
	/*! 
	 * \class operaException
	 * \brief operaException(string Message, operaErrorCode Errorcode)
	 * \brief Basic Exception class constructor with an additional informatonal message and error code.
	 * \param Message
	 * \param Errorcode
	 * \return none
	 */
	operaException(string Message, operaErrorCode Errorcode);
	/*! 
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
	operaException(string Message, operaErrorCode Errorcode, string Filename, string Function, int Line);
	~operaException() throw() {}

	/*
	 * getters and setters
	 */
	
	/*! 
	 * operaException::getFormattedMessage()
	 * \brief Return a string containing a fully formatted error message.
	 * \return string
	 */
	string getFormattedMessage();
	
	/*! 
	 * operaException::setErrorCode(const operaErrorCode errcode)
	 * \brief Set the error code part of an error message.
	 * \return void
	 */
	void setErrorCode(const operaErrorCode errcode);
	/*! 
	 * string operaException::getErrorCode()
	 * \brief Return the error code part of an error message.
	 * \param errcode const operaErrorCode
	 * \return string
	 */
	operaErrorCode getErrorCode();
	
	/*! 
	 * string operaException::getMessage()
	 * \brief Return the string part of an error message.
	 * \param m - message string
	 * \return string
	 */
	void setMessage(string m);
	/*! 
	 * string operaException::getMessage()
	 * \brief Return the string part of an error message.
	 * \return string
	 */
	string getMessage();
	
	/*! 
	 * operaException::setLine((int l)
	 * \brief Set the line number part of an error message.
	 * \param l int
	 * \return void
	 */
	void setLine(int l);
	/*! 
	 * string operaException::getLine()
	 * \brief Return the line number part of an error message.
	 * \return int
	 */
	int getLine();
	
	/*! 
	 * operaException::setFunction(string func)
	 * \brief Set the function name part of an error message.
	 * \param func string
	 * \return void
	 */
	void setFunction(string func);
	/*! 
	 * string operaException::getFunction()
	 * \brief Return the function name part of an error message.
	 * \return string
	 */
	string getFunction();
	
	/*! 
	 * operaException::setFile(string func)
	 * \brief Set the file name part of an error message.
	 * \param func string
	 * \return void
	 */
	void setFile(string f);
	/*! 
	 * string operaException::getFile()
	 * \brief Return the file name part of an error message.
	 * \return string
	 */
	string getFile();
	
};



#endif /* OPERAEXCEPTION_H_ */
