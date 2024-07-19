/* *****************************************************************************
Copyright (c) 2012-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

#include "DL_Interpreter.h"
#include <iostream>

DL_Interpreter  *ops_TheActiveInterpreter = 0;

DL_Interpreter::DL_Interpreter()
{
    // does nothing
}

DL_Interpreter::~DL_Interpreter()
{
    // does nothing
}

int 
DL_Interpreter::addCommand(const char *, Command &)
{
    return -1;
}

int 
DL_Interpreter::removeCommand(const char *)
{
    return -1;
}

int 
DL_Interpreter::getNumRemainingInputArgs(void)
{
    return -1;
}

int 
DL_Interpreter::getInt(int *, int numArgs)
{
    return -1;
}

int 
DL_Interpreter::getDouble(double *, int numArgs)
{
    return -1;
}

int DL_Interpreter::getDoubleList(int* size, Vector* data)
{
    return -1;
}

const char*
DL_Interpreter::getString()
{
    return 0;
}

const char*
DL_Interpreter::getStringFromAll(char* buffer, int len)
{
    return 0;
}

int 
DL_Interpreter::getStringCopy(char **stringPtr)
{
    return -1;
}

int DL_Interpreter::evalDoubleStringExpression(const char* theExpression, double& current_val)
{
	return -1;
}

void
DL_Interpreter::resetInput(int cArg)
{
    // does nothing
}

int
DL_Interpreter::setInt(int *, int numArgs, bool scalar)
{
    return -1;
}

int DL_Interpreter::setInt(std::vector<std::vector<int>> &data) {
    return -1;
}

int DL_Interpreter::setInt(std::map<const char*, int>& data) {
    return -1;
}

int DL_Interpreter::setInt(std::map<const char*, std::vector<int>>& data) {
    return -1;
}

int
DL_Interpreter::setDouble(double *, int numArgs, bool scalar)
{
    return -1;
}

int
DL_Interpreter::setDouble(std::vector<std::vector<double>>& data)
{
    return -1;
}

int DL_Interpreter::setDouble(std::map<const char*, double>& data) {
    return -1;
}

int DL_Interpreter::setDouble(std::map<const char*, std::vector<double>>& data) {
    return -1;
}

int
DL_Interpreter::setString(const char*)
{
    return -1;
}

int
DL_Interpreter::setString(std::vector<const char*>& data)
{
    return -1;
}

int DL_Interpreter::setString(std::map<const char*, const char*>& data) {
    return -1;
}

int
DL_Interpreter::setString(std::vector<std::vector<const char*>>& data)
{
    return -1;
}

int DL_Interpreter::setString(std::map<const char*, std::vector<const char*>>& data) {
    return -1;
}

int
DL_Interpreter::runCommand(const char*)
{
    return -1;
}
