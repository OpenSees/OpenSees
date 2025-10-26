/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** ****************************************************************** */

// Written: Gustavo A. Araújo R. (gaaraujor@gmail.com)
//
// Description: Generic dictionary type for interpreter output
//
// This provides a type-safe way to return mixed-type dictionaries from
// OpenSees commands to both TCL and Python interpreters.
//
// Future consideration: This could be extended to support fully generic
// key-value pairs (e.g., map<DictValue, DictValue>) to allow non-string
// keys such as integer indices. The current implementation uses string
// keys to simplify the API.
//
// If extending to variant keys (map<DictValue, DictValue>), the following
// functions would need to be updated to handle non-string keys:
//   - PythonWrapper::setGenericDictOutput()  (SRC/interpreter/PythonWrapper.cpp)
//   - TclWrapper::setGenericDictOutput()     (SRC/interpreter/TclWrapper.cpp)
//   - OPS_SetGenericDict()                   (SRC/api/elementAPI_TCL.cpp)
// These currently use PyDict_SetItemString() and Tcl_NewStringObj() which
// assume string keys. For variant keys, you'd need PyDict_SetItem() and
// convert keys using the same std::visit pattern used for values.

#ifndef GenericDict_h
#define GenericDict_h

#include <string>
#include <map>
#include <vector>
#include <variant>

// Forward declarations
class Vector;
class ID;

// Type-safe variant for dictionary values
using DictValue = std::variant<
    int,
    double,
    std::vector<int>,
    std::vector<double>,
    std::string,
    std::vector<std::string>
>;

// Generic dictionary type for interpreter output
// Keys are strings, values can be any type in the DictValue variant
using GenericDict = std::map<std::string, DictValue>;

#endif

