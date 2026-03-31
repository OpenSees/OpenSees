/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Jose Abell (UANDES, github.com/jaabell)
// Massimo Petracca - ASDEA Software, Italy (2022)
//
// Utility header file to interface with the Eigen library
//

#ifndef EigenAPI_h
#define EigenAPI_h

#include <OPS_Globals.h>
#include <Matrix.h>
#include <Vector.h>

// MACROS that should be defined before including ANY eigen file

/*
The macro EIGEN_DONT_ALIGN is deprecated (it is a synonim of EIGEN_MAX_ALIGN_BYTES=0).
It disables automatic allignment of fixed size arrays.
It is necessary to do so to avoid problems using those types as member of structures.
More information here:
https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
*/

#define EIGEN_MAX_ALIGN_BYTES 0
#define EIGEN_MAX_CPP_VER 14


#include "Eigen/Dense"

// namespace EigenAPI {

#include "converters.h"
#include "typedefs.h"
#include "operations.h"

// }

#endif // EigenAPI_h