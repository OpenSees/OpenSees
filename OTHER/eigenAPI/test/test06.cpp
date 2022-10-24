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
// Test interface between eigen and opensees matrices.
//
// testing: copyToVectorReference

//
//
//
// types and operations
//
//


#include <iostream>
#include "../EigenAPI.h"



int main()
{
	VoigtVector sigma = {1, 1, 1, 1, 1, 1};
	std::cout << "sigma = " << sigma << std::endl;
	std::cout << "sigma.trace() = " << sigma.trace() << std::endl;
	std::cout << "sigma.deviator() = " << sigma.deviator() << std::endl;
}