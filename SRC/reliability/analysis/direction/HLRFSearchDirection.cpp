/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.6 $
// $Date: 2007-07-11 23:51:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/HLRFSearchDirection.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <HLRFSearchDirection.h>
#include <SearchDirection.h>
#include <Vector.h>


HLRFSearchDirection::HLRFSearchDirection()
:SearchDirection(), searchDirection(1)
{
}

HLRFSearchDirection::~HLRFSearchDirection()
{
}




const Vector&
HLRFSearchDirection::getSearchDirection()
{
	return searchDirection;
}



int
HLRFSearchDirection::computeSearchDirection(int stepNumber,
					    const Vector &u, 
					    double gFunctionValue, 
					    const Vector &gradientInStandardNormalSpace )
{

	// Compute the norm of the gradient
	double normOfGradient = gradientInStandardNormalSpace.Norm();


	// Check that the norm is not zero
	if (normOfGradient == 0.0) {
		opserr << "HLRFSearchDirection::computeSearchDirection() - " << endln
			<< " the norm of the gradient is zero. " << endln;
		return -1;
	}

	
	// Compute the alpha-vector
	//Vector alpha = gradientInStandardNormalSpace * ( (-1) / normOfGradient );
	Vector alpha(gradientInStandardNormalSpace);
	alpha *= -1.0/normOfGradient;


	// Compute the direction vector
	double alpha_times_u = alpha ^ u ;
	//Vector direction = alpha * ( gFunctionValue / normOfGradient + alpha_times_u ) - u;

	//searchDirection = direction;

	searchDirection = alpha;
	searchDirection.addVector(gFunctionValue / normOfGradient + alpha_times_u, u, -1.0);

	return 0;
}

