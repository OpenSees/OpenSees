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
                                                                        
// $Revision: 1.5 $
// $Date: 2007-07-11 23:51:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/HLRFSearchDirection.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef HLRFSearchDirection_h
#define HLRFSearchDirection_h

#include <SearchDirection.h>
#include <Vector.h>

class HLRFSearchDirection : public SearchDirection
{

public:
	HLRFSearchDirection();
	~HLRFSearchDirection();

	int computeSearchDirection(int stepNumber,
				   const Vector &passed_u, 
				   double passed_gFunctionValue, 
				   const Vector &passedGradientInStandardNormalSpace);
	const Vector &getSearchDirection();

protected:

private:
	Vector searchDirection;

};

#endif
