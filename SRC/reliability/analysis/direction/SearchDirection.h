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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-03-04 00:39:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/SearchDirection.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef SearchDirection_h
#define SearchDirection_h

#include <Vector.h>

class SearchDirection
{

public:
	SearchDirection();
	virtual ~SearchDirection();

	virtual int computeSearchDirection(	
			int stepNumber,
			Vector passed_u, 
			double passed_gFunctionValue, 
			Vector passedGradientInStandardNormalSpace) =0;
	virtual Vector getSearchDirection() =0;

protected:

private:

};

#endif
