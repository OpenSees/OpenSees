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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-10-27 23:45:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/PolakHeSearchDirectionAndMeritFunction.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef PolakHeSearchDirectionAndMeritFunction_h
#define PolakHeSearchDirectionAndMeritFunction_h

#include <SearchDirection.h>
#include <MeritFunctionCheck.h>
#include <Vector.h>

class PolakHeSearchDirectionAndMeritFunction : public SearchDirection, public MeritFunctionCheck
{

public:
	PolakHeSearchDirectionAndMeritFunction(double gamma, double delta);
	~PolakHeSearchDirectionAndMeritFunction();

	int computeSearchDirection(	int stepNumber, 
								Vector passed_u, 
								double passed_gFunctionValue, 
								Vector passedGradientInStandardNormalSpace);
	Vector getSearchDirection();

	int	check(Vector u_old, 
			  double g_old, 
			  Vector grad_G_old, 
			  double stepSize,
			  Vector stepDirection,
			  double g_new);
	double getMeritFunctionValue(Vector u, double g, Vector grad_G);
	int updateMeritParameters(Vector u, double g, Vector grad_G);

	int setAlpha(double alpha);

protected:

private:
	Vector searchDirection;
	double thetaFunction;
	double alpha;
	double gamma;
	double delta;
};

#endif
