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
// $Date: 2008-02-29 19:47:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/meritFunction/CriteriaReductionMeritFunctionCheck.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef CriteriaReductionMeritFunctionCheck_h
#define CriteriaReductionMeritFunctionCheck_h

#include <MeritFunctionCheck.h>
#include <ReliabilityConvergenceCheck.h>

class CriteriaReductionMeritFunctionCheck : public MeritFunctionCheck
{

public:
	CriteriaReductionMeritFunctionCheck(ReliabilityConvergenceCheck *theReliabilityConvergenceCheck);
	~CriteriaReductionMeritFunctionCheck();

	/*int check(const Vector &u_old, 
		  double g_old, 
		  const Vector &grad_G_old, 
		  double stepSize,
		  const Vector &stepDirection,
		  double g_new, 
		  const Vector &grad_G_new,);
	double getMeritFunctionValue(const Vector &u, double g,
				     const Vector &grad_G);
	int updateMeritParameters(const Vector &u, double g,
				  const Vector &grad_G);*/

 /////S added by K Fujimura /////
	int	check(const Vector &u_old, 
			  double g_old, 
			  const Vector &grad_G_old, 
			  double stepSize,
			  const Vector &stepDirection,
			  double g_new, 
			  const Vector &grad_G_new,
			  int reschk=0);
	double getMeritFunctionValue(const Vector &u, double g, const Vector &grad_G);
	int updateMeritParameters(const Vector &u, double g, const Vector &grad_G,
		int reschk=0);
/////E added by K Fujimura /////

protected:

private:
	ReliabilityConvergenceCheck *theReliabilityConvergenceCheck;

};

#endif
