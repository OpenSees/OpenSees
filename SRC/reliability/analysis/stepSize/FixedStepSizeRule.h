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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-02-29 19:47:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/stepSize/FixedStepSizeRule.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef FixedStepSizeRule_h
#define FixedStepSizeRule_h

#include <StepSizeRule.h>

class FixedStepSizeRule : public StepSizeRule
{

public:
	FixedStepSizeRule(double stepSize);
	~FixedStepSizeRule();
	/////S added by K Fujimura /////
	/*int		computeStepSize(Vector u, Vector grad_G, double G, Vector d, int stepNumber,
		int reschk=0);*/
	int		computeStepSize(const Vector &u, const Vector &grad_G, double G, const Vector &d, int stepNumber, int reschk=0);
	/////E added by K Fujimura /////
	double	getStepSize();
	double	getInitialStepSize();
	double getGFunValue();
	int getNumReductions();	///// added by K Fujimura /////


protected:

private:
	double stepSize;
	double gFunValue;
	int numReducntion;	///// added by K Fujimura /////

};

#endif
