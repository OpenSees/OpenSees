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
// $Date: 2008-02-29 19:47:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/convergenceCheck/ReliabilityConvergenceCheck.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef ReliabilityConvergenceCheck_h
#define ReliabilityConvergenceCheck_h

#include <Vector.h>

class ReliabilityConvergenceCheck
{

public:
	ReliabilityConvergenceCheck();
	virtual ~ReliabilityConvergenceCheck();

	virtual int check(const Vector &u, double g, const Vector &gradG) = 0;
	virtual int getNumberOfCriteria() = 0;
	virtual double getCriteriaValue(int whichCriteria) = 0;
	virtual int setScaleValue(double scaleValue) = 0;

	//S Modified by K Fujimura 10/10/2004
	virtual double getScaleValue() = 0;
	virtual void Scalefix(bool) = 0;
	virtual bool getScfix()=0;
	virtual double getCheck1()=0;
	virtual double getCheck2()=0;
	virtual int	checkG(double g) = 0;
	//E Modified by K Fujimura 10/10/2004
protected:

private:

};

#endif
