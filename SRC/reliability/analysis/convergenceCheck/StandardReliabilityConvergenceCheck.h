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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-03-04 00:38:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/convergenceCheck/StandardReliabilityConvergenceCheck.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef StandardReliabilityConvergenceCheck_h
#define StandardReliabilityConvergenceCheck_h

#include <ReliabilityConvergenceCheck.h>

#include <fstream>
using std::ofstream;

class StandardReliabilityConvergenceCheck : public ReliabilityConvergenceCheck
{

public:
	StandardReliabilityConvergenceCheck(double e1, double e2, int printFlag);
	~StandardReliabilityConvergenceCheck();

	int	check(Vector u, double g, Vector gradG);
	int getNumberOfCriteria();
	double getCriteriaValue(int whichCriteria);
	int setScaleValue(double scaleValue);

protected:

private:
	double e1, e2;
	double criterium1, criterium2;
	double scaleValue;
	int printFlag;
};

#endif
