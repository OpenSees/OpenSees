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
// $Date: 2008-02-29 19:47:19 $
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
	StandardReliabilityConvergenceCheck(double e1, double e2, double scaleValue, int printFlag);
	~StandardReliabilityConvergenceCheck();

	int check(const Vector &u, double g, const Vector &gradG);
	int getNumberOfCriteria();
	double getCriteriaValue(int whichCriteria);
	int setScaleValue(double scaleValue);
/////////////////////////////////////
//S  Modified by K Fujimura 10/10/2004
/////////////////////////////////////
	double getScaleValue(){ return scaleValue; }
	void Scalefix(bool);
	bool getScfix(){ return fixscale; }
	double getCheck1();
	double getCheck2();
	int	checkG(double g);
/////////////////////////////////////
//E  Modified by K Fujimura 10/10/2004
/////////////////////////////////////
protected:

private:  //these variables are not in K.F.
	double e1, e2;
	double criterium1, criterium2;
	double scaleValue;
	int printFlag;
	ofstream logfile;
	bool fixscale;
};

#endif
