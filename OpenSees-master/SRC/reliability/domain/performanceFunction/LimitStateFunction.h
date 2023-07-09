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
                                                                        
// $Revision: 1.15 $
// $Date: 2010-06-10 20:14:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunction.h,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#ifndef LimitStateFunction_h
#define LimitStateFunction_h

#include <PerformanceFunction.h>
//#include <FunctionEvaluator.h>
#include <map>
#include <string>
using namespace std;

class LimitStateFunction : public PerformanceFunction
{

public:
	LimitStateFunction(int tag, const char *expression);
	~LimitStateFunction();

	// Methods to access LSF
	void Print(OPS_Stream &s, int flag = 0);
	const char *getExpression(void);
	
	// Methods to add/remove/get gradient of LSF
	int addGradientExpression(const char *expression, int rvTag);
	int removeGradientExpression(int rvTag);
	const char* getGradientExpression(int rvTag);

protected:

private:
	char *theExpression;

	// STL map for analytic gradients of this LSF
	map<int, string> mapOfGradientExpressions;

};

#endif
