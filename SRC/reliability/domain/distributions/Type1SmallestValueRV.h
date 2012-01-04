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
// $Date: 2007-02-17 21:27:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type1SmallestValueRV.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef Type1SmallestValueRV_h
#define Type1SmallestValueRV_h

#include <RandomVariable.h>

class Type1SmallestValueRV : public RandomVariable
{

public:
	Type1SmallestValueRV(int tag, double mean, double stdv);
	Type1SmallestValueRV(int tag, const Vector &parameters);
	~Type1SmallestValueRV();
	
	// pure virtual defining variable type and properties
	const char* getType();
	double getMean();
	double getStdv();
	const Vector &getParameters();
	int setParameters(double mean, double stdv);
	
	// RV functionality
	double getPDFvalue(double rvValue);
	double getCDFvalue(double rvValue);
	double getInverseCDFvalue(double rvValue); 
	
	// starting point methods
	int setStartValue(double newVal) {startValue = newVal; return 0;}
	double getStartValue() {return startValue;}
	
	// other
	void Print(OPS_Stream &s, int flag = 0);


protected:

private:
	double startValue;
	double u;
	double alpha;

};

#endif
