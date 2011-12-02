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
// $Date: 2003-10-27 23:04:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/modulatingFunction/TrapezoidalModulatingFunction.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef TrapezoidalModulatingFunction_h
#define TrapezoidalModulatingFunction_h

#include <ModulatingFunction.h>
#include <Filter.h>

class TrapezoidalModulatingFunction : public ModulatingFunction
{

public:
	TrapezoidalModulatingFunction(int tag,
		                          Filter *theFilter,
								  double t1, 
								  double t2, 
								  double t3, 
								  double t4,
								  double amplitude);
	~TrapezoidalModulatingFunction();

	void Print(OPS_Stream &s, int flag =0);

	double getAmplitude(double time);
	double getMaxAmplitude();
	Filter *getFilter();

protected:

private:
	Filter *theFilter;
	double t1;
	double t2;
	double t3;
	double t4;
	double amplitude;
};

#endif
