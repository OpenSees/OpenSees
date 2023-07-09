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
**   Optimization module developed by:                                **
**   Quan Gu  (qgu@ucsd.edu)                                          **
**   Joel Conte (jpconte@ucsd.edu)                                    **
**   Philip Gill (pgill@ucsd.edu)                                     **
** ****************************************************************** */


//
// Written by Quan Gu (qgu@ucsd.edu)    March 2010
//


#ifndef ObjectiveFunction_h
#define ObjectiveFunction_h

#include <OptimizationDomainComponent.h>
#include <OptimizationDomain.h>
#include <Vector.h>
#include <tcl.h>

//#include <fstream>

class ObjectiveFunction : public OptimizationDomainComponent
{

public:
	ObjectiveFunction(				   int passedTag, 
									   OptimizationDomain * passedOptimizationDomain, 
									   Tcl_Interp *passedTclInterp, 
									   bool passedIsGradProvided,
									   Vector * passedLinearAddition,
									   char * passedTclFileName,
									   char * passedName, 
									   char * passedGradientName, 
									   double passedLowerBound, 
									   double passedUpperBound, 
									   double passedMultiplier, 									   
									   double passedState		   
									   );
	~ObjectiveFunction();

	double getValue();
	double getLowerBound();
	double getUpperBound();
	double getMultiplier();
	double getState();
	Vector * getGradientPtr();
	Vector * getLinearAdditionPtr();

	int update();
	bool isGradientProvided(){return isGradProvided;};
	int getNumberOfEvaluations() {return numOfComputation;};

	void Print(OPS_Stream &s, int flag){};

	int sendSelf(int,Channel &){ return 0;}
	int recvSelf(int,Channel &,FEM_ObjectBroker &) {return 0;}

private:
	
	double value;
	double lowerBound;
	double upperBound;
	double multiplier;
	double state;

	char name[25];
	char gradientName[25];

	char TclFileName[50];

	Vector * gradient;
 
	Tcl_Interp * theTclInterp;
	
	bool isGradProvided; // true: yes, false: not needed. 
	int numOfComputation;
	Vector * linearAddition;
};


#endif
