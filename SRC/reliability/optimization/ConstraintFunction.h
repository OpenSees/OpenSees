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
**   Philip Gill (pgill@ucsd.edu)                                     **                                 **
** ****************************************************************** */
                                                                        

//
// Written by Quan Gu (qgu@ucsd.edu)
//

#ifndef ConstraintFunction_h
#define ConstraintFunction_h

#include <ReliabilityDomainComponent.h>
#include <ReliabilityDomain.h>
#include <Vector.h>
#include <Matrix.h>
#include <tcl.h>


class ConstraintFunction : public ReliabilityDomainComponent
{

public:
	ConstraintFunction(int passedTag, int passedNumberOfConstraint,
									   ReliabilityDomain * passedReliabilityDomain, 
									   Tcl_Interp *passedTclInterp, 
									   bool passedIsGradProvided, 
									   Matrix * passedLinearAdd,
									   char * passedTclFileName,
									   char * passedName, 
									   char * passedGradientName, 
									   
									   Vector * passedLowerBound, 
									   Vector * passedUpperBound, 
									   Vector * passedMultiplier, 									   
									   Vector * passedState
									   );

	~ConstraintFunction();
	int getNumOfConstraintFunctions();
	Vector * getValuePtr();
	Vector * getLowerBoundPtr();
	Vector * getUpperBoundPtr();
	Vector * getMultiplierPtr();
	Vector * getStatePtr();

	Matrix * getGradientPtr();
	Matrix * getLinearAdditionPtr();

//	int computeValueAndGradient();
	int update();
	
	bool isGradientProvided(){return isGradProvided;};
	void Print(OPS_Stream &s, int flag){};

private:
	int numOfConstraintFunction;
	Vector * value;
	Vector * lowerBound;
	Vector * upperBound;
	Vector * multiplier;
	Vector * state;

	char name[25];
	char * gradientName ;


	char TclFileName[50];

	Matrix * gradient ;

	ReliabilityDomain *theReliabilityDomain;
	Tcl_Interp * theTclInterp;
	
	bool isGradProvided; // true: yes, false: not needed. 
	int numOfComputation;
	Matrix * linearAddition ;  //A

};


#endif
