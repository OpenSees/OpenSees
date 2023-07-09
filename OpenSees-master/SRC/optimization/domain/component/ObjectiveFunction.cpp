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


#include <ObjectiveFunction.h>
#include <OptimizationDomainComponent.h>
#include <iostream>
#include <string.h>
#include <classTags.h> 
#include <OPS_Globals.h>

# define OBJECTIVE_FUNCTION 1001  //should move to classTags.h

ObjectiveFunction::ObjectiveFunction(int passedTag, 
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
									   )
:OptimizationDomainComponent(passedTag, OBJECTIVE_FUNCTION)
{
//	gradient =0;
//	linearAddition=0;
	value=0;
	numOfComputation=0;	

	theOptimizationDomain=passedOptimizationDomain;
	theTclInterp = passedTclInterp;
	
	isGradProvided=passedIsGradProvided;

	strcpy(TclFileName,passedTclFileName);
	strcpy(name,passedName);


	if (passedGradientName !=0) strcpy(gradientName,passedGradientName);
	
	lowerBound = passedLowerBound;
	upperBound = passedUpperBound;

	multiplier = passedMultiplier;
	state = passedState;
	
	int numOfDesignVariables = theOptimizationDomain->getNumberOfDesignVariables();
	
	if (isGradProvided)
		gradient= new Vector(numOfDesignVariables);
	else gradient=0;

	if (passedLinearAddition !=0)  {
		linearAddition=new Vector(numOfDesignVariables);
		(*linearAddition).addVector(0.0,*passedLinearAddition,1.0);
	}
	else linearAddition=0;



//	update();





};

ObjectiveFunction::~ObjectiveFunction()
{
	if (gradient !=0) delete [] gradient;
	if (linearAddition !=0) delete [] linearAddition;
};


	


double ObjectiveFunction::getValue() { return value;};

double ObjectiveFunction::getLowerBound(){ return lowerBound;};

double ObjectiveFunction::getUpperBound(){ return upperBound;};

double ObjectiveFunction::getMultiplier(){ return multiplier;};

double ObjectiveFunction::getState(){ return state;};

Vector * ObjectiveFunction::getGradientPtr(){
	return gradient;
	};

Vector * ObjectiveFunction::getLinearAdditionPtr(){ 
	return linearAddition;
	};


int ObjectiveFunction::update(){
	/* check and 
	    1. run userinput tcl file 
		2, to get value Vector
		4, to get Graident Vector (if isGradientPrivided)
	*/

	if(Tcl_EvalFile(theTclInterp, TclFileName) !=TCL_OK){
		opserr<<"the file "<<TclFileName<<" can not be run!"<<endln;
		exit(-1);
	}
		

/*	if (Tcl_GetDouble(theTclInterp, name, &value) != TCL_OK) {
			opserr << "ERROR: invalid input: value \n";
			return TCL_ERROR;
	}
*/
	const char * myStr;
	myStr = Tcl_GetVar(theTclInterp, name,TCL_GLOBAL_ONLY );
	value = atof(myStr);	
	
//	opserr<<"ObjectiveFunction::update theValue is"<<value<<endln;

	if(isGradientProvided()){
		int numOfDvs = theOptimizationDomain->getNumberOfDesignVariables();
		char index[5];
		const char * myValue;
		for(int i=0; i<numOfDvs;i++){
			sprintf(index,"%d",i+1);  // index from 1
			 myValue =  Tcl_GetVar2(theTclInterp, gradientName,index,TCL_GLOBAL_ONLY);
			 (*gradient)(i) = atof(myValue);
			
			
		}//for
	}//if


	numOfComputation++;
	 return 0;
	};

