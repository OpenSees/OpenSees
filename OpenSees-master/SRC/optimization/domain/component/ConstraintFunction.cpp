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


#include <ConstraintFunction.h>
#include <OptimizationDomainComponent.h>
#include <iostream>
#include <string.h>
#include <classTags.h> 
#include <OPS_Globals.h>

# define CONSTRAINT_FUNCTION 10045650  //should move to classTags.h

ConstraintFunction::ConstraintFunction(int passedTag, int passedNumberOfConstraint,
									   OptimizationDomain * passedOptimizationDomain, 
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
									   )
:OptimizationDomainComponent(passedTag, CONSTRAINT_FUNCTION)
{
	theOptimizationDomain=passedOptimizationDomain;
	theTclInterp = passedTclInterp;


	numOfComputation=0;

//	int numOfDV=theDomain->getNumberOfDesignVariables();

	numOfConstraintFunction = passedNumberOfConstraint;
	isGradProvided=passedIsGradProvided;

	strcpy(TclFileName,passedTclFileName);
	strcpy(name,passedName);

	if (passedGradientName !=0){
		gradientName = new char[25];
		strcpy(gradientName,passedGradientName);
	}
	else gradientName=0;

	
	
	int numOfDesignVariables = theOptimizationDomain->getNumberOfDesignVariables();

	value = new Vector(numOfConstraintFunction);
	(*value).Zero();

	lowerBound = new Vector(numOfConstraintFunction);
	if (passedLowerBound !=0) (*lowerBound).addVector(0,(*passedLowerBound),1.0);
	else { 
		for(int i=0;i<numOfConstraintFunction; i++)
			(*lowerBound)(i)=-1e20;
	}//else 

	upperBound = new Vector(numOfConstraintFunction);
	if (passedUpperBound !=0) (*upperBound).addVector(0,(*passedUpperBound),1.0);
	else { 
		for(int i=0;i<numOfConstraintFunction; i++)
			(*upperBound)(i)= 1e20;
	}//else 

	
	if (passedMultiplier !=0) {
		multiplier = new Vector(numOfConstraintFunction);
		(*multiplier).addVector(0,(*passedMultiplier),1.0);
	}
	else multiplier = 0;


	if(passedState !=0){
		state = new Vector(numOfConstraintFunction);
		(*state).addVector(0,(*passedState),1.0);
	}
	else state=0;


    if (isGradProvided) 
		gradient= new Matrix(numOfConstraintFunction,numOfDesignVariables);
	else gradient=0;


	
	if (passedLinearAdd !=0){
		linearAddition=new Matrix(numOfConstraintFunction,numOfDesignVariables);
		(*linearAddition).addMatrix(0,(*passedLinearAdd),1.0);
	}
	else 
		linearAddition=0;
	

//	update();


}

ConstraintFunction::~ConstraintFunction()
{
if (value !=0) delete  value;
if (lowerBound !=0) delete  lowerBound;
if (upperBound !=0) delete  upperBound;
if (multiplier !=0) delete  multiplier;
if (state !=0) delete  state;
if (gradient !=0) delete  gradient;
if (linearAddition !=0) delete  linearAddition;
}


	

int ConstraintFunction::getNumOfConstraintFunctions(){return numOfConstraintFunction;};


Vector * ConstraintFunction::getValuePtr() { return value;};

Vector * ConstraintFunction::getLowerBoundPtr(){ return lowerBound;};

Vector * ConstraintFunction::getUpperBoundPtr(){ return upperBound;};

Vector * ConstraintFunction::getMultiplierPtr(){ return multiplier;};

Vector * ConstraintFunction::getStatePtr(){ return state;};

Matrix * ConstraintFunction::getGradientPtr(){
	return gradient;
	};

Matrix * ConstraintFunction::getLinearAdditionPtr(){ 
	return linearAddition;
	};

int ConstraintFunction::update(){
	/*  1, check and run userinput tcl file 
		2, to get value Vector
		4, to get Graident Matrix (if isGradientPrivided)
    */ 

	if(Tcl_EvalFile(theTclInterp, TclFileName) !=TCL_OK){
		opserr<<"ConstraintFunction: the file "<<TclFileName<<" can not be run!"<<endln;
		exit(-1);
	}
	
	
	
	const char *  theValue;	
	char index[10];

	for(int i=0; i<numOfConstraintFunction; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(theTclInterp, name,index,TCL_GLOBAL_ONLY );
				(*value)(i) = atof(theValue);
//	opserr<<"ConstraintFunction::update theValue is"<<(*value)(i)<<endln;
	};


  	
	if (isGradientProvided()){
	
		char temp[5];
		int numOfDVs = theOptimizationDomain->getNumberOfDesignVariables();
		for(int i=0;i<numOfConstraintFunction; i++){
			for(int j=0; j<numOfDVs; j++){
			       sprintf(temp,"%d",i+1);   // begin with 1
				   strcpy(index,temp);
				   sprintf(temp,"%d",j+1);   // begin with 1
				   strcat(index,",");
				   strcat(index,temp);
					
				   if(Tcl_GetVar2(theTclInterp, gradientName,index,TCL_GLOBAL_ONLY ) ==NULL){
					   opserr<<"constraintFunction: theGradient can not be read index G("<<index<<")"<<endln;
		//             exit(-1);				 
					   }

					 theValue = Tcl_GetVar2(theTclInterp, gradientName,index,TCL_GLOBAL_ONLY );
				   (*gradient)(i,j) = atof(theValue);
			}; //for
		} //for
	}

	



	numOfComputation++;
	 return 0;
	};
