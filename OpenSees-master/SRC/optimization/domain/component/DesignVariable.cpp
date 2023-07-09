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

 

#include "DesignVariable.h"
#include "DesignVariablePositioner.h"

#define DESIGN_VARIABLE 10007567

#include <math.h>

DesignVariable::DesignVariable(int passedTag, 
			 char *passedName,
			 double passedValue,
			 double passedUpperBound,
			 double passedLowerBound,
			 Tcl_Interp *passedTclInterp, 
			 OptimizationDomain *passedOptimizationDomain,
			 
			 double passedXMultiplier

			 			 )
:OptimizationDomainComponent(passedTag,DESIGN_VARIABLE)

{
	value = passedValue;
	if (value !=0)
		scale = fabs(value);  // 2006  Dec..
	else scale =1.0;

	//--- note: all scaling things are done in the upper level, not in this class 
	//(i.e., done by ones who handle the designVariable, like snoptanalysis class. )

	oldValue = value;
    strcpy(name,passedName);
	upperBound = passedUpperBound;
	lowerBound = passedLowerBound;

	xMultiplier = passedXMultiplier; // 0

	theTclInterp = passedTclInterp;
	theOptimizationDomain =  passedOptimizationDomain;
	
	numOfMyDVPositioners = 0;

	maxNumOfMyDVPositioners =256;   // could enlarge if needed.
	myDVPositionerTags = new int [maxNumOfMyDVPositioners]; 




};

DesignVariable::~DesignVariable()
{

	if (myDVPositionerTags ==0) delete[] myDVPositionerTags;


};



double DesignVariable::getUpperBound(){return upperBound;};

double DesignVariable::getLowerBound(){return lowerBound;};


char * DesignVariable::getName(){return name;};

double DesignVariable::getValue(){return value;};

double DesignVariable::getOldValue(){return oldValue;};

double DesignVariable::getXMultiplier(){return xMultiplier;};


int DesignVariable::setXMultiplier(double newXMultiplier){
	xMultiplier=newXMultiplier;
	return 0;
};

int DesignVariable::setValue(double newValue){
	value = newValue;
	return 0;
};

int DesignVariable::setOldValue(double newValue){
	oldValue=newValue;
	return 0;

};

int DesignVariable::setUpperBound(double newBound){

	upperBound=newBound;
	return 0;

};

int DesignVariable::setLowerBound(double newBound){

	lowerBound=newBound;
	return 0;

};

int DesignVariable::setName(char * newName){

	strcpy(name,newName);
	return 0;

}


int DesignVariable::addDVPositioner(int tag){

	if (numOfMyDVPositioners>=maxNumOfMyDVPositioners){  //need enlarge

		int * theNewDVPositionerTags = new int [maxNumOfMyDVPositioners*2];
			
		for (int i=0; i<maxNumOfMyDVPositioners;i++)
				theNewDVPositionerTags[i] = myDVPositionerTags[i];
			
		int * temp =myDVPositionerTags;

		myDVPositionerTags = theNewDVPositionerTags; 
			
		delete [] temp;
		opserr <<"enlarge maxNumOfMyDVPositioners by " <<maxNumOfMyDVPositioners*2<<endln;
		maxNumOfMyDVPositioners =maxNumOfMyDVPositioners*2;
	}

	
	this->myDVPositionerTags[numOfMyDVPositioners] = tag;		
	numOfMyDVPositioners++;
		
	
	return 0;

};

int DesignVariable::getDVPositionerTag(int i){ // i begins from 0;
	return myDVPositionerTags[i];
}



int DesignVariable::update(double newX){

//  1.--- in tcl, set value

	  char tclAssignment[50];
	  strcpy (tclAssignment, "set ");
	  strcat (tclAssignment, name);
	  strcat (tclAssignment, " ");
      char theValueString[40];

	  sprintf(theValueString," %20.14e", newX);
	  strcat (tclAssignment, theValueString);

   
	  if (Tcl_GetVar(theTclInterp, name, TCL_GLOBAL_ONLY) ==NULL)
	   {
			opserr<<"Fatal::the variable with name: "<<name <<" is missed!"<<endln;   
			exit(-1);
	   }

	   if (Tcl_Eval(theTclInterp,tclAssignment) !=TCL_OK ){
   			opserr<<"Fatal,DesignVariable::update, can not set varuable with name: "<<name <<"in tcl command!"<<endln;  
			opserr<<"command is:"<<tclAssignment<<endln;
			exit(-1);
	   }

// 2.--- change value in this class
       this->setValue(newX);

// 3.--- update Parameter in structure system if parameter is in structure system .


	   if (numOfMyDVPositioners > 0)  { // this means structure is involved
														
	//		char theRevertToStartCommand[10] = "reset";
	//		Tcl_Eval( theTclInterp, theRevertToStartCommand );

			DesignVariablePositioner * theDesignVariablePositioner;
			int dvTag;
			for (int i=0 ; i< numOfMyDVPositioners ; i++ )  {
				dvTag = myDVPositionerTags[i];
				theDesignVariablePositioner = theOptimizationDomain->getDesignVariablePositionerPtr(dvTag);
    			theDesignVariablePositioner->update(newX);
			}

	   } //if
  

	   return 0;
}
