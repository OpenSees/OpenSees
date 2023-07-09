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



#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>
#include <Domain.h>
#include <Parameter.h>
#include <TclOptimizationBuilder.h>
#include <OptimizationDomain.h>

 
///////////#include <SnoptProblem.h>  

#if _HAVESNOPT
#include <SnoptAnalysis.h>  
#endif

#include <DesignVariable.h>
#include <DesignVariablePositioner.h>
#include <ConstraintFunction.h>
#include <ObjectiveFunction.h>

#ifdef _HAVESNOPT
SNOPTAnalysis * theSNOPTAnalysis=0;
#endif

static Domain *theStructuralDomain = 0;
OptimizationDomain *theOptimizationDomain = 0;

 
int TclOptimizationModelBuilder_addDesignVariable(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv);
int TclOptimizationModelBuilder_addDesignVariablePositioner(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclOptimizationModelBuilder_addObjectiveFunction(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv);
int TclOptimizationModelBuilder_addConstraintFunction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclOptimizationModelBuilder_runSNOPTAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);


 
// constructor: the constructor will add certain commands to the interpreter

TclOptimizationBuilder::TclOptimizationBuilder( Domain &passedSructuralDomain, Tcl_Interp *interp)
{

   theInterp = interp;
   theStructuralDomain	= &passedSructuralDomain;            
   theOptimizationDomain	= new OptimizationDomain();


  // call Tcl_CreateCommand for class specific commands
   //Tcl_CreateCommand(interp, "wipe",                  &wipeModel,                           (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
   Tcl_CreateCommand(interp, "designVariable",(Tcl_CmdProc *)TclOptimizationModelBuilder_addDesignVariable,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "designVariablePositioner",(Tcl_CmdProc *)TclOptimizationModelBuilder_addDesignVariablePositioner,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "constraintFunction", (Tcl_CmdProc *)TclOptimizationModelBuilder_addConstraintFunction,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "objectiveFunction", (Tcl_CmdProc *)TclOptimizationModelBuilder_addObjectiveFunction,(ClientData)NULL, NULL);
   Tcl_CreateCommand(interp, "runSNOPTAnalysis", (Tcl_CmdProc *)TclOptimizationModelBuilder_runSNOPTAnalysis,(ClientData)NULL, NULL); 
 
 
  



} 

TclOptimizationBuilder::~TclOptimizationBuilder()
{


#ifdef _HAVESNOPT
  if (theSNOPTAnalysis != 0)
    delete theSNOPTAnalysis;
 ////// if (theSNOPT != 0)
 //////   delete theSNOPT;

  theSNOPTAnalysis = 0;
 // theSNOPT = 0;

#endif

  Tcl_DeleteCommand(theInterp, "designVariable");
  Tcl_DeleteCommand(theInterp, "designVariablePositioner");
  Tcl_DeleteCommand(theInterp, "constraintFunction");
  Tcl_DeleteCommand(theInterp, "objectiveFunction");
  Tcl_DeleteCommand(theInterp, "runSNOPTAnalysis");
  


}

 

OptimizationDomain *
TclOptimizationBuilder::getOptimizationDomain()
{
	return theOptimizationDomain;
}
 



// command: designVariable  1 -name E  < -startPt $E   -lowerBound [expr $E*0.8] -upperBound [expr $E*1.2] >
  
int 
TclOptimizationModelBuilder_addDesignVariable(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv)
{


	
 DesignVariable *theDesignVariable = 0;
  int tag;
  char name[20]="";
  char valueString[20]="";
  double value=0;
  double lowerBound=0;
  double upperBound=0;
  int numberOfArguments = argc;

  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 3) {
		opserr << "ERROR: invalid number of arguments to designVariable command \n";
		return TCL_ERROR;
  }


  // GET TAG NUMBER
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }


	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-name") == 0)||(strcmp(argv[argvCounter],"-Name") == 0)) {
			argvCounter++;
			strcpy(name,argv[argvCounter]);
			argvCounter++;
		}// if

		else if (strcmp(argv[argvCounter],"-startPt") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &value) != TCL_OK) {
			opserr << "ERROR: invalid input: startPt \n";
			return TCL_ERROR;
			}
			strcpy(valueString,argv[argvCounter]);
			lowerBound=value;
			upperBound=value;
				argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-lowerBound") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &lowerBound) != TCL_OK) {
			opserr << "ERROR: invalid input: startPt \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-upperBound") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &upperBound) != TCL_OK) {
			opserr << "ERROR: invalid input: startPt \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if
		
		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}

	};  // while


// in tcl, run 'set name value'
  char tclAssignment[50];
  strcpy (tclAssignment, "set ");
  strcat (tclAssignment, name);
  strcat (tclAssignment, " ");

   
 // _gcvt( value, 7, buffer );
   strcat (tclAssignment, valueString);
  
   
   if (Tcl_GetVar(interp, name, TCL_GLOBAL_ONLY) !=NULL)
   {
		opserr<<"Fatal::the variable with name: "<<name <<" is already in system, please use another name!"<<endln;   
		exit(-1);
   }

   if (Tcl_Eval(interp,tclAssignment) !=TCL_OK ){
   		opserr<<"Fatal::can not set varuable with name: "<<name <<"in tcl command!"<<endln;   
		exit(-1);
   }

// here tag ;
  theDesignVariable = new DesignVariable(tag, 
			 name,
			 value,
			 upperBound,
			 lowerBound,
			 interp,
			 theOptimizationDomain,
			 0);

  if (theDesignVariable == 0) {
		opserr << "ERROR: could not create random theDesignVariable number " << tag << endln;
		return TCL_ERROR;
	  }
  

  // ADD THE OBJECT TO THE DOMAIN
  if (theOptimizationDomain->addDesignVariable(theDesignVariable) == false) {
	opserr << "ERROR: failed to add theDesignVariable variable to the domain (wrong number of arguments?)\n";
	opserr << "theDesignVariable variable: " << tag << endln;
	delete theDesignVariable; // otherwise memory leak
	return TCL_ERROR;
  }
 
  return TCL_OK;
}
					   
 


// command designVariablePositioner 1   -dvNum 1 -element 1     -material E  < -parameter 1>
int 
TclOptimizationModelBuilder_addDesignVariablePositioner(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	
	DesignVariablePositioner *theDesignVariablePositioner = 0;
	int tag;
	int dvNumber;
	int tagOfObject;
	DomainComponent *theObject;
	int argvCounter = 1;


	// READ THE TAG NUMBER
	if (Tcl_GetInt(interp, argv[argvCounter++], &tag) != TCL_OK) {
		opserr << "ERROR: Invalid input tag to design variable positioner." << endln;
		return TCL_ERROR;
	}


	if (strcmp(argv[argvCounter],"-dvNum") == 0) {
		argvCounter++;
		
		// READ THE RANDOM VARIABLE NUMBER
		if (Tcl_GetInt(interp, argv[argvCounter++], &dvNumber) != TCL_OK) {
			opserr << "ERROR: invalid input: dvNumber \n";
			return TCL_ERROR;
		}

		// CHECK THAT THE RANDOM VARIABLE ACTUALLY EXISTS
		DesignVariable *theDesignVariable = 0;
		theDesignVariable = theOptimizationDomain->getDesignVariablePtr(dvNumber);
		if (theDesignVariable == 0){
			opserr << "ERROR:: A non-existing design variable number " << dvNumber << " is being positioned in the model " << endln;
			return TCL_ERROR;
		}
	}
	else {
		opserr << "ERROR: Illegal design variable specification in  " << endln
			<< " design variable positioner command. " << endln;
		return TCL_ERROR;
	}
	

	const char **data = new const char *[argc-argvCounter-2];
	int ii,jj;
	for (ii=argvCounter+2, jj=0; ii<argc; ii++, jj++)
	  data[jj] = argv[ii];

	// IF UNCERTAIN *ELEMENT* PROPERTY
	if (strcmp(argv[argvCounter],"-element") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			argvCounter++;
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}

		theObject = (DomainComponent *)theStructuralDomain->getElement(tagOfObject);

		theDesignVariablePositioner = new DesignVariablePositioner(tag,
									   theOptimizationDomain,
									   dvNumber,
									   theObject,
									   data,
									   argc-argvCounter);


	//	int dvnumber = theDesignVariablePositioner->getDVNumber();
	}

	// IF UNCERTAIN *LOAD*
	else if (strcmp(argv[argvCounter],"-loadPattern") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getLoadPattern(tagOfObject);


		theDesignVariablePositioner = new DesignVariablePositioner(tag,
									   theOptimizationDomain,
									   dvNumber,
									   theObject,
									   data,
									   argc-argvCounter);
	}

	// IF UNCERTAIN *NODE* PROPERTY
	else if (strcmp(argv[argvCounter],"-node") == 0) {
		argvCounter++;

		if (Tcl_GetInt(interp, argv[argvCounter++], &tagOfObject) != TCL_OK) {
			opserr << "ERROR: invalid input: tagOfObject \n";
			return TCL_ERROR;
		}
		theObject = (DomainComponent *)theStructuralDomain->getNode(tagOfObject);

		theDesignVariablePositioner = new DesignVariablePositioner(tag,
							           theOptimizationDomain,
									   dvNumber,
									   theObject,
									   data,
									   argc-argvCounter);
	}
	else if (strcmp(argv[argvCounter],"-parameter") == 0) {
	  argvCounter++;
	  int paramTag;
	  if (Tcl_GetInt(interp, argv[argvCounter++], &paramTag) != TCL_OK) {
	    opserr << "ERROR: invalid input in positioner: parameter tag \n";
	    return TCL_ERROR;
	  }

	  Parameter *theParameter = theStructuralDomain->getParameter(paramTag);

	  if (theParameter == 0) {
	    opserr << "ERROR: parameter with tag " << paramTag 
		   << " not found in structural domain\n";
	    return TCL_ERROR;
	  }
	  else {
	    theDesignVariablePositioner =
	      new DesignVariablePositioner(tag, theOptimizationDomain, dvNumber, *theParameter);
	  }
	}
	else {
		opserr << "ERROR: Unknown parameter in designVariablePositioner" << endln;
		return TCL_ERROR;
	}

	delete [] data;

	// ADD THE RANDOMVARIABLEPOSITIONER TO THE DOMAIN
	if (theOptimizationDomain->addDesignVariablePositioner(theDesignVariablePositioner) == false) {
		opserr << "ERROR: failed to add random variable positioner number " << tag << " to the domain." << endln;
		delete theDesignVariablePositioner; // otherwise memory leak
		return TCL_ERROR;
	}

	return TCL_OK;

}




// command: objectiveFunction 1  -name F  -GradientName G -tclFile objective.tcl  ; # -lowerBound -1.e20 -upperBound 1.e20  - multiplier $E -state 0 -linearAdd A ;
int 
TclOptimizationModelBuilder_addObjectiveFunction(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv)
{
 
  ObjectiveFunction *theObjectiveFunction = 0;
  int tag;
  char name[25] = "" ;
  char * gradientName = 0;
  char tclFileName[35] = "";
  double lowerBound= -1.e20;;
  double upperBound=1.e20;
  double multiplier = 0; 
  double state = 0; 
  int numberOfArguments = argc;

  Vector *linearAdd=0;

  bool isGradProvided=false;

  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 7) {
		opserr << "ERROR: invalid number of arguments to designVariable command "<<endln;
	    opserr <<"command: objectiveFunction 1  -name F  -GradientName G  -tclFile objective.tcl  ; # -lowerBound -1.e20 -upperBound 1.e20 - multiplier $E -state 0 -linearAdd A"<<endln;
		return TCL_ERROR;
  }


  // GET TAG NUMBER
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }


	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-name") == 0)||(strcmp(argv[argvCounter],"-Name") == 0)) {

			argvCounter++;
			strcpy(name,argv[argvCounter]);

			if ((Tcl_GetVar(interp, name, TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, name,"1",TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, name,"1,1",TCL_GLOBAL_ONLY ) != NULL)){
				opserr<<"Fatal: objectiveFunction Name "<< name << " is already been used, please change another name"<<endln; 
				exit(-1);
			}
				
			
			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-GradientName") == 0)||(strcmp(argv[argvCounter],"-gradientName") == 0)) {
			
			argvCounter++;
			
			gradientName = new char[30];

			strcpy(gradientName,argv[argvCounter]);
			
			if ((Tcl_GetVar(interp, gradientName, TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, gradientName,"1",TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, gradientName,"1,1",TCL_GLOBAL_ONLY ) != NULL)){
				opserr<<"Fatal: objectiveFunction gradient Name "<< gradientName << " is already been used, please change another name"<<endln; 
				exit(-1);
			}
				
			isGradProvided = true;
			argvCounter++;
		}// else if
		
		else if ((strcmp(argv[argvCounter],"-tclFile") == 0)||(strcmp(argv[argvCounter],"-TclFile") == 0)) {
			argvCounter++;
			strcpy(tclFileName,argv[argvCounter]);
			argvCounter++;
		}// else if
		
		else if (strcmp(argv[argvCounter],"-lowerBound") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &lowerBound) != TCL_OK) {
			opserr << "ERROR: invalid input: lowerBound \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-upperBound") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &upperBound) != TCL_OK) {
			opserr << "ERROR: invalid input: upperBound \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if

		else if (strcmp(argv[argvCounter],"-multiplier") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &multiplier) != TCL_OK) {
			opserr << "ERROR: invalid input: multiplier \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if		

		else if (strcmp(argv[argvCounter],"-state") == 0) {
			argvCounter++;
			if (Tcl_GetDouble(interp, argv[argvCounter], &state) != TCL_OK) {
			opserr << "ERROR: invalid input: state \n";
			return TCL_ERROR;
			}
				argvCounter++;
		}// else if		

		else if (strcmp(argv[argvCounter],"-linearAdd") == 0) {
			char linearAddName[20];

			argvCounter++;
			strcpy(linearAddName,argv[argvCounter]);


         if( Tcl_GetVar2(interp, linearAddName,"1",TCL_GLOBAL_ONLY ) != NULL)  {

			int numOfDV = theOptimizationDomain ->getNumberOfDesignVariables();
			linearAdd = new Vector(numOfDV);

			const char *  theValue;	
			char index[5];
			for(int i=0; i<numOfDV; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, linearAddName,index,TCL_GLOBAL_ONLY );
				(*linearAdd)(i) = atof(theValue);
			};
			
		 }  //if 
		  else {
			opserr<<"warning: the linearAdd with name "<<linearAddName<< " does not exit"<<endln;
		 } // else

				argvCounter++;
		}// else if	"-linearAdd"

		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}

	};  // while


// here tag ;
  theObjectiveFunction = new ObjectiveFunction(   tag,
												  theOptimizationDomain,
												  interp,
												  isGradProvided,
											      linearAdd,
												  tclFileName,
												  name, 
												  gradientName, 
												  lowerBound, 
												  upperBound, 
												  multiplier, 									   
												  state		   
												  );




  if (theObjectiveFunction == 0) {
		opserr << "ERROR: could not create random theObjectiveFunction "<< endln;
		return TCL_ERROR;
	  }
 // release memory 
	if (gradientName !=0) delete gradientName;
	if (linearAdd !=0)	 delete linearAdd;


  // ADD THE OBJECT TO THE DOMAIN
  if (theOptimizationDomain->addObjectiveFunction(theObjectiveFunction) == false) {
	opserr << "ERROR: failed to add theObjectiveFunction  to the domain (wrong number of arguments?)\n";

	delete theObjectiveFunction;  
 
	return TCL_ERROR;
  }

  return TCL_OK;
}




int 
TclOptimizationModelBuilder_addConstraintFunction(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv){
 
	
  ConstraintFunction *theConstraintFunction = 0;
  int tag;
  int numberOfConstraints;

  char name[25] = "" ;
  char * gradientName = 0;
  char tclFileName[35] = "";

  Vector * lowerBound = 0;
  Vector * upperBound = 0;
  Vector * multiplier = 0; 
  Vector * state = 0; 
  Matrix *linearAdd=0;

  int numberOfArguments = argc;
  bool isGradProvided=false;





  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 6) {
		opserr << "ERROR: invalid number of arguments to designVariable command "<<endln;
	    opserr <<"command: constraintFunction 1  -name F1  -GradientName G1 -lowerBound  lBound1 -upperBound uBound1 -tclFile constraint1.tcl -multiplier M1 -state S1 -linearAdd A1 "<<endln;
		return TCL_ERROR;
  }


  // GET TAG NUMBER
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "ERROR: invalid input: tag \n";
	return TCL_ERROR;
  }


	// Loop through arguments
	int argvCounter = 2;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-name") == 0)||(strcmp(argv[argvCounter],"-Name") == 0)) {

			argvCounter++;
			strcpy(name,argv[argvCounter]);

			if ((Tcl_GetVar(interp, name, TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, name,"1",TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, name,"1,1",TCL_GLOBAL_ONLY ) != NULL)){
				opserr<<"Fatal: ConstraintFunction Name "<< name << " is already been used, please change another name"<<endln; 
				exit(-1);
			}
				
			
			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-GradientName") == 0)||(strcmp(argv[argvCounter],"-gradientName") == 0)) {
			
			argvCounter++;
			
			gradientName = new char[30];

			strcpy(gradientName,argv[argvCounter]);
			
			if ((Tcl_GetVar(interp, gradientName, TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, gradientName,"1",TCL_GLOBAL_ONLY ) != NULL)||(Tcl_GetVar2(interp, gradientName,"1,1",TCL_GLOBAL_ONLY ) != NULL)){
				opserr<<"Fatal: ConstraintFunction gradient Name "<< gradientName << " is already been used, please change another name"<<endln; 
				exit(-1);
			}
				
			isGradProvided = true;
			argvCounter++;
		}// else if
		
		else if ((strcmp(argv[argvCounter],"-tclFile") == 0)||(strcmp(argv[argvCounter],"-TclFile") == 0)) {
			argvCounter++;
			strcpy(tclFileName,argv[argvCounter]);
			argvCounter++;
		}// else if
		

		else if (strcmp(argv[argvCounter],"-lowerBound") == 0) {
			char lowerBoundName[20];
			argvCounter++;
			strcpy(lowerBoundName,argv[argvCounter]);

			// get Number of constraintFunction
			
			int ii=1;
			char index[5];
			sprintf(index,"%d",ii);
			while ( Tcl_GetVar2(interp, lowerBoundName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: lowerBound error!"<<endln; exit(-1);}
//			int numOfDV = theStructuralDomain ->getNumberOfDesignVariables();
			
			lowerBound = new Vector(numberOfConstraints);

			const char *  theValue;	
			
			for(int i=0; i<numberOfConstraints; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, lowerBoundName,index,TCL_GLOBAL_ONLY );
				(*lowerBound)(i) = atof(theValue);
			};

				argvCounter++;
		}// else if "-lowerBound"


		else if (strcmp(argv[argvCounter],"-upperBound") == 0) {
			char upperBoundName[20];
			argvCounter++;
			strcpy(upperBoundName,argv[argvCounter]);

			// get Number of constraintFunction
			
			int ii=1;
			char index[5];
			sprintf(index,"%d",ii);
			while ( Tcl_GetVar2(interp, upperBoundName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: upperBound error!"<<endln; exit(-1);}
//			int numOfDV = theStructuralDomain ->getNumberOfDesignVariables();
			
			upperBound = new Vector(numberOfConstraints);

			const char *  theValue;	
			
			for(int i=0; i<numberOfConstraints; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, upperBoundName,index,TCL_GLOBAL_ONLY );
				(*upperBound)(i) = atof(theValue);
			};

				argvCounter++;
		}// else if "-upperBound"


		else if (strcmp(argv[argvCounter],"-state") == 0) {
			char stateName[20];
			argvCounter++;
			strcpy(stateName,argv[argvCounter]);

			// get Number of constraintFunction
			
			int ii=1;
			char index[5];
			sprintf(index,"%d",ii);
			while ( Tcl_GetVar2(interp, stateName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: state error!"<<endln; exit(-1);}
//			int numOfDV = theStructuralDomain ->getNumberOfDesignVariables();
			
			state = new Vector(numberOfConstraints);

			const char *  theValue;	
			
			for(int i=0; i<numberOfConstraints; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, stateName,index,TCL_GLOBAL_ONLY );
				(*state)(i) = atof(theValue);
			};

				argvCounter++;
		}// else if "-state"



		else if (strcmp(argv[argvCounter],"-multiplier") == 0) {
			char multiplierName[20];
			argvCounter++;
			strcpy(multiplierName,argv[argvCounter]);

			// get Number of constraintFunction
			
			int ii=1;
			char index[5];
			sprintf(index,"%d",ii);
			while ( Tcl_GetVar2(interp, multiplierName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: lmultiplier error!"<<endln; exit(-1);}
//			int numOfDV = theStructuralDomain ->getNumberOfDesignVariables();
			
			multiplier = new Vector(numberOfConstraints);

			const char *  theValue;	
			
			for(int i=0; i<numberOfConstraints; i++){
				sprintf(index,"%d",i+1);   // begin with 1
				theValue = Tcl_GetVar2(interp, multiplierName,index,TCL_GLOBAL_ONLY );
				(*multiplier)(i) = atof(theValue);
			};

				argvCounter++;
		}// else if "-multiplier"

		
		
		else if (strcmp(argv[argvCounter],"-linearAdd") == 0) {
			char linearAddName[20];

			argvCounter++;
			strcpy(linearAddName,argv[argvCounter]);

// get Number of constraintFunction
		
			int ii=1;
			char index[10];
			sprintf(index,"%d",ii);
			strcat(index,",1");

			while ( Tcl_GetVar2(interp, linearAddName,index,TCL_GLOBAL_ONLY ) != NULL)  {
				ii++;
				sprintf(index,"%d",ii);
				strcat(index,",1");
			}
			
			numberOfConstraints = ii-1;
			if (numberOfConstraints ==0){opserr<<"Fatal: linearAdd error!"<<endln; exit(-1);}

			int numOfDVs = theOptimizationDomain ->getNumberOfDesignVariables();

  
         if( Tcl_GetVar2(interp, linearAddName,"1,1",TCL_GLOBAL_ONLY ) != NULL)  {

			linearAdd = new Matrix(numberOfConstraints,numOfDVs);

			const char *  theValue;	
			char temp[5];

			for(int i=0;i<numberOfConstraints; i++){
				for(int j=0; j<numOfDVs; j++){

			       sprintf(temp,"%d",i+1);   // begin with 1
				   strcpy(index,temp);
				   sprintf(temp,"%d",j+1);   // begin with 1
				   strcat(index,",");
				   strcat(index,temp);
					
				   theValue = Tcl_GetVar2(interp, linearAddName,index,TCL_GLOBAL_ONLY );
				   (*linearAdd)(i,j) = atof(theValue);
				}; //for
			} //for


		 }  //if 
		  else {
			opserr<<"warning: the linearAdd with name "<<linearAddName<< " does not exit"<<endln;
		 } // else

				argvCounter++;
		}// else if	"-linearAdd"



		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}

	};  // while


// here tag ;

  theConstraintFunction = new ConstraintFunction(tag, 
									   numberOfConstraints,
									   theOptimizationDomain, 
									   interp, 

									   isGradProvided, 
									   linearAdd,
									   tclFileName,
									   name, 
									   gradientName, 
									   
									   lowerBound, 
									   upperBound, 
									   multiplier, 									   
									   state
									   ); 




  if (theConstraintFunction == 0) {
		opserr << "ERROR: could not create random theConstraintFunction "<< endln;
		return TCL_ERROR;
	  }
  
// delete 

if (gradientName !=0) delete gradientName;
if (linearAdd !=0)	 delete linearAdd;
if (lowerBound !=0) delete lowerBound;
if (upperBound !=0) delete upperBound;
if (multiplier !=0) delete multiplier;
if (state !=0) delete state;



  // ADD THE OBJECT TO THE DOMAIN
  if (theOptimizationDomain->addConstraintFunction(theConstraintFunction) == false) {
	opserr << "ERROR: failed to add theConstraintFunction  to the domain (wrong number of arguments?)\n";

	delete theConstraintFunction; // otherwise memory leak
 


	return TCL_ERROR;
  }



 
  return TCL_OK;

}






// command: runSNOPTAnalysis -maxNumIter 100 -printOptPointX OptX.out -tclFileToRun tclFileToRun.tcl -printFlag 1
int 
TclOptimizationModelBuilder_runSNOPTAnalysis(ClientData clientData,Tcl_Interp *interp,int argc,TCL_Char **argv)
{

 

  int maxNumberOfIterations=100;
  int printFlag=0;
  char fileNamePrint[25]="";
  char probType[25]="SNOPTAnalysis";
  char * tclFileName=0;

  int numberOfArguments = argc;
  
  // CHECK THAT AT LEAST ENOUGH ARGUMENTS ARE GIVEN
  if (numberOfArguments < 3) {
		opserr << "ERROR: invalid number of arguments to designVariable command "<<endln;
	    opserr <<"command: runSNOPTAnalysis -maxNumIter 100 -printOptPointX OptX.out -tclFileToRun tclFileToRun.tcl"<<endln;
		return TCL_ERROR;
  }


	// Loop through arguments
	int argvCounter = 1;
	while (argc > argvCounter) {
		if ((strcmp(argv[argvCounter],"-maxNumIter") == 0)||(strcmp(argv[argvCounter],"-maxnumiter") == 0)) {

			argvCounter++;
			
			if (Tcl_GetInt(interp, argv[argvCounter], &maxNumberOfIterations) != TCL_OK) {
			opserr << "ERROR: invalid input: maxNumberOfIterations \n";
			return TCL_ERROR;
			}

			argvCounter++;
		}// if

		else if ((strcmp(argv[argvCounter],"-printFlag") == 0)||(strcmp(argv[argvCounter],"-printflag") == 0)) {

			argvCounter++;
			
			if (Tcl_GetInt(interp, argv[argvCounter], &printFlag) != TCL_OK) {
			opserr << "ERROR: invalid input: printFlag \n";
			return TCL_ERROR;
			}

			argvCounter++;
		}// if
		else if ((strcmp(argv[argvCounter],"-printOptPointX") == 0)||(strcmp(argv[argvCounter],"-printoptpointx") == 0)) {
			
			argvCounter++;
			
//  		gradientName = new char[30];

			strcpy(fileNamePrint,argv[argvCounter]);
			
			argvCounter++;
		}// else if
		
		else if ((strcmp(argv[argvCounter],"-tclFileToRun") == 0)||(strcmp(argv[argvCounter],"-TclFileToRun") == 0)) {
			argvCounter++;
			tclFileName = new char[30];
			strcpy(tclFileName,argv[argvCounter]);
			argvCounter++;
		}// else if
		
	
		else {
			opserr<<"warning: unknown command: "<<argv[argvCounter]<<endln;
			argvCounter++;
		
		}

	};  // while


// here tag ;
#ifdef _HAVESNOPT
  theSNOPTAnalysis = new SNOPTAnalysis(maxNumberOfIterations, 
					printFlag,
					fileNamePrint,
					probType, 
					theOptimizationDomain, 
					interp,
					tclFileName
					);




  if (theSNOPTAnalysis == 0) {
		opserr << "ERROR: could not create random theSNOPTAnalysis "<< endln;
		return TCL_ERROR;
	  }


	theSNOPTAnalysis->runOptAnalysis(theOptimizationDomain);
#else
	opserr << "FATAL ERROR: SNOPT NOT LINKED\n";
	exit(-1);
#endif
	if (tclFileName !=0) delete tclFileName;


  return TCL_OK;
}


 	