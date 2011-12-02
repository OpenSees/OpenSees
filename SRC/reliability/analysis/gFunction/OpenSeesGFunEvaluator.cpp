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
                                                                        
// $Revision: 1.18 $
// $Date: 2007-04-30 20:03:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/OpenSeesGFunEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <OpenSeesGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
#include <RandomVariablePositionerIter.h>

#include <Node.h>
#include <Element.h>
#include <ElementResponse.h>

#include <DummyStream.h>

#include <tcl.h>
#include <string.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;

OpenSeesGFunEvaluator::OpenSeesGFunEvaluator(Tcl_Interp *passedTclInterp,
					     ReliabilityDomain *passedReliabilityDomain,
					     Domain *passedOpenSeesDomain,
					     TCL_Char *passedFileName)
  :GFunEvaluator(passedTclInterp, passedReliabilityDomain),
   theOpenSeesDomain(passedOpenSeesDomain)
{
	strcpy(fileName,passedFileName);
	nsteps = 0;
	dt = 0.0;
	// (here the user has provided a file with the analysis commmands)

	createTclVariables();

}

OpenSeesGFunEvaluator::OpenSeesGFunEvaluator(Tcl_Interp *passedTclInterp,
					     ReliabilityDomain *passedReliabilityDomain,
					     Domain *passedOpenSeesDomain,
					     int p_nsteps, double p_dt)
  :GFunEvaluator(passedTclInterp, passedReliabilityDomain),
   theOpenSeesDomain(passedOpenSeesDomain)
{
	fileName[0] = '\0';
	nsteps = p_nsteps;
	dt = p_dt;
	// (here the user has specified number of steps and possibly dt)

	createTclVariables();
}

OpenSeesGFunEvaluator::~OpenSeesGFunEvaluator()
{
  
}


int
OpenSeesGFunEvaluator::runGFunAnalysis(const Vector &x)
{
	// Zero out the response in the structural domain to make ready for next analysis
	char theRevertToStartCommand[10] = "reset";
	if (Tcl_Eval(theTclInterp, theRevertToStartCommand) == TCL_ERROR) {
	  opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_Eval for the reset command" << endln;
	  return -1;
	}


	// Put random variables into the structural domain according to the RandomVariablePositioners
	int rvNumber;
	RandomVariablePositionerIter &rvPosIter =
	  theReliabilityDomain->getRandomVariablePositioners();
	RandomVariablePositioner *theRVPos;
	while ((theRVPos = rvPosIter()) != 0) {
	  rvNumber = theRVPos->getRvNumber();
	  theRVPos->update(x(rvNumber-1));
	}


	// Run the structural analysis according to user specified scheme
	double result = 0;
	if (dt==0.0 && nsteps==0 && !strcmp(fileName,"0") ) {
		// Run up to max time in fFuncs
		opserr << "OpenSeesGFunEvaluator: The option -runToMaxTimeInGFun " << endln
			<< " is not yet implemented." << endln;
	}
	else if (dt==0.0 && nsteps==0) {
		// Read commands from file and execute them
		// (The reason for doing it like this is that we want to check
		// the flag returned by the analysis to see if it converged). 
		char theAnalyzeCommand[30];
		sprintf(theAnalyzeCommand,"[source %s]",fileName);
		if (Tcl_ExprDouble( theTclInterp, theAnalyzeCommand, &result) == TCL_ERROR) {
		  opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_ExprDouble for the analyze file command" << endln;
		  return -1;
		}

	}
	else {
		// User has given "nsteps" and possibly "dt"
		char theAnalyzeCommand[30];
		if (dt == 0.0) {
			sprintf(theAnalyzeCommand,"[analyze %d]",nsteps);
			
			if (Tcl_ExprDouble( theTclInterp, theAnalyzeCommand, &result) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_ExprDouble for the analyze command" << endln;
			  return -1;
			}
		}
		else {
			sprintf(theAnalyzeCommand,"[analyze %d %10.5f]", nsteps, dt);
			if (Tcl_ExprDouble( theTclInterp, theAnalyzeCommand, &result) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_ExprDouble for the analyze command" << endln;
			  return -1;
			}
			
		}
	
	}

	return ((int)result);
}




int
OpenSeesGFunEvaluator::tokenizeSpecials(TCL_Char *theExpression)
{
  // Set value of OpenSees finite element response quantities 
  // appearing in the limit-state function in the Tcl domain
  char tempchar[1000];
  double fileValue = 0.0;
  
  //opserr << "OSGFE:: " << theExpression << endln;
  
  
  char separators[5] = "}{";
  char lsf_forTokenizing[1000];
  strcpy(lsf_forTokenizing,theExpression);
  char lsf_expression[1000] = "";
  char *dollarSign = "$";
  char *underscore = "_";
  char *tokenPtr = strtok( lsf_forTokenizing, separators);

  while ( tokenPtr != NULL ) {
    
    strcpy(tempchar,tokenPtr);
    
    //opserr << "OSGFE::tok " << tokenPtr << endln;
    
    // If a nodal velocity is detected
    if ( strncmp(tokenPtr, "ud", 2) == 0) {
      
      // Get node number and dof number
      int nodeNumber, direction;
      sscanf(tempchar,"ud_%i_%i", &nodeNumber, &direction);
      
      // Assign value to the displacement quantity
      char tclAssignment[100];
      sprintf(tclAssignment,"set ud_%d_%d [nodeVel %d %d ]",nodeNumber,direction,nodeNumber,direction);
      if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
	opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_Eval for tokenizeSpecials" << endln;
	return -1;
      }
    }
    // If a nodal displacement is detected
    else if ( strncmp(tokenPtr, "u", 1) == 0) {
      
      // Get node number and dof number
      int nodeNumber, direction;
      sscanf(tempchar,"u_%i_%i", &nodeNumber, &direction);
      
      // Assign value to the displacement quantity
      char tclAssignment[100];
      sprintf(tclAssignment,"set u_%d_%d [nodeDisp %d %d ]",nodeNumber,direction,nodeNumber,direction);
      if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
	opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_Eval for tokenizeSpecials" << endln;
	return -1;
      }
      
    }
    else if ( strncmp(tokenPtr, "rec",3) == 0) {
      
      // Determine file name 
      int lineNum = 0;
      int colNum = 0;
      char variableName[100];
      if ( strncmp(tokenPtr, "rec_node",8) == 0) {
	rec_nodeTclVariable(tempchar, variableName);
      }
      else if ( strncmp(tokenPtr, "rec_element",11) == 0) {
	rec_elementTclVariable(tempchar, variableName);
      }
      
    }
    
    tokenPtr = strtok( NULL, separators);
  }
  
  return 0;
}


void
OpenSeesGFunEvaluator::setNsteps(int p_nsteps)
{
	nsteps = p_nsteps;
}

double
OpenSeesGFunEvaluator::getDt()
{
	return dt;
}





int
OpenSeesGFunEvaluator::createTclVariables()
{
  // Download active limit-state function
  int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
  LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
  char *theExpression = theLimitStateFunction->getExpression();

  // Initial declarations
  char separators[5] = "}{";
  char lsf_forTokenizing[500];
  strcpy(lsf_forTokenizing,theExpression);

  // Go through the limit-state function
  char *tokenPtr = strtok( lsf_forTokenizing, separators);
  while ( tokenPtr != NULL ) {

    if ( strncmp(tokenPtr, "rec_node",8) == 0 ||
	 strncmp(tokenPtr, "rec_element",11) == 0 ) {
      if (Tcl_SetVar(theTclInterp, tokenPtr, "0.0", TCL_LIST_ELEMENT) == NULL) {
	opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_SetVar for createTclVariables" << endln;
	return -1;
      }
    }
    
    tokenPtr = strtok( NULL, separators);
  }

  return 0;
}




int
OpenSeesGFunEvaluator::removeTclVariables()
{
  // Download active limit-state function
  int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
  LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
  char *theExpression = theLimitStateFunction->getExpression();

  // Initial declarations
  char separators[5] = "}{";
  char lsf_forTokenizing[500];
  strcpy(lsf_forTokenizing,theExpression);

  // Go through the limit-state function
  char *tokenPtr = strtok( lsf_forTokenizing, separators);
  while ( tokenPtr != NULL ) {

    if ( strncmp(tokenPtr, "rec_node",8) == 0 ||
	 strncmp(tokenPtr, "rec_element",11) == 0 ) {
      if (Tcl_UnsetVar(theTclInterp, tokenPtr, TCL_LIST_ELEMENT) == TCL_ERROR) {
	opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_UnsetVar for removeTclVariables" << endln;
	return -1;
      }
    }
    
    tokenPtr = strtok( NULL, separators);
  }

  return 0;
}

int
OpenSeesGFunEvaluator::rec_nodeTclVariable(char *tempchar, char *variableName)
{
  // Get node and dof number etc.
  int nodeNumber, direction;
  char dispOrWhat[10];
  
  char restString[100];
  sscanf(tempchar,"rec_node_%s",restString);
  
  if (strncmp(restString, "disp",4) == 0) {
    strcpy(dispOrWhat,"disp");
    sscanf(restString,"disp_%i_%i", &nodeNumber, &direction);
  }
  else if (strncmp(restString, "vel",3) == 0) {
    strcpy(dispOrWhat,"vel");
    sscanf(restString,"vel_%i_%i", &nodeNumber, &direction);
  }
  else if (strncmp(restString, "accel",5) == 0) {
    strcpy(dispOrWhat,"accel");
    sscanf(restString,"accel_%i_%i", &nodeNumber, &direction);
  }
  else {
    opserr << "ERROR in syntax of limit-state function with recorder." << endln;
  }

  Node *theNode = theOpenSeesDomain->getNode(nodeNumber);
  if (theNode == 0) {
    opserr << "OpenSeesGFunEvaluator -- node with tag " << nodeNumber
	   << " not found in OpenSees Domain" << endln;
    return 0;
  }

  double gFunValue = 0.0;
    
  if ( strcmp(dispOrWhat, "disp") == 0) {
    const Vector &theDisp = theNode->getDisp();
    gFunValue = theDisp(direction-1);
  }
  else if ( strcmp(dispOrWhat, "vel") == 0) {
    const Vector &theDisp = theNode->getVel();
    gFunValue = theDisp(direction-1);
  }
  else if ( strcmp(dispOrWhat, "accel") == 0) {
    const Vector &theDisp = theNode->getAccel();
    gFunValue = theDisp(direction-1);
  }
  else {
    opserr << "ERROR in syntax of limit-state function with recorder." << endln;
  }

  // Determine variable name 
  char tempString[100];
  sprintf(tempString, "rec_node_%s_%d_%d", dispOrWhat, nodeNumber, direction);
  
  char gString[80];
  sprintf(gString, "%25.20e", gFunValue);
  if (Tcl_SetVar(theTclInterp, tempString, gString, TCL_LIST_ELEMENT) == NULL) {
    opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_SetVar for rec_nodeTclVariable" << endln;
    return -1;
  }

  return 0;
}

int
OpenSeesGFunEvaluator::rec_elementTclVariable(char *tempchar, char *variableName)
{
  // Initial declarations
  char restString[100];
  
  // Start obtaining information about element number etc. 
  int eleNumber;
  sscanf(tempchar,"rec_element_%i_%s", &eleNumber, restString);

  Element *theElement = theOpenSeesDomain->getElement(eleNumber);
  if (theElement == 0) {
    opserr << "OpenSeesGFunEvaluator -- element with tag " << eleNumber
	   << " not found in OpenSees Domain" << endln;
    return 0;
  }
  
  int rowNumber; // index in to vector containing element response
  const int argvLength = 20;
  const int argcMax = 10;
  char workspace[argcMax*argvLength];
  char *argv[argcMax];
  for (int i = 0; i < argcMax; i++)
    argv[i] = &workspace[i*argvLength];

  
  int argc = 0;

  if ( strncmp(restString, "globalForce",11) == 0) {
    sscanf(restString,"globalForce_%i", &rowNumber);
    strcpy(argv[0], "globalForce");
    sprintf(argv[1], "%d", rowNumber);
    argc = 2;
  }
  else if ( strncmp(restString, "localForce",10) == 0) {
    sscanf(restString,"localForce_%i", &rowNumber);
    strcpy(argv[0], "localForce");
    sprintf(argv[1], "%d", rowNumber);
    argc = 2;
  }
  else if ( strncmp(restString, "plasticDeformation",18) == 0) {
    sscanf(restString,"plasticDeformation_%i", &rowNumber);
    strcpy(argv[0], "plasticDeformation");
    sprintf(argv[1], "%d", rowNumber);
    argc = 2;
  }
  else if ( strncmp(restString, "section",7) == 0) {
    int sectionNumber;
    sscanf(restString,"section_%i_%s", &sectionNumber, restString);
    strcpy(argv[0], "section");
    sprintf(argv[1], "%d", sectionNumber);
    if ( strncmp(restString, "force",5) == 0) {
      sscanf(restString,"force_%i", &rowNumber);
      strcpy(argv[2], "force");
      sprintf(argv[3], "%d", rowNumber);
      argc = 4;
    }
    else if ( strncmp(restString, "deformation",11) == 0) {
      sscanf(restString,"deformation_%i", &rowNumber);
      strcpy(argv[2], "deformation");
      sprintf(argv[3], "%d", rowNumber);
      argc = 4;
    }
    else if ( strncmp(restString, "fiber",5) == 0) {
      int ya, yb, za, zb;
      sscanf(restString,"fiber_%i_%i_%i_%i_%s", &ya, &yb, &za, &zb, restString);
      strcpy(argv[2],"fiber");
      sprintf(argv[3],"%d.%d", ya, yb);
      sprintf(argv[4],"%d.%d", za, zb);
      if (strcmp(restString,"stress")==0) {
	strcpy(argv[5],"stress");
	rowNumber = 1;
	argc = 6;
      }
      else if (strcmp(restString,"strain")==0) {
	strcpy(argv[5],"strain");
	rowNumber = 2;
	argc = 6;
      }
      else {
	opserr << "ERROR in syntax of limit-state function for element response quantity." << endln;
      }
    }
  }

  DummyStream theHandler;
  Response *theResponse;
  theResponse = theElement->setResponse((const char **)argv, argc, theHandler);
  
  double gFunValue = 0.0;
  if (theResponse != 0) {
    theResponse->getResponse();
    Information &eleInfo = theResponse->getInformation();
    gFunValue = eleInfo.theVector->operator()(rowNumber-1); // C-index
    delete theResponse;
  }
  
  char gString[80];
  sprintf(gString, "%25.20e", gFunValue);
  if (Tcl_SetVar(theTclInterp, tempchar, gString, TCL_LIST_ELEMENT) == NULL) {
    opserr << "ERROR OpenSeesGFunEvaluator -- error in Tcl_SetVar for rec_elementTclVariable" << endln;
    return -1;
  }

  return 0;
}
