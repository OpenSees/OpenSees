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
                                                                        
// $Revision: 1.9 $
// $Date: 2010-09-13 21:39:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/BasicGFunEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <BasicGFunEvaluator.h>
#include <Vector.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>

#include <DummyStream.h>
#include <ElementResponse.h>
#include <Node.h>
#include <Element.h>

#include <tcl.h>
#include <string.h>


BasicGFunEvaluator::BasicGFunEvaluator(Tcl_Interp *passedTclInterp, 
				       ReliabilityDomain *passedReliabilityDomain,
				       Domain *passedOpenSeesDomain)

  :GFunEvaluator(), theTclInterp(passedTclInterp), theReliabilityDomain(passedReliabilityDomain), theOpenSeesDomain(passedOpenSeesDomain), g(0.0)
{
}

BasicGFunEvaluator::~BasicGFunEvaluator()
{
}

int
BasicGFunEvaluator::evaluateG(const Vector &x) 
{
  g = 0.0;

  // Set random variable values in Tcl namespace
  if (this->setTclRandomVariables(x) != 0) {
    opserr << "ERROR BasicGFunEvaluator::evaluateG -- error in setTclRandomVariables" << endln;
    return -1;
  }

  // "Download" limit-state function from reliability domain
  int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
  LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
  
  // Get the limit-state function expression
  const char *theExpression = theLimitStateFunction->getExpression();

  if (Tcl_ExprDouble( theTclInterp, theExpression, &g) != TCL_OK) {
    opserr << "BasicGFunEvaluator::evaluateGnoRecompute -- expression \"" << theExpression;
    opserr << "\" caused error:" << endln << theTclInterp->result << endln;
    return -1;
  }

  numberOfEvaluations++;

  return 0;
}

double
BasicGFunEvaluator::getG(void)
{
  return g;
}

int
BasicGFunEvaluator::setTclRandomVariables(const Vector &x)
{
  char theIndex[80];
  double xval;
  RandomVariable *theRV;
	
  // Set values of random variables in the Tcl interpreter
  int nrv = theReliabilityDomain->getNumberOfRandomVariables();

  int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();

  for (int i = 0; i < nrv; i++) {
    theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(i);
    int rvTag = theRV->getTag();

    xval = x(i);

    // put in x(1) format
    sprintf(theIndex,"%d",rvTag);
    if (Tcl_SetVar2Ex(theTclInterp,"xrv",theIndex,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG) == NULL) {
      opserr << "ERROR GFunEvaluator -- error in setTclRandomVariables xrv" << endln;
      opserr << theTclInterp->result << endln;
      return -1;
    }
    
    // put in x(1,lsfTag) format (useful for reporting design point)
    sprintf(theIndex,"%d,%d",rvTag,lsf);
    if (Tcl_SetVar2Ex(theTclInterp,"xrv",theIndex,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG) == NULL) {
      opserr << "ERROR GFunEvaluator -- error in setTclRandomVariables xrv" << endln;
      opserr << theTclInterp->result << endln;
      return -1;
    }
    
    // for legacy reasons, also put random variables in x_1 format
    sprintf(theIndex,"x_%d",rvTag);
    if (Tcl_SetVar2Ex(theTclInterp,theIndex,NULL,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG) == NULL) {
      opserr << "ERROR GFunEvaluator -- error in setTclRandomVariables x" << endln;
      opserr << theTclInterp->result << endln;
      return -1;
    }
  }

  return 0;
}

int
BasicGFunEvaluator::tokenizeSpecials(TCL_Char *theExpression, Tcl_Obj *paramList)
{
	// No need to set additional quantities here
	// (the basic random variables are set in the base class)

	return 0;
}

int
BasicGFunEvaluator::uParse(char *tempchar, int *node, int *dirn, char* disp, char* varName, char* arrName)
{
  // this is deprecated but keep in there for now
	
  // match regular expressions of known lsf parameters
  char* pattern = {"(ud*)\\(([0-9]+),([0-9]+)\\)"};
  const char* first;
  const char* last;
  char par_result[30];
	
  // save instances of parameters
  Tcl_RegExp regexp = Tcl_RegExpCompile(theTclInterp, pattern);
  if (regexp == NULL) {
    opserr << "GFunEvaluator::uParse ERROR compiling regular expression -- " << pattern << endln;
    opserr << theTclInterp->result << endln;
  }
	
  if ( Tcl_RegExpExec(theTclInterp, regexp, tempchar, tempchar) == 1 ) {
    // match found
		
    // get var name
    Tcl_RegExpRange(regexp, 1, &first, &last);
    if (first) {
      strncpy(par_result,first,last-first);
      par_result[last-first] = '\0';
      strcpy(varName, par_result);
      
      if ( strncmp(par_result, "udd", 3) == 0 )
	strcpy(disp,"accel");
      
      else if ( strncmp(par_result, "ud", 2) == 0 )
	strcpy(disp,"vel");
      
      else if ( strncmp(par_result, "u", 1) == 0 )
	strcpy(disp,"disp");
      
    }
    
    // get node number
    Tcl_RegExpRange(regexp, 2, &first, &last);
    if (first) {
      strncpy(par_result,first,last-first);
      par_result[last-first] = '\0';
      *node = atoi(par_result);
    }
    
    // get direction
    Tcl_RegExpRange(regexp, 3, &first, &last);
    if (first) {
      strncpy(par_result,first,last-first);
      par_result[last-first] = '\0';
      *dirn = atoi(par_result);
    }
    
    sprintf(arrName,"%d,%d",*node,*dirn);
  } 
  
  return 0;
}

int
BasicGFunEvaluator::nodeParse(char *tempchar, int *node, int *dirn, char* disp, char* varName, char* arrName)
{
  // match regular expressions of known lsf parameters
  char* pattern = {"((rec_)?node)\\(([0-9]+),([0-9]+),([a-z]+)\\)"};
  const char* first;
  const char* last;
  char par_result[30];
  
  // save instances of parameters
  Tcl_RegExp regexp = Tcl_RegExpCompile(theTclInterp, pattern);
  if (regexp == NULL) {
    opserr << "GFunEvaluator::nodeParse ERROR compiling regular expression -- " << pattern << endln;
    opserr << theTclInterp->result << endln;
  }
  
  if ( Tcl_RegExpExec(theTclInterp, regexp, tempchar, tempchar) == 1 ) {
    // match found
    
    // get var name
    Tcl_RegExpRange(regexp, 1, &first, &last);
    if (first) {
      strncpy(par_result,first,last-first);
      par_result[last-first] = '\0';
      strcpy(varName, par_result);
    }
    
    // get node number
    Tcl_RegExpRange(regexp, 3, &first, &last);
    if (first) {
      strncpy(par_result,first,last-first);
      par_result[last-first] = '\0';
      *node = atoi(par_result);
    }
    
    // get direction
    Tcl_RegExpRange(regexp, 4, &first, &last);
    if (first) {
      strncpy(par_result,first,last-first);
      par_result[last-first] = '\0';
      *dirn = atoi(par_result);
    }
    
    // get response type
    Tcl_RegExpRange(regexp, 5, &first, &last);
    if (first) {
      strncpy(par_result,first,last-first);
      par_result[last-first] = '\0';
      strcpy(disp, par_result);
    }
    
    sprintf(arrName,"%d,%d,%s",*node,*dirn,disp);
  } 
  
  return 0;
}

int
BasicGFunEvaluator::elementParse(char *tempchar, int *element, char* varName, char* eleArgs)
{
  int args = 0;
  
  // Start obtaining information about element number etc. 
  if ( strncmp(tempchar, "element", 7) == 0 ) {
    args = sscanf(tempchar,"element(%i,%s)",element,eleArgs);
    if (args == 2)
      strcpy(varName,"element");
  }
  else if ( strncmp(tempchar, "rec_element", 11) == 0 ) {
    args = sscanf(tempchar,"rec_element(%i,%s)",element,eleArgs);
    if (args == 2)
      strcpy(varName,"rec_element");
  }
  
  int strbrk = strcspn(eleArgs,")");
  eleArgs[strbrk] = '\0';	
  
  return 0;
}

int
BasicGFunEvaluator::nodeTclVariable(int nodeNumber, int direction, char* dispOrWhat, char* varName, char* arrName)
{
  // now obtain the information directly from domain without creating recorder files
  Node *theNode = theOpenSeesDomain->getNode(nodeNumber);
  if (theNode == 0) {
    opserr << "GFunEvaluator::nodeTclVariable -- node with tag " << nodeNumber
	   << " not found in OpenSees Domain" << endln;
    return 0;
  }
  
  double gFunValue = 0.0;
  
  if ( strncmp(dispOrWhat, "disp", 4) == 0) {
    const Vector &theDisp = theNode->getDisp();
    gFunValue = theDisp(direction-1);
  }
  else if ( strncmp(dispOrWhat, "vel", 3) == 0) {
    const Vector &theDisp = theNode->getVel();
    gFunValue = theDisp(direction-1);
  }
  else if ( strncmp(dispOrWhat, "accel", 5) == 0) {
    const Vector &theDisp = theNode->getAccel();
    gFunValue = theDisp(direction-1);
  }
  else {
    opserr << "ERROR GFunEvaluator::nodeTclVariable in syntax (" << dispOrWhat << ") of limit-state function " 
	   << "with node recorder." << endln;
  }
  
  // set Tcl value
  if (Tcl_SetVar2Ex(theTclInterp,varName,arrName,Tcl_NewDoubleObj(gFunValue),TCL_LEAVE_ERR_MSG) == NULL) {
    opserr << "ERROR  GFunEvaluator::nodeTclVariable -- SetVar error" << endln;
    opserr << theTclInterp->result << endln;
    return -1;
  }
  
  return 0;
}

int
BasicGFunEvaluator::elementTclVariable(int eleNumber, char* varName, char* inString)
{
  // now obtain the information directly from domain without creating recorder files  
  Element *theElement = theOpenSeesDomain->getElement(eleNumber);
  if (theElement == 0) {
    opserr << "GFunEvaluator::elementTclVariable -- element with tag " << eleNumber
	   << " not found in OpenSees Domain" << endln;
    return 0;
  }
  
  int rowNumber; // index into vector containing element response
  const int argvLength = 20;
  const int argcMax = 10;
  char workspace[argcMax*argvLength];
  char *argv[argcMax];
  for (int i = 0; i < argcMax; i++)
    argv[i] = &workspace[i*argvLength];
  
  int argc = 0;
  char restString[100];
  strcpy(restString, inString);
  
  if ( strncmp(restString, "section",7) == 0) {
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
      if (strncmp(restString,"stress",6)==0) {
	strcpy(argv[5],"stress");
	rowNumber = 1;
	argc = 6;
      }
      else if (strncmp(restString,"strain",6)==0) {
	strcpy(argv[5],"strain");
	rowNumber = 2;
	argc = 6;
      }
      else {
	opserr << "ERROR in syntax of limit-state function for element response quantity." << endln;
      }
    }
  }
  else {
    // should work for any user input response quantity type (as long as element supports it)
    char responseType[30] = "";
    char* tokstr = strtok(restString,"_");
    strcpy(responseType,tokstr);
    
    tokstr = strtok(NULL,"_");
    if (tokstr != NULL)
      rowNumber = atoi(tokstr);
    
    strcpy(argv[0], responseType);
    sprintf(argv[1], "%d", rowNumber);
    argc = 2;
  }
  
  char arrName[100];
  sprintf(arrName,"%i,%s",eleNumber,inString);
  opserr << "var = " << varName << " and array = " << arrName << endln;
  
  // add the element get and set response
  DummyStream theHandler;
  Response *theResponse;
  theResponse = theElement->setResponse((const char **)argv, argc, theHandler);
  
  double gFunValue = 0.0;
  if (theResponse != 0) {
    theResponse->getResponse();
    Information &eleInfo = theResponse->getInformation();
    const Vector &eleData = eleInfo.getData();
    if (eleData.Size() > 0)
      gFunValue = eleData(rowNumber-1); // C-index
    delete theResponse;
  }
  
  // set tcl variable
  if (Tcl_SetVar2Ex(theTclInterp,varName,arrName,Tcl_NewDoubleObj(gFunValue),TCL_LEAVE_ERR_MSG) == NULL) {
    opserr << "GFunEvaluator::elementTclVariable -- SetVar error" << endln;
    opserr << theTclInterp->result << endln;
    return -1;
  }
  
  return 0;
}

/*
int 
GFunEvaluator::setTclRandomVariables(const Vector &x)
{
  char theIndex[20];
  double xval;
  RandomVariable *theRV;
  
  // Set values of random variables in the Tcl interpreter
  int nrv = theReliabilityDomain->getNumberOfRandomVariables();
  //RandomVariableIter theRViter = theReliabilityDomain->getRandomVariables();
  //while ((theRV = theRViter()) != 0) {
  for (int i = 0; i < nrv; i++) {
    theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(i);
    int rvTag = theRV->getTag();
    sprintf(theIndex,"%d",rvTag);
    //xval = x( theReliabilityDomain->getRandomVariableIndex(rvTag) );
    xval = x( i );
    
    if (Tcl_SetVar2Ex(theTclInterp,"xrv",theIndex,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG) == NULL) {
      opserr << "ERROR GFunEvaluator -- error in setTclRandomVariables xrv" << endln;
      opserr << theTclInterp->result << endln;
      return -1;
    }
    
    // for legacy reasons, also put random variables in x_1 format
    sprintf(theIndex,"x_%d",rvTag);
    if (Tcl_SetVar2Ex(theTclInterp,theIndex,NULL,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG) == NULL) {
      opserr << "ERROR GFunEvaluator -- error in setTclRandomVariables x" << endln;
      opserr << theTclInterp->result << endln;
      return -1;
    }
  }
  
  return 0;
}
*/
