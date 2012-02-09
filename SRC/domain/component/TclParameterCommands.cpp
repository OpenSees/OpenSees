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
** ****************************************************************** */

// $Revision: 1.6 $
// $Date: 2010-06-09 17:31:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/TclParameterCommands.cpp,v $

#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>
#include <Domain.h>

#include <TclModelBuilder.h>

#include <Parameter.h>

#ifdef _RELIABILITY

#include <RVParameter.h>
#include <ReliabilityDomain.h>

extern ReliabilityDomain *theReliabilityDomain;

#endif

int
TclModelBuilderParameterCommand(ClientData clientData, Tcl_Interp *interp,
				int argc, TCL_Char **argv,
				Domain *theTclDomain,
				TclModelBuilder *theTclBuilder)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }
  
  // check at least two arguments so don't segemnt fault on strcmp
  if (argc < 2) {
    opserr << "WARNING need to specify a parameter tag\n";
    opserr << "Want: parameter tag <specific parameter args> .. see manual for valid parameter types and arguments\n";
    return TCL_ERROR;
  }

  // Figure out which parameter we are dealing with
  int paramTag;
  if (Tcl_GetInt(interp, argv[1], &paramTag) != TCL_OK) {

    return TCL_ERROR;    
  }

  Parameter *theParameter = theTclDomain->getParameter(paramTag);
  
  // First, check special case of a blank parameter
  if (argc == 2 && strcmp(argv[0],"parameter") == 0) {
    Parameter *newParameter = new Parameter(paramTag, 0, 0, 0);
	
    theTclDomain->addParameter(newParameter);

    char buffer[40];
    sprintf(buffer, "%d", paramTag);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
  }

  // First, check special case of a blank parameter
  if (argc == 3 && strcmp(argv[0],"parameter") == 0) {
    Parameter *newParameter = new Parameter(paramTag, 0, 0, 0);
	
    double value;
    if (Tcl_GetDouble(interp,argv[2],&value) != TCL_OK)
      return TCL_ERROR;

    newParameter->setValue(value);

    theTclDomain->addParameter(newParameter);

    char buffer[40];
    sprintf(buffer, "%d", paramTag);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
  }

#ifdef _RELIABILITY
  if (strcmp(argv[0],"parameter") == 0 && strcmp(argv[2],"randomVariable") == 0) {
    int rvTag;
    if (Tcl_GetInt(interp, argv[3], &rvTag) != TCL_OK) {
      return TCL_ERROR;    
    }

    if (theReliabilityDomain == 0) {
      opserr << "ERROR parameter " << paramTag << " -- reliability domain has not been created" << endln;      
    }

    RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtr(rvTag);
    if (theRV == 0) {
      opserr << "ERROR parameter " << paramTag << " -- random variable with tag " << rvTag << " not defined" << endln;
      return TCL_ERROR;
    }

    Parameter *newParameter = new RVParameter(paramTag, theRV);

    theTclDomain->addParameter(newParameter);

    return TCL_OK;
  }
#endif

  // Now handle the parameter according to which command is invoked
  if (strcmp(argv[0],"parameter") == 0 || strcmp(argv[0],"addToParameter") == 0) {

    int argStart = 0;

    DomainComponent *theObject;

    if (strstr(argv[2],"element") != 0) {

      if (argc < 4) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return TCL_ERROR;
      }

      int eleTag;
      if (Tcl_GetInt(interp, argv[3], &eleTag) != TCL_OK) {
	opserr << "WARNING parameter -- invalid element tag\n";
	return TCL_ERROR;    
      }

      // Retrieve element from domain
      theObject = (DomainComponent *) theTclDomain->getElement(eleTag);

      argStart = 4;
    }
    else if (strstr(argv[2],"node") != 0) {
      if (argc < 4) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return TCL_ERROR;
      }

      int nodeTag;
      if (Tcl_GetInt(interp, argv[3], &nodeTag) != TCL_OK) {
	opserr << "WARNING parameter -- invalid node tag\n";
	return TCL_ERROR;    
      }

      // Retrieve element from domain
      theObject = (DomainComponent *) theTclDomain->getNode(nodeTag);

      argStart = 4;
    }
    else if (strstr(argv[2],"loadPattern") != 0) {
      if (argc < 4) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return TCL_ERROR;
      }

      int loadTag;
      if (Tcl_GetInt(interp, argv[3], &loadTag) != TCL_OK) {
	opserr << "WARNING parameter -- invalid load pattern tag\n";
	return TCL_ERROR;    
      }

      // Retrieve element from domain
      theObject = (DomainComponent *) theTclDomain->getLoadPattern(loadTag);

      argStart = 4;
    }
    else {
      opserr << "WARNING - unable to assign parameter to object of type "
	     << argv[2] << '\n';
      return TCL_ERROR;
    }

    ///////////////////////////////////

    // Create new parameter
    if (strcmp(argv[0],"parameter") == 0) {
      
      if (theParameter != 0) {
	opserr << "WARNING parameter -- parameter with tag " << paramTag
	       << " already exists in domain\n";
	return TCL_ERROR;
      }
      else {
	Parameter *newParameter = new Parameter(paramTag, theObject,
						(const char **)&argv[argStart],
						argc-argStart);
	
	theTclDomain->addParameter(newParameter);
	char buffer[40];
	sprintf(buffer, "%d", paramTag);
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);
      }
    }
    // Add to an existing parameter
    if (strcmp(argv[0],"addToParameter") == 0) {
      
      if (theParameter == 0) {
	opserr << "WARNING addToParameter -- parameter with tag " << paramTag
	       << " not found in domain\n";
	return TCL_ERROR;
      }
      else {
	theParameter->addComponent(theObject, (const char **)&argv[argStart], argc-argStart);
      }
    }

    return TCL_OK;
  }
  
  if (strcmp(argv[0],"updateParameter") == 0) {
    
    // Cannot update a parameter that is not present
    if (theParameter == 0) {
      opserr << "WARNING updateParameter -- parameter with tag " << paramTag
	     << " not found in domain\n";
      return TCL_ERROR;
    }
    
    double newValue;
    if (Tcl_GetDouble(interp, argv[2], &newValue) != TCL_OK) {
      opserr << "WARNING updateParameter -- invalid parameter value\n";
      return TCL_ERROR;
    }

    theParameter->update(newValue);

  }

  return TCL_OK;
}
