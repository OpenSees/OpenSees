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
#include <ElementParameter.h>

#include <RVParameter.h>
#include <NodeResponseParameter.h>
#include <LoadFactorParameter.h>

#ifdef _RELIABILITY

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
  int eleTag = -1;
  bool isele = false;

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


  if (argc >= 6 && strcmp(argv[0],"parameter") == 0 && strcmp(argv[2],"node") == 0 
      && strcmp(argv[4],"disp") == 0) {

      int nodeTag;
      if (Tcl_GetInt(interp, argv[3], &nodeTag) != TCL_OK) {
	return TCL_ERROR;    
      }
      Node *theNode = theTclDomain->getNode(nodeTag);

      int dof;
      if (Tcl_GetInt(interp, argv[5], &dof) != TCL_OK) {
	return TCL_ERROR;    
      }

      Parameter *newParameter = new NodeResponseParameter(paramTag, theNode, Disp, dof);

      theTclDomain->addParameter(newParameter);
      
      char buffer[40];
      sprintf(buffer, "%d", paramTag);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);
      
      return TCL_OK;
  }

  if (argc >= 5 && strcmp(argv[0],"parameter") == 0 && strcmp(argv[2],"pattern") == 0 
      && strcmp(argv[4],"lambda") == 0) {

      int patternTag;
      if (Tcl_GetInt(interp, argv[3], &patternTag) != TCL_OK) {
	return TCL_ERROR;    
      }
      LoadPattern *thePattern = theTclDomain->getLoadPattern(patternTag);

      Parameter *newParameter = new LoadFactorParameter(paramTag, thePattern);

      theTclDomain->addParameter(newParameter);
      
      char buffer[40];
      sprintf(buffer, "%d", paramTag);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);
      
      return TCL_OK;
  }

  RandomVariable *theRV = 0;

  // Now handle the parameter according to which command is invoked
  if (strcmp(argv[0],"parameter") == 0 || strcmp(argv[0],"addToParameter") == 0) {

    DomainComponent *theObject;

    if (strstr(argv[2],"randomVariable") != 0) {
#ifdef _RELIABILITY
      int rvTag;
      if (Tcl_GetInt(interp, argv[3], &rvTag) != TCL_OK) {
	return TCL_ERROR;    
      }
      
      if (theReliabilityDomain == 0) {
	opserr << "ERROR parameter " << paramTag << " -- reliability domain has not been created" << endln;      
      }
      
      theRV = theReliabilityDomain->getRandomVariablePtr(rvTag);
      if (theRV == 0) {
	opserr << "ERROR parameter " << paramTag << " -- random variable with tag " << rvTag << " not defined" << endln;
	return TCL_ERROR;
      }
#endif
    }

    int argStart = (theRV) ? 4 : 2;

    if (argc > argStart && strstr(argv[argStart],"element") != 0) {

      if (argc < 4) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[argStart+1], &eleTag) != TCL_OK) {
	opserr << "WARNING parameter -- invalid element tag\n";
	return TCL_ERROR;    
      }
      isele = true;
	 
      // Retrieve element from domain
      //  FMK theObject = (DomainComponent *) theTclDomain->getElement(eleTag);
      theObject = (DomainComponent *) theTclDomain->getElement(eleTag);

      argStart = (theRV) ? 6 : 4;
    }
    else if (argc > argStart && strstr(argv[argStart],"node") != 0) {
      if (argc < 4) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return TCL_ERROR;
      }

      int nodeTag;
      if (Tcl_GetInt(interp, argv[argStart+1], &nodeTag) != TCL_OK) {
	opserr << "WARNING parameter -- invalid node tag\n";
	return TCL_ERROR;    
      }

      // Retrieve element from domain
      theObject = (DomainComponent *) theTclDomain->getNode(nodeTag);

      argStart = (theRV) ? 6 : 4;
    }
    else if (argc > argStart && strstr(argv[argStart],"loadPattern") != 0) {
      if (argc < 4) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return TCL_ERROR;
      }

      int loadTag;
      if (Tcl_GetInt(interp, argv[argStart+1], &loadTag) != TCL_OK) {
	opserr << "WARNING parameter -- invalid load pattern tag\n";
	return TCL_ERROR;    
      }

      // Retrieve element from domain
      theObject = (DomainComponent *) theTclDomain->getLoadPattern(loadTag);

      argStart = (theRV) ? 6 : 4;
    }
    else if (argc > argStart) {
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

      Parameter *newParameter;
      if (argc > argStart) {
		if (isele == false) {
			newParameter = new Parameter(paramTag, theObject,
				       (const char **)&argv[argStart],
				       argc-argStart);
		} else { 
			newParameter = new ElementParameter(paramTag, eleTag,
					      (const char **)&argv[argStart],
					      argc-argStart);
		}
	  }  else
			newParameter = new Parameter(paramTag, 0, 0, 0);

      if (theRV != 0) {
		RVParameter *newRVParameter = new RVParameter(paramTag, theRV, newParameter);
		theTclDomain->addParameter(newRVParameter);
	  }
      else {
		theTclDomain->addParameter(newParameter);
      }

      char buffer[40];
      sprintf(buffer, "%d", paramTag);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);

	  }
    // Add to an existing parameter
    if (strcmp(argv[0],"addToParameter") == 0) {
     
      if (theParameter == 0) {
	       opserr << "WARNING addToParameter -- parameter with tag " << paramTag
	       << " not found in domain\n";
	       return TCL_ERROR;
      }
      else {
	if (isele == false) 
	  theParameter->addComponent(theObject, (const char **)&argv[argStart], argc-argStart);
	else {
	  theObject = (DomainComponent *) theTclDomain->getElement(eleTag);
	  theParameter->addComponent(theObject, (const char **)&argv[argStart], argc-argStart);	  
	  // Sorry, Frank, had to change this -- MHS
	  //theParameter->addComponent(eleTag, (const char **)&argv[argStart], argc-argStart);	  
	}
      }
    }

    return TCL_OK;
  }
  
  if (strcmp(argv[0],"updateParameter") == 0) {
    
    // Cannot update a parameter that is not present
    if (theParameter == 0) {
      opserr << "WARNING updateParameter -- parameter with tag " << paramTag
	     << " not found in domain\n";
    //  return TCL_ERROR;
    }
    
    double newValue;
    if (Tcl_GetDouble(interp, argv[2], &newValue) != TCL_OK) {
      opserr << "WARNING updateParameter -- invalid parameter value\n";
      return TCL_ERROR;
    }

    //    theParameter->update(newValue);
    theTclDomain->updateParameter(paramTag, newValue);

  }

  return TCL_OK;
}
