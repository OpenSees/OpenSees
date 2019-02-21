/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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

// Written: Minjie

// Description: parameter commands

#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>
#include <Domain.h>
#include <elementAPI.h>

#include <Parameter.h>
#include <ElementParameter.h>
#include <ParameterIter.h>

#include <RVParameter.h>
#include <NodeResponseParameter.h>
#include <LoadFactorParameter.h>
#include <ElementStateParameter.h>

#include <vector>

#ifdef _RELIABILITY

// #include <ReliabilityDomain.h>

// extern ReliabilityDomain *theReliabilityDomain;

#endif

int
OPS_Parameter()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // check at least two arguments so don't segemnt fault on strcmp
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING need to specify a parameter tag\n";
	opserr << "Want: parameter tag <specific parameter args> .. see manual for valid parameter types and arguments\n";
	return -1;
    }

    // Figure out which parameter we are dealing with
    int paramTag;
    int num = 1;
    if (OPS_GetIntInput(&num, &paramTag) < 0) {
	opserr << "WARING: parameter - failed to get parameter tag\n";
	return -1;
    }

    // check existence
    Parameter *theParameter = theDomain->getParameter(paramTag);
    int eleTag = -1;
    bool isele = false;
    if (theParameter != 0) {
	opserr << "WARNNG: parameter "<<paramTag<<" already exists\n";
	return -1;
    }

    // First, check special case of a blank parameter
    if (OPS_GetNumRemainingInputArgs() == 0) {
	Parameter *newParameter = new Parameter(paramTag, 0, 0, 0);

	if (theDomain->addParameter(newParameter) == false) {
	    opserr << "WARNING: failed to add parameter\n";
	    return -1;
	}

	if (OPS_SetIntOutput(&num, &paramTag) < 0) {
	    opserr << "WARING: parameter - failed to set parameter tag\n";
	    return -1;
	}

	return 0;
    }

    // Second, check special case of one parameter
    if (OPS_GetNumRemainingInputArgs() == 1) {
	Parameter *newParameter = new Parameter(paramTag, 0, 0, 0);

	double value;
	if (OPS_GetDoubleInput(&num, &value) < 0) {
	    opserr << "WARNING: failed to get paramber value\n";
	    return -1;
	}

	newParameter->setValue(value);

	if (theDomain->addParameter(newParameter) == false) {
	    opserr << "WARNING: failed to add parameter\n";
	    return -1;
	}

	if (OPS_SetIntOutput(&num, &paramTag) < 0) {
	    opserr << "WARING: parameter - failed to set parameter tag\n";
	    return -1;
	}

	return 0;
    }

    // check special case of node disp
    if (OPS_GetNumRemainingInputArgs() >= 4) {

	const char* type1 = OPS_GetString();
	if (strcmp(type1,"node") == 0) {

	    int nodeTag;
	    if (OPS_GetIntInput(&num, &nodeTag) < 0) {
		opserr << "WARNING: failed to get node tag\n";
		return -1;
	    }
	    Node *theNode = theDomain->getNode(nodeTag);
	    if (theNode == 0) {
		opserr << "WARNING: node "<<nodeTag<<" not exist\n";
		return -1;
	    }

	    const char* type2 = OPS_GetString();
	    if (strcmp(type2,"disp") == 0) {
		int dof;
		if (OPS_GetIntInput(&num, &dof) < 0) {
		    opserr << "WARNING: failed to get disp dof\n";
		    return -1;
		}

		Parameter *newParameter = new NodeResponseParameter(paramTag, theNode, Disp, dof);

		if (theDomain->addParameter(newParameter) == false) {
		    opserr << "WARNING: failed to add parameter\n";
		    return -1;
		}

		if (OPS_SetIntOutput(&num, &paramTag) < 0) {
		    opserr << "WARING: parameter - failed to set parameter tag\n";
		    return -1;
		}

		return 0;
	    } else {

		OPS_ResetCurrentInputArg(-3);

	    }
	} else {
	    OPS_ResetCurrentInputArg(-1);
	}


    }

    // check special case of pattern load factor
    if (OPS_GetNumRemainingInputArgs() >= 3) {

	const char* type1 = OPS_GetString();
	if (strcmp(type1,"pattern") == 0) {

	    int patternTag;
	    if (OPS_GetIntInput(&num, &patternTag) < 0) {
		opserr << "WARNING: failed to get pattern tag\n";
		return -1;
	    }
	    LoadPattern *thePattern = theDomain->getLoadPattern(patternTag);
	    if (thePattern == 0) {
		opserr << "WARNING: pattern "<<patternTag<<" not exists\n";
		return -1;
	    }

	    const char* type2 = OPS_GetString();
	    if (strcmp(type2,"lambda") == 0) {


		Parameter *newParameter = new LoadFactorParameter(paramTag, thePattern);

		if (theDomain->addParameter(newParameter) == false) {
		    opserr << "WARNING: failed to add parameter\n";
		    return -1;
		}

		if (OPS_SetIntOutput(&num, &paramTag) < 0) {
		    opserr << "WARING: parameter - failed to set parameter tag\n";
		    return -1;
		}

		return 0;
	    } else {
		OPS_ResetCurrentInputArg(-3);
	    }
	} else {
	    OPS_ResetCurrentInputArg(-1);
	}
    }


    // Now handle the parameter
    RandomVariable *theRV = 0;
    DomainComponent *theObject;

    const char* type = OPS_GetString();

    if (strcmp(type,"randomVariable") == 0) {
#ifdef _RELIABILITY
	// int rvTag;
	// if (Tcl_GetInt(interp, argv[3], &rvTag) != TCL_OK) {
	// 	return TCL_ERROR;
	// }

	// if (theReliabilityDomain == 0) {
	// 	opserr << "ERROR parameter " << paramTag << " -- reliability domain has not been created" << endln;
	// }

	// theRV = theReliabilityDomain->getRandomVariablePtr(rvTag);
	// if (theRV == 0) {
	// 	opserr << "ERROR parameter " << paramTag << " -- random variable with tag " << rvTag << " not defined" << endln;
	// 	return TCL_ERROR;
	// }
#endif
    }

    if (theRV != 0) {
	type = OPS_GetString();
    }

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return -1;
    }

    if (strcmp(type,"element") == 0) {

	if (OPS_GetIntInput(&num, &eleTag) < 0) {
	    opserr << "WARNING parameter -- invalid element tag\n";
	    return -1;
	}
	isele = true;

	// Retrieve element from domain
	theObject = (DomainComponent *) theDomain->getElement(eleTag);

    } else if (strcmp(type,"node") == 0) {

	int nodeTag;
	if (OPS_GetIntInput(&num, &nodeTag) < 0) {
	    opserr << "WARNING parameter -- invalid node tag\n";
	    return -1;
	}

	// Retrieve node from domain
	theObject = (DomainComponent *) theDomain->getNode(nodeTag);

    } else if (strcmp(type,"loadPattern") == 0) {

	int loadTag;
	if (OPS_GetIntInput(&num, &loadTag) < 0) {
	    opserr << "WARNING parameter -- invalid load pattern tag\n";
	    return -1;
	}

	// Retrieve load pattern from domain
	theObject = (DomainComponent *) theDomain->getLoadPattern(loadTag);

    } else {
	opserr << "WARNING - unable to assign parameter to object of type "
	       << type << '\n';
	return -1;
    }

    ///////////////////////////////////

    // Create new parameter
    Parameter *newParameter;
    int numr = OPS_GetNumRemainingInputArgs();
    if (numr > 0) {

	char** argv = new char*[numr];

	for (int i=0; i<numr; ++i) {

	    argv[i] = new char[128];

	    double value;
	    int val;
	    if (OPS_GetIntInput(&num, &val) == 0) {

		// convert to string
		snprintf(argv[i], 128, "%d", val);

	    } else {

		// back one
		OPS_ResetCurrentInputArg(-1);

		if (OPS_GetDoubleInput(&num, &value) == 0) {

		    // convert to string
		    snprintf(argv[i], 128, "%.20f", value);

		} else {

		    // back one
		    OPS_ResetCurrentInputArg(-1);

		    // copy string
		    strcpy(argv[i], OPS_GetString());
		}
	    }

	}

	if (isele == false) {
	    newParameter = new Parameter(paramTag, theObject,
					 (const char**)argv, numr);
	} else {
	    newParameter = new ElementParameter(paramTag, eleTag,
						(const char**)argv, numr);
	}

	for (int i=0; i<numr; ++i) {

	    delete [] argv[i];
	}
	delete [] argv;

    } else {
	newParameter = new Parameter(paramTag, 0, 0, 0);
    }

    if (theRV != 0) {
	RVParameter *newRVParameter = new RVParameter(paramTag, theRV, newParameter);
	theDomain->addParameter(newRVParameter);
    } else {
	theDomain->addParameter(newParameter);
    }

    if (OPS_SetIntOutput(&num, &paramTag) < 0) {
	opserr << "WARING: parameter - failed to set parameter tag\n";
	return -1;
    }

    return 0;
}

int
OPS_addToParameter()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // check at least two arguments so don't segemnt fault on strcmp
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING need to specify a parameter tag\n";
	opserr << "Want: addToParameter tag <specific parameter args> .. see manual for valid parameter types and arguments\n";
	return -1;
    }

    // Figure out which parameter we are dealing with
    int paramTag;
    int num = 1;
    if (OPS_GetIntInput(&num, &paramTag) < 0) {
	opserr << "WARING: parameter - failed to get parameter tag\n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return -1;
    }

    // check existence
    Parameter *theParameter = theDomain->getParameter(paramTag);
    if (theParameter == 0) {
	opserr << "WARNNG: parameter "<<paramTag<<" not exists\n";
	return -1;
    }

    // Now handle the parameter
    DomainComponent *theObject;
    const char* type = OPS_GetString();

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return -1;
    }

    if (strcmp(type,"element") == 0) {

	int eleTag;
	if (OPS_GetIntInput(&num, &eleTag) < 0) {
	    opserr << "WARNING parameter -- invalid element tag\n";
	    return -1;
	}

	// Retrieve element from domain
	theObject = (DomainComponent *) theDomain->getElement(eleTag);

    } else if (strcmp(type,"node") == 0) {

	int nodeTag;
	if (OPS_GetIntInput(&num, &nodeTag) < 0) {
	    opserr << "WARNING parameter -- invalid node tag\n";
	    return -1;
	}

	// Retrieve node from domain
	theObject = (DomainComponent *) theDomain->getNode(nodeTag);

    } else if (strcmp(type,"loadPattern") == 0) {

	int loadTag;
	if (OPS_GetIntInput(&num, &loadTag) < 0) {
	    opserr << "WARNING parameter -- invalid load pattern tag\n";
	    return -1;
	}

	// Retrieve load pattern from domain
	theObject = (DomainComponent *) theDomain->getLoadPattern(loadTag);

    } else {
	opserr << "WARNING - unable to assign parameter to object of type "
	       << type << '\n';
	return -1;
    }

    ///////////////////////////////////

    // Add to an existing parameter
    int numr = OPS_GetNumRemainingInputArgs();
    if (numr > 0) {

	char** argv = new char*[numr];

	for (int i=0; i<numr; ++i) {

	    argv[i] = new char[128];


	    double value;
	    int val;
	    if (OPS_GetIntInput(&num, &val) == 0) {

		// convert to string
		snprintf(argv[i], 128, "%d", val);

	    } else {

		// back one
		OPS_ResetCurrentInputArg(-1);

		if (OPS_GetDoubleInput(&num, &value) == 0) {

		    // convert to string
		    snprintf(argv[i], 128, "%.20f", value);

		} else {

		    // back one
		    OPS_ResetCurrentInputArg(-1);

		    // copy string
		    strcpy(argv[i], OPS_GetString());
		}
	    }
	}

	theParameter->addComponent(theObject,(const char **)argv,numr);

	for (int i=0; i<numr; ++i) {

	    delete [] argv[i];
	}
	delete [] argv;

    }

    if (OPS_SetIntOutput(&num, &paramTag) < 0) {
	opserr << "WARING: parameter - failed to set parameter tag\n";
	return -1;
    }

    return 0;
}

int
OPS_updateParameter()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // check at least two arguments so don't segemnt fault on strcmp
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING need to specify a parameter tag\n";
	opserr << "Want: updateParameter tag <specific parameter args> .. see manual for valid parameter types and arguments\n";
	return -1;
    }

    // Figure out which parameter we are dealing with
    int paramTag;
    int num = 1;
    if (OPS_GetIntInput(&num, &paramTag) < 0) {
	opserr << "WARING: parameter - failed to get parameter tag\n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING parameter -- insufficient number of arguments for parameter with tag " << paramTag << '\n';
	return -1;
    }

    // check existence
    Parameter *theParameter = theDomain->getParameter(paramTag);
    if (theParameter == 0) {
	opserr << "WARNNG: parameter "<<paramTag<<" not exists\n";
	return -1;
    }

    // update values
    double newValue;
    if (OPS_GetDouble(&num, &newValue) < 0) {
	opserr << "WARNING updateParameter -- invalid parameter value\n";
	return -1;
    }

    theDomain->updateParameter(paramTag, newValue);

    if (OPS_SetIntOutput(&num, &paramTag) < 0) {
	opserr << "WARING: parameter - failed to set parameter tag\n";
	return -1;
    }

    return 0;
}


int OPS_getParamTags()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Parameter *theParam;
    ParameterIter &paramIter = theDomain->getParameters();

    std::vector<int> tags;
    while ((theParam = paramIter()) != 0) {
	tags.push_back(theParam->getTag());
    }

    if (tags.empty()) return 0;

    int size = (int)tags.size();
    int* data = &tags[0];

    if (OPS_SetIntOutput(&size, data) < 0) {
	opserr << "WARNING failed to set outputs\n";
	return -1;
    }

    return 0;

}

int OPS_getParamValue()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "Insufficient arguments to getParamValue" << endln;
	return -1;
    }

    int paramTag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &paramTag) < 0) {
	opserr << "WARNING getParamValue -- could not read paramTag \n";
	return -1;
    }

    Parameter *theParam = theDomain->getParameter(paramTag);
    if (theParam == 0) {
	opserr << "WARNING parameter "<<paramTag<<" is not found\n";
	return -1;
    }

    double value = theParam->getValue();

    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_setParameter()
{
    double newValue = 0.0;
    ID eleIDs(0, 32);
    int numEle = 0;
    int flag = 0;

    const char* opt = OPS_GetString();

    int numdata = 1;
    if (strcmp(opt,"-val") == 0) {
	if (OPS_GetDoubleInput(&numdata, &newValue) < 0) {
	    opserr << "WARNING: failed to get paramber value\n";
	    return -1;
	}
    } else {
	opserr << "WARNING setParameter:  -val not found \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() == 0) return 0;

    opt = OPS_GetString();

    if ((strcmp(opt,"-ele") == 0) ||
	(strcmp(opt,"-eles") == 0) ||
	(strcmp(opt,"-element") == 0)) {

	    //
	    // read in a list of ele until end of command or other flag
	    //

	    int eleTag;
	    while (OPS_GetNumRemainingInputArgs() > 0) {
		if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
		    // back on arg
		    OPS_ResetCurrentInputArg(-1);
		    break;
		} else {
		    eleIDs[numEle] = eleTag;
		    numEle++;
		}
	    }

	    if (numEle > 0)
		flag = 1;

    } else if (strcmp(opt,"-eleRange") == 0) {

	flag = 2;

	// ensure no segmentation fault if user messes up
	if (OPS_GetNumRemainingInputArgs() < 2) {
	    opserr << "WARNING recorder Element .. -eleRange start? end?  .. - no ele tags specified\n";
	    return -1;
	}

	//
	// read in start and end tags of two elements & add set [start,end]
	//

	int start, end;
	if (OPS_GetIntInput(&numdata, &start) < 0) {
	    opserr << "WARNING recorder Element -eleRange start? end? - invalid start\n";
	    return -1;
	}
	if (OPS_GetIntInput(&numdata, &end) < 0) {
	    opserr << "WARNING recorder Element -eleRange start? end? - invalid end\n ";
	    return -1;
	}
	if (start > end) {
	    int swap = end;
	    end = start;
	    start = swap;
	}
	eleIDs[0] = start;
	eleIDs[1] = end;
    }

    // all left args
    std::vector<const char*> argv;
    while (OPS_GetNumRemainingInputArgs() > 0) {
	opt = OPS_GetString();
	if (strcmp(opt, "Invalid String Inpu!") == 0) {
	    opserr << opt <<"\n";
	    return -1;
	}
	argv.push_back(opt);
    }

    if (argv.empty()) return 0;

    ElementStateParameter theParameter(newValue, &argv[0], (int)argv.size(), flag, &eleIDs);

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return 0;
    theDomain->addParameter(&theParameter);

    return 0;
}
