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

#include <Domain.h>
#include <Element.h>
#include <ElementParameter.h>
#include <ElementStateParameter.h>
#include <LoadFactorParameter.h>
#include <NodeResponseParameter.h>
#include <OPS_Stream.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <RVParameter.h>
#include <elementAPI.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#ifdef _RELIABILITY

#include <ReliabilityDomain.h>

ReliabilityDomain *OPS_GetReliabilityDomain();

#endif

int OPS_Parameter() {
  Domain *theDomain = OPS_GetDomain();
  #ifdef _RELIABILITY
    ReliabilityDomain *theReliabilityDomain = OPS_GetReliabilityDomain();
  #endif
  if (theDomain == 0) return -1;

  // check at least two arguments so don't segment fault on strcmp
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING need to specify a parameter tag\n";
    opserr << "Want: parameter tag <specific parameter args> .. see "
              "manual for valid parameter types and arguments\n";
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
  if (theParameter != 0) {
    opserr << "WARNNG: parameter " << paramTag << " already exists\n";
    return -1;
  }

  // First, check special case of a blank parameter
  if (OPS_GetNumRemainingInputArgs() == 0) {
    Parameter *newParameter = new Parameter(paramTag, 0, 0, 0);

    if (theDomain->addParameter(newParameter) == false) {
      opserr << "WARNING: failed to add parameter\n";
      return -1;
    }

    if (OPS_SetIntOutput(&num, &paramTag, true) < 0) {
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

    if (OPS_SetIntOutput(&num, &paramTag, true) < 0) {
      opserr << "WARING: parameter - failed to set parameter tag\n";
      return -1;
    }

    return 0;
  }

  // loop through all other parameters
  Node *node = 0;
  Element *element = 0;
  LoadPattern *pattern = 0;
  RandomVariable *theRV = 0;
  DomainComponent *theObject = 0;
  std::vector<const char *> argv;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();

    if (strcmp(type, "node") == 0) {
      // node object
      if (OPS_GetNumRemainingInputArgs() == 0) {
        opserr << "WARNING: need node tag\n";
        return -1;
      }
      if (theObject != 0) {
        opserr << "WARNING: another object is already set\n";
        return -1;
      }
      int tag;
      if (OPS_GetIntInput(&num, &tag) < 0) {
        opserr << "WARNING: failed to get node tag\n";
        return -1;
      }
      node = theDomain->getNode(tag);
      if (node == 0) {
        opserr << "WARNING: node " << tag << " not exist\n";
        return -1;
      }
      theObject = (DomainComponent *)node;

    } else if (strcmp(type, "pattern") == 0 ||
               strcmp(type, "loadPattern") == 0) {
      // pattern object
      if (OPS_GetNumRemainingInputArgs() == 0) {
        opserr << "WARNING: need load pattern tag\n";
        return -1;
      }
      if (theObject != 0) {
        opserr << "WARNING: another object is already set\n";
        return -1;
      }
      int tag;
      if (OPS_GetIntInput(&num, &tag) < 0) {
        opserr << "WARNING: failed to get pattern tag\n";
        return -1;
      }
      pattern = theDomain->getLoadPattern(tag);
      if (pattern == 0) {
        opserr << "WARNING: pattern " << tag << " not exists\n";
        return -1;
      }
      theObject = (DomainComponent *)pattern;

    } else if (strcmp(type, "element") == 0) {
      // element object
      if (OPS_GetNumRemainingInputArgs() == 0) {
        opserr << "WARNING: need element tag\n";
        return -1;
      }
      if (theObject != 0) {
        opserr << "WARNING: another object is already set\n";
        return -1;
      }
      int tag;
      if (OPS_GetIntInput(&num, &tag) < 0) {
        opserr << "WARNING parameter -- invalid element tag\n";
        return -1;
      }

      element = theDomain->getElement(tag);
      if (element == 0) {
        opserr << "WARNING: element " << tag << " not exists\n";
        return -1;
      }
      theObject = (DomainComponent *)element;

    } else if (strcmp(type, "randomVariable") == 0) {
#ifdef _RELIABILITY
      // associate with random variable
      if (OPS_GetNumRemainingInputArgs() == 0) {
        opserr << "WARNING: need random variable tag\n";
        return -1;
      }
      int rvTag;
      if (OPS_GetIntInput(&num, &rvTag) < 0) {
        opserr << "WARNING parameter -- invalid rv tag\n";
        return -1;
      }

      if (theReliabilityDomain == 0) {
        opserr << "WARNING: the reliability domain is not defined\n";
        return -1;
      }
      theRV = theReliabilityDomain->getRandomVariablePtr(rvTag);
      if (theRV == 0) {
        opserr << "ERROR parameter " << paramTag
               << " -- random variable with tag " << rvTag
               << " not defined\n";
        return -1;
      }
#endif

    } else if (node != 0 && strcmp(type, "disp") == 0) {
      // special node disp object
      if (OPS_GetNumRemainingInputArgs() == 0) {
        opserr << "WARNING: need disp dof\n";
        return -1;
      }

      int dof;
      if (OPS_GetIntInput(&num, &dof) < 0) {
        opserr << "WARNING: failed to get disp dof\n";
        return -1;
      }

      Parameter *newParameter =
          new NodeResponseParameter(paramTag, node, Disp, dof);

      if (theDomain->addParameter(newParameter) == false) {
        opserr << "WARNING: failed to add parameter\n";
        return -1;
      }

      if (OPS_SetIntOutput(&num, &paramTag, true) < 0) {
        opserr << "WARING: parameter - failed to set parameter tag\n";
        return -1;
      }

      return 0;

    } else if (pattern != 0 && strcmp(type, "lambda") == 0) {
      // special pattern lambda object
      Parameter *newParameter = new LoadFactorParameter(paramTag, pattern);

      if (theDomain->addParameter(newParameter) == false) {
        opserr << "WARNING: failed to add parameter\n";
        return -1;
      }

      if (OPS_SetIntOutput(&num, &paramTag, true) < 0) {
        opserr << "WARING: parameter - failed to set parameter tag\n";
        return -1;
      }

      return 0;

    } else if (theObject != 0) {
      // put the arg back
      OPS_ResetCurrentInputArg(-1);

      // collect parameter object data
      char *buffer = new char[128];
      if (OPS_GetStringFromAll(buffer, 128) == 0) {
        opserr << "WARNING: failed to read string from all\n";
        delete[] buffer;
        for (int i = 0; i < argv.size(); ++i) {
          delete[] argv[i];
        }
        return -1;
      }
      argv.push_back(buffer);

    } else {
      opserr << "WARING: parameter - inputs must be one of node, element, "
                "pattern, randomVariable\n";
      return -1;
    }
  }

  // Create new parameter
  Parameter *newParameter;
  if (element != 0) {
    newParameter = new ElementParameter(paramTag, element->getTag(),
                                        &argv[0], argv.size());
  } else if (theObject != 0) {
    newParameter =
        new Parameter(paramTag, theObject, &argv[0], argv.size());
  } else {
    newParameter = new Parameter(paramTag, 0, 0, 0);
  }

  // clean up argv
  for (int i = 0; i < argv.size(); ++i) {
    delete[] argv[i];
  }

  // add parameter
  if (theRV != 0) {
    RVParameter *newRVParameter =
        new RVParameter(paramTag, theRV, newParameter);
    theDomain->addParameter(newRVParameter);
  } else {
    theDomain->addParameter(newParameter);
  }

  if (OPS_SetIntOutput(&num, &paramTag, true) < 0) {
    opserr << "WARING: parameter - failed to set parameter tag\n";
    return -1;
  }

  return 0;
}

int OPS_addToParameter() {
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0) return -1;

  // check at least two arguments so don't segemnt fault on strcmp
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING need to specify a parameter tag\n";
    opserr << "Want: addToParameter tag <specific parameter args> .. see "
              "manual for valid parameter types and arguments\n";
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
    opserr << "WARNING parameter -- insufficient number of arguments for "
              "parameter with tag "
           << paramTag << '\n';
    return -1;
  }

  // check existence
  Parameter *theParameter = theDomain->getParameter(paramTag);
  if (theParameter == 0) {
    opserr << "WARNNG: parameter " << paramTag << " not exists\n";
    return -1;
  }

  // Now handle the parameter
  DomainComponent *theObject;
  const char *type = OPS_GetString();

  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING parameter -- insufficient number of arguments for "
              "parameter with tag "
           << paramTag << '\n';
    return -1;
  }

  if (strcmp(type, "element") == 0) {
    int eleTag;
    if (OPS_GetIntInput(&num, &eleTag) < 0) {
      opserr << "WARNING parameter -- invalid element tag\n";
      return -1;
    }

    // Retrieve element from domain
    theObject = (DomainComponent *)theDomain->getElement(eleTag);

  } else if (strcmp(type, "node") == 0) {
    int nodeTag;
    if (OPS_GetIntInput(&num, &nodeTag) < 0) {
      opserr << "WARNING parameter -- invalid node tag\n";
      return -1;
    }

    // Retrieve node from domain
    theObject = (DomainComponent *)theDomain->getNode(nodeTag);

  } else if (strcmp(type, "loadPattern") == 0) {
    int loadTag;
    if (OPS_GetIntInput(&num, &loadTag) < 0) {
      opserr << "WARNING parameter -- invalid load pattern tag\n";
      return -1;
    }

    // Retrieve load pattern from domain
    theObject = (DomainComponent *)theDomain->getLoadPattern(loadTag);

  } else {
    opserr << "WARNING - unable to assign parameter to object of type "
           << type << '\n';
    return -1;
  }

  ///////////////////////////////////

  // Add to an existing parameter
  int numr = OPS_GetNumRemainingInputArgs();
  if (numr > 0) {
    char **argv = new char *[numr];
    char buffer[128];

    for (int i = 0; i < numr; ++i) {
      argv[i] = new char[128];

      strcpy(argv[i], OPS_GetStringFromAll(buffer, 128));
    }

    theParameter->addComponent(theObject, (const char **)argv, numr);

    for (int i = 0; i < numr; ++i) {
      delete[] argv[i];
    }
    delete[] argv;
  }

  if (OPS_SetIntOutput(&num, &paramTag, true) < 0) {
    opserr << "WARING: parameter - failed to set parameter tag\n";
    return -1;
  }

  return 0;
}

int OPS_updateParameter() {
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0) return -1;

  // check at least two arguments so don't segemnt fault on strcmp
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING need to specify a parameter tag\n";
    opserr << "Want: updateParameter tag <specific parameter args> .. see "
              "manual for valid parameter types and arguments\n";
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
    opserr << "WARNING parameter -- insufficient number of arguments for "
              "parameter with tag "
           << paramTag << '\n';
    return -1;
  }

  // check existence
  Parameter *theParameter = theDomain->getParameter(paramTag);
  if (theParameter == 0) {
    opserr << "WARNNG: parameter " << paramTag << " not exists\n";
    return -1;
  }

  // update values
  double newValue;
  if (OPS_GetDouble(&num, &newValue) < 0) {
    opserr << "WARNING updateParameter -- invalid parameter value\n";
    return -1;
  }

  theDomain->updateParameter(paramTag, newValue);

  if (OPS_SetIntOutput(&num, &paramTag, true) < 0) {
    opserr << "WARING: parameter - failed to set parameter tag\n";
    return -1;
  }

  return 0;
}

int OPS_getParamTags() {
  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0) return -1;

  Parameter *theParam;
  ParameterIter &paramIter = theDomain->getParameters();

  std::vector<int> tags;
  while ((theParam = paramIter()) != 0) {
    tags.push_back(theParam->getTag());
  }

  int size = 0;
  int *data = 0;
  if (!tags.empty()) {
    size = (int)tags.size();
    data = &tags[0];
  }

  if (OPS_SetIntOutput(&size, data, false) < 0) {
    opserr << "WARNING failed to set outputs\n";
    return -1;
  }

  return 0;
}

int OPS_getParamValue() {
  Domain *theDomain = OPS_GetDomain();
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
    opserr << "WARNING parameter " << paramTag << " is not found\n";
    return -1;
  }

  double value = theParam->getValue();

  if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
    opserr << "WARNING failed to set output\n";
    return -1;
  }

  return 0;
}

int OPS_setParameter() {
  double newValue = 0.0;
  ID eleIDs(0, 32);
  int numEle = 0;
  int flag = 0;

  const char *opt = OPS_GetString();

  int numdata = 1;
  if (strcmp(opt, "-val") == 0) {
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

  if ((strcmp(opt, "-ele") == 0) || (strcmp(opt, "-eles") == 0) ||
      (strcmp(opt, "-element") == 0)) {
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

    if (numEle > 0) flag = 1;

  } else if (strcmp(opt, "-eleRange") == 0) {
    flag = 2;

    // ensure no segmentation fault if user messes up
    if (OPS_GetNumRemainingInputArgs() < 2) {
      opserr << "WARNING recorder Element .. -eleRange start? end?  .. - "
                "no ele tags specified\n";
      return -1;
    }

    //
    // read in start and end tags of two elements & add set [start,end]
    //

    int start, end;
    if (OPS_GetIntInput(&numdata, &start) < 0) {
      opserr << "WARNING recorder Element -eleRange start? end? - invalid "
                "start\n";
      return -1;
    }
    if (OPS_GetIntInput(&numdata, &end) < 0) {
      opserr << "WARNING recorder Element -eleRange start? end? - invalid "
                "end\n ";
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
  std::vector<const char *> argv;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    opt = OPS_GetString();
    if (strcmp(opt, "Invalid String Inpu!") == 0) {
      opserr << opt << "\n";
      return -1;
    }
    argv.push_back(opt);
  }

  if (argv.empty()) return 0;

  ElementStateParameter theParameter(newValue, &argv[0], (int)argv.size(),
                                     flag, &eleIDs);

  Domain *theDomain = OPS_GetDomain();
  if (theDomain == 0) return 0;
  theDomain->addParameter(&theParameter);

  return 0;
}
