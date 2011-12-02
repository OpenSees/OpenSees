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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-02-02 23:41:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/SimulationInformation.cpp,v $
//
// Description: This file contains the class definition for SimulationInformation.
// SimulationInformation is a stopwatch.
//
// What: "@(#) SimulationInformation.h, revA"



#include<SimulationInformation.h>
#include <OPS_Globals.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <tcl.h>

static int numSimulationInformation = 0;
static SimulationInformation *theLastSimulationInformation = 0;

SimulationInformation::SimulationInformation() 
  :title(0), description(0), contactName(0),
   lengthUnit(0), forceUnit(0), timeUnit(0)
{
  strcpy(startTime," ");
  strcpy(endTime," ");

  if (numSimulationInformation == 0) {
    numSimulationInformation++;
  }

  theLastSimulationInformation = this;
}

SimulationInformation::~SimulationInformation()
{ 
  if (title != 0)
    delete [] title;
  if (description != 0)
    delete [] description;
  if (contactName != 0)
    delete [] contactName;
  if ( lengthUnit != 0)
    delete [] lengthUnit;
  if ( forceUnit!= 0)
    delete [] forceUnit;
  if ( timeUnit != 0)
    delete [] timeUnit;

  if (theLastSimulationInformation == this)
    theLastSimulationInformation = 0;
}


int
SimulationInformation::start(void)
{
  inputFiles.clear();
  outputFiles.clear();
  paramNames.clear();
  paramValues.clear();
  analysisTypes.clear();
  modelTypes.clear();
  elementTypes.clear();
  materialTypes.clear();

  // now set the start time
  
  time_t timeT;
  if (time(&timeT) != 0) {
#ifdef _WIN32
    const char *sTime = ctime(&timeT);
    strcpy(startTime, sTime);
#else
    ctime_r(&timeT, &startTime[0]);
#endif

  }
  
  return 0;
}

int
SimulationInformation::end(void)
{
  time_t timeT;
  if (time(&timeT) != 0) {
#ifdef _WIN32
    const char *eTime = ctime(&timeT);
    strcpy(endTime, eTime);
#else
    ctime_r(&timeT, &endTime[0]);
#endif
  }
  
  return 0;
}

int 
SimulationInformation::setTitle(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (title != 0)
    delete [] title;

  // create new space & copy string
  title = new char[strlen(name)+1];
  if (title == 0)
    return -1;
  strcpy(title, name);

  return 0;
}

int 
SimulationInformation::setDescription(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (description != 0)
    delete [] description;

  // create new space & copy string
  description = new char[strlen(name)+1];
  if (description == 0)
    return -1;
  strcpy(description, name);

  return 0;
}

int 
SimulationInformation::setContact(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (contactName != 0)
    delete [] contactName;

  // create new space & copy string
  contactName = new char[strlen(name)+1];
  if (contactName == 0)
    return -1;
  strcpy(contactName, name);

  return 0;
}

int 
SimulationInformation::setLengthUnit(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (lengthUnit != 0)
    delete [] lengthUnit;

  // create new space & copy string
  lengthUnit = new char[strlen(name)+1];
  if (lengthUnit == 0)
    return -1;
  strcpy(lengthUnit, name);

  return 0;
}

int 
SimulationInformation::setForceUnit(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (forceUnit != 0)
    delete [] forceUnit;

  // create new space & copy string
  forceUnit = new char[strlen(name)+1];
  if (forceUnit == 0)
    return -1;
  strcpy(forceUnit, name);

  return 0;
}

int 
SimulationInformation::setTimeUnit(const char *name)
{
  // check for valid i/p
  if (name == 0)
    return -1;

  // delete the old
  if (timeUnit != 0)
    delete [] timeUnit;

  // create new space & copy string
  timeUnit = new char[strlen(name)+1];
  if (timeUnit == 0)
    return -1;
  strcpy(timeUnit, name);

  return 0;
}



int 
SimulationInformation::addInputFile(const char *fileName)
{

  if (strstr(fileName,"history.tcl") != 0)
    return 0;

  return inputFiles.addString(fileName);
}

int 
SimulationInformation::addOutputFile(const char *fileName)
{
  return outputFiles.addString(fileName);
}



int 
SimulationInformation::addParameter(const char *name, const char *value)
{
  // check for valid name & value
  if (name == 0 || value == 0) 
    return -1;
  
  // create new array to hold pointers and copy pointers there
  paramNames.addString(name);
  paramValues.addString(value);
  
  return 0;
}

int 
SimulationInformation::addModelType(const char *theType)
{
  return modelTypes.addString(theType);
}

int 
SimulationInformation::addAnalysisType(const char *theType)
{
  return analysisTypes.addString(theType);
}

int 
SimulationInformation::addLoadingType(const char *theType)
{
  return loadingTypes.addString(theType);
}

int 
SimulationInformation::addElementType(const char *theType)
{
  return elementTypes.addString(theType);
}

int 
SimulationInformation::addMaterialType(const char *theType)
{
  return materialTypes.addString(theType);
}


// TclSimulationInformation_defaultUnits()
// to define basic units. the following is based on code provided by S. Mazzoni

int
TclSimulationInformation_defaultUnits(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theLastSimulationInformation == 0 || argc < 7)
    return -1;

  const char *length = 0;
  const char *force = 0;
  const char *time = 0;
  
  int count = 1;
  while (count < 7) {
    if ((strcmp(argv[count],"-force") == 0) || (strcmp(argv[count],"-Force") == 0) 
	|| (strcmp(argv[count],"-FORCE") == 0)) {
      force = argv[count+1];
    } else if ((strcmp(argv[count],"-time") == 0) || (strcmp(argv[count],"-Time") == 0) 
	       || (strcmp(argv[count],"-TIME") == 0)) {
      time = argv[count+1];
    } else if ((strcmp(argv[count],"-length") == 0) || (strcmp(argv[count],"-Length") == 0) 
	       || (strcmp(argv[count],"-LENGTH") == 0)) {
      length = argv[count+1];
    } else {
      opserr << "units - unrecognized unit: " << argv[count] << " want: units -Force type? -Length type? - Time type\n";
      return -1;
    }
    count += 2;
  }

  if (length == 0 || force == 0 || time == 0) {
    opserr << "defaultUnits - missing a unit type want: units -Force type? -Length type? - Time type\n";
    return -1;
  }

  double in, ft, mm, cm, m;
  double lb, kip, n, kn;
  double sec, msec;
  

  if ((strcmp(length,"in") == 0) || (strcmp(length,"inch") == 0)) {
    in = 1.0;
  } else if ((strcmp(length,"ft") == 0) || (strcmp(length,"feet") == 0)) {
    in = 1.0 / 12.0;
  } else if ((strcmp(length,"mm") == 0)) {
    in = 25.4;
  } else if ((strcmp(length,"cm") == 0)) {
    in = 2.54;
  } else if ((strcmp(length,"m") == 0)) {
    in = 0.0254;
  } else {
    in = 1.0;
    opserr << "defaultUnits - unknown length type, valid options: in, ft, mm, cm, m\n";
    return TCL_ERROR;
  }

  if ((strcmp(force,"lb") == 0) || (strcmp(force,"lbs") == 0)) {
    lb = 1.0;
  } else if ((strcmp(force,"kip") == 0) || (strcmp(force,"kips") == 0)) {
    lb = 0.001;
  } else if ((strcmp(force,"N") == 0)) {
    lb = 4.4482216152605;
  } else if ((strcmp(force,"kN") == 0) || (strcmp(force,"KN") == 0) || (strcmp(force,"kn") == 0)) {
    lb = 0.0044482216152605;
  } else {
    lb = 1.0;
    opserr << "defaultUnits - unknown force type, valid options: lb, kip, N, kN\n";
    return TCL_ERROR;
  }

  if ((strcmp(time,"sec") == 0) || (strcmp(time,"sec") == 0)) {
    sec = 1.0;
  } else if ((strcmp(time,"msec") == 0) || (strcmp(time,"mSec") == 0)) {
    sec = 1000.0;
  } else {
    sec = 1.0;
    opserr << "defaultUnits - unknown time type, valid options: sec, msec\n";
    return TCL_ERROR;
  }

  ft = in * 12.0;
  mm = in / 25.44;
  cm = in / 2.54;
  m  = in / 0.0254;

  kip = lb / 0.001;
  n =   lb / 4.4482216152605;
  kn  = lb / 0.0044482216152605;

  msec = sec * 0.001;

  char string[50];


  sprintf(string,"set in %.18e", in);   Tcl_Eval(interp, string);
  sprintf(string,"set inch %.18e", in);   Tcl_Eval(interp, string);
  sprintf(string,"set ft %.18e", ft);   Tcl_Eval(interp, string);
  sprintf(string,"set mm %.18e", mm);   Tcl_Eval(interp, string);
  sprintf(string,"set cm %.18e", cm);   Tcl_Eval(interp, string);
  sprintf(string,"set m  %.18e", m);   Tcl_Eval(interp, string);
  sprintf(string,"set meter  %.18e", m);   Tcl_Eval(interp, string);

  sprintf(string,"set lb %.18e", lb);   Tcl_Eval(interp, string);
  sprintf(string,"set lbf %.18e", lb);   Tcl_Eval(interp, string);
  sprintf(string,"set kip %.18e", kip);   Tcl_Eval(interp, string);
  sprintf(string,"set N %.18e", n);   Tcl_Eval(interp, string);
  sprintf(string,"set kN %.18e", kn);   Tcl_Eval(interp, string);
  sprintf(string,"set Newton %.18e", n);   Tcl_Eval(interp, string);
  sprintf(string,"set kNewton %.18e", kn);   Tcl_Eval(interp, string);

  sprintf(string,"set sec %.18e", sec);   Tcl_Eval(interp, string);
  sprintf(string,"set msec %.18e", msec);   Tcl_Eval(interp, string);

  double g = 32.174049*ft/(sec*sec);
  sprintf(string,"set g %.18e", g);   Tcl_Eval(interp, string);
  sprintf(string,"set Pa %.18e",n/(m*m));   Tcl_Eval(interp, string);
  sprintf(string,"set MPa %.18e",1e6*n/(m*m));   Tcl_Eval(interp, string);
  sprintf(string,"set ksi %.18e",kip/(in*in));   Tcl_Eval(interp, string);
  sprintf(string,"set psi %.18e",lb/(in*in));   Tcl_Eval(interp, string);
  sprintf(string,"set pcf %.18e",lb/(ft*ft*ft));   Tcl_Eval(interp, string);
  sprintf(string,"set psf %.18e",lb/(ft*ft));   Tcl_Eval(interp, string);
  sprintf(string,"set in2 %.18e",in*in);   Tcl_Eval(interp, string);
  sprintf(string,"set m2 %.18e", m*m);   Tcl_Eval(interp, string);
  sprintf(string,"set mm2 %.18e",mm*mm);   Tcl_Eval(interp, string);
  sprintf(string,"set cm2 %.18e",cm*cm);   Tcl_Eval(interp, string);
  sprintf(string,"set in4 %.18e",in*in*in*in);   Tcl_Eval(interp, string);
  sprintf(string,"set mm4 %.18e",mm*mm*mm*mm);   Tcl_Eval(interp, string);
  sprintf(string,"set cm4 %.18e",cm*cm*cm*cm);   Tcl_Eval(interp, string);
  sprintf(string,"set m4 %.18e",m*m*m*m);   Tcl_Eval(interp, string);
  sprintf(string,"set pi %.18e",2.0*asin(1.0));   Tcl_Eval(interp, string);
  sprintf(string,"set PI %.18e",2.0*asin(1.0));   Tcl_Eval(interp, string);

  int res = theLastSimulationInformation->setLengthUnit(length);
  res += theLastSimulationInformation->setTimeUnit(time);
  res += theLastSimulationInformation->setForceUnit(force);

  return res;
}


int
TclSimulationInformation_neesMetaData(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theLastSimulationInformation == 0 || argc < 2)
    return -1;

  int count = 1;
  while (count < argc) {
    if ((strcmp(argv[count],"-title") == 0) || (strcmp(argv[count],"-Title") == 0) 
	|| (strcmp(argv[count],"-TITLE") == 0)) {
      if (count+1 < argc) {
	theLastSimulationInformation->setTitle(argv[count+1]);	
	count += 2;
      }
    } else if ((strcmp(argv[count],"-contact") == 0) || (strcmp(argv[count],"-Contact") == 0) 
	       || (strcmp(argv[count],"-CONTACT") == 0)) {
      if (count+1 < argc) {
	theLastSimulationInformation->setContact(argv[count+1]);	
	count += 2;
      }
    } else if ((strcmp(argv[count],"-description") == 0) || (strcmp(argv[count],"-Description") == 0) 
	       || (strcmp(argv[count],"-DESCRIPTION") == 0)) {
      if (count+1 < argc) {
	theLastSimulationInformation->setDescription(argv[count+1]);	
	count += 2;
      }
    } else if ((strcmp(argv[count],"-modelType") == 0) || (strcmp(argv[count],"-ModelType") == 0) 
	       || (strcmp(argv[count],"-MODELTYPE") == 0)) {
      if (count+1 < argc) {
	theLastSimulationInformation->addModelType(argv[count+1]);
	count += 2;
      }
    } else if ((strcmp(argv[count],"-analysisType") == 0) || (strcmp(argv[count],"-AnalysisType") == 0) 
	       || (strcmp(argv[count],"-ANALYSISTYPE") == 0)) {
      if (count+1 < argc) {
	theLastSimulationInformation->addAnalysisType(argv[count+1]);
	count += 2;
      }
    } else if ((strcmp(argv[count],"-elementType") == 0) || (strcmp(argv[count],"-ElementType") == 0) 
	       || (strcmp(argv[count],"-ELEMENTTYPE") == 0)) {
      if (count+1 < argc) {
	theLastSimulationInformation->addElementType(argv[count+1]);
	count += 2;
      }
    } else if ((strcmp(argv[count],"-materialType") == 0) || (strcmp(argv[count],"-MaterialType") == 0) 
	       || (strcmp(argv[count],"-MATERIALTYPE") == 0)) {
      if (count+1 < argc) {
	theLastSimulationInformation->addMaterialType(argv[count+1]);
	count += 2;
      }
    } else if ((strcmp(argv[count],"-loadingType") == 0) || (strcmp(argv[count],"-LoadingType") == 0) 
	       || (strcmp(argv[count],"-LOADINGTYPE") == 0)) {
      if (count+1 < argc) {
	theLastSimulationInformation->addLoadingType(argv[count+1]);
	count += 2;
      }
    } else {
      opserr << "WARNING unknown arg type: " << argv[count] << endln;
      count++;
    }
  }
  return TCL_OK;
}


int 
SimulationInformation::addTclInformationCommands(Tcl_Interp *interp)
{
  Tcl_CreateCommand(interp, "neesMetaData", TclSimulationInformation_neesMetaData,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "defaultUnits", TclSimulationInformation_defaultUnits,(ClientData)NULL, NULL);

  return 0;
}


void 
SimulationInformation::Print(OPS_Stream &s) const
{
  char version[10];
  strcpy(version,OPS_VERSION);
  
  s.tag("Central");
  s.tag("SimulationRun");

  if (title != 0)
    s.tag("title",title);

  if (description != 0)
    s.tag("description",description);

  if (contactName != 0)
    s.tag("contact",contactName);

  if (lengthUnit != 0)
    s.tag("lengthUnit", lengthUnit);

  if (forceUnit != 0)
    s.tag("lengthUnit", forceUnit);

  if (timeUnit != 0)
    s.tag("timeUnit", forceUnit);


  // need anotherTime to get rid of /n
  char *c = (char *)strchr(startTime,'\n');
  if (c != 0)
    strcpy(c,"");
  s.tag("startDate",startTime);

  c = (char *)strchr(endTime,'\n');
  if (c != 0)
    strcpy(c,"");
  s.tag("endDate",endTime);


  int numStrings = inputFiles.getNumStrings();
  if (numStrings != 0)
    s.tag("MainInputFile",inputFiles.getString(0));

  for (int i=1; i<numStrings; i++) 
    s.tag("InputFile",inputFiles.getString(i));    

  numStrings = outputFiles.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("OutputFile",outputFiles.getString(i));    

  numStrings = modelTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationModelType",modelTypes.getString(i));    

  numStrings = analysisTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationAnalysisType",analysisTypes.getString(i));    

  numStrings = loadingTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationLoadingType",loadingTypes.getString(i));    

  numStrings = elementTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationElementType",elementTypes.getString(i));    

  numStrings = materialTypes.getNumStrings();
  for (int i=0; i<numStrings; i++) 
    s.tag("SimulationMaterialType",materialTypes.getString(i));    

  s.tag("Software");
  s.tag("program","OpenSees");
  s.tag("version",version);
  s.endTag(); // Software


  s.tag("ComputerResource");
  s.tag("OS","Linux");
  s.tag("machine","local");
  s.endTag(); // Computer

  s.endTag(); // SimulationRun
  s.endTag(); // Central


  /*
  s << "Input Files:\n";
  for (int i=0; i<numFilesRead; i++)
    s << "  " << filesRead[i] << "\n";
  s << "\nOutput Files:\n";
  for (int k=0; k<numFilesWritten; k++)
    s << "  " << filesWritten[k] << "\n";
  s << "\nParameters:\n";
  for (int j=0; j<numParameters; j++)
    s << "  " << paramNames[j] << " " << paramValues[j] << "\n";
  s << endln;
  */

}    

OPS_Stream &operator<<(OPS_Stream &s, const SimulationInformation &E)
{
  E.Print(s);
  return s;
}


