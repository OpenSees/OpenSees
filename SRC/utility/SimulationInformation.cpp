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
                                                                        
// $Revision: 1.19 $
// $Date: 2008-12-18 22:46:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/SimulationInformation.cpp,v $
//
// Description: This file contains the class definition for SimulationInformation.
// SimulationInformation is a stopwatch.
//
// What: "@(#) SimulationInformation.h, revA"



#include<SimulationInformation.h>
#include <FileIter.h>
#include <OPS_Globals.h>
#include <string.h>
#include <time.h>
#include <math.h>
//#include <tcl.h>
#include <stdlib.h>


static int numSimulationInformation = 0;
static SimulationInformation *theLastSimulationInformation = 0;

SimulationInformation::SimulationInformation() 
  :title(0), description(0), contactName(0),
   lengthUnit(0), forceUnit(0), timeUnit(0),
   temperatureUnit(0)
{
  strcpy(startTime," ");
  strcpy(endTime," ");

  if (numSimulationInformation == 0) {
    numSimulationInformation++;
  }

  theLastSimulationInformation = this;

  theFiles = new File("", "", true);

  this->start();
}

//#include <XmlFileStream.h>

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
  if (temperatureUnit != 0)
      delete[] temperatureUnit;

  if (theLastSimulationInformation == this)
    theLastSimulationInformation = 0;

  delete theFiles;
}


int
SimulationInformation::start(void)
{
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

  strcpy(endTime," ");
  numInputFiles = 0;
  
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
SimulationInformation::setTemperatureUnit(const char *name)
{
    // check for valid i/p
    if (name == 0)
        return -1;

    // delete the old
    if (temperatureUnit != 0)
        delete[] temperatureUnit;

    // create new space & copy string
    temperatureUnit = new char[strlen(name) + 1];
    if (temperatureUnit == 0)
        return -1;
    strcpy(temperatureUnit, name);

    return 0;
}

int 
SimulationInformation::addInputFile(const char *fileName, const char *path)
{

  // windows throws this one in if nothing provided on the cmd line
 if (strstr(fileName,"history.tcl") != 0)
    return 0;

 
  if (numInputFiles == 0)
    theFiles->addFile(fileName, path, "Main Input File");
  else
    theFiles->addFile(fileName, path, "Input File");

  numInputFiles++;

  return 0;
}



int 
SimulationInformation::addOutputFile(const char *fileName, const char *path)
{
  theFiles->addFile(fileName, path, "Output File");

  return 0;
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


void 
PrintFiles(OPS_Stream &s, File *theFile) 
{

  if (theFile == 0)
    return;
  
  const char *fileName = theFile->getName();

  if (theFile->isDir() == true) {
    if (fileName != 0) {
      s.tag("Directory");
      s.attr("name", fileName);

    }

    FileIter theDirFiles = theFile->getFiles();
    File *theDirFile;
    while ((theDirFile = theDirFiles()) != 0)
      PrintFiles(s, theDirFile);
  } else {
    s.tag("File");
    s.attr("name", fileName);
    }

  s.endTag();
}

void 
SimulationInformation::Print(OPS_Stream &s, int flag) const
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s.tag("Central");
        s.tag("SimulationRun");

        if (title != 0)
            s.tag("title", title);

        if (description != 0)
            s.tag("description", description);

        if (contactName != 0)
            s.tag("contact", contactName);

        if (lengthUnit != 0)
            s.tag("lengthUnit", lengthUnit);

        if (forceUnit != 0)
            s.tag("forceUnit", forceUnit);

        if (timeUnit != 0)
            s.tag("timeUnit", timeUnit);

        if (temperatureUnit != 0)
            s.tag("temperatureUnit", temperatureUnit);

        // need anotherTime to get rid of /n
        char *c = (char *)strchr(startTime, '\n');
        if (c != 0)
            strcpy(c, "");
        s.tag("startDate", startTime);

        c = (char *)strchr(endTime, '\n');
        if (c != 0)
            strcpy(c, "");
        s.tag("endDate", endTime);

        int numStrings;

        numStrings = modelTypes.getNumStrings();
        for (int i = 0; i < numStrings; i++)
            s.tag("SimulationModelType", modelTypes.getString(i));

        numStrings = analysisTypes.getNumStrings();
        for (int i = 0; i < numStrings; i++)
            s.tag("SimulationAnalysisType", analysisTypes.getString(i));

        numStrings = loadingTypes.getNumStrings();
        for (int i = 0; i < numStrings; i++)
            s.tag("SimulationLoadingType", loadingTypes.getString(i));

        numStrings = elementTypes.getNumStrings();
        for (int i = 0; i < numStrings; i++)
            s.tag("SimulationElementType", elementTypes.getString(i));

        numStrings = materialTypes.getNumStrings();
        for (int i = 0; i < numStrings; i++)
            s.tag("SimulationMaterialType", materialTypes.getString(i));

        s.tag("Software");
        s.tag("program", "OpenSees");
        //  s.tag("version",version);
        s.endTag(); // Software


        s.tag("ComputerResource");
        s.tag("OS", "Linux");
        s.tag("machine", "local");
        s.endTag(); // Computer

        s.tag("Files");
        PrintFiles(s, theFiles);
        s.endTag(); // Files

        s.endTag(); // SimulationRun
        s.endTag(); // Central
    }

    else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "{\n";
        s << "\"StructuralAnalysisModel\": {\n";

        // BIM model name
        if (title != 0)
            s << "\t\"BIM\": \"" << title << "\",\n";
        else
            s << "\t\"BIM\": \"unknown\",\n";

        // description of model
        if (description != 0)
            s << "\t\"description\": \"" << description << "\",\n";
        else
            s << "\t\"description\": \"\",\n";

        // engineer/modeler contact information
        if (contactName != 0)
            s << "\t\"engineer\": \"" << contactName << "\",\n";
        else
            s << "\t\"engineer\": \"\",\n";

        // units
        s << "\t\"units\": {\n";
        if (forceUnit != 0)
            s << "\t\t\"force\": \"" << forceUnit << "\",\n";
        else
            s << "\t\t\"force\": \"\",\n";
        if (lengthUnit != 0)
            s << "\t\t\"length\": \"" << lengthUnit << "\",\n";
        else
            s << "\t\t\"length\": \"\",\n";
        if (timeUnit != 0)
            s << "\t\t\"time\": \"" << timeUnit << "\",\n";
        else
            s << "\t\t\"time\": \"\",\n";
        if (temperatureUnit != 0)
            s << "\t\t\"temperature\": \"" << temperatureUnit << "\"\n";
        else
            s << "\t\t\"temperature\": \"\"\n";
        s << "\t},\n";

        // model types
        int numStrings = modelTypes.getNumStrings();
        if (numStrings > 0) {
            s << "\t\"modelTypes\": {\n";
            for (int i = 0; i < numStrings-1; i++)
                s << "\t\t\"type\": \"" << modelTypes.getString(i) << "\",\n";
            s << "\t\t\"type\": \"" << modelTypes.getString(numStrings - 1) << "\"\n";
            s << "\t},\n";
        }

        // analysis types
        numStrings = analysisTypes.getNumStrings();
        if (numStrings > 0) {
            s << "\t\"analysisTypes\": {\n";
            for (int i = 0; i < numStrings - 1; i++)
                s << "\t\t\"type\": \"" << analysisTypes.getString(i) << "\",\n";
            s << "\t\t\"type\": \"" << analysisTypes.getString(numStrings - 1) << "\"\n";
            s << "\t},\n";
        }

        // loading types
        numStrings = loadingTypes.getNumStrings();
        if (numStrings > 0) {
            s << "\t\"loadingTypes\": {\n";
            for (int i = 0; i < numStrings - 1; i++)
                s << "\t\t\"type\": \"" << loadingTypes.getString(i) << "\",\n";
            s << "\t\t\"type\": \"" << loadingTypes.getString(numStrings - 1) << "\"\n";
            s << "\t},\n";
        }

        // element types
        numStrings = elementTypes.getNumStrings();
        if (numStrings > 0) {
            s << "\t\"elementTypes\": {\n";
            for (int i = 0; i < numStrings - 1; i++)
                s << "\t\t\"type\": \"" << elementTypes.getString(i) << "\",\n";
            s << "\t\t\"type\": \"" << elementTypes.getString(numStrings - 1) << "\"\n";
            s << "\t},\n";
        }

        // material types
        numStrings = materialTypes.getNumStrings();
        if (numStrings > 0) {
            s << "\t\"materialTypes\": {\n";
            for (int i = 0; i < numStrings - 1; i++)
                s << "\t\t\"type\": \"" << materialTypes.getString(i) << "\",\n";
            s << "\t\t\"type\": \"" << materialTypes.getString(numStrings - 1) << "\"\n";
            s << "\t},\n";
        }
    }
}



OPS_Stream &operator<<(OPS_Stream &s, const SimulationInformation &E)
{
  E.Print(s);
  return s;
}


