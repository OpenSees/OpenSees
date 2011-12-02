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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-12-18 22:46:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/SimulationInformation.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 04/06
//
// Description: This file contains the class definition for SimulationInformation.
// SimulationInformation is used to store information about a simulation; this
// includes start and end times; program version, files opened for reading, files
// opened for writing, and all parameters used (specified with pset or -par option
// to program)
//
// What: "@(#) SimulationInformation.h, revA"

#ifndef SimulationInformation_h
#define SimulationInformation_h

#include <OPS_Globals.h>
#include <StringContainer.h>
#include <File.h>

class SimulationInformation
{
  public:
  SimulationInformation();    
  ~SimulationInformation();

  int start(void);
  int end(void);

  int setTitle(const char *name);
  int setDescription(const char *name);
  int setContact(const char *name);
  int setLengthUnit(const char *name);
  int setForceUnit(const char *name);
  int setTimeUnit(const char *name);

  int addInputFile(const char *filename, const char *path);
  int addOutputFile(const char *filename, const char *path);
  int addDir(const char *path);

  int addParameter(const char *name, const char *value);
  int addModelType(const char *type);
  int addAnalysisType(const char *type);
  int addLoadingType(const char *type);
  int addElementType(const char *type);
  int addMaterialType(const char *type);

  /*  int addTclInformationCommands(Tcl_Interp *interp); */

  void Print(OPS_Stream &s) const;   
  friend OPS_Stream &operator<<(OPS_Stream &s, const SimulationInformation &E);    
  
  int neesUpload(const char *username, const char *passwd, int projID, int expID);

 protected:
  
 private:
  char *title;
  char *description;
  char *contactName;

  char *lengthUnit;
  char *forceUnit;
  char *timeUnit;

  char startTime[30];
  char endTime[30];

  File *theFiles;
  int numInputFiles;

  StringContainer paramNames;
  StringContainer paramValues;
  StringContainer analysisTypes;
  StringContainer modelTypes;
  StringContainer loadingTypes;
  StringContainer elementTypes;
  StringContainer materialTypes;
};


#endif
