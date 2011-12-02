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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-04-28 17:54:28 $
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

class SimulationInformation
{
  public:
  SimulationInformation();    
  ~SimulationInformation();
  int start(void);
  int end(void);
  int addReadFile(const char *);
  int addWriteFile(const char *);
  int addParameter(const char *name, const char *value);

  
  void Print(OPS_Stream &s) const;   
  friend OPS_Stream &operator<<(OPS_Stream &s, const SimulationInformation &E);    
  
 protected:
  
 private:
  char startTime[30];
  char endTime[30];
  char **filesRead;
  char **filesWritten;
  char **paramNames;
  char **paramValues;
  int numFilesWritten;
  int numFilesRead;
  int numParameters;
};


#endif
