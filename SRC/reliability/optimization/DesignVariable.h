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
**   Optimization module developed by:                                   **
**   Quan Gu  (qgu@ucsd.edu)                          **
**   Joel Conte (jpconte@ucsd.edu)                       **
**   Philip Gill (pgill@ucsd.edu)                                                                 **
** ****************************************************************** */
                                                                        

//
// Written by Quan Gu (qgu@ucsd.edu)

#ifndef DesignVariable_h
#define DesignVariable_h

#include <ReliabilityDomainComponent.h>
#include <iostream>
#include <string.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <ReliabilityDomain.h>
#include <tcl.h>

class DesignVariable : public ReliabilityDomainComponent  
{

public:
  DesignVariable(int passedTag, 
		 char *passedName,
		 double passedValue,
		 double passedUpperBound,
		 double passedLowerBound,
		 Tcl_Interp *passedTclInterp, 
		 ReliabilityDomain *passedReliabilityDomain,
		 double passedXMultiplier);
  ~DesignVariable();
  
  //---- useful interfaces -------------------
  
  double getUpperBound();
  double getLowerBound();
  char * getName();
  double getValue();
  double getOldValue();
  double getXMultiplier();
  
  int setName(char * newName);
  int setXMultiplier(double newXMultiplier);
  int setValue(double newValue);
  int setOldValue(double newValue);
  int setUpperBound(double newBound);
  int setLowerBound(double newBound);
  int update(double newX);
  
  // --- deal with DVPositioner --
  int addDVPositioner(int tag);
  int getNumOFDVPositioners(){return numOfMyDVPositioners;};
  int getDVPositionerTag(int i);
  
  double getScale(){return scale;};
  
  void Print(OPS_Stream &s, int flag =0){return;};  
  
 private:
  Tcl_Interp *theTclInterp;
  ReliabilityDomain *theReliabilityDomain;
  double upperBound;
  double lowerBound;
  char name[20];
  
  
  double value;
  double oldValue; // save for xstate
  double xMultiplier; // for xmul
  
  
  // --- deal with DVPositioner --
  int * myDVPositionerTags;
  int numOfMyDVPositioners;
  int maxNumOfMyDVPositioners;
  double scale;
};

#endif
