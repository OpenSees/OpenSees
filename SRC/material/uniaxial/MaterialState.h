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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: August 2000
//
// Description: This file contains the interface for MaterialState,
// which is the base for classes which model states of hysteretic
// degradation.

#ifndef MaterialState_h
#define MaterialState_h

#include <DomainComponent.h>
#include <MovableObject.h>

class MaterialState : public TaggedObject, public MovableObject
{
 public:
  MaterialState(int tag, int classTag);    
  virtual ~MaterialState();

  virtual int setVariable(const char *argv);
  virtual int getVariable(int variableID, double &info);
  
  virtual int setParameter(const char **argv, int argc,
			   Information &eleInformation);
  virtual int updateParameter(int responseID, Information &eleInformation);  
 protected:
  
 private:

};

#endif
