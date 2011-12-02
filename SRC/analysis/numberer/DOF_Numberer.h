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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/numberer/DOF_Numberer.h,v $
                                                                        
                                                                        
// File: ~/analysis/numberer/DOF_Numberer.h
// 
// Written: fmk 
// Created: 9/96
// Revision: A
//
// Description: This file contains the class definition for DOF_Numberer.
// DOF_Numberer is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// What: "@(#) DOF_Numberer.h, revA"

#ifndef DOF_Numberer_h
#define DOF_Numberer_h

#include <MovableObject.h>

class AnalysisModel;
class GraphNumberer;
class FEM_ObjectBroker;
class ID;

class DOF_Numberer: public MovableObject
{
  public:
    DOF_Numberer(int classTag);
    DOF_Numberer(GraphNumberer &theGraphNumberer);    
    DOF_Numberer();    
    virtual ~DOF_Numberer();

    virtual void setLinks(AnalysisModel &theModel);
    

    // pure virtual functions
    virtual int numberDOF(int lastDOF_Group = -1);
    virtual int numberDOF(ID &lastDOF_Groups);    

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    

  protected:
    AnalysisModel *getAnalysisModelPtr(void) const;
    GraphNumberer *getGraphNumbererPtr(void) const;
    
  private:
    AnalysisModel *theAnalysisModel;
    GraphNumberer *theGraphNumberer;
};

#endif

