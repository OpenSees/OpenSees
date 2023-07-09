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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:01:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/PartitionedQuick2dFrame.h,v $
                                                                        
                                                                        
// File: ~/modelbuilder/PartitionedQuick2dFrame.h
// 
// Written: fmk 
// Created: Mon Sept 15 14:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for PartitionedQuick2dFrame.
// A PartitionedQuick2dFrame creates the plane frame models.

//
// What: "@(#) ModelBuilder.h, revA"

#ifndef PartitionedQuick2dFrame_h
#define PartitionedQuick2dFrame_h

#include <PartitionedModelBuilder.h>
#include <iOPS_Stream.h>

class Element;
class Node;
class SP_Constraint;
class MP_Constraint;
class Subdomain;

class PartitionedQuick2dFrame : public PartitionedModelBuilder
{
  public:
    PartitionedQuick2dFrame(PartitionedDomain &theDomain, int numX, int numY, int eleType);
    PartitionedQuick2dFrame(Subdomain &theSubdomain); // for actorsubdomain
    ~PartitionedQuick2dFrame();    

    int buildInterface(int numSubdomains);
    int buildSubdomain(int partition, int numPartitions, Subdomain &theSubdomain);
    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    
  private:
    int numX, numY, eleType;
    int numEle, numNode, numSPs, numMPs, numLCs;
};

#endif
