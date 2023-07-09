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
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomainSubIter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/partitioned/PartitionedDomainSubIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// PartitionedDomainSubIter. PartitionedDomainSubIter is an iter for returning
// the Subdomains of an object of class PartitionedDomain. 
//
// What: "@(#) PartitionedDomainSubIter.h, revA"

#ifndef PartitionedDomainSubIter_h
#define PartitionedDomainSubIter_h

#include <SubdomainIter.h>

class TaggedObjectStorage;
class TaggedObjectIter;

class PartitionedDomainSubIter: public SubdomainIter
{
  public:
    PartitionedDomainSubIter(TaggedObjectStorage *theStorage);
    virtual ~PartitionedDomainSubIter();

    virtual void reset(void);
    virtual Subdomain *operator()(void);
    
  private:
    TaggedObjectIter &myIter;
};

#endif





