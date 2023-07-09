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
// $Date: 2009-08-25 22:10:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomainEleIter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/partitioned/PartitionedDomainEleIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// PartitionedDomainEleIter. PartitionedDomainEleIter is an iter 
// for returning the elements of an object of class
// PartitionedDomain.
//
// What: "@(#) PartitionedDomainEleIter.h, revA"


#ifndef PartitionedDomainEleIter_h
#define PartitionedDomainEleIter_h

#include <ElementIter.h>
#include <SingleDomEleIter.h>


class PartitionedDomain;
class ArrayOfTaggedObjectsIter;
class Subdomain;
class Domain;

class PartitionedDomainEleIter: public ElementIter
{
  public:
    PartitionedDomainEleIter(PartitionedDomain *theDomain);
    virtual ~PartitionedDomainEleIter();

    virtual void reset(void);
    virtual Element *operator()(void);
    
  private:
    SingleDomEleIter	        *mainEleIter;
    ArrayOfTaggedObjectsIter 	*subdomainIter;
    ElementIter 	     	*currentIter;
    Subdomain 	 	     	*currentSubdomain;    
    PartitionedDomain	        *thePartitionedDomain;
    bool			mainDomain;
};

#endif









