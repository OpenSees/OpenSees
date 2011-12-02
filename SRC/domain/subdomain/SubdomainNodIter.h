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
// $Date: 2009-08-25 22:09:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/SubdomainNodIter.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// SubdomainNodIter. SubdomainNodIter is an iter 
// for returning the Nodes of an object of class
// Subdomain.
//
// What: "@(#) SubdomainNodIter.h, revA"


#ifndef SubdomainNodIter_h
#define SubdomainNodIter_h

#include <NodeIter.h>

class Subdomain;
class ArrayOfTaggedObjectsIter;
class Subdomain;

class SubdomainNodIter: public NodeIter
{
  public:
    SubdomainNodIter(Subdomain &theSubdomain);
    virtual ~SubdomainNodIter();

    virtual void reset(void);
    virtual Node *operator()(void);
    
  private:
    NodeIter   *currentIter;
    Subdomain  *theSubdomain;    
    bool 	external;
};

#endif









