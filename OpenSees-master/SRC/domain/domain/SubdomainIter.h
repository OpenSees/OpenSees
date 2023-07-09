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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/SubdomainIter.h,v $
                                                                        
                                                                        
// File: ~/domain/domain/SubdomainIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for SubdomainIter.
// SubdomainIter is an abstract base class. A  SubdomainIter object is an iter 
// for returning the subdomains of an object of class PartitionedDomain. SubdomainIters 
// must be written for each subclass of PartitionedDomain, hence the abstract base class.


#ifndef SubdomainIter_h
#define SubdomainIter_h

class Subdomain;

class SubdomainIter
{
  public:
    SubdomainIter() {};
    virtual ~SubdomainIter() {};

    virtual Subdomain *operator()(void) =0;
    
  protected:

  private:    

};

#endif

