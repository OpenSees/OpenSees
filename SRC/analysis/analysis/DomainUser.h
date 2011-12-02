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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/DomainUser.h,v $
                                                                        
                                                                        
// File: ~/analysis/analysis/DomainUser.h
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for DomainUser.
// DomainUser is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// What: "@(#) DomainUser.h, revA"

#ifndef DomainUser_h
#define DomainUser_h

class Domain;

class DomainUser
{
  public:
    DomainUser(Domain &theDomain);
    virtual ~DomainUser();

    // pure virtual functions
    virtual void domainChanged(void) =0;
    
  protected:
    Domain &getDomain(void) const;
    
  private:
    Domain &myDomain;
};

#endif

