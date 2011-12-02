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
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/ModelBuilder.h,v $
                                                                        
                                                                        
// File: ~/modelbuilder/ModelBuilder.h
// 
// Written: fmk 
// Created: Mon Sept 15 14:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for ModelBuilder.
// ModelBuilder is an abstract base class, i.e. no objects of it's
// type can be created. A ModelBuilder creates the discritization of 
// the structure.
//
// What: "@(#) ModelBuilder.h, revA"

#ifndef ModelBuilder_h
#define ModelBuilder_h

class Domain;

// #include <Domain.h>

class ModelBuilder
{
  public:
    ModelBuilder(Domain &theDomain);
    virtual ~ModelBuilder();
    
    virtual int buildFE_Model(void) = 0;
    
  protected:
    Domain *getDomainPtr(void) const;
    
  private:
    Domain *myDomain;
};

#endif

