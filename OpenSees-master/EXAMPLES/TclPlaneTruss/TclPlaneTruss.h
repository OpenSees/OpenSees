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
// $Date: 2003-02-19 21:49:00 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/TclPlaneTruss/TclPlaneTruss.h,v $

                                                                        
// File: ~/modelbuilder/TclPlaneTruss.h
// 
// Written: fmk 
// Created: Mon Sept 15 14:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for TclPlaneTruss.
// A TclPlaneTruss creates the plane frame models.

//
// What: "@(#) ModelBuilder.h, revA"

#ifndef TclPlaneTruss_h
#define TclPlaneTruss_h

#include <ModelBuilder.h>

class Element;
class Node;
class SP_Constraint;
class MP_Constraint;
class UniaxialMaterial;
class TaggedObjectStorage;

#include <tcl.h>
#include <tk.h>

class TclPlaneTruss : public ModelBuilder
{
  public:
    TclPlaneTruss(Domain &theDomain,Tcl_Interp *interp);
    ~TclPlaneTruss();    

    int buildFE_Model(void);
    
    int addMaterial(UniaxialMaterial &theMaterial);
    UniaxialMaterial *getMaterial(int tag);
    Domain *getDomain;

  protected:

  private:
    TaggedObjectStorage *theMaterials;
};

#endif







