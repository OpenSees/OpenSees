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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-07-24 18:28:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclUniaxialMaterialTester.h,v $
                                                                        
// File: ~/modelbuilder/tcl/TclUniaxialMaterialTester.h
// 
// Written: fmk 
// Created: 03/01
// Revision: A
//
// Description: This file contains the class definition for TclUniaxialMaterialTester.
// A TclUniaxialMaterialTester adds the commands to create and test uniaxial materials
//
// What: "@(#) TclUniaxialMaterialTester.h, revA"

#ifndef TclUniaxialMaterialTester_h
#define TclUniaxialMaterialTester_h

#include <TclModelBuilder.h>

class SectionForceDeformation;
class SectionRepres;
class UniaxialMaterial;
class NDMaterial;
class TaggedObjectStorage;

class CrdTransf2d;
class CrdTransf3d;

#include <tcl.h>
#include <tk.h>

class TclUniaxialMaterialTester : public TclModelBuilder
{
  public:
    TclUniaxialMaterialTester(Domain &theDomain,Tcl_Interp *interp, int count=1);
    ~TclUniaxialMaterialTester();    

  protected:

  private:
    Tcl_Interp *theInterp;
};

#endif







