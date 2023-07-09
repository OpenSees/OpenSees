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
                                                                        
// $Revision: 1.4 $
// $Date: 2007-06-26 20:13:24 $
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

#ifndef TclPlaneStressMaterialTester_h
#define TclPlaneStressMaterialTester_h

#include <TclModelBuilder.h>
#include <tcl.h>

class TclPlaneStressMaterialTester : public TclModelBuilder
{
  public:
    TclPlaneStressMaterialTester(Domain &theDomain,Tcl_Interp *interp, int count=1);
    ~TclPlaneStressMaterialTester();    

  protected:

  private:
    Tcl_Interp *theInterp;
};

#endif







