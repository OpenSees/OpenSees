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
                                                                        
// $Revision: 1.1 $
// $Date: 2008-12-09 21:24:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclSectionTester.h,v $

// Written: fmk 
// Created: 10/08
// Revision: A

#ifndef TclSectionTester_h
#define TclSectionTester_h

#include <TclModelBuilder.h>

class SectionForceDeformation;

#include <tcl.h>

class TclSectionTester : public TclModelBuilder
{
  public:
  TclSectionTester(Domain &theDomain,Tcl_Interp *interp, int count=1);
  ~TclSectionTester();    
  
  protected:
  
 private:
  Tcl_Interp *theInterp;
};

#endif







