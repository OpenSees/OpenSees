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
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/StaticIntegrator.h,v $
                                                                        
                                                                        
#ifndef StaticIntegrator_h
#define StaticIntegrator_h

// File: ~/analysis/integrator/StaticIntegrator.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for StaticIntegrator.
// StaticIntegrator is an algorithmic class for setting up the finite element
// equations for a static analysis and for Incrementing the nodal displacements
// with the values in the soln vector to the LinearSOE object. 
//
// What: "@(#) StaticIntegrator.h, revA"

#include <IncrementalIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;


class StaticIntegrator : public IncrementalIntegrator
{
  public:
    StaticIntegrator(int classTag);    

    virtual ~StaticIntegrator();

    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    virtual int formEleTangent(FE_Element *theEle);
    virtual int formEleResidual(FE_Element *theEle);
    virtual int formNodTangent(DOF_Group *theDof);        
    virtual int formNodUnbalance(DOF_Group *theDof);    
   virtual int formEleTangentSensitivity(FE_Element *theEle,int gradNumber); 
   
   virtual int newStep(void) =0;    

  protected:

 
  private:
};

#endif

