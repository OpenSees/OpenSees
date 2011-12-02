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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/TransientIntegrator.h,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/TransientIntegrator.h
// 
// Written: fmk 
// Created: Tue Sept 17 15:54:47: 1996
// Revision: A
//
// Description: This file contains the class definition for TransientIntegrator.
// TransientIntegrator is an algorithmic class for setting up the finite element
// equations for a static analysis and for Incrementing the nodal displacements
// with the values in the soln vector to the LinearSOE object. 
//
// What: "@(#) TransientIntegrator.h, revA"

#ifndef TransientIntegrator_h
#define TransientIntegrator_h

#include <IncrementalIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class DOF_Group;
class Vector;

class TransientIntegrator : public IncrementalIntegrator
{
  public:
    TransientIntegrator(int classTag);
    virtual ~TransientIntegrator();

    virtual int formTangent(void);    
    virtual int formEleResidual(FE_Element *theEle);
    virtual int formNodUnbalance(DOF_Group *theDof);    

    virtual int initialize(void) {return 0;};
    virtual int newStep(double deltaT) =0;
    
  protected:
    
  private:
};

#endif

