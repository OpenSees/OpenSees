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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-11-19 23:38:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/penalty/PenaltySP_FE.h,v $
                                                                        
                                                                        
#ifndef PenaltySP_FE_h
#define PenaltySP_FE_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for PenaltySP_FE.
// PenaltySP_FE is a subclass of FE_Element which handles SP_Constraints
// using the penalty method.
//
// What: "@(#) PenaltySP_FE.h, revA"

#include <FE_Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class Element;
class Integrator;
class AnalysisModel;
class Domain;
class SP_Constraint;
class Node;

class PenaltySP_FE: public FE_Element
{
  public:
    PenaltySP_FE(int tag, Domain &theDomain, SP_Constraint &theSP, double alpha=1.0e8);    
    virtual ~PenaltySP_FE();    

    // public methods
    virtual int  setID(void);
    virtual const Matrix &getTangent(Integrator *theIntegrator);
    virtual const Vector &getResidual(Integrator *theIntegrator);
    virtual const Vector &getTangForce(const Vector &x, double fact = 1.0);

    virtual const Vector &getK_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getKi_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getC_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getM_Force(const Vector &x, double fact = 1.0);    

  protected:
    
  private:
    double alpha;
    SP_Constraint *theSP;
    Node *theNode;
    static Matrix tang;
    static Vector resid;
};

#endif


