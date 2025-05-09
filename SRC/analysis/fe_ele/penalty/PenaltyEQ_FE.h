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
                                                                        
// $Revision: 1.0 $
// $Date: 2025-05-09$
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/penalty/PenaltyEQ_FE.h,v $
                                                                        
                                                                        
#ifndef PenaltyEQ_FE_h
#define PenaltyEQ_FE_h

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 05/2020
// Revision: A
//
// Description: This file contains the class definition for PenaltyEQ_FE.
// PenaltyEQ_FE is a subclass of FE_Element which handles EQ_Constraints
// using the penalty method.
//
// What: "@(#) PenaltyEQ_FE.h, revA"

#include <FE_Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class Element;
class Integrator;
class AnalysisModel;
class Domain;
class EQ_Constraint;
class Node;

class PenaltyEQ_FE: public FE_Element
{
  public:
    PenaltyEQ_FE(int tag, Domain &theDomain, EQ_Constraint &theEQ, double alpha);
    virtual ~PenaltyEQ_FE();    

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
    void determineTangent(void);
    
    EQ_Constraint *theEQ;
    Node *theConstrainedNode;
    Node **theRetainedNode;    

    Matrix *tang;
    Vector *resid;
    Vector *C;
    double alpha;
	
    
};

#endif


