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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/lagrange/LagrangeEQ_FE.h,v $
                                                                        
                                                                        
// File: ~/analysis/fe_ele/lagrange/LagrangeEQ_FE.h
// 
// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 05/2020
// Revision: A
//
// Description: This file contains the class definition for LagrangeEQ_FE.
// LagrangeEQ_FE is a subclass of FE_Element which handles EQ_Constraints
// using the Lagrange method.
//
// What: "@(#) LagrangeEQ_FE.h, revA"


#ifndef LagrangeEQ_FE_h
#define LagrangeEQ_FE_h

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
class DOF_Group;

class LagrangeEQ_FE: public FE_Element
{
  public:
    LagrangeEQ_FE(int tag, Domain &theDomain, EQ_Constraint &theEQ, 
		  DOF_Group &theDofGrp, double alpha = 1.0);
    virtual ~LagrangeEQ_FE();    

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
    void determineTangent(void);
    
    EQ_Constraint *theEQ;
    Node *theConstrainedNode;
    Node **theRetainedNode;    

    DOF_Group *theDofGroup;
    Matrix *tang;
    Vector *resid;
};

#endif


