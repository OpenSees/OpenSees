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
// $Date: 2004-10-12 21:55:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/lagrange/LagrangeSP_FE.h,v $
                                                                        
                                                                        
#ifndef LagrangeSP_FE_h
#define LagrangeSP_FE_h

// Written: fmk 
// Created: 02/99
// Revision: A
//
// Description: This file contains the class definition for LagrangeSP_FE.
// LagrangeSP_FE is a subclass of FE_Element which handles SP_Constraints
// using the Lagrange method.
//
// What: "@(#) LagrangeSP_FE.h, revA"

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
class DOF_Group;

class LagrangeSP_FE: public FE_Element
{
  public:
    LagrangeSP_FE(Domain &theDomain, SP_Constraint &theSP, 
		  DOF_Group &theDofGrp, double alpha = 1.0);
    virtual ~LagrangeSP_FE();    

    // public methods
    virtual int  setID(void);
    virtual const Matrix &getTangent(Integrator *theIntegrator);    
    virtual const Vector &getResidual(Integrator *theIntegrator);    
    virtual const Vector &getTangForce(const Vector &x, double fact = 1.0);

    virtual const Vector &getK_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getC_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getM_Force(const Vector &x, double fact = 1.0);        
  protected:
    
  private:
    double alpha;
    Matrix *tang;
    Vector *resid;
    
    SP_Constraint *theSP;
    Node *theNode;    

    DOF_Group *theDofGroup;
};

#endif


