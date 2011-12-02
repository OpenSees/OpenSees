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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/lagrange/LagrangeMP_FE.h,v $
                                                                        
                                                                        
// File: ~/analysis/fe_ele/lagrange/LagrangeMP_FE.h
// 
// Written: fmk 
// Created: 02/99
// Revision: A
//
// Description: This file contains the class definition for LagrangeMP_FE.
// LagrangeMP_FE is a subclass of FE_Element which handles MP_Constraints
// using the Lagrange method.
//
// What: "@(#) LagrangeMP_FE.h, revA"


#ifndef LagrangeMP_FE_h
#define LagrangeMP_FE_h

#include <FE_Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class Element;
class Integrator;
class AnalysisModel;
class Domain;
class MP_Constraint;
class Node;
class DOF_Group;

class LagrangeMP_FE: public FE_Element
{
  public:
    LagrangeMP_FE(Domain &theDomain, MP_Constraint &theMP, 
		  DOF_Group &theDofGrp, double alpha = 1.0);
    virtual ~LagrangeMP_FE();    

    // public methods
    virtual int  setID(void);
    virtual const Matrix &getTangent(Integrator *theIntegrator);    
    virtual const Vector &getResidual(Integrator *theIntegrator);    
    virtual const Vector &getTangForce(const Vector &x, double fact = 1.0);
    
  protected:
    
  private:
    double alpha;
    void determineTangent(void);
    
    MP_Constraint *theMP;
    Node *theConstrainedNode;
    Node *theRetainedNode;    

    DOF_Group *theDofGroup;
    Matrix *tang;
    Vector *resid;
};

#endif


