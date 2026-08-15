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

class PenaltySP_FE final : public FE_Element
{
  public:
    PenaltySP_FE(int tag, Domain &theDomain, SP_Constraint &theSP, double alpha=1.0e8);    
    ~PenaltySP_FE() override;

    // public methods
    int  setID(void) override;
    const Matrix &getTangent(Integrator *theIntegrator) override;
    const Vector &getResidual(Integrator *theIntegrator) override;
    const Vector &getTangForce(const Vector &x, double fact = 1.0) override;

    const Vector &getK_Force(const Vector &x, double fact = 1.0) override;
    const Vector &getKi_Force(const Vector &x, double fact = 1.0) override;
    const Vector &getC_Force(const Vector &x, double fact = 1.0) override;
    const Vector &getM_Force(const Vector &x, double fact = 1.0) override;

    // methods to allow integrator to build tangent
    void  zeroTangent(void) override;
    void  addKtToTang(double fact = 1.0) override;
    void  addKiToTang(double fact = 1.0) override;
    void  addCtoTang(double fact = 1.0) override;
    void  addMtoTang(double fact = 1.0) override;

    // methods to allow integrator to build residual
    void  zeroResidual(void) override;
    void  addRtoResidual(double fact = 1.0) override;
    void  addRIncInertiaToResidual(double fact = 1.0) override;
    void  addM_Force(const Vector &accel, double fact = 1.0) override;
    void  addD_Force(const Vector &vel, double fact = 1.0) override;

  private:
    double alpha;
    SP_Constraint *theSP;
    Node *theNode;
    static Matrix tang;
    static Vector resid;
};

#endif
