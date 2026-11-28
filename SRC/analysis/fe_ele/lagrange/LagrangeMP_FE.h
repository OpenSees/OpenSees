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
// $Date: 2006-02-08 20:20:00 $
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

class LagrangeMP_FE final : public FE_Element
{
public:
    LagrangeMP_FE(int tag, Domain &theDomain, MP_Constraint &theMP,
                  DOF_Group &theDofGrp, double alpha = 1.0);
    ~LagrangeMP_FE() override;

    int setID(void) override;
    const Matrix &getTangent(Integrator *theIntegrator) override;
    const Vector &getResidual(Integrator *theIntegrator) override;
    const Vector &getTangForce(const Vector &x, double fact = 1.0) override;

    const Vector &getK_Force(const Vector &x, double fact = 1.0) override;
    const Vector &getKi_Force(const Vector &x, double fact = 1.0) override;
    const Vector &getC_Force(const Vector &x, double fact = 1.0) override;
    const Vector &getM_Force(const Vector &x, double fact = 1.0) override;

    void zeroTangent(void) override;
    void addKtToTang(double fact = 1.0) override;
    void addKiToTang(double fact = 1.0) override;
    void addCtoTang(double fact = 1.0) override;
    void addMtoTang(double fact = 1.0) override;

    void zeroResidual(void) override;
    void addRtoResidual(double fact = 1.0) override;
    void addRIncInertiaToResidual(double fact = 1.0) override;
    void addM_Force(const Vector &accel, double fact = 1.0) override;
    void addD_Force(const Vector &vel, double fact = 1.0) override;

private:
    void determineTangent(void);

    const Matrix &getStaticTangent(void)
    {
        if (timeVarying)
            this->determineTangent();
        return *tang;
    }

    MP_Constraint *theMP;
    Node *theConstrainedNode;
    Node *theRetainedNode;

    DOF_Group *theDofGroup;
    Matrix *tang;     // unscaled static coupling
    Matrix *sysTang;  // integrator-assembled system tangent
    Vector *resid;
    const bool timeVarying;
    bool urLoaded;  // true once upper-right C^T has been copied into sysTang
    double alpha;
    int numU;  // number of displacement dofs before lambda block
};

#endif
