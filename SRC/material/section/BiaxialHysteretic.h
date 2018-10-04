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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
                                                                        
// File: ~/section/BiaxialHysteretic.h
//
// Written: MHS
// Created: Jun 2000
// Revision: A
//
// Description: This file contains the class definition for 
// BiaxialHysteretic.  BiaxialHysteretic decorates an MP
// section (couple bending and axial) with an uncoupled shear
// relation.
//
// What: "@(#) BiaxialHysteretic.h, revA"

#ifndef BiaxialHysteretic_h
#define BiaxialHysteretic_h

#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>

#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <vector>

class BiaxialHysteretic : public SectionForceDeformation
{
public:
    BiaxialHysteretic();
    BiaxialHysteretic(int tag, double k, double fc, double fn,
		      double alp, double als, double eta, double r0,
		      double rp, double rs, double rc, double rn,
		      double Rs, double sig, double lmbda,
		      int code1=SECTION_RESPONSE_MZ, int code2=SECTION_RESPONSE_MY);
    ~BiaxialHysteretic();

    const char *getClassType(void) const {return "BiaxialHysteretic";};

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);
    const Matrix &getInitialFlexibility(void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
 
    SectionForceDeformation *getCopy(void);
    const ID &getType (void);
    int getOrder (void) const;

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    const Vector & getSectionDeformationSensitivity(int gradIndex);
    const Matrix & getSectionTangentSensitivity(int gradIndex);
    const Matrix & getInitialTangentSensitivity(int gradIndex);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const Vector &getdedh(void); // MHS hack
    // AddingSensitivity:END ///////////////////////////////////////////

private:

    // update strenghs and stiffness of 3 springs
    void updateSprings();

    // update loading state of a spring
    int updateLoadingState(int a);

    // update zero force point for a spring
    int updateZeroForcePoint(int a);

    // update total force for a spring
    int updateForce(int a);

    // update tangent stiffness for a spring
    int updateTangent(int a);

    // update energy
    void updateEnergy();

    // sign function
    static int sign(double v);

    // single variable newton iteration
    int newton(double& x, double du1, double F1, 
	       double tol=1e-8, int maxiter = 50);

    // nonlinear spring1 function
    double spring1(double du1, double F1, double F);
    double dspring1(double du1, double F1, double F);
    
private:

    // strength and stiffness of 3 springs
    double Fh, kh, Fp, kp, ku;

    // energy
    double Et, Eh;

    // deterioration
    double r0, rp, rs, rc, rn;

    // strength, stiffness
    double fn, fc, k, alp, als, eta;

    // pinching
    double Rs, sig, lmbda; 

    // zero force point for spring 1 and 2, [spring1, spring2]
    Vector ufx, ufy;

    // previous and current displacements [x, y]
    Vector ui, u;

    // distance to zero force point and total force at step i [spring1, spring2]
    Vector Li, Fi;

    // distance to zero force point and total force at step i+1 [spring1, spring2]
    Vector L, F;

    // Total spring forces in x and y directions [Fx, Fy]
    Vector sF;

    // current du values of spring 1 and 2, [[du0,du1],[du0,du1]]
    std::vector<Vector> du;

    // loading state for spring1 and 2, [spring1 and spring2]
    ID loading, loadingprev;

    // maximum inelastic displacement of the system in + and - directions
    // [positive, negative]
    Vector uxmax, uymax;

    // tangent stiffness matrix
    Matrix Kt;

    static double sqrtpi;
    static double sqrttwo;

    ID code;

    int otherDbTag;

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;

    Vector dedh; // MHS hack
// AddingSensitivity:END ///////////////////////////////////////////

};

#endif
