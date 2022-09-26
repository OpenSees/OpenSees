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
                                                                        
// $Revision: 1.17 $
// $Date: 2009-10-01 23:04:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ASDCoupledHinge3D.h,v $
                                                                        
                                                                        
// File: ~/section/ASDCoupledHinge3D.h
//
// Written: Diego Talledo
// Created: Jun 2021
// Revision: A
//
// Description: This file contains the class definition for 
// ASDCoupledHinge3D.  ASDCoupledHinge3D decorates a PMM
// section (couple bending and axial load) with an uncoupled shear
// relation.
//
// What: "@(#) ASDCoupledHinge3D.h, revA"

#ifndef ASDCoupledHinge3D_h
#define ASDCoupledHinge3D_h

#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>

#include <Vector.h>
#include <Matrix.h>
#include <string>
#include <Parameter.h>

#define STRENGTH_DOMAIN_UNDEFINED 0
#define STRENGTH_DOMAIN_SIMPLIFIED 1
#define STRENGTH_DOMAIN_BY_POINTS 2

#define ASD_HINGE_NUM_TANG

class ASDCoupledHinge3D;
class ASDCoupledHinge3DDomainData {
    friend class ASDCoupledHinge3D;
public:
    ASDCoupledHinge3DDomainData() = default;
    ASDCoupledHinge3DDomainData(int nN, int nTheta, int nData);

    double getValue(int i, int j, int k);

    void setValue(int i, int j, int k, double val);

    void print(void);

    int getMyMzForNAndDirection(double N, double theta, double& My, double& Mz);

    void getRangeN(double& Nmin, double& Nmax);

private:
    int size = 0;
    int numberAxial = 0;
    int numberTheta = 0;
    int numberData = 0;
    Vector theVector;

};

class ASDCoupledHinge3D : public SectionForceDeformation
{
public:
    ASDCoupledHinge3D(); 

    ASDCoupledHinge3D(int tag, 
        UniaxialMaterial* theTorsionMaterial, 
        UniaxialMaterial* theAxialMaterial,
        UniaxialMaterial* theShearYMaterial,
        UniaxialMaterial* theShearZMaterial, 
        UniaxialMaterial* theMomentYMaterial, 
        UniaxialMaterial* theMomentZMaterial,
        const ASDCoupledHinge3DDomainData &theStrengthDomain,
        const std::string &theRawInitialStiffnessExpressionY, 
        const std::string &theRawInitialStiffnessExpressionZ, 
        const std::string &theRawThetaPExpressionY, 
        const std::string &theRawThetaPExpressionZ,
        const std::string &theRawThetaPCExpressionY, 
        const std::string &theRawThetaPCExpressionZ, 
        double a_s_i);
    ~ASDCoupledHinge3D();

    const char *getClassType(void) const {return "ASDCoupledHinge3D";};

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);
    const Matrix &getSectionFlexibility(void);
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

    int getVariable(const char *, Information &);

    virtual Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    virtual int getResponse(int responseID, Information& info);

private:

    void updateLaws(void);
    void resetStrengthDomain(double& My_u_p, double& My_u_n, double& Mz_u_p, double& Mz_u_n);
    void setUncoupledStrengthDomainforAxial(const double N, double& My_u_p, double& My_u_n, double& Mz_u_p, double& Mz_u_n);
    void setupParameters();

    // materials
    UniaxialMaterial* axialMaterial;  // Material pointer to uniaxial spring (elastic for axial load)
    UniaxialMaterial* MyMaterial;   // Material pointer to uniaxial spring (for now Pinching4 for My)
    UniaxialMaterial* MzMaterial;   // Material pointer to uniaxial spring (for now Pinching4 for Mz)
    UniaxialMaterial* torsionMaterial; // Material pointer to uniaxial spring (elastic for torion)
    UniaxialMaterial* VyMaterial; // Material pointer to uniaxial spring (elastic for shear Vy load)
    UniaxialMaterial* VzMaterial; // Material pointer to uniaxial spring (elastic for shear Vz load)

    // settings
    double tolN = 1e-5;
    double tolM = 1e-5;
    double MmaxAbs = 0.0;
    double errExplicit = 0.0;

    ASDCoupledHinge3DDomainData strengthDomain; // strength domain is point My. Mmax = My * a_s

    // strings
    std::string rawInitialStiffnessExpressionY;
    std::string rawInitialStiffnessExpressionZ;
    std::string rawThetaPExpressionY;
    std::string rawThetaPExpressionZ;
    std::string rawThetaPCExpressionY; 
    std::string rawThetaPCExpressionZ;
    double a_s = 1.001;

    // parameters for updating backbone during analysis
    // My
    // positive
    Parameter par_f1p_Y;
    Parameter par_d1p_Y;
    Parameter par_f2p_Y;
    Parameter par_d2p_Y;
    Parameter par_f3p_Y;
    Parameter par_d3p_Y;
    Parameter par_f4p_Y;
    Parameter par_d4p_Y;
    // negative
    Parameter par_f1n_Y;
    Parameter par_d1n_Y;
    Parameter par_f2n_Y;
    Parameter par_d2n_Y;
    Parameter par_f3n_Y;
    Parameter par_d3n_Y;
    Parameter par_f4n_Y;
    Parameter par_d4n_Y;
    // Mz
    // positive
    Parameter par_f1p_Z;
    Parameter par_d1p_Z;
    Parameter par_f2p_Z;
    Parameter par_d2p_Z;
    Parameter par_f3p_Z;
    Parameter par_d3p_Z;
    Parameter par_f4p_Z;
    Parameter par_d4p_Z;
    // negative
    Parameter par_f1n_Z;
    Parameter par_d1n_Z;
    Parameter par_f2n_Z;
    Parameter par_d2n_Z;
    Parameter par_f3n_Z;
    Parameter par_d3n_Z;
    Parameter par_f4n_Z;
    Parameter par_d4n_Z;

    // working memory
   
    int otherDbTag = 0;

#ifdef ASD_HINGE_NUM_TANG
    Vector m_num_tang = Vector(6);
#endif // ASD_HINGE_NUM_TANG


};

#endif
