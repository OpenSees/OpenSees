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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/CoupledSection3D.h,v $
                                                                        
                                                                        
// File: ~/section/CoupledSection3D.h
//
// Written: MHS
// Created: Jun 2000
// Revision: A
//
// Description: This file contains the class definition for 
// CoupledSection3D.  CoupledSection3D decorates an MP
// section (couple bending and axial) with an uncoupled shear
// relation.
//
// What: "@(#) CoupledSection3D.h, revA"

#ifndef CoupledSection3D_h
#define CoupledSection3D_h

#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>

#include <Vector.h>
#include <Matrix.h>
#include <string>
#include <Parameter.h>

#define STRENGTH_DOMAIN_UNDEFINED 0
#define STRENGTH_DOMAIN_SIMPLIFIED 1
#define STRENGTH_DOMAIN_BY_POINTS 2


// Massimo: una cosa del genere dicevi? Dopo casomai lo implementiamo così.

class DomainData {

public:
    DomainData();
    DomainData(int nN, int nTheta, int nData);
    ~DomainData();

    double getValue(int i, int j, int k);

    void setValue(int i, int j, int k, double val);

    void print(void);

    DomainData* getCopy(void);

    int getMyMzForNAndDirection(double N, double theta, double& My, double& Mz);

    void getRangeN(double& Nmin, double& Nmax);

    static const double pi;

private:
    int size;
    int numberAxial;
    int numberTheta;
    int numberData;
    Vector* theVector;

};

class CoupledSection3D : public SectionForceDeformation
{
public:
    CoupledSection3D(); 

        CoupledSection3D(int tag, UniaxialMaterial* theTorsionMaterial, UniaxialMaterial* theAxialMaterial, UniaxialMaterial* theShearYMaterial, UniaxialMaterial* theShearZMaterial, UniaxialMaterial* theMomentYMaterial, UniaxialMaterial* theMomentZMaterial,
            DomainData* ultDomain, std::string theRawInitialStiffnessExpressionY, std::string theRawInitialStiffnessExpressionZ, std::string theRawThetaPExpressionY, std::string theRawThetaPExpressionZ,
            std::string theRawThetaPCExpressionY, std::string theRawThetaPCExpressionZ, double a_s_i);
    ~CoupledSection3D();

    const char *getClassType(void) const {return "CoupledSection3D";};

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

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	//by SAJalali
	int getResponse(int responseID, Information &info);

    void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *, Information &);

    //// AddingSensitivity:BEGIN //////////////////////////////////////////
    //int setParameter(const char **argv, int argc, Parameter &param);
    //const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    //const Vector & getSectionDeformationSensitivity(int gradIndex);
    //const Matrix & getSectionTangentSensitivity(int gradIndex);
    //const Matrix & getInitialTangentSensitivity(int gradIndex);
    //int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    //const Vector &getdedh(void); // MHS hack
    //// AddingSensitivity:END ///////////////////////////////////////////

protected:
    
private:

    void resetStrengthDomain(double& My_u_p, double& My_u_n, double& Mz_u_p, double& Mz_u_n);
    void setUncoupledStrengthDomainforAxial(const double N, double& My_u_p, double& My_u_n, double& Mz_u_p, double& Mz_u_n);

    UniaxialMaterial* axialMaterial;  // Material pointer to uniaxial spring (elastic for axial load)
    UniaxialMaterial* MyMaterial;   // Material pointer to uniaxial spring (for now Pinching4 for My)
    UniaxialMaterial* MzMaterial;   // Material pointer to uniaxial spring (for now Pinching4 for Mz)
    UniaxialMaterial* torsionMaterial; // Material pointer to uniaxial spring (elastic for torion)
    UniaxialMaterial* VyMaterial; // Material pointer to uniaxial spring (elastic for shear Vy load)
    UniaxialMaterial* VzMaterial; // Material pointer to uniaxial spring (elastic for shear Vz load)

    double tolN = 1e-5;
    double tolM = 1e-5;

    ID *matCodes;
    int numMats;
    
    Vector *e;    // Storage for section deformations
    Vector *s;    // Storage for stress resultants
    Matrix *ks;   // Storage for section stiffness
    Matrix *fs;   // Storage for section flexibility
    ID     *theCode;     // Storage for section type information

    DomainData* ultimateDomain;

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
   
    int otherDbTag;

    static double workArea[];
    static int codeArea[];

//// AddingSensitivity:BEGIN //////////////////////////////////////////
//    Vector dedh; // MHS hack
//// AddingSensitivity:END ///////////////////////////////////////////

};

#endif
