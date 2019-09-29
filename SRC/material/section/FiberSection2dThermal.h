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

// $Revision: 1.1 $
// $Date: 2011-07-18 10:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FiberSection2dThermal.h,v $

//Modified by Jian Zhang, [University of Edinburgh]
//Modified by Panagiotis Kotsovinos, [University of Edinburgh]

// Description: This file contains the class definition for FiberSection2dThermal
// FiberSection2dThermal provides the abstraction of a 2d beam section discretized by fibers.
// The section stiffness and stress resultants are obtained by summing fiber contributions.
// Also the thermal stress are integrated through section.

#ifndef FiberSection2dThermal_h
#define FiberSection2dThermal_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>


class UniaxialMaterial;
class Fiber;
class Response;
class SectionIntegration;

class FiberSection2dThermal : public SectionForceDeformation
{
  public:
    FiberSection2dThermal();
    FiberSection2dThermal(int tag, int numFibers, Fiber **fibers);
    FiberSection2dThermal(int tag, int num);
    FiberSection2dThermal(int tag, int numFibers, UniaxialMaterial **mats,
			  SectionIntegration &si);
    ~FiberSection2dThermal();

    const char *getClassType(void) const {return "FiberSection2dThermal";};
    //virtual int   setTrialSectionDeformation(const Vector &deforms, const Vector&); //PK changed 18 to 27
    virtual int   setTrialSectionDeformation(const Vector &deforms);
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);
    const Vector &getTemperatureStress(const Vector &tData);

	const Vector&  getThermalElong(void);

    int   commitState(void);
    int   revertToLastCommit(void);
    int   revertToStart(void);

    SectionForceDeformation *getCopy(void);
    const ID &getType (void);
    int getOrder (void) const;

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel,
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc,
			  OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    int addFiber(Fiber &theFiber);
    //double getMaterialPara(void);


    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    const Vector& getStressResultantSensitivity(int gradIndex,
						bool conditional);
    const Vector& getSectionDeformationSensitivity(int gradIndex);
    const Matrix& getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(const Vector& sectionDeformationGradient,
			  int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

    const Vector& determineFiberTemperature(const Vector& , double );  //Added by Liming (UoE)

  protected:

    //  private:
    int numFibers, sizeFibers;       // number of fibers in the section
    UniaxialMaterial **theMaterials; // array of pointers to materials
    double   *matData;               // data for the materials [yloc and area]
    double   kData[4];               // data for ks matrix
    double   sData[2];               // data for s vector

    double QzBar, ABar, yBar;       // Section centroid

    SectionIntegration *sectionIntegr;

    static ID code;

    Vector e;          // trial section deformations
    Vector eCommit;    // committed section deformations
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness
    Vector DataMixed;

    double   sTData[2];   //Data for section resisting force due to thermal load
    Vector  *sT;  //  Pointer to sTData
    double *Fiber_Tangent;
    double *Fiber_ElongP;
    Vector AverageThermalElong;
	//Basiclly this data stores the last committed fiber tangent for calculating thermal foreces
// AddingSensitivity:BEGIN //////////////////////////////////////////
    Vector dedh; // MHS hack
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
