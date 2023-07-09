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

// Written: MHS
// Created: 2012
//
// Description: This file contains the class definition for 
// NDFiberSectionWarping2d.h. NDFiberSectionWarping2d provides the abstraction of a 
// 2d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef NDFiberSectionWarping2d_h
#define NDFiberSectionWarping2d_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class NDMaterial;
class Fiber;
class Response;
class SectionIntegration;

class NDFiberSectionWarping2d : public SectionForceDeformation
{
  public:
    NDFiberSectionWarping2d(); 
    NDFiberSectionWarping2d(int tag, int numFibers, Fiber **fibers, double a = 1.0);
  NDFiberSectionWarping2d(int tag, int num, double a = 1.0);
  NDFiberSectionWarping2d(int tag, int numFibers, NDMaterial **mats,
			  SectionIntegration &si, double a = 1.0);
    ~NDFiberSectionWarping2d();

    const char *getClassType(void) const {return "NDFiberSectionWarping2d";};

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);
 
    SectionForceDeformation *getCopy(void);
    const ID &getType (void);
    int getOrder (void) const;
	double getRho(void);

    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
	    
    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    int addFiber(Fiber &theFiber);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);
    const Vector& getStressResultantSensitivity(int gradIndex,
						bool conditional);
    const Vector& getSectionDeformationSensitivity(int gradIndex);
    const Matrix& getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(const Vector& sectionDeformationGradient,
			  int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
    //  private:
  int numFibers,sizeFibers;                   // number of fibers in the section
    NDMaterial **theMaterials; // array of pointers to materials
    double   *matData;               // data for the materials [yloc and area]
    double   kData[25];               // data for ks matrix 
    double   sData[5];               // data for s vector 
    
    double yBar;       // Section centroid
    double alpha;      // Shear shape factor
	
	double yBarZero;
	double DeltaYbar;

    SectionIntegration *sectionIntegr;

    static ID code;

    Vector e;          // trial section deformations
    Vector eCommit;    // committed section deformations 
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness


// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    Vector dedh; // MHS hack
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
