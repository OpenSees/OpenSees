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
                                                                        
// $Revision: 1.10 $
// $Date: 2003-03-04 00:48:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection2d.h,v $
                                                                        
// Written: fmk
// Created: 04/01
//
// Description: This file contains the class definition for 
// FiberSection2d.h. FiberSection2d provides the abstraction of a 
// 2d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef FiberSection2d_h
#define FiberSection2d_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;
class Fiber;
class Response;

class FiberSection2d : public SectionForceDeformation
{
  public:
    FiberSection2d(); 
    FiberSection2d(int tag, int numFibers, Fiber **fibers); 
    ~FiberSection2d();

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
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
	    
    Response *setResponse(const char **argv, int argc, Information &info);
    int getResponse(int responseID, Information &info);

    int addFiber(Fiber &theFiber);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int   setParameter(const char **argv, int argc, Information &info);
    int   updateParameter(int parameterID, Information &info);
    int   activateParameter(int parameterID);
    const Vector & getStressResultantSensitivity(int gradNumber, bool conditional);
    const Vector & getSectionDeformationSensitivity(int gradNumber);
    const Matrix & getSectionTangentSensitivity(int gradNumber);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    int numFibers;                   // number of fibers in the section
    UniaxialMaterial **theMaterials; // array of pointers to materials
    double   *matData;               // data for the materials [yloc and area]
    double   kData[4];               // data for ks matrix 
    double   sData[2];               // data for s vector 
    
    double yBar;       // Section centroid
  
    static ID code;

    Vector e;          // trial section deformations 
    Vector eCommit;    // committed section deformations 
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
