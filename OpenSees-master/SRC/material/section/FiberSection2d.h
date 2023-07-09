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
                                                                        
// $Revision: 1.14 $
// $Date: 2008-08-26 16:47:26 $
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
class SectionIntegration;

class FiberSection2d : public SectionForceDeformation
{
  public:
    FiberSection2d(); 
    FiberSection2d(int tag, int numFibers, Fiber **fibers, bool compCentroid=true);
    FiberSection2d(int tag, int numFibers, bool compCentroid=true);
    FiberSection2d(int tag, int numFibers, UniaxialMaterial **mats,
		   SectionIntegration &si, bool compCentroid=true);
    ~FiberSection2d();

    const char *getClassType(void) const {return "FiberSection2d";};

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
	    
    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    int addFiber(Fiber &theFiber);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    const Vector& getStressResultantSensitivity(int gradIndex,
						bool conditional);
    const Vector& getSectionDeformationSensitivity(int gradIndex);
    const Matrix& getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(const Vector& sectionDeformationGradient,
			  int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////
	//by SAJalali
	double getEnergy() const;

  protected:
    
    //  private:
    int numFibers, sizeFibers;       // number of fibers in the section
    UniaxialMaterial **theMaterials; // array of pointers to materials
    double   *matData;               // data for the materials [yloc and area]
    double   kData[4];               // data for ks matrix 
    double   sData[2];               // data for s vector 
    
    double QzBar, ABar, yBar;       // Section centroid
    bool computeCentroid;
      
    SectionIntegration *sectionIntegr;

    static ID code;

    Vector e;          // trial section deformations 
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness

// AddingSensitivity:BEGIN //////////////////////////////////////////
    Vector dedh; // MHS hack
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
