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
// $Date: 2008-08-26 16:47:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionAsym3d.h,v $
                                                                        
// Written: fmk
// Created: 04/01
//
// Description: This file contains the class definition for 
// FiberSectionAsym3d.h. FiberSectionAsym3d provides the abstraction of a 
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

// Modified by: Xinlong Du and Jerome F. Hajjar, Northeastern University, USA; Year 2019
// Description: Modified FiberSection3d.h (from version 3.0.0 on 11/5/2018)
//              to include shear center coordinates and high-order longitudinal strain terms.
// References:
// Du, X., & Hajjar, J. (2021). Three-dimensional nonlinear displacement-based beam element
// for members with angle and tee sections. Engineering Structures, 239, 112239.
// Du, X., & Hajjar, J. F. (2021). Three-dimensional nonlinear mixed 6-DOF beam element 
// for thin-walled members. Thin-Walled Structures, 164, 107817.

#ifndef FiberSectionAsym3d_h
#define FiberSectionAsym3d_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;
class Fiber;
class Response;
class SectionIntegration;

class FiberSectionAsym3d : public SectionForceDeformation
{
  public:
    FiberSectionAsym3d(); 
    FiberSectionAsym3d(int tag, int numFibers, Fiber **fibers,             //Xinlong
		   UniaxialMaterial *torsion = 0, double ys = 0.0, double zs = 0.0);                             //Xinlong
    FiberSectionAsym3d(int tag, int numFibers, UniaxialMaterial *torsion = 0, double ys=0.0, double zs=0.0); //Xinlong
    FiberSectionAsym3d(int tag, int numFibers, UniaxialMaterial **mats,
		   SectionIntegration &si, UniaxialMaterial *torsion = 0, double ys=0.0, double zs=0.0);         //Xinlong
    ~FiberSectionAsym3d();

    const char *getClassType(void) const {return "FiberSectionAsym3d";};

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

    const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    const Matrix & getSectionTangentSensitivity(int gradIndex);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const Vector & getSectionDeformationSensitivity(int gradIndex);
    // AddingSensitivity:END ///////////////////////////////////////////



  protected:
    
  private:
    int numFibers, sizeFibers;       // number of fibers in the section
    UniaxialMaterial **theMaterials; // array of pointers to materials
    double   *matData;               // data for the materials [yloc, zloc, area]
    double   kData[25];              // data for ks matrix       Xinlong 
    double   sData[5];               // data for s vector        Xinlong

    double QzBar, QyBar, Abar;
    double yBar;       // Section centroid
    double zBar;
	double ys; //Xinlong: y coord of shear center relative to centroid
	double zs; //Xinlong: z coord of shear center relative to centroid
  
    SectionIntegration *sectionIntegr;

    static ID code;

    Vector e;          // trial section deformations 
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness

    UniaxialMaterial *theTorsion;
};

#endif
