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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-08-26 16:48:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ParallelSection.h,v $
                                                                        
                                                                        
// File: ~/section/ParallelSection.h
//
// Written: MHS
// Created: Jun 2000
// Revision: A
//
// Description: This file contains the class definition for 
// ParallelSection.  ParallelSection decorates an MP
// section (couple bending and axial) with an uncoupled shear
// relation.
//
// What: "@(#) ParallelSection.h, revA"

#ifndef ParallelSection_h
#define ParallelSection_h

#include <SectionForceDeformation.h>
#include <UniaxialMaterial.h>

#include <Vector.h>
#include <Matrix.h>

class ParallelSection : public SectionForceDeformation
{
  public:
    ParallelSection(); 

    ParallelSection(int tag, int numSections, 
		    SectionForceDeformation **theSections);

    ~ParallelSection();

    const char *getClassType(void) const {return "ParallelSection";};

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
 
    void Print(OPS_Stream &s, int flag =0);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    const Vector & getStressResultantSensitivity(int gradIndex, bool 
conditional);
    const Vector & getSectionDeformationSensitivity(int gradIndex);
    const Matrix & getSectionTangentSensitivity(int gradIndex);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const Vector &getdedh(void); // MHS hack
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:

    int numSections;
    SectionForceDeformation **theSections;

    Vector *e;    // Storage for section deformations
    Vector *s;    // Storage for stress resultants
    Matrix *ks;   // Storage for section stiffness
    Matrix *fs;   // Storage for section flexibility

    int order;
    ID     *theCode;     // Storage for section type information

  void setUpCode(void);
   
    int otherDbTag;

    static double workArea[];
    static int codeArea[];

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;

    Vector dedh; // MHS hack
// AddingSensitivity:END ///////////////////////////////////////////

};

#endif

