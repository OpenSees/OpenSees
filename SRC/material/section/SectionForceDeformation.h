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
// $Date: 2008-08-26 16:45:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/SectionForceDeformation.h,v $
                                                                        
                                                                        
#ifndef SectionForceDeformation_h
#define SectionForceDeformation_h

// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for SectionForceDeformation.
// SectionForceDeformation is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) SectionForceDeformation.h, revA"

#include <Material.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Information;
class Response;

#define MAX_SECTION_RESPONSE_ID 10000

#define SECTION_RESPONSE_NONE   -1

#define SECTION_RESPONSE_MZ		1
#define SECTION_RESPONSE_P		2
#define SECTION_RESPONSE_VY		3
#define SECTION_RESPONSE_MY		4
#define SECTION_RESPONSE_VZ		5
#define SECTION_RESPONSE_T		6	
#define SECTION_RESPONSE_R		7	
#define SECTION_RESPONSE_Q		8	
#define SECTION_RESPONSE_B              9 // Bi-moment (FiberSectionWarping3d)
#define SECTION_RESPONSE_W             10 // (FiberSectionWarping3d)

// section responses for shells
#define SECTION_RESPONSE_FXX 11 // membrane xx
#define SECTION_RESPONSE_FYY 12 // membrane yy
#define SECTION_RESPONSE_FXY 13 // membrane xy
#define SECTION_RESPONSE_MXX 14 // bending xx
#define SECTION_RESPONSE_MYY 15 // bending yy
#define SECTION_RESPONSE_MXY 16 // bending xy
#define SECTION_RESPONSE_VXZ 17 // bending yy
#define SECTION_RESPONSE_VYZ 18 // bending xy


class SectionForceDeformation : public Material
{
 public:
  SectionForceDeformation (int tag, int classTag);
  SectionForceDeformation ();
  virtual ~SectionForceDeformation ();
  
  virtual int setTrialSectionDeformation (const Vector&) = 0;
  virtual const Vector &getSectionDeformation (void) = 0;
  
  virtual const Vector &getStressResultant (void) = 0;
  virtual const Matrix &getSectionTangent (void) = 0;
  virtual const Matrix &getInitialTangent (void) = 0;
  virtual const Matrix &getSectionFlexibility (void);
  virtual const Matrix &getInitialFlexibility (void);
  
  virtual double getRho(void);
  
  virtual int commitState (void) = 0;
  virtual int revertToLastCommit (void) = 0;
  virtual int revertToStart (void) = 0;
  
  virtual SectionForceDeformation *getCopy (void) = 0;
  virtual const ID &getType (void) = 0;
  virtual int getOrder (void) const = 0;
  
  virtual Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  virtual int getResponse(int responseID, Information &info);

  virtual int getResponseSensitivity(int responseID, int gradIndex,
				     Information &info);
  
  

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  virtual const Vector &getStressResultantSensitivity(int gradIndex,
						      bool conditional);
  virtual const Vector &getSectionDeformationSensitivity(int gradIndex);
  virtual const Matrix &getSectionTangentSensitivity(int gradIndex);
  virtual const Matrix &getSectionFlexibilitySensitivity(int gradIndex);
  virtual const Matrix &getInitialTangentSensitivity(int gradIndex);
  virtual const Matrix &getInitialFlexibilitySensitivity(int gradIndex);
  virtual double getRhoSensitivity(int gradIndex);
  virtual int commitSensitivity(const Vector& sectionDeformationGradient,
				int gradIndex, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////
  
  //--- Adding Thermal Materials:[BEGIN]   by UoE OpenSees Group ----//  
  virtual int setTrialSectionDeformation(const Vector&, const Vector &); //JZ
  virtual const Vector &getTemperatureStress(const Vector &tData);//27 is for 'FireLoadPattern'
  //--- Adding Thermal Functions:[END]   by UoE OpenSees Group ----//

  virtual const Vector& getThermalElong(void);
  virtual double getEnergy() const { return 0; };		//by SAJalali

 protected:
  Matrix *fDefault;	// Default flexibility matrix
  Vector *sDefault;
  
 private:

};

extern bool OPS_addSectionForceDeformation(SectionForceDeformation *newComponent);
extern SectionForceDeformation *OPS_getSectionForceDeformation(int tag);
extern bool OPS_removeSectionForceDeformation(int tag);
extern void OPS_clearAllSectionForceDeformation(void);
extern void OPS_printSectionForceDeformation(OPS_Stream &s, int flag=0);

#endif
