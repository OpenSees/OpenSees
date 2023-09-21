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

// $Revision: 1.6 $
// $Date: 2007-02-02 01:18:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionGJThermal.h,v $

// Written: fmk
// Created: 04/01
//
// Description: This file contains the class definition for
// FiberSectionGJThermal.h. FiberSectionGJThermal provides the abstraction of a
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef FiberSectionGJThermal_h
#define FiberSectionGJThermal_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;
class Fiber;
class Response;

class FiberSectionGJThermal : public SectionForceDeformation
{
 public:
  FiberSectionGJThermal();
  FiberSectionGJThermal(int tag, int numFibers, Fiber **fibers, double GJ = 1.0e10);
    FiberSectionGJThermal(int tag, int numFibers, double GJ = 1.0e10);
  ~FiberSectionGJThermal();

    const char *getClassType(void) const {return "FiberSectionGJThermal";};

  int  setTrialSectionDeformation(const Vector &deforms);

  const Vector &getSectionDeformation(void);

  const Vector &getTemperatureStress(const Vector& dataMixed); //UoE for obtaining thermal stress
  const Vector& getThermalElong(void);                       //Added by Liming (UoE)
  double determineFiberTemperature(const Vector& , double , double);   //Added by Liming (UoE)


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

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &info);

  int addFiber(Fiber &theFiber);

  int setParameter(const char **argv, int argc, Parameter &param);

 protected:

 private:
    int numFibers, sizeFibers;                   // number of fibers in the section
  UniaxialMaterial **theMaterials; // array of pointers to materials
  double *matData;               // data for the materials [yloc and area]
  double kData[16];               // data for ks matrix
  double sData[4];               // data for s vector

    double QzBar, QyBar, ABar;
  double yBar;       // Section centroid
  double zBar;

  Vector e;          // trial section deformations
  Vector eCommit;    // committed section deformations

  static ID code;
  Vector *s;         // section resisting forces
  Matrix *ks;        // section stiffness

  double GJ;

	Vector sT;  // JZ  section resisting forces, caused by the temperature
	Vector dataMixed;
	double *Fiber_ElongP;
	Vector AverageThermalElong;
};

#endif
