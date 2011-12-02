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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-12-09 19:50:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/WrapperUniaxialMaterial.h,v $

// Written: fmk                                                                         
                                                                        
#ifndef WrapperUniaxialMaterial_h
#define WrapperUniaxialMaterial_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

#include <elementAPI.h>

class WrapperUniaxialMaterial : public UniaxialMaterial
{
 public:
  // constructors  
  WrapperUniaxialMaterial(const char *functName, matObject *theMat);
  WrapperUniaxialMaterial();

  // destructor
  ~WrapperUniaxialMaterial();

  int setTrialStrain (double strain, double strainRate = 0.0);
  int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
  double getStrain (void);
  double getStress (void);
  double getTangent (void);
  double getInitialTangent (void);
  double getDampTangent (void);
  
  int commitState (void);
  int revertToLastCommit (void);    
  int revertToStart (void);        
  
  UniaxialMaterial *getCopy (void);


  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);
	       
 protected:
  
 private:
  char *funcName;
  matObject *theMat;

  double strain;
  double stress;
  double tangent;
  double initTangent;
};

#endif




