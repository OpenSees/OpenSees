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

// Written: fmk                                                                         
                                                                        
#ifndef WrapperNDMaterial_h
#define WrapperNDMaterial_h

#include <NDMaterial.h>
#include <Matrix.h>

#include <elementAPI.h>

class WrapperNDMaterial : public NDMaterial
{
 public:
  // constructors  
  WrapperNDMaterial(const char *functName, matObject *theMat, int matType);
  WrapperNDMaterial();

  // destructor
  ~WrapperNDMaterial();

  int setTrialStrain (const Vector &strain);
  const Vector & getStrain (void);
  const Vector & getStress (void);
  const Matrix & getTangent (void);
  const Matrix & getInitialTangent (void);
  const Vector & getDampTangent (void);
  
  int commitState (void);
  int revertToLastCommit (void);    
  int revertToStart (void);        


  virtual const char *getType(void) const;  
  NDMaterial *getCopy (void);
  NDMaterial *getCopy (const char *code);

  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);
	       
 protected:
  
 private:
  char *funcName;
  matObject *theMat;
  int matType;

  int dataSize;
  double *data;
  Vector *strain;
  Vector *stress;
  Matrix *tangent;
  Matrix *initTangent;
};

#endif




