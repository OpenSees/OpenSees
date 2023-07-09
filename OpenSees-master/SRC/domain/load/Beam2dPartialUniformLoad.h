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
                                                                        
#ifndef Beam2dPartialUniformLoad_h
#define Beam2dPartialUniformLoad_h

// Written: fmk 

// Purpose: This file contains the class definition for Beam2dPartialUniformLoad.

#include <ElementalLoad.h>

class Beam2dPartialUniformLoad : public ElementalLoad
{
 public:
  Beam2dPartialUniformLoad(int tag, double wTrans, double wAxial, int eleTag);
  Beam2dPartialUniformLoad(int tag, double wTrans, double wAxial, double aL, double bL, int eleTag);
  Beam2dPartialUniformLoad(int tag, double wTransA, double wTransB, double wAxialA, double wAxialB, double aL, double bL, int eleTag);
  Beam2dPartialUniformLoad();    
  ~Beam2dPartialUniformLoad();

  const Vector &getData(int &type, double loadFactor);

  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int paramID);
  
  const Vector &getSensitivityData(int gradNumber);
  
 protected:
  
 private:
  double wTrans_a;
  double wTrans_b;
  double wAxial_a;
  double wAxial_b;
  double aOverL;
  double bOverL;
  static Vector data;
  
  int parameterID;
};

#endif
