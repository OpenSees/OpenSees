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
                                                                        
#ifndef Beam3dPartialUniformLoad_h
#define Beam3dPartialUniformLoad_h

// Written: fmk 

// Purpose: This file contains the class definition for Beam3dPartialUniformLoad.

#include <ElementalLoad.h>

class Beam3dPartialUniformLoad : public ElementalLoad
{
 public:
  Beam3dPartialUniformLoad(int tag, double wTransya, double wTransza, double wAxiala, double aL, double bL, double wTransyb, double wTranszb, double wAxialb, int eleTag);
  Beam3dPartialUniformLoad();    
  ~Beam3dPartialUniformLoad();

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
  double wTransya;
  double wTransza;
  double wAxiala;
  double aOverL;
  double bOverL;
  double wTransyb;
  double wTranszb;
  double wAxialb;
  static Vector data;
  
  int parameterID;
};

#endif
