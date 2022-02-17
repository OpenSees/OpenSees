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

#ifndef CementedSoil_h
#define CementedSoil_h

#include <HystereticBackbone.h>
#include <Vector.h>

/**
 * @brief Cemented Soils - the Evans and Duncan (1982) SILT model at
 * http://www.findapile.com/p-y-curves/p-y-curves-models
 *
 */
class CementedSoil : public HystereticBackbone {
 public:
  CementedSoil(int tag, double pM, double pU, double Kpy,
               double z, double b);
  CementedSoil();
  ~CementedSoil();

  double getStress(double strain);
  double getTangent(double strain);
  double getEnergy(double strain);

  double getYieldStrain(void);

  HystereticBackbone *getCopy(void);

  void Print(OPS_Stream &s, int flag = 0);

  int setVariable(char *argv);
  int getVariable(int varID, double &theValue);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
               FEM_ObjectBroker &theBroker);

 protected:
 private:
  double pm;
  double pu;
  double kpy;
  double depth;
  double diameter;
};

#endif
