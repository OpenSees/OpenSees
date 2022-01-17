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

#ifndef ReeseStiffClayAboveWS_h
#define ReeseStiffClayAboveWS_h

#include <HystereticBackbone.h>
#include <Vector.h>

/**
 * @brief Response of Stiff Clay above the water surface
 * (https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/PB94108305.xhtml)
 * page 336
 *
 */
class ReeseStiffClayAboveWS : public HystereticBackbone {
 public:
  ReeseStiffClayAboveWS(int tag, double pu, double y50);
  ReeseStiffClayAboveWS();
  ~ReeseStiffClayAboveWS();

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
  double pu;
  double y50;
  double hl;
};

#endif
