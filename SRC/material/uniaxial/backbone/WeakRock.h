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

#ifndef WeakRock_h
#define WeakRock_h

#include <HystereticBackbone.h>
#include <Vector.h>

/**
 * @brief WeakRock - on page 56 of this manual
 * https://www.dropbox.com/s/h3s9qvvz3pb4tlf/Weak%20Rock%20Reference%20Book.pdf?dl=0
 *
 */
class WeakRock : public HystereticBackbone {
 public:
  WeakRock(int tag, double kir, double pur, double yrm);
  WeakRock();
  ~WeakRock();

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
  double Kir;
  double pur;
  double yrm;
};

#endif
