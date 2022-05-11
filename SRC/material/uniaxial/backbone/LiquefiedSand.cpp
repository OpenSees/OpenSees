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

#include <Channel.h>
#include <LiquefiedSand.h>
#include <Vector.h>
#include <elementAPI.h>
#include <math.h>

/**
 * @brief OPS_API function to create a LiquefiedSand object
 * @param command LiquefiedSand tag X D
 * @param tag backbone tag
 * @param X the depth
 * @param D the pile diameter between 0.3m and 2.6m
 * @param kN define the unit kilo Newton
 * @param m define the unit meter
 *
 * @return void*
 */
void *OPS_LiquefiedSand() {
  // check inputs
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING: need hystereticBackbone LiquefiedSand "
           << "tag X D kN m\n";
  }

  // get tag
  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING: invalid tag for hystereticBackbone "
              "LiquefiedSand\n";
    return 0;
  }

  // get X and D, kN, m
  double data[4];
  numData = 4;
  if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
    opserr << "WARNING: invalid data for hystereticBackbone "
              "LiquefiedSand\n";
    return 0;
  }

  // check data
  if (data[0] < 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "LiquefiedSand -- X < 0\n";
    return 0;
  }
  if (data[2] < 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "LiquefiedSand -- kN < 0\n";
    return 0;
  }
  if (data[3] < 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "LiquefiedSand -- m < 0\n";
    return 0;
  }
  double m = data[3];
  if (data[1] < 0.3 * m) {
    opserr << "WARNING: hystereticBackbone "
              "LiquefiedSand -- D < 0.3 m\n";
    return 0;
  }
  if (data[1] > 2.6 * m) {
    opserr << "WARNING: hystereticBackbone "
              "LiquefiedSand -- D > 2.6 m\n";
    return 0;
  }

  // create object
  LiquefiedSand *theBackbone =
      new LiquefiedSand(tag, data[0], data[1], data[2], data[3]);

  return theBackbone;
}

/**
 * @brief Construct a new LiquefiedSand object
 *
 * @param tag backbone tag
 * @param X the depth
 * @param D the pile diameter between 0.3m and 2.6m
 * @param kN define the unit kilo Newton
 * @param m define the unit meter
 */
LiquefiedSand::LiquefiedSand(int tag, double x, double d, double kn,
                             double m)
    : HystereticBackbone(tag, BACKBONE_TAG_LiquefiedSand),
      X(x),
      D(d),
      kN(kn),
      meter(m),
      yu(0.15 * m) {}

/**
 * @brief Default Construct
 *
 */
LiquefiedSand::LiquefiedSand()
    : HystereticBackbone(0, BACKBONE_TAG_LiquefiedSand),
      X(0.0),
      D(0.0),
      kN(1.0),
      meter(1.0),
      yu(0.15) {}

LiquefiedSand::~LiquefiedSand() {}

/**
 * @brief get tangent of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double LiquefiedSand::getTangent(double strain) {
  double A = 3e-7 * pow(X + 1, 6.05);
  double B = 2.8 * pow(X + 1, 0.11);
  double C = 2.85 * pow(X + 1, -0.41);
  double Pd = 3.81 * log(D) + 5.6;

  double k0 = Pd * A * B * C * pow(B * 0.001 * yu, C - 1);

  if (strain < 0.001 * yu) {
    return k0;
  }

  if (strain < yu) {
    return Pd * A * B * C * pow(B * strain, C - 1);
  }

  return 0.001 * k0;
}

/**
 * @brief get stress of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double LiquefiedSand::getStress(double strain) {
  double A = 3e-7 * pow(X + 1, 6.05);
  double B = 2.8 * pow(X + 1, 0.11);
  double C = 2.85 * pow(X + 1, -0.41);
  double Pd = 3.81 * log(D) + 5.6;

  if (strain < 0.001 * yu) {
    double k0 = Pd * A * B * C * pow(B * 0.001 * yu, C - 1);
    return k0 * strain;
  }

  if (strain < yu) {
    double P03m = A * pow(B * strain, C);
    if (P03m > 15 * kN / meter) {
      opserr << "WARNING: P0.3m > 15 kN/m\n";
    }
    return Pd * P03m;
  }

  double P03m = A * pow(B * yu, C);
  if (P03m > 15 * kN / meter) {
    opserr << "WARNING: P0.3m > 15 kN/m\n";
  }
  return Pd * P03m;
}

double LiquefiedSand::getEnergy(double strain) { return 0.0; }

double LiquefiedSand::getYieldStrain(void) { return 0.0; }

HystereticBackbone *LiquefiedSand::getCopy(void) {
  LiquefiedSand *theCopy =
      new LiquefiedSand(this->getTag(), X, D, kN, meter);

  return theCopy;
}

void LiquefiedSand::Print(OPS_Stream &s, int flag) {
  s << "LiquefiedSand, tag: " << this->getTag() << endln;
  s << "\tX: " << X << endln;
  s << "\tD: " << D << endln;
  s << "\tkN: " << kN << endln;
  s << "\tmeter: " << meter << endln;
}

int LiquefiedSand::setVariable(char *argv) { return -1; }

int LiquefiedSand::getVariable(int varID, double &theValue) { return -1; }

int LiquefiedSand::sendSelf(int commitTag, Channel &theChannel) {
  int res = 0;

  static Vector data(5);

  data(0) = this->getTag();
  data(1) = X;
  data(2) = D;
  data(3) = kN;
  data(4) = meter;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "LiquefiedSand::sendSelf -- could not send Vector" << endln;

    return res;
  }

  return res;
}

int LiquefiedSand::recvSelf(int commitTag, Channel &theChannel,
                            FEM_ObjectBroker &theBroker) {
  int res = 0;

  static Vector data(5);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "LiquefiedSand::recvSelf -- could not receive Vector"
           << endln;

    return res;
  }

  this->setTag(int(data(0)));
  X = data(1);
  D = data(2);
  kN = data(3);
  meter = data(4);

  return res;
}
