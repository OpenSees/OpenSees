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

#include <CementedSoil.h>
#include <Channel.h>
#include <Vector.h>
#include <elementAPI.h>
#include <math.h>

/**
 * @brief OPS_API function to create a CementedSoil object
 * @param command CementedSoil tag pM pU Kpy z b
 * @param tag backbone tag
 * @param pM the soil resistance at the end of parabolic portion
 * @param pU the ultimate soil resistance
 * @param Kpy initial modulus of subgrade reaction
 * @param z depth below ground surface
 * @param b Diameter (width) of the pile
 *
 *
 * @return void*
 */
void *OPS_CementedSoil() {
  // check inputs
  if (OPS_GetNumRemainingInputArgs() < 6) {
    opserr << "WARNING: need hystereticBackbone CementedSoil "
           << "tag pM pU Kpy z b\n";
  }

  // get tag
  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING: invalid tag for hystereticBackbone "
              "CementedSoil\n";
    return 0;
  }

  // get pM pU Kpy z b
  double data[5];
  numData = 5;
  if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
    opserr << "WARNING: invalid data for hystereticBackbone "
              "CementedSoil\n";
    return 0;
  }

  // check data
  if (data[0] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "CementedSoil -- pM <= 0\n";
    return 0;
  }
  if (data[1] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "CementedSoil -- pU <= 0\n";
    return 0;
  }
  if (data[2] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "CementedSoil -- Kpy <= 0\n";
    return 0;
  }
  if (data[3] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "CementedSoil -- z <= 0\n";
    return 0;
  }
  if (data[4] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "CementedSoil -- b <= 0\n";
    return 0;
  }

  // create object
  CementedSoil *theBackbone =
      new CementedSoil(tag, data[0], data[1], data[2], data[3], data[4]);

  return theBackbone;
}

/**
 * @brief Construct a new CementedSoil object
 *
 * @param tag backbone tag
 * @param pM the soil resistance at the end of parabolic portion
 * @param pU the ultimate soil resistance
 * @param Kpy initial modulus of subgrade reaction
 * @param z depth below ground surface
 * @param b Diameter (width) of the pile
 */
CementedSoil::CementedSoil(int tag, double pM, double pU, double Kpy,
                           double z, double b)
    : HystereticBackbone(tag, BACKBONE_TAG_CementedSoil),
      pm(pM),
      pu(pU),
      kpy(Kpy),
      depth(z),
      diameter(b) {}

/**
 * @brief Default Construct
 *
 */
CementedSoil::CementedSoil()
    : HystereticBackbone(0, BACKBONE_TAG_CementedSoil),
      pm(0.0),
      pu(0.0),
      kpy(0.0),
      depth(0.0),
      diameter(0.0) {}

CementedSoil::~CementedSoil() {}

/**
 * @brief get tangent of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double CementedSoil::getTangent(double strain) {
  double ym = diameter / 60.0;
  double yu = diameter * 3.0 / 80.0;
  double m = (pu - pm) / (yu - ym);
  double n = pm / (m * ym);
  double C = pm / pow(ym, 1.0 / n);
  double yk = pow(C / (kpy * depth), n / (n - 1));

  if (strain < yk) {
    return kpy * depth;
  }

  if (strain < ym) {
    return C / n * pow(strain, (1 - n) / n);
  }

  if (strain < yu) {
    return m;
  }

  return 0.001 * kpy * depth;
}

/**
 * @brief get stress of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double CementedSoil::getStress(double strain) {
  double ym = diameter / 60.0;
  double yu = diameter * 3.0 / 80.0;
  double m = (pu - pm) / (yu - ym);
  double n = pm / (m * ym);
  double C = pm / pow(ym, 1.0 / n);
  double yk = pow(C / (kpy * depth), n / (n - 1));

  if (strain < yk) {
    return kpy * depth * strain;
  }

  if (strain < ym) {
    return C * pow(strain, 1.0 / n);
  }

  if (strain < yu) {
    return pm + m * (strain - ym);
  }

  return pu;
}

double CementedSoil::getEnergy(double strain) { return 0.0; }

double CementedSoil::getYieldStrain(void) { return 0.0; }

HystereticBackbone *CementedSoil::getCopy(void) {
  CementedSoil *theCopy =
      new CementedSoil(this->getTag(), pm, pu, kpy, depth, diameter);

  return theCopy;
}

void CementedSoil::Print(OPS_Stream &s, int flag) {
  s << "CementedSoil, tag: " << this->getTag() << endln;
  s << "\tpM: " << pm << endln;
  s << "\tpU: " << pu << endln;
  s << "\tKpy: " << kpy << endln;
  s << "\tz: " << depth << endln;
  s << "\tb: " << diameter << endln;
}

int CementedSoil::setVariable(char *argv) { return -1; }

int CementedSoil::getVariable(int varID, double &theValue) { return -1; }

int CementedSoil::sendSelf(int commitTag, Channel &theChannel) {
  int res = 0;

  static Vector data(6);

  data(0) = this->getTag();
  data(1) = pm;
  data(2) = pu;
  data(3) = kpy;
  data(4) = depth;
  data(5) = diameter;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CementedSoil::sendSelf -- could not send Vector" << endln;

    return res;
  }

  return res;
}

int CementedSoil::recvSelf(int commitTag, Channel &theChannel,
                           FEM_ObjectBroker &theBroker) {
  int res = 0;

  static Vector data(6);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CementedSoil::recvSelf -- could not receive Vector"
           << endln;

    return res;
  }

  this->setTag(int(data(0)));
  pm = data(1);
  pu = data(2);
  kpy = data(3);
  depth = data(4);
  diameter = data(5);

  return res;
}
