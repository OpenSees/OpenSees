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
#include <ReeseStiffClayAboveWS.h>
#include <Vector.h>
#include <elementAPI.h>
#include <math.h>

/**
 * @brief OPS_API function to create a ReeseStiffClayAboveWS object
 * @param command ReeseStiffClayAboveWS tag pu y50
 * @param tag backbone tag
 * @param pu the ultimate soil resistance per unit length of shaft
 * @param y50 the deflection at one-half the ultimate soil resistance
 * @return void*
 */
void *OPS_ReeseStiffClayAboveWS() {
  // check inputs
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING: need hystereticBackbone ReeseStiffClayAboveWS "
           << "tag pu y50\n";
  }

  // get tag
  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING: invalid tag for hystereticBackbone "
              "ReeseStiffClayAboveWS\n";
    return 0;
  }

  // get pu y50
  double data[2];
  numData = 2;
  if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
    opserr << "WARNING: invalid data for hystereticBackbone "
              "ReeseStiffClayAboveWS\n";
    return 0;
  }

  // check data
  if (data[0] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "ReeseStiffClayAboveWS -- pu <= 0\n";
    return 0;
  }
  if (data[1] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "ReeseStiffClayAboveWS -- y50 <= 0\n";
    return 0;
  }

  // create object
  ReeseStiffClayAboveWS *theBackbone =
      new ReeseStiffClayAboveWS(tag, data[0], data[1]);

  return theBackbone;
}

/**
 * @brief Construct a new ReeseStiffClayAboveWS object
 *
 * @param tag backbone tag
 * @param pu the ultimate soil resistance per unit length of shaft
 * @param y50 the deflection at one-half the ultimate soil resistance
 */
ReeseStiffClayAboveWS::ReeseStiffClayAboveWS(int tag, double pu,
                                             double y50)
    : HystereticBackbone(tag, BACKBONE_TAG_ReeseStiffClayAboveWS),
      pu(pu),
      y50(y50),
      hl(0.001) {}

/**
 * @brief Default Construct
 *
 */
ReeseStiffClayAboveWS::ReeseStiffClayAboveWS()
    : HystereticBackbone(0, BACKBONE_TAG_ReeseStiffClayAboveWS),
      pu(0.0),
      y50(0.0),
      hl(0.001) {}

ReeseStiffClayAboveWS::~ReeseStiffClayAboveWS() {}

/**
 * @brief get tangent of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double ReeseStiffClayAboveWS::getTangent(double strain) {
  double yhl = hl * y50;
  double k0 = getStress(yhl) / yhl;
  if (strain < yhl) {
    return k0;
  }

  if (strain > 16.0 * y50) {
    return hl * k0;
  }

  return pu * pow(strain / y50, -0.75) / 8.0 / y50;
}

/**
 * @brief get stress of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double ReeseStiffClayAboveWS::getStress(double strain) {
  double yhl = hl * y50;
  if (strain < yhl) {
    return strain * getStress(yhl) / yhl;
  }

  if (strain > 16.0 * y50) {
    return pu;
  }

  return 0.5 * pu * pow(strain / y50, 0.25);
}

double ReeseStiffClayAboveWS::getEnergy(double strain) { return 0.0; }

double ReeseStiffClayAboveWS::getYieldStrain(void) { return 0.0; }

HystereticBackbone *ReeseStiffClayAboveWS::getCopy(void) {
  ReeseStiffClayAboveWS *theCopy =
      new ReeseStiffClayAboveWS(this->getTag(), pu, y50);

  return theCopy;
}

void ReeseStiffClayAboveWS::Print(OPS_Stream &s, int flag) {
  s << "ReeseStiffClayAboveWS, tag: " << this->getTag() << endln;
  s << "\tpu: " << pu << endln;
  s << "\ty50: " << y50 << endln;
}

int ReeseStiffClayAboveWS::setVariable(char *argv) { return -1; }

int ReeseStiffClayAboveWS::getVariable(int varID, double &theValue) {
  return -1;
}

int ReeseStiffClayAboveWS::sendSelf(int commitTag, Channel &theChannel) {
  int res = 0;

  static Vector data(3);

  data(0) = this->getTag();
  data(1) = pu;
  data(2) = y50;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ReeseStiffClayAboveWS::sendSelf -- could not send Vector"
           << endln;

    return res;
  }

  return res;
}

int ReeseStiffClayAboveWS::recvSelf(int commitTag, Channel &theChannel,
                                    FEM_ObjectBroker &theBroker) {
  int res = 0;

  static Vector data(3);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ReeseStiffClayAboveWS::recvSelf -- could not receive Vector"
           << endln;

    return res;
  }

  this->setTag(int(data(0)));
  pu = data(1);
  y50 = data(2);

  return res;
}
