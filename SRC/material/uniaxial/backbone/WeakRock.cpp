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
#include <Vector.h>
#include <WeakRock.h>
#include <elementAPI.h>
#include <math.h>

/**
 * @brief OPS_API function to create a WeakRock object
 * @param command WeakRock tag Kir pur yrm
 * @param tag backbone tag
 * @param Kir initial slope of the curve
 * @param pur the rock mass ultimate resistance
 * @param yrm yrm = krm*B, B is the shaft diameter, krm is a constant
ranging from 0.0005 to 0.00005 that serves to establish the overall
stiffness of the curve.
 *
 *
 * @return void*
 */
void *OPS_WeakRock() {
  // check inputs
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING: need hystereticBackbone WeakRock "
           << "tag Kir pur yrm\n";
  }

  // get tag
  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING: invalid tag for hystereticBackbone "
              "WeakRock\n";
    return 0;
  }

  // get Kir pur yrm
  double data[3];
  numData = 3;
  if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
    opserr << "WARNING: invalid data for hystereticBackbone "
              "WeakRock\n";
    return 0;
  }

  // check data
  if (data[0] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "WeakRock -- Kir <= 0\n";
    return 0;
  }
  if (data[1] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "WeakRock -- pur <= 0\n";
    return 0;
  }
  if (data[2] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "WeakRock -- yrm <= 0\n";
    return 0;
  }

  // create object
  WeakRock *theBackbone = new WeakRock(tag, data[0], data[1], data[2]);

  return theBackbone;
}

/**
 * @brief Construct a new WeakRock object
 *
 * @param tag backbone tag
 * @param Kir initial slope of the curve
 * @param pur the rock mass ultimate resistance
 * @param yrm yrm = krm*B, B is the shaft diameter, krm is a constant
ranging from 0.0005 to 0.00005 that serves to establish the overall
stiffness of the curve.
 */
WeakRock::WeakRock(int tag, double k, double p, double y)
    : HystereticBackbone(tag, BACKBONE_TAG_WeakRock),
      Kir(k),
      pur(p),
      yrm(y) {}

/**
 * @brief Default Construct
 *
 */
WeakRock::WeakRock()
    : HystereticBackbone(0, BACKBONE_TAG_WeakRock),
      Kir(0.0),
      pur(0.0),
      yrm(0.0) {}

WeakRock::~WeakRock() {}

/**
 * @brief get tangent of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double WeakRock::getTangent(double strain) {
  double yA = pow(pur / (2 * pow(yrm, 0.25) * Kir), 1.333);
  double yu = 16.0 * yrm;

  if (strain < yA) {
    return Kir;
  }

  if (strain < yu) {
    return pur / (8 * yrm) * pow(strain / yrm, -0.75);
  }

  return 0.001 * Kir;
}

/**
 * @brief get stress of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double WeakRock::getStress(double strain) {
  double yA = pow(pur / (2 * pow(yrm, 0.25) * Kir), 1.333);
  double yu = 16.0 * yrm;

  if (strain < yA) {
    return Kir * strain;
  }

  if (strain < yu) {
    return pur / 2 * pow(strain / yrm, 0.25);
  }

  return pur;
}

double WeakRock::getEnergy(double strain) { return 0.0; }

double WeakRock::getYieldStrain(void) { return 0.0; }

HystereticBackbone *WeakRock::getCopy(void) {
  WeakRock *theCopy =
      new WeakRock(this->getTag(), Kir, pur, yrm);

  return theCopy;
}

void WeakRock::Print(OPS_Stream &s, int flag) {
  s << "WeakRock, tag: " << this->getTag() << endln;
  s << "\tKir: " << Kir << endln;
  s << "\tpur: " << pur << endln;
  s << "\tyrm: " << yrm << endln;
}

int WeakRock::setVariable(char *argv) { return -1; }

int WeakRock::getVariable(int varID, double &theValue) { return -1; }

int WeakRock::sendSelf(int commitTag, Channel &theChannel) {
  int res = 0;

  static Vector data(4);

  data(0) = this->getTag();
  data(1) = Kir;
  data(2) = pur;
  data(3) = yrm;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "WeakRock::sendSelf -- could not send Vector" << endln;

    return res;
  }

  return res;
}

int WeakRock::recvSelf(int commitTag, Channel &theChannel,
                       FEM_ObjectBroker &theBroker) {
  int res = 0;

  static Vector data(4);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "WeakRock::recvSelf -- could not receive Vector" << endln;

    return res;
  }

  this->setTag(int(data(0)));
  Kir = data(1);
  pur = data(2);
  yrm = data(3);

  return res;
}
