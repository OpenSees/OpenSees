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
#include <VuggyLimestone.h>
#include <elementAPI.h>
#include <math.h>

/**
 * @brief OPS_API function to create a VuggyLimestone object
 * @param command VuggyLimestone tag b su
 * @param tag backbone tag
 * @param b the pile diameter
 * @param su the shear strength of the rock
 * @return void*
 */
void *OPS_VuggyLimestone() {
  // check inputs
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "WARNING: need hystereticBackbone VuggyLimestone "
           << "tag b su\n";
  }

  // get tag
  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING: invalid tag for hystereticBackbone "
              "VuggyLimestone\n";
    return 0;
  }

  // get b su
  double data[2];
  numData = 2;
  if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
    opserr << "WARNING: invalid data for hystereticBackbone "
              "VuggyLimestone\n";
    return 0;
  }

  // check data
  if (data[0] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "VuggyLimestone -- b <= 0\n";
    return 0;
  }
  if (data[1] <= 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "VuggyLimestone -- su <= 0\n";
    return 0;
  }

  // create object
  VuggyLimestone *theBackbone = new VuggyLimestone(tag, data[0], data[1]);

  return theBackbone;
}

/**
 * @brief Construct a new VuggyLimestone object
 *
 * @param tag backbone tag
 * @param b the pile diameter
 * @param su the shear strength of the rock
 */
VuggyLimestone::VuggyLimestone(int tag, double b, double su)
    : HystereticBackbone(tag, BACKBONE_TAG_VuggyLimestone),
      diameter(b),
      shearStrength(su) {}

/**
 * @brief Default Construct
 *
 */
VuggyLimestone::VuggyLimestone()
    : HystereticBackbone(0, BACKBONE_TAG_VuggyLimestone),
      diameter(0.0),
      shearStrength(0.0) {}

VuggyLimestone::~VuggyLimestone() {}

/**
 * @brief get tangent of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double VuggyLimestone::getTangent(double strain) {
  if (strain <= 0.0004 * diameter) {
    return 2000.0 * shearStrength;
  } else if (strain <= 0.0024 * diameter) {
    return 100.0 * shearStrength;
  }
  return shearStrength;
}

/**
 * @brief get stress of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double VuggyLimestone::getStress(double strain) {
  if (strain <= 0.0004 * diameter) {
    return 2000.0 * shearStrength * strain;
  } else if (strain <= 0.0024 * diameter) {
    return 0.8 * diameter * shearStrength +
           100.0 * shearStrength * (strain - 0.0004 * diameter);
  }
  return 0.0;
}

double VuggyLimestone::getEnergy(double strain) { return 0.0; }

double VuggyLimestone::getYieldStrain(void) { return 0.0; }

HystereticBackbone *VuggyLimestone::getCopy(void) {
  VuggyLimestone *theCopy =
      new VuggyLimestone(this->getTag(), diameter, shearStrength);

  return theCopy;
}

void VuggyLimestone::Print(OPS_Stream &s, int flag) {
  s << "VuggyLimestone, tag: " << this->getTag() << endln;
  s << "\tb: " << diameter << endln;
  s << "\tsu: " << shearStrength << endln;
}

int VuggyLimestone::setVariable(char *argv) { return -1; }

int VuggyLimestone::getVariable(int varID, double &theValue) { return -1; }

int VuggyLimestone::sendSelf(int commitTag, Channel &theChannel) {
  int res = 0;

  static Vector data(3);

  data(0) = this->getTag();
  data(1) = diameter;
  data(2) = shearStrength;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "VuggyLimestone::sendSelf -- could not send Vector" << endln;

    return res;
  }

  return res;
}

int VuggyLimestone::recvSelf(int commitTag, Channel &theChannel,
                             FEM_ObjectBroker &theBroker) {
  int res = 0;

  static Vector data(3);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "VuggyLimestone::recvSelf -- could not receive Vector"
           << endln;

    return res;
  }

  this->setTag(int(data(0)));
  diameter = data(1);
  shearStrength = data(2);

  return res;
}
