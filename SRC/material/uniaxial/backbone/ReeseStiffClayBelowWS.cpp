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
#include <ReeseStiffClayBelowWS.h>
#include <Vector.h>
#include <elementAPI.h>
#include <math.h>

/**
 * @brief OPS_API function to create a ReeseStiffClayBelowWS object
 * @param command ReeseStiffClayBelowWS tag Esi y50 As Pc
 * @param tag backbone tag
 * @param Esi Esi = ks*x, where ks is the k value for static load in
 * $lb/in^3$ and x is the depth.
 * @param y50 y50 = e50*b, where e50 is from results of laboratory tests
 * or, in the absence of laboratory tests, from Table 3.3, and b is the
 * pile diameter.
 * @param As the particular nondimensional depth from Fig. 3.7
 * @param Pc the ultimate soil resistance per unit length of pile, using
 * the smaller of the values given by eqs. 3.14 and 3.15.
 * @return void*
 */
void *OPS_ReeseStiffClayBelowWS() {
  // check inputs
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING: need hystereticBackbone ReeseStiffClayBelowWS "
           << "tag Esi y50 As Pc\n";
  }

  // get tag
  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING: invalid tag for hystereticBackbone "
              "ReeseStiffClayBelowWS\n";
    return 0;
  }

  // get Esi, y50, As, Pc
  double data[4];
  numData = 4;
  if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
    opserr << "WARNING: invalid data for hystereticBackbone "
              "ReeseStiffClayBelowWS\n";
    return 0;
  }

  // check data
  if (data[0] < 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "ReeseStiffClayBelowWS -- Esi < 0\n";
    return 0;
  }
  if (data[1] < 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "ReeseStiffClayBelowWS -- y50 < 0\n";
    return 0;
  }
  if (data[2] < 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "ReeseStiffClayBelowWS -- As < 0\n";
    return 0;
  }
  if (data[3] < 0.0) {
    opserr << "WARNING: hystereticBackbone "
              "ReeseStiffClayBelowWS -- Pc < 0\n";
    return 0;
  }

  // create object
  ReeseStiffClayBelowWS *theBackbone =
      new ReeseStiffClayBelowWS(tag, data[0], data[1], data[2], data[3]);

  return theBackbone;
}

/**
 * @brief Constructor
 *
 * @param tag backbone tag
 * @param Esi Esi = ks*x, where ks is the k value for static load in
 * $lb/in^3$ and x is the depth.
 * @param y50 y50 = e50*b, where e50 is from results of laboratory tests
 * or, in the absence of laboratory tests, from Table 3.3, and b is the
 * pile diameter.
 * @param As the particular nondimensional depth from Fig. 3.7
 * @param Pc the ultimate soil resistance per unit length of pile, using
 * the smaller of the values given by eqs. 3.14 and 3.15.
 */
ReeseStiffClayBelowWS::ReeseStiffClayBelowWS(int tag, double esi, double y,
                                             double as, double pc)
    : HystereticBackbone(tag, BACKBONE_TAG_ReeseStiffClayBelowWS),
      Esi(esi),
      y50(y),
      As(as),
      Pc(pc) {
  if (Esi < 0.0)
    opserr << "ReeseStiffClayBelowWS::ReeseStiffClayBelowWS -- Esi < 0"
           << endln;

  if (y50 < 0.0)
    opserr << "ReeseStiffClayBelowWS::ReeseStiffClayBelowWS -- y50 < 0"
           << endln;

  if (As < 0.0)
    opserr << "ReeseStiffClayBelowWS::ReeseStiffClayBelowWS -- As < 0"
           << endln;

  if (Pc < 0.0)
    opserr << "ReeseStiffClayBelowWS::ReeseStiffClayBelowWS -- Pc < 0"
           << endln;
}

/**
 * @brief Default Construct
 *
 */
ReeseStiffClayBelowWS::ReeseStiffClayBelowWS()
    : HystereticBackbone(0, BACKBONE_TAG_ReeseStiffClayBelowWS),
      Esi(0.0),
      y50(0.0),
      As(0.0),
      Pc(0.0) {}

ReeseStiffClayBelowWS::~ReeseStiffClayBelowWS() {}

/**
 * @brief get tangent of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double ReeseStiffClayBelowWS::getTangent(double strain) {
  strain = fabs(strain);

  double y1 = 0.25 * Pc * Pc / (y50 * Esi * Esi);
  double y2 = As * y50;

  if (strain <= y1)
    return Esi;

  else if (strain <= y2 && strain > y1)
    return 0.25 * Pc / y50 * pow(strain / y50, -0.5);

  else if (strain <= 6 * y2 && strain > y2)
    return 0.25 * Pc / y50 * pow(strain / y50, -0.5) -
           0.06875 * Pc / y2 * pow((strain - y2) / y2, 0.25);

  else if (strain <= 18 * y2 && strain > 6 * y2)
    return -0.0625 * Pc / y50;

  else
    return 0.001 * Esi;
}

/**
 * @brief get stress of the backbone function
 *
 * @param strain the strain value to be inquried
 * @return double
 */
double ReeseStiffClayBelowWS::getStress(double strain) {
  int signStrain = (strain > 0.0) ? 1 : -1;
  strain = signStrain * strain;

  double y1 = 0.25 * Pc * Pc / (y50 * Esi * Esi);
  double y2 = As * y50;

  double stress = 0.0;

  if (strain <= y1)
    stress = Esi * strain;

  else if (strain <= y2 && strain > y1)
    stress = 0.5 * Pc * pow(strain / y50, 0.5);

  else if (strain <= 6 * y2 && strain > y2)
    stress = 0.5 * Pc * pow(strain / y50, 0.5) -
             0.055 * Pc * pow((strain - y2) / y2, 1.25);

  else if (strain <= 18 * y2 && strain > 6 * y2)
    stress = 0.5 * Pc * pow(6 * As, 0.5) - 0.411 * Pc -
             0.0625 / y50 * Pc * (strain - 6 * y2);

  else if (strain > 18 * y2)
    stress = Pc * (1.225 * sqrt(As) - 0.75 * As - 0.411);

  return signStrain * stress;
}

double ReeseStiffClayBelowWS::getEnergy(double strain) { return 0.0; }

double ReeseStiffClayBelowWS::getYieldStrain(void) { return 0.0; }

HystereticBackbone *ReeseStiffClayBelowWS::getCopy(void) {
  ReeseStiffClayBelowWS *theCopy =
      new ReeseStiffClayBelowWS(this->getTag(), Esi, y50, As, Pc);

  return theCopy;
}

void ReeseStiffClayBelowWS::Print(OPS_Stream &s, int flag) {
  s << "ReeseStiffClayBelowWS, tag: " << this->getTag() << endln;
  s << "\tEsi: " << Esi << endln;
  s << "\ty50: " << y50 << endln;
  s << "\tAs: " << As << endln;
  s << "\tPc: " << Pc << endln;
}

int ReeseStiffClayBelowWS::setVariable(char *argv) { return -1; }

int ReeseStiffClayBelowWS::getVariable(int varID, double &theValue) {
  return -1;
}

int ReeseStiffClayBelowWS::sendSelf(int commitTag, Channel &theChannel) {
  int res = 0;

  static Vector data(5);

  data(0) = this->getTag();
  data(1) = Esi;
  data(2) = y50;
  data(3) = As;
  data(4) = Pc;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ReeseStiffClayBelowWS::sendSelf -- could not send Vector"
           << endln;

    return res;
  }

  return res;
}

int ReeseStiffClayBelowWS::recvSelf(int commitTag, Channel &theChannel,
                                    FEM_ObjectBroker &theBroker) {
  int res = 0;

  static Vector data(5);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ReeseStiffClayBelowWS::recvSelf -- could not receive Vector"
           << endln;

    return res;
  }

  this->setTag(int(data(0)));
  Esi = data(1);
  y50 = data(2);
  As = data(3);
  Pc = data(4);
  return res;
}
