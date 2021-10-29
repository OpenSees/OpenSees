#include <stdlib.h>
#include <string.h>

#include <DegradingUniaxialWrapper.hh>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export

OPS_Export void *
OPS_DegradingUniaxialWrapper(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain = 1.0e16;
  int iData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial DegradingUniaxialWrapper $tag "
              "$otherTag <-min $minStrain> <-max $maxStrain>"
           << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial DegradingUniaxialWrapper $tag "
              "$otherTag"
           << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial MinMax tag: "
           << iData[0] << endln;
    return 0;
  }

  argc = OPS_GetNumRemainingInputArgs();
  while (argc > 1) {
    // char argvLoc[10];
    const char *argvLoc = OPS_GetString();
    /*    if (OPS_GetString(argvLoc, 10) != 0) {
      opserr << "WARNING invalid string option uniaxialMaterial MinMax tag: " <<
    iData[0] << endln; return 0;
    }
    */
    numData = 1;

    if ((strcmp(argvLoc, "-min") == 0) || (strcmp(argvLoc, "-Min") == 0) ||
        (strcmp(argvLoc, "-MIN") == 0)) {
      if (OPS_GetDouble(&numData, &minStrain) != 0) {
        opserr << "WARNING invalid min value  uniaxialMaterial MinMax tag: "
               << iData[0] << endln;
        return 0;
      }
    } else if ((strcmp(argvLoc, "-max") == 0) ||
               (strcmp(argvLoc, "-Max") == 0) ||
               (strcmp(argvLoc, "-MAX") == 0)) {
      if (OPS_GetDouble(&numData, &maxStrain) != 0) {
        opserr << "WARNING invalid min value  uniaxialMaterial MinMax tag: "
               << iData[0] << endln;
        return 0;
      }
    } else {
      opserr << "WARNING invalid option:" << argvLoc
             << " uniaxialMaterial MinMax tag: " << iData[0] << endln;
      return 0;
    }

    argc = OPS_GetNumRemainingInputArgs();
  }

  // Parsing was successful, allocate the material
  theMaterial = new DegradingUniaxialWrapper(iData[0], *theOtherMaterial,
                                             minStrain, maxStrain);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "DegradingUniaxialWrapper\n";
    return 0;
  }

  return theMaterial;
}

DegradingUniaxialWrapper::DegradingUniaxialWrapper(int tag,
                                                   UniaxialMaterial &material,
                                                   double min, double max)
    : UniaxialMaterial(tag, MAT_TAG_MinMax), theMaterial(0), minStrain(min),
      maxStrain(max), Tfailed(false), Cfailed(false)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr << "DegradingUniaxialWrapper::DegradingUniaxialWrapper -- failed to "
              "get copy of material\n";
    exit(-1);
  }
}

DegradingUniaxialWrapper::DegradingUniaxialWrapper()
    : UniaxialMaterial(0, MAT_TAG_MinMax), theMaterial(0), minStrain(0.0),
      maxStrain(0.0), Tfailed(false), Cfailed(false)
{
}

DegradingUniaxialWrapper::~DegradingUniaxialWrapper()
{
  if (theMaterial)
    delete theMaterial;
}

int
DegradingUniaxialWrapper::setTrialStrain(double strain, double strainRate)
{
  if (Cfailed)
    return 0;

  if (strain >= maxStrain || strain <= minStrain) {
    Tfailed = true;
    return 0;
  } else {
    Tfailed = false;
    return theMaterial->setTrialStrain(strain, strainRate);
  }
}

int
DegradingUniaxialWrapper::setTrialStrain(double strain, double temp,
                                         double strainRate)
{
  if (Cfailed)
    return 0;

  if (strain >= maxStrain || strain <= minStrain) {
    Tfailed = true;
    return 0;
  } else {
    Tfailed = false;
    return theMaterial->setTrialStrain(strain, temp, strainRate);
  }
}

double
DegradingUniaxialWrapper::getStress(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getStress();
}

double
DegradingUniaxialWrapper::getTangent(void)
{
  if (Tfailed)
    // return 0.0;
    return 1.0e-8 * theMaterial->getInitialTangent();
  else
    return theMaterial->getTangent();
}

double
DegradingUniaxialWrapper::getDampTangent(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}

double
DegradingUniaxialWrapper::getStrain(void)
{
  return theMaterial->getStrain();
}

double
DegradingUniaxialWrapper::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int
DegradingUniaxialWrapper::commitState(void)
{
  Cfailed = Tfailed;

  // Check if failed at current step
  if (Tfailed)
    return 0;
  else
    return theMaterial->commitState();
}

int
DegradingUniaxialWrapper::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int
DegradingUniaxialWrapper::revertToStart(void)
{
  Cfailed = false;
  Tfailed = false;

  return theMaterial->revertToStart();
}

UniaxialMaterial *
DegradingUniaxialWrapper::getCopy(void)
{
  DegradingUniaxialWrapper *theCopy = new DegradingUniaxialWrapper(
      this->getTag(), *theMaterial, minStrain, maxStrain);

  theCopy->Cfailed = Cfailed;
  theCopy->Tfailed = Tfailed;

  return theCopy;
}

int
DegradingUniaxialWrapper::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "DegradingUniaxialWrapper::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(3);
  dataVec(0) = minStrain;
  dataVec(1) = maxStrain;
  if (Cfailed == true)
    dataVec(2) = 1.0;
  else
    dataVec(2) = 0.0;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr
        << "DegradingUniaxialWrapper::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "DegradingUniaxialWrapper::sendSelf() - failed to send the "
              "Material\n";
    return -3;
  }

  return 0;
}

int
DegradingUniaxialWrapper::recvSelf(int cTag, Channel &theChannel,
                                   FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "DegradingUniaxialWrapper::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "DegradingUniaxialWrapper::recvSelf() - failed to create "
                "Material with classTag "
             << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(3);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr
        << "DegradingUniaxialWrapper::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  minStrain = dataVec(0);
  maxStrain = dataVec(1);

  if (dataVec(2) == 1.0)
    Cfailed = true;
  else
    Cfailed = false;

  Tfailed = Cfailed;

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "DegradingUniaxialWrapper::recvSelf() - failed to get the "
              "Material\n";
    return -4;
  }
  return 0;
}

void
DegradingUniaxialWrapper::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "DegradingUniaxialWrapper, tag: " << this->getTag() << endln;
    s << "  material: " << theMaterial->getTag() << endln;
    s << "  min strain: " << minStrain << endln;
    s << "  max strain: " << maxStrain << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"DegradingUniaxialWrapper\", ";
    s << "\"material\": \"" << theMaterial->getTag() << "\", ";
    s << "\"epsMin\": " << minStrain << ", ";
    s << "\"epsMax\": " << maxStrain << "}";
  }
}

int
DegradingUniaxialWrapper::setParameter(const char **argv, int argc,
                                       Parameter &param)
{
  //
  // I suppose epsMin and epsMax could be parameters, but not for now -- MHS
  //
  return theMaterial->setParameter(argv, argc, param);
}

int
DegradingUniaxialWrapper::updateParameter(int parameterID, Information &info)
{
  return 0;
}

double
DegradingUniaxialWrapper::getStressSensitivity(int gradIndex, bool conditional)
{
  if (Cfailed)
    return 0.0;
  else
    return theMaterial->getStressSensitivity(gradIndex, conditional);
}

double
DegradingUniaxialWrapper::getStrainSensitivity(int gradIndex)
{
  return theMaterial->getStrainSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getInitialTangentSensitivity(int gradIndex)
{
  return theMaterial->getInitialTangentSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getDampTangentSensitivity(int gradIndex)
{
  return theMaterial->getDampTangentSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getRhoSensitivity(int gradIndex)
{
  return theMaterial->getRhoSensitivity(gradIndex);
}

int
DegradingUniaxialWrapper::commitSensitivity(double strainGradient,
                                            int gradIndex, int numGrads)
{
  if (Cfailed)
    return 0;
  else
    return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
}
