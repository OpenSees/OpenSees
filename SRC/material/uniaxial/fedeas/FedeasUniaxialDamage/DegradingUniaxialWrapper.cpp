#include <functional> // std::hash
#include <string>

#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

#ifdef OPS_USE_RUNTIME
#  include <InputAPI.h>
#else
#  include <elementAPI.h>
#endif
#include "DegradingUniaxialWrapper.h"
#include "FedeasAPI.h"

#define MATERIAL_NAME "FedeasUniaxialDamage"
static const std::hash<std::string> hasher;
static const int MatTag = hasher(MATERIAL_NAME);

UniaxialMaterial *
TclModelBuilder_addFedeasUniaxialDamage(Tcl_Interp* interp, void *cd, int argc,
            TCL_Char **argv)
{
  // Pointer to the uniaxial material that will be returned
  DegradingUniaxialWrapper *theMaterial = nullptr;
  UniaxialMaterial *theWrappedMaterial  = nullptr;
  int tags[2];

  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial " MATERIAL_NAME " $tag "
              "$wrapTag <-damage $damageTag>"
           << endln;
    return nullptr;
  }

  // Get wrapper tag
  if (Tcl_GetInt(interp, argv[2], &tags[0]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return nullptr;
  }

  // Get base tag
  if (Tcl_GetInt(interp, argv[3], &tags[1]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return nullptr;
  }

  // Get instance of base material
#ifdef OPS_USE_RUNTIME
  theWrappedMaterial = G3_getUniaxialMaterialInstance(G3_getRuntime(interp), tags[1]);
#else
  theWrappedMaterial = OPS_GetUniaxialMaterial(tags[1]);
#endif
  if (theWrappedMaterial == 0) {
    opserr << "WARNING unable to retrieve uniaxialMaterial with tag"
              " tag: "
           << tags[1] << endln;
    return nullptr;
  }

  int argn = 4;
  const char *dmgtag = 0;
  double Ccd = 0.5;
  StateOperator *wrapper = new StateOperator;
  while (argn < argc) {
    const char *param = argv[argn];

    if ((strcmp(param, "-damage") == 0) || 
        (strcmp(param, "-dmg") == 0) ||
        (strcmp(param, "-DMG") == 0)) {
      *wrapper = *(StateOperator *)Tcl_GetAssocData(interp, 
                                                    "fedeas::damage::UniaxialDamage", NULL);

      wrapper->call(wrapper, interp, ISW_CREATE, argc - argn, &argv[++argn], 0, 0, 0, 0, 0);
      wrapper->call(wrapper, interp, ISW_MALLOC, 0, 0, 0, 0, 0, 0, 0);

    } else if ((strcmp(param, "-couple") == 0) ||
               (strcmp(param, "-ccd") == 0) || 
               (strcmp(param, "-Ccd") == 0)) {
      Ccd = std::stod(argv[++argn]);

    } else {
      break;
    }
    argn++;
  }

  // Parsing was successful, allocate the material
  theMaterial =
      new DegradingUniaxialWrapper(tags[0], *theWrappedMaterial, wrapper);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type " MATERIAL_NAME
           << endln;
    return nullptr;
  }
  theMaterial->setCoupling(Ccd);
  return theMaterial;
}

DegradingUniaxialWrapper::DegradingUniaxialWrapper(int tag,
                                                   UniaxialMaterial &material,
                                                   StateOperator *wrapper)
    : UniaxialMaterial(tag, MatTag), theMaterial(0), m_stress(0.0)
{
  theMaterial = material.getCopy();
  m_tangent = theMaterial->getInitialTangent();
  this->degrade = wrapper;

  if (theMaterial == 0) {
    opserr << MATERIAL_NAME "::DegradingUniaxialWrapper -- failed to "
              "get copy of material\n";
  }
}

DegradingUniaxialWrapper::DegradingUniaxialWrapper()
    : UniaxialMaterial(0, MatTag), theMaterial(0)
{
}

DegradingUniaxialWrapper::~DegradingUniaxialWrapper()
{
  if (theMaterial)
    delete theMaterial;
}

int
DegradingUniaxialWrapper::setCoupling(double coupling_param)
{
  return 0;
}

int
DegradingUniaxialWrapper::setTrialStrain(double strain, double temp,
                                         double strainRate)
{

  double pastStrain = theMaterial->getStrain();
  theMaterial->setTrialStrain(strain, temp, strainRate);
  double trialStrain = theMaterial->getStrain();

  double strain_incr = trialStrain - pastStrain;

  if (degrade) { //  && abs(strain_incr) > m_rate_tol){

    double inputs[] = {
        strain,
        theMaterial->getStress(),
        theMaterial->getTangent(),
    };
    double outputs[2];

    this->degrade->call(degrade, 0, ISW_UPDATE, 0, 0, 3, inputs, 2, outputs, 0);

    this->m_stress = outputs[0];
    this->m_tangent = outputs[1];

  } else {
    this->m_stress = theMaterial->getStress();
    this->m_tangent = theMaterial->getTangent();
  }
  return 0;
}

int
DegradingUniaxialWrapper::setTrialStrain(double trialStrain, double strainRate)
{
  return this->setTrialStrain(trialStrain, 0.0, strainRate);
}

double
DegradingUniaxialWrapper::getTangent(void)
{

  if (degrade)
    return this->m_tangent;
  else
    return theMaterial->getTangent();
}

double
DegradingUniaxialWrapper::getStress(void)
{
  if (degrade)
    return m_stress;
  else
    return theMaterial->getStress();
}

double
DegradingUniaxialWrapper::getInitialTangent(void)
{
  return theMaterial->getInitialTangent();
}

double
DegradingUniaxialWrapper::getDampTangent(void)
{
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
  degrade->call(degrade, 0, ISW_COMMIT, 0, 0, 0, 0, 0, 0, 0);
  return theMaterial->commitState();
}

int
DegradingUniaxialWrapper::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int
DegradingUniaxialWrapper::revertToStart(void)
{
  return theMaterial->revertToStart();
}

UniaxialMaterial *
DegradingUniaxialWrapper::getCopy(void)
{
  // TODO!! call copy routine for degradation
  StateOperator *damage = new StateOperator;
  *damage = *degrade;
  damage->call(damage, 0, ISW_MALLOC, 0, 0, 0, 0, 0, 0, 0);
  DegradingUniaxialWrapper *theCopy =
      new DegradingUniaxialWrapper(this->getTag(), *theMaterial, damage);

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
    opserr << MATERIAL_NAME "::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(3);
  // dataVec(0) = minStrain;
  // dataVec(1) = maxStrain;
  dataVec(2) = 0.0;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr
        << MATERIAL_NAME "::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << MATERIAL_NAME "::sendSelf() - failed to send the "
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
    opserr << MATERIAL_NAME "::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << MATERIAL_NAME "::recvSelf() - failed to create "
                "Material with classTag "
             << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(3);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr
        << MATERIAL_NAME "::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << MATERIAL_NAME "::recvSelf() - failed to get the "
              "Material\n";
    return -4;
  }
  return 0;
}

void
DegradingUniaxialWrapper::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << MATERIAL_NAME ", tag: " << this->getTag() << endln;
    s << "  material: " << theMaterial->getTag() << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"" MATERIAL_NAME "\", ";
    s << "\"material\": \"" << theMaterial->getTag() << "\", ";
  }
}

int
DegradingUniaxialWrapper::setParameter(const char **argv, int argc,
                                       Parameter &param)
{
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
  return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
}
